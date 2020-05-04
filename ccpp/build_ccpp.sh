#!/bin/bash

set +x
set -eu

# List of valid/tested machines
VALID_MACHINES=( wcoss_cray wcoss_dell_p3 gaea.intel jet.intel \
                 hera.intel hera.gnu orion.intel \
                 cheyenne.intel cheyenne.intel-impi cheyenne.gnu cheyenne.pgi endeavor.intel \
                 stampede.intel supermuc_phase2.intel macosx.gnu \
                 linux.intel linux.gnu linux.pgi )

###################################################################################################

function usage   {
  echo "Usage: "
  echo "build_ccpp.sh MACHINE_ID CCPP_DIR CCPP_MK [ 'MAKE_OPT' ] [ clean_before ] [ clean_after ]"
  echo "    Where: MACHINE      [required] can be : ${VALID_MACHINES[@]}"
  echo "           CCPP_DIR     [required] is the target installation directory for CCPP"
  echo "           CCPP_MK      [required] is the location/name of the CCPP ESMF makefile fragment"
  echo "           MAKE_OPT     [optional] can be any of the NEMSfv3gfs MAKE_OPT options,"
  echo "                                   enclosed in a single string; used:"
  echo "                                   SION=Y/N        (default N)"
  echo "                                   DEBUG=Y/N       (default N)"
  echo "                                   REPRO=Y/N       (default N)"
  echo "                                   OPENMP=Y/N      (default Y)"
  echo "                                   32BIT=Y/N       (default N, affects dynamics/fast physics only)"
  echo "                                   SUITES=ABC,XYZ  (comma-separated list of CCPP suites; "
  echo "                                                    corresponding filenames: suite_ABC.xml. ...)"
  echo "                                   MULTI_GASES=Y/N (default N)"
  echo "           clean_before [optional] can be 'YES' (default) or 'NO'"
  echo "           clean_after  [optional] can be 'YES' (default) or 'NO'"
  exit 1
}

function checkvalid {
# Ensure value ($2) of variable ($1) is contained in list of validvalues ($3)
  if [ $# -lt 3 ]; then
    echo $FUNCNAME requires at least 3 arguments: stopping
    exit 1
  fi

  var_name=$1 && shift
  input_val=$1 && shift
  valid_vars=($@)

  for x in ${valid_vars[@]}; do
    if [ "$input_val" == "$x" ]; then
      echo "${var_name}=${input_val} is valid."
      return
    fi
  done
  echo "ERROR: ${var_name}=${input_val} is invalid."
  usage
  exit 1
}

function trim {
    local var="$1"
    # remove leading whitespace characters
    var="${var#"${var%%[![:space:]]*}"}"
    # remove trailing whitespace characters
    var="${var%"${var##*[![:space:]]}"}"   
    echo -n "$var"
}

###################################################################################################

# Check and process command line arguments

if [[ $# -lt 2 ]]; then usage; fi

readonly MACHINE_ID=$1
readonly CCPP_DIR=$2
readonly CCPP_MK=$3
readonly MAKE_OPT=${4:-}
readonly clean_before=${5:-YES}
readonly clean_after=${6:-YES}

checkvalid MACHINE_ID $MACHINE_ID ${VALID_MACHINES[@]}

# Set compilers for cmake
source ./set_compilers.sh

# Obtain ESMF_LIB from ESMFMKFILE's ESMF_LIBSDIR entry
readonly ESMF_LIB=$(cat $ESMFMKFILE | grep -E '^ESMF_LIBSDIR=.+' | cut -d = -f 2)
echo "Obtained ESMF_LIB=${ESMF_LIB} from ${ESMFMKFILE}"

# Account for inconsistencies in HPC modules: if environment variable
# NETCDF is undefined, try to set from NETCDF_DIR, NETCDF_ROOT, ...
if [[ "${MACHINE_ID}" == "wcoss_cray" || "${MACHINE_ID}" == "gaea.intel" ]]; then
  NETCDF=${NETCDF:-${NETCDF_DIR}}
fi

# Generate CCPP cmake flags from MAKE_OPT
CCPP_CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=${CCPP_DIR} -DNETCDF_DIR=${NETCDF} -DMPI=ON"
CCPP_MAKE_FLAGS=""
if [[ "${MAKE_OPT}" == *"SION=Y"* ]]; then
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DSIONLIB=${SIONLIB}"
fi
if [[ "${MAKE_OPT}" == *"DEBUG=Y"* ]]; then
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DCMAKE_BUILD_TYPE=Debug"
  CCPP_MAKE_FLAGS="${CCPP_MAKE_FLAGS} VERBOSE=1"
elif [[ "${MAKE_OPT}" == *"REPRO=Y"* ]]; then
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DCMAKE_BUILD_TYPE=Bitforbit"
  CCPP_MAKE_FLAGS="${CCPP_MAKE_FLAGS} VERBOSE=1"
else
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DCMAKE_BUILD_TYPE=Release"
  # Don't use the AVX512 flags yet on hera
  #if [[ "${MACHINE_ID}" == "jet.intel" || "${MACHINE_ID}" == "hera.intel" ]]; then
  if [[ "${MACHINE_ID}" == "jet.intel" ]]; then
    CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DSIMDMULTIARCH=ON"
  fi
  CCPP_MAKE_FLAGS="${CCPP_MAKE_FLAGS} VERBOSE=1"
fi
if [[ "${MAKE_OPT}" == *"OPENMP=N"* ]]; then
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DOPENMP=OFF"
else
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DOPENMP=ON"
fi
if [[ "${MAKE_OPT}" == *"32BIT=Y"* ]]; then
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DDYN32=ON"
else
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DDYN32=OFF"
fi
if [[ "${MAKE_OPT}" == *"MULTI_GASES=Y"* ]]; then
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DMULTI_GASES=ON"
else
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DMULTI_GASES=OFF"
fi

# Flag to cmake that modern Intel compilers are used
CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DLEGACY_INTEL=OFF"

# Generate additional CCPP cmake flags depending on machine / compiler
if [[ "${MACHINE_ID}" == "macosx.gnu" ]]; then
  # Intel MKL (for FFTW)
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DMKL_DIR=${MKL_DIR}"
  # ESMF (for DGEMM) - replace with MKL version in the future?
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DESMF_LIB_DIR=${ESMF_LIB}"
  # netCDF (needed when linking ESMF)
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DNETCDF_DIR=${NETCDF}"
elif [[ "${MACHINE_ID}" == "linux.gnu" ]]; then
  # ESMF (for DGEMM) - replace with MKL version in the future?
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DESMF_LIB_DIR=${ESMF_LIB}"
  # netCDF (needed when linking ESMF)
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DNETCDF_DIR=${NETCDF}"
elif [[ "${MACHINE_ID}" == "gaea.intel" || "${MACHINE_ID}" == "wcoss_cray" ]]; then
  # Fix broken libxml2 installation on gaea by
  # using own version (not known to the system)
  if [[ "${MACHINE_ID}" == "gaea.intel" ]]; then
    CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DLIBXML2_LIB_DIR=${LIBXML2_LIB_DIR} -DLIBXML2_INCLUDE_DIR=${LIBXML2_INCLUDE_DIR}"
  fi
fi

CCPP_CMAKE_FLAGS=$(trim "${CCPP_CMAKE_FLAGS}")
CCPP_MAKE_FLAGS=$(trim "${CCPP_MAKE_FLAGS}")

# Build and install CCPP

echo "Building CCPP with options '${CCPP_CMAKE_FLAGS}' ..."
PATH_CCPP=${PWD}
PATH_CCPP_BUILD=${PWD}/build
PATH_CCPP_INC=${PWD}/include
PATH_CCPP_LIB=${PWD}/lib

if [ $clean_before = YES ]; then
    rm -fr ${PATH_CCPP_BUILD}
    rm -fr ${PATH_CCPP_INC}
    rm -fr ${PATH_CCPP_LIB}
    rm -f ${CCPP_MK}
fi
mkdir -p ${PATH_CCPP_BUILD}
cd ${PATH_CCPP_BUILD}
cmake ${CCPP_CMAKE_FLAGS} ${PATH_CCPP}
make ${CCPP_MAKE_FLAGS}
make ${CCPP_MAKE_FLAGS} install

# Generate ESMF makefile fragment

set -u

# Set linker flags
CCPP_LINK_OBJS="-L${PATH_CCPP_LIB} -lccpp -lccppphys"

echo "ESMF_DEP_INCPATH=${PATH_CCPP_INC} ${PATH_CCPP_BUILD}/physics" > ${CCPP_MK}
echo "ESMF_DEP_LINK_OBJS=${CCPP_LINK_OBJS}" >> ${CCPP_MK}

if [ $clean_after = YES ]; then
    rm -fr ${PATH_CCPP_BUILD}
fi

