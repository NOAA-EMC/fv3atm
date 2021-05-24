
# fv3atm

This is the NOAA EMC version of the FV3 atmospheric model dynamical
core, originally from GFDL.

## Authors

Not sure what to list here...

## Prerequisites

This package requires the following packages:
 - [netcdf-c Library](https://github.com/Unidata/netcdf-c)
 - [netcdf-fortran Library](https://github.com/Unidata/netcdf-fortran)
 - [ESMF](https://github.com/esmf-org/esmf)
 - [Jasper](https://github.com/jasper-software/jasper)

## Building and Installing

```
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/path/to/install ..
make -j2
make install
```


## References

## Disclaimer

The United States Department of Commerce (DOC) GitHub project code is
provided on an "as is" basis and the user assumes responsibility for
its use. DOC has relinquished control of the information and no longer
has responsibility to protect the integrity, confidentiality, or
availability of the information. Any claims against the Department of
Commerce stemming from the use of its GitHub project will be governed
by all applicable Federal law. Any reference to specific commercial
products, processes, or services by service mark, trademark,
manufacturer, or otherwise, does not constitute or imply their
endorsement, recommendation or favoring by the Department of
Commerce. The Department of Commerce seal and logo, or the seal and
logo of a DOC bureau, shall not be used in any manner to imply
endorsement of any commercial product or activity by DOC or the United
States Government.

