#undef _atm_src_dycore
#undef _atm_src_gfs
#undef _atm_src_array
#undef _atm_src_ptr
#undef _atm_src_trac
#undef _atm_dst_ptr
#undef _atm_src_cnv
#undef _atm_copy_lev
#undef _atm_copy_3d

          !--- JEDI fields

#define _atm_dst_ptr datar83d
#define _atm_copy_3d

          case ('u')
#define _atm_src_ptr u
#include "atm_data_copy.h"

          case ('v')
#define _atm_src_ptr v
#include "atm_data_copy.h"

          case ('ua')
#define _atm_src_ptr ua
#include "atm_data_copy.h"

          case ('va')
#define _atm_src_ptr va
#include "atm_data_copy.h"

          case ('t')
#define _atm_src_ptr pt
#include "atm_data_copy.h"

          case ('delp')
#define _atm_src_ptr delp
#include "atm_data_copy.h"

#define _atm_src_trac sphum
          case ('sphum')
            sphum = get_tracer_index(MODEL_ATMOS, 'sphum')
#define _atm_src_ptr q
#include "atm_data_copy.h"

#define _atm_src_trac ice_wat
          case ('ice_wat')
            ice_wat = get_tracer_index(MODEL_ATMOS, 'ice_wat')
#define _atm_src_ptr q
#include "atm_data_copy.h"

#define _atm_src_trac liq_wat
          case ('liq_wat')
            liq_wat = get_tracer_index(MODEL_ATMOS, 'liq_wat')
#define _atm_src_ptr q
#include "atm_data_copy.h"

#define _atm_src_trac o3mr
          case ('o3mr')
            o3mr = get_tracer_index(MODEL_ATMOS, 'o3mr')
#define _atm_src_ptr q
#include "atm_data_copy.h"

#undef _atm_copy_3d
#undef _atm_dst_ptr
#define _atm_dst_ptr datar82d

          case ('phis')
#define _atm_src_ptr phis
#include "atm_data_copy.h"

          case ('u_srf')
#define _atm_src_ptr u_srf
#include "atm_data_copy.h"

          case ('v_srf')
#define _atm_src_ptr v_srf
#include "atm_data_copy.h"

          case ('slmsk')
#define _atm_src_gfs
#define _atm_src_ptr sfcprop%slmsk
#include "gfs_data_copy.h"

          case ('weasd')
#define _atm_src_ptr sfcprop%weasd
#include "gfs_data_copy.h"

          case ('tsea')
#define _atm_src_ptr sfcprop%tsfco
#include "gfs_data_copy.h"

          case ('vtype')
#define _atm_src_ptr sfcprop%vtype
#include "gfs_data_copy.h"

          case ('stype')
#define _atm_src_ptr sfcprop%stype
#include "gfs_data_copy.h"

          case ('vfrac')
#define _atm_src_ptr sfcprop%vfrac
#include "gfs_data_copy.h"

#undef  _atm_dst_ptr
#define _atm_dst_ptr datar83d

          case ('stc')
#define _atm_src_ptr sfcprop%stc
#define _atm_copy_lev
#include "gfs_data_copy.h"

          case ('smc')
#define _atm_src_ptr sfcprop%smc
#define _atm_copy_lev
#include "gfs_data_copy.h"

#undef  _atm_dst_ptr
#define _atm_dst_ptr datar82d

          case ('snwdph')
#define _atm_src_ptr sfcprop%snowd
#include "gfs_data_copy.h"

          case ('f10m')
#define _atm_src_ptr sfcprop%f10m
#include "gfs_data_copy.h"

          case ('zorl')
#define _atm_src_ptr sfcprop%zorl
#include "gfs_data_copy.h"

          case ('t2m')
#define _atm_src_ptr sfcprop%t2m
#include "gfs_data_copy.h"
