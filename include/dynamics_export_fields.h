#undef _atm_src_gfs
#undef _atm_src_dycore
#undef _atm_src_array
#undef _atm_src_cnv
#undef _atm_src_ptr
#undef _atm_src_trac
#undef _atm_dst_ptr
#undef _atm_copy_lev
#undef _atm_copy_3d

          !--- Dycore quantities

#define _atm_src_dycore
#define _atm_dst_ptr datar82d

          ! bottom layer temperature (t)
          case('inst_temp_height_lowest')
#define _atm_src_ptr coupling%t_bot
#include "gfs_data_copy.h"

          ! bottom layer specific humidity (q)
          !    !    ! CHECK if tracer 1 is for specific humidity     !    !    !
          case('inst_spec_humid_height_lowest')
#define _atm_src_ptr coupling%tr_bot
#define _atm_src_trac 1
#include "gfs_data_copy.h"

          ! bottom layer zonal wind (u)
          case('inst_zonal_wind_height_lowest')
#define _atm_src_ptr coupling%u_bot
#include "gfs_data_copy.h"

          ! bottom layer meridionalw wind (v)
          case('inst_merid_wind_height_lowest')
#define _atm_src_ptr coupling%v_bot
#include "gfs_data_copy.h"

          ! bottom layer pressure (p)
          case('inst_pres_height_lowest')
#define _atm_src_ptr coupling%p_bot
#include "gfs_data_copy.h"

          ! bottom layer height (z)
          case('inst_height_lowest')
#define _atm_src_ptr coupling%z_bot
#include "gfs_data_copy.h"
