#undef _atm_src_dycore
#undef _atm_src_gfs
#undef _atm_src_array
#undef _atm_src_ptr
#undef _atm_src_trac
#undef _atm_dst_ptr
#undef _atm_src_cnv
#undef _atm_copy_lev

          !--- Instantaneous quantities

#define _atm_src_gfs
#define _atm_dst_ptr datar82d

          ! Instantaneous u wind (m/s) 10 m above ground
          case ('inst_zonal_wind_height10m')
#define _atm_src_ptr coupling%u10mi_cpl
#include "gfs_data_copy.h"

          ! Instantaneous v wind (m/s) 10 m above ground
          case ('inst_merid_wind_height10m')
#define _atm_src_ptr coupling%v10mi_cpl
#include "gfs_data_copy.h"

          ! Instantaneous Zonal compt of momentum flux (N/m**2)
          case ('inst_zonal_moment_flx')
#define _atm_src_ptr coupling%dusfci_cpl
#include "gfs_data_copy.h"

          ! Instantaneous Merid compt of momentum flux (N/m**2)
          case ('inst_merid_moment_flx')
#define _atm_src_ptr coupling%dvsfci_cpl
#include "gfs_data_copy.h"

          ! Instantaneous Sensible heat flux (W/m**2)
          case ('inst_sensi_heat_flx')
#define _atm_src_ptr coupling%dtsfci_cpl
#include "gfs_data_copy.h"

          ! Instantaneous Latent heat flux (W/m**2)
          case ('inst_laten_heat_flx')
#define _atm_src_ptr coupling%dqsfci_cpl
#include "gfs_data_copy.h"

          ! Instantaneous Downward long wave radiation flux (W/m**2)
          case ('inst_down_lw_flx')
#define _atm_src_ptr coupling%dlwsfci_cpl
#include "gfs_data_copy.h"

          ! Instantaneous Downward solar radiation flux (W/m**2)
          case ('inst_down_sw_flx')
#define _atm_src_ptr coupling%dswsfci_cpl
#include "gfs_data_copy.h"

          ! Instantaneous Temperature (K) 2 m above ground
          case ('inst_temp_height2m')
#define _atm_src_ptr coupling%t2mi_cpl
#include "gfs_data_copy.h"

          ! Instantaneous Specific humidity (kg/kg) 2 m above ground
          case ('inst_spec_humid_height2m')
#define _atm_src_ptr coupling%q2mi_cpl
#include "gfs_data_copy.h"

          ! Instantaneous Temperature (K) at surface
          case ('inst_temp_height_surface')
#define _atm_src_ptr coupling%tsfci_cpl
#include "gfs_data_copy.h"

          ! Instantaneous Pressure (Pa) land and sea surface
          case ('inst_pres_height_surface')
#define _atm_src_ptr coupling%psurfi_cpl
#include "gfs_data_copy.h"

          ! Instantaneous Surface height (m)
          case ('inst_surface_height')
#define _atm_src_ptr coupling%oro_cpl
#include "gfs_data_copy.h"

          ! Instantaneous NET long wave radiation flux (W/m**2)
          case ('inst_net_lw_flx')
#define _atm_src_ptr coupling%nlwsfci_cpl
#include "gfs_data_copy.h"

          ! Instantaneous NET solar radiation flux over the ocean (W/m**2)
          case ('inst_net_sw_flx')
#define _atm_src_ptr coupling%nswsfci_cpl
#include "gfs_data_copy.h"

          ! Instantaneous sfc downward nir direct flux (W/m**2)
          case ('inst_down_sw_ir_dir_flx')
#define _atm_src_ptr coupling%dnirbmi_cpl
#include "gfs_data_copy.h"

          ! Instantaneous sfc downward nir diffused flux (W/m**2)
          case ('inst_down_sw_ir_dif_flx')
#define _atm_src_ptr coupling%dnirdfi_cpl
#include "gfs_data_copy.h"

          ! Instantaneous sfc downward uv+vis direct flux (W/m**2)
          case ('inst_down_sw_vis_dir_flx')
#define _atm_src_ptr coupling%dvisbmi_cpl
#include "gfs_data_copy.h"

          ! Instantaneous sfc downward uv+vis diffused flux (W/m**2)
          case ('inst_down_sw_vis_dif_flx')
#define _atm_src_ptr coupling%dvisdfi_cpl
#include "gfs_data_copy.h"

          ! Instantaneous net sfc nir direct flux (W/m**2)
          case ('inst_net_sw_ir_dir_flx')
#define _atm_src_ptr coupling%nnirbmi_cpl
#include "gfs_data_copy.h"

          ! Instantaneous net sfc nir diffused flux (W/m**2)
          case ('inst_net_sw_ir_dif_flx')
#define _atm_src_ptr coupling%nnirdfi_cpl
#include "gfs_data_copy.h"

          ! Instantaneous net sfc uv+vis direct flux (W/m**2)
          case ('inst_net_sw_vis_dir_flx')
#define _atm_src_ptr coupling%nvisbmi_cpl
#include "gfs_data_copy.h"

          ! Instantaneous net sfc uv+vis diffused flux (W/m**2)
          case ('inst_net_sw_vis_dif_flx')
#define _atm_src_ptr coupling%nvisdfi_cpl
#include "gfs_data_copy.h"

          ! Land/Sea mask (sea:0,land:1)
          case ('inst_land_sea_mask')
#define _atm_src_ptr coupling%slmsk_cpl
#include "gfs_data_copy.h"

          !--- Mean quantities

#define _atm_src_cnv rtime

          ! MEAN Zonal compt of momentum flux (N/m**2)
          case ('mean_zonal_moment_flx_atm')
#define _atm_src_ptr coupling%dusfc_cpl
#include "gfs_data_copy.h"

          ! MEAN Merid compt of momentum flux (N/m**2)
          case ('mean_merid_moment_flx_atm')
#define _atm_src_ptr coupling%dvsfc_cpl
#include "gfs_data_copy.h"

          ! MEAN Sensible heat flux (W/m**2)
          case ('mean_sensi_heat_flx')
#define _atm_src_ptr coupling%dtsfc_cpl
#include "gfs_data_copy.h"

          ! MEAN Latent heat flux (W/m**2)
          case ('mean_laten_heat_flx')
#define _atm_src_ptr coupling%dqsfc_cpl
#include "gfs_data_copy.h"

          ! MEAN Downward LW heat flux (W/m**2)
          case ('mean_down_lw_flx')
#define _atm_src_ptr coupling%dlwsfc_cpl
#include "gfs_data_copy.h"

          ! MEAN Downward SW heat flux (W/m**2)
          case ('mean_down_sw_flx')
#define _atm_src_ptr coupling%dswsfc_cpl
#include "gfs_data_copy.h"

          ! MEAN NET long wave radiation flux (W/m**2)
          case ('mean_net_lw_flx')
#define _atm_src_ptr coupling%nlwsfc_cpl
#include "gfs_data_copy.h"

          ! MEAN NET solar radiation flux over the ocean (W/m**2)
          case ('mean_net_sw_flx')
#define _atm_src_ptr coupling%nswsfc_cpl
#include "gfs_data_copy.h"

          ! MEAN sfc downward nir direct flux (W/m**2)
          case ('mean_down_sw_ir_dir_flx')
#define _atm_src_ptr coupling%dnirbm_cpl
#include "gfs_data_copy.h"

          ! MEAN sfc downward nir diffused flux (W/m**2)
          case ('mean_down_sw_ir_dif_flx')
#define _atm_src_ptr coupling%dnirdf_cpl
#include "gfs_data_copy.h"

          ! MEAN sfc downward uv+vis direct flux (W/m**2)
          case ('mean_down_sw_vis_dir_flx')
#define _atm_src_ptr coupling%dvisbm_cpl
#include "gfs_data_copy.h"

          ! MEAN sfc downward uv+vis diffused flux (W/m**2)
          case ('mean_down_sw_vis_dif_flx')
#define _atm_src_ptr coupling%dvisdf_cpl
#include "gfs_data_copy.h"

          ! MEAN NET sfc nir direct flux (W/m**2)
          case ('mean_net_sw_ir_dir_flx')
#define _atm_src_ptr coupling%nnirbm_cpl
#include "gfs_data_copy.h"

          ! MEAN NET sfc nir diffused flux (W/m**2)
          case ('mean_net_sw_ir_dif_flx')
#define _atm_src_ptr coupling%nnirdf_cpl
#include "gfs_data_copy.h"

          ! MEAN NET sfc uv+vis direct flux (W/m**2)
          case ('mean_net_sw_vis_dir_flx')
#define _atm_src_ptr coupling%nvisbm_cpl
#include "gfs_data_copy.h"

          ! MEAN NET sfc uv+vis diffused flux (W/m**2)
          case ('mean_net_sw_vis_dif_flx')
#define _atm_src_ptr coupling%nvisdf_cpl
#include "gfs_data_copy.h"

#undef  _atm_src_cnv
#define _atm_src_cnv rtimek

          ! MEAN precipitation rate (kg/m2/s)
          case ('mean_prec_rate')
#define _atm_src_ptr coupling%rain_cpl
#include "gfs_data_copy.h"

          ! MEAN snow precipitation rate (kg/m2/s)
          case ('mean_fprec_rate')
#define _atm_src_ptr coupling%snow_cpl
#include "gfs_data_copy.h"

          ! oceanfrac used by atm to calculate fluxes
          case ('openwater_frac_in_atm')
#undef  _atm_src_cnv
#define _atm_src_cnv (one-GFS_Data(nb)%sfcprop%fice(ix))
#define _atm_src_ptr sfcprop%oceanfrac
          if (associated(GFS_Data(nsb)%sfcprop%fice)) then
#include "gfs_data_copy.h"
          endif
