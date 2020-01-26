module module_write_netcdf_parallel

  use esmf
  implicit none
  private
  public write_netcdf_parallel

  contains

!----------------------------------------------------------------------------------------
  subroutine write_netcdf_parallel(fieldbundle, wrtfb, filename, mpi_comm, mype, im, jm, rc)
!
    type(ESMF_FieldBundle), intent(in) :: fieldbundle
    type(ESMF_FieldBundle), intent(in) :: wrtfb
    character(*), intent(in)           :: filename
    integer, intent(in)                :: mpi_comm
    integer, intent(in)                :: mype
    integer, intent(in)                :: im, jm
    integer, optional,intent(out)      :: rc

    print *,'in stub write_netcdf_parallel - model not built with parallel netcdf support, return'

  end subroutine write_netcdf_parallel

 
end module module_write_netcdf_parallel
