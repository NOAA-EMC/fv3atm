#undef __atm_data_src
#ifdef _atm_src_ptr
#  define __atm_data_src Atm(mygrid)%_atm_src_ptr
#endif

#undef __atm_data_dst
#ifdef _atm_dst_ptr
#  define __atm_data_dst _atm_dst_ptr
#endif 

#undef  __ik
#undef  __beg_k_loop
#undef  __end_k_loop
#ifdef _atm_copy_3d
#  define __ik ,k
#  define __beg_k_loop do k=1,nk
#  define __end_k_loop enddo
#else
#  define __ik
#  define __beg_k_loop
#  define __end_k_loop
#endif

#undef __lev
#ifdef _atm_copy_lev
#  define __lev ,:
#else
#  define __lev
#endif

#undef __trac
#undef __itrac
#undef __beg_check_trac
#undef __end_check_trac
#ifdef _atm_src_trac
#  define __trac _atm_src_trac
#  define __itrac ,__trac
#  define __beg_check_trac if (__trac > 0) then
#  define __end_check_trac endif
#else
#  define __trac
#  define __itrac
#  define __beg_check_trac
#  define __end_check_trac
#endif

#if (defined(__atm_data_dst) && defined(__atm_data_src))
            __beg_check_trac
!$omp parallel do default(shared) private(i,ib,ix,j,jb __ik)
            __beg_k_loop
              do jb=jsc,jec
                i = 0
                j = jb - jsc + 1
                do ib=isc,iec
                  i = i + 1
                  __atm_data_dst(i, j __ik __lev) = __atm_data_src(i, j __ik __itrac __lev)
                enddo
              enddo
            __end_k_loop
            __end_check_trac
#endif

#undef _atm_src_ptr
#undef _atm_src_trac
