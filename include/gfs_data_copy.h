#undef  __atm_data_type
#if defined(_atm_src_gfs)
#  define __atm_data_type GFS_data
#elif defined(_atm_src_dycore)
#  define __atm_data_type DYCORE_data
#else
#  define __atm_data_type GFS_data
#endif

#undef  __atm_data_src
#undef  __atm_data_zero
#undef  __beg_check_ptr
#undef  __end_check_ptr
#if defined(_atm_src_ptr)
#  define __atm_data_src(N) __atm_data_type(N)%_atm_src_ptr
#  define __beg_check_ptr if (associated(__atm_data_src(nsb))) then
#  define __end_check_ptr endif
#  ifdef _atm_src_dycore
#    define __atm_data_zero
#  endif
#elif defined(_atm_src_array)
#  define __atm_data_src(N) __atm_data_type(N)%_atm_src_array
#  define __beg_check_ptr
#  define __end_check_ptr
#endif 

#undef __atm_data_dst
#ifdef _atm_dst_ptr
#  define __atm_data_dst _atm_dst_ptr
#endif 

#undef __atm_data_cnv
#ifdef _atm_src_cnv
#  define __atm_data_cnv *_atm_src_cnv
#else
#  define __atm_data_cnv
#endif

#undef __lev
#ifdef _atm_copy_lev
#  define __lev ,:
#else
#  define __lev
#endif

#undef __itrac
#ifdef _atm_src_trac
#  define __itrac ,_atm_src_trac
#else
#  define __itrac
#endif

#if (defined(__atm_data_dst) && defined(__atm_data_src))
            __beg_check_ptr
!$omp parallel do default(shared) private(i,ib,ix,j,jb,nb)
              do jb=jsc,jec
                i = 0
                j = jb - jsc + 1
                do ib=isc,iec
                  nb = Atm_block%blkno(ib,jb)
                  ix = Atm_block%ixp(ib,jb)
                  i = i + 1
                  __atm_data_dst(i,j __lev) = __atm_data_src(nb)(ix __lev __itrac)__atm_data_cnv
                enddo
              enddo
#  ifdef __atm_data_zero
            else
!$omp parallel do default(shared) private(i,ib,j,jb)
              do jb=jsc,jec
                i = 0
                j = jb - jsc + 1
                do ib=isc,iec
                  i = i + 1
                  __atm_data_dst(i,j __lev) = zero
                enddo
              enddo
#  endif
            __end_check_ptr
#endif

#undef _atm_src_array
#undef _atm_src_ptr
#undef _atm_src_trac
#undef _atm_copy_lev
