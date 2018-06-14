
MODULE module_mp_wsm6
!
!
   REAL, PARAMETER, PRIVATE :: dtcldcr     = 120.
   REAL, PARAMETER, PRIVATE :: n0r = 8.e6
   REAL, PARAMETER, PRIVATE :: n0g = 4.e6
   REAL, PARAMETER, PRIVATE :: avtr = 841.9
   REAL, PARAMETER, PRIVATE :: bvtr = 0.8
   REAL, PARAMETER, PRIVATE :: r0 = .8e-5 ! 8 microm  in contrast to 10 micro m
   REAL, PARAMETER, PRIVATE :: peaut = .55   ! collection efficiency
   REAL, PARAMETER, PRIVATE :: xncr = 3.e8   ! maritime cloud in contrast to 3.e8 in tc80
   REAL, PARAMETER, PRIVATE :: xmyu = 1.718e-5 ! the dynamic viscosity kgm-1s-1
   REAL, PARAMETER, PRIVATE :: avts = 11.72
   REAL, PARAMETER, PRIVATE :: bvts = .41
   REAL, PARAMETER, PRIVATE :: avtg = 330.
   REAL, PARAMETER, PRIVATE :: bvtg = 0.8
   REAL, PARAMETER, PRIVATE :: deng = 500.
   REAL, PARAMETER, PRIVATE :: n0smax =  1.e11 ! t=-90C unlimited
   REAL, PARAMETER, PRIVATE :: lamdarmax = 8.e4
   REAL, PARAMETER, PRIVATE :: lamdasmax = 1.e5
   REAL, PARAMETER, PRIVATE :: lamdagmax = 6.e4
   REAL, PARAMETER, PRIVATE :: betai = .6
   REAL, PARAMETER, PRIVATE :: xn0 = 1.e-2
   REAL, PARAMETER, PRIVATE :: dicon = 11.9
   REAL, PARAMETER, PRIVATE :: di0 = 12.9e-6
   REAL, PARAMETER, PRIVATE :: dimax = 500.e-6
   REAL, PARAMETER, PRIVATE :: n0s = 2.e6             ! temperature dependent n0s
   REAL, PARAMETER, PRIVATE :: alpha = .12        ! .122 exponen factor for n0s
   REAL, PARAMETER, PRIVATE :: pfrz1 = 100.
   REAL, PARAMETER, PRIVATE :: pfrz2 = 0.66
   REAL, PARAMETER, PRIVATE :: qcrmin = 1.e-9
   REAL, PARAMETER, PRIVATE :: t40c = 233.16
   REAL, PARAMETER, PRIVATE :: eacrc = 1.0
   REAL, PARAMETER, PRIVATE :: dens  =  100.0
   REAL, PARAMETER, PRIVATE :: qs0   =  6.e-4   ! pgaut
   REAL, SAVE ::                                     &
             qc0, qck1,bvtr1,bvtr2,bvtr3,bvtr4,g1pbr,&
             g3pbr,g4pbr,g5pbro2,pvtr,eacrr,pacrr,   &
             bvtr6,g6pbr,                            &
             precr1,precr2,xm0,xmmax,roqimax,bvts1,  &
             bvts2,bvts3,bvts4,g1pbs,g3pbs,g4pbs,    &
             g5pbso2,pvts,pacrs,precs1,precs2,pidn0r,&
             pidn0s,xlv1,pacrc,                      &
             bvtg1,bvtg2,bvtg3,bvtg4,g1pbg,          &
             g3pbg,g4pbg,g5pbgo2,pvtg,pacrg,         &
             precg1,precg2,pidn0g,                   &
             rslopermax,rslopesmax,rslopegmax,       &
             rsloperbmax,rslopesbmax,rslopegbmax,    &
             rsloper2max,rslopes2max,rslopeg2max,    &
             rsloper3max,rslopes3max,rslopeg3max
CONTAINS
!===================================================================
!
  SUBROUTINE wsm6(th, q, qc, qr, qi, qs, qg                        &
                 ,den, pii, p, delz                                &
                 ,delt,g, cpd, cpv, rd, rv, t0c                    &
                 ,ep1, ep2, qmin                                   &
                 ,XLS, XLV0, XLF0, den0, denr                      &
                 ,cliq,cice,psat                                   &
                 ,rain, rainncv                                    &
                 ,snow, snowncv                                    &
                 ,graupel, graupelncv                              &
                 ,sr                                               &
                 ,ims,ime, jms,jme, lm                             &
                 ,d_ss,mprates                                     &
                                                                   )
!-------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------
!
!  This code is a 6-class GRAUPEL phase microphyiscs scheme (WSM6) of the WRF
!  Single-Moment MicroPhyiscs (WSMMP). The WSMMP assumes that ice nuclei
!  number concentration is a function of temperature, and seperate assumption
!  is developed, in which ice crystal number concentration is a function
!  of ice amount. A theoretical background of the ice-microphysics and related
!  processes in the WSMMPs are described in Hong et al. (2004).
!  All production terms in the WSM6 scheme are described in Hong and Lim (2006).
!  All units are in m.k.s. and source/sink terms in kgkg-1s-1.
!
!  WSM6 cloud scheme
!
!  Coded by Song-You Hong and Jeong-Ock Jade Lim (Yonsei Univ.)
!           Summer 2003
!
!  Implemented by Song-You Hong (Yonsei Univ.) and Jimy Dudhia (NCAR)
!           Summer 2004
!
!  Reference) Hong, Dudhia, Chen (HDC, 2004) Mon. Wea. Rev. 
!             Hong and Lim (HL, 2006) J. Korean Meteor. Soc. 
!             Lin, Farley, Orville (LFO, 1983) J. Appl. Meteor.
!             Rutledge, Hobbs (RH83, 1983) J. Atmos. Sci.
!             Rutledge, Hobbs (RH84, 1984) J. Atmos. Sci.
!
  INTEGER,      INTENT(IN   )    ::   ims,ime, jms,jme, lm
  REAL, DIMENSION( ims:ime , jms:jme, lm ),                       &
        INTENT(INOUT) ::                                          &
                                                             th,  &
                                                              q,  &
                                                              qc, &
                                                              qi, &
                                                              qr, &
                                                              qs, &
                                                              qg
  REAL, DIMENSION( ims:ime , jms:jme, lm ),                       &
        INTENT(IN   ) ::                                          &
                                                             den, &
                                                             pii, &
                                                               p, &
                                                            delz
  REAL, INTENT(IN   ) ::                                    delt, &
                                                               g, &
                                                              rd, &
                                                              rv, &
                                                             t0c, &
                                                            den0, &
                                                             cpd, &
                                                             cpv, &
                                                             ep1, &
                                                             ep2, &
                                                            qmin, &
                                                             XLS, &
                                                            XLV0, &
                                                            XLF0, &
                                                            cliq, &
                                                            cice, &
                                                            psat, &
                                                            denr
  REAL, DIMENSION( ims:ime , jms:jme ),                           &
        INTENT(INOUT) ::                                    rain, &
                                                         rainncv, &
                                                              sr

  REAL, DIMENSION( ims:ime , jms:jme ), OPTIONAL,                &
        INTENT(INOUT) ::                                    snow, &
                                                         snowncv

  REAL, DIMENSION( ims:ime , jms:jme ), OPTIONAL,                &
        INTENT(INOUT) ::                                graupel, &
                                                        graupelncv

     INTEGER :: d_ss
     REAL,               DIMENSION(ims:ime, jms:jme,lm,d_ss) ::  &
                          mprates 
! LOCAL VAR
  REAL, DIMENSION( ims:ime , lm ) ::   t
  REAL, DIMENSION( ims:ime , lm ) ::   q2,den2,p2,delz2
  REAL, DIMENSION( ims:ime , lm, 2 ) ::   qci
  REAL, DIMENSION( ims:ime , lm, 3 ) ::   qrs

  INTEGER ::               i,j,k,kflip 
  REAL, DIMENSION( ims:ime , lm ) ::                             &
                  pigen2d, psmlt2d                               &
                 ,pgfrz2d, pgacw2d, psacw2d                      &
                 ,pimlt2d, pihmf2d, pihtf2d                      &
                 ,praut2d, psaut2d, pgaut2d, pracw2d, psevp2d    &
                 ,pgacr2d, psaci2d, pgmlt2d, pgevp2d             &
                 ,pgaci2d, pseml2d, pgeml2d                      &
                 ,praci2d_s,praci2d_g                            &
                 ,piacr2d_s,piacr2d_g,pracs2d,psacr2d_s          &
                 ,psacr2d_g,psdep2d_d,psdep2d_s                  &
                 ,pgdep2d_d,pgdep2d_s,pidep2d_d,pidep2d_s        &
                 ,pcondc0,pconde0,prevp2d_e                      &
                 ,prevp2d_c,vt2s2d,vt2r2d,vt2i2d,vt2g2d
!-------------------------------------------------------------------
! IN NEMS k index is from top to bottom, while from bot to top in wsm6
      DO j=jms,jme
         DO k=1,lm
            kflip = lm-k+1
         DO i=ims,ime
            t(i,k)=th(i,j,kflip)*pii(i,j,kflip)
            qci(i,k,1) = qc(i,j,kflip)
            qci(i,k,2) = qi(i,j,kflip)
            qrs(i,k,1) = qr(i,j,kflip)
            qrs(i,k,2) = qs(i,j,kflip)
            qrs(i,k,3) = qg(i,j,kflip)
            q2(i,k) = q(i,j,kflip)
            p2(i,k) = p(i,j,kflip)
            delz2(i,k) = delz(i,j,kflip)
            den2(i,k) = den(i,j,kflip)            
          prevp2d_e(i,k)=0.
          prevp2d_c(i,k)=0.
          pigen2d(i,k)=0.
          psmlt2d(i,k)=0.
          pgfrz2d(i,k)=0.
          pgacw2d(i,k)=0.
          psacw2d(i,k)=0.
          pimlt2d(i,k)=0.
          pihmf2d(i,k)=0.
          pihtf2d(i,k)=0.
          praut2d(i,k)=0.
          psaut2d(i,k)=0.
          pgaut2d(i,k)=0.
          pracw2d(i,k)=0.
          pgevp2d(i,k)=0.
          pgacr2d(i,k)=0.
          psaci2d(i,k)=0.
          pgmlt2d(i,k)=0.
          pgaci2d(i,k)=0.
          pgeml2d(i,k)=0.
          psevp2d(i,k)=0.
          pseml2d(i,k)=0.
          praci2d_s(i,k)=0.
          praci2d_g(i,k)=0.
          piacr2d_s(i,k)=0.
          piacr2d_g(i,k)=0.
          pracs2d(i,k)=0.
          psacr2d_s(i,k)=0.
          psacr2d_g(i,k)=0.
          psdep2d_d(i,k)=0.
          psdep2d_s(i,k)=0.
          pgdep2d_d(i,k)=0.
          pgdep2d_s(i,k)=0.
          pidep2d_d(i,k)=0.
          pidep2d_s(i,k)=0.
          pcondc0(i,k)=0.
          pconde0(i,k)=0.
          vt2s2d(i,k)=0.
          vt2r2d(i,k)=0.
          vt2i2d(i,k)=0.
          vt2g2d(i,k)=0.
         ENDDO
         ENDDO
            
        CALL wsm62D(t, q2, qci, qrs                                &
                    ,den2                                          &
                    ,p2, delz2                                     &
                    ,delt,g, cpd, cpv, rd, rv, t0c                 &
                    ,ep1, ep2, qmin                                &
                    ,XLS, XLV0, XLF0, den0, denr                   &
                    ,cliq,cice,psat                                &
                    ,j                                             &
                    ,rain(ims,j),rainncv(ims,j)                    &
                    ,sr(ims,j)                                     &
                    ,ims,ime, jms,jme, lm                          &
                    ,snow(ims,j),snowncv(ims,j)                    &
                    ,graupel(ims,j),graupelncv(ims,j)              &
                 ,pigen2d,   psmlt2d                               &
                 ,pgfrz2d,   pgacw2d,   psacw2d,   pimlt2d,    pihmf2d    &
                 ,pihtf2d,   praut2d,   psaut2d,   pgaut2d,    pracw2d    &
                 ,pgevp2d,   pgacr2d,   psaci2d,   pgmlt2d,    pgaci2d    &
                 ,pgeml2d,   psevp2d,   pseml2d,   praci2d_s              &
                 ,praci2d_g, piacr2d_s, piacr2d_g, pracs2d,    prevp2d_e  &
                 ,psacr2d_s, psacr2d_g, psdep2d_d, prevp2d_c,  pconde0    &
                 ,psdep2d_s, pgdep2d_d, pidep2d_d, pidep2d_s,  pcondc0    &
                 ,pgdep2d_s, vt2s2d,    vt2r2d,    vt2i2d,     vt2g2d     &
                                                                   )

         DO K=1,lm
            kflip = lm-k+1
         DO I=ims,ime
            th(i,j,k)=t(i,kflip)/pii(i,j,k)
            qc(i,j,k) = qci(i,kflip,1)
            qi(i,j,k) = qci(i,kflip,2)
            qr(i,j,k) = qrs(i,kflip,1)
            qs(i,j,k) = qrs(i,kflip,2)
            qg(i,j,k) = qrs(i,kflip,3)
             q(i,j,k) = q2(i,kflip)
!---convert 2D source/sink terms to one 4D array
!---d_ss is the number of source/sink terms.  When d_ss is 1
!---only 1 source/sink term is used
         IF(D_SS.EQ.1)THEN
           mprates(I,J,K,1) = 0. 
         ELSE
           mprates(I,J,K,1) = mprates(I,J,K,1) + prevp2d_e(i,kflip)
           mprates(I,J,K,2) = mprates(I,J,K,2) + prevp2d_c(i,kflip)
           mprates(I,J,K,3) = mprates(I,J,K,3) + pigen2d(i,kflip)
           mprates(I,J,K,4) = mprates(I,J,K,4) + psmlt2d(i,kflip)
           mprates(I,J,K,5) = mprates(I,J,K,5) + pgfrz2d(i,kflip)
           mprates(I,J,K,6) = mprates(I,J,K,6) + pgacw2d(i,kflip)
           mprates(I,J,K,7) = mprates(I,J,K,7) + psacw2d(i,kflip)
           mprates(I,J,K,8) = mprates(I,J,K,8) + pimlt2d(i,kflip)
           mprates(I,J,K,9) = mprates(I,J,K,9) + pihmf2d(i,kflip)
           mprates(I,J,K,10) = mprates(I,J,K,10) + pihtf2d(i,kflip)
           mprates(I,J,K,11) = mprates(I,J,K,11) + praut2d(i,kflip)
           mprates(I,J,K,12) = mprates(I,J,K,12) + psaut2d(i,kflip)
           mprates(I,J,K,13) = mprates(I,J,K,13) + pgaut2d(i,kflip)
           mprates(I,J,K,14) = mprates(I,J,K,14) + pracw2d(i,kflip)
           mprates(I,J,K,15) = mprates(I,J,K,15) + pgacr2d(i,kflip)
           mprates(I,J,K,16) = mprates(I,J,K,16) + psaci2d(i,kflip)
           mprates(I,J,K,17) = mprates(I,J,K,17) + pgmlt2d(i,kflip)
           mprates(I,J,K,18) = mprates(I,J,K,18) + pgaci2d(i,kflip)
           mprates(I,J,K,19) = mprates(I,J,K,19) + pseml2d(i,kflip)
           mprates(I,J,K,20) = mprates(I,J,K,20) + pgeml2d(i,kflip)
           mprates(I,J,K,21) = mprates(I,J,K,21) + psevp2d(i,kflip)
           mprates(I,J,K,22) = mprates(I,J,K,22) + pgevp2d(i,kflip)
           mprates(I,J,K,23) = mprates(I,J,K,23) + praci2d_s(i,kflip)
           mprates(I,J,K,24) = mprates(I,J,K,24) + praci2d_g(i,kflip)
           mprates(I,J,K,25) = mprates(I,J,K,25) + piacr2d_s(i,kflip)
           mprates(I,J,K,26) = mprates(I,J,K,26) + piacr2d_g(i,kflip)
           mprates(I,J,K,27) = mprates(I,J,K,27) + pracs2d(i,kflip)
           mprates(I,J,K,28) = mprates(I,J,K,28) + psacr2d_s(i,kflip)
           mprates(I,J,K,29) = mprates(I,J,K,29) + psacr2d_g(i,kflip)
           mprates(I,J,K,30) = mprates(I,J,K,30) + psdep2d_d(i,kflip)
           mprates(I,J,K,31) = mprates(I,J,K,31) + psdep2d_s(i,kflip)
           mprates(I,J,K,32) = mprates(I,J,K,32) + pgdep2d_d(i,kflip)
           mprates(I,J,K,33) = mprates(I,J,K,33) + pgdep2d_s(i,kflip)
           mprates(I,J,K,34) = mprates(I,J,K,34) + pidep2d_d(i,kflip)
           mprates(I,J,K,35) = mprates(I,J,K,35) + pidep2d_s(i,kflip)
           mprates(I,J,K,36) = mprates(I,J,K,36) + pcondc0(i,kflip)
           mprates(I,J,K,37) = mprates(I,J,K,37) + pconde0(i,kflip)
           mprates(I,J,K,38) = vt2s2d(i,kflip)
           mprates(I,J,K,39) = vt2r2d(i,kflip)
           mprates(I,J,K,40) = vt2i2d(i,kflip)
           mprates(I,J,K,41) = vt2g2d(i,kflip)
         ENDIF
         ENDDO
         ENDDO
      ENDDO
  END SUBROUTINE wsm6
!===================================================================
!
  SUBROUTINE wsm62D(t, q, qci, qrs, den, p, delz                  &
                   ,delt,g, cpd, cpv, rd, rv, t0c                 &
                   ,ep1, ep2, qmin                                &
                   ,XLS, XLV0, XLF0, den0, denr                   &
                   ,cliq,cice,psat                                &
                   ,lat                                           &
                   ,rain,rainncv                                  &
                   ,sr                                            &
                   ,ims,ime, jms,jme, lm                          &
                   ,snow,snowncv                                  &
                   ,graupel,graupelncv                            &
                 ,pigen2d,   psmlt2d                              &
                 ,pgfrz2d,   pgacw2d,   psacw2d,   pimlt2d,    pihmf2d    &
                 ,pihtf2d,   praut2d,   psaut2d,   pgaut2d,    pracw2d    &
                 ,pgevp2d,   pgacr2d,   psaci2d,   pgmlt2d,    pgaci2d    &
                 ,pgeml2d,   psevp2d,   pseml2d,   praci2d_s              &
                 ,praci2d_g, piacr2d_s, piacr2d_g, pracs2d,    prevp2d_e  &
                 ,psacr2d_s, psacr2d_g, psdep2d_d, prevp2d_c,  pconde0    &
                 ,psdep2d_s, pgdep2d_d, pidep2d_d, pidep2d_s,  pcondc0    &
                 ,pgdep2d_s, vt2s2d,    vt2r2d,    vt2i2d,     vt2g2d     &
                                                                  )
!-------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------
  INTEGER,      INTENT(IN   )    ::   ims,ime, jms,jme, lm ,lat
  REAL, DIMENSION(ims:ime,lm),   INTENT(INOUT) :: t
  REAL, DIMENSION(ims:ime,lm,2), INTENT(INOUT) :: qci
  REAL, DIMENSION(ims:ime,lm,3), INTENT(INOUT) :: qrs
  REAL, DIMENSION(ims:ime,lm),   INTENT(INOUT) :: q
  REAL, DIMENSION(ims:ime,lm),   INTENT(IN   ) ::                 &
                                                             den, &
                                                               p, &
                                                            delz
  REAL, INTENT(IN   ) ::                                    delt, &
                                                               g, &
                                                             cpd, &
                                                             cpv, &
                                                             t0c, &
                                                            den0, &
                                                              rd, &
                                                              rv, &
                                                             ep1, &
                                                             ep2, &
                                                            qmin, &
                                                             XLS, &
                                                            XLV0, &
                                                            XLF0, &
                                                            cliq, &
                                                            cice, &
                                                            psat, &
                                                            denr
  REAL, DIMENSION( ims:ime ), INTENT(INOUT) ::              rain, &
                                                         rainncv, &
                                                              sr
  REAL, DIMENSION( ims:ime ),     OPTIONAL,                       &
        INTENT(INOUT) ::                                    snow, &
                                                         snowncv

  REAL, DIMENSION( ims:ime ),     OPTIONAL,                       &
        INTENT(INOUT) ::                                 graupel, &
                                                      graupelncv
! LOCAL VAR
  REAL, DIMENSION( ims:ime , lm , 3) ::                           &
        rh, qs, rslope, rslope2, rslope3, rslopeb,                &
        falk, fall, work1
  REAL, DIMENSION( ims:ime , lm ) ::                              &
              falkc, work1c, work2c, fallc
  REAL, DIMENSION( ims:ime , lm) ::                               &
        prevp, psdep, pgdep, praut, psaut, pgaut,                 &
        pracw, psacw, pgacw, pgacr, pgacs, psaci, pgmlt, praci,   &
        piacr, pracs, psacr, pgaci, pseml, pgeml      
  REAL, DIMENSION( ims:ime , lm ) ::                              &
        pigen, pidep, pcond, xl, cpm, work2, psmlt, psevp, denfac,&
        xni, pgevp,n0sfac
  REAL, DIMENSION( ims:ime , lm ) ::                              &
                  pigen2d                                         &
                 ,psmlt2d, pgfrz2d, pgacw2d, psacw2d              &
                 ,pimlt2d, pihmf2d, pihtf2d                       &
                 ,praut2d, psaut2d, pgaut2d, pracw2d, psevp2d     &
                 ,pgacr2d, psaci2d, pgmlt2d, pgevp2d              &
                 ,pgaci2d, pseml2d, pgeml2d                       &
                 ,pihtf,pgfrz,pimlt,pihmf                         &
                 ,praci2d_s                                       &
                 ,praci2d_g,piacr2d_s,piacr2d_g,fdelta2,fdelta3   &
                 ,pracs2d,  psacr2d_s,psacr2d_g                   &
                 ,psdep2d_d,psdep2d_s                             &
                 ,pgdep2d_d,pgdep2d_s,pidep2d_d,pidep2d_s         &
                 ,pcondc0,pconde0,prevp2d_e                       &
                 ,prevp2d_c,vt2s2d,vt2r2d,vt2i2d,vt2g2d
! variables for optimization
  REAL, DIMENSION( ims:ime )           :: tvec1
  REAL :: temp
  INTEGER, DIMENSION( ims:ime ) :: mstep, numdt
  LOGICAL, DIMENSION( ims:ime ) :: flgcld
  REAL  ::  pi,                                                   &
            cpmcal, xlcal, lamdar, lamdas, lamdag, diffus,        &
            viscos, xka, venfac, conden, diffac,                  &
            x, y, z, a, b, c, d, e,                               &
            qdt, holdrr, holdrs, holdrg, supcol, pvt,             &
            coeres, supsat, dtcld, xmi, eacrs, satdt,             &
            qimax, diameter, xni0, roqi0,                         &
            fallsum, fallsum_qsi, fallsum_qg,                     &
             vt2i,vt2r,vt2s,vt2g,acrfac,egs,egi,     &
            xlwork2, factor, source, value,              &
            xlf, pfrzdtc, pfrzdtr, supice, alpha2, delta2, delta3  
  REAL  :: holdc, holdci
  INTEGER :: i, j, k, mstepmax,                                   &
            iprt, latd, lond, loop, loops, ifsat, n
! Temporaries used for inlining fpvs function
  REAL  :: dldti, xb, xai, tr, xbi, xa, hvap, cvap, hsub, dldt, ttp
!
  real :: ppp1, c88   ! for test
  REAL, DIMENSION(ims:ime) :: rain_dt,snow_dt,graupel_dt     ! rain, snow and graupel at a big time step (i.e.
                                                             !, loop*dtcld)
!=================================================================
!   compute internal functions
!
      cpmcal(x) = cpd*(1.-max(x,qmin))+max(x,qmin)*cpv
      xlcal(x) = xlv0-xlv1*(x-t0c)
!----------------------------------------------------------------
!     size distributions: (x=mixing ratio, y=air density):
!     valid for mixing ratio > 1.e-9 kg/kg.
!
! Optimizatin : A**B => exp(log(A)*(B)) 
      lamdar(x,y)=   sqrt(sqrt(pidn0r/(x*y)))      ! (pidn0r/(x*y))**.25
      lamdas(x,y,z)= sqrt(sqrt(pidn0s*z/(x*y)))    ! (pidn0s*z/(x*y))**.25
      lamdag(x,y)=   sqrt(sqrt(pidn0g/(x*y)))      ! (pidn0g/(x*y))**.25
!
!----------------------------------------------------------------
!     diffus: diffusion coefficient of the water vapor
!     viscos: kinematic viscosity(m2s-1)
!
      diffus(x,y) = 8.794e-5 * exp(log(x)*(1.81)) / y        ! 8.794e-5*x**1.81/y
      viscos(x,y) = 1.496e-6 * (x*sqrt(x)) /(x+120.)/y  ! 1.496e-6*x**1.5/(x+120.)/y
      xka(x,y) = 1.414e3*viscos(x,y)*y
      diffac(a,b,c,d,e) = d*a*a/(xka(c,d)*rv*c*c)+1./(e*diffus(c,b))
      venfac(a,b,c) = exp(log((viscos(b,c)/diffus(b,a)))*((.3333333)))    &
                     /sqrt(viscos(b,c))*sqrt(sqrt(den0/c))
      conden(a,b,c,d,e) = (max(b,qmin)-c)/(1.+d*d/(rv*e)*c/(a*a))
!
      pi = 4. * atan(1.)
      c88 =1.0 
!
!
!----------------------------------------------------------------
!     paddint 0 for negative values generated by dynamics
!
      do k = 1, lm
        do i = ims, ime
          qci(i,k,1) = max(qci(i,k,1),0.0)
          qrs(i,k,1) = max(qrs(i,k,1),0.0)
          qci(i,k,2) = max(qci(i,k,2),0.0)
          qrs(i,k,2) = max(qrs(i,k,2),0.0)
          qrs(i,k,3) = max(qrs(i,k,3),0.0)
        enddo
      enddo
!
!----------------------------------------------------------------
!     latent heat for phase changes and heat capacity. neglect the
!     changes during microphysical process calculation
!     emanuel(1994)
!
      do k = 1, lm
        do i = ims, ime
          cpm(i,k) = cpmcal(q(i,k))
          xl(i,k) = xlcal(t(i,k))
        enddo
      enddo
!
!----------------------------------------------------------------
!     compute the minor time steps.
!
      loops = max(nint(delt/dtcldcr),1)
      dtcld = delt/loops
      if(delt.le.dtcldcr) dtcld = delt

      do i = ims, ime
         rain_dt(i) = 0.0
         snow_dt(i) = 0.0
         graupel_dt(i) = 0.0
      enddo
!
      do loop = 1,loops
!
!----------------------------------------------------------------
!     initialize the large scale variables
!
      do i = ims, ime
        mstep(i) = 1
        flgcld(i) = .true.
      enddo
!
      do k = 1, lm
        CALL VREC( tvec1(ims), den(ims,k), ime-ims+1)
        do i = ims, ime
          tvec1(i) = tvec1(i)*den0
        enddo
        CALL VSQRT( denfac(ims,k), tvec1(ims), ime-ims+1)
      enddo
!
! Inline expansion for fpvs
!         qs(i,k,1) = fpvs(t(i,k),0,rd,rv,cpv,cliq,cice,xlv0,xls,psat,t0c)
!         qs(i,k,2) = fpvs(t(i,k),1,rd,rv,cpv,cliq,cice,xlv0,xls,psat,t0c)
      hsub = xls
      hvap = xlv0
      cvap = cpv
      ttp=t0c+0.01
      dldt=cvap-cliq
      xa=-dldt/rv
      xb=xa+hvap/(rv*ttp)
      dldti=cvap-cice
      xai=-dldti/rv
      xbi=xai+hsub/(rv*ttp)
      do k = 1, lm
        do i = ims, ime
          tr=ttp/t(i,k)
          qs(i,k,1)=psat*exp(log(tr)*(xa))*exp(xb*(1.-tr))
          qs(i,k,1) = ep2 * qs(i,k,1) / (p(i,k) - qs(i,k,1))
          qs(i,k,1) = max(qs(i,k,1),qmin)
          rh(i,k,1) = max(q(i,k) / qs(i,k,1),qmin)
          tr=ttp/t(i,k)
          if(t(i,k).lt.ttp) then
            qs(i,k,2)=psat*exp(log(tr)*(xai))*exp(xbi*(1.-tr))
          else
            qs(i,k,2)=psat*exp(log(tr)*(xa))*exp(xb*(1.-tr))
          endif
          qs(i,k,2) = ep2 * qs(i,k,2) / (p(i,k) - qs(i,k,2))
          qs(i,k,2) = max(qs(i,k,2),qmin)
          rh(i,k,2) = max(q(i,k) / qs(i,k,2),qmin)
        enddo
      enddo
!
!----------------------------------------------------------------
!     initialize the variables for microphysical physics
!
!
      do k = 1, lm
        do i = ims, ime
          prevp(i,k) = 0.
          psdep(i,k) = 0.
          pgdep(i,k) = 0.
          praut(i,k) = 0.
          psaut(i,k) = 0.
          pgaut(i,k) = 0.
          pracw(i,k) = 0.
          praci(i,k) = 0.
          piacr(i,k) = 0.
          psaci(i,k) = 0.
          psacw(i,k) = 0.
          pracs(i,k) = 0.
          psacr(i,k) = 0.
          pgacw(i,k) = 0.
          pgaci(i,k) = 0.
          pgacr(i,k) = 0.
          pgacs(i,k) = 0.
          pigen(i,k) = 0.
          pidep(i,k) = 0.
          pcond(i,k) = 0.
          pimlt(i,k) = 0.
          pihmf(i,k) = 0.
          pihtf(i,k) = 0.
          pgfrz(i,k) = 0.
          psmlt(i,k) = 0.
          pgmlt(i,k) = 0.
          pseml(i,k) = 0.
          pgeml(i,k) = 0.
          psevp(i,k) = 0.
          pgevp(i,k) = 0.
          falk(i,k,1) = 0.
          falk(i,k,2) = 0.
          falk(i,k,3) = 0.
          fall(i,k,1) = 0.
          fall(i,k,2) = 0.
          fall(i,k,3) = 0.
          fallc(i,k) = 0.
          falkc(i,k) = 0.
          fdelta3(i,k)=0.
          fdelta2(i,k)=0.
          xni(i,k) = 1.e3
        enddo
      enddo
!
!----------------------------------------------------------------
!     compute the fallout term:
!     first, vertical terminal velosity for minor loops
!
      do k = 1, lm
        do i = ims, ime
          supcol = t0c-t(i,k)
!---------------------------------------------------------------
! n0s: Intercept parameter for snow [m-4] [HDC 6]
!---------------------------------------------------------------
          n0sfac(i,k) = max(min(exp(alpha*supcol),n0smax/n0s),1.)
          if(qrs(i,k,1).le.qcrmin)then
            rslope(i,k,1) = rslopermax
            rslopeb(i,k,1) = rsloperbmax
            rslope2(i,k,1) = rsloper2max
            rslope3(i,k,1) = rsloper3max
          else
            rslope(i,k,1) = 1./lamdar(qrs(i,k,1),den(i,k))
            rslopeb(i,k,1) = rslope(i,k,1)**bvtr
            rslope2(i,k,1) = rslope(i,k,1)*rslope(i,k,1)
            rslope3(i,k,1) = rslope2(i,k,1)*rslope(i,k,1)
          endif
          if(qrs(i,k,2).le.qcrmin)then
            rslope(i,k,2) = rslopesmax
            rslopeb(i,k,2) = rslopesbmax
            rslope2(i,k,2) = rslopes2max
            rslope3(i,k,2) = rslopes3max
          else
            rslope(i,k,2) = 1./lamdas(qrs(i,k,2),den(i,k),n0sfac(i,k))
            rslopeb(i,k,2) = rslope(i,k,2)**bvts
            rslope2(i,k,2) = rslope(i,k,2)*rslope(i,k,2)
            rslope3(i,k,2) = rslope2(i,k,2)*rslope(i,k,2)
          endif
          if(qrs(i,k,3).le.qcrmin)then
            rslope(i,k,3) = rslopegmax
            rslopeb(i,k,3) = rslopegbmax
            rslope2(i,k,3) = rslopeg2max
            rslope3(i,k,3) = rslopeg3max
          else
            rslope(i,k,3) = 1./lamdag(qrs(i,k,3),den(i,k))
            rslopeb(i,k,3) = rslope(i,k,3)**bvtg
            rslope2(i,k,3) = rslope(i,k,3)*rslope(i,k,3)
            rslope3(i,k,3) = rslope2(i,k,3)*rslope(i,k,3)
          endif
!-------------------------------------------------------------
! Ni: ice crystal number concentraiton   [HDC 5c]
!-------------------------------------------------------------
!         xni(i,k) = min(max(5.38e7*(den(i,k)                           &
!                   *max(qci(i,k,2),qmin))**0.75,1.e3),1.e6)
          temp = (den(i,k)*max(qci(i,k,2),qmin))
          temp = sqrt(sqrt(temp*temp*temp))
          xni(i,k) = min(max(5.38e7*temp,1.e3),1.e6)
        enddo
      enddo
!
      mstepmax = 1
      numdt = 1
      do k = lm, 1, -1
        do i = ims, ime
          work1(i,k,1) = pvtr*rslopeb(i,k,1)*denfac(i,k)/delz(i,k)
          work1(i,k,2) = pvts*rslopeb(i,k,2)*denfac(i,k)/delz(i,k)
          work1(i,k,3) = pvtg*rslopeb(i,k,3)*denfac(i,k)/delz(i,k)
          numdt(i) = max(nint(max(work1(i,k,1),work1(i,k,2),work1(i,k,3)) &
                    *dtcld+.5),1)
          if(numdt(i).ge.mstep(i)) mstep(i) = numdt(i)
        enddo
      enddo
      do i = ims, ime
        if(mstepmax.le.mstep(i)) mstepmax = mstep(i)
      enddo
!
      do n = 1, mstepmax
        k = lm
        do i = ims, ime
          if(n.le.mstep(i)) then
              falk(i,k,1) = den(i,k)*qrs(i,k,1)*work1(i,k,1)/mstep(i)
              falk(i,k,2) = den(i,k)*qrs(i,k,2)*work1(i,k,2)/mstep(i)
              falk(i,k,3) = den(i,k)*qrs(i,k,3)*work1(i,k,3)/mstep(i)
              fall(i,k,1) = fall(i,k,1)+falk(i,k,1)
              fall(i,k,2) = fall(i,k,2)+falk(i,k,2)
              fall(i,k,3) = fall(i,k,3)+falk(i,k,3)
              qrs(i,k,1) = max(qrs(i,k,1)-falk(i,k,1)*dtcld/den(i,k),0.)
              qrs(i,k,2) = max(qrs(i,k,2)-falk(i,k,2)*dtcld/den(i,k),0.)
              qrs(i,k,3) = max(qrs(i,k,3)-falk(i,k,3)*dtcld/den(i,k),0.)
            endif
          enddo
        do k = lm-1, 1, -1
          do i = ims, ime
            if(n.le.mstep(i)) then
              falk(i,k,1) = den(i,k)*qrs(i,k,1)*work1(i,k,1)/mstep(i)
              falk(i,k,2) = den(i,k)*qrs(i,k,2)*work1(i,k,2)/mstep(i)
              falk(i,k,3) = den(i,k)*qrs(i,k,3)*work1(i,k,3)/mstep(i)
              fall(i,k,1) = fall(i,k,1)+falk(i,k,1)
              fall(i,k,2) = fall(i,k,2)+falk(i,k,2)
              fall(i,k,3) = fall(i,k,3)+falk(i,k,3)
              qrs(i,k,1) = max(qrs(i,k,1)-(falk(i,k,1)-falk(i,k+1,1)    &
                          *delz(i,k+1)/delz(i,k))*dtcld/den(i,k),0.)
              qrs(i,k,2) = max(qrs(i,k,2)-(falk(i,k,2)-falk(i,k+1,2)    &
                          *delz(i,k+1)/delz(i,k))*dtcld/den(i,k),0.)
              qrs(i,k,3) = max(qrs(i,k,3)-(falk(i,k,3)-falk(i,k+1,3)    &
                          *delz(i,k+1)/delz(i,k))*dtcld/den(i,k),0.)
            endif
          enddo
        enddo
        do k = lm, 1, -1
          do i = ims, ime
            if(n.le.mstep(i).and.t(i,k).gt.t0c) then
!---------------------------------------------------------------
! psmlt: melting of snow [HL A33] [RH83 A25]
!       (T>T0: S->R)
!---------------------------------------------------------------
              xlf = xlf0
              work2(i,k) = venfac(p(i,k),t(i,k),den(i,k))
              if(qrs(i,k,2).gt.0.) then
                coeres = rslope2(i,k,2)*sqrt(rslope(i,k,2)*rslopeb(i,k,2))
                psmlt(i,k) = xka(t(i,k),den(i,k))/xlf*(t0c-t(i,k))*pi/2. &
                           *n0sfac(i,k)*(precs1*rslope2(i,k,2)          &
                           +precs2*work2(i,k)*coeres)
                psmlt(i,k) = min(max(psmlt(i,k)*dtcld/mstep(i),           &
                            -qrs(i,k,2)/mstep(i)),0.)
                qrs(i,k,2) = qrs(i,k,2) + psmlt(i,k)
                qrs(i,k,1) = qrs(i,k,1) - psmlt(i,k)
                t(i,k) = t(i,k) + xlf/cpm(i,k)*psmlt(i,k)*c88
              endif
!---------------------------------------------------------------
! pgmlt: melting of graupel [HL A23]  [LFO 47]
!       (T>T0: G->R)
!---------------------------------------------------------------
              if(qrs(i,k,3).gt.0.) then
                coeres = rslope2(i,k,3)*sqrt(rslope(i,k,3)*rslopeb(i,k,3))
                pgmlt(i,k) = xka(t(i,k),den(i,k))/xlf                    &
                           *(t0c-t(i,k))*(precg1*rslope2(i,k,3)         &
                           +precg2*work2(i,k)*coeres)
                pgmlt(i,k) = min(max(pgmlt(i,k)*dtcld/mstep(i),           &
                            -qrs(i,k,3)/mstep(i)),0.)
                qrs(i,k,3) = qrs(i,k,3) + pgmlt(i,k)
                qrs(i,k,1) = qrs(i,k,1) - pgmlt(i,k)
                t(i,k) = t(i,k) + xlf/cpm(i,k)*pgmlt(i,k)*c88
              endif
            endif
          enddo
        enddo
      enddo
!---------------------------------------------------------------
! Vice [ms-1] : fallout of ice crystal [HDC 5a]
!---------------------------------------------------------------
      mstepmax = 1
      mstep = 1
      numdt = 1
      do k = lm, 1, -1
        do i = ims, ime
          if(qci(i,k,2).le.0.) then
            work2c(i,k) = 0.
          else
            xmi = den(i,k)*qci(i,k,2)/xni(i,k)
!           diameter  = min(dicon * sqrt(xmi),dimax)
            diameter  = max(min(dicon * sqrt(xmi),dimax), 1.e-25)
            work1c(i,k) = 1.49e4*diameter**1.31
            work2c(i,k) = work1c(i,k)/delz(i,k)
          endif
          numdt(i) = max(nint(work2c(i,k)*dtcld+.5),1)
          if(numdt(i).ge.mstep(i)) mstep(i) = numdt(i)
        enddo
      enddo
      do i = ims, ime
        if(mstepmax.le.mstep(i)) mstepmax = mstep(i)
      enddo
!
      do n = 1, mstepmax
        k = lm
        do i = ims, ime
          if(n.le.mstep(i)) then
            falkc(i,k) = den(i,k)*qci(i,k,2)*work2c(i,k)/mstep(i)
            holdc = falkc(i,k)
            fallc(i,k) = fallc(i,k)+falkc(i,k)
            holdci = qci(i,k,2)
            qci(i,k,2) = max(qci(i,k,2)-falkc(i,k)*dtcld/den(i,k),0.)
          endif
        enddo
        do k = lm-1, 1, -1
          do i = ims, ime
            if(n.le.mstep(i)) then
              falkc(i,k) = den(i,k)*qci(i,k,2)*work2c(i,k)/mstep(i)
              holdc = falkc(i,k)
              fallc(i,k) = fallc(i,k)+falkc(i,k)
              holdci = qci(i,k,2)
              qci(i,k,2) = max(qci(i,k,2)-(falkc(i,k)-falkc(i,k+1)      &
                          *delz(i,k+1)/delz(i,k))*dtcld/den(i,k),0.)
            endif
          enddo
        enddo
      enddo
!
!----------------------------------------------------------------
!      rain (unit is mm/sec;kgm-2s-1: /1000*delt ===> m)==> mm for wrf
!
      do i = ims, ime
        fallsum = fall(i,1,1)+fall(i,1,2)+fall(i,1,3)+fallc(i,1)
        fallsum_qsi = fall(i,1,2)+fallc(i,1)
        fallsum_qg = fall(i,1,3)
        rainncv(i) = 0.
        if(fallsum.gt.0.) then
          rainncv(i) = fallsum*delz(i,1)/denr*dtcld*1000.
          rain(i) = fallsum*delz(i,1)/denr*dtcld*1000. + rain(i)
          rain_dt(i) = rain_dt(i) + rainncv(i)                         ! W. Wang 2/26/10
        endif
        IF ( PRESENT (snowncv) .AND. PRESENT (snow)) THEN
        snowncv(i) = 0.
        if(fallsum_qsi.gt.0.) then
          snowncv(i) = fallsum_qsi*delz(i,1)/denr*dtcld*1000.
          snow(i) = fallsum_qsi*delz(i,1)/denr*dtcld*1000. + snow(i)
          snow_dt(i) = snow_dt(i) + snowncv(i)                         ! W. wang 2/26/10
        endif
        ENDIF
        IF ( PRESENT (graupelncv) .AND. PRESENT (graupel)) THEN
        graupelncv(i) = 0.
        if(fallsum_qg.gt.0.) then
          graupelncv(i) = fallsum_qg*delz(i,1)/denr*dtcld*1000.
          graupel(i) = fallsum_qg*delz(i,1)/denr*dtcld*1000. + graupel(i)
          graupel_dt(i) = graupel_dt(i) + graupelncv(i)                !W. wang 2/26/10
        endif
        ENDIF
        sr(i) = 0.
        if(fallsum.gt.0.)sr(i)=(fallsum_qsi*delz(i,1)/denr*dtcld*1000. + &
                                fallsum_qg*delz(i,1)/denr*dtcld*1000.)/(rainncv(i)+1.e-12)
      enddo
!
!---------------------------------------------------------------
! pimlt: instantaneous melting of cloud ice [HL A47] [RH83 A28]
!       (T>T0: I->C)
!---------------------------------------------------------------
      do k = 1, lm
        do i = ims, ime
          supcol = t0c-t(i,k)
          xlf = xls-xl(i,k)
          if(supcol.lt.0.) xlf = xlf0
          if(supcol.lt.0.and.qci(i,k,2).gt.0.) then
!
!aligo
!  pimlt was added as a diagnostic array not in any other version of
!  the WSM6 scheme.
!  pimlt,pihmf,pihtf,and pgfrz are the added diagnostics.  Declared
!  in wsm62D and initialized there as well 
! 
!not instantaneous
            pimlt(i,k)=qci(i,k,2)
            pimlt(i,k) = max(qmin,pimlt(i,k))
!aligo
            qci(i,k,1) = qci(i,k,1) + qci(i,k,2)
            t(i,k) = t(i,k) - xlf/cpm(i,k)*qci(i,k,2)*c88
            qci(i,k,2) = 0.
          endif
!---------------------------------------------------------------
! pihmf: homogeneous freezing of cloud water below -40c [HL A45]
!        (T<-40C: C->I)
!---------------------------------------------------------------
          if(supcol.gt.40..and.qci(i,k,1).gt.0.) then
!aligo
            pihmf(i,k) = qci(i,k,1)
            pihmf(i,k) = max(qmin,pihmf(i,k))
!aligo
            qci(i,k,2) = qci(i,k,2) + qci(i,k,1)
            t(i,k) = t(i,k) + xlf/cpm(i,k)*qci(i,k,1)*c88
            qci(i,k,1) = 0.
          endif
!---------------------------------------------------------------
! pihtf: heterogeneous freezing of cloud water [HL A44]
!        (T0>T>-40C: C->I)
!---------------------------------------------------------------
          if(supcol.gt.0..and.qci(i,k,1).gt.qmin) then
!           pfrzdtc = min(pfrz1*(exp(pfrz2*supcol)-1.)                  &
!              *den(i,k)/denr/xncr*qci(i,k,1)**2*dtcld,qci(i,k,1))
            pfrzdtc = min(pfrz1*(exp(pfrz2*supcol)-1.)                  &
            *den(i,k)/denr/xncr*qci(i,k,1)*qci(i,k,1)*dtcld,qci(i,k,1))
!aligo
            pihtf(i,k) = pfrzdtc
            pihtf(i,k) = max(qmin,pihtf(i,k))
!aligo
            qci(i,k,2) = qci(i,k,2) + pfrzdtc
            t(i,k) = t(i,k) + xlf/cpm(i,k)*pfrzdtc*c88
            qci(i,k,1) = qci(i,k,1)-pfrzdtc
          endif
!---------------------------------------------------------------
! pgfrz: freezing of rain water [HL A20] [LFO 45]
!        (T<T0, R->G)
!---------------------------------------------------------------
          if(supcol.gt.0..and.qrs(i,k,1).gt.0.) then
!           pfrzdtr = min(20.*pi**2*pfrz1*n0r*denr/den(i,k)             &
!                 *(exp(pfrz2*supcol)-1.)*rslope3(i,k,1)**2             &
!                 *rslope(i,k,1)*dtcld,qrs(i,k,1))
            temp = rslope3(i,k,1)
   !i         temp = temp*temp*rslope(i,k,1)
         !   pfrzdtr = min(20.*(pi*pi)*pfrz1*n0r*denr/den(i,k)           &
         !         *(exp(pfrz2*supcol)-1.)*temp*dtcld,                   &
         !         qrs(i,k,1))
                  pfrzdtr = pfrz2*supcol
                  pfrzdtr = exp(pfrzdtr)
                  pfrzdtr = (pfrzdtr -1.0)*temp
                   
                  ppp1= temp*rslope(i,k,1)*20.*(pi*pi)*pfrz1*n0r*denr/den(i,k)
                  pfrzdtr = pfrzdtr*ppp1
              pfrzdtr = pfrzdtr*dtcld
               pfrzdtr = min( pfrzdtr, qrs(i,k,1) ) 
!aligo
            pgfrz(i,k) = pfrzdtr
            pgfrz(i,k) = max(qmin,pgfrz(i,k))
!aligo
            qrs(i,k,3) = qrs(i,k,3) + pfrzdtr
            t(i,k) = t(i,k) + xlf/cpm(i,k)*pfrzdtr*c88
            qrs(i,k,1) = qrs(i,k,1)-pfrzdtr
          endif
        enddo
      enddo
!
!
!----------------------------------------------------------------
!     rsloper: reverse of the slope parameter of the rain(m)
!     xka:    thermal conductivity of air(jm-1s-1k-1)
!     work1:  the thermodynamic term in the denominator associated with
!             heat conduction and vapor diffusion
!             (ry88, y93, h85)
!     work2: parameter associated with the ventilation effects(y93)
!
      do k = 1, lm
        do i = ims, ime
          if(qrs(i,k,1).le.qcrmin)then
            rslope(i,k,1) = rslopermax
            rslopeb(i,k,1) = rsloperbmax
            rslope2(i,k,1) = rsloper2max
            rslope3(i,k,1) = rsloper3max
          else
            rslope(i,k,1) = 1./lamdar(qrs(i,k,1),den(i,k))
            rslopeb(i,k,1) = rslope(i,k,1)**bvtr
            rslope2(i,k,1) = rslope(i,k,1)*rslope(i,k,1)
            rslope3(i,k,1) = rslope2(i,k,1)*rslope(i,k,1)
          endif
          if(qrs(i,k,2).le.qcrmin)then
            rslope(i,k,2) = rslopesmax
            rslopeb(i,k,2) = rslopesbmax
            rslope2(i,k,2) = rslopes2max
            rslope3(i,k,2) = rslopes3max
          else
            rslope(i,k,2) = 1./lamdas(qrs(i,k,2),den(i,k),n0sfac(i,k))
            rslopeb(i,k,2) = rslope(i,k,2)**bvts
            rslope2(i,k,2) = rslope(i,k,2)*rslope(i,k,2)
            rslope3(i,k,2) = rslope2(i,k,2)*rslope(i,k,2)
          endif
          if(qrs(i,k,3).le.qcrmin)then
            rslope(i,k,3) = rslopegmax
            rslopeb(i,k,3) = rslopegbmax
            rslope2(i,k,3) = rslopeg2max
            rslope3(i,k,3) = rslopeg3max
          else
            rslope(i,k,3) = 1./lamdag(qrs(i,k,3),den(i,k))
            rslopeb(i,k,3) = rslope(i,k,3)**bvtg
            rslope2(i,k,3) = rslope(i,k,3)*rslope(i,k,3)
            rslope3(i,k,3) = rslope2(i,k,3)*rslope(i,k,3)
          endif
        enddo
      enddo
!
      do k = 1, lm
        do i = ims, ime
          work1(i,k,1) = diffac(xl(i,k),p(i,k),t(i,k),den(i,k),qs(i,k,1))
          work1(i,k,2) = diffac(xls,p(i,k),t(i,k),den(i,k),qs(i,k,2))
          work2(i,k) = venfac(p(i,k),t(i,k),den(i,k))
        enddo
      enddo
!
!===============================================================
!
! warm rain processes
!
! - follows the processes in RH83 and LFO except for autoconcersion
!
!===============================================================
!
      do k = 1, lm
        do i = ims, ime
          supsat = max(q(i,k),qmin)-qs(i,k,1)
          satdt = supsat/dtcld
!---------------------------------------------------------------
! praut: auto conversion rate from cloud to rain [HDC 16]
!        (C->R)
!---------------------------------------------------------------
          if(qci(i,k,1).gt.qc0) then
            praut(i,k) = qck1*qci(i,k,1)**(7./3.)
            praut(i,k) = min(praut(i,k),qci(i,k,1)/dtcld)
          endif
!---------------------------------------------------------------
! pracw: accretion of cloud water by rain [HL A40] [LFO 51]
!        (C->R)
!---------------------------------------------------------------
          if(qrs(i,k,1).gt.qcrmin.and.qci(i,k,1).gt.qmin) then
            pracw(i,k) = min(pacrr*rslope3(i,k,1)*rslopeb(i,k,1)        &
                        *qci(i,k,1)*denfac(i,k),qci(i,k,1)/dtcld)
          endif
!---------------------------------------------------------------
! prevp: evaporation/condensation rate of rain [HDC 14]
!        (V->R or R->V)
!---------------------------------------------------------------
          if(qrs(i,k,1).gt.0.) then
            coeres = rslope2(i,k,1)*sqrt(rslope(i,k,1)*rslopeb(i,k,1))
            prevp(i,k) = (rh(i,k,1)-1.)*(precr1*rslope2(i,k,1)         &
                         +precr2*work2(i,k)*coeres)/work1(i,k,1)
            if(prevp(i,k).lt.0.) then
              prevp(i,k) = max(prevp(i,k),-qrs(i,k,1)/dtcld)
              prevp(i,k) = max(prevp(i,k),satdt/2)
            else
              prevp(i,k) = min(prevp(i,k),satdt/2)
            endif
          endif
        enddo
      enddo
!
!===============================================================
!
! cold rain processes
!
! - follows the revised ice microphysics processes in HDC
! - the processes same as in RH83 and RH84  and LFO behave
!   following ice crystal hapits defined in HDC, inclduing
!   intercept parameter for snow (n0s), ice crystal number
!   concentration (ni), ice nuclei number concentration
!   (n0i), ice diameter (d)
!
!===============================================================
!
      do k = 1, lm
        do i = ims, ime
          supcol = t0c-t(i,k)
          supsat = max(q(i,k),qmin)-qs(i,k,2)
          satdt = supsat/dtcld
          ifsat = 0
!-------------------------------------------------------------
! Ni: ice crystal number concentraiton   [HDC 5c]
!-------------------------------------------------------------
!         xni(i,k) = min(max(5.38e7*(den(i,k)                           &
!                      *max(qci(i,k,2),qmin))**0.75,1.e3),1.e6)
          temp = (den(i,k)*max(qci(i,k,2),qmin))
          temp = sqrt(sqrt(temp*temp*temp))
          xni(i,k) = min(max(5.38e7*temp,1.e3),1.e6)
          eacrs = exp(0.07*(-supcol))
!
          xmi = den(i,k)*qci(i,k,2)/xni(i,k)
          diameter  = min(dicon * sqrt(xmi),dimax)
          vt2i = 1.49e4*diameter**1.31
          vt2r=pvtr*rslopeb(i,k,1)*denfac(i,k)
          vt2s=pvts*rslopeb(i,k,2)*denfac(i,k)
          vt2g=pvtg*rslopeb(i,k,3)*denfac(i,k)
!aligo
          vt2i2d(i,k)=vt2i
          vt2r2d(i,k)=vt2r
          vt2s2d(i,k)=vt2s
          vt2g2d(i,k)=vt2g
!aligo
          if(supcol.gt.0.and.qci(i,k,2).gt.qmin) then
            if(qrs(i,k,1).gt.qcrmin) then
!-------------------------------------------------------------
! praci: Accretion of cloud ice by rain [HL A15] [LFO 25]
!        (T<T0: I->R)
!-------------------------------------------------------------
              acrfac = 2.*rslope3(i,k,1)+2.*diameter*rslope2(i,k,1)     &
                      +diameter**2*rslope(i,k,1)
              praci(i,k) = pi*qci(i,k,2)*n0r*abs(vt2r-vt2i)*acrfac/4.
              praci(i,k) = min(praci(i,k),qci(i,k,2)/dtcld)
!-------------------------------------------------------------
! piacr: Accretion of rain by cloud ice [HL A19] [LFO 26]
!        (T<T0: R->S or R->G)
!-------------------------------------------------------------
              piacr(i,k) = pi**2*avtr*n0r*denr*xni(i,k)*denfac(i,k)     &
                          *g6pbr*rslope3(i,k,1)*rslope3(i,k,1)          &
                          *rslopeb(i,k,1)/24./den(i,k)
              piacr(i,k) = min(piacr(i,k),qrs(i,k,1)/dtcld)
            endif
!-------------------------------------------------------------
! psaci: Accretion of cloud ice by snow [HDC 10]
!        (T<T0: I->S)
!-------------------------------------------------------------
            if(qrs(i,k,2).gt.qcrmin) then
              acrfac = 2.*rslope3(i,k,2)+2.*diameter*rslope2(i,k,2)     &
                      +diameter**2*rslope(i,k,2)
              psaci(i,k) = pi*qci(i,k,2)*eacrs*n0s*n0sfac(i,k)          &
                          *abs(vt2s-vt2i)*acrfac/4.
              psaci(i,k) = min(psaci(i,k),qci(i,k,2)/dtcld)
            endif
!-------------------------------------------------------------
! pgaci: Accretion of cloud ice by graupel [HL A17] [LFO 41]
!        (T<T0: I->G)
!-------------------------------------------------------------
            if(qrs(i,k,3).gt.qcrmin) then
              egi = exp(0.07*(-supcol))
              acrfac = 2.*rslope3(i,k,3)+2.*diameter*rslope2(i,k,3)     &
                      +diameter**2*rslope(i,k,3)
              pgaci(i,k) = pi*egi*qci(i,k,2)*n0g*abs(vt2g-vt2i)*acrfac/4.
              pgaci(i,k) = min(pgaci(i,k),qci(i,k,2)/dtcld)
            endif
          endif
!-------------------------------------------------------------
! psacw: Accretion of cloud water by snow  [HL A7] [LFO 24]
!        (T<T0: C->G, and T>=T0: C->R)
!-------------------------------------------------------------
          if(qrs(i,k,2).gt.qcrmin.and.qci(i,k,1).gt.qmin) then
            psacw(i,k) = min(pacrc*n0sfac(i,k)*rslope3(i,k,2)           &
                        *rslopeb(i,k,2)*qci(i,k,1)*denfac(i,k)          &
                        ,qci(i,k,1)/dtcld)
          endif
!-------------------------------------------------------------
! pgacw: Accretion of cloud water by graupel [HL A6] [LFO 40]
!        (T<T0: C->G, and T>=T0: C->R)
!-------------------------------------------------------------
          if(qrs(i,k,3).gt.qcrmin.and.qci(i,k,1).gt.qmin) then
            pgacw(i,k) = min(pacrg*rslope3(i,k,3)*rslopeb(i,k,3)        &
                        *qci(i,k,1)*denfac(i,k),qci(i,k,1)/dtcld)
          endif
!-------------------------------------------------------------
! pracs: Accretion of snow by rain [HL A11] [LFO 27]
!         (T<T0: S->G)
!-------------------------------------------------------------
          if(qrs(i,k,2).gt.qcrmin.and.qrs(i,k,1).gt.qcrmin) then
            if(supcol.gt.0) then
              acrfac = 5.*rslope3(i,k,2)*rslope3(i,k,2)*rslope(i,k,1)   &
                      +2.*rslope3(i,k,2)*rslope2(i,k,2)*rslope2(i,k,1)  &
                      +.5*rslope2(i,k,2)*rslope2(i,k,2)*rslope3(i,k,1)
              pracs(i,k) = pi**2*n0r*n0s*n0sfac(i,k)*abs(vt2r-vt2s)     &
                          *(dens/den(i,k))*acrfac
              pracs(i,k) = min(pracs(i,k),qrs(i,k,2)/dtcld)
            endif
!-------------------------------------------------------------
! psacr: Accretion of rain by snow [HL A10] [LFO 28]
!         (T<T0:R->S or R->G) (T>=T0: enhance melting of snow)
!-------------------------------------------------------------
            acrfac = 5.*rslope3(i,k,1)*rslope3(i,k,1)*rslope(i,k,2)     &
                    +2.*rslope3(i,k,1)*rslope2(i,k,1)*rslope2(i,k,2)    &
                    +.5*rslope2(i,k,1)*rslope2(i,k,1)*rslope3(i,k,2)
            psacr(i,k) = pi**2*n0r*n0s*n0sfac(i,k)*abs(vt2s-vt2r)       &
                        *(denr/den(i,k))*acrfac
            psacr(i,k) = min(psacr(i,k),qrs(i,k,1)/dtcld)
          endif
!-------------------------------------------------------------
! pgacr: Accretion of rain by graupel [HL A12] [LFO 42]
!         (T<T0: R->G) (T>=T0: enhance melting of graupel)
!-------------------------------------------------------------
          if(qrs(i,k,3).gt.qcrmin.and.qrs(i,k,1).gt.qcrmin) then
            acrfac = 5.*rslope3(i,k,1)*rslope3(i,k,1)*rslope(i,k,3)     &
                    +2.*rslope3(i,k,1)*rslope2(i,k,1)*rslope2(i,k,3)    &
                    +.5*rslope2(i,k,1)*rslope2(i,k,1)*rslope3(i,k,3)
            pgacr(i,k) = pi**2*n0r*n0g*abs(vt2g-vt2r)*(denr/den(i,k))   &
                        *acrfac
            pgacr(i,k) = min(pgacr(i,k),qrs(i,k,1)/dtcld)
          endif
!
!-------------------------------------------------------------
! pgacs: Accretion of snow by graupel [HL A13] [LFO 29]
!        (S->G)
!-------------------------------------------------------------
          if(qrs(i,k,3).gt.qcrmin.and.qrs(i,k,2).gt.qcrmin) then
            acrfac = 5.*rslope3(i,k,2)*rslope3(i,k,2)*rslope(i,k,3)     &
                    +2.*rslope3(i,k,2)*rslope2(i,k,2)*rslope2(i,k,3)    &
                    +.5*rslope2(i,k,2)*rslope2(i,k,2)*rslope3(i,k,3)
            if(supcol.gt.0) then
              egs = exp(-0.09*supcol)
            else
              egs = 1.
            endif
            pgacs(i,k) = pi**2*egs*n0s*n0sfac(i,k)*n0g*abs(vt2g-vt2s)   &
                        *(dens/den(i,k))*acrfac
            pgacs(i,k) = min(pgacs(i,k),qrs(i,k,2)/dtcld)
          endif
          if(supcol.le.0) then
            xlf = xlf0
!-------------------------------------------------------------
! pseml: Enhanced melting of snow by accretion of water [HL A34]
!        (T>=T0: S->R)
!-------------------------------------------------------------
            if(qrs(i,k,2).gt.0.)                                        &
              pseml(i,k) = min(max(cliq*supcol*(psacw(i,k)+psacr(i,k))  &
                          /xlf,-qrs(i,k,2)/dtcld),0.)
!-------------------------------------------------------------
! pgeml: Enhanced melting of graupel by accretion of water [HL A24] [RH84 A21-A22]
!        (T>=T0: G->R)
!-------------------------------------------------------------
            if(qrs(i,k,3).gt.0.)                                        &
              pgeml(i,k) = min(max(cliq*supcol*(pgacw(i,k)+pgacr(i,k))  &
                          /xlf,-qrs(i,k,3)/dtcld),0.)
          endif
          if(supcol.gt.0) then
!-------------------------------------------------------------
! pidep: Deposition/Sublimation rate of ice [HDC 9]
!       (T<T0: V->I or I->V)
!-------------------------------------------------------------
            if(qci(i,k,2).gt.0.and.ifsat.ne.1) then
              pidep(i,k) = 4.*diameter*xni(i,k)*(rh(i,k,2)-1.)/work1(i,k,2)
              supice = satdt-prevp(i,k)
              if(pidep(i,k).lt.0.) then
                pidep(i,k) = max(max(pidep(i,k),satdt/2),supice)
                pidep(i,k) = max(pidep(i,k),-qci(i,k,2)/dtcld)
              else
                pidep(i,k) = min(min(pidep(i,k),satdt/2),supice)
              endif
              if(abs(prevp(i,k)+pidep(i,k)).ge.abs(satdt)) ifsat = 1
            endif
!-------------------------------------------------------------
! psdep: deposition/sublimation rate of snow [HDC 14]
!        (T<T0: V->S or S->V)
!-------------------------------------------------------------
            if(qrs(i,k,2).gt.0..and.ifsat.ne.1) then
              coeres = rslope2(i,k,2)*sqrt(rslope(i,k,2)*rslopeb(i,k,2))
              psdep(i,k) = (rh(i,k,2)-1.)*n0sfac(i,k)*(precs1          &
                           *rslope2(i,k,2)+precs2*work2(i,k)            &
                           *coeres)/work1(i,k,2)
              supice = satdt-prevp(i,k)-pidep(i,k)
              if(psdep(i,k).lt.0.) then
                psdep(i,k) = max(psdep(i,k),-qrs(i,k,2)/dtcld)
                psdep(i,k) = max(max(psdep(i,k),satdt/2),supice)
              else
                psdep(i,k) = min(min(psdep(i,k),satdt/2),supice)
              endif
              if(abs(prevp(i,k)+pidep(i,k)+psdep(i,k)).ge.abs(satdt))  &
                ifsat = 1
            endif
!-------------------------------------------------------------
! pgdep: deposition/sublimation rate of graupel [HL A21] [LFO 46]
!        (T<T0: V->G or G->V)
!-------------------------------------------------------------
            if(qrs(i,k,3).gt.0..and.ifsat.ne.1) then
              coeres = rslope2(i,k,3)*sqrt(rslope(i,k,3)*rslopeb(i,k,3))
              pgdep(i,k) = (rh(i,k,2)-1.)*(precg1*rslope2(i,k,3)       &
                              +precg2*work2(i,k)*coeres)/work1(i,k,2)
              supice = satdt-prevp(i,k)-pidep(i,k)-psdep(i,k)
              if(pgdep(i,k).lt.0.) then
                pgdep(i,k) = max(pgdep(i,k),-qrs(i,k,3)/dtcld)
                pgdep(i,k) = max(max(pgdep(i,k),satdt/2),supice)
              else
                pgdep(i,k) = min(min(pgdep(i,k),satdt/2),supice)
              endif
              if(abs(prevp(i,k)+pidep(i,k)+psdep(i,k)+pgdep(i,k)).ge. &
                abs(satdt)) ifsat = 1
            endif
!-------------------------------------------------------------
! pigen: generation(nucleation) of ice from vapor [HL 50] [HDC 7-8]
!       (T<T0: V->I)
!-------------------------------------------------------------
            if(supsat.gt.0.and.ifsat.ne.1) then
              supice = satdt-prevp(i,k)-pidep(i,k)-psdep(i,k)-pgdep(i,k)
              xni0 = 1.e3*exp(0.1*supcol)
              roqi0 = 4.92e-11*xni0**1.33
              pigen(i,k) = max(0.,(roqi0/den(i,k)-max(qci(i,k,2),0.))    &
                         /dtcld)
              pigen(i,k) = min(min(pigen(i,k),satdt),supice)
            endif
!
!-------------------------------------------------------------
! psaut: conversion(aggregation) of ice to snow [HDC 12]
!        (T<T0: I->S)
!-------------------------------------------------------------
            if(qci(i,k,2).gt.0.) then
              qimax = roqimax/den(i,k)
              psaut(i,k) = max(0.,(qci(i,k,2)-qimax)/dtcld)
            endif
!
!-------------------------------------------------------------
! pgaut: conversion(aggregation) of snow to graupel [HL A4] [LFO 37]
!        (T<T0: S->G)
!-------------------------------------------------------------
            if(qrs(i,k,2).gt.0.) then
              alpha2 = 1.e-3*exp(0.09*(-supcol))
              pgaut(i,k) = min(max(0.,alpha2*(qrs(i,k,2)-qs0))         &
                           ,qrs(i,k,2)/dtcld)
            endif
          endif
!
!-------------------------------------------------------------
! psevp: Evaporation of melting snow [HL A35] [RH83 A27]
!       (T>=T0: S->V)
!-------------------------------------------------------------
          if(supcol.lt.0.) then
            if(qrs(i,k,2).gt.0..and.rh(i,k,1).lt.1.) then
              coeres = rslope2(i,k,2)*sqrt(rslope(i,k,2)*rslopeb(i,k,2))
              psevp(i,k) = (rh(i,k,1)-1.)*n0sfac(i,k)*(precs1            &
                           *rslope2(i,k,2)+precs2*work2(i,k)            &
                           *coeres)/work1(i,k,1)
              psevp(i,k) = min(max(psevp(i,k),-qrs(i,k,2)/dtcld),0.)
            endif
!-------------------------------------------------------------
! pgevp: Evaporation of melting graupel [HL A25] [RH84 A19]
!       (T>=T0: G->V)
!-------------------------------------------------------------
            if(qrs(i,k,3).gt.0..and.rh(i,k,1).lt.1.) then
              coeres = rslope2(i,k,3)*sqrt(rslope(i,k,3)*rslopeb(i,k,3))
              pgevp(i,k) = (rh(i,k,1)-1.)*(precg1*rslope2(i,k,3)         &
                         +precg2*work2(i,k)*coeres)/work1(i,k,1)
              pgevp(i,k) = min(max(pgevp(i,k),-qrs(i,k,3)/dtcld),0.)
            endif
          endif
        enddo
      enddo
!
!
!----------------------------------------------------------------
!     check mass conservation of generation terms and feedback to the
!     large scale
!
      do k = 1, lm
        do i = ims, ime
!
          delta2=0.
          delta3=0.
          if(qrs(i,k,1).lt.1.e-4.and.qrs(i,k,2).lt.1.e-4) delta2=1.
          if(qrs(i,k,1).lt.1.e-4) delta3=1.
          if(t(i,k).le.t0c) then
!
!     cloud water
!
            value = max(qmin,qci(i,k,1))
            source = (praut(i,k)+pracw(i,k)+psacw(i,k)+pgacw(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              praut(i,k) = praut(i,k)*factor
              pracw(i,k) = pracw(i,k)*factor
              psacw(i,k) = psacw(i,k)*factor
              pgacw(i,k) = pgacw(i,k)*factor
            endif
!
!     cloud ice
!
            value = max(qmin,qci(i,k,2))
            source = (psaut(i,k)-pigen(i,k)-pidep(i,k)+praci(i,k)        &
                    +psaci(i,k)+pgaci(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              psaut(i,k) = psaut(i,k)*factor
              pigen(i,k) = pigen(i,k)*factor
              pidep(i,k) = pidep(i,k)*factor
              praci(i,k) = praci(i,k)*factor
              psaci(i,k) = psaci(i,k)*factor
              pgaci(i,k) = pgaci(i,k)*factor
            endif
!
!     rain
!
            value = max(qmin,qrs(i,k,1))
            source = (-praut(i,k)-prevp(i,k)-pracw(i,k)+piacr(i,k)    &
                    +psacr(i,k)+pgacr(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              praut(i,k) = praut(i,k)*factor
              prevp(i,k) = prevp(i,k)*factor
              pracw(i,k) = pracw(i,k)*factor
              piacr(i,k) = piacr(i,k)*factor
              psacr(i,k) = psacr(i,k)*factor
              pgacr(i,k) = pgacr(i,k)*factor
            endif
!
!     snow
!
            value = max(qmin,qrs(i,k,2))
            source = -(psdep(i,k)+psaut(i,k)-pgaut(i,k)+psacw(i,k)      &
                     +piacr(i,k)*delta3+praci(i,k)*delta3               &
                     -pracs(i,k)*(1.-delta2)+psacr(i,k)*delta2          &
                     +psaci(i,k)-pgacs(i,k) )*dtcld
            if (source.gt.value) then
              factor = value/source
              psdep(i,k) = psdep(i,k)*factor
              psaut(i,k) = psaut(i,k)*factor
              pgaut(i,k) = pgaut(i,k)*factor
              psacw(i,k) = psacw(i,k)*factor
              piacr(i,k) = piacr(i,k)*factor
              praci(i,k) = praci(i,k)*factor
              psaci(i,k) = psaci(i,k)*factor
              pracs(i,k) = pracs(i,k)*factor
              psacr(i,k) = psacr(i,k)*factor
              pgacs(i,k) = pgacs(i,k)*factor
            endif
!
!     graupel
!
            value = max(qmin,qrs(i,k,3))
            source = -(pgdep(i,k)+pgaut(i,k)               &
                     +piacr(i,k)*(1.-delta3)+praci(i,k)*(1.-delta3)     &
                     +psacr(i,k)*(1.-delta2)+pracs(i,k)*(1.-delta2)     &
                     +pgaci(i,k)+pgacw(i,k)+pgacr(i,k)+pgacs(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              pgdep(i,k) = pgdep(i,k)*factor
              pgaut(i,k) = pgaut(i,k)*factor
              piacr(i,k) = piacr(i,k)*factor
              praci(i,k) = praci(i,k)*factor
              psacr(i,k) = psacr(i,k)*factor
              pracs(i,k) = pracs(i,k)*factor
              pgacw(i,k) = pgacw(i,k)*factor
              pgaci(i,k) = pgaci(i,k)*factor
              pgacr(i,k) = pgacr(i,k)*factor
              pgacs(i,k) = pgacs(i,k)*factor
            endif
!
            work2(i,k)=-(prevp(i,k)+psdep(i,k)+pgdep(i,k)+pigen(i,k)   &
                       +pidep(i,k))
!     update
            q(i,k) = q(i,k)+work2(i,k)*dtcld
            qci(i,k,1) = max(qci(i,k,1)-(praut(i,k)+pracw(i,k)          &
                           +psacw(i,k)+pgacw(i,k))*dtcld,0.)
            qrs(i,k,1) = max(qrs(i,k,1)+(praut(i,k)+pracw(i,k)          &
                           +prevp(i,k)-piacr(i,k)-pgacr(i,k)            &
                           -psacr(i,k))*dtcld,0.)
            qci(i,k,2) = max(qci(i,k,2)-(psaut(i,k)+praci(i,k)          &
                           +psaci(i,k)+pgaci(i,k)-pigen(i,k)-pidep(i,k))   &
                           *dtcld,0.)
            qrs(i,k,2) = max(qrs(i,k,2)+(psdep(i,k)+psaut(i,k)+psacw(i,k)  &
                           -pgaut(i,k)+piacr(i,k)*delta3                &
                           +praci(i,k)*delta3+psaci(i,k)-pgacs(i,k)      &
                           -pracs(i,k)*(1.-delta2)+psacr(i,k)*delta2)    &
                           *dtcld,0.)
            qrs(i,k,3) = max(qrs(i,k,3)+(pgdep(i,k)+pgaut(i,k)         &
                           +piacr(i,k)*(1.-delta3)            &
                           +praci(i,k)*(1.-delta3)+psacr(i,k)*(1.-delta2)&
                           +pracs(i,k)*(1.-delta2)+pgaci(i,k)+pgacw(i,k) &
                           +pgacr(i,k)+pgacs(i,k))*dtcld,0.)
            xlf = xls-xl(i,k)
            xlwork2 = -xls*(psdep(i,k)+pgdep(i,k)+pidep(i,k)+pigen(i,k)) &
                      -xl(i,k)*prevp(i,k)-xlf*(piacr(i,k)+psacw(i,k)    &
                      +pgacw(i,k)+pgacr(i,k)+psacr(i,k))
            t(i,k) = t(i,k)-xlwork2/cpm(i,k)*dtcld*c88
          else
!
!     cloud water
!
            value = max(qmin,qci(i,k,1))
            source=(praut(i,k)+pracw(i,k)+psacw(i,k)+pgacw(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              praut(i,k) = praut(i,k)*factor
              pracw(i,k) = pracw(i,k)*factor
              psacw(i,k) = psacw(i,k)*factor
              pgacw(i,k) = pgacw(i,k)*factor
            endif
!
!     rain
!
            value = max(qmin,qrs(i,k,1))
            source = (-psacw(i,k)-praut(i,k)+pseml(i,k)+pgeml(i,k)     &
                     -pracw(i,k)-pgacw(i,k)-prevp(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              praut(i,k) = praut(i,k)*factor
              prevp(i,k) = prevp(i,k)*factor
              pgacw(i,k) = pgacw(i,k)*factor
              pracw(i,k) = pracw(i,k)*factor
              psacw(i,k) = psacw(i,k)*factor
              pseml(i,k) = pseml(i,k)*factor
              pgeml(i,k) = pgeml(i,k)*factor
            endif
!
!     snow
!
            value = max(qcrmin,qrs(i,k,2))
            source=(pgacs(i,k)-pseml(i,k)-psevp(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              pgacs(i,k) = pgacs(i,k)*factor
              psevp(i,k) = psevp(i,k)*factor
              pseml(i,k) = pseml(i,k)*factor
            endif
!
!     graupel
!
            value = max(qcrmin,qrs(i,k,3))
            source=-(pgacs(i,k)+pgevp(i,k)+pgeml(i,k))*dtcld
            if (source.gt.value) then
              factor = value/source
              pgacs(i,k) = pgacs(i,k)*factor
              pgevp(i,k) = pgevp(i,k)*factor
              pgeml(i,k) = pgeml(i,k)*factor
            endif
            work2(i,k)=-(prevp(i,k)+psevp(i,k)+pgevp(i,k))
!     update
            q(i,k) = q(i,k)+work2(i,k)*dtcld
            qci(i,k,1) = max(qci(i,k,1)-(praut(i,k)+pracw(i,k)         &
                    +psacw(i,k)+pgacw(i,k))*dtcld,0.)
            qrs(i,k,1) = max(qrs(i,k,1)+(praut(i,k)+pracw(i,k)         &
                    +prevp(i,k)+psacw(i,k)+pgacw(i,k)-pseml(i,k)       &
                    -pgeml(i,k))*dtcld,0.)
            qrs(i,k,2) = max(qrs(i,k,2)+(psevp(i,k)-pgacs(i,k)           &
                    +pseml(i,k))*dtcld,0.)
            qrs(i,k,3) = max(qrs(i,k,3)+(pgacs(i,k)+pgevp(i,k)           &
                    +pgeml(i,k))*dtcld,0.)
            xlf = xls-xl(i,k)
            xlwork2 = -xl(i,k)*(prevp(i,k)+psevp(i,k)+pgevp(i,k))        &
                      -xlf*(pseml(i,k)+pgeml(i,k))
            t(i,k) = t(i,k)-xlwork2/cpm(i,k)*dtcld*c88
          endif
! save delta2 and delta3 
!
          fdelta2(i,k) = delta2
          fdelta3(i,k) = delta3
        enddo
      enddo
!
! Inline expansion for fpvs
!         qs(i,k,1) = fpvs(t(i,k),0,rd,rv,cpv,cliq,cice,xlv0,xls,psat,t0c)
!         qs(i,k,2) = fpvs(t(i,k),1,rd,rv,cpv,cliq,cice,xlv0,xls,psat,t0c)
      hsub = xls
      hvap = xlv0
      cvap = cpv
      ttp=t0c+0.01
      dldt=cvap-cliq
      xa=-dldt/rv
      xb=xa+hvap/(rv*ttp)
      dldti=cvap-cice
      xai=-dldti/rv
      xbi=xai+hsub/(rv*ttp)
      do k = 1, lm
        do i = ims, ime
          tr=ttp/t(i,k)

          qs(i,k,1)=psat*exp(log(tr)*(xa))*exp(xb*(1.-tr))
          qs(i,k,1) = ep2 * qs(i,k,1) / (p(i,k) - qs(i,k,1))
          qs(i,k,1) = max(qs(i,k,1),qmin)
          tr=ttp/t(i,k)
          if(t(i,k).lt.ttp) then
            qs(i,k,2)=psat*exp(log(tr)*(xai))*exp(xbi*(1.-tr))
          else
            qs(i,k,2)=psat*exp(log(tr)*(xa))*exp(xb*(1.-tr))
          endif
          qs(i,k,2) = ep2 * qs(i,k,2) / (p(i,k) - qs(i,k,2))
          qs(i,k,2) = max(qs(i,k,2),qmin)
        enddo
      enddo
!
!----------------------------------------------------------------
!  pcond: condensational/evaporational rate of cloud water [HL A46] [RH83 A6]
!     if there exists additional water vapor condensated/if
!     evaporation of cloud water is not enough to remove subsaturation
!
      do k = 1, lm
        do i = ims, ime
          work1(i,k,1) = conden(t(i,k),q(i,k),qs(i,k,1),xl(i,k),cpm(i,k))
          work2(i,k) = qci(i,k,1)+work1(i,k,1)
          pcond(i,k) = min(max(work1(i,k,1)/dtcld,0.),max(q(i,k),0.)/dtcld)
          if(qci(i,k,1).gt.0..and.work1(i,k,1).lt.0.)                   &
            pcond(i,k) = max(work1(i,k,1),-qci(i,k,1))/dtcld
          q(i,k) = q(i,k)-pcond(i,k)*dtcld
          qci(i,k,1) = max(qci(i,k,1)+pcond(i,k)*dtcld,0.)
          t(i,k) = t(i,k)+pcond(i,k)*xl(i,k)/cpm(i,k)*dtcld*c88
        enddo
      enddo
      do k = 1, lm
        do i = ims, ime

            psmlt2d(i,k)= psmlt2d(i,k) + psmlt(i,k)
            pgmlt2d(i,k)= pgmlt2d(i,k) + pgmlt(i,k)
            pimlt2d(i,k)= pimlt2d(i,k) + pimlt(i,k)
            pihmf2d(i,k)= pihmf2d(i,k) + pihmf(i,k)
            pihtf2d(i,k)= pihtf2d(i,k) + pihtf(i,k)
            pgfrz2d(i,k)= pgfrz2d(i,k) + pgfrz(i,k)
            praut2d(i,k)= praut2d(i,k) + praut(i,k)*dtcld
            pracw2d(i,k)= pracw2d(i,k) + pracw(i,k)*dtcld
            praci2d_s(i,k)=praci2d_s(i,k)+praci(i,k)*dtcld*fdelta3(i,k)
            praci2d_g(i,k)=praci2d_g(i,k)+praci(i,k)*dtcld*(1-fdelta3(i,k))
            piacr2d_s(i,k)=piacr2d_s(i,k)+piacr(i,k)*dtcld*fdelta3(i,k)
            piacr2d_g(i,k)=piacr2d_g(i,k)+piacr(i,k)*dtcld*(1-fdelta3(i,k))
            psaci2d(i,k)= psaci2d(i,k) + psaci(i,k)*dtcld
            pgaci2d(i,k)= pgaci2d(i,k) + pgaci(i,k)*dtcld
            psacw2d(i,k)= psacw2d(i,k) + psacw(i,k)*dtcld
            pgacw2d(i,k)= pgacw2d(i,k) + pgacw(i,k)*dtcld
            pracs2d(i,k)= pracs2d(i,k) + pracs(i,k)*dtcld*(1-fdelta2(i,k))
            psacr2d_s(i,k)=psacr2d_s(i,k)+psacr(i,k)*dtcld*fdelta2(i,k)
            psacr2d_g(i,k)=psacr2d_g(i,k)+psacr(i,k)*dtcld*(1-fdelta2(i,k))
            pgacr2d(i,k)= pgacr2d(i,k) + pgacr(i,k)*dtcld
            pseml2d(i,k)= pseml2d(i,k) + pseml(i,k)*dtcld
            pgeml2d(i,k)= pgeml2d(i,k) + pgeml(i,k)*dtcld
            pigen2d(i,k)= pigen2d(i,k) + pigen(i,k)*dtcld
            psaut2d(i,k)= psaut2d(i,k) + psaut(i,k)*dtcld
            pgaut2d(i,k)= pgaut2d(i,k) + pgaut(i,k)*dtcld
            psevp2d(i,k)= psevp2d(i,k) + psevp(i,k)*dtcld
            pgevp2d(i,k)= pgevp2d(i,k) + pgevp(i,k)*dtcld
            if(prevp(i,k).lt.0)then
              prevp2d_e(i,k)=prevp2d_e(i,k)-prevp(i,k)*dtcld
            elseif(prevp(i,k).gt.0)then
              prevp2d_c(i,k)=prevp2d_c(i,k)+prevp(i,k)*dtcld
            endif

            if(psdep(i,k).lt.0)then
              psdep2d_s(i,k)=psdep2d_s(i,k)-psdep(i,k)*dtcld
            elseif(psdep(i,k).gt.0)then
              psdep2d_d(i,k)=psdep2d_d(i,k)+psdep(i,k)*dtcld
            endif
            if(pgdep(i,k).lt.0)then
              pgdep2d_s(i,k)=pgdep2d_s(i,k)-pgdep(i,k)*dtcld
            elseif(pgdep(i,k).gt.0)then
              pgdep2d_d(i,k)=pgdep2d_d(i,k)+pgdep(i,k)*dtcld
            endif

            if(pidep(i,k).lt.0)then
              pidep2d_s(i,k)=pidep2d_s(i,k)-pidep(i,k)*dtcld
            elseif(pidep(i,k).gt.0)then
              pidep2d_d(i,k)=pidep2d_d(i,k)+pidep(i,k)*dtcld
            endif


!          ------------------------------------------------------
!           Separate condensation and evaporation of cloud water
!           to be consistent with the Detailed mp (12/16/2008 ki)
!          ------------------------------------------------------

!
!                  condensation/evaporation of cloud drops
!
           if (pcond(i,k) .ge. 0.0) then      ! > 0 for deposition
              pcondc0(i,k) = pcondc0(i,k) + pcond(i,k) * dtcld
           else                               ! < 0 for evaporation of cloud drops
              pconde0(i,k) = pconde0(i,k) - pcond(i,k) * dtcld
           endif
      enddo
       enddo
!
!
!----------------------------------------------------------------
!     padding for small values
!
      do k = 1, lm
        do i = ims, ime
          if(qci(i,k,1).le.qmin) qci(i,k,1) = 0.0
          if(qci(i,k,2).le.qmin) qci(i,k,2) = 0.0
        enddo
      enddo
      enddo                  ! big loops
! W.wang 2/26/10 , still use rainncv,snowncv,graupelncv to OUT
      do i = ims, ime
         rainncv(i) = rain_dt(i)
         snowncv(i) = snow_dt(i)
         graupelncv(i) = graupel_dt(i)
         sr(i) = 0.0
         if(rainncv(i) .gt. 0) sr(i) = (snowncv(i)+graupelncv(i))/(rainncv(i)+1.0e-12)
      enddo

  END SUBROUTINE wsm62d
! ...................................................................
      REAL FUNCTION rgmma(x)
!-------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------
!     rgmma function:  use infinite product form
      REAL :: euler
      PARAMETER (euler=0.577215664901532)
      REAL :: x, y
      INTEGER :: i
      if(x.eq.1.)then
        rgmma=0.
          else
        rgmma=x*exp(euler*x)
        do i=1,10000
          y=float(i)
          rgmma=rgmma*(1.000+x/y)*exp(-x/y)
        enddo
        rgmma=1./rgmma
      endif
      END FUNCTION rgmma
!
!--------------------------------------------------------------------------
      REAL FUNCTION fpvs(t,ice,rd,rv,cvap,cliq,cice,hvap,hsub,psat,t0c)
!--------------------------------------------------------------------------
      IMPLICIT NONE
!--------------------------------------------------------------------------
      REAL t,rd,rv,cvap,cliq,cice,hvap,hsub,psat,t0c,dldt,xa,xb,dldti,  &
           xai,xbi,ttp,tr
      INTEGER ice
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ttp=t0c+0.01
      dldt=cvap-cliq
      xa=-dldt/rv
      xb=xa+hvap/(rv*ttp)
      dldti=cvap-cice
      xai=-dldti/rv
      xbi=xai+hsub/(rv*ttp)
      tr=ttp/t
      if(t.lt.ttp.and.ice.eq.1) then
        fpvs=psat*(tr**xai)*exp(xbi*(1.-tr))
      else
        fpvs=psat*(tr**xa)*exp(xb*(1.-tr))
      endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END FUNCTION fpvs
!-------------------------------------------------------------------
  SUBROUTINE wsm6init(den0,denr,dens,cl,cpv)
!-------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------
!.... constants which may not be tunable
   REAL, INTENT(IN) :: den0,denr,dens,cl,cpv
   REAL :: pi
!
   pi = 4.*atan(1.)
   xlv1 = cl-cpv
!
   qc0  = 4./3.*pi*denr*r0**3*xncr/den0  ! 0.419e-3 -- .61e-3
   qck1 = .104*9.8*peaut/(xncr*denr)**(1./3.)/xmyu*den0**(4./3.) ! 7.03
!
   bvtr1 = 1.+bvtr
   bvtr2 = 2.5+.5*bvtr
   bvtr3 = 3.+bvtr
   bvtr4 = 4.+bvtr
   bvtr6 = 6.+bvtr
   g1pbr = rgmma(bvtr1)
   g3pbr = rgmma(bvtr3)
   g4pbr = rgmma(bvtr4)            ! 17.837825
   g6pbr = rgmma(bvtr6)
   g5pbro2 = rgmma(bvtr2)          ! 1.8273
   pvtr = avtr*g4pbr/6.
   eacrr = 1.0
   pacrr = pi*n0r*avtr*g3pbr*.25*eacrr
   precr1 = 2.*pi*n0r*.78
   precr2 = 2.*pi*n0r*.31*avtr**.5*g5pbro2
   xm0  = (di0/dicon)**2
   xmmax = (dimax/dicon)**2
   roqimax = 2.08e22*dimax**8
!
   bvts1 = 1.+bvts
   bvts2 = 2.5+.5*bvts
   bvts3 = 3.+bvts
   bvts4 = 4.+bvts
   g1pbs = rgmma(bvts1)    !.8875
   g3pbs = rgmma(bvts3)
   g4pbs = rgmma(bvts4)    ! 12.0786
   g5pbso2 = rgmma(bvts2)
   pvts = avts*g4pbs/6.
   pacrs = pi*n0s*avts*g3pbs*.25
   precs1 = 4.*n0s*.65
   precs2 = 4.*n0s*.44*avts**.5*g5pbso2
   pidn0r =  pi*denr*n0r
   pidn0s =  pi*dens*n0s
!
   pacrc = pi*n0s*avts*g3pbs*.25*eacrc
!
   bvtg1 = 1.+bvtg
   bvtg2 = 2.5+.5*bvtg
   bvtg3 = 3.+bvtg
   bvtg4 = 4.+bvtg
   g1pbg = rgmma(bvtg1)
   g3pbg = rgmma(bvtg3)
   g4pbg = rgmma(bvtg4)
   pacrg = pi*n0g*avtg*g3pbg*.25
   g5pbgo2 = rgmma(bvtg2)
   pvtg = avtg*g4pbg/6.
   precg1 = 2.*pi*n0g*.78
   precg2 = 2.*pi*n0g*.31*avtg**.5*g5pbgo2
   pidn0g =  pi*deng*n0g
!
   rslopermax = 1./lamdarmax
   rslopesmax = 1./lamdasmax
   rslopegmax = 1./lamdagmax
   rsloperbmax = rslopermax ** bvtr
   rslopesbmax = rslopesmax ** bvts
   rslopegbmax = rslopegmax ** bvtg
   rsloper2max = rslopermax * rslopermax
   rslopes2max = rslopesmax * rslopesmax
   rslopeg2max = rslopegmax * rslopegmax
   rsloper3max = rsloper2max * rslopermax
   rslopes3max = rslopes2max * rslopesmax
   rslopeg3max = rslopeg2max * rslopegmax
!
  END SUBROUTINE wsm6init

          subroutine vrec(y,x,n)   ! copied from wrf frame/libmassv.F
          real x(*),y(*)
          do j=1,n
           y(j)=1.0/x(j)
          enddo
           end subroutine vrec
 
          subroutine vsqrt(y,x,n)
          real x(*),y(*)
          do j=1,n
           y(j)=sqrt(x(j))
          enddo
          end subroutine vsqrt



END MODULE module_mp_wsm6
