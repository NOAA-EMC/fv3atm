cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE solar_heat(np,nps,O,O2,N2,HO,HO2,HN2,effeuv,effuv,     
     & F107, COSPASS,sheat,sh1,sh2)
!-------------------------------------------------------------------------
! calculate solar heating from Tim Fuller-Rowell
!-------------------------------------------------------------------------
!c  **
!c  calculates solar heating from EUV and SRC wavelengths
!c  assumes a heating efficiency profile on pressure levels
!c  code was written in SI units
!c  Input:
!c  O atomic oxygen number density profile m-3
!c  O2 molecular oxygen number density profile m-3
!c  N2 molecular nitrogen number density profile m-3
!c  HO atomic oxygen scale height profile m
!c  HO2 molecular oxygen scale height profile m
!c  HN2 molecular nitrogen scale height profile m
!c  solar flux F10.7
!c  COSPASS cosine of solar zenith angle
!c  Output:
!c  SHEAT heating rate profile J/m3
!---------------------------------------------------------------------------
      implicit none
      integer, intent(in) :: np    ! number of pressure levels
      integer, intent(in) :: nps    !pressure index to start  
      real, intent(in)    :: o(np),o2(np),n2(np) ! number density/m3
      real, intent(in)    :: ho(np),ho2(np),hn2(np) ! scale height(m)
      real, intent(in)    :: effeuv(np),effuv(np) !heating efficiency 
      real, intent(in)    :: f107    !f10.7cm 
      real, intent(in)    :: cospass !cos zenith angle
      real, intent(out)   :: sheat(np),sh1(np),sh2(np)   !J/m3 heating rate
      real SO(np),SO2(np),SN2(np),                                      
     &SFL(57),CSAO(57),CSAO2(57),CSAN2(57),CSIO(57),                    
     &CSIO2(57),CSIN2(57),SFH(57),SF(57),PAEUV(np,65),SFUV(8),          
     &A(8),UVXS(8),RLAM(65)
      real coschi,rnight,seco,seco2,secn2,wo,wo2,wn2,tau,tauo,tauo2,    
     & taun2,pcc
      integer i,j,jj
!c  **
!c  number of pressure levels to process for solar heating
!c  pressure levels defined by pressure(n)=5.2285*exp(1-n)
!c  **
!c  wavelength/energy conversion SI E=hc/lamda
      PCC=1.985E-25
!c  **
!C  WAVELENGTHS Angstroms
!C  **
      DATA RLAM/18.6,19.0,21.6,21.8,22.1,28.5,28.8,29.5,30.0,           
     &30.4,33.7,41.0,43.8,44.0,44.2,45.7,46.4,46.7,47.9,49.2,           
     &75.,125.,175.,225.,256.3,284.15,275.,303.31,303.78,               
     &325.,368.07,375.,425.,465.22,475.,525.,554.37,584.33,             
     &575.,609.76,629.73,625.,675.,730.36,725.,765.15,770.41,           
     &789.36,775.,825.,875.,925.,977.62,975.,1025.72,1031.91,           
     &1025.,1387.5,1425.,1475.,1525.,1575.,1625.,1675.,1725./
!C  **
!C  REVISED FLUXES BY TORR AND TORR 85 JGR 90 6675
!C  WITH ADDITIONAL VALUES 1 TO 20 FOR WAVELENGTHS BELOW 50A.
!C  **
           DATA SFL/
     &.0001,.0001,.0003,.0001,.0003,.0005,.0025,.0022,.0012,            
     &.0006,.0011,.0006,.0021,.0008,.0009,.0005,.0027,.0052,            
     &.0059,.0043,                                                      
     &.38,.13,1.84,.92,.27,.1,.84,.24,6.,.87,.74,.21,.39,.18,           
     &.31,.51,.80,1.58,.48,.45,1.5,.17,.22,.39,.17,.2,.24,              
     &.79,.87,1.93,4.43,4.22,5.96,1.79,4.38,3.18,3.64/
!c  **
           DATA SFH/                                                    
     &.0016,.0053,.0048,.0016,.0048,.0072,.0211,.0186,.0024,            
     &.0104,.0158,.0073,.0130,.0097,.0109,.0061,.0168,.0107,            
     &.0121,.0267,                                                      
     &1.37,.468,5.7,7.14,1.08,5.72,12.16,4.69,14.39,6.83,1.53,          
     &2.54,1.53,.736,1.82,1.64,1.52,4.3,1.048,2.48,3.87,1.37,           
     &.539,.746,.429,.439,1.19,1.514,2.454,4.85,12.219,9.85,            
     &10.217,4.078,11.85,6.1,6.09/           
!c  **
!c  UV FLUX IN SRC AND O2 X-SECTIONS FROM M.R.TORR ET AL 
!c  JGR 1980 6063
!c  **
      DATA A/9.73,17.93,27.38,51.57,70.99,97.4,205.,374.24/
      DATA UVXS/1.2E-17,1.5E-17,1.3E-17,1.0E-17,6.0E-18,3.4E-18,        
     &1.5E-18,5.0E-19/
!C  **
!C  REVISED VALUES FROM SAMSON AND PAREEK 85
!C  **
           DATA CSAO/                                                   
     &.34,.36,.5,.51,.52,.05,.05,.06,.06,.06,.08,.13,.15,               
     &.15,.16,.17,.18,.18,.19,.21,                                      
     &0.7,1.7,3.0,5.1,6.2,7.3,7.0,7.7,7.7,8.5,10.,10.,                  
     &11.21,11.25,11.64,11.91,12.13,12.17,11.9,12.23,12.22,             
     &12.21,10.04,11.35,8.0,4.18,4.18,4.28,4.23,4.38,4.18,2.12,         
     &0.,0.,0.,0.,0./
!c  **
           DATA CSAO2/                                                  
     &.69,.72,.99,1.01,1.05,.10,.11,.11,.12,.12,.16,.26,.3,.31,         
     &.31,.34,.35,.36,.38,.41,                                          
     &1.18,3.61,7.27,10.5,12.8,14.8,13.65,15.98,16.,17.19,18.40,        
     &18.17,19.39,20.4,21.59,24.06,25.59,22.0,25.04,26.1,25.8,          
     &26.02,26.27,25.,29.05,21.98,25.18,26.66,27.09,20.87,9.85,         
     &15.54,4.0,16.53,1.6,1.0,1.1/
!c  **
           DATA CSAN2/                                                  
     &.44,.47,.65,.67,.69,1.13,1.13,1.12,1.11,1.10,.1,.16,.19,          
     &.19,.19,.21,.22,.22,.24,.25,                                      
     &.6,2.32,5.4,8.15,9.65,10.6,10.8,11.58,11.6,                       
     &14.6,18.0,17.51,21.07,21.8,                                       
     &21.85,24.53,24.69,23.2,22.38,23.1,23.2,23.22,29.75,26.3,          
     &30.94,35.46,26.88,19.26,30.71,15.05,46.63,16.99,.7,               
     &36.16,0.,0.,0./
!c  **
! VAY-2016 SECO corrections
!
      COSCHI=COSPASS
      rnight=1.0
      if(coschi.lt.0.07)then
        rnight=1.e-6
        SECO=1./0.07      ! was 1.0 in 2013
      else
        rnight=1.0
        SECO=1./COSCHI         
      end if
        SECO2=SECO
        SECN2=SECO

      do j=1,57
        SF(j)=1.e9*((SFH(j)-SFL(j))*F107/172.-0.413               
     &        *SFH(j)+1.413*SFL(j))
        if(sf(j).lt.0.0)sf(j)=0.0
      enddo
      do  j=1,8
        SFUV(j)=A(j)*1.E9*(0.00086*F107+0.94)
      enddo
      do i=nps,np
!c  **
!        SECO=1./COSCHI
!        SECO2=SECO
!        SECN2=SECO
        WO=O(i)*HO(i)*SECO*1.e-4
        WO2=O2(i)*HO2(i)*SECO2*1.e-4
        WN2=N2(i)*HN2(i)*SECN2*1.e-4
!c  **
!c  loop over all wavelengths bands
!c  **
        sheat(i)=0.0
        sh1(i)=0.
        sh2(i)=0.
        do j=1,57
          TAUO=CSAO(j)*WO*1.e-18
          TAUO2=CSAO2(j)*WO2*1.e-18
          TAUN2=CSAN2(j)*WN2*1.e-18
          TAU=TAUO+TAUO2+TAUN2
          PAEUV(i,j)=SF(j)*EXP(-TAU)*(CSAO(j)*O(i)+                     
     &    CSAO2(j)*O2(i)+CSAN2(j)*N2(i))*PCC*rnight*1.e-8/RLAM(j)
          sheat(i)=sheat(i)+paeuv(i,j)*effeuv(i)
          sh1(i)=sh1(i)+paeuv(i,j)*effeuv(i)
        enddo
!c  **
!c  add SRC channels
!c  ** 
        do j=58,65
          JJ=j-57
          TAU=UVXS(JJ)*WO2
          PAEUV(i,j)=SFUV(JJ)*EXP(-TAU)*UVXS(JJ)*O2(i)*PCC*rnight       
     &    /RLAM(j)*1.e10
          sheat(i)=sheat(i)+paeuv(i,j)*effuv(i)
          sh2(i)=sh2(i)+paeuv(i,j)*effuv(i)
        enddo
      enddo
      if(nps.ge.2) then
        do i=1,nps-1
          sheat(i)=0.
          sh1(i)=0.
          sh2(i)=0.
        enddo
      endif
      RETURN
      END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE COOLNO1(np,nps,T,O,NO,QNO)                   
!-------------------------------------------------------------------------
! calculate NO cooling from Tim Fuller-Rowell
!-------------------------------------------------------------------------
!c  **
!c  input:
!c  T temperature profile K
!c  O atomic oxygen number density profile m-3
!c  NO nitric oxide number density profile m-3
!c  output:
!c  QNO: NO cooling rate J/m-3
!c  **
      implicit none
      integer, intent(in):: np           !numer of pressure levels
      integer, intent(in):: nps          ! pressure index to start
      real, intent(in)   :: O(np),NO(np) !number density/m3        
      real, intent(in)   :: T(np)        !temp (K)        
      real, intent(out)  :: QNO(np)         
      real K10,HV,A10,BZ,A1,A2,A3,OM1,OM,G        
      integer i
      K10=3.6E-17                                                       
      HV=3.726E-20                                                      
      A10=13.3                                                          
      G=1.0                                                             
      BZ=1.38E-23                                                       
      A2=5.4E-6*(1./(EXP(HV/BZ/5800.)-1.))
      A3=0.5*EXP(-HV/BZ/247.5)
      do i=nps,np                                                   
        OM1=K10*O(i)          
        OM=OM1/(OM1+A10)
        A1=EXP(-HV/BZ/T(i)) 
        QNO(i)=HV*NO(i)*OM*A10*G*(A1-A2-A3)
      enddo
      if(nps.ge.2) then
        do i=1,nps-1                                                   
          QNO(i)=0.
        enddo
      endif
      return                                                            
      end
!
      subroutine z15toz(ain,levs, Z, aout,down)
! interpolate 15 pressure levels (from Tim's grid) to
! idea pressure grid pr(levs)
      use idea_composition, only : pr=> pr_idea
      implicit none
      integer, parameter :: np=15      !number of pressure levels of input
      integer, intent(in) :: levs      !number of pressure levels of output 
      real,    intent(in) :: Z(levs)   !model -log(pressure) grid
      real,    intent(in) :: ain(np)   !input field in 15 pressure grid
      real,    intent(in) :: down      !field value below 1.0376Pa => surface
      real,    intent(out):: aout(levs)!output in levs pressure grid
!local variable
      real p15(np),z15(np),dz
      integer kref,k,i
!
      do k=1,np
        p15(k)=1.0376*exp(1.-k)
        z15(k)=-log(p15(k))
      enddo
!
! grids
! interpolation
!
      do k=1,levs
        do i=1,np-1
          if(z(k).ge.z15(i).and.z(k).le.z15(i+1)) then
          dz=(z(k)-z15(i))/(z15(i+1)-z15(i))*(ain(i+1)-ain(i))
          aout(k)=ain(i) +dz
          endif
        enddo   
        if(z(k).lt.z15(1))  aout(k)=down    ! zero
        if(z(k).gt.z15(15)) aout(k)=ain(15)
    
      enddo
      return
      end
      subroutine z63toz(ain, levs, Z, aout,down)
! interpolate 63 pressure levels (from Tim's grid) to
! idea pressure grid pr(levs)
      use idea_composition, only : pr=> pr_idea
      implicit none
      integer, parameter  :: np=63     !number of pressure levels of input
      integer, intent(in) :: levs      !number of pressure levels of output 
      real,    intent(in) :: Z(levs)   ! model grid
      real,    intent(in) :: ain(np)   !input field in 63 pressure grid
      real,    intent(in) :: down      !field value under 6.9Pa
      real,    intent(out):: aout(levs)!output in levs pressure grid
!local variable
      real p63(np),z63(np),dz
      integer kref,k,i
!
      DATA p63/6.90775528,  6.57442194, 
     &   6.24108859,  5.90775525,  5.57442191,
     &   5.24108856,  4.90775522,  4.57442188,
     &   4.24108853,  3.90775519,  3.57442185,
     &   3.2410885,   2.90775516,  2.57442182,
     &   2.24108847,  1.90775513,  1.57442179,
     &   1.24108844,  0.9077551,   0.574421757,
     &   0.241088414, -0.0922449296,-0.425578273,
     &   -0.758911616,-1.09224496,  -1.4255783, 
     &   -1.75891165, -2.09224499,  -2.42557833,
     &   -2.75891168, -3.09224502,  -3.42557836,
     &   -3.75891171, -4.09224505,  -4.42557839,
     &   -4.75891174, -5.09224508,  -5.42557842,
     &   -5.75891177, -6.09224511,  -6.42557845,
     &   -6.75891179, -7.09224514,  -7.42557848,
     &   -7.75891182, -8.09224517,  -8.42557851,
     &   -8.75891185, -9.0922452,   -9.42557854,
     &   -9.75891188, -10.0922452,  -10.4255786,
     &   -10.7589119, -11.0922453,  -11.4255786,
     &   -11.7589119, -12.0922453,  -12.4255786,
     &   -12.758912,  -13.0922453,  -13.4255787,
     &   -13.758912/
      do k=1,np
        z63(k)=-p63(k)
      enddo
!
      do k=1,levs
        do i=1,np-1
          if(z(k).ge.z63(i).and.z(k).le.z63(i+1)) then
          dz=(z(k)-z63(i))/(z63(i+1)-z63(i))*(ain(i+1)-ain(i))
           aout(k)=dz+ain(i)
          endif
        enddo
        if(z(k).lt.z63(1))  aout(k)=down
        if(z(k).gt.z63(63)) aout(k)=ain(63)
      
      enddo
      return
      end
!
