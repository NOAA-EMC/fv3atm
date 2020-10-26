!-----------------------------------------------------------------------
!
                        module module_constants
!
!-----------------------------------------------------------------------
!
      use module_kinds
!
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!***  Physical constants
!-----------------------------------------------------------------------
real(kind=kfpt),parameter:: &
 a=6376000. &                ! radius of earth
,a2=17.2693882 &             ! saturation spec. humidity formula coeff.
,a3=273.15 &                 ! saturation spec. humidity formula coeff.
,a4=35.86 &                  ! saturation spec. humidity formula coeff.
,cp=1004.6 &                 ! spec. heat for dry air at constant pressure
,elwv=2.501e6 &              ! latent heat, liquid/vapor
,eliv=2.850e6 &              ! latent heat, ice/vapor
,eliwv=2.683e6 &             ! latent heat, mix ice/water - vapor
,elivw=2.72e6 &              ! another one
,eliw=3.50e5 &               ! latent heat of fusion from WRF
,epsilon=1.e-15 &
!! ,epsq2=0.02 &                ! floor value for 2tke
,epsq=1.e-12 &               ! floor value for specific humidity (kg/kg)
,g=9.8060226 &               ! gravity
,pi=3.141592653589793 &      ! ludolf number
,pihf=0.5*pi &               ! pi/2
,pq0=379.90516 &             ! water vapor pressure for tetens formula
,r=287.04 &                  ! gas constant for dry air
,r_d=287.04 &                ! gas constant for dry air
,r_v=461.6 &                 ! gas constant for water vapor
,cv=cp-r_d &
,cpv=4.*r_v &
,ep_1=r_v/r_d-1. &
,ep_2=r/r_v &
,cice=2106. &
,cliq=4190. &
,psat=610.78 &
,rhoair0=1.28 &
,rhowater=1000. &            ! density of water (kg/m3)
,rhosnow=100. &
,p608=r_v/r-1. &             ! factor for water vapor in virtual temperature
,svpt0=273.15 &
,ti0=271.15 &                ! water temperature below sea ice
,tiw=273.15 &                ! melting point
,twom=.00014584 &            ! 2 x angular velocity of earth
,cappa=r/cp &                ! kappa
,dtr=pi/180. &               ! factor converting degrees to radians
,rlag=14.8125 &       
,stbolt=5.67051E-8 &         ! Stefan-Boltzmann constant
,dbzmin=-20.  &              ! Minimum radar reflectivity (dBZ)
,xlf=3.50E5 &
,xlv=2.5E6
!-----------------------------------------------------------------------
!***  Soil layers
!-----------------------------------------------------------------------
integer(kind=kint),parameter :: &
 kmsc=6 &                    ! max # of LSM layers
,ksnoc=1 &                   ! number of snow layers in LSM scheme
,nosnoc=kmsc-ksnoc &         ! # of soil layers without sno
,nwetc=nosnoc-1              ! # of soil layers with soil moisture
!
                       end module module_constants
!
!-----------------------------------------------------------------------
