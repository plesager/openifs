SUBROUTINE TM5M7_DRYDEP &
 !---input
 & ( KIDIA , KFDIA   , KLON, KLEV, KTRAC, KAERO   , &
 &   PTSPHY, PTENCI,  PT      , PAP , PAPH, PWIND, &
 &   PLSM,PCI, PAERUST, PZ0M,PRHCL,PGEOH, PDZ,&
 &   SSHF, SLHF, &
 !---prognostic fields
 &   RW_MODE, DENS_MODE, PCEN, &
 !---output
 &   PFLUXAER, PTENC, PVDA)

!**** *TM5M7_DRYDEP* -  ROUTINE FOR PARAMETRIZATION OF TM5M7 AEROSOL DEPOSITION
!
! !DESCRIPTION:  This module calculates deposition of aerosol particles
!
! It is assumed that the distribution is described by nmodes log-normal
! distributions
!
! Each mode has a number and mass and a sigma_lognormal. Number and mass are
! related and the mean aerosol radius can thus be calculated for each mode.
! 
! mass=pi*4./3.*radius**3.*number*exp(9./2.*ln2s) /1e18*density, with:
!
!     density = density of aerosol type
!     ln2s    = (alog(sigma_lognormal(mode)))**2
!     mass    = mode mass
!     number  = mode number
!
!    mode1 : accumulation
!    mode2 : coarse
!    mode3 : super coarse (ss) coarse
!
! For each mode a separate fall velocity is calculated according to the mass
! and number mean radii. Water take-up by seasalt particles is taken into
! account. This changes the density, radius, and sigma of the distribution.
!
! Also included is the deposition calculation. based on a lookup table
! calculated for a reference aerosol density (e.g. 1800 kg/m3) and a number of
! radii. This deposition curve is convoluted with the number/volume
! distribution of the aerosols.
!

! Again, for SS the water takeup is accounted for, and the effects on density,
! sigma and radius are calculated. The density has effect on the impaction
! term is the depotion calculation. This can be modeled by a shift in the
! radius. Thus the radii of the lookup table are adapted for density
! differences when impaction becomes important.

!**   INTERFACE.
!     ----------
!          *TM5M7_DRYDEP* IS CALLED FROM *TM5M7_PHY2*.

! INPUTS:
! -------
! PTSPHY                : TIMESTEP                  (s)
! PTENCI(KLON,KLEV,KTRAC)    : INPUT TENDENCIES          (xx kg s-1)
! PAP     (KLON,KLEV)   : LEVEL PRESSURE            (Pa)
! PAPH    (KLON,KLEV+1) : HALF-LEVEL PRESSURE       (Pa)
! PWIND   (KLON)        : 10-meter wind speed
! PCI     (KLON)        : FRACTION OF SEA ICE
! PLSM    (KLON)        : LAND-SEA MASK   
! PAERUST  KLON)        : Friction velocity 
! PZ0M                  : roughness length for momentum              (m)
! PRHCL   (KLON)        : Rel. humidity at surface level (0-1)
! PT      (KLON)        : SURFACE LEVEL TEMPERATURE         (K)
! PCEN    (KLON,KLEV,KTRAC)   : CONCENTRATION OF TRACERS  (xx kg-1)

! OUTPUTS:
! --------
! PTENC  (KLON,KLEV,KTRAC)    : TENDENCY                  (xx kg-1 s-1)
! PVDA   (KLON,NACTAERO)      : dry deposition velocities (m/sec)

!     EXTERNALS.
!     ----------
!          NONE

!     MODIFICATIONS.
!     -------------

!     SWITCHES.
!     --------

!     MODEL PARAMETERS
!     ----------------

!-----------------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMLUN, ONLY : NULERR
USE YOEAERSRC ,ONLY : YREAERSRC
USE YOEAERSNK ,ONLY : YREAERSNK
USE YOMCST    ,ONLY : RG, RD,RNAVO,RPI
USE YOM_YGFL  ,ONLY : YGFL 
USE TM5M7_DATA, ONLY : NMOD, MODE_NM_SED, mode_start, mode_tracers_sed, XMAIR,MODAL_DATA, &
    & SIGMA_LOGNORMAL, NRDEP, LUR, DENSITY_REF
USE TM5M7_EMIS_DATA, ONLY : VKARMAN    

IMPLICIT NONE

!-----------------------------------------------------------------------

!*       0.1  ARGUMENTS
!             ---------

!---input fields
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KTRAC
INTEGER(KIND=JPIM),INTENT(IN)    :: KAERO(YGFL%NAERO)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTENCI(KLON,KLEV,KTRAC) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON) ! Temp
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV)! Mid-lev pres 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWIND(KLON),PLSM(KLON),PCI(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAERUST(KLON),PZ0M(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRHCL(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOH(KLON,0:KLEV),PDZ(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: SSHF(KLON),SLHF(KLON) ! sensible and latent heat fluxes (W/m2)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1) ! interface pres
TYPE(MODAL_DATA)  ,INTENT(IN)    :: RW_MODE(NMOD)
TYPE(MODAL_DATA)  ,INTENT(IN)    :: DENS_MODE(NMOD)

REAL(KIND=JPRB)   ,INTENT(IN)    :: PCEN(KLON,KLEV,KTRAC)
!---output fields
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFLUXAER(KLON,YGFL%NACTAERO)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENC(KLON,KLEV,KTRAC)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVDA(KLON,YGFL%NACTAERO) !eehol: correct indices

!---local fields
REAL(KIND=JPRB)     :: TEMP, PB, DP, ZVIS, RHO_AIR, TO_PASCAL,RADIUS,VT
REAL(KIND=JPRB)     :: UM,ZR, VD_SEA,CUNNING,DC,DENSAER,RELAX,SC
REAL(KIND=JPRB)     :: UST_LAND, ST_LAND,VB_LAND,VI_LAND, VKD_LAND, VD_LAND
REAL(KIND=JPRB)     :: FREESEA, UST_SEA, ST_SEA, RE, VBSEA,VISEA, VKDACCSEA, VKD_SEA

REAL(KIND=JPRB)     :: ZRAERO(KLON)
REAL(KIND=JPRB)     :: USTAR_SEA(KLON), SR_SEA(KLON), ALPHA(KLON), BUBBLE(KLON),ALPHAE(KLON)
REAL(KIND=JPRB)     :: Y0,YRA, ZRA,OBUK, BUOY, TSTV, Z1TSPHY, ZAERI, ZAERO

INTEGER(KIND=JPIM)  :: JK, JL, JN, JNN, MODE, IMODE, INMODE, IRDEP, ITN, NRD, JAER

REAL(KIND=JPRB)      :: m_to_pa  !factor from m/s --> Pa/s 

REAL(KIND=JPRB), PARAMETER :: dynvisc=1.789e-4*2. ! g cm-1 s-1 CHECKED FD OK, there is temp. dependence Perry p. 3-248.
				                  ! unit is g cm-1 s+1 FD 
				                  ! checked with Seinfeld---> factor 2. came out (diameter ? radius)
REAL(KIND=JPRB), PARAMETER :: cl=0.066*1e-4	  ! mean free path [cm] (particle size also in cm)
REAL(KIND=JPRB), PARAMETER :: bc= 1.38e-16	  ! boltzman constant [g cm-2 s-1 K-1] (1.38e-23 J deg-1) =>binas
REAL(KIND=JPRB), PARAMETER :: kappa=1.  	  ! shapefactor
REAL(KIND=JPRB), PARAMETER :: visc=0.15 	  ! KINEMATIC molecular viscocity [cm2 s-1] 
                                                  ! this is also function of temperature FD


REAL(KIND=JPRB),PARAMETER    ::  eff=0.5		 ! parameters needed for sea/bubble formation
REAL(KIND=JPRB),PARAMETER    ::  rdrop=0.005	 ! cm
REAL(KIND=JPRB),PARAMETER    ::  zdrop=10.0	 ! cm
REAL(KIND=JPRB),PARAMETER    ::  eps=0.6		 ! parameter related to bubble formation
REAL(KIND=JPRB)		     ::  qdrop,ZS,phi,alpha1,vk1,vk2,alpharat   ! auxiliary parameters for bubble formation
REAL(KIND=JPRB)		     ::  sr_help         ! surface roughnes consistent with 10 m wind.


      ! Href: constant reference height for calculations of raero
      real, parameter :: Href=30.     
      real,parameter  :: rhoCp          = 1231.0
      real,parameter  :: rhoLv          = 3013000.0


REAL(KIND=JPRB),DIMENSION(KLON,NRDEP)  :: ZVDA

REAL(KIND=JPRB)    :: LNS (NMOD)

REAL(KIND=JPRB),DIMENSION(:),POINTER   :: VS  



REAL(KIND=JPRB),  DIMENSION(nrdep)	     :: d_aer	      ! diameter vd loopup table (um)
REAL(KIND=JPRB),  DIMENSION(nrdep)	     :: nnumb,nvolume ! number and volume distribution
REAL(KIND=JPRB),  DIMENSION(nrdep)	     :: vdi	      ! for the integration
REAL(KIND=JPRB),  DIMENSION(nrdep)	     :: vdi_def       ! for the integration
REAL(KIND=JPRB),  DIMENSION(nrdep)	     :: lure	      ! effective loo
REAL(KIND=JPRB),  DIMENSION(nrdep)	     :: ddpi	      ! integration bin-sizes
REAL(KIND=JPRB),  DIMENSION(nrdep+1)	     :: dlogdp,ddp    ! integration edges
REAL(KIND=JPRB),  DIMENSION(nrdep)	     :: logdp,logde   ! log(diameter)
    
REAL(KIND=JPRB)       :: SIGMA, LNSIGMA, DENSITY, DPG, NTOT
INTEGER(KIND=JPIM)    :: NSHIFT, IR1, IR

REAL(KIND=JPRB)    :: Z1RG, ZALPHA, ZHGT, ZZ0M,ZAERUST,ZDP(KLON,KLEV)

! Arrays to collect deposition velocities
TYPE(MODAL_DATA), DIMENSION(NMOD), TARGET :: vn_deposition
TYPE(MODAL_DATA), DIMENSION(NMOD), TARGET :: vm_deposition


REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TM5M7_DRYDEP',0,ZHOOK_HANDLE)

! Initialize..
PVDA(1:KLON,1:KLEV)=0._JPRB

do JK=1,KLEV
   ZDP(KIDIA:KFDIA,JK)=PAPH(KIDIA:KFDIA,JK)-PAPH(KIDIA:KFDIA,JK-1)
end do
!-- Compute aerodynamic resistance.
! VH - for now use fixed resistance. Particul
! VH - there is an issue with SSHF and SLHF (and in turn buoy), 
! VH - which does not seem to lead to reasonable numbers.
! VH - needs to be double-checked if these quantities are already computed before call to aerini_layer
!  DO JL=KIDIA,KFDIA
!    ZZ0M=MAX(1e-9_JPRB,PZ0M(JL))
!    ZAERUST=MAX(1e-9_JPRB,PAERUST(JL))
!    buoy = -SSHF(JL) / rhoCp &
!      &     -0.61 * PT(JL) * SLHF(JL)/rhoLv
!    tstv=-buoy/ZAERUST 
!    obuk=1E-6_JPRB
!    IF( abs(tstv) .gt. 1.E-25_JPRB .AND. abs(tstv) .LT. 1E25_JPRB ) THEN 
!       obuk = ZAERUST*ZAERUST*PT(JL)/(tstv*RG*vKarman)
!    ENDIF
!
!    IF( obuk > 0. ) THEN !   stable conditions
!       ZRA = 0.74*( alog(Href/ZZ0M) + 6.4*(Href-ZZ0M)/obuk )&
!         &   / (vKarman*ZAERUST) 
!    ELSE              !   unstable
!       y0 = sqrt(1.-9.*ZZ0M/obuk)+1.       
!       yra = sqrt(1.-9.*Href/obuk)+1. 
!       ZRA = 0.74*(alog(Href/ZZ0M)+2.*(alog(y0)-alog(yra)))/ &
!         &   (vKarman*ZAERUST) 
!    ENDIF
!    IF (ZRA .GT. 0._JPRB) THEN
!      ZRA=ZRA
!    ELSE
!!VH DEBUG
!    WRITE(NULERR,*)"DDEP DEBUG: ",ZRA, PZ0M(JL), ZZ0M, OBUK
!    WRITE(NULERR,*)"DDEB DEBUG-B:",TSTV,buoy,SSHF(JL),SLHF(JL)
!!VH END DEBUG
!
!    ENDIF 
!
!    ZRAERO(JL)=max(10._JPRB,min(ZRA,1e10_JPRB))  
!  ENDDO
! VH end. 
! VH For now take aerodynamic resistance as from 'AER' model.
    DO JL=KIDIA,KFDIA
      ZRAERO(JL)=LOG(PDZ(JL)/PZ0M(JL))/(VKARMAN*PAERUST(JL))
    ENDDO


    DO JL=KIDIA,KFDIA
            ! SEA:
            ! surface roughness  from Charnock equation
            ! friction velocity from surface stress
            !
            if(PLSM(JL) < 0.99) then 
               !VH ustar_sea(JL)=sqrt(sstr(JL))
               !VH sr_sea(JL)=alfa_charnock1*v_charnock/ustar_sea(JL) +  &
               !VH            alfa_charnock2*sstr(JL)/RG 
               !VH replace with something (hopefully) appropriate... (check?!)
               ustar_sea(JL)=PAERUST(JL)
               sr_sea(JL)=PZ0M(JL)
            else 
               ustar_sea(JL)=0.0
               sr_sea(JL)=0.0
            endif
    
	    
           freesea = MAX(0.,1.-PLSM(JL)-PCI(JL))
           if (freesea>0.01) then   !

              ! bubble bursting effect,see Hummelshoj, equation 10
              ! relationship by Wu (1988), note that Hummelshoj has not 
              ! considered the cunningham factor which yields a different 
              ! vb curve, with smaller values for small particles
              ! according to LG indeed 10 m windspeed (instead of 1 m windspeed in use in ECHAM5
              
              alpha(JL)=MIN(1.0,MAX(1.e-10,1.7e-6*PWIND(JL)**3.75))   ! set maximum!
              
              qdrop=5.*(100.*alpha(JL))   !  100 is the flux of bubbles per cm^2/s (see Monohan(1988))
              
              bubble(JL)=((100.*ustar_sea(JL))**2.)/(100.*PWIND(JL)) +  &
                   eff*(2.*rpi*rdrop**2.)*(2.*zdrop)*(qdrop/alpha(JL))
              
              !--- Correction of particle radius for particle growth close to the 
              !    surface according to Fitzgerald, 1975, the relative humidity over 
              !    the ocean is restricted to 98% (0.98) due to the salinity

               ZS=MIN(0.98,PRHCL(JL))

              !fd23012004 beta(JL)=EXP((0.00077*ZS)/(1.009-ZS)) THIS term was present in ECHAM code !
              !; max. value reached for this parameter is 1.04 and is ignored here.

              phi=1.058-((0.0155*(ZS-0.97))/(1.02-ZS**1.4))
              alpha1=1.2*EXP((0.066*ZS)/(phi-ZS))   
              vk1=10.2-23.7*ZS+14.5*ZS**2. 
              vk2=-6.7+15.5*ZS-9.2*ZS**2. 
              alpharat=1.-vk1*(1.-eps)-vk2*(1.-eps**2)
              alphae(JL)=alpharat*alpha1

              !--- Over land no correction for large humidity close to the surface:

            else ! land surface
              alpha(JL)=0.
              bubble(JL)=0.
              alphae(JL)=1.
            endif
     ENDDO	    



 ! look up table for different aerosol radii:
 do irdep = 1, nrdep

    do JL=KIDIA,KFDIA

      zr = 2.0e-4 * lur(irdep)   ! diameter in cm !
      vd_land=0.
      vd_sea=0.


      um=PWIND(JL)
      ! xland = PLSM(JL)   ![]fraction . xland replaced with PLSM(JL)

      !--- Cunningham factor: 

      cunning=1.+(cl/(alphae(JL)*zr))*(2.514+0.800*EXP(-0.55*zr/cl))

      !--- Diffusivity:

      dc=(bc* PT(JL)*cunning)/(3.*RPI*dynvisc*(alphae(JL)*zr)) ! [cm2/s] FD

      ! Relaxation:
      ! relax represents characteristic relaxation time scale [Seinfeld, p. 319]
      !      [     kg m-3 => g cm-3    ]
      densaer = density_ref   ! reference density (e.q. 1800 kg/m3)
      relax=cunning*densaer*1.E-3*((alphae(JL)*zr)**2. ) &
	 &   /(18.*dynvisc*kappa)


      ! Sedimentation is calculated operator split in the subroutine sedimentation:
      !
      ! sedspeed=(((((alphae(JL)*zr(JL)))**2.)* &
      ! 	 densaer(jl,klev,jmod,jrow)*1.E-3*grav*cunning)/(18.*dynvisc)) ! note grav should be in cm
      !    [	 kg m-3 => g cm-3    ]

      ! Calculation of schmidt 

      sc      =visc/dc  ![cm2 s-1]/[cm2 s-1] dimensionless

      if (PLSM(JL).gt.0.01) then
	 !    note that in ECHAM there is a difference between vegetaton and snow/bare soil

	 ust_land=PAERUST(JL)

	 !
	 !    calculation of stokes numbers

	 st_land  =max(0.1,(relax*(100.*ust_land)**2.)/visc)

	 !--- Calculation of the dry deposition velocity
	 !    See paper slinn and slinn, 1980, vd is related to d**2/3 
	 !    over land, whereas over sea there is accounted for slip
	 !    vb_ represents the contribution in vd of the brownian diffusion [cm s-1]
	 !    and vi_ represents the impaction [cm s-1]
	 ! 

	 vb_land   =(1./vkarman)*((ust_land/um)**2)*100.*um*(sc**(-2./3.))	 ![cm s-1]
	 vi_land   =(1./vkarman)*((ust_land/um)**2)*100.*um*(10.**(-3./st_land)) ![cm s-1]
	 vkd_land  =(vb_land+vi_land)*1e-2					 ![m s-1]
	 vd_land=1./(ZRAERO(JL)+1./vkd_land)
	 !	if (vkd_land.gt.1.) write(*,999) &
	 !	  'after vb',i,j,'jtype',jtype,'jmod',jmod,'vb',vb_land,'vi',vi_land,&
	 !	     'um10',um,'ust_land',ust_land,'dc',dc,'relax',relax,'schmidt',sc,'stokes',st_land
      endif ! PLSM.gt.0
      if (PLSM(JL) < 0.99) then
	 !--- Over sea:
	 !    Brownian diffusion for rough elements, see Hummelshoj
	 !    re is the reynolds stress:

	 ust_sea=ustar_sea(JL)
	 st_sea  =max(0.1,(relax*(100.*ust_sea)**2.)/visc)
	 re  =(100.*ust_sea*100.*sr_sea(JL))/visc     ! [cm/s]*[cm]/[cm2/s]
	 vbsea     =(1./3.)*100.*ust_sea*((sc**(-0.5))*re**(-0.5))
	 visea     =100.*ust_sea*10.**(-3./st_sea)
	 vkdaccsea =vbsea+visea
	 vkd_sea   =((1.-alpha(JL))*vkdaccsea+alpha(JL)*bubble(JL))*1e-2 ! [m s-1]
	 vd_sea=1./(ZRAERO(JL)+1./vkd_sea) ! [m s-1]
         !    if (vkd_sea.gt.1.)  write(*,999) 'sea',i,j,'jtype',jtype,'jmod',jmod,'vb',vbsea,&
         !    'vi',visea,'vkd_sea',vkd_sea,'alpha*bubble',alpha(JL)*bubble(JL),'alpha',alpha(JL)

     endif ! PLSM.lt.0.99

     ZVDA(JL,IRDEP) = min(0.1,(1.-PLSM(JL))*vd_sea + PLSM(JL)*vd_land)  ! [m s-1] limit to 10 cm/s

   enddo !JL
enddo  ! loop over nrdep






! -----------------------------------------
! Part 2: Apply aerosol deposition velocity
! -----------------------------------------





m_to_pa =  7.24e16*RG*xmair*1e3/RNAVO  !factor from m/s --> Pa/s 

!- PCEN      in unit of kg kg-1 (mass mixing ratio)
!- ZRHO      in unit of kg m-3
!- v[n|m]_deposition     in unit of m s-1    speed


DO IMODE=1,NMOD
  ALLOCATE(vn_deposition(IMODE)%surf(KIDIA:KFDIA))
  ALLOCATE(vm_deposition(IMODE)%surf(KIDIA:KFDIA))

  vn_deposition(IMODE)%surf(KIDIA:KFDIA)=0.0_JPRB
  vm_deposition(IMODE)%surf(KIDIA:KFDIA)=0.0_JPRB
ENDDO


    do mode =1,nmod
       lns(mode)      = log(sigma_lognormal(mode))
    enddo

    ! calculate the binsizes (um) around the radii of the pre-calculated vd's
    d_aer(1:NRDEP)    = 2.0*lur(1:NRDEP)      ! diameter (um)
    logdp(1:NRDEP)    = log10(d_aer(1:NRDEP)) ! log(diameter)


   do JL=KIDIA, KFDIA


          temp      = PT(JL)                                ! at surface to temp box
          pb        = PAPH(JL,KLEV+1)                            ! pressure at bottom of box (Pa)
          dp        = PAPH(JL,KLEV+1)-PAPH(JL,KLEV)              ! layer thickness
          ! to_pascal = m_to_pa*dt*pb/temp                         ! convert from m/s ---> Pa/timestep

          ! do IRDEP=1,nrdep
          !    vdi_def(IRDEP) = VDA(JL,IRDEP)
          ! enddo

          M7MODES: do mode = 1, nmod

             vt = 0.0_JPRB

             itn = mode_start(mode)   ! position of number tracer

             ! compute radius:
             radius = rw_mode(mode)%d2(JL,KLEV)

             ! initial deposition velocities for increasing radia:
             vdi(1:NRDEP) = ZVDA(JL,1:NRDEP)

             sigma   = sigma_lognormal(mode)
             lnsigma = log(sigma)
             density = dens_mode(mode)%d2(JL,KLEV)

             !if(okdebug) then 
             !   if(radius > tiny(radius)) then
             !      r_mean(mode) = r_mean(mode) + radius
             !      nr(mode)     = nr(mode) + 1
             !      r_max(mode)  = max(r_max(mode), radius)
             !   endif
             !endif

             RADENS: if (radius > 1e-11 .and. density > 1e-2) then

                ! account for density different than density_ref of the look-up table (lur --> vdi):
                lure(:) = lur(:)
                logde(:) = logdp(:)
                do ir = 2, nrdep
                   if(vdi(ir) > vdi(ir-1)) exit  ! for bigger r's :  impaction dominates (density effects)
                   if ( ir == nrdep ) exit   ! trap upper boundary
                enddo
                do ir1 = ir, nrdep
                   lure(ir1) = lur(ir1)*sqrt(density_ref/density)
                   logde(ir1) = log10(2*lure(ir1))
                enddo

                ! compress look-up table such that radii are increasing monotonic:
                nshift = 0
                ir1 = ir
                do
                   if ( logde(ir1) > logde(ir-1) ) exit
                   nshift = nshift + 1
                   ir = ir -1
                   if(ir == 1) exit
                enddo
                nrd = nrdep - nshift
                if (nshift > 0) then 
                   do ir1 = ir, nrd
                      logde(ir1) = logde(ir1+nshift)
                      lure(ir1) = lure(ir1+nshift)
                      vdi(ir1) = vdi(ir1+nshift)
                   enddo
                endif

                ! do the integration of the shifted lookup table:
                dlogdp(1) = -3.0 
                do ir=2,nrd 
                   dlogdp(ir) = 0.5*(logde(ir-1)+logde(ir))   ! take middle of the log scale
                enddo
                dlogdp(nrd+1) = 3.0   ! 1000 um
                ddp(1:nrd+1)  = 10**dlogdp(1:nrd+1)
                ddpi(1:nrd) = ddp(2:nrd+1)-ddp(1:nrd)   ! integration intervals (um)
                d_aer(1:nrd) = 2.0*lure(1:nrd)

                ! perform convolution with log-normal distribution:
                dpg = 2*radius*1e6   ! diameter in um
                ! In TM5: ntot=rm(JL,1,itn). Double-check if this is consistent?!
                ntot = PCEN(JL,KLEV,KAERO(itn))

                ! calculate the distribution (number and mass) over the deposition bins:
                if(ntot > 1.0 .and. radius > tiny(radius) ) then ! you need some aerosol!
                   do JNN=1,nrd
                      nnumb(JNN) = ntot/(sqrt(2.*RPI)*d_aer(JNN)*lnsigma)*exp(-(log(d_aer(JNN))-log(dpg))**2/(2*lnsigma**2))
                      nvolume(JNN) = nnumb(JNN)*(RPI/6.0)*d_aer(JNN)**3
                   enddo
                   vt = sum(nnumb(1:nrd)*ddpi(1:nrd)*vdi(1:nrd))/sum(nnumb(1:nrd)*ddpi(1:nrd))
                else   
                   vt = 0.0
                endif

                ! vn_deposition_mean%surf(JL,mode) = vn_deposition_mean%surf(JL, mode) + vt
                ! vn_deposition(mode)%surf(JL) = min(to_pascal*vt,ndp*dp)   ! in Pa/timestep downwards
                vn_deposition(mode)%surf(JL) = max(0._JPRB,vt)  ! keep units in (presumably ) m/sec

                !if(okdebug) then 
                !   if ( vt > tiny(vt) ) then
                !      vd_mean(mode,1) = vd_mean(mode,1) + vt
                !      vd_max(mode,1)  = max(vd_max(mode,1) , vt)
                !      nvd(mode,1)     = nvd(mode,1) + 1
                !   endif
                !endif

                ! for mass:
                if(ntot > 1.0 .and. radius > tiny(radius) ) then ! you need some aerosol!
                   vt = sum(nvolume(1:nrd)*ddpi(1:nrd)*vdi(1:nrd))/sum(nvolume(1:nrd)*ddpi(1:nrd))
                else
                   vt = 0.0
                endif
                ! vm_deposition_mean%surf(JL, mode) = vm_deposition_mean%surf(JL, mode) + vt
                ! vm_deposition(mode)%surf(JL) = min(to_pascal*vt,ndp*dp)   ! in Pa/timestep downwards
                vm_deposition(mode)%surf(JL) = max(0._JPRB, vt)   ! in m/sec


                !if(okdebug) then 
                !   if ( vt > tiny(vt) ) then
                !      vd_mean(mode,2) = vd_mean(mode,2) + vt
                !      vd_max(mode,2)  = max(vd_max(mode,2) , vt)
                !      nvd(mode,2)     = nvd(mode,2) + 1
                !   endif
                !endif   ! 

             else
                vm_deposition(mode)%surf(JL) = 0.0
                vn_deposition(mode)%surf(JL) = 0.0
             endif RADENS

          end do M7MODES

      enddo ! JL



! Compute corresponding tendencies
Z1RG = 1.0_JPRB/RG
Z1TSPHY = 1.0_JPRB/PTSPHY

! ================= 
! Loop over tracers
! ================= 
do mode =1,nmod
   !do inmode=0,mode_nm(mode)
   do inmode=0,mode_nm_sed(mode)

      JN = mode_tracers_sed(inmode,mode)
      JAER=KAERO(JN)

      !------------- reset
      NULLIFY(VS)
      if (inmode == 0) then    ! number or mass tracer
         vs => vn_deposition(mode)%surf 
      else
         vs => vm_deposition(mode)%surf
      endif
      
      !output velocity field for diagnostics purposes..
      PVDA(KIDIA:KFDIA,JAER)=VS(KIDIA:KFDIA)      

      DO JL =KIDIA,KFDIA
         ZAERO = PCEN(JL,KLEV,JAER) + PTSPHY * PTENCI(JL,KLEV,JAER)

!      using the analytical solution (Flemming et al., 2011, D_GRG_4.6)
!      The tendency in the lowest layer is modified, but the surface flux remains 
!      untouched. 

        ZHGT= PGEOH(JL,KLEV-1) * Z1RG
        ZALPHA= PTSPHY* VS(JL) /ZHGT
!-- using Euler forward
!        ZAERI = ZAERO * (1.0_JPRB - ZALPHA)
!-- using Euler backward
!        ZAERI = ZAERO * (1.0_JPRB + ZALPHA)
!-- using Euler centered
!        ZAERI = ZAERO * ((1.0_JPRB - ZALPHA)/(1.0_JPRB + ZALPHA))
!-- using the analytical solution (Flemming et al., 2011, D_GRG_4.6)
        ZAERI = ZAERO * EXP(-1.0_JPRB * ZALPHA)
        PTENC(JL,KLEV,JAER)= PTENCI(JL,KLEV,JAER) &
         &                  + (ZAERI-ZAERO) * Z1TSPHY
        !vh check indices.. PFLUXAER(JL,JAER)= (ZAERO - ZAERI)*Z1TSPHY * PDP(JL,KLEV) * Z1RG
        ! No update of surface flux:
         PFLUXAER(JL,JAER)= (ZAERO - ZAERI)*Z1TSPHY * ZDP(JL,KLEV) * Z1RG!PFLUXAER(JL,JAER)
      ENDDO
      
   enddo ! loop over tracers in mode
enddo ! loop over modes



DO IMODE=1,NMOD
  IF(ASSOCIATED(vn_deposition(IMODE)%surf)) DEALLOCATE(vn_deposition(IMODE)%surf)
  IF(ASSOCIATED(vm_deposition(IMODE)%surf)) DEALLOCATE(vm_deposition(IMODE)%surf)
ENDDO



!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TM5M7_DRYDEP',1,ZHOOK_HANDLE)

  

END SUBROUTINE TM5M7_DRYDEP
