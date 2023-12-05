SUBROUTINE TM5M7_SEDIMENT &
 !---input
 & ( KIDIA , KFDIA   , KLON, KLEV, KTRAC, KAERO   , &
 &   PTSPHY, PT      , PAP , PAPH, &
 !---prognostic fields
 &   RW_MODE, DENS_MODE, PCEN, &
 !---output
 &   PFLUXAER, PTENC)  

!**** *TM5M7_SEDIMENT* -  ROUTINE FOR PARAMETRIZATION OF TM5M7 AEROSOL SEDIMENTATION
!
! !DESCRIPTION:  This module calculates sedimentation of aerosol particles. 
!                Also the deposition at the surface is calculated here since
!                 it uses similar routines. 
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
!     ln2s    = (log(sigma_lognormal(mode)))**2
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
!          *TM5M7_SEDIMENT* IS CALLED FROM *TM5M7_PHY2*.

! INPUTS:
! -------
! PTSPHY                : TIMESTEP                  (s)
! PAP     (KLON,KLEV)   : LEVEL PRESSURE            (Pa)
! PAPH    (KLON,KLEV+1) : HALF-LEVEL PRESSURE       (Pa)
! PCI     (KLON)        : FRACTION OF SEA ICE
! PLSM    (KLON)        : LAND-SEA MASK   
! PT      (KLON,KLEV)   : LEVEL TEMPERATURE         (K)
! PCEN    (KLON,KLEV,KTRAC)   : CONCENTRATION OF TRACERS  (xx kg-1)

! OUTPUTS:
! --------
! PTENC  (KLON,KLEV,KTRAC)    : TENDENCY                  (xx kg-1 s-1)

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
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOEAERSRC ,ONLY : YREAERSRC
USE YOEAERSNK ,ONLY : YREAERSNK
USE YOMCST    ,ONLY : RG, RD,RNAVO
USE YOM_YGFL  ,ONLY : YGFL 
USE TM5M7_DATA, ONLY : NMOD, MODE_NM_SED, mode_tracers_sed, XMAIR,MODAL_DATA, &
    & SIGMA_LOGNORMAL

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
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV) ! Temp
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV)! Mid-lev pres 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1) ! interface pres
!REAL(KIND=JPRB)  ,INTENT(IN)    :: PAERO(KLON,KLEV) 
!REAL(KIND=JPRB)  ,INTENT(IN)    :: PTAERI(KLON,KLEV) 
TYPE(MODAL_DATA)  ,INTENT(IN)    :: RW_MODE(NMOD)
TYPE(MODAL_DATA)  ,INTENT(IN)    :: DENS_MODE(NMOD)

REAL(KIND=JPRB)   ,INTENT(IN)    :: PCEN(KLON,KLEV,KTRAC)
!---output fields
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFLUXAER(KLON,YGFL%NACTAERO)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENC(KLON,KLEV,KTRAC)


!---local fields
REAL(KIND=JPRB)     :: TEMP, PB, DP, ZVIS, RHO_AIR, TO_PASCAL,RADIUS,DPG,VT
REAL(KIND=JPRB)     :: ZSEDFLX(KLON), ZAERONWM1(KLON), ZAERI(KLON,KLEV), &
       &               ZSOLAERS, ZSOLAERB, ZGDP, ZDTGDP, ZKK, & 
       &               SIGMA, LNSIGMA,DENSITY,SLINNFAC, &
       &               ZRHO, ZAERONW
INTEGER(KIND=JPIM)  :: JK, JL, JN, MODE, IMODE, INMODE
INTEGER(KIND=JPIM), PARAMETER           :: ndp = 3 ! limit of the velocity to ndp times the layer thickness
                                                   ! iteration will account for a fall length through multiple layers

REAL(KIND=JPRB)      :: m_to_pa  !factor from m/s --> Pa/s 

REAL(KIND=JPRB),DIMENSION(:,:),POINTER   :: VS  

! Arrays to collect sedimentation velocities
TYPE(MODAL_DATA), DIMENSION(NMOD), TARGET :: VN_SEDIMENTATION
TYPE(MODAL_DATA), DIMENSION(NMOD), TARGET :: VM_SEDIMENTATION


REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TM5M7_SEDIMENT',0,ZHOOK_HANDLE)


m_to_pa =  7.24e16*RG*xmair*1e3/RNAVO  !factor from m/s --> Pa/s 

!- PAERO     in unit of xx kg-1 (mixing ratio)
!- ZAERONW   in unit of xx kg-1
!- ZSOLAERS, ZSOLAERB in xx kg-1  as m s-2/DeltaP xx m-2 = m s-2/(kg m-1 s-2) xx m-2
!- ZSEDFLX   in unit of xx m-2
!- ZRHO      in unit of kg m-3
!- v[n|m]_sedimentation     in unit of m s-1    speed


DO IMODE=1,NMOD
  ALLOCATE(VN_SEDIMENTATION(IMODE)%d2(KIDIA:KFDIA,KLEV))
  ALLOCATE(VM_SEDIMENTATION(IMODE)%d2(KIDIA:KFDIA,KLEV))

  VN_SEDIMENTATION(IMODE)%d2(KIDIA:KFDIA,1:KLEV)=0.0_JPRB
  VM_SEDIMENTATION(IMODE)%d2(KIDIA:KFDIA,1:KLEV)=0.0_JPRB
ENDDO

! Loop over levels
do JK=2, KLEV


   do JL=KIDIA, KFDIA

         if(JK == KLEV) then 
            temp = PT(JL,KLEV)                   ! at surface to temp box
         else
            temp = 0.5*(PT(JL,JK)+PT(JL,JK+1))   ! interpolate to bottom of box
         endif

         pb           = PAPH(JL,JK)                    ! pressure at bottom of box (Pa)
         dp           = PAPH(JL,JK)-PAPH(JL,JK-1)      ! layer thickness
         zvis           = 0.0000014963*temp*sqrt(temp)/(temp+120.)        ! viscosity  (kg/ms)
         rho_air   = 7.24e16*pb/temp  * xmair*1e3/RNAVO         ! kg/m3
         to_pascal = m_to_pa*pb/temp                         ! convert from m/s ---> Pa/sec

         M7MODES: do mode = 1, nmod

            radius = rw_mode(mode)%d2(JL,JK)

            sigma = sigma_lognormal(mode)
            lnsigma = log(sigma)
            density = dens_mode(mode)%d2(JL,JK)

            slinnfac = exp(2.0*lnsigma*lnsigma)


            ! for number: take first mode (Seinfeld, 1986, pg 286)
            dpg = radius*2.0*exp(lnsigma*lnsigma)   ! diameter in m  
            vt =  TM5M7_fall1(pb,Dpg,zvis,temp,rho_air,density)


            !VH vn_sedimentation(mode)%d2(JL,JK) = min(to_pascal*vt,ndp*dp)   ! in Pa/sec downwards
            vn_sedimentation(mode)%d2(JL,JK) = vt     ! in m/sec downwards

            !VH if(JK == KLEV) then
            !VH   vn_sedimentation_mean%surf(JL,mode) = &
            !VH        vn_sedimentation_mean%surf(JL,mode) + vt
            !VH endif

            ! for mass: the third moment
            dpg = radius*2.0*exp(3*lnsigma*lnsigma)  ! diameter in m
            vt =  TM5M7_fall1(pb,Dpg,zvis,temp,rho_air,density)


            !VH vm_sedimentation(mode)%d2(JL,JK) = min(to_pascal*vt*slinnfac,ndp*dp)  ! multiply with slinnfac
            vm_sedimentation(mode)%d2(JL,JK) = vt*slinnfac  ! multiply with slinnfac

            !VH if(JK == KLEV) then
            !VH   vm_sedimentation_mean%surf(JL,mode) = &
            !VH        vm_sedimentation_mean%surf(JL,mode) + vt*slinnfac
            !VH endif

         enddo M7MODES
      enddo ! JL
enddo ! LEVS



! Compute corresponding tendencies
! ================= 
! Loop over tracers
! ================= 
do mode =1,nmod
   !do inmode=0,mode_nm(mode)
   do inmode=0,mode_nm_sed(mode)

      JN = mode_tracers_sed(inmode,mode)

      !------------- reset
      NULLIFY(VS)
      if (inmode == 0) then    ! number or mass tracer
         vs => vn_sedimentation(mode)%d2 
      else
         vs => vm_sedimentation(mode)%d2 
      endif


      !--initialisations of variables carried out from one layer to the next layer
      !--actually not needed if (JK>1) test is on
        ZSEDFLX(KIDIA:KFDIA)=0.0_JPRB
        ZAERONWM1(KIDIA:KFDIA)=0.0_JPRB
      
      DO JK=1, KLEV
      
        DO JL=KIDIA, KFDIA
      !--initialisations
          ZSOLAERS=0.0
          ZSOLAERB=0.0
          ZGDP=RG/(PAPH(JL,JK+1)-PAPH(JL,JK))
          ZDTGDP=PTSPHY*ZGDP
      
      !- Starting point is input concentration field
      
         ZAERI(JL,JK) = PCEN(JL,JK,KAERO(JN))
      
      ! (Or: apply intermediate tendencies??)
      !    ZAERI(JL,JK) = PCEN(JL,JK,KAERO(JN)) + PTSPHY * PTAERI(JL,JK,KAERO(JN)) 
      !  with PTAERI containing the tendencies from dry deposition
      
      ! source from above
          IF (JK>1) THEN 
            ZSEDFLX(JL)=ZSEDFLX(JL)*ZAERONWM1(JL)  
            ZSOLAERS=ZSOLAERS+ZSEDFLX(JL)*ZDTGDP
          ENDIF
      
      ! sink to next layer
          ZRHO=PAP(JL,JK)/(RD*PT(JL,JK))
          ZSEDFLX(JL)=VS(JL,JK)*ZRHO
          ZSOLAERB=ZSOLAERB+ZDTGDP*ZSEDFLX(JL)
      
      !---implicit solver
          ZAERONW=(ZAERI(JL,JK)+ZSOLAERS)/(1.0_JPRB+ZSOLAERB)
      
      !---new time-step AER variable needed for next layer
          ZAERONWM1(JL)=ZAERONW
      
      !---tendency in unit of xx kg-1 s-1
          PTENC(JL,JK,KAERO(JN))=PTENC(JL,JK,KAERO(JN))+(ZAERONW-ZAERI(JL,JK))/PTSPHY
      !
        ENDDO
      ENDDO
      
      !---sedimentation flux to the surface
      !---ZAERONWM1 now contains the surface concentration at the new timestep
      !---PFLUXAER in unit of xx m-2 s-1 
      DO JL=KIDIA,KFDIA 
        ZRHO=PAP(JL,KLEV)/(RD*PT(JL,KLEV))
        PFLUXAER(JL,JN)=ZRHO*ZAERONWM1(JL)*VS(JL,KLEV)
      ENDDO
      
   enddo ! loop over tracers in mode
enddo ! loop over modes



DO IMODE=1,NMOD
  DEALLOCATE(VN_SEDIMENTATION(IMODE)%d2)
  DEALLOCATE(VM_SEDIMENTATION(IMODE)%d2)
ENDDO



!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TM5M7_SEDIMENT',1,ZHOOK_HANDLE)
CONTAINS



  ! !IROUTINE:    TM5M7_FALL1
  !
  ! !DESCRIPTION: function to calculate the fall velocity from the particle 
  !               diameter, as a function of density, temperature and pressure.
  !\\
  !\\
  ! !INTERFACE:
  !
  REAL(KIND=JPRB) FUNCTION TM5M7_FALL1( press, zmd, zvis, t, zdensair, densaer_p) result(vt)
    !
    ! !INPUT PARAMETERS:
    !
    REAL(KIND=JPRB), INTENT(IN)      :: press      ! pressure in Pa (or bar?)
    REAL(KIND=JPRB), INTENT(IN)      :: zmd        ! median radius of aerosol (m)
    REAL(KIND=JPRB), INTENT(IN)      :: zvis       ! viscosity  (kg/(sm))
    REAL(KIND=JPRB), INTENT(IN)      :: t          ! temperature (K)
    REAL(KIND=JPRB), INTENT(IN)      :: zdensair   ! density air (kg/m3)
    REAL(KIND=JPRB), INTENT(IN)      :: densaer_p  ! density aerosol particles (kg/m3)
    !
    ! !REVISION HISTORY: 
    !    2 Apr 2010 - P. Le Sager - 
    !
    ! !REMARKS:
    ! Called in Sedimentation_Apply, only if m7 used.
    !
    !EOP
    !------------------------------------------------------------------------------
    !BOC

    REAL(KIND=JPRB) :: zlair
    
    ! ----------- start
    if (zmd > tiny(zmd)) then 
       vt=2./9.*(densaer_p-zdensair) * RG/zvis*(zmd/2.)**2.

       zlair=0.066*(1.01325e5/press)*t/293.15*1e-6

       !--- With cunnigham slip- flow correction (S&P, Equation 8.34):

       vt = vt * (1.+ 1.257*zlair/zmd*2. +  0.4*zlair/zmd*2. *   exp(-1.1/(zlair/zmd*2.)) )
    else
       vt = 0.0
    endif

  END FUNCTION TM5M7_FALL1
  

END SUBROUTINE TM5M7_SEDIMENT
