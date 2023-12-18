SUBROUTINE m7_simple_sulfur_DRYDEP &
!---inputs
 & ( YDMODEL, KIDIA  , KFDIA , KLON  , KLEV  , &
 &   PAERO  , PFAERI , &
 &   PDP    , PGEOH, PRHO  , PTAERI, PTSPHY, &
 &   PSO2DD, PGELAM,  &
!-- outputs
 &   PFAERO , PTAERO, PFDRYD)


!**** *m7_simple_sulfur_DRYDEP* -  ROUTINE FOR PARAMETRIZATION OF DRY DEPOSITION
  
  !      Modified from Dry dep routine for AER aerosols (aer_drydep.f90) for use in 
  !      OpenIFS/AC dry deposition of so2 when using simple sulfur chemistry
  !      Original by
  !      Jean-Jacques Morcrette 
  !      following O.Boucher's formulation for LMD-Z
  
  !      Modifications Tommi Bergman (FMI)
  !       
  !      Dry deposition is (simply) represented by a modification of the 
  !      instantaneous surface flux by what comes down from layer just 
  !      above the surface
  
!**   INTERFACE.
!     ----------
!          *M7_simple_sulfur_DRYDEP* IS CALLED FROM *TM5M7*.

! INPUTS:
! -------
! PFAERI(KLON,NACTAERO)      : INPUT SURFACE FLUX        (xx m-2)
! PTAERI(KLON,KLEV,NACTAERO) : INPUT TENDENCIES          (xx kg s-1)

! OUTPUTS:
! --------
! PFAERO(KLON,NACTAERO)      : SURFACE FLUX              (xx m-2)
! PTAERO(KLON,KLEV,NACTAERO) : UPDATED TENDENCIES        (xx kg s-1)
! PFDRYD(KLON)             : DIAGNOSTIC DRY DEPOSITION AT THE SURFACE (xx m-2)

!     EXTERNALS.
!     ----------
!          NONE

!     MODIFICATIONS.
!     -------------
!          JJMorcrette 20110725 maximum deposition speed
!          SRémy       20160309 SO2 dry deposition velocity from SUMO (same as
!          CHEM)
!          SRémy       20160830 dry deposition velocities (except SO2) computed
!          following Zhang et al 2001.
!          TBergman    20230207 remove most of the code which is not needed for so2

!     SWITCHES.
!     --------

!     MODEL PARAMETERS
!     ----------------

!-----------------------------------------------------------------------
USE TYPE_MODEL         , ONLY : MODEL
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

!USE YOEAERATM, ONLY : YREAERATM
!USE YOEAERSNK, ONLY : YREAERSNK
!USE YOEAERSRC ,ONLY : YREAERSRC
!USE YOM_YGFL , ONLY : YGFL
USE YOMCST   , ONLY : RG, RPI
!USE YOMRIP  , ONLY : YRRIP
!USE YOMLUN   , ONLY : NULOUT
USE MO_TRACDEF, ONLY: ntrac, trlist
USE MO_HAM, ONLY: nclass, naerocomp, sizeclass, nccndiag, subm_ngasspec

IMPLICIT NONE

!-----------------------------------------------------------------------

!*       0.1  ARGUMENTS
!             ---------

!---input fields
TYPE(MODEL)       ,INTENT(IN)    :: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
!INTEGER(KIND=JPIM),INTENT(IN)    :: KTDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
!INTEGER(KIND=JPIM),INTENT(IN)    :: KSTEP 

!REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSM(KLON), PCI(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAERO(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NACTAERO), PTAERI(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NACTAERO)
!REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,0:KLEV)
!REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDP(KLON,KLEV), PGEOH(KLON,0:KLEV), PRHO(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFAERI(KLON,YDMODEL%YRML_GCONF%YGFL%NACTAERO)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSO2DD(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGELAM(KLON)
!REAL(KIND=JPRB)   ,INTENT(IN)    :: PAERUST(KLON)
!REAL(KIND=JPRB)   ,INTENT(IN)    :: PZ0M(KLON)
!REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV)
!REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ(KLON)
!REAL(KIND=JPRB)   ,INTENT(IN)    :: PDZ(KLON)

!---output fields
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFAERO(KLON,YDMODEL%YRML_GCONF%YGFL%NACTAERO)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTAERO(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NACTAERO)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFDRYD(KLON)!,YGFL%NACTAERO)

!--- local variables
INTEGER(KIND=JPIM) :: JL, JAER, JAERO
REAL(KIND=JPRB) :: ZAERO, ZFROC, ZOCEA, ZLAND, ZLICE, ZSICE, ZRHO
REAL(KIND=JPRB) :: Z1RG , Z1TSPHY, ZHGT, ZALPHA, ZAERI
REAL(KIND=JPRB) :: ZVDEP,ZVDEP2(YDMODEL%YRML_GCONF%YGFL%NACTAERO)

REAL(KIND=JPRB) :: ZRHOP
REAL(KIND=JPRB) :: ZWETD
REAL(KIND=JPRB) :: ZSIGMA

!!$!* Taken from J.J.Morcrette
!!$REAL(KIND=JPRB), PARAMETER    :: ZR_OM = 0.13E-6 !m
!!$REAL(KIND=JPRB), PARAMETER    :: ZR_BC = 0.04E-6 !m
!!$REAL(KIND=JPRB), PARAMETER    :: ZR_SO4 = 0.9E-6   !m SO4 dry particle radius,Martin et al., 2003
!!$REAL(KIND=JPRB), PARAMETER    :: ZR_AM = 0.35E-6   !m Wang et al ACP 2014
!!$REAL(KIND=JPRB), DIMENSION(2),PARAMETER    :: ZR_NI =(/0.35E-6,1.5E-6/)   !m Wang et al ACP 2014
!!$!* Growth factors corresponding to RH table as given in RRHTAB, according to
!!$!J.J. Morcrette
!!$REAL(KIND=JPRB), DIMENSION(12), PARAMETER    :: ZRH_GROWTH_SO4= &
!!$         & (/1.00,1.00,1.00,1.00,1.169,1.220,1.282,1.363,1.485,1.581,1.732,2.085/)
!!$! According to Chin et al., AMS 2002
!!$REAL(KIND=JPRB), DIMENSION(12), PARAMETER    :: ZRH_GROWTH_BC= &
!!$         & (/1.00,1.00,1.00,1.00,1.00,1.000,1.000,1.000,1.200,1.300,1.400,1.500/)
!!$REAL(KIND=JPRB), DIMENSION(12), PARAMETER    :: ZRH_GROWTH_OM= &
!!$         & (/1.00,1.00,1.00,1.00,1.169,1.200,1.300,1.400,1.500,1.550,1.600,1.800/)
!!$! Nitrate : Svenningsson et al ACP 2006
!!$REAL(KIND=JPRB), DIMENSION(12), PARAMETER    :: ZRH_GROWTH_NI= &
!!$         & (/1.00,1.00,1.00,1.00,1.100,1.200,1.250,1.300,1.350,1.500,1.700,2.100/)
!!$REAL(KIND=JPRB), DIMENSION(12), PARAMETER    :: ZRH_GROWTH_AM= &
!!$         & (/1.00,1.00,1.00,1.00,1.169,1.220,1.282,1.363,1.485,1.581,1.732,2.085/)

!!$REAL(KIND=JPRB), PARAMETER    :: ZRHO_OM=1800  ! kg/m^3
!!$REAL(KIND=JPRB), PARAMETER    :: ZRHO_BC=1000  ! kg/m^3
!!$REAL(KIND=JPRB), PARAMETER    :: ZRHO_SO4=1760  ! kg/m^3
!!$REAL(KIND=JPRB), PARAMETER    :: ZRHO_H2O=1000  ! kg/m^3 (water)
!!$REAL(KIND=JPRB), PARAMETER    :: ZRHO_AM=1760  ! kg/m^3 (water)
!!$REAL(KIND=JPRB), DIMENSION(2), PARAMETER    :: ZRHO_NI=(/1730,1400/)  ! kg/m^3 (water)



REAL(KIND=JPRB) :: ZMAXVDRY(KLON), ZDZ(KLON),ZDVMAX,ZHOURLT, ZSCALE, ZVFRAC
REAL(KIND=JPRB) :: ZQSAT(KLON,KLEV)
REAL(KIND=JPRB) :: ZRHCL(KLON)
INTEGER(KIND=JPIM) :: IAER, ICAER, ITAER, ITYP, IBIN, JTAB
INTEGER(KIND=JPIM) :: INAER(18),IRH(KLON)



LOGICAL :: LLPRINT


REAL(KIND=JPRB) :: ZHOOK_HANDLE

!#include "aer_drydepvel.intfb.h"!
!#include "satur.intfb.h"!

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('M7_SIMPLE_SULFUR_DRYDEP',0,ZHOOK_HANDLE) 
ASSOCIATE(YGFL=>YDMODEL%YRML_GCONF%YGFL, YREAERSNK=>YDMODEL%YRML_PHY_AER%YREAERSNK, YDRIP=>YDMODEL%YRML_GCONF%YRRIP)
ASSOCIATE(NDRYDEP=>YREAERSNK%NDRYDEP, NACTAERO=>YGFL%NACTAERO)


! & RVDPLIC=>YREAERSNK%RVDPLIC, RVDPLND=>YREAERSNK%RVDPLND, &
! & RVDPOCE=>YREAERSNK%RVDPOCE, RVDPSIC=>YREAERSNK%RVDPSIC, &
! & RRHO_DD=>YREAERSNK%RRHO_DD, RRHO_SS=>YREAERSNK%RRHO_SS, &
! & RMMD_DD=>YREAERSNK%RMMD_DD, RMMD_SS=>YREAERSNK%RMMD_SS, &
! & RSSDENS_RHTAB=>YREAERSNK%RSSDENS_RHTAB,                 &
! & RSSGROWTH_RHTAB=>YREAERSNK%RSSGROWTH_RHTAB,             &
! & RRHTAB=>YREAERSNK%RRHTAB,                               &
! & NACTAERO=>YGFL%NACTAERO, &
! & NINDDDEP=>YREAERATM%NINDDDEP, &
! & LDRYDEPVEL_DYN => YREAERSNK%LDRYDEPVEL_DYN, &
! & NTDDEP=>YREAERATM%NTDDEP,YAERO=>YGFL%YAERO, &
! & NMAXTAER=>YREAERSRC%NMAXTAER, NTYPAER=>YREAERSRC%NTYPAER)

!- PFAERI   in unit of xx m-2 s-1 (surface flux)
!- PFAERO    in unit of xx m-2 s-1 (surface flux)
!- N.B. Surface emission fluxes are negative upward, so contribution of 
!  dry deposition is to make the surface flux less negative 

LLPRINT=.FALSE.
!DO JL=1,NSTPDBG
!  IF (KSTEP == KSTPDBG(JL)) THEN
!    LLPRINT=.TRUE.
!  ENDIF
!ENDDO

!!$ICAER=0
!!$DO JAER=1,NMAXTAER
!!$  IF (NTYPAER(JAER) /= 0) THEN
!!$    ITAER=NTYPAER(JAER)
!!$    DO IAER=1,ITAER
!!$      ICAER=ICAER+1
!!$      INAER(ICAER)=JAER*10+IAER
!!$    ENDDO
!!$  ENDIF
!!$ENDDO


Z1RG   = 1.0_JPRB/RG
Z1TSPHY= 1.0_JPRB/PTSPHY
PTAERO(KIDIA:KFDIA,1:KLEV,1:YGFL%NACTAERO)=PTAERI(KIDIA:KFDIA,1:KLEV,1:YGFL%NACTAERO)
PFAERO(KIDIA:KFDIA,1:YGFL%NACTAERO)= PFAERI(KIDIA:KFDIA,1:YGFL%NACTAERO)
ZVDEP=0._JPRB
!!$ZVDEP2(:)=0._JPRB

!!$DO JL=KIDIA,KFDIA
!!$  ZDZ(JL)=(PAPH(JL,KLEV)-PAPH(JL,KLEV-1)) / (RG*PRHO(JL,KLEV))  
!!$  ZMAXVDRY(JL)= ZDZ(JL) / PTSPHY
!!$ENDDO

! limit max der dep velocity based on 30 m box height (it is 10-15 m)
ZDVMAX=30.0_JPRB/PTSPHY

!CALL SATUR (KIDIA, KFDIA, KLON  , KTDIA, KLEV, PAP, PT, ZQSAT, 2 )
PFDRYD(:)=0.0_JPRB
!IRH=12
!write(3434,*)NDRYDEP
!NDRYDEP=1
DO JAERO=1,subm_ngasspec!NACTAERO!NTDDEP
  !JAER=NINDDDEP(JDDEP)
   !write(3535,*)JAERO,TRIM(trlist%ti(JAERO)%basename)
   !IF (TRIM(YAERO(JAERO)%CNAME)=='SO2')THEN
   IF(TRIM(trlist%ti(JAERO)%basename)=='SO2')THEN
      JAER=JAERO
  ELSE
     CYCLE
  END IF
!  ITYP=INAER(JAER)/10
!  IBIN=INAER(JAER)-ITYP*10
  DO JL=KIDIA,KFDIA

!!$    ZRHCL(JL)=PQ(JL)/ZQSAT(JL,KLEV)
    ZRHO = PRHO(JL,KLEV)
!!$    DO JTAB=1,12
!!$      IF (ZRHCL(JL)*100._JPRB > RRHTAB(JTAB)) THEN
!!$        IRH(JL)=JTAB
!!$      ENDIF
!!$    ENDDO
    ZAERO = PAERO(JL,KLEV,JAER) + PTSPHY * PTAERI(JL,KLEV,JAER)
    !write(3434,*)PAERO(JL,KLEV,JAER) , PTSPHY * PTAERI(JL,KLEV,JAER)
    !IF (ITYP == 5 .AND. IBIN == 2) THEN
!     IF (.true.) THEN
       ! use SUMO dry dep velocity for SO2
       !ZVDEP(JAER)=PSO2DD(JL)
       ZVDEP=PSO2DD(JL)
       ZHOURLT  =  (RPI - YDRIP%RWSOVR) - PGELAM(JL)
       ! Difference w.r.t to longitude of sza max
       ZSCALE = 1.0 + COS(ZHOURLT) * 0.7_JPRB
       !ZVDEP(JAER) = ZSCALE * MIN(ZDVMAX,ZSCALE*ZVDEP(JAER))
       ZVDEP = ZSCALE * MIN(ZDVMAX,ZSCALE*ZVDEP)
!!$    ELSE
!!$      ZLAND= PLSM(JL)
!!$      ZOCEA= 1._JPRB-PLSM(JL)
!!$      ZFROC= ZOCEA*(1._JPRB-PCI(JL))
!!$      ZSICE= ZOCEA*PCI(JL)
!!$      ZLICE= 0._JPRB
!!$
!!$      ZVDEP2(JAER)= ZLAND * RVDPLND(JAER) + ZFROC * RVDPOCE(JAER) + &
!!$        &    ZSICE * RVDPSIC(JAER) + ZLICE * RVDPLIC(JAER)
!!$      ZVDEP2(JAER)= MIN(ZVDEP2(JAER),ZMAXVDRY(JL))
!!$
!!$      ! set RHOP and WETD for each aerosol type
!!$      ZSIGMA=2.0_JPRB
!!$      IF (ITYP == 1) THEN
          ! for land, adjust size
!!$          IF (PLSM(JL) > 0.5_JPRB) THEN
!!$            ZRHOP=RSSDENS_RHTAB(12)
!!$            ZWETD=1.E-6_JPRB*RMMD_SS(IBIN)*RSSGROWTH_RHTAB(12)/RSSGROWTH_RHTAB(9)
!!$          ELSE
!!$            ZRHOP=RSSDENS_RHTAB(IRH(JL))
!!$            ZWETD=1.E-6_JPRB*RMMD_SS(IBIN)*RSSGROWTH_RHTAB(IRH(JL))/RSSGROWTH_RHTAB(9)
!!$          ENDIF
!!$      ELSEIF (ITYP == 2) THEN
!!$          ZRHOP=RRHO_DD(IBIN)
!!$          ZWETD=1.E-6_JPRB*RMMD_DD(IBIN)
!!$        ELSEIF (ITYP == 3) THEN
!!$          ZWETD=2._JPRB*ZR_OM*ZRH_GROWTH_OM(IRH(JL))
!!$          ZVFRAC = 1.0 / ZRH_GROWTH_OM(IRH(JL))**3
!!$          ZRHOP = ZRHO_H2O*(1.0-ZVFRAC) + ZVFRAC*ZRHO_OM
!!$        ELSEIF (ITYP == 4) THEN
!!$          ZVFRAC = 1.0 / ZRH_GROWTH_BC(IRH(JL))**3
!!$          ZRHOP = ZRHO_H2O*(1.0-ZVFRAC) + ZVFRAC*ZRHO_BC
!!$          ZWETD=2._JPRB*ZR_BC*ZRH_GROWTH_BC(IRH(JL))
!!$        ELSEIF (ITYP == 5) THEN ! SO4 only
!!$          ZVFRAC = 1.0 / ZRH_GROWTH_SO4(IRH(JL))**3
!!$          ZRHOP = ZRHO_H2O*(1.0-ZVFRAC) + ZVFRAC*ZRHO_SO4
!!$          ZWETD=2._JPRB*ZR_SO4*ZRH_GROWTH_SO4(IRH(JL))
!!$        ELSEIF (ITYP == 6) THEN ! nitrate
!!$          ZVFRAC = 1.0 / ZRH_GROWTH_NI(IRH(JL))**3
!!$          ZRHOP = ZRHO_H2O*(1.0-ZVFRAC) + ZVFRAC*ZRHO_NI(IBIN)
!!$          ZWETD=2._JPRB*ZR_NI(IBIN)*ZRH_GROWTH_NI(IRH(JL))
!!$        ELSEIF (ITYP == 7) THEN ! Ammonium
!!$          ZVFRAC = 1.0 / ZRH_GROWTH_AM(IRH(JL))**3
!!$          ZRHOP = ZRHO_H2O*(1.0-ZVFRAC) + ZVFRAC*ZRHO_AM
!!$          ZWETD=2._JPRB*ZR_AM*ZRH_GROWTH_AM(IRH(JL))
!!$        ENDIF
!!$
!!$        IF (LDRYDEPVEL_DYN) THEN
!!$          ! compute deposition velocity following Zhang et al 2001
!!$          ! CALL AER_DRYDEPVEL(ZRHOP,ZWETD,ZSIGMA,PZ0M(JL),PCI(JL),PAERUST(JL),PDZ(JL),PT(JL,KLEV),ZRHO,ZVDEP(JAER))
!!$          CALL ABOR1("dynamic dry dep for AER not yet supported in OIFS!")
!!$
!!$           ZVDEP=MIN(0.1_JPRB,ZVDEP)
!!$        ELSE
!!$           ZVDEP(JAER)=ZVDEP2(JAER)
!!$        ENDIF

! RVDPxxx (from su_aerp, derived from LMDZ in m s-1)
!   formula is consistent: (xx m-2 s-1) - (m s-1) * (xx kg-1) * (kg m-3)
!   PFAERO is then in xx m-2 s-1
!    ENDIF

    IF (NDRYDEP == 1) THEN

!-- only the surface flux is diminished of the equivalent effect of the dry deposition
!   but the vertical distribution of tendencies is untouched 

      PFAERO(JL,JAER)=PFAERI(JL,JAER) &
        & + ZVDEP * ZAERO * ZRHO
      !PFDRYD(JL,JAER)=0._JPRB
      PFDRYD(JL)=0._JPRB

    ELSE

!-- Alternate formulation
!    using the analytical solution (Flemming et al., 2011, D_GRG_4.6)
!   The tendency in the lowest layer is modified, but the surface flux remains 
!    untouched. 

      ZHGT= PGEOH(JL,KLEV-1) * Z1RG
      ZALPHA= PTSPHY* ZVDEP/ZHGT
!-- using Euler forward
!      ZAERI = ZAERO * (1.0_JPRB - ZALPHA)
!-- using Euler backward
!      ZAERI = ZAERO * (1.0_JPRB + ZALPHA)
!-- using Euler centered
!      ZAERI = ZAERO * ((1.0_JPRB - ZALPHA)/(1.0_JPRB + ZALPHA))
!-- using the analytical solution (Flemming et al., 2011, D_GRG_4.6)
      ZAERI = ZAERO * EXP(-1.0_JPRB * ZALPHA)
!      PTAERO(JL,KLEV,JAER)= PTAERI(JL,KLEV,JAER) &
!       &                  + (ZAERI-PAERO(JL,KLEV,JAER)) * Z1TSPHY
      PTAERO(JL,KLEV,JAER)= PTAERO(JL,KLEV,JAER) &
       &                  + (ZAERI-ZAERO) * Z1TSPHY
      !write(3434,*)ZAERO , ZAERI
      !write(3434,*)Z1TSPHY
      !write(3434,*)PDP(JL,KLEV)
      !write(3434,*)Z1RG
      !if (((ZAERO - ZAERI))>1e-20_JPRB)THEN
         PFDRYD(JL)= (ZAERO - ZAERI)*Z1TSPHY * PDP(JL,KLEV) * Z1RG
      !ELSE
      !   PFDRYD(JL)= 0.0_JPRB
      !END if
         PFAERO(JL,JAER)= PFAERI(JL,JAER)

   ENDIF
!   write(3434,*)NDRYDEP,JL,PFAERO(JL,JAER),PFDRYD(JL)
  ENDDO
ENDDO

!IF (LLPRINT) THEN
!  WRITE(UNIT=NULOUT,FMT='(1x,''DRYDEP'',I5,9E12.5)') KSTEP,PRHO(KIDIA,KLEV),&
!    & (PFAERI(KIDIA,JAER),PFAERO(KIDIA,JAER),ZVDEP(JAER),PAERO(KIDIA,KLEV,JAER),JAER=3,YGFL%NACTAERO,3)
!ENDIF

!-----------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('M7_SIMPLE_SULFUR_DRYDEP',1,ZHOOK_HANDLE)
END SUBROUTINE M7_SIMPLE_SULFUR_DRYDEP
