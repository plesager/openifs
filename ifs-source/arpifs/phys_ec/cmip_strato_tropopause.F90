! (C) Copyright 2026- FMI.
! (C) Copyright 2026- KNMI.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
SUBROUTINE CMIP_STRATO_TROPOPAUSE &
& (KIDIA, KFDIA, KLON, KLEV, PAPF, PAPRS, PT, PTH, KTROP)

!     PURPOSE.
!     --------
!     Tropopause level used together with the prescribed CMIP
!     stratospheric aerosol forcing.
!
!     In the radiation scheme, the prescribed stratospheric aerosols
!     are applied at the levels JK < KTROP, while the interactive
!     tropospheric aerosols (M7, AER) are only used at the levels
!     JK >= KTROP. The same routine is used in AER_PHY3_LAYER on the
!     model grid, so that the aerosol column optical depth diagnostics
!     only include the levels where the interactive aerosols are used.
!
!     The tropopause follows the WMO (1975) definition: the lowest
!     level at which the lapse rate decreases to 2 K/km or less,
!     provided that the average lapse rate between this level and all
!     higher levels within 2 km does not exceed 2 K/km (see Reichler
!     et al., 2003). Only levels between 75 and 550 hPa are
!     considered. If the criterion fails, the tropopause is set at
!     100 hPa.

!     INTERFACE.
!     ----------
!     CALL *CMIP_STRATO_TROPOPAUSE* FROM *CMIP_STRATO_AERO_PROCESS*
!     AND *AER_PHY3_LAYER*

!     AUTHOR.
!     -------
!     M.Ménégoz for EC-EARTH 2017
!     A.Laakso for EC-EARTH4/OIFS 2022
!     Moved out of CMIP_STRATO_AERO_PROCESS, 2026

!---------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RD, RG

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KLON
INTEGER(KIND=JPIM),INTENT(IN) :: KLEV
REAL(KIND=JPRB)   ,INTENT(IN) :: PAPF(KLON,KLEV)    ! full-level pressure (Pa)
REAL(KIND=JPRB)   ,INTENT(IN) :: PAPRS(KLON,KLEV+1) ! half-level pressure (Pa), top to surface
REAL(KIND=JPRB)   ,INTENT(IN) :: PT(KLON,KLEV)      ! full-level temperature (K)
REAL(KIND=JPRB)   ,INTENT(IN) :: PTH(KLON,KLEV+1)   ! half-level temperature (K), top to surface
INTEGER(KIND=JPIM),INTENT(OUT):: KTROP(KLON)        ! highest model level in the troposphere

REAL(KIND=JPRB)    :: ZSIG, ZLAPSE_RATE, ZLAPSE_RATE_MEAN, ZDZ
REAL(KIND=JPRB)    :: ZTRPAUS(KLON), ZPAPHD(KLON)
INTEGER(KIND=JPIM) :: JK, JL, JLEV
LOGICAL            :: LLTEST(KLON)
REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('CMIP_STRATO_TROPOPAUSE',0,ZHOOK_HANDLE)

DO JL=KIDIA,KFDIA
  ZTRPAUS(JL)=0.1_JPRB
  ZPAPHD(JL)=1.0_JPRB/PAPRS(JL,KLEV+1)
ENDDO

LLTEST(:)=.FALSE.

DO JK=KLEV,1,-1
  DO JL=KIDIA,KFDIA
    ! As in cloudsc.F90, the pressure level is normalized by the
    ! surface pressure before checking that the level is located on
    ! the tropopause altitude.
    ZSIG=PAPF(JL,JK)*ZPAPHD(JL)
    IF (PAPF(JL,JK)>7500._JPRB.AND.PAPF(JL,JK)<55000._JPRB) THEN
      ! The lapse rate is computed for each layer (here we use the gas
      ! equation and the hydrostatic approximation):
      ZLAPSE_RATE=(PTH(JL,JK+1)-PTH(JL,JK))/(PAPRS(JL,JK+1)-PAPRS(JL,JK))*PAPF(JL,JK)/PT(JL,JK)*RG/RD
      IF (ZLAPSE_RATE <= 0.002_JPRB) THEN
        ZDZ=0
        JLEV=JK
        DO WHILE (ZDZ <= 2000_JPRB) ! Compute the delta(P) corresponding to 2km
          ZDZ=ZDZ+(PAPF(JL,JLEV)-PAPF(JL,JLEV-1))/PAPRS(JL,JLEV)*PTH(JL,JLEV)*RD/RG
          JLEV=JLEV-1
        END DO
        ZLAPSE_RATE_MEAN=(PT(JL,JK)-PT(JL,JLEV))/ZDZ ! Lapse rate over 2 km
        IF (ZLAPSE_RATE_MEAN <= 0.002_JPRB) THEN
          ZTRPAUS(JL)=MAX(ZSIG,ZTRPAUS(JL)) ! We keep the lowest value that fits the two conditions
          LLTEST(JL)=.TRUE.
        ENDIF
      ENDIF
    ENDIF
  ENDDO
ENDDO

! Highest level below the tropopause.
! Use default tropopause at 100 hPa in case the above tropopause criteria fails
KTROP(KIDIA:KFDIA)=KLEV
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZSIG=PAPF(JL,JK)*ZPAPHD(JL)
    IF ((LLTEST(JL) .AND. ZSIG>ZTRPAUS(JL)) .OR. &
    &   (.NOT.LLTEST(JL) .AND. PAPF(JL,JK)>10000._JPRB)) THEN
      KTROP(JL)=MIN(KTROP(JL),JK)
    ENDIF
  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('CMIP_STRATO_TROPOPAUSE',1,ZHOOK_HANDLE)

END SUBROUTINE CMIP_STRATO_TROPOPAUSE
