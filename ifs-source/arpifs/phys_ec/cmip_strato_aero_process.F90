SUBROUTINE CMIP_STRATO_AERO_PROCESS                      &
&(YDMODEL, KIDIA , KFDIA , KLON , KLEV,                           &
& KRINT , KSHIFT,                                         &
& KINDAT,KSTADD,                                          &
& PAPF, PAPRS , ZPT, ZPTH, PGELAM,                 &
& STRATAOD_M, STRATAAOD_M, STRATREFAOD_M, STRATAAOD_LW_M, &
& STRATAOD , STRATAAOD , STRATREFAOD, STRATAAOD_LW, TROPOPAUSE )

!     PURPOSE.
!     --------

!     INTERFACE.
!     ----------
!     CALL *CMIP_STRATO_AERO_PROCESS* FROM *RADLSWR*

!     AUTHOR.
!     -------
!     M.Ménégoz for EC-EARTH 2017
!     A.Laakso for EC-EARTH4/OIFS 2022 
!
!     Time interpolation of the stratospheric aerosol CMIP dataset
!     Remove the aerosols located below the tropopause and above 10
!     hPa. (WMO criteria based on the vertical lapse rate temperature)

!     MODIFICATIONS.
!     -------------- 

!---------------------------------------------------------------------

USE YOEAEROP ,ONLY : STRATO_CMIP_NLAT,STRATO_CMIP_NALT,STRATO_CMIP_NTIME, &
                & STRATO_CMIP_LAT,STRATO_CMIP_SB,STRATO_CMIP_TB,&
                & STRATO_CMIP_NSB,STRATO_CMIP_NTB,STRATO_CMIP_NMONTH

USE YOMCT0   , ONLY : LNF
USE YOMCT2   , ONLY : NSTAR2
USE YOMCST   , ONLY : RDAY
USE YOMLUN   , ONLY : NULOUT
USE PARKIND1 , ONLY : JPIM, JPIB, JPRB, JPRD !JBRD needed for fcttim.func
USE YOMHOOK,    ONLY: LHOOK, DR_HOOK, JPHOOK
USE TYPE_MODEL, ONLY: MODEL

! For the tropopause criteria:
USE YOMCST   , ONLY : RD, RG

IMPLICIT NONE

!     -----------------------------------------------------------------
!       ARGUMENTS.
!      -------------
TYPE(MODEL),       INTENT(IN) :: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN) :: KLON 
INTEGER(KIND=JPIM),INTENT(IN) :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN) :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN) :: KRINT !Alaak: not used 
INTEGER(KIND=JPIM),INTENT(IN) :: KSHIFT !Alaak: not used 
INTEGER(KIND=JPIM),INTENT(IN) :: KINDAT 
INTEGER(KIND=JPIM),INTENT(IN) :: KSTADD
REAL(KIND=JPRB)   ,INTENT(IN) :: PAPF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN) :: PAPRS(KLON,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(IN) :: ZPT(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN) :: ZPTH(KLON,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(IN) :: PGELAM(KLON)
REAL(KIND=JPRB)   ,INTENT(IN) :: STRATAOD_M(KLON,KLEV,STRATO_CMIP_NSB,STRATO_CMIP_NMONTH)
REAL(KIND=JPRB)   ,INTENT(IN) :: STRATAAOD_M(KLON,KLEV,STRATO_CMIP_NSB,STRATO_CMIP_NMONTH)
REAL(KIND=JPRB)   ,INTENT(IN) :: STRATREFAOD_M(KLON,KLEV,STRATO_CMIP_NSB,STRATO_CMIP_NMONTH)
REAL(KIND=JPRB)   ,INTENT(IN) :: STRATAAOD_LW_M(KLON,KLEV,STRATO_CMIP_NTB,STRATO_CMIP_NMONTH)
REAL(KIND=JPRB)   ,INTENT(OUT):: STRATAOD(KLON,KLEV,STRATO_CMIP_NSB)
REAL(KIND=JPRB)   ,INTENT(OUT):: STRATAAOD(KLON,KLEV,STRATO_CMIP_NSB)
REAL(KIND=JPRB)   ,INTENT(OUT):: STRATREFAOD(KLON,KLEV,STRATO_CMIP_NSB)
REAL(KIND=JPRB)   ,INTENT(OUT):: STRATAAOD_LW(KLON,KLEV,STRATO_CMIP_NTB)
INTEGER(KIND=JPIM),INTENT(OUT):: TROPOPAUSE(KLON)

!     -----------------------------------------------------------------
!       LOCAL ARRAYS.
!      -------------

INTEGER(KIND=JPIB) :: ITIME, IZT
INTEGER(KIND=JPIM) :: ISTADD
INTEGER(KIND=JPIM) :: IMV1, IMV2, IYR1, IYR2, INYR, INDY, INMN, IT1, IT2, IMP1, IMP2
INTEGER(KIND=JPIM) :: ILMO(12)
INTEGER(KIND=JPIM) :: IDY0,IMN0,IYR0,IYNR
INTEGER(KIND=JPIM) :: JK,JL
REAL(KIND=JPRB)    :: ZW1, ZW2, ZSIG
REAL(KIND=JPRB)    :: ZTRPAUS(KLON), ZPAPHD(KLON)
REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE

! For the Reichler criteria (tropopause):
REAL(KIND=JPRB)    :: LAPSE_RATE, LAPSE_RATE_MEAN, DZ
INTEGER(KIND=JPIM) :: JLEV
LOGICAL            :: LTEST(KLON)
INTEGER(KIND=JPIM) :: RVOLCDEC

!---included functions
#include "fcttim.func.h" 
#include "updcal.intfb.h"

IF (LHOOK) CALL DR_HOOK('CMIP_STRATO_AERO_PROCESS',0,ZHOOK_HANDLE)

WRITE(NULOUT,*)'Entering Subroutine CMIP_STRATO_AERO_PROCESS'

ASSOCIATE(&
  &  YDPHYAER     => YDMODEL%YRML_PHY_AER,  &
  &LTWOTL=>YDMODEL%YRML_DYN%YRDYNA%LTWOTL,&
  &YDRIP=>YDMODEL%YRML_GCONF%YRRIP &
  )

ASSOCIATE(&
    & YDAERSTRAT        => YDPHYAER%YREAEROSTRAT &
    )
ASSOCIATE(&
     STRATO_CMIP_ALT=>YDAERSTRAT%STRATO_CMIP_ALT&
    )

STRATAOD    (KIDIA:KFDIA,:,:) = 0.0_JPRB
STRATAAOD   (KIDIA:KFDIA,:,:) = 0.0_JPRB
STRATREFAOD (KIDIA:KFDIA,:,:) = 0.0_JPRB
STRATAAOD_LW(KIDIA:KFDIA,:,:) = 0.0_JPRB

! Time interpolation
! Assume monthly means valid at 0Z on the 16th of each month
! and interpolate linearly

! prepare monthly weights
IF(.NOT.LNF.AND.KSTADD == 0) THEN
  ! IN CASE OF RESTART:
  IF (LTWOTL) THEN
    IZT=NINT(YDRIP%TSTEP*(REAL(NSTAR2,JPRB)+0.5_JPRB), KIND=JPIB)
   ELSE
     ITIME=NINT(YDRIP%TSTEP, KIND=JPIB)
     IZT=ITIME*INT(NSTAR2, KIND=JPIB)
  ENDIF
  ISTADD=INT(IZT/NINT(RDAY, KIND=JPIB), KIND=JPIM)
 ELSE
  ISTADD=KSTADD ! NUMBER OF DAYS SINCE START OF THE MODEL
ENDIF

IYR0=NCCAA(KINDAT)
IMN0=NMM(KINDAT)
IDY0=NDD(KINDAT)

CALL UPDCAL(IDY0,IMN0,IYR0, ISTADD, INDY,INMN,INYR, ILMO, NULOUT)

! Assume monthly means valid at 0Z on the 16th of each month
! and interpolate linearly (daily values)

IF(INDY >= 16) THEN
  ! Month to consider to apply the interpolation
  IMP1=2
  IMP2=3
  !
  IMV1=INMN
  IMV2=1+MOD(INMN,12)
  IYR1=INYR
  IF(IMV1 /= 12) THEN
    IYR2=INYR
   ELSE
    IYR2=INYR+1
  ENDIF
  IT1=16
  IT2=16+ILMO(IMV1)
 ELSE
  ! Month to consider to apply the interpolation
  IMP1=1
  IMP2=2
  !
  IMV1=1+MOD(INMN+10,12)
  IMV2=INMN
  IF(IMV1 /= 12) THEN
    IYR1=INYR
   ELSE
    IYR1=INYR-1
  ENDIF
  IYR2=IYNR
  IT1=16-ILMO(IMV1)
  IT2=16
ENDIF

! Weights
ZW1=REAL(IT2-INDY,JPRB)/REAL(IT2-IT1,JPRB)
ZW2=1.0_JPRB-ZW1

! Interpolation and application of the RVOLCDEC factor (2-year
! exponential decay if NVOLCSTOP different from 99999999)
! --ALaakso: no exponential decay, RVOLCDEC is set to 1
RVOLCDEC=1
STRATAOD    (KIDIA:KFDIA,:,:) = RVOLCDEC*(ZW1*STRATAOD_M    (KIDIA:KFDIA,:,:,IMP1) + ZW2*STRATAOD_M    (KIDIA:KFDIA,:,:,IMP2))
STRATAAOD   (KIDIA:KFDIA,:,:) = RVOLCDEC*(ZW1*STRATAAOD_M   (KIDIA:KFDIA,:,:,IMP1) + ZW2*STRATAAOD_M   (KIDIA:KFDIA,:,:,IMP2))
STRATREFAOD (KIDIA:KFDIA,:,:) = RVOLCDEC*(ZW1*STRATREFAOD_M (KIDIA:KFDIA,:,:,IMP1) + ZW2*STRATREFAOD_M (KIDIA:KFDIA,:,:,IMP2))
STRATAAOD_LW(KIDIA:KFDIA,:,:) = RVOLCDEC*(ZW1*STRATAAOD_LW_M(KIDIA:KFDIA,:,:,IMP1) + ZW2*STRATAAOD_LW_M(KIDIA:KFDIA,:,:,IMP2))

! Computing the tropopause height
DO JL=KIDIA,KFDIA
  ZTRPAUS(JL)=0.1_JPRB
  ZPAPHD(JL)=1.0_JPRB/PAPRS(JL,KLEV+1)
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! According to the WMO (1975), the tropopause is defined as the lowest
! level at which the lapse-rate decreases to 2°C/km or less, provided
! that the average lapse-rate between this level and all higher levels
! within 2 km does not exceed 2°C/km.(see Reichler et al.,2003)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
LTEST(:)=.FALSE.

DO JK=KLEV,1,-1
  DO JL=KIDIA,KFDIA
    ! As in cloudsc.F90, the pressure level is normalized by the
    ! surface pressure before checking that the level is located on
    ! the tropopause altitude.
    ZSIG=PAPF(JL,JK)*ZPAPHD(JL)
    IF (PAPF(JL,JK)>7500._JPRB.AND.PAPF(JL,JK)<55000._JPRB) THEN
      ! The lapse rate is computed for each layer (here we use the gas
      ! equation and the hydrostatic approximation):
      LAPSE_RATE=(ZPTH(JL,JK+1)-ZPTH(JL,JK))/(PAPRS(JL,JK+1)-PAPRS(JL,JK))*PAPF(JL,JK)/ZPT(JL,JK)*RG/RD
      IF (LAPSE_RATE <= 0.002_JPRB) THEN
        DZ=0
        JLEV=JK
        DO WHILE (DZ <= 2000_JPRB) ! Compute the delta(P) corresponding to 2km
          DZ=DZ+(PAPF(JL,JLEV)-PAPF(JL,JLEV-1))/PAPRS(JL,JLEV)*ZPTH(JL,JLEV)*RD/RG
          JLEV=JLEV-1
        END DO
        LAPSE_RATE_MEAN=(ZPT(JL,JK)-ZPT(JL,JLEV))/DZ ! Lapse rate over 2 km
        IF (LAPSE_RATE_MEAN <= 0.002_JPRB) THEN
          ZTRPAUS(JL)=MAX(ZSIG,ZTRPAUS(JL)) ! We kep the lowest value that fit the two conditions
          LTEST(JL)=.TRUE.
        ENDIF
      ENDIF
    ENDIF
  ENDDO
ENDDO

! Removing the stratospheric forcing that is located below the tropopause 
! Use default tropopause at 100 hPa in case the above tropopause criteria fails
TROPOPAUSE(KIDIA:KFDIA)=KLEV
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZSIG=PAPF(JL,JK)*ZPAPHD(JL)
    IF ((LTEST(JL) .AND. ZSIG>ZTRPAUS(JL)) .OR. &
    &   (.NOT.LTEST(JL) .AND. PAPF(JL,JK)>10000._JPRB)) THEN
      STRATAOD(JL,JK,:)=0.0_JPRB
      STRATAAOD(JL,JK,:)=0.0_JPRB  
      STRATREFAOD(JL,JK,:)=0.0_JPRB
      STRATAAOD_LW(JL,JK,:)=0.0_JPRB
      TROPOPAUSE(JL)=MIN(TROPOPAUSE(JL),JK)
    ENDIF
  ENDDO
ENDDO

END ASSOCIATE
END ASSOCIATE
END ASSOCIATE
WRITE(NULOUT,*)'Leaving subroutine CMIP_STRATO_AERO_PROCESS'
IF (LHOOK) CALL DR_HOOK('CMIP_STRATO_AERO_PROCESS',1,ZHOOK_HANDLE)

END SUBROUTINE CMIP_STRATO_AERO_PROCESS
