SUBROUTINE HAMM7_DIAG_PM &
 &( YDMODEL, KIDIA  , KFDIA  , KLON   , KLEV , NAERO, &
 &  PAEROP,PAEPM1,PAEPM25,PAEPM10,PRHO, PDRYRADIUS, PWETRADIUS, PRHOP) 

! ╭────────────────────────────────────────────────────────────────────────────╮
! │                                                                            │
! │ Purpose :                                                                  │
! │ -------                                                                    │
! │  *hamm7_diag_pm* -                                                         │
! │                                                                            │
! │                                                                            │
! │ Interface :                                                                │
! │ ---------                                                                  │
! │   *HAMM7_DIAG_PM* is called from HAMM7_INTERFACE                           │
! │                                                                            │
! │                                                                            │
! │ Input :                                                                    │
! │ -----                                                                      │
! │  Dimensions                                                                │
! │  Aerosol mass mixing                                                       │
! │  aerosol mode sizes                                                        │
! │  aerosol densitities                                                       │
! │                                                                            │
! │ Output :                                                                   │
! │ ------                                                                     │
! │    PAEPM1   PM1 consentration                                              │
! │    PAEPM25  PM2.5 consentration                                            │
! │    PAEPM10  PM10 consentration                                             │
! │                                                                            │
! │ Method :                                                                   │
! │ ------                                                                     │
! │   Calcluate particulate concentration for given                            │
! │   aerosol size limits. Takes into account that observations.               │
! │   select particles with aerodynmic diameter while mass is calculated       │
! │   from dry mass.                                                           │
! │                                                                            │
! │ Reference :                                                                │
! │ ---------                                                                  │
! │                                                                            │
! │ Author :                                                                   │
! │ -------                                                                    │
! │  Orginal version:                                                          │
! │  Tommi Bergman (FMI)                                                       │
! │                                                                            │
! │ Modifications :                                                            │
! │ -------------                                                              │
! │                                                                            │
! ╰────────────────────────────────────────────────────────────────────────────╯
USE TYPE_MODEL,      ONLY : MODEL
USE MO_HAM,          ONLY : sigma_fine, sigma_coarse, nsol, nclass
USE PARKIND1,        ONLY : JPIM, JPRB
USE YOMHOOK,         ONLY : LHOOK, DR_HOOK, JPHOOK
USE mo_ham_m7ctl,    ONLY : cmedr2mmedr

IMPLICIT NONE

TYPE(MODEL)       , INTENT(IN)  :: YDMODEL                 ! For finding tracer indices
INTEGER(KIND=JPIM), INTENT(IN)  :: KLON 
INTEGER(KIND=JPIM), INTENT(IN)  :: KIDIA 
INTEGER(KIND=JPIM), INTENT(IN)  :: KFDIA                                           
INTEGER(KIND=JPIM), INTENT(IN)  :: KLEV                                            
INTEGER(KIND=JPIM), INTENT(IN)  :: NAERO                   ! Number of active aerosol species       
REAL(KIND=JPRB),    INTENT(IN)  :: PAEROP(KLON,KLEV,NAERO) ! aerosol mass mixing ratios

REAL(KIND=JPRB),    INTENT(OUT) :: PAEPM1(KLON),PAEPM25(KLON),PAEPM10(KLON)  ! PM output variables
REAL(KIND=JPRB),    INTENT(IN)  :: PDRYRADIUS(KLON,KLEV,NSOL)                ! Dry radii of particles (only for soluble classes)
REAL(KIND=JPRB),    INTENT(IN)  :: PWETRADIUS(KLON,KLEV,NCLASS)              ! Wet radii of particles (actually dry for insoluble)
REAL(KIND=JPRB),    INTENT(IN)  :: PRHOP(KLON,KLEV,NCLASS)                   ! Particle density
REAL(KIND=JPRB),    INTENT(IN)  :: PRHO(KLON,KLEV)                           ! Air density Kg(air)/M3 for unit conversion

!*      0.1 LOCAL VARIABLES

INTEGER(KIND=JPIM) :: JAER,  JL

! inidices for modes in local ZMASSDRYDIAMETER (also in pwetradius and dryradius (1-4))

INTEGER(KIND=JPIM), PARAMETER   :: I_NS = 1    ! Nucleation mode
!INTEGER(KIND=JPIM), PARAMETER   :: I_KS = 2    ! Aitken soluble (not needed at the moment)
INTEGER(KIND=JPIM), PARAMETER   :: I_AS = 3    ! accumulation soluble
INTEGER(KIND=JPIM), PARAMETER   :: I_CS = 4    ! coarse soluble
INTEGER(KIND=JPIM), PARAMETER   :: I_KI = 5    ! Aitken insoluble
INTEGER(KIND=JPIM), PARAMETER   :: I_AI = 6    ! accumulation insoluble
INTEGER(KIND=JPIM), PARAMETER   :: I_CI = 7    ! coarse insoluble

INTEGER(KIND=JPIM)              :: ICLASS      ! loop index
INTEGER(KIND=JPIM)              :: IPM         ! loop index

REAL(KIND=JPRB)                 :: ZPM_FRACTION(3)                    ! Fraction of mass in the PM range temporary for a mode
REAL(KIND=JPRB)                 :: ZPM_LIMIT(3)                       ! diameter limits for the three PM
REAL(KIND=JPRB)                 :: ZDRYRADIUS(KLON,KLEV,NCLASS)       ! count median dry radii of the M7 classes (modes)
REAL(KIND=JPRB)                 :: ZMASSDRYDIAMETER(KLON,KLEV,NCLASS) ! mass median dry diameters 

REAL(KIND=JPRB)                 :: ZHR     ! Temporary variable
REAL(KIND=JPRB)                 :: ZTMP    ! Temporary variable

INTRINSIC ERF

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('HAMM7_DIAG_PM',0,ZHOOK_HANDLE)
ASSOCIATE(  YAERO_NL         => YDMODEL%YRML_GCONF%YGFL%YAERO_NL )

!*   0.5 Initialisations

! output
PAEPM1(:)  = 0.0_JPRB
PAEPM25(:) = 0.0_JPRB
PAEPM10(:) = 0.0_JPRB

! local variables
ZDRYRADIUS(:,:,:)       = 0.0_JPRB
ZMASSDRYDIAMETER(:,:,:) = 0.0_JPRB

! Size limits
ZPM_LIMIT(1)    = 1.0e-6_JPRB ! sizelimit for PM1
ZPM_LIMIT(2)    = 2.5e-6_JPRB ! sizelimit for PM2.5
ZPM_LIMIT(3)    = 1.0e-5_JPRB ! sizelimit for PM10
ZPM_FRACTION(:) = 0.0_JPRB

! dry radii into one array
ZDRYRADIUS(KIDIA:KFDIA, 1:KLEV, I_NS:I_CS) = PDRYRADIUS(KIDIA:KFDIA, 1:KLEV, I_NS:I_CS) ! Dry radius defined only for soluble modes
ZDRYRADIUS(KIDIA:KFDIA, 1:KLEV, I_KI:I_CI) = PWETRADIUS(KIDIA:KFDIA, 1:KLEV, I_KI:I_CI) ! for Insoluble WET=DRY

! change count median radii to mass median diameter
do iclass=1,nclass
  ZMASSDRYDIAMETER(KIDIA:KFDIA, 1:KLEV, iclass) = 2.0_JPRB * ZDRYRADIUS(KIDIA:KFDIA, 1:KLEV, iclass) * cmedr2mmedr(iclass)
ENDDO

! determine PMs
DO JAER=1,NAERO
  IF ((INDEX(YAERO_NL(JAER)%CNAME ,'NS')>0 ) .or. (INDEX(YAERO_NL(JAER)%CNAME ,'KS')>0) .or.&
    & (INDEX(YAERO_NL(JAER)%CNAME ,'KI')>0 )  ) THEN
    ! exclude Number concentrations e.g. "xxxAS_N"
    IF(INDEX(YAERO_NL(JAER)%CNAME ,'_N') ==0 ) THEN
      DO JL=KIDIA,KFDIA
        PAEPM1(JL)  = PAEPM1(JL)  + PAEROP(JL,KLEV,JAER) * PRHO(JL,KLEV)
        PAEPM25(JL) = PAEPM25(JL) + PAEROP(JL,KLEV,JAER) * PRHO(JL,KLEV)
        PAEPM10(JL) = PAEPM10(JL) + PAEROP(JL,KLEV,JAER) * PRHO(JL,KLEV)
      ENDDO
    ENDIF
  ELSEIF((INDEX(YAERO_NL(JAER)%CNAME ,'AS')>0 ) )  THEN
    ! exclude Number concentrations e.g. "xxxAS_N"
    IF(INDEX(YAERO_NL(JAER)%CNAME ,'_N') ==0 ) THEN ! exclude number mixing ratios
      DO JL=KIDIA,KFDIA    
        if ((ZMASSDRYDIAMETER(JL,KLEV,I_AS))>1e-20_JPRB) then
          ZHR = 0.5 * SQRT(2.0)
          DO IPM=1,3
            ! the limit is multiplied to account for aerodynamic diameter (in contrast to geometric mean diameter)
            ZTMP = ( log(ZPM_LIMIT(IPM)*((1000./prhop(JL,KLEV,I_AS))**0.5)) &
            & - log(ZMASSDRYDIAMETER(JL,KLEV,I_AS)  ) ) / log(sigma_fine)
            ZPM_FRACTION(IPM) = 0.5 + 0.5 * ERF(ZTMP * ZHR)
          ENDDO
          PAEPM1(JL)  = PAEPM1(JL)  + ZPM_FRACTION(1) * PAEROP(JL,KLEV,JAER) * PRHO(JL,KLEV)
          PAEPM25(JL) = PAEPM25(JL) + ZPM_FRACTION(2) * PAEROP(JL,KLEV,JAER) * PRHO(JL,KLEV)
          PAEPM10(JL) = PAEPM10(JL) + ZPM_FRACTION(3) * PAEROP(JL,KLEV,JAER) * PRHO(JL,KLEV)    
        ENDIF
      ENDDO
    ENDIF
  ELSEIF((INDEX(YAERO_NL(JAER)%CNAME ,'AI')>0 ) )  THEN
    ! exclude Number concentrations e.g. "xxxAS_N"
    IF(INDEX(YAERO_NL(JAER)%CNAME ,'_N') ==0 ) THEN ! exclude number mixing ratios
      DO JL=KIDIA,KFDIA
        if ((ZMASSDRYDIAMETER(JL,KLEV,I_AI))>1e-20_JPRB) then
          ZHR = 0.5 * SQRT(2.0)
          DO IPM=1,3
            ! the limit is multiplied to account for aerodynamic diameter (in contrast to geometric mean diameter)
            ZTMP = ( log(ZPM_LIMIT(IPM)*((1000./prhop(JL,KLEV,I_AI))**0.5)) &
            & - log(ZMASSDRYDIAMETER(JL,KLEV,I_AI)  ) ) / log(sigma_fine)
            ZPM_FRACTION(IPM) = 0.5 + 0.5 * ERF(ZTMP * ZHR)
          ENDDO
          PAEPM1(JL)  = PAEPM1(JL)  + ZPM_FRACTION(1) * PAEROP(JL,KLEV,JAER) * PRHO(JL,KLEV)
          PAEPM25(JL) = PAEPM25(JL) + ZPM_FRACTION(2) * PAEROP(JL,KLEV,JAER) * PRHO(JL,KLEV)
          PAEPM10(JL) = PAEPM10(JL) + ZPM_FRACTION(3) * PAEROP(JL,KLEV,JAER) * PRHO(JL,KLEV)    
        ENDIF
      ENDDO
    ENDIF
  ELSEIF((INDEX(YAERO_NL(JAER)%CNAME ,'CS')>0 ) )  THEN
    ! exclude Number concentrations e.g. "xxxAS_N"
    IF(INDEX(YAERO_NL(JAER)%CNAME ,'_N') ==0 ) THEN ! exclude number mixing ratios
      DO JL=KIDIA,KFDIA
        if ((ZMASSDRYDIAMETER(JL,KLEV,I_CS))>1e-20_JPRB) then
          ZHR = 0.5 * SQRT(2.0)
          DO IPM=1,3
            ! the limit is multiplied to account for aerodynamic diameter (in contrast to geometric mean diameter)
            ZTMP = ( log(ZPM_LIMIT(IPM)*((1000./prhop(JL,KLEV,I_CS))**0.5)) &
            & - log( ZMASSDRYDIAMETER(JL,KLEV,I_CS) ) ) / log(sigma_coarse)
            ZPM_FRACTION(IPM) = 0.5 + 0.5 * ERF(ZTMP * ZHR)
          ENDDO
          PAEPM1(JL)  = PAEPM1(JL)  + ZPM_FRACTION(1) * PAEROP(JL,KLEV,JAER) * PRHO(JL,KLEV)
          PAEPM25(JL) = PAEPM25(JL) + ZPM_FRACTION(2) * PAEROP(JL,KLEV,JAER) * PRHO(JL,KLEV)
          PAEPM10(JL) = PAEPM10(JL) + ZPM_FRACTION(3) * PAEROP(JL,KLEV,JAER) * PRHO(JL,KLEV)    
        ENDIF
      ENDDO
    ENDIF
  ELSEIF (INDEX(YAERO_NL(JAER)%CNAME ,'CI')>0) THEN
    ! exclude Number concentrations e.g. "xxxAS_N"
    IF(INDEX(YAERO_NL(JAER)%CNAME ,'_N') ==0 ) THEN ! exclude number mixing ratios
      DO JL=KIDIA,KFDIA
        if ((ZMASSDRYDIAMETER(JL,KLEV,I_CI))>1e-20_JPRB) then
          ZHR = 0.5 * SQRT(2.0)
          DO IPM=1,3
            ! the limit is multiplied to account for aerodynamic diameter (in contrast to geometric mean diameter)
            ZTMP = ( log(ZPM_LIMIT(IPM)*((1000./prhop(JL,KLEV,I_CI))**0.5)) &
            & - log( ZMASSDRYDIAMETER(JL,KLEV,I_CI) ) ) / log(sigma_coarse)
            ZPM_FRACTION(IPM) = 0.5 + 0.5 * ERF(ZTMP * ZHR)
          ENDDO
          PAEPM1(JL)  = PAEPM1(JL)  + ZPM_FRACTION(1) * PAEROP(JL,KLEV,JAER) * PRHO(JL,KLEV)
          PAEPM25(JL) = PAEPM25(JL) + ZPM_FRACTION(2) * PAEROP(JL,KLEV,JAER) * PRHO(JL,KLEV)
          PAEPM10(JL) = PAEPM10(JL) + ZPM_FRACTION(3) * PAEROP(JL,KLEV,JAER) * PRHO(JL,KLEV)
        ENDIF
      ENDDO

    ENDIF
  ENDIF
ENDDO



END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('HAMM7_DIAG_PM',1,ZHOOK_HANDLE)
END SUBROUTINE HAMM7_DIAG_PM
