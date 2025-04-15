  MODULE YOE_AER_ACTIV

  !---inherited functions, types, variables and constants 
  USE PARKIND1,            ONLY: JPIM,JPRB
  !USE PHY_DIAG_MOD,        ONLY: T_DIAG

  IMPLICIT NONE

  !---public member functions
  PUBLIC AER_ACTIV 
  !PUBLIC GET_CDNC_FACTOR
  !PUBLIC SETUP_ACI_DIAG

  !---private member functions
  !PRIVATE AER_ACTIV_FOUNTOUKIS_NENES
  !PRIVATE AER_ACTIV_MORALES_NENES
  PRIVATE AER_ACTIV_MORALES_NENES_FULL
  !PRIVATE AER_ACTIV_MORALES_NENES_FULL_OLDPDF
  !PRIVATE AER_ACTIV_ABDULRAZZAK_GHAN
  !PRIVATE AER_ACTIV_MENON
  !PRIVATE SINTEGRAL                    
  !PRIVATE SINTEGRAL_MN
  !PRIVATE DIAGNOSE_AEROSOL_MASS
  PRIVATE GET_HAMM7_AERO_PROP
  !PRIVATE LIQ_CLOUD_RE 
  PRIVATE ICE_CLOUD_PROP
  !PRIVATE PDF_UPDRAFT

  !---module types, variables and constants

  ! surface tension of pure water at 273.15K [J m-2]
  !REAL(KIND=JPRB), PARAMETER, PUBLIC :: PPSURFTEN = 0.075_JPRB   

  ! standard deviation of updraft PDF
  !REAL(KIND=JPRB), PARAMETER, PUBLIC :: PPDF_SIGMA = 0.8_JPRB

  ! minimum CDNC value
  REAL(KIND=JPRB), PARAMETER, PUBLIC :: PPMINCDNC = 1._JPRB

  ! default CDNC value
  REAL(KIND=JPRB), PARAMETER, PUBLIC :: PPDEFCDNC = 1._JPRB !1.E-14_JPRB !10._JPRB

  ! Diagnostics
  !TYPE(T_DIAG), POINTER, PUBLIC :: D_CDNC
  !TYPE(T_DIAG), POINTER, PUBLIC :: D_REFF
  !TYPE(T_DIAG), POINTER, PUBLIC :: D_LIQCLDT

  ! default value of effective radius (value when gridbox is cloud-free) [um]
  REAL(KIND=JPRB), PARAMETER, PRIVATE :: PPREFFL_DEF = 4.0_JPRB
 
CONTAINS
  SUBROUTINE AER_ACTIV(KIDIA,   KFDIA,  KTDIA,   KLON,    KLEV,   KSTGLO, &
                     !&  KLEVX,   KFLDX,  KFLDX2,                   &
                     &  PAPH,    PAP,    PT,      PQ,      PQSAT,  &
                     &  PVERVEL, PA,     PL,      PI,              &
                     &  PLSM,    PGELAM,   PGEMU, & !PSLON,   PGEMU,  &
                     &  PGFL, YDMODEL, PCDNCACT, PICNC, PREFFL, PREFFI, PSMAX, PDRYRSOLU, PXTM1, KTRAC, PSIGMA_W) !,    PEXTRA, PEXTR2)
   
   ! *AER_ACTIV* is the interface to the cloud droplet activation scheme. 
   !  Four schemes are available, depending on the aerosol scheme used.
   !  1. Menon et al. for TM5 aerosols or Tegen climatological aerosols
   !  2. Abdul-Razzak & Ghan (2000)
   !  3. Fountoukis & Nenes (2005)
   !  4. Morales Betancourt & Nenes (2014)
   !  5. Morales Betancourt & Nenes (2014) with adsorption activation (Kumar et al., 2011), PDF integration with Gaussian-Legendre quadrature, explicit CCN spectra
   !  6. Same as 5. but with PDF sampling as in 2-4
   !  Schemes 2-6 currently require TM5 aerosols,
   !  either via interactive coupling to TM5 or a prescribed pre-industrial climatology
      
   !---inherited functions, types, variables and constants
   USE YOMHOOK,             ONLY: LHOOK, DR_HOOK
   USE YOMCST,              ONLY: RD, RPI, RG
   !  USE YOERAD,              ONLY: LNEWAER, LCMIP5, LTM5AER, LCMIP6_PI_AEROSOLS, LMAC2SPACI
   !  USE YOECLDP,             ONLY: NCLOUDACT, JP_ACT_FOUNTOUKIS_NENES, JP_ACT_ABDULRAZZAK_GHAN, & 
   !                               & JP_ACT_MENON, JP_ACT_MORALES_NENES, JP_ACT_MORALES_NENES_FULL, &
   !                               & JP_ACT_MORALES_NENES_FULL_OLDPDF, &
   !                               & LAERICESED, LAERICEAUTO, LACI_DIAG, NACTPDF,   &
   !                               & RLMIN, RAMIN, RTHOMO
   !USE YOECLDP,             ONLY: LAERICESED, LAERICEAUTO,   &
   !                             & RLMIN, RAMIN, RTHOMO
   !USE MO_ACTIV,            ONLY: nw !eehol: to replace NACTPDF
   !USE YOE_AERO_M7_DATA,    ONLY: NSOL
   USE TM5M7_DATA,          ONLY: NSOL
   !USE YOMCT3,              ONLY: NSTEP
   !USE YOMCT0,              ONLY: NFRPOS
   !USE YOMDYN,              ONLY: TSTEP
   !USE YOM_YGFL,            ONLY: YGFL, YCDNC, YICNC, YRE_LIQ, YRE_ICE
   !USE YOE_PI_AERO,         ONLY: LPI_AERO_UPDATED
   USE TYPE_MODEL,          ONLY: MODEL

   IMPLICIT NONE
  
   !---included functions from header files
#include "abor1.intfb.h"
   !#include "ice_effective_radius.intfb.h"

   !---subroutine interface
   !
   ! -   Input arguments
   !     ---------------
   ! KIDIA    : start of horizontal loop
   ! KFDIA    : end   of horizontal loop
   ! KLON     : horizontal dimension
   ! KLEV     : vertical dimension
   ! KSTGLO   : offset of horizontal block in coupling arrays
   ! KLEVX    : vertical dimension, 3-D diagnostic fields
   ! KFLDX    : number of 3-D diagnostic fields
   ! KFLDX2   : number of 2-D diagnostic fields
   ! PAPH     : pressure at half levels [Pa]
   ! PAP      : pressure at full levels [Pa]
   ! PT       : temperature [K]
   ! PQ       : specific humidity [kg/kg]
   ! PQSAT    : saturation specific humidity [kg/kg]
   ! PVERVEL  : large-scale vertical veloctiy [Pa/s]
   ! PA       : cloud fraction
   ! PL       : liquid water content   
   ! PI       : ice content [kg/kg]
   ! PLSM     : land-sea mask
   ! PGELAM   : Longitude
   ! PCLON    : Cosine of longitude
   ! PSLON    : Sine of longitude
   ! PGEMU    : Sine of latitude

   ! -   output arguments
   !     ----------------
   ! The following are written to the appropriate PGFL:
   ! CDNC   : cloud droplet number concentration
   ! ICNC   : ice crystal number concentration 
   ! PRE_LIQ : liquid drop effective radius
   ! PRE_ICE : ice crystal effective radius
   ! PEXTRA  : 3-D diagnostics
   ! PEXTR2  : 2-D diagnostics

   TYPE(MODEL),        INTENT(IN) :: YDMODEL 
   INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA                  ! beginning of horizontal block
   INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA                  ! beginning of horizontal block
   INTEGER(KIND=JPIM), INTENT(IN) :: KTDIA                  ! beginning of vertical block
   INTEGER(KIND=JPIM), INTENT(IN) :: KLON                   ! horizontal dimension
   INTEGER(KIND=JPIM), INTENT(IN) :: KLEV                   ! number of model vertical levels
   INTEGER(KIND=JPIM), INTENT(IN) :: KSTGLO                 ! offset of horizontal block in coupling arrays
   INTEGER(KIND=JPIM), INTENT(IN) :: KTRAC                  ! number of tracers 

   !INTEGER(KIND=JPIM), INTENT(IN) :: KLEVX
   !INTEGER(KIND=JPIM), INTENT(IN) :: KFLDX
   !INTEGER(KIND=JPIM), INTENT(IN) :: KFLDX2

   REAL(KIND=JPRB), INTENT(IN)    :: PAPH(KLON,KLEV+1)      ! half level pressure
   REAL(KIND=JPRB), INTENT(IN)    :: PAP(KLON,KLEV)         ! full level pressure
   REAL(KIND=JPRB), INTENT(IN)    :: PT(KLON,KLEV)          ! temperature
   REAL(KIND=JPRB), INTENT(IN)    :: PQ(KLON,KLEV)          ! specific humidity
   REAL(KIND=JPRB), INTENT(IN)    :: PQSAT(KLON,KLEV)       ! saturation specific humidity
   REAL(KIND=JPRB), INTENT(IN)    :: PVERVEL(KLON,KLEV)     ! vertical velocity

   REAL(KIND=JPRB), INTENT(IN)    :: PA(KLON,KLEV)          ! cloud fraction
   REAL(KIND=JPRB), INTENT(IN)    :: PL(KLON,KLEV)          ! cloud liquid water
   REAL(KIND=JPRB), INTENT(IN)    :: PI(KLON,KLEV)          ! cloud ice

   REAL(KIND=JPRB), INTENT(IN)    :: PLSM(KLON)             ! land-sea mask
   REAL(KIND=JPRB), INTENT(IN)    :: PGELAM(KLON)           ! longitude
   !REAL(KIND=JPRB), INTENT(IN)    :: PCLON(KLON)            ! cosine of longitude
   !REAL(KIND=JPRB), INTENT(IN)    :: PSLON(KLON)            ! sine of longitude
   REAL(KIND=JPRB), INTENT(IN)    :: PGEMU(KLON)            ! sine of latitude
   REAL(KIND=JPRB), INTENT(IN)    :: PDRYRSOLU(KLON,KLEV,NSOL) ! rdry of soluble modes [m]
   REAL(KIND=JPRB), INTENT(IN)    :: PXTM1(KLON,KLEV,KTRAC) ! tracer mixing ratios
   REAL(KIND=JPRB), INTENT(IN)    :: PSIGMA_W(KLON,KLEV)    ! sigma_w

   REAL(KIND=JPRB), INTENT(INOUT) :: PGFL(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NDIM) !YGFL%NDIM)
   REAL(KIND=JPRB), INTENT(INOUT) :: PCDNCACT(KLON,KLEV) ! cloud droplet number concentration [#/cm-3]
   REAL(KIND=JPRB), INTENT(INOUT) :: PICNC(KLON,KLEV) ! ice crystal number concentration [#/cm-3]
   REAL(KIND=JPRB), INTENT(INOUT) :: PREFFL(KLON,KLEV) ! liquid droplet effective radius [um]
   REAL(KIND=JPRB), INTENT(INOUT) :: PREFFI(KLON,KLEV) ! ice effective radius [um]
   REAL(KIND=JPRB), INTENT(INOUT) :: PSMAX(KLON,KLEV) ! maximum supersaturation [%]
   
   !---extra diagnostics
   !REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEXTRA(KLON,KLEVX,KFLDX) 
   !REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEXTR2(KLON,KFLDX2) 

   !---local variables
   REAL(KIND=JPRB) :: ZRHO(KLON,KLEV)                       ! air density [kg/m3]
   REAL(KIND=JPRB) :: ZSO4MASS(KLON,KLEV,NSOL)              ! modal so4 mass mixing ratio [kg/kg]
   REAL(KIND=JPRB) :: ZBCMASS(KLON,KLEV,NSOL)               ! modal black carbon mass mixing ratio [kg/kg]
   REAL(KIND=JPRB) :: ZOMMASS(KLON,KLEV,NSOL)               ! modal organic matter mass mixing ratio [kg/kg]
   REAL(KIND=JPRB) :: ZSSMASS(KLON,KLEV,NSOL)               ! modal sea salt mass mixing ratio [kg/kg]
   REAL(KIND=JPRB) :: ZDUMASS(KLON,KLEV,NSOL)               ! modal dust mass mixing ratio [kg/kg]
   REAL(KIND=JPRB) :: ZNO3MASS(KLON,KLEV)                   ! nitrate mass mixing ratio [kg/kg]
   REAL(KIND=JPRB) :: ZMSAMASS(KLON,KLEV)                   ! MSA mass mixing ratio [kg/kg]
   REAL(KIND=JPRB) :: ZAERONUM(KLON,KLEV,NSOL)              ! aerosol number mixing ratio [#/kg]
   REAL(KIND=JPRB) :: ZDRYRSOL(KLON,KLEV,NSOL)              ! dry count median radius of soluble modes
   REAL(KIND=JPRB) :: ZDRYRSOLOLD(KLON,KLEV,NSOL)           ! dry count median radius of soluble modes eehol: not calculated here but in HAM
   REAL(KIND=JPRB) :: ZSO4BULK(KLON,KLEV)                   ! bulk so4 mass mixing ratio [kg/kg]
   REAL(KIND=JPRB) :: ZBCBULK(KLON,KLEV)                    ! bulk black carbon mass mixing ratio [kg/kg]
   REAL(KIND=JPRB) :: ZOMBULK(KLON,KLEV)                    ! bulk organic matter mass mixing ratio [kg/kg]
   REAL(KIND=JPRB) :: ZSSBULK(KLON,KLEV)                    ! bulk sea salt mass mixing ratio [kg/kg]
   REAL(KIND=JPRB) :: ZDUBULK(KLON,KLEV)                    ! bulk dust mass mixing ratio [kg/kg]
   !REAL(KIND=JPRB) :: ZW(KLON,KLEV,nw)                 ! updraft speed [m/s]
   !REAL(KIND=JPRB) :: ZWPDF(KLON,KLEV,nw)              ! updraft probability
   REAL(KIND=JPRB) :: ZHOOK_HANDLE
   REAL(KIND=JPRB) :: ZTMP(KLON,KLEV)                       ! interim for diagnostics

   REAL(KIND=JPRB) :: ZCDNC(KLON,KLEV)                      ! cloud droplet number concentration [#/cm-3]
   REAL(KIND=JPRB) :: ZICNC(KLON,KLEV)                      ! ice crystal number concentration [#/cm-3]
   REAL(KIND=JPRB) :: ZRE_LIQ(KLON,KLEV)                    ! liquid droplet effective radius [um]
   REAL(KIND=JPRB) :: ZRE_ICE(KLON,KLEV)                    ! ice crystal effective radius [um]
   REAL(KIND=JPRB) :: ZSMAX(KLON,KLEV)                    ! maximum supersaturation [%]

   !variables for liquid droplet eff rad calculations
   REAL(KIND=JPRB) :: ZQLWC(KLON,KLEV) !tresholded LWC [kg/kg]
   REAL(KIND=JPRB) :: RCLDMAX=5.E-3_JPRB !max cloud water

   LOGICAL :: LLIQCLD(KLON,KLEV)                            ! true if liquid cloud is present
   LOGICAL :: LLIQCLDD(KLON,KLEV)                           ! true if liquid cloud is present (for activation calculations)
   LOGICAL :: LICECLD(KLON,KLEV)                            ! true if ice cloud is present

   LOGICAL :: LL1
   LOGICAL :: LBULK, LMODE                                  ! fetch TM5 aerosols as bulk mass / per-mode mass

   LOGICAL :: LCALCINCLOUD = .TRUE.                        ! eehol: logical to calculate activation only in-cloud (T) or everywhere (F)
   REAL(KIND=JPRB) :: ZEPS                                  ! epsilon(1.)
   INTEGER(KIND=JPIM) :: ITOP                               ! highest level for water cloud
   INTEGER(KIND=JPIM) :: JK, JL, JMOD                       ! loop indices
   !INTEGER(KIND=JPIM) :: IX                                 ! index to PEXTRA diagnostic

   REAL(KIND=JPRB), PARAMETER :: ZEPSEC  = 1.E-14_JPRB      ! taken from cloudsc.F90
   REAL(KIND=JPRB) :: ZTMPA

   REAL(KIND=JPRB) :: ZMAC2SP_CDNC_FACTOR(KLON)

   !---executable procedure
   IF (LHOOK) CALL DR_HOOK('YOE_AER_ACTIV.AER_ACTIV', 0, ZHOOK_HANDLE)
   ASSOCIATE(YGFL=>YDMODEL%YRML_GCONF%YGFL,YDECLDP=>YDMODEL%YRML_PHY_EC%YRECLDP)
   ASSOCIATE(YCDNC=>YGFL%YCDNC, YICNC=>YGFL%YICNC, YRE_LIQ=>YGFL%YRE_LIQ, YRE_ICE=>YGFL%YRE_ICE, &
      & LAERICESED=>YDECLDP%LAERICESED, LAERICEAUTO=>YDECLDP%LAERICEAUTO, &
      & RLMIN=>YDECLDP%RLMIN, RAMIN=>YDECLDP%RAMIN, RTHOMO=>YDECLDP%RTHOMO, RNICE=>YDECLDP%RNICE)

   !---air density
   DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
         ZRHO(JL,JK)=PAP(JL,JK)/(RD*PT(JL,JK))
         IF (ZRHO(JL,JK) .LE. 1.E-12_JPRB) THEN
            ZRHO(JL,JK)=1.E-12_JPRB
         END IF
      END DO
   END DO

   !---background CDNC and Reff and ICNC and init ZSMAX
   ZSMAX(KIDIA:KFDIA,1:KLEV) = 0._JPRB
   ZCDNC(KIDIA:KFDIA,1:KLEV) = PPDEFCDNC
   ZRE_LIQ(KIDIA:KFDIA,1:KLEV) = PPREFFL_DEF
   ZICNC(KIDIA:KFDIA,1:KLEV) = RNICE
   ZRE_ICE(KIDIA:KFDIA,1:KLEV) = 80._JPRB*0.64952_JPRB
  !PGFL(KIDIA:KFDIA,:,YCDNC%MP9_PH) = PPDEFCDNC
  !PGFL(KIDIA:KFDIA,:,YRE_LIQ%MP9_PH) = PPREFFL_DEF
  !PGFL(KIDIA:KFDIA,:,YICNC%MP9_PH) = ZICNC(KIDIA:KFDIA,:)
  !PGFL(KIDIA:KFDIA,:,YRE_ICE%MP9_PH) = PREFFI(KIDIA:KFDIA,:)

  !---get aerosols from TM5
  !IF (LTM5AER .OR. LCMIP6_PI_AEROSOLS) THEN
   LMODE = .TRUE.
   LBULK = (LAERICESED .OR. LAERICEAUTO) ! .OR. NCLOUDACT==JP_ACT_MENON)
   CALL GET_HAMM7_AERO_PROP(KIDIA, KFDIA, KLON, KTDIA, KLEV, KSTGLO, LMODE, LBULK, &
                        & PAPH,     PGELAM,  PGEMU, PXTM1, KTRAC, &
                        & ZDRYRSOLOLD, ZAERONUM, &
                        & ZSO4MASS, ZBCMASS, ZOMMASS, ZSSMASS, ZDUMASS, &
                        & ZSO4BULK, ZBCBULK, ZOMBULK, ZSSBULK, ZDUBULK, &
                        & ZNO3MASS, ZMSAMASS)
   !END IF

   ZDRYRSOL(KIDIA:KFDIA,1:KLEV,:) = PDRYRSOLU(KIDIA:KFDIA,1:KLEV,:)

  !---find cells where there is cloud (same criteria as cloudsc.F90 when it uses CDNC)

   DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
         ZTMPA = 1.0_JPRB/MAX(PA(JL,JK),ZEPSEC)
         LLIQCLD(JL,JK) = ( PL(JL,JK)*ZTMPA  ) > ZEPSEC
         LICECLD(JL,JK) = ( PI(JL,JK)*ZTMPA  ) > ZEPSEC ! logical for ice cloud
      END DO
   END DO

   IF (LCALCINCLOUD) THEN
      !---find highest model level where there is cloud
      DO JK=1,KLEV
         IF (ANY(LLIQCLD(KIDIA:KFDIA,JK))) EXIT
      END DO
      ITOP=JK
      LLIQCLDD(KIDIA:KFDIA,1:KLEV) = LLIQCLD(KIDIA:KFDIA,1:KLEV)
   ELSE
      ITOP=KTDIA
      LLIQCLDD(KIDIA:KFDIA,1:KLEV) = .TRUE.
   END IF

   IF (ITOP.LE.KLEV) THEN
!     
!     SELECT CASE(NCLOUDACT)
!     CASE(JP_ACT_MENON) 
!
!        IF (.NOT. (LTM5AER .OR. LCMIP6_PI_AEROSOLS .OR. LCMIP5 .OR. LNEWAER)) &
!           & CALL ABOR1('YOE_AER_ACTIV: No supported aerosol scheme for Menon cloud activation')
!
!        !---Tegen or CMIP5 aerosols: map from optical properties to bulk aerosol mass
!        IF ((.NOT.LTM5AER) .AND. (.NOT. LCMIP6_PI_AEROSOLS)) THEN
!           CALL DIAGNOSE_AEROSOL_MASS(KIDIA, KFDIA, KLON, KLEV, PT, PQ, PQSAT, &
!                                    & PAPH, PAP, PGELAM, PGEMU, PCLON, PSLON, &
!                                    & ZSO4BULK, ZBCBULK, ZOMBULK, ZSSBULK, ZDUBULK)
!        END IF
!       
!        CALL AER_ACTIV_MENON (KIDIA, KFDIA, KLON, ITOP, KLEV, PT, ZRHO, PLSM, ZSO4BULK, ZSSBULK, ZOMBULK, &
!                            & PGFL(:,:,YCDNC%MP9_PH))
!
!
!     CASE(JP_ACT_FOUNTOUKIS_NENES)
!        
!        IF ((.NOT.LTM5AER) .AND. (.NOT. LCMIP6_PI_AEROSOLS)) &
!         & CALL ABOR1('YOE_AER_ACTIV:Fountoukis and Nenes scheme requires LTM5AER=T or LMIP6_PI_AEROSOLS=T')
!
!        CALL PDF_UPDRAFT(KIDIA, KFDIA, KLON, ITOP, KLEV, NACTPDF, ZRHO, PLSM, PVERVEL, ZW, ZWPDF)
!
!        CALL AER_ACTIV_FOUNTOUKIS_NENES(KIDIA, KFDIA, ITOP, KLON, KLEV, NACTPDF, LLIQCLD, PT, PAP, ZRHO, & 
!                                      & PQ,  ZW, ZWPDF, ZSO4MASS, ZBCMASS, ZOMMASS, ZSSMASS, &
!                                      & ZDUMASS, ZNO3MASS, ZMSAMASS, ZAERONUM, ZDRYRSOL, &
!                                      & PGFL(:,:,YCDNC%MP9_PH), KFLDX, PEXTRA)
!
!     CASE(JP_ACT_MORALES_NENES)
!
!        IF ((.NOT.LTM5AER) .AND. (.NOT. LCMIP6_PI_AEROSOLS)) &
!         & CALL ABOR1('YOE_AER_ACTIV:Morales and Nenes scheme requires LTM5AER=T or LMIP6_PI_AEROSOLS=T')
!
!        CALL PDF_UPDRAFT(KIDIA, KFDIA, KLON, ITOP, KLEV, NACTPDF, ZRHO, PLSM, PVERVEL, ZW, ZWPDF)
!
!        CALL AER_ACTIV_MORALES_NENES(KIDIA, KFDIA, ITOP, KLON, KLEV, NACTPDF, LLIQCLD, PT, PAP, ZRHO, &
!                                      & PQ,  ZW, ZWPDF, ZSO4MASS, ZBCMASS, ZOMMASS, ZSSMASS, &
!                                      & ZDUMASS, ZNO3MASS, ZMSAMASS, ZAERONUM, ZDRYRSOL, &
!                                      & PGFL(:,:,YCDNC%MP9_PH), KFLDX, PEXTRA)
!
!     CASE(JP_ACT_MORALES_NENES_FULL)
!
!        IF ((.NOT.LTM5AER) .AND. (.NOT. LCMIP6_PI_AEROSOLS)) &
!         & CALL ABOR1('YOE_AER_ACTIV:Morales and Nenes full scheme requires LTM5AER=T or LMIP6_PI_AEROSOLS=T')

        ! The new scheme either approximates the integral over the updraft velocity PDF
        ! by Gauss-Legendre quadrature or uses a PPDEFsingle characteristic velocity.
        ! PDF_UPDRAFT, NACTPDF, ZW and ZWPDF are not used.
      CALL AER_ACTIV_MORALES_NENES_FULL(KIDIA, KFDIA, ITOP, KLON, KLEV, LLIQCLDD, PT, PAP, ZRHO, &
                                      & PVERVEL, ZSO4MASS, ZBCMASS, ZOMMASS, ZSSMASS, &
                                      & ZDUMASS, ZNO3MASS, ZMSAMASS, ZAERONUM, ZDRYRSOL, &
                                      & ZCDNC, ZSMAX, PGEMU, PSIGMA_W) !PSLON, PGEMU)                                    
                                      !& PGFL(:,:,YCDNC%MP9_PH), ZSMAX, PGEMU, PSIGMA_W) !PSLON, PGEMU)
                                      !& PGFL(:,:,YCDNC%MP9_PH), KFLDX, PEXTRA, PSLON, PGEMU)
!
!     CASE(JP_ACT_MORALES_NENES_FULL_OLDPDF)
!
!        IF ((.NOT.LTM5AER) .AND. (.NOT. LCMIP6_PI_AEROSOLS)) &
!         & CALL ABOR1('YOE_AER_ACTIV:Morales and Nenes full scheme with old PDF sampling requires LTM5AER=T or LMIP6_PI_AEROSOLS=T')
!
!        CALL PDF_UPDRAFT(KIDIA, KFDIA, KLON, ITOP, KLEV, NACTPDF, ZRHO, PLSM, PVERVEL, ZW, ZWPDF)
!
!        CALL AER_ACTIV_MORALES_NENES_FULL_OLDPDF(KIDIA, KFDIA, ITOP, KLON, KLEV, NACTPDF, LLIQCLD, PT, PAP, ZRHO, &
!                                      & ZW, ZWPDF, ZSO4MASS, ZBCMASS, ZOMMASS, ZSSMASS, &
!                                      & ZDUMASS, ZNO3MASS, ZMSAMASS, ZAERONUM, ZDRYRSOL, &
!                                      & PGFL(:,:,YCDNC%MP9_PH), KFLDX, PEXTRA, PSLON, PGEMU)
!
!     CASE(JP_ACT_ABDULRAZZAK_GHAN)
!
!        IF ((.NOT.LTM5AER) .AND. (.NOT. LCMIP6_PI_AEROSOLS)) &
!          & CALL ABOR1('YOE_AER_ACTIV:Abdul-Razzak and Ghan scheme requires LTM5AER=T or LCMIP6_PI_AEROSOLS=T')
!
!        CALL PDF_UPDRAFT(KIDIA, KFDIA, KLON, ITOP, KLEV, NACTPDF, ZRHO, PLSM, PVERVEL, ZW, ZWPDF)
!
!        CALL AER_ACTIV_ABDULRAZZAK_GHAN(KIDIA, KFDIA, ITOP, KLON, KLEV, NACTPDF, LLIQCLD, PT, PAP, ZRHO, & 
!                                      & PQ,  ZW, ZWPDF, ZSO4MASS, ZBCMASS, ZOMMASS, ZSSMASS, &
!                                      & ZDUMASS, ZNO3MASS, ZMSAMASS, ZAERONUM, ZDRYRSOL, &
!                                      & PGFL(:,:,YCDNC%MP9_PH), KFLDX, PEXTRA)           
!     END SELECT

      !---get CDNC_FACTOR from MAC2SP if needed
      ! IF (LMAC2SPACI .AND. .NOT.LTM5AER) THEN
      !    CALL GET_CDNC_FACTOR(KIDIA,KFDIA,KLON,PGEMU,PGELAM,ZMAC2SP_CDNC_FACTOR)
      ! ELSE
      !    ZMAC2SP_CDNC_FACTOR(:)=1._JPRB
      ! END IF

      ZMAC2SP_CDNC_FACTOR(:)=1._JPRB

      !---limit CDNC to min PPMINCDNC, set default value for CDNC outside clouds
      DO JK=KTDIA,KLEV
         ZCDNC(KIDIA:KFDIA,JK)=MAX(ZCDNC(KIDIA:KFDIA,JK)*ZMAC2SP_CDNC_FACTOR(KIDIA:KFDIA),PPMINCDNC)
         !ZCDNC(KIDIA:KFDIA,JK)=MERGE( &
         !& MAX(ZCDNC(KIDIA:KFDIA,JK)*ZMAC2SP_CDNC_FACTOR(KIDIA:KFDIA),PPMINCDNC), &
         !& PPDEFCDNC, LLIQCLD(KIDIA:KFDIA,JK) )
      !   PGFL(KIDIA:KFDIA,JK,YCDNC%MP9_PH)=MERGE( &
      !   & MAX(PGFL(KIDIA:KFDIA,JK,YCDNC%MP9_PH)*ZMAC2SP_CDNC_FACTOR(KIDIA:KFDIA),PPMINCDNC), &
      !   & PPDEFCDNC, LLIQCLD(KIDIA:KFDIA,JK) )
      END DO

      !---cloud liquid water: droplet effective radius is computed in radlswr now
      !CALL LIQ_CLOUD_RE(KIDIA, KFDIA, KLON, ITOP, KLEV, LLIQCLD, PL, PA, ZRHO, PGFL)

      ! liquid effective radius                                                                                                                         
      DO JK=1,KLEV
         DO JL=KIDIA,KFDIA                                                                
            ZQLWC(JL,JK) = PL(JL,JK)/MAX(PA(JL,JK),1.E-10_JPRB) ! calculate lwc                                                                         
            ZQLWC(JL,JK) = MIN(MAX(ZQLWC(JL,JK),0.0_JPRB),RCLDMAX) ! treshold lwc                                                                         
            ! effective radius calculated similarly as in radlswr.F90                                                                                     
            ! 2.387e-10 is 3/(4*pi*rho_liq*10^6)  [10^6 for N in right units]                                                                             
            !ZRE_LIQ(JL,JK) = 1.E+06_JPRB*(2.387e-10_JPRB*ZRHO(JL,JK)*ZQLWC(JL,JK)/(MAX(PGFL(JL,JK,YCDNC%MP9_PH),PPMINCDNC)))**0.333_JPRB ! calculate effective radius in um (use minimum value for CDNC if CDNC is small)
            ZRE_LIQ(JL,JK) = 1.E+06_JPRB*(2.387e-10_JPRB*ZRHO(JL,JK)*ZQLWC(JL,JK)/(MAX(ZCDNC(JL,JK),PPMINCDNC)))**0.333_JPRB ! calculate effective radius in um (use minimum value for CDNC if CDNC is small)
         END DO
      END DO
      ZRE_LIQ(KIDIA:KFDIA,1:KLEV) = MERGE(ZRE_LIQ(KIDIA:KFDIA,1:KLEV),4._JPRB,LLIQCLD(KIDIA:KFDIA,1:KLEV)) 
   END IF

   !---cloud ice: ICNC and effective radius for ice crystals 
   IF (LAERICESED .OR. LAERICEAUTO) THEN
      CALL ICE_CLOUD_PROP(KIDIA, KFDIA, KLON, KLEV, PT, ZRHO, PI, PA, PAP, &
                       &  PQSAT, ZSO4BULK, ZBCBULK, ZDUBULK, PGFL, YDMODEL, ZRE_ICE, ZICNC)        
   END IF
      
   !eehol: diagnostics:
   !--CDNC
   !PCDNCACT(KIDIA:KFDIA,:) = 1.0E6_JPRB*PGFL(KIDIA:KFDIA,:,YCDNC%MP9_PH) !eehol: output CDNC [#/m3]
   PCDNCACT(KIDIA:KFDIA,1:KLEV) = 1.0E6_JPRB*MAX(ZCDNC(KIDIA:KFDIA,1:KLEV),PPMINCDNC) !eehol: output CDNC [#/m3]

   !--ICNC
   !PICNC(KIDIA:KFDIA,:) = PGFL(KIDIA:KFDIA,:,YICNC%MP9_PH) !eehol: output ICNC [#/cm3]
   PICNC(KIDIA:KFDIA,1:KLEV) = ZICNC(KIDIA:KFDIA,1:KLEV) !eehol: output ICNC [#/cm3]
   
   !--Liq eff rad
   !PGFL(KIDIA:KFDIA,1:KLEV,YRE_LIQ%MP9_PH) = 1.0E-06_JPRB * ZRE_LIQ(KIDIA:KFDIA,1:KLEV) !save liq eff rad in meters
   PREFFL(KIDIA:KFDIA,1:KLEV) = MERGE(ZRE_LIQ(KIDIA:KFDIA,1:KLEV),4._JPRB,LLIQCLD(KIDIA:KFDIA,1:KLEV)) !eehol: output liq eff rad [um]

   !--Ice eff rad
   ! only if there is ice cloud else minimum value
   PREFFI(KIDIA:KFDIA,1:KLEV) = MERGE(ZRE_ICE(KIDIA:KFDIA,1:KLEV),20._JPRB,LICECLD(KIDIA:KFDIA,1:KLEV))
   !PREFFI(KIDIA:KFDIA,:) = 1.E-6_JPRB*PGFL(KIDIA:KFDIA,:,YRE_ICE%MP9_PH) !eehol: output ice eff rad [um]
   ! only if there is ice cloud else minimum value
   !PREFFI(KIDIA:KFDIA,1:KLEV) = MERGE(PREFFI(KIDIA:KFDIA,1:KLEV), 20._JPRB, LICECLD(KIDIA:KFDIA,1:KLEV))
   ! add effective radii to PGFL fields
   !PGFL(KIDIA:KFDIA,1:KLEV,YRE_ICE%MP9_PH) = 1.0E-06_JPRB * PREFFI(KIDIA:KFDIA,1:KLEV) ! convert um to meters and save to PGFL fields

   !--Maximum supersaturation
   PSMAX(KIDIA:KFDIA,1:KLEV) = ZSMAX(KIDIA:KFDIA,1:KLEV) !eehol: output maximum supersaturation [%]
   
!  IF (LACI_DIAG) THEN
!     !---2D diagnostics: AOD - removed. Can be re-implemented in RADLSWR if needed.
!
!     !---3D diagnostics: accumulated CDNC, ICNC, RE_liq, Re_ice, 
!     !   liquid cloud time, ice cloud time (ice diagnostics disabled for now)
!     !DO JK=1,KLEV
!     !   DO JL=KIDIA,KFDIA
!     !      LICECLD(JL,JK) = (PI(JL,JK) > RLMIN)
!     !   END DO
!     !END DO
!
!     IF (MOD(NSTEP,NFRPOS) == 0) THEN    ! if first timestep after output step
!        PEXTRA(KIDIA:KFDIA,:,D_CDNC%IXTRA) = 0._JPRB           
!        PEXTRA(KIDIA:KFDIA,:,D_LIQCLDT%IXTRA) = 0._JPRB           
!     END IF
!
!     !--Reff (liq)
!     !  moved to RADLSWR
!
!     !--CDNC
!     ZTMP(KIDIA:KFDIA,:) = PEXTRA(KIDIA:KFDIA,:,D_CDNC%IXTRA) + PGFL(KIDIA:KFDIA,:,YCDNC%MP9_PH)*TSTEP
!     PEXTRA(KIDIA:KFDIA,:,D_CDNC%IXTRA) = MERGE(ZTMP(KIDIA:KFDIA,:), &
!                                & PEXTRA(KIDIA:KFDIA,:,D_CDNC%IXTRA), LLIQCLD(KIDIA:KFDIA,:))
!     !--cloud time (liq)
!     ZTMP(KIDIA:KFDIA,:) = PEXTRA(KIDIA:KFDIA,:,D_LIQCLDT%IXTRA) + TSTEP
!     PEXTRA(KIDIA:KFDIA,:,D_LIQCLDT%IXTRA) = MERGE(ZTMP(KIDIA:KFDIA,:), &
!                                & PEXTRA(KIDIA:KFDIA,:,D_LIQCLDT%IXTRA), LLIQCLD(KIDIA:KFDIA,:))
!
!     !--ICNC (not in CMIP6)
!     ! IF (LAERICESED .OR. LAERICEAUTO) THEN
!     !    ZTMP(KIDIA:KFDIA,:) = PEXTRA(KIDIA:KFDIA,:,IX+3) + PGFL(KIDIA:KFDIA,:,YICNC%MP9_PH)*TSTEP
!     !   PEXTRA(KIDIA:KFDIA,:,IX+3) = MERGE(ZTMP(KIDIA:KFDIA,:), &
!     !                              & PEXTRA(KIDIA:KFDIA,:,IX+3), LICECLD(KIDIA:KFDIA,:))
!     ! END IF
!
!     !--Reff (ice) (not in CMIP6)
!     ! IF (LAERICESED .OR. LAERICEAUTO) THEN
!     !    ZTMP(KIDIA:KFDIA,:) = PEXTRA(KIDIA:KFDIA,:,IX+4) + PGFL(KIDIA:KFDIA,:,YRE_ICE%MP9_PH)*TSTEP
!     !    PEXTRA(KIDIA:KFDIA,:,IX+4) = MERGE(ZTMP(KIDIA:KFDIA,:), &
!     !                            & PEXTRA(KIDIA:KFDIA,:,IX+4), LICECLD(KIDIA:KFDIA,:))
!     ! END IF
!
!     !--cloud time (ice) (not in CMIP6)
!     ! IF (LAERICESED .OR. LAERICEAUTO) THEN
!     !    ZTMP(KIDIA:KFDIA,:) = PEXTRA(KIDIA:KFDIA,:,IX+5) + TSTEP
!     !    PEXTRA(KIDIA:KFDIA,:,IX+5) = MERGE(ZTMP(KIDIA:KFDIA,:), &
!     !                            & PEXTRA(KIDIA:KFDIA,:,IX+5), LICECLD(KIDIA:KFDIA,:))
!     ! END IF
!
!  END IF

   END ASSOCIATE
   END ASSOCIATE

   IF (LHOOK) CALL DR_HOOK('YOE_AER_ACTIV.AER_ACTIV', 1, ZHOOK_HANDLE)

  END SUBROUTINE AER_ACTIV


!  SUBROUTINE AER_ACTIV_FOUNTOUKIS_NENES(KIDIA, KFDIA, KTDIA, KLON, KLEV, KPDF, LCLOUD, PT, PAP, PRHO, & 
!       & PQ,  PW, PWPDF, PSO4MASS, PBCMASS, POMMASS, PSSMASS, &
!       & PDUMASS,PNO3MASS, PMSAMASS, PAERONUM, PRDRY, PCDNC, KFLDX, PEXTRA)
!
!
!    ! *aer_activ_fountoukis_nenes* calculates the number of activated aerosol 
!    !              particles from the aerosol size-distribution,
!    !              composition and ambient supersaturation
!    !
!    ! Author:
!    ! -------
!    ! Sami Romakkaniemi, FMI
!    ! Philip Stier, University of Oxford
!    ! Declan O'Donnell,  FMI
!    ! Twan van Noije, KNMI
!    !
!    ! References:
!    ! -----------
!    ! Abdul-Razzak et al., JGR, 103, D6, 6123-6131, 1998
!    ! Abdul-Razzak and Ghan, JGR, 105, D5, 6837-6844, 2000 (referred to as ARG)
!    ! Pruppbacher and Klett, Kluewer Ac. Pub., 1997
!    ! Ghan et al., JAMES 3, M10001, 2011
!    ! Nenes and Seinfeld, JGR, 108, D14, 4415, 2003 (referred to as NS)
!    ! Fountoukis and Nenes, JGR, 110, D11212, 2005 (referred to as FN)
!    ! Morales Betancourt and Nenes, GMD, 7, 2345-2357, 2014 (referred to as MN)
!    ! Seinfeld and Pandis, Atmospheric Chemistry and Physics, Second Edition (referred to as SP)
!    
!
!    USE YOMCST,              ONLY: R, RV, RPI, RCPD, RG, RLVTT, RMV, RMD, RTT , RLSTT
!    USE YOECLDP,             ONLY: RTHOMO, PPRHO_WAT
!    USE YOMHOOK,             ONLY: LHOOK, DR_HOOK
!    USE YOMLUN,              ONLY: NULOUT
!    USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
!         & R4IES    ,R5LES    ,R5IES    ,R5ALVCP  ,R5ALSCP  ,&
!         & RALVDCP  ,RALSDCP  ,RTWAT    ,&
!         & RTICE    ,RTICECU  ,&
!         & RTWAT_RTICE_R      ,RTWAT_RTICECU_R,&
!         & RKOOP1   ,RKOOP2
!    USE YOE_AERO_M7_DATA,    ONLY: NMOD, NSOL, SIGMALN, CMR2RAM, &
!         & DH2SO4, DBC, DOC, DNACL, DDUST, &
!         & DNA2SO4, DNH4NO3, DMSA, NH4NO3_FACTOR, &
!         & PPKAPPA_H2SO4, PPKAPPA_NACL, PPKAPPA_NA2SO4, &
!         & PPKAPPA_BC, PPKAPPA_OC, PPKAPPA_DU, &
!         & PPKAPPA_NH4NO3, PPKAPPA_MSA, &
!         & WSO4, WH2SO4, WNACL, WNA2SO4, &
!         & WH2O, WDAIR
!
!    IMPLICIT NONE
!
!    !---included functions from header files
!#include "fcttre.h"
!
!    !---subroutine interface
!    INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA 
!    INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA 
!    INTEGER(KIND=JPIM), INTENT(IN) :: KTDIA 
!    INTEGER(KIND=JPIM), INTENT(IN) :: KLON 
!    INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
!    INTEGER(KIND=JPIM), INTENT(IN) :: KPDF
!    INTEGER(KIND=JPIM), INTENT(IN) :: KFLDX
!
!    LOGICAL, INTENT(IN)            :: LCLOUD(KLON,KLEV)
!
!    REAL(KIND=JPRB), INTENT(IN)    :: PT(KLON,KLEV)
!    REAL(KIND=JPRB), INTENT(IN)    :: PAP(KLON,KLEV)
!    REAL(KIND=JPRB), INTENT(IN)    :: PRHO(KLON,KLEV)
!    REAL(KIND=JPRB), INTENT(IN)    :: PQ(KLON,KLEV)
!    REAL(KIND=JPRB), INTENT(IN)    :: PW(KLON,KLEV,KPDF)
!    REAL(KIND=JPRB), INTENT(IN)    :: PWPDF(KLON,KLEV,KPDF)
!    REAL(KIND=JPRB), INTENT(IN)    :: PSO4MASS(KLON,KLEV,NSOL) ! [KG(SO4)/KG(AIR)]
!    REAL(KIND=JPRB), INTENT(IN)    :: PBCMASS(KLON,KLEV,NSOL)  ! [KG(BC)/KG(AIR)]
!    REAL(KIND=JPRB), INTENT(IN)    :: POMMASS(KLON,KLEV,NSOL)  ! [KG(OM)/KG(AIR)]
!    REAL(KIND=JPRB), INTENT(IN)    :: PSSMASS(KLON,KLEV,NSOL)  ! [KG(SS)/KG(AIR)]
!    REAL(KIND=JPRB), INTENT(IN)    :: PDUMASS(KLON,KLEV,NSOL)  ! [KG(DU)/KG(AIR)]
!    REAL(KIND=JPRB), INTENT(IN)    :: PNO3MASS(KLON,KLEV)      ! [KG(NO3)/KG(AIR)]
!    REAL(KIND=JPRB), INTENT(IN)    :: PMSAMASS(KLON,KLEV)      ! [KG(MSA)/KG(AIR)]
!    REAL(KIND=JPRB), INTENT(IN)    :: PAERONUM(KLON,KLEV,NSOL)
!    REAL(KIND=JPRB), INTENT(IN)    :: PRDRY(KLON,KLEV,NSOL)
!    REAL(KIND=JPRB), INTENT(IN)    :: PEXTRA(KLON,KLEV,KFLDX)
!
!    !   Output:
!    REAL(KIND=JPRB), INTENT(INOUT) :: PCDNC(KLON,KLEV)
!
!    !---local data
!    INTEGER(KIND=JPIM) :: JL, JK, JMOD, JW
!
!    REAL(KIND=JPRB)    :: ZN(KLON,KLEV,NSOL)      ! aerosol number concentration for each mode [m-3]
!    REAL(KIND=JPRB)    :: ZSM(KLON,KLEV,NSOL)     ! critical supersaturation for activating particles
!                                                  ! with the mode number median radius
!    REAL(KIND=JPRB)    :: ZVOL(KLON)              ! total dry particle volume
!    REAL(KIND=JPRB)    :: ZKAPPA(KLON)            ! volume-weighted kappa
!    REAL(KIND=JPRB)    :: ZSMAX(KLON,KLEV,KPDF)   ! maximum supersaturation 
!    REAL(KIND=JPRB)    :: ZESW(KLON,KLEV)         ! saturation water vapour pressure
!    REAL(KIND=JPRB)    :: ZDIF(KLON,KLEV)         ! diffusivity
!    REAL(KIND=JPRB)    :: ZK(KLON,KLEV)           ! thermal conductivity
!    REAL(KIND=JPRB)    :: ZA(KLON,KLEV)           ! Kelvin coefficient
!    REAL(KIND=JPRB)    :: ZALPHA(KLON,KLEV)       ! Intermediate term in supersaturation calulation
!    REAL(KIND=JPRB)    :: ZGAMMA(KLON,KLEV)       ! Intermediate term in supersaturation calulation
!    REAL(KIND=JPRB)    :: ZGROWTH(KLON,KLEV)      ! Growth coefficient
!    REAL(KIND=JPRB)    :: ZTERM1(KLON)            ! Intermediate term in growth coefficient calulation
!    REAL(KIND=JPRB)    :: ZTERM2(KLON)            ! Intermediate term in growth coefficient calulation
!    REAL(KIND=JPRB)    :: ZTERM3(KLON)            ! Intermediate term in growth coefficient calulation
!
!    REAL(KIND=JPRB)    :: ZAMW                    ! molecular weight of water [kg mol-1]
!    REAL(KIND=JPRB)    :: ZAMD                    ! molecular weight of dry air [kg mol-1]
!
!    REAL(KIND=JPRB)    :: ZKA(KLON), ZKV(KLON)    ! Intermediate terms in supersaturation calulation
!    REAL(KIND=JPRB)    :: ZXV(KLON), ZB(KLON)     ! Intermediate terms in supersaturation calulation
!
!    REAL(KIND=JPRB)    :: ZDIFMOD(KLON)           ! Modified diffusivity
!    REAL(KIND=JPRB)    :: ZNACT                   ! Intermediate values of the activated number concentration (#/m3)
!    REAL(KIND=JPRB)    :: ZNACT_WSUM(KLON,KLEV)   ! Weighted sum of activated number concentration
!    REAL(KIND=JPRB)    :: ZPDF_NORM(KLON,KLEV)    ! Normalization factor for ZNACT_WSUM
!
!    REAL(KIND=JPRB)    :: ZSSMASS(KLON)           ! Sea salt MMR
!    REAL(KIND=JPRB)    :: ZDUMASS(KLON)           ! Dust MMR
!    REAL(KIND=JPRB)    :: ZNO3MASS(KLON)          ! Nitrate MMR
!    REAL(KIND=JPRB)    :: ZMSAMASS(KLON)          ! MSA MMR
!
!    REAL(KIND=JPRB)    :: NSO4(KLON), NH2SO4(KLON) ! Particle numbers [kmol/kg air]
!    REAL(KIND=JPRB)    :: NNACL(KLON), NNA(KLON), NCL(KLON), NNA2SO4(KLON) 
!
!    ! Per-mode constants
!    REAL(KIND=JPRB) :: ZMODECST1(NSOL), ZMODECST2(NSOL), ZMODECST3(NSOL), &
!              & ZMODECST4(NSOL), ZMODECST5(NSOL)
!
!    ! Intermediate values for the Fountoukis & Nenes iterative scheme
!    REAL(KIND=JPRB) :: ZCF1, ZCF2, ZCF3
!    REAL(KIND=JPRB) :: ZVALUE1
!    REAL(KIND=JPRB) :: ZVALUE2
!    REAL(KIND=JPRB) :: ZVALUE3
!    REAL(KIND=JPRB) :: ZINT1
!    REAL(KIND=JPRB) :: ZINT2
!    REAL(KIND=JPRB) :: ZSMAXTEMP1, ZSMAXTEMP2, ZSMAXTEMP3
!
!    ! Control of iterative loop:
!    ! According to Ghan et al. (2011) the FN scheme
!    ! takes about 30 interations to converge;
!    ! so safer to increase the maximum here.
!    ! xxx to be tested
!    INTEGER(KIND=JPIM), PARAMETER :: NMAXITER = 30
!    INTEGER(KIND=JPIM) :: NITERATIONS 
!    LOGICAL  :: LCONVERGED 
!
!    ! Miscellaneous
!    REAL(KIND=JPRB)    :: ZEPS
!    REAL(KIND=JPRB)    :: Z4PIOVER3, ZSQRT2
!    REAL(KIND=JPRB)    :: ZHOOK_HANDLE
!    REAL(KIND=JPRB), PARAMETER :: PPEPSSEC = 1.E-25_JPRB  ! used to avoid division by 0
!
!    ! mass accomodation coefficient
!    ! between 0.1 and 1.0
!    ! Raatikainen et al., 2013 PNAS
!    REAL(KIND=JPRB), PARAMETER :: PPALPHA_C = 1.E-1_JPRB
!  
!    ! Upper and lower size bounds (m) used for calculating the average water vapor diffusivity
!    REAL(KIND=JPRB), PARAMETER :: DPBIG = 5.E-6_JPRB
!    REAL(KIND=JPRB) :: DPLOW
!
!    !--- executable procedure
!    IF (LHOOK) CALL DR_HOOK('YOE_AER_ACTIV.AER_ACTIV_FOUNTOUKIS_NENES',0,ZHOOK_HANDLE)
!
!    !--- 0) Initializations:
!
!    ZSMAX(KIDIA:KFDIA,KTDIA:KLEV,:) = 0._JPRB
!    ZSM(KIDIA:KFDIA,KTDIA:KLEV,:) = 0._JPRB
!    ZNACT_WSUM(KIDIA:KFDIA,KTDIA:KLEV) = 0._JPRB
!    ZPDF_NORM(KIDIA:KFDIA,KTDIA:KLEV) = 0._JPRB
!    PCDNC(KIDIA:KFDIA,KTDIA:KLEV) = 0._JPRB
!
!    ZEPS=EPSILON(1._JPRB)
!
!    !--- Conversions to SI units [g mol-1 to kg mol-1]:    
!    ZAMW=WH2O*1.E-3_JPRB
!    ZAMD=WDAIR*1.E-3_JPRB
!
!    !---miscellaneous
!    Z4PIOVER3 = 4._JPRB*RPI/3._JPRB
!    ZSQRT2 = SQRT(2._JPRB)
!    
!    ! FN, Eq. (24) converted to m.
!    DPLOW = MIN(0.207683E-6_JPRB * PPALPHA_C**(-0.33048_JPRB), DPBIG)
!
!    !---grid-point calculations
!    DO JK=KTDIA,KLEV
!       DO JL=KIDIA,KFDIA
!       IF (LCLOUD(JL,JK)) THEN
!          !---Kelvin (curvature) coefficient, usually denoted by capital A.
!          !   Here we use the definition of SP (introduced on p. 770),
!          !   which is adopted in NS, FN, MN,
!          !   and include an extra factor of 2 compared to the definition of Abdul-Razzak et al. (Eq. 5).
!          ZA(JL,JK) = 4._JPRB * ZAMW * PPSURFTEN / (R * PPRHO_WAT * PT(JL,JK))
!
!          !--- Abdul-Razzak et al. (1998) (Eq. 11):
!          ZALPHA(JL,JK) = (RG*ZAMW*RLVTT) / (RCPD*R*PT(JL,JK)*PT(JL,JK)) - &
!                        & (RG*ZAMD) / (R*PT(JL,JK))
!
!          ! Saturation water vapour pressure:
!          ZESW(JL,JK) = FOEEWM(PT(JL,JK))
!
!          !--- Following the definitions of NS, FN and MN,
!          !    gamma is defined as a dimensionless coefficient.
!          !    We use Eq. (12) from Abdul-Razzak et al. (1998) and multiply it with the air density:
!          ZGAMMA(JL,JK) = ( (R*PT(JL,JK)) / (ZESW(JL,JK)*ZAMW) +  &
!                        &   (ZAMW*RLVTT*RLVTT) / (RCPD*PAP(JL,JK)*ZAMD*PT(JL,JK)) ) * PRHO(JL,JK)
!
!          !--- Diffusivity of water vapour in air (P&K, 13.3) [m2 s-1]:
!
!          ZDIF(JL,JK)=0.211_JPRB * (PT(JL,JK)/RTT)**1.94_JPRB * (101325._JPRB/PAP(JL,JK)) *1.E-4_JPRB
!
!          !--- modified diffusivity 
!          !    Average mode-independent value using FN, Eq. (23)
!          ZB(JL) = (2._JPRB*ZDIF(JL,JK)/PPALPHA_C)*SQRT(2._JPRB*RPI*ZAMW/(R*PT(JL,JK)))
!
!          !--- For any reasonable value of PPALPHA_C: DPBIG > DPLOW
!          ZDIFMOD(JL) = ZDIF(JL,JK)*(1._JPRB-(ZB(JL)/(DPBIG-DPLOW))*LOG((DPBIG+ZB(JL))/(DPLOW+ZB(JL))))
!
!          !--- Thermal conductivity zk (P&K, 13.18) [cal cm-1 s-1 K-1]:
!
!          ! Mole fraction of water:
!
!          ZXV(JL) = PQ(JL,JK)*(ZAMD/ZAMW)
!
!          ZKA(JL) = (5.69_JPRB+0.017_JPRB*(PT(JL,JK)-273.15_JPRB))*1.E-5_JPRB
!
!          ZKV(JL) = (3.78_JPRB+0.020_JPRB*(PT(JL,JK)-273.15_JPRB))*1.E-5_JPRB
!
!          ! Moist air, convert to [J m-1 s-1 K-1]:
!
!          ZK(JL,JK) = ZKA(JL)*(1._JPRB-(1.17_JPRB-1.02_JPRB*ZKV(JL)/ZKA(JL))*ZXV(JL)) &
!                    & * 4.1868_JPRB*1.E2_JPRB
!
!          !--- growth coefficient due to gas kinetic effects:
!          !--- NS, Eq. (15)
!
!          ZTERM1(JL) = (PPRHO_WAT*R*PT(JL,JK)) / (ZESW(JL,JK)*ZDIFMOD(JL)*ZAMW)
!
!          !--- Note that no size dependence is introduced in the thermal conductivity  
!          !    See FN, p. 5
!          ZTERM2(JL) = (RLVTT*PPRHO_WAT) / (ZK(JL,JK)*PT(JL,JK))
!
!          ZTERM3(JL) = (RLVTT*ZAMW) / (R*PT(JL,JK))-1._JPRB
!
!          !--- Note that the expression for G in NS (Eq. 15) has an additional factor of 4
!          !    compared to Abdul-Razzak et al. (1998), Eq. (16).
!          !    FN follow the definition from NS,
!          !    but have omitted the factor of 4 from their Eq. (12). This is a typo.
!          !    See also SP (Eq. 17.70) and MN (Eq. A3).
!          ZGROWTH(JL,JK) = 4._JPRB / (ZTERM1(JL) + ZTERM2(JL) * ZTERM3(JL))
!
!       END IF
!       END DO
!    END DO
!
!    !---per-mode calculations:
!    !   soluble mode number and critical supersaturation, ignore nucleation mode (mode 1)
!    DO JMOD=2, NSOL
!       DO JK=KTDIA,KLEV
!          DO JL=KIDIA,KFDIA
!          IF (LCLOUD(JL,JK)) THEN
!          
!             !---total volume per mode [m^3 / kg(air)], used for kappa calculation
!             ZVOL(JL) = Z4PIOVER3 * PAERONUM(JL,JK,JMOD) *   &
!                      & (CMR2RAM(JMOD)*PRDRY(JL,JK,JMOD))**3 
!
!             !--- Number per unit volume [# m-3] for each mode:
!             ZN(JL,JK,JMOD) = PAERONUM(JL,JK,JMOD)*PRHO(JL,JK)
!
!             !---sea salt and dust do not exist in mode 2:
!             ZSSMASS(JL) = MERGE(0._JPRB, PSSMASS(JL,JK,JMOD), JMOD==2)
!             ZDUMASS(JL) = MERGE(0._JPRB, PDUMASS(JL,JK,JMOD), JMOD==2)
!
!             !---ammonium-nitrate and MSA do not exit in mode 2 and 4:
!             ZNO3MASS(JL) = MERGE(0._JPRB, PNO3MASS(JL,JK), JMOD==2 .OR. JMOD==4)
!             ZMSAMASS(JL) = MERGE(0._JPRB, PMSAMASS(JL,JK), JMOD==2 .OR. JMOD==4)
!
!             NNA(JL) = ZSSMASS(JL) / WNACL
!             NCL(JL) = NNA(JL)
!             NSO4(JL) = PSO4MASS(JL,JK,JMOD) / WSO4
!             NNA2SO4(JL) = MIN(NNA(JL)/2._JPRB, NSO4(JL))
!             NNA(JL) = NNA(JL) - 2._JPRB*NNA2SO4(JL)
!             NNACL(JL) = MIN(NCL(JL), NNA(JL))
!             NCL(JL) = NNACL(JL)
!             NH2SO4(JL) = NSO4(JL) - NNA2SO4(JL)
!
!             !---mode kappa = volume-weighted sum of component kappa's
!             ZKAPPA(JL) = ( (PPKAPPA_NACL * NNACL(JL) * WNACL / (DNACL*1.E3_JPRB)) + &
!                  & (PPKAPPA_NA2SO4 * NNA2SO4(JL) * WNA2SO4 / (DNA2SO4*1.E3_JPRB)) + &
!                  & (PPKAPPA_H2SO4 * NH2SO4(JL) * WH2SO4 / (DH2SO4*1.E3_JPRB)) + &
!                  & (PPKAPPA_BC * PBCMASS(JL,JK,JMOD) / (DBC*1.E3_JPRB))       + &
!                  & (PPKAPPA_OC * POMMASS(JL,JK,JMOD) / (DOC*1.E3_JPRB))       + &
!                  & (PPKAPPA_DU * ZDUMASS(JL) / (DDUST*1.E3_JPRB)) + &
!                  & (PPKAPPA_NH4NO3 * ZNO3MASS(JL) * NH4NO3_FACTOR / (DNH4NO3*1.E3_JPRB)) + &
!                  & (PPKAPPA_MSA * ZMSAMASS(JL) / (DMSA*1.E3_JPRB)) )   / &
!                  & ZVOL(JL)
!
!             !---defensive step: minimum kappa to avoid divide by zero errors
!             ZKAPPA(JL) = MERGE(ZKAPPA(JL), 0.04_JPRB, ZKAPPA(JL) > 0.04_JPRB )
!
!             !---eqn. (2) from Ghan et al (2011)
!             !   an addition factor 0.5**3 has been included in the nominator
!             !   because the Kelvin coefficient used here (ZA) 
!             !   is twice that in Ghan et al. (ZKELV)
!             ZSM(JL,JK,JMOD) = SQRT(0.5_JPRB * ZA(JL,JK)**3 / &
!                                 & (27._JPRB * ZKAPPA(JL) * PRDRY(JL,JK,JMOD)**3) )
!
!          END IF
!          END DO
!       END DO
!    END DO
!
!
!    !--- 2) Calculate maximum supersaturation:
!
!    DO JMOD=2,NSOL
!
!       !---some per-mode constants
!       ZMODECST1(JMOD) = 3.0_JPRB*SIGMALN(JMOD)/ZSQRT2      ! final term in FN, Eq. (18)
!       ZMODECST2(JMOD) = EXP(4.5_JPRB*SIGMALN(JMOD)**2)     ! factor in FN, Eq. (18)
!       
!       ZMODECST3(JMOD) = ZMODECST1(JMOD)/2._JPRB            ! final term in FN, Eq. (19)
!       ZMODECST4(JMOD) = EXP(1.125_JPRB*SIGMALN(JMOD)**2)   ! factor in FN, Eq. (19)
!       
!       ZMODECST5(JMOD) = 3.0_JPRB*ZSQRT2*SIGMALN(JMOD)      ! denominator in FN, Eqs. (8) and (20)
!
!    END DO
!
!    DO JW=1,KPDF
!       DO JK=KTDIA, KLEV
!          DO JL=KIDIA,KFDIA
!             IF (LCLOUD(JL,JK)) THEN
!                IF (ZALPHA(JL,JK) > ZEPS .AND. PW(JL,JK,JW) > ZEPS) THEN
!                   ZSMAXTEMP1 = 1.0E-5_JPRB ! min cloud supersat.
!                   ZSMAXTEMP2 = 0.1_JPRB    ! max cloud supersat.
!                   LCONVERGED = .FALSE.
!                   NITERATIONS = 0
!
!                   ZCF1 = 0.5_JPRB*SQRT(ZGROWTH(JL,JK)/(ZALPHA(JL,JK)*PW(JL,JK,JW)))
!                   ZCF2 = ZA(JL,JK)/3.0_JPRB
!                   ! When gamma is defined as a dimensionless coefficient as above,
!                   ! ZCF3 should include a division by the air density.
!                   ! See MN, Eq. 5.
!                   ! This is an error in NS (Eq. 32) and FN (Eq. 10).
!                   ! It is can be traced back to the definition of the rate of water condensation in NS, Eq. (11),
!                   ! where W is expressed in units kg(liquid water)/m3(air) instead of kg(liquid water)/kg(air).
!                   ! See also SP, Eqs. (17.73) and (17.79).
!                   ZCF3 = 0.5_JPRB*RPI*ZGAMMA(JL,JK)*PPRHO_WAT*ZGROWTH(JL,JK) &
!                        & /ZALPHA(JL,JK)/PW(JL,JK,JW)/PRHO(JL,JK)
!
!                   CALL SINTEGRAL(ZA(JL,JK),    ZALPHA(JL,JK), ZGAMMA(JL,JK), ZGROWTH(JL,JK), &
!                               &  PW(JL,JK,JW), ZN(JL,JK,:),   ZSM(JL,JK,:),  ZSMAXTEMP1,     &
!                               &  ZMODECST1(:), ZMODECST2(:),  ZMODECST3(:),  ZMODECST4(:),   &
!                               &  ZMODECST5(:), ZINT1,         ZINT2)
!                      
!                   ZVALUE1 = (ZINT1*ZCF1 + ZINT2*ZCF2)*ZCF3*ZSMAXTEMP1 - 1.0_JPRB
!                         
!                   !DIR$ INLINE
!                   CALL SINTEGRAL(ZA(JL,JK),    ZALPHA(JL,JK), ZGAMMA(JL,JK), ZGROWTH(JL,JK), &
!                               &  PW(JL,JK,JW), ZN(JL,JK,:),   ZSM(JL,JK,:),  ZSMAXTEMP2,     &
!                               &  ZMODECST1(:), ZMODECST2(:),  ZMODECST3(:),  ZMODECST4(:),   &
!                               &  ZMODECST5(:), ZINT1,         ZINT2)
!                   ZVALUE2 = (ZINT1*ZCF1 + ZINT2*ZCF2) * ZCF3*ZSMAXTEMP2 - 1.0_JPRB
!                
!                   DO WHILE ((.NOT.LCONVERGED) .AND. NITERATIONS < NMAXITER)
!                      NITERATIONS = NITERATIONS + 1 
!                
!                      ZSMAXTEMP3 = 0.5_JPRB*(ZSMAXTEMP1+ZSMAXTEMP2)
!                      !DIR$ INLINE
!                      CALL SINTEGRAL(ZA(JL,JK),    ZALPHA(JL,JK), ZGAMMA(JL,JK), ZGROWTH(JL,JK), &
!                                  &  PW(JL,JK,JW), ZN(JL,JK,:),   ZSM(JL,JK,:),  ZSMAXTEMP3,     &
!                                  &  ZMODECST1(:), ZMODECST2(:),  ZMODECST3(:),  ZMODECST4(:),   &
!                                  &  ZMODECST5(:), ZINT1,         ZINT2)
!                      ZVALUE3 = (ZINT1*ZCF1 + ZINT2*ZCF2) * ZCF3*ZSMAXTEMP3 - 1.0_JPRB
!
!                      IF (SIGN(1.0_JPRB,ZVALUE1)*SIGN(1.0_JPRB,ZVALUE3) <= 0.0_JPRB) THEN
!                         ZVALUE2 = ZVALUE3
!                         ZSMAXTEMP2 = ZSMAXTEMP3
!                      ELSE
!                         ZVALUE1 = ZVALUE3
!                         ZSMAXTEMP1 = ZSMAXTEMP3
!                      END IF
!
!                      !xxx 1e-3, i.e. 0.1%, should be enough
!                      IF (ABS(ZSMAXTEMP2-ZSMAXTEMP1) <= 1.0E-5_JPRB*ZSMAXTEMP1) THEN
!
!                         ZSMAX(JL,JK,JW) = 0.5_JPRB*(ZSMAXTEMP1+ZSMAXTEMP2)
!                         LCONVERGED = .TRUE.  
!
!                         ZNACT = 0._JPRB 
!                         DO JMOD=2,NSOL 
!                            IF (ZSM(JL,JK,JMOD) > ZEPS) THEN
!                               ZNACT = ZNACT + 0.5_JPRB*ZN(JL,JK,JMOD)*ERFC(2._JPRB*(LOG(ZSM(JL,JK,JMOD)/ZSMAX(JL,JK,JW)))/ZMODECST5(JMOD))
!                            END IF
!                         END DO
!
!                         !---Sum up the total number of activated particles, weighted by the updraft PDF [m-3]:
!                         ZNACT_WSUM(JL,JK) = ZNACT_WSUM(JL,JK) + ZNACT*PWPDF(JL,JK,JW)
!
!                         ! Including the normalization here effectively means that conditions that don't yield a converged solution
!                         ! are discarded when calculating the mean CDNC over the updraft PDF.
!                         ZPDF_NORM(JL,JK) = ZPDF_NORM(JL,JK) + PWPDF(JL,JK,JW)
!
!                      END IF
!                      
!                   END DO    ! do while
!
!                   IF (.NOT.LCONVERGED) THEN
!                      WRITE(NULOUT,*) 'WARNING: Fountoukis and Nenes scheme not converged for updraft velocity', PW(JL,JK,JW)
!                   ENDIF
!
!                ELSE
!
!                   ! Set contribution to activation to zero when conditions are not met
!                   ZPDF_NORM(JL,JK) = ZPDF_NORM(JL,JK) + PWPDF(JL,JK,JW)
!
!                END IF          ! END IF safe wrt divide by 0, LOG(0)...
!             END IF 
!          END DO !jl
!       END DO !jk
!    END DO
!
!    DO JK=KTDIA,KLEV
!       DO JL=KIDIA,KFDIA
!          !---normalize the total number of activated particles over the PDF, and convert to [# cm-3]
!          IF (LCLOUD(JL,JK) .AND. ZPDF_NORM(JL,JK) > ZEPS ) THEN
!             PCDNC(JL,JK) = 1.E-6_JPRB * ZNACT_WSUM(JL,JK) / ZPDF_NORM(JL,JK)
!          END IF
!       END DO
!    END DO
!
!    IF (LHOOK) CALL DR_HOOK('YOE_AER_ACTIV.AER_ACTIV_FOUNTOUKIS_NENES',1,ZHOOK_HANDLE)
!
!  END SUBROUTINE AER_ACTIV_FOUNTOUKIS_NENES
!
!  SUBROUTINE SINTEGRAL(PA, PALPHA, PGAMMA, PGROWTH, PW, PN, PSM, &
!                       PSMAX, PUP1, PF1, PUP2, PF2, PUPD, PINTEG1, PINTEG2)
!
!    USE YOE_AERO_M7_DATA, ONLY: NSOL
!
!    IMPLICIT NONE
!
!    REAL(KIND=JPRB), INTENT(IN) :: PN(NSOL), PSM(NSOL)
!    REAL(KIND=JPRB), INTENT(IN) :: PUP1(NSOL), PF1(NSOL), PUP2(NSOL), PF2(NSOL), PUPD(NSOL)
!    REAL(KIND=JPRB), INTENT(IN) :: PA, PALPHA, PGAMMA, PGROWTH, PW, PSMAX  
!
!    REAL(KIND=JPRB), INTENT(OUT) :: PINTEG1, PINTEG2
!
!    REAL(KIND=JPRB) :: ZUPART, ZUPM, ZUMAX, ZUMM, ZUPP
!
!    REAL(KIND=JPRB) :: ZDELTA, ZRATIO, ZSPART2
!
!    REAL(KIND=JPRB) :: ZEPS
!
!    INTEGER(KIND=JPIM) :: JMOD
!
!    ! note: this subroutine is intended for inline compilation and therefore does not call dr_hook
!
!    ZEPS=EPSILON(1._JPRB)
!
!    ZDELTA = 1.0_JPRB - 16.0_JPRB/(9.0_JPRB*PGROWTH)*PALPHA*PW*(PA/PSMAX**2)**2
!
!    !    ZDELTA = ZSMAX**4 - 16.0_JPRB/(9.0_JPRB*ZGROWTH)*ZALPHA*ZW*ZA**2
!    IF (ZDELTA <= 0.0_JPRB) THEN
!       ZRATIO = (2.0E7_JPRB/3.0_JPRB)*PA*PSMAX**(-0.3824)
!       IF (ZRATIO > 1.0_JPRB) THEN
!          ZRATIO = 1.0_JPRB
!       END IF
!       ZSPART2 = PSMAX*ZRATIO 
!    ELSE
!       ZSPART2 = 0.5_JPRB*(1.0_JPRB + SQRT(ZDELTA)) !max root
!       ZSPART2 = SQRT(ZSPART2)*PSMAX
!    ENDIF
!
!    !** Calculate integrals
!
!    PINTEG1 = 0._JPRB
!    PINTEG2 = 0._JPRB
!    DO JMOD=2,NSOL
!       IF (PSM(JMOD) > ZEPS) THEN
!          ! FN, Eq. (20):
!          ZUPART        = 2.0_JPRB*LOG(PSM(JMOD)/ZSPART2)/PUPD(JMOD)
! 
!          ZUMAX         = 2.0_JPRB*LOG(PSM(JMOD)/PSMAX)/ PUPD(JMOD)
!  
!          ! argument to last erfc in FN, Eq. (18)
!          ZUPP          = ZUPART + PUP1(JMOD)
!
!          ! argument to first erf in FN, Eq. (19)
!          ZUPM          = ZUPART - PUP2(JMOD)
!
!          ! argument to second erf in FN, Eq. (19)
!          ZUMM          = ZUMAX - PUP2(JMOD)
!
!          ! intergral I1: FN, Eq. (18) without the factor 0.5 sqrt(G/aV)
!          PINTEG1 = PINTEG1 + PN(JMOD)*PSMAX*(ERFC(ZUPART) - &
!                  & 0.5_JPRB*(PSM(JMOD)/PSMAX)**2*PF1(JMOD)*ERFC(ZUPP))
!
!          ! integral I2: FN, Eq. (19) without the factor A/3
!          PINTEG2 = PINTEG2 + (PF2(JMOD) * PN(JMOD)/PSM(JMOD))*(ERF(ZUPM) - ERF(ZUMM))
!       ENDIF
!    ENDDO
!
!
!  END SUBROUTINE SINTEGRAL


!  SUBROUTINE AER_ACTIV_MORALES_NENES(KIDIA, KFDIA, KTDIA, KLON, KLEV, KPDF, LCLOUD, PT, PAP, PRHO, & 
!       & PQ,  PW, PWPDF, PSO4MASS, PBCMASS, POMMASS, PSSMASS, &
!       & PDUMASS,PNO3MASS, PMSAMASS, PAERONUM, PRDRY, PCDNC, KFLDX, PEXTRA)
!
!
!    ! *aer_activ_morales_nenes* calculates the number of activated aerosol 
!    !              particles from the aerosol size-distribution,
!    !              composition and ambient supersaturation
!    !
!    ! Author:
!    ! -------
!    ! Twan van Noije, KNMI
!    !
!    ! References:
!    ! -----------
!    ! Abdul-Razzak et al., JGR, 103, D6, 6123-6131, 1998
!    ! Abdul-Razzak and Ghan, JGR, 105, D5, 6837-6844, 2000 (referred to as ARG)
!    ! Pruppbacher and Klett, Kluewer Ac. Pub., 1997
!    ! Ghan et al., JAMES 3, M10001, 2011
!    ! Nenes and Seinfeld, JGR, 108, D14, 4415, 2003 (referred to as NS)
!    ! Fountoukis and Nenes, JGR, 110, D11212, 2005 (referred to as FN)
!    ! Morales Betancourt and Nenes, GMD, 7, 2345-2357, 2014 (referred to as MN)
!    ! Seinfeld and Pandis, Atmospheric Chemistry and Physics, Second Edition (referred to as SP)
!
!    USE YOMCST,              ONLY: R, RV, RPI, RCPD, RG, RLVTT, RMV, RMD, RTT , RLSTT
!    USE YOECLDP,             ONLY: RTHOMO, PPRHO_WAT
!    USE YOMHOOK,             ONLY: LHOOK, DR_HOOK
!    USE YOMLUN,              ONLY: NULOUT
!    USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
!         & R4IES    ,R5LES    ,R5IES    ,R5ALVCP  ,R5ALSCP  ,&
!         & RALVDCP  ,RALSDCP  ,RTWAT    ,&
!         & RTICE    ,RTICECU  ,&
!         & RTWAT_RTICE_R      ,RTWAT_RTICECU_R,&
!         & RKOOP1   ,RKOOP2
!    USE YOE_AERO_M7_DATA,    ONLY: NMOD, NSOL, SIGMALN, CMR2RAM, &
!         & DH2SO4, DBC, DOC, DNACL, DDUST, &
!         & DNA2SO4, DNH4NO3, DMSA, NH4NO3_FACTOR, &
!         & PPKAPPA_H2SO4, PPKAPPA_NACL, PPKAPPA_NA2SO4, &
!         & PPKAPPA_BC, PPKAPPA_OC, PPKAPPA_DU, &
!         & PPKAPPA_NH4NO3, PPKAPPA_MSA, &
!         & WSO4, WH2SO4, WNACL, WNA2SO4, &
!         & WH2O, WDAIR
!
!    IMPLICIT NONE
!
!    !---included functions from header files
!#include "fcttre.h"
!
!    !---subroutine interface
!    INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA 
!    INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA 
!    INTEGER(KIND=JPIM), INTENT(IN) :: KTDIA 
!    INTEGER(KIND=JPIM), INTENT(IN) :: KLON 
!    INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
!    INTEGER(KIND=JPIM), INTENT(IN) :: KPDF
!    INTEGER(KIND=JPIM), INTENT(IN) :: KFLDX
!
!    LOGICAL, INTENT(IN)            :: LCLOUD(KLON,KLEV)
!
!    REAL(KIND=JPRB), INTENT(IN)    :: PT(KLON,KLEV)
!    REAL(KIND=JPRB), INTENT(IN)    :: PAP(KLON,KLEV)
!    REAL(KIND=JPRB), INTENT(IN)    :: PRHO(KLON,KLEV)
!    REAL(KIND=JPRB), INTENT(IN)    :: PQ(KLON,KLEV)
!    REAL(KIND=JPRB), INTENT(IN)    :: PW(KLON,KLEV,KPDF)
!    REAL(KIND=JPRB), INTENT(IN)    :: PWPDF(KLON,KLEV,KPDF)
!    REAL(KIND=JPRB), INTENT(IN)    :: PSO4MASS(KLON,KLEV,NSOL) ! [KG(SO4)/KG(AIR)]
!    REAL(KIND=JPRB), INTENT(IN)    :: PBCMASS(KLON,KLEV,NSOL)  ! [KG(BC)/KG(AIR)]
!    REAL(KIND=JPRB), INTENT(IN)    :: POMMASS(KLON,KLEV,NSOL)  ! [KG(OM)/KG(AIR)]
!    REAL(KIND=JPRB), INTENT(IN)    :: PSSMASS(KLON,KLEV,NSOL)  ! [KG(SS)/KG(AIR)]
!    REAL(KIND=JPRB), INTENT(IN)    :: PDUMASS(KLON,KLEV,NSOL)  ! [KG(DU)/KG(AIR)]
!    REAL(KIND=JPRB), INTENT(IN)    :: PNO3MASS(KLON,KLEV)      ! [KG(NO3)/KG(AIR)]
!    REAL(KIND=JPRB), INTENT(IN)    :: PMSAMASS(KLON,KLEV)      ! [KG(MSA)/KG(AIR)]
!    REAL(KIND=JPRB), INTENT(IN)    :: PAERONUM(KLON,KLEV,NSOL)
!    REAL(KIND=JPRB), INTENT(IN)    :: PRDRY(KLON,KLEV,NSOL)
!    REAL(KIND=JPRB), INTENT(IN)    :: PEXTRA(KLON,KLEV,KFLDX)
!
!    !   Output:
!    REAL(KIND=JPRB), INTENT(INOUT) :: PCDNC(KLON,KLEV)
!
!    !---local data
!    INTEGER(KIND=JPIM) :: JL, JK, JMOD, JW
!
!    REAL(KIND=JPRB)    :: ZN(KLON,KLEV,NSOL)      ! aerosol number concentration for each mode [m-3]
!    REAL(KIND=JPRB)    :: ZSM(KLON,KLEV,NSOL)     ! critical supersaturation for activating particles
!                                                  ! with the mode number median radius
!    REAL(KIND=JPRB)    :: ZVOL(KLON)              ! total dry particle volume
!    REAL(KIND=JPRB)    :: ZKAPPA(KLON)            ! volume-weighted kappa
!    REAL(KIND=JPRB)    :: ZSMAX(KLON,KLEV,KPDF)   ! maximum supersaturation 
!    REAL(KIND=JPRB)    :: ZESW(KLON,KLEV)         ! saturation water vapour pressure
!    REAL(KIND=JPRB)    :: ZDIF(KLON,KLEV)         ! diffusivity
!    REAL(KIND=JPRB)    :: ZK(KLON,KLEV)           ! thermal conductivity
!    REAL(KIND=JPRB)    :: ZA(KLON,KLEV)           ! Kelvin coefficient
!    REAL(KIND=JPRB)    :: ZALPHA(KLON,KLEV)       ! Intermediate term in supersaturation calulation
!    REAL(KIND=JPRB)    :: ZGAMMA(KLON,KLEV)       ! Intermediate term in supersaturation calulation
!    REAL(KIND=JPRB)    :: ZGROWTH(KLON,KLEV)      ! Growth coefficient
!    REAL(KIND=JPRB)    :: ZTERM1(KLON)            ! Intermediate term in growth coefficient calulation
!    REAL(KIND=JPRB)    :: ZTERM2(KLON)            ! Intermediate term in growth coefficient calulation
!    REAL(KIND=JPRB)    :: ZTERM3(KLON)            ! Intermediate term in growth coefficient calulation
!
!    REAL(KIND=JPRB)    :: ZAMW                    ! molecular weight of water [kg mol-1]
!    REAL(KIND=JPRB)    :: ZAMD                    ! molecular weight of dry air [kg mol-1]
!
!    REAL(KIND=JPRB)    :: ZKA(KLON), ZKV(KLON)    ! Intermediate terms in supersaturation calulation
!    REAL(KIND=JPRB)    :: ZXV(KLON), ZB(KLON)     ! Intermediate terms in supersaturation calulation
!
!    REAL(KIND=JPRB)    :: ZDIFMOD(KLON)           ! Modified diffusivity
!    REAL(KIND=JPRB)    :: ZNACT                   ! Intermediate values of the activated number concentration (#/m3)
!    REAL(KIND=JPRB)    :: ZNACT_WSUM(KLON,KLEV)   ! Weighted sum of activated number concentration
!    REAL(KIND=JPRB)    :: ZPDF_NORM(KLON,KLEV)    ! Normalization factor for ZNACT_WSUM
!
!    REAL(KIND=JPRB)    :: ZSSMASS(KLON)           ! Sea salt MMR
!    REAL(KIND=JPRB)    :: ZDUMASS(KLON)           ! Dust MMR
!    REAL(KIND=JPRB)    :: ZNO3MASS(KLON)          ! Nitrate MMR
!    REAL(KIND=JPRB)    :: ZMSAMASS(KLON)          ! MSA MMR
!
!    REAL(KIND=JPRB)    :: NSO4(KLON), NH2SO4(KLON) ! Particle numbers [kmol/kg air]
!    REAL(KIND=JPRB)    :: NNACL(KLON), NNA(KLON), NCL(KLON), NNA2SO4(KLON) 
!
!    ! Per-mode constants
!    REAL(KIND=JPRB) :: ZMODECST1(NSOL), ZMODECST2(NSOL), ZMODECST3(NSOL), &
!              & ZMODECST4(NSOL), ZMODECST5(NSOL)
!
!    ! Intermediate values for the Morales & Nenes iterative scheme
!    REAL(KIND=JPRB) :: ZCF1, ZCF2, ZCF3
!    REAL(KIND=JPRB) :: ZVALUE1
!    REAL(KIND=JPRB) :: ZVALUE2
!    REAL(KIND=JPRB) :: ZVALUE3
!    REAL(KIND=JPRB) :: ZINT1
!    REAL(KIND=JPRB) :: ZINT2
!    REAL(KIND=JPRB) :: ZSMAXTEMP1, ZSMAXTEMP2, ZSMAXTEMP3
!
!    ! Control of iterative loop:
!    ! According to Ghan et al. (2011) the FN scheme 
!    ! takes about 30 interations to converge;
!    ! so safer to increase the maximum here.
!    ! xxx to be tested
!    ! xxx should be enough, rarely exceeds 30
!    ! define convergence in terms of CDNC 1E-3, i.e. 0.1%
!    INTEGER(KIND=JPIM), PARAMETER :: NMAXITER = 30
!    INTEGER(KIND=JPIM) :: NITERATIONS 
!    LOGICAL  :: LCONVERGED 
!
!    ! Miscellaneous
!    REAL(KIND=JPRB)    :: ZEPS
!    REAL(KIND=JPRB)    :: Z4PIOVER3, ZSQRT2
!    REAL(KIND=JPRB)    :: ZHOOK_HANDLE
!    REAL(KIND=JPRB), PARAMETER :: PPEPSSEC = 1.E-25_JPRB  ! used to avoid division by 0
!
!    ! mass accomodation coefficient
!    ! between 0.1 and 1.0
!    ! Raatikainen et al., 2013 PNAS
!    REAL(KIND=JPRB), PARAMETER :: PPALPHA_C = 1.E-1_JPRB
!  
!    ! Upper and lower size bounds (m) used for calculating the average water vapor diffusivity
!    REAL(KIND=JPRB), PARAMETER :: DPBIG = 5.E-6_JPRB
!    REAL(KIND=JPRB) :: DPLOW
!
!    !--- executable procedure
!    IF (LHOOK) CALL DR_HOOK('YOE_AER_ACTIV.AER_ACTIV_MORALES_NENES',0,ZHOOK_HANDLE)
!
!    !--- 0) Initializations:
!
!    ZSMAX(KIDIA:KFDIA,KTDIA:KLEV,:) = 0._JPRB
!    ZSM(KIDIA:KFDIA,KTDIA:KLEV,:) = 0._JPRB
!    ZNACT_WSUM(KIDIA:KFDIA,KTDIA:KLEV) = 0._JPRB
!    ZPDF_NORM(KIDIA:KFDIA,KTDIA:KLEV) = 0._JPRB
!    PCDNC(KIDIA:KFDIA,KTDIA:KLEV) = 0._JPRB
!
!    ZEPS=EPSILON(1._JPRB)
!
!    !--- Conversions to SI units [g mol-1 to kg mol-1]:    
!    ZAMW=WH2O*1.E-3_JPRB
!    ZAMD=WDAIR*1.E-3_JPRB
!
!    !---miscellaneous
!    Z4PIOVER3 = 4._JPRB*RPI/3._JPRB
!    ZSQRT2 = SQRT(2._JPRB)
!
!    ! FN, Eq. (24) converted to m.
!    DPLOW = MIN(0.207683E-6_JPRB * PPALPHA_C**(-0.33048_JPRB), DPBIG)
!
!    !---grid-point calculations
!    DO JK=KTDIA,KLEV
!       DO JL=KIDIA,KFDIA
!       IF (LCLOUD(JL,JK)) THEN
!          !---Kelvin (curvature) coefficient, usually denoted by capital A.
!          !   Here we use the definition of SP (introduced on p. 770),
!          !   which is adopted in NS, FN, MN,
!          !   and include an extra factor of 2 compared to the definition of Abdul-Razzak et al. (Eq. 5).
!          ZA(JL,JK) = 4._JPRB * ZAMW * PPSURFTEN / (R * PPRHO_WAT * PT(JL,JK))
!
!          !--- Abdul-Razzak et al. (1998) (Eq. 11):
!
!          ZALPHA(JL,JK) = (RG*ZAMW*RLVTT) / (RCPD*R*PT(JL,JK)*PT(JL,JK)) - &
!                        & (RG*ZAMD) / (R*PT(JL,JK))
!
!          ! Saturation water vapour pressure:
!          ZESW(JL,JK) = FOEEWM(PT(JL,JK))
!
!          !--- Following the definitions of NS, FN and MN,
!          !    gamma is defined as a dimensionless coefficient.
!          !    We use Eq. (12) from Abdul-Razzak et al. (1998) and multiply it with the air density:
!          ZGAMMA(JL,JK) = ( (R*PT(JL,JK)) / (ZESW(JL,JK)*ZAMW) +  &
!                        & (ZAMW*RLVTT*RLVTT) / (RCPD*PAP(JL,JK)*ZAMD*PT(JL,JK)) ) * PRHO(JL,JK)
!
!          !--- Diffusivity of water vapour in air (P&K, 13.3) [m2 s-1]:
!
!          ZDIF(JL,JK)=0.211_JPRB * (PT(JL,JK)/RTT)**1.94_JPRB * (101325._JPRB/PAP(JL,JK)) *1.E-4_JPRB
!
!          !--- modified diffusivity 
!          !    Average mode-independent value using FN, Eq. (23)
!          ZB(JL) = (2._JPRB*ZDIF(JL,JK)/PPALPHA_C)*SQRT(2._JPRB*RPI*ZAMW/(R*PT(JL,JK)))
!
!          !--- For any reasonable value of PPALPHA_C: DPBIG > DPLOW
!          ZDIFMOD(JL) = ZDIF(JL,JK)*(1._JPRB-(ZB(JL)/(DPBIG-DPLOW))*LOG((DPBIG+ZB(JL))/(DPLOW+ZB(JL))))
!
!          !--- Thermal conductivity zk (P&K, 13.18) [cal cm-1 s-1 K-1]:
!
!          ! Mole fraction of water:
!
!          ZXV(JL) = PQ(JL,JK)*(ZAMD/ZAMW)
!
!          ZKA(JL) = (5.69_JPRB+0.017_JPRB*(PT(JL,JK)-273.15_JPRB))*1.E-5_JPRB
!
!          ZKV(JL) = (3.78_JPRB+0.020_JPRB*(PT(JL,JK)-273.15_JPRB))*1.E-5_JPRB
!
!          ! Moist air, convert to [J m-1 s-1 K-1]:
!
!          ZK(JL,JK) = ZKA(JL)*(1._JPRB-(1.17_JPRB-1.02_JPRB*ZKV(JL)/ZKA(JL))*ZXV(JL)) &
!                    & * 4.1868_JPRB*1.E2_JPRB
!
!          !--- growth coefficient due to gas kinetic effects:
!          !--- NS, Eq. (15)
!
!          ZTERM1(JL) = (PPRHO_WAT*R*PT(JL,JK)) / (ZESW(JL,JK)*ZDIFMOD(JL)*ZAMW)
!
!          !--- Note that no size dependence is introduced in the thermal conductivity  
!          !    See FN, p. 5
!          ZTERM2(JL) = (RLVTT*PPRHO_WAT) / (ZK(JL,JK)*PT(JL,JK))
!
!          ZTERM3(JL) = (RLVTT*ZAMW) / (R*PT(JL,JK))-1._JPRB
!
!          !--- Note that the expression for G in NS (Eq. 15) has an additional factor of 4
!          !    compared to Abdul-Razzak et al. (1998), Eq. (16).
!          !    FN follow the definition from NS,
!          !    but have omitted the factor of 4 from their Eq. (12). This is a typo.
!          !    See also SP (Eq. 17.70) MN (Eq. A3).
!          ZGROWTH(JL,JK) = 4._JPRB / (ZTERM1(JL) + ZTERM2(JL) * ZTERM3(JL))
!
!       END IF
!       END DO
!    END DO
!
!    !---per-mode calculations:
!    !   soluble mode number and critical supersaturation, ignore nucleation mode (mode 1)
!    DO JMOD=2, NSOL
!       DO JK=KTDIA,KLEV
!          DO JL=KIDIA,KFDIA
!          IF (LCLOUD(JL,JK)) THEN
!          
!             !---total volume per mode [m-3 / kg(air)], used for kappa calculation
!             ZVOL(JL) = Z4PIOVER3 * PAERONUM(JL,JK,JMOD) *   &
!                      & (CMR2RAM(JMOD)*PRDRY(JL,JK,JMOD))**3 
!
!             !--- Number per unit volume [# m-3] for each mode:
!             ZN(JL,JK,JMOD) = PAERONUM(JL,JK,JMOD)*PRHO(JL,JK)
!
!             !---sea salt and dust do not exist in mode 2:
!             ZSSMASS(JL) = MERGE(0._JPRB, PSSMASS(JL,JK,JMOD), JMOD==2)
!             ZDUMASS(JL) = MERGE(0._JPRB, PDUMASS(JL,JK,JMOD), JMOD==2)
!
!             !---ammonium-nitrate and MSA do not exit in mode 2 and 4:
!             ZNO3MASS(JL) = MERGE(0._JPRB, PNO3MASS(JL,JK), JMOD==2 .OR. JMOD==4)
!             ZMSAMASS(JL) = MERGE(0._JPRB, PMSAMASS(JL,JK), JMOD==2 .OR. JMOD==4)
!
!             NNA(JL) = ZSSMASS(JL) / WNACL
!             NCL(JL) = NNA(JL)
!             NSO4(JL) = PSO4MASS(JL,JK,JMOD) / WSO4
!             NNA2SO4(JL) = MIN(NNA(JL)/2._JPRB, NSO4(JL))
!             NNA(JL) = NNA(JL) - 2._JPRB*NNA2SO4(JL)
!             NNACL(JL) = MIN(NCL(JL), NNA(JL))
!             NCL(JL) = NNACL(JL)
!             NH2SO4(JL) = NSO4(JL) - NNA2SO4(JL)
!
!             !---mode kappa = volume-weighted sum of component kappa's
!             ZKAPPA(JL) = ( (PPKAPPA_NACL * NNACL(JL) * WNACL / (DNACL*1.E3_JPRB)) + &
!                  & (PPKAPPA_NA2SO4 * NNA2SO4(JL) * WNA2SO4 / (DNA2SO4*1.E3_JPRB)) + &
!                  & (PPKAPPA_H2SO4 * NH2SO4(JL) * WH2SO4 / (DH2SO4*1.E3_JPRB)) + &
!                  & (PPKAPPA_BC * PBCMASS(JL,JK,JMOD) / (DBC*1.E3_JPRB))       + &
!                  & (PPKAPPA_OC * POMMASS(JL,JK,JMOD) / (DOC*1.E3_JPRB))       + &
!                  & (PPKAPPA_DU * ZDUMASS(JL) / (DDUST*1.E3_JPRB)) + &
!                  & (PPKAPPA_NH4NO3 * ZNO3MASS(JL) * NH4NO3_FACTOR / (DNH4NO3*1.E3_JPRB)) + &
!                  & (PPKAPPA_MSA * ZMSAMASS(JL) / (DMSA*1.E3_JPRB)) )   / &
!                  & ZVOL(JL)
!
!             !---defensive step: minimum kappa to avoid divide by zero errors
!             ZKAPPA(JL) = MERGE(ZKAPPA(JL), 0.04_JPRB, ZKAPPA(JL) > 0.04_JPRB )
!
!             !---eqn. (2) from Ghan et al (2011)
!             !   an addition factor 0.5**3 has been included in the nominator
!             !   because the Kelvin coefficient used here (ZA) 
!             !   is twice that in Ghan et al. (ZKELV)
!             ZSM(JL,JK,JMOD) = SQRT(0.5_JPRB * ZA(JL,JK)**3 / &
!                                 & (27._JPRB * ZKAPPA(JL) * PRDRY(JL,JK,JMOD)**3) )
!
!          END IF
!          END DO
!       END DO
!    END DO
!
!
!    !--- 2) Calculate maximum supersaturation:
!
!    DO JMOD=2,NSOL
!
!       !---some per-mode constants
!       ZMODECST1(JMOD) = 3.0_JPRB*SIGMALN(JMOD)/ZSQRT2      ! final term in FN, Eq. (18)
!       ZMODECST2(JMOD) = EXP(4.5_JPRB*SIGMALN(JMOD)**2)     ! factor in FN, Eq. (18)
!       
!       ZMODECST3(JMOD) = ZMODECST1(JMOD)/2._JPRB            ! final term in FN, Eq. (19)
!       ZMODECST4(JMOD) = EXP(1.125_JPRB*SIGMALN(JMOD)**2)   ! factor in FN, Eq. (19)
!       
!       ZMODECST5(JMOD) = 3.0_JPRB*ZSQRT2*SIGMALN(JMOD)      ! denominator in FN, Eqs. (8) and (20)
!
!    END DO
!
!    DO JW=1,KPDF
!       DO JK=KTDIA, KLEV
!          DO JL=KIDIA,KFDIA
!             IF (LCLOUD(JL,JK)) THEN
!                IF (ZALPHA(JL,JK) > ZEPS .AND. PW(JL,JK,JW) > ZEPS) THEN
!                   ZSMAXTEMP1 = 1.0E-5_JPRB ! min cloud supersat.
!                   ZSMAXTEMP2 = 0.1_JPRB    ! max cloud supersat.
!                   LCONVERGED = .FALSE.
!                   NITERATIONS = 0
!
!                   ZCF1 = 0.5_JPRB*SQRT(ZGROWTH(JL,JK)/(ZALPHA(JL,JK)*PW(JL,JK,JW)))
!                   ZCF2 = ZA(JL,JK)/3.0_JPRB
!                   ! When gamma is defined as a dimensionless coefficient as above,
!                   ! ZCF3 should include a division by the air density.
!                   ! See MN, Eq. 5.
!                   ! This is an error in NS (Eq. 32) and FN (Eq. 10).
!                   ! It is can be traced back to the definition of the rate of water condensation in NS, Eq. (11),
!                   ! where W is expressed in units kg(liquid water)/m3(air) instead of kg(liquid water)/kg(air).
!                   ! See also SP, Eqs. (17.73) and (17.79).
!                   ZCF3 = 0.5_JPRB*RPI*ZGAMMA(JL,JK)*PPRHO_WAT*ZGROWTH(JL,JK) &
!                        & /ZALPHA(JL,JK)/PW(JL,JK,JW)/PRHO(JL,JK)
!
!                   !DIR$ INLINE
!                   CALL SINTEGRAL_MN(ZA(JL,JK),    ZALPHA(JL,JK), ZGAMMA(JL,JK), ZGROWTH(JL,JK), &
!                                  &  PW(JL,JK,JW), ZN(JL,JK,:),   ZSM(JL,JK,:),  ZSMAXTEMP1,     &
!                                  &  ZMODECST1(:), ZMODECST2(:),  ZMODECST3(:),  ZMODECST4(:),   &
!                                  &  ZMODECST5(:), ZINT1,         ZINT2)
!                   ZVALUE1 = (ZINT1*ZCF1 + ZINT2*ZCF2) * ZCF3*ZSMAXTEMP1 - 1.0_JPRB
!                         
!                   !DIR$ INLINE
!                   CALL SINTEGRAL_MN(ZA(JL,JK),    ZALPHA(JL,JK), ZGAMMA(JL,JK), ZGROWTH(JL,JK), &
!                                  &  PW(JL,JK,JW), ZN(JL,JK,:),   ZSM(JL,JK,:),  ZSMAXTEMP2,     &
!                                  &  ZMODECST1(:), ZMODECST2(:),  ZMODECST3(:),  ZMODECST4(:),   &
!                                  &  ZMODECST5(:), ZINT1,         ZINT2)
!                   ZVALUE2 = (ZINT1*ZCF1 + ZINT2*ZCF2) * ZCF3*ZSMAXTEMP2 - 1.0_JPRB
!                
!                   DO WHILE ((.NOT.LCONVERGED) .AND. NITERATIONS < NMAXITER)
!                      NITERATIONS = NITERATIONS + 1 
!             
!                      ZSMAXTEMP3 = 0.5_JPRB*(ZSMAXTEMP1+ZSMAXTEMP2)
!                      !DIR$ INLINE
!                      CALL SINTEGRAL_MN(ZA(JL,JK),    ZALPHA(JL,JK), ZGAMMA(JL,JK), ZGROWTH(JL,JK), &
!                                     &  PW(JL,JK,JW), ZN(JL,JK,:),   ZSM(JL,JK,:),  ZSMAXTEMP3,     &
!                                     &  ZMODECST1(:), ZMODECST2(:),  ZMODECST3(:),  ZMODECST4(:),   &
!                                     &  ZMODECST5(:), ZINT1,         ZINT2)
!                      ZVALUE3 = (ZINT1*ZCF1 + ZINT2*ZCF2) * ZCF3*ZSMAXTEMP3 - 1.0_JPRB
!
!                      IF (SIGN(1.0_JPRB,ZVALUE1)*SIGN(1.0_JPRB,ZVALUE3) <= 0.0_JPRB) THEN
!                         ZVALUE2 = ZVALUE3
!                         ZSMAXTEMP2 = ZSMAXTEMP3
!                      ELSE
!                         ZVALUE1 = ZVALUE3
!                         ZSMAXTEMP1 = ZSMAXTEMP3
!                      END IF
!
!                      IF (ABS(ZSMAXTEMP2-ZSMAXTEMP1) <= 1.0E-5_JPRB*ZSMAXTEMP1) THEN
!
!                         ZSMAX(JL,JK,JW) = 0.5_JPRB*(ZSMAXTEMP1+ZSMAXTEMP2)
!                         LCONVERGED = .TRUE.
!
!                         ZNACT = 0._JPRB
!                         DO JMOD=2,NSOL
!                            IF (ZSM(JL,JK,JMOD) > ZEPS) THEN
!                               ZNACT = ZNACT + 0.5_JPRB*ZN(JL,JK,JMOD)*ERFC(2._JPRB*(LOG(ZSM(JL,JK,JMOD)/ZSMAX(JL,JK,JW)))/ZMODECST5(JMOD))
!                            ENDIF
!                         END DO
!
!                         !---Sum up the total number of activated particles, weighted by the updraft PDF [m-3]:
!                         ZNACT_WSUM(JL,JK) = ZNACT_WSUM(JL,JK) + ZNACT*PWPDF(JL,JK,JW)
!
!                         ! Including the normalization here effectively means that conditions that don't yield a converged solution
!                         ! are discarded when calculating the mean CDNC over the updraft PDF.
!                         ZPDF_NORM(JL,JK) = ZPDF_NORM(JL,JK) + PWPDF(JL,JK,JW)
!
!                      END IF
!                      
!                   END DO    ! do while
!
!                   IF (.NOT.LCONVERGED) THEN
!                      WRITE(NULOUT,*) 'WARNING: Morales and Nenes scheme not converged for updraft velocity ', PW(JL,JK,JW)
!                   ENDIF
!
!                ELSE
!
!                   ! Set contribution to activation to zero when conditions are not met
!                   ZPDF_NORM(JL,JK) = ZPDF_NORM(JL,JK) + PWPDF(JL,JK,JW)
!
!                END IF          ! END IF safe wrt divide by 0, LOG(0)...
!             END IF
!          END DO !jl
!       END DO !jk
!    END DO
!
!    DO JK=KTDIA,KLEV
!       DO JL=KIDIA,KFDIA
!          !---normalize the total number of activated particles over the PDF, and convert to [# cm-3]
!          IF (LCLOUD(JL,JK) .AND. ZPDF_NORM(JL,JK) > ZEPS ) THEN
!             PCDNC(JL,JK) = 1.E-6_JPRB * ZNACT_WSUM(JL,JK) / ZPDF_NORM(JL,JK)
!          END IF
!       END DO
!    END DO
!
!    IF (LHOOK) CALL DR_HOOK('YOE_AER_ACTIV.AER_ACTIV_MORALES_NENES',1,ZHOOK_HANDLE)
!
!  END SUBROUTINE AER_ACTIV_MORALES_NENES 


!  SUBROUTINE SINTEGRAL_MN(PA, PALPHA, PGAMMA, PGROWTH, PW, PN, PSM, &
!                          PSMAX, PUP1, PF1, PUP2, PF2, PUPD, PINTEG1, PINTEG2)
!
!    USE YOE_AERO_M7_DATA, ONLY: NSOL
!
!    IMPLICIT NONE
!
!    REAL(KIND=JPRB), INTENT(IN) :: PN(NSOL), PSM(NSOL)
!    REAL(KIND=JPRB), INTENT(IN) :: PUP1(NSOL), PF1(NSOL), PUP2(NSOL), PF2(NSOL), PUPD(NSOL)
!    REAL(KIND=JPRB), INTENT(IN) :: PA, PALPHA, PGAMMA, PGROWTH, PW, PSMAX
!    REAL(KIND=JPRB), INTENT(OUT) :: PINTEG1, PINTEG2
!
!    REAL(KIND=JPRB) :: ZUPARTPLUS, ZUPARTMIN, ZUMAX
!
!    REAL(KIND=JPRB) :: ZDELTA, ZRATIO, ZSPARTPLUS, ZSPARTMIN
!
!    REAL(KIND=JPRB) :: XIC4, XIC4_EXP
!    REAL(KIND=JPRB) :: PINTEG1_FACTOR, PINTEG1_MODE, PINTEG2_MODE
!
!    REAL(KIND=JPRB) :: ZEPS
!
!    ! Empirically derived exponent in MN, Eq. (14)
!    REAL(KIND=JPRB), PARAMETER :: EXP_VALUE = -0.3824_JPRB
!
!    INTEGER(KIND=JPIM) :: JMOD
!
!    ! note: this subroutine is intended for inline compilation and therefore does not call dr_hook
!
!    ZEPS=EPSILON(1._JPRB)
!
!    ! xi_c to power 4
!    XIC4 = 16.0_JPRB*PALPHA*PW*(PA**2)/(9.0_JPRB*PGROWTH)
!    ZDELTA = 1.0_JPRB - XIC4/PSMAX**4
!
!    IF (ZDELTA <= 0.0_JPRB) THEN
!       ! MN, Eq. (14)
!       XIC4_EXP = EXP_VALUE/4._JPRB
!       ZRATIO = (2.0E7_JPRB/3.0_JPRB)*PA*(PSMAX**EXP_VALUE-XIC4**XIC4_EXP) + 1._JPRB/SQRT(2._JPRB)
!       IF (ZRATIO > 1.0_JPRB) THEN
!          ZRATIO = 1.0_JPRB
!       END IF
!       ZSPARTPLUS = PSMAX*ZRATIO 
!       ZSPARTMIN = ZSPARTPLUS
!    ELSE
!       ! MN, Eq. (10)
!       ZSPARTPLUS = 0.5_JPRB*(1.0_JPRB + SQRT(ZDELTA))
!       ZSPARTPLUS = SQRT(ZSPARTPLUS)*PSMAX
!       ZSPARTMIN = 0.5_JPRB*(1.0_JPRB - SQRT(ZDELTA))
!       ZSPARTMIN = SQRT(ZSPARTMIN)*PSMAX
!    ENDIF
!
!    !** Calculate integrals
!
!    PINTEG1 = 0._JPRB
!    PINTEG2 = 0._JPRB
!    DO JMOD=2,NSOL
!       IF ( PSM(JMOD) > ZEPS) THEN
!
!          ! FN, Eq. (20);
!          ! u_part for max root:
!          ZUPARTPLUS    = 2.0_JPRB*LOG(PSM(JMOD)/ZSPARTPLUS)/PUPD(JMOD)
!          ! u_part for min root:
!          ZUPARTMIN     = 2.0_JPRB*LOG(PSM(JMOD)/ZSPARTMIN)/PUPD(JMOD)
!          ! u_max:
!          ZUMAX         = 2.0_JPRB*LOG(PSM(JMOD)/PSMAX)/ PUPD(JMOD)
!
!          ! The first integral describes the contribution from "population I"  to I(0,smax),
!          ! i.e. the term [I1(0,sp+) - I1(0,sp-)] in MN, Eq. (15).
!          ! The integral I1 is given in FN, Eq. (18),
!          ! but the factor 0.5 sqrt(G/aV) is not included here.
!          IF (ZDELTA > 0.0_JPRB) THEN
!             ! Contribution I1(0,sp+),
!             ! i.e. the term between curly brackets in FN, Eq. (18)
!             PINTEG1_FACTOR = 0.5_JPRB*(PSM(JMOD)/PSMAX)**2*PF1(JMOD)
!             PINTEG1_MODE = ERFC(ZUPARTPLUS) - PINTEG1_FACTOR * ERFC(ZUPARTPLUS+PUP1(JMOD))
!
!             ! Subtract I1(0,sp1-):
!             PINTEG1_MODE = PINTEG1_MODE - (ERFC(ZUPARTMIN) - PINTEG1_FACTOR * ERFC(ZUPARTMIN+PUP1(JMOD)))
!
!             ! Include factor N_i * S_max from FN, Eq. (18):
!             PINTEG1_MODE = PN(JMOD)*PSMAX*PINTEG1_MODE
! 
!             ! Add mode contribution
!             PINTEG1 = PINTEG1 + PINTEG1_MODE
!          ENDIF
!
!          ! The second integral describes the contributions from populations II and III,
!          ! i.e. the last and first terms in MN, Eq. (15).
!          ! The two contributions are given by Eq. (19) of FN,
!          ! but the factor A/3 is not included here.
!
!          ! Contribution from population II, i.e. I2(sp+,s_max)
!          PINTEG2_MODE = ERF(ZUPARTPLUS-PUP2(JMOD)) - ERF(ZUMAX-PUP2(JMOD))
!
!          ! Add contribution from population III, i.e. I2(0,sp-),
!          ! which includes an additional factor 1/sqrt(3)
!          PINTEG2_MODE = PINTEG2_MODE + (1._JPRB/SQRT(3._JPRB)) * ERFC(ZUPARTMIN-PUP2(JMOD))
!
!          ! Include prefactors in Eq. (19) of FN, but not the factor A/3.
!          PINTEG2_MODE = PINTEG2_MODE*PF2(JMOD)*PN(JMOD)/PSM(JMOD)
!  
!          ! Add mode contribution
!          PINTEG2 = PINTEG2 + PINTEG2_MODE
!       ENDIF
!    ENDDO
!
!  END SUBROUTINE SINTEGRAL_MN


  SUBROUTINE AER_ACTIV_MORALES_NENES_FULL(KIDIA, KFDIA, KTDIA, KLON, KLEV, LCLOUD, PT, PAP, PRHO, & 
       & PVERVEL, PSO4MASS, PBCMASS, POMMASS, PSSMASS, &
       & PDUMASS,PNO3MASS, PMSAMASS, PAERONUM, PRDRY, PCDNC, PSMAX, PGEMU, PSIGMA_W) !, KFLDX, PEXTRA, PSLON, PGEMU)


    ! *aer_activ_morales_nenes_full* calculates the number of activated aerosol 
    !              particles from the aerosol size-distribution,
    !              composition and ambient supersaturation
    !
    ! Author:
    ! -------
    ! Twan van Noije, KNMI
    ! Thanos Nenes, EPFL/FORTH
    !
    ! References:
    ! -----------
    ! Abdul-Razzak et al., JGR, 103, D6, 6123-6131, 1998
    ! Abdul-Razzak and Ghan, JGR, 105, D5, 6837-6844, 2000 (referred to as ARG)
    ! Pruppbacher and Klett, Kluewer Ac. Pub., 1997
    ! Ghan et al., JAMES 3, M10001, 2011
    ! Nenes and Seinfeld, JGR, 108, D14, 4415, 2003 (referred to as NS)
    ! Fountoukis and Nenes, JGR, 110, D11212, 2005 (referred to as FN)
    ! Morales Betancourt and Nenes, GMD, 7, 2345-2357, 2014 (referred to as MN)
    ! Kumar et al., ACP, 9, 2517-2532, 2009
    ! Seinfeld and Pandis, Atmospheric Chemistry and Physics, Second Edition (referred to as SP)
    ! Morales and Nenes, JGR, D18220, 2010

    USE YOMCST,              ONLY: RG, RPI
    USE YOMHOOK,             ONLY: LHOOK, DR_HOOK
    USE YOMLUN,              ONLY: NULOUT
    USE TM5M7_DATA,          ONLY: NMOD, NSOL, DDUST, DNACL, &
                                 & DOC, DBC, DH2SO4, DNA2SO4, DNH4NO3, DMSA, &
                                 & NH4NO3_FACTOR, Kap_su,Kap_pom,Kap_soa,    &
                                 & Kap_bc,Kap_ss,Kap_du,Kap_na2so4,Kap_msa,    &
                                 & Kap_no3, WSO4, WH2SO4, WNACL, WNA2SO4,    &
                                 & WH2O, WDAIR
    USE MO_HAM_M7CTL,        ONLY: CMR2RAM, SIGMA, SIGMALN       
   !  USE YOE_AERO_M7_DATA,    ONLY: NMOD, NSOL, SIGMA, SIGMALN, CMR2RAM, &
   !       & DH2SO4, DBC, DOC, DNACL, DDUST, &
   !       & DNA2SO4, DNH4NO3, DMSA, NH4NO3_FACTOR, &
   !       & PPKAPPA_H2SO4, PPKAPPA_NACL, PPKAPPA_NA2SO4, &
   !       & PPKAPPA_BC, PPKAPPA_OC, PPKAPPA_DU, &
   !       & PPKAPPA_NH4NO3, PPKAPPA_MSA, &
   !       & WSO4, WH2SO4, WNACL, WNA2SO4, &
   !       & WH2O, WDAIR

    IMPLICIT NONE

    !---subroutine interface
    INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA 
    INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA 
    INTEGER(KIND=JPIM), INTENT(IN) :: KTDIA 
    INTEGER(KIND=JPIM), INTENT(IN) :: KLON 
    INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
    !INTEGER(KIND=JPIM), INTENT(IN) :: KFLDX

    LOGICAL, INTENT(IN)            :: LCLOUD(KLON,KLEV)

    REAL(KIND=JPRB), INTENT(IN)    :: PT(KLON,KLEV)
    REAL(KIND=JPRB), INTENT(IN)    :: PAP(KLON,KLEV)
    REAL(KIND=JPRB), INTENT(IN)    :: PRHO(KLON,KLEV)
    REAL(KIND=JPRB), INTENT(IN)    :: PVERVEL(KLON,KLEV)
    REAL(KIND=JPRB), INTENT(IN)    :: PSO4MASS(KLON,KLEV,NSOL) ! [KG(SO4)/KG(AIR)]
    REAL(KIND=JPRB), INTENT(IN)    :: PBCMASS(KLON,KLEV,NSOL)  ! [KG(BC)/KG(AIR)]
    REAL(KIND=JPRB), INTENT(IN)    :: POMMASS(KLON,KLEV,NSOL)  ! [KG(OM)/KG(AIR)]
    REAL(KIND=JPRB), INTENT(IN)    :: PSSMASS(KLON,KLEV,NSOL)  ! [KG(SS)/KG(AIR)]
    REAL(KIND=JPRB), INTENT(IN)    :: PDUMASS(KLON,KLEV,NSOL)  ! [KG(DU)/KG(AIR)]
    REAL(KIND=JPRB), INTENT(IN)    :: PNO3MASS(KLON,KLEV)      ! [KG(NO3)/KG(AIR)]
    REAL(KIND=JPRB), INTENT(IN)    :: PMSAMASS(KLON,KLEV)      ! [KG(MSA)/KG(AIR)]
    REAL(KIND=JPRB), INTENT(IN)    :: PAERONUM(KLON,KLEV,NSOL)
    REAL(KIND=JPRB), INTENT(IN)    :: PRDRY(KLON,KLEV,NSOL)
    !REAL(KIND=JPRB), INTENT(IN)    :: PEXTRA(KLON,KLEV,KFLDX)
    !REAL(KIND=JPRB), INTENT(IN)    :: PSLON(KLON), PGEMU(KLON)
    REAL(KIND=JPRB), INTENT(IN)    :: PGEMU(KLON)
    REAL(KIND=JPRB), INTENT(IN)    :: PSIGMA_W(KLON,KLEV) !eehol: input sigma_w from outside [m/s]

    !   Output:
    REAL(KIND=JPRB), INTENT(INOUT) :: PCDNC(KLON,KLEV) ! # cm-3
    REAL(KIND=JPRB), INTENT(INOUT) :: PSMAX(KLON,KLEV) ! maximum supersaturation in %

    !---local data
    INTEGER(KIND=JPIM) :: JL, JK, JMOD, JW

    REAL(KIND=JPRB)    :: ZN(KLON,KLEV,NSOL)      ! aerosol number concentration for each mode [m-3]
                                                  ! with the mode number median radius
    REAL(KIND=JPRB)    :: ZVOL(KLON)              ! total dry particle volume
    REAL(KIND=JPRB)    :: ZKAPPA(KLON,KLEV,NSOL)  ! volume-weighted kappa
    REAL(KIND=JPRB)    :: ZWLARGE(KLON,KLEV)      ! large-scale velocity (m/s)
    !REAL(KIND=JPRB)    :: PSMAX(KLON,KLEV)        ! maximum supersaturation in %

    REAL(KIND=JPRB)    :: ZSSMASS(KLON)           ! Sea salt MMR
    REAL(KIND=JPRB)    :: ZDUMASS(KLON)           ! Dust MMR
    REAL(KIND=JPRB)    :: ZNO3MASS(KLON)          ! Nitrate MMR
    REAL(KIND=JPRB)    :: ZMSAMASS(KLON)          ! MSA MMR

    REAL(KIND=JPRB)    :: NSO4(KLON), NH2SO4(KLON) ! Particle numbers [kmol/kg air]
    REAL(KIND=JPRB)    :: NNACL(KLON), NNA(KLON), NCL(KLON), NNA2SO4(KLON) 

    ! Miscellaneous
    REAL(KIND=JPRB)    :: ZEPS
    REAL(KIND=JPRB)    :: Z4PIOVER3
    REAL(KIND=JPRB)    :: ZHOOK_HANDLE

    ! Variables to interface with Nenes routines
    ! These are declared double precision as they are in those routines.
    ! Can be changed to IFS style.
    DOUBLE PRECISION TPI(NSOL-1), DPGI(NSOL-1), SIGI(NSOL-1), AKKI(NSOL-1), &
                   & TPARC, PPARC, WPARC, SG(NSOL-1), NACT, SMAX

    DOUBLE PRECISION, PARAMETER :: A = 2.25 ! Default FHH adsorption parameters (in the case of FHH-AT)
    DOUBLE PRECISION, PARAMETER :: B = 1.20 ! See Kumar et al., (2011) ACP
    DOUBLE PRECISION, PARAMETER :: ACCOM = 1.0 ! Accommodation coefficient

    ! Standard deviation of the updraft velocity distribution (m/s)
    ! For the moment it is set to a constant value of 0.8 m/s,
    ! as we did for the other activation schemes.
    ! In reality it depends on turbulence characteristics.
    ! Several parameterizations have been proposed 
    ! (see e.g., Hoose et al., 2010; Zheng et al., GRL, 2016)
    ! A common approach is to use TKE or, alternatively, 
    ! the vertical diffusion coefficient 
    ! (see module VDIFLCZ in sinvect directory).
    !DOUBLE PRECISION, PARAMETER :: SIGW = 0.6_JPRB
    DOUBLE PRECISION SIGW !eehol: add sigma_w to not be a parameter but rather an input value

    ! Logical switch to use a single characteristic velocity
    ! instead of Gauss-Legendre quadrature.
    LOGICAL, PARAMETER :: CHAR_VELOCITY = .FALSE.

    INTEGER MODEI(NSOL-1)

    !--- executable procedure
    IF (LHOOK) CALL DR_HOOK('YOE_AER_ACTIV.AER_ACTIV_MORALES_NENES_FULL',0,ZHOOK_HANDLE)

    !--- 0) Initializations:

    !PCDNC(KIDIA:KFDIA,KTDIA:KLEV) = 10._JPRB
    PSMAX(KIDIA:KFDIA,KTDIA:KLEV) = 0._JPRB

    ZEPS=EPSILON(1._JPRB)

    ZWLARGE(KIDIA:KFDIA,KTDIA:KLEV) = -1._JPRB* PVERVEL(KIDIA:KFDIA,KTDIA:KLEV) / &
                                   &  (RG*PRHO(KIDIA:KFDIA,KTDIA:KLEV))

    !---miscellaneous
    Z4PIOVER3 = 4._JPRB*RPI/3._JPRB

    !---per-mode calculations:
    !   soluble mode number and critical supersaturation, ignore nucleation mode (mode 1)
    DO JMOD=2, NSOL
       DO JK=KTDIA,KLEV
          DO JL=KIDIA,KFDIA
            IF (LCLOUD(JL,JK)) THEN
            
               !---total volume per mode [m-3 / kg(air)], used for kappa calculation
               ZVOL(JL) = Z4PIOVER3 * PAERONUM(JL,JK,JMOD) *   &
                        & (CMR2RAM(JMOD)*PRDRY(JL,JK,JMOD))**3 

               !--- Number per unit volume [# m-3] for each mode:
               ZN(JL,JK,JMOD) = PAERONUM(JL,JK,JMOD)*PRHO(JL,JK)

               !---sea salt and dust do not exist in mode 2:
               ZSSMASS(JL) = MERGE(0._JPRB, PSSMASS(JL,JK,JMOD), JMOD==2)
               ZDUMASS(JL) = MERGE(0._JPRB, PDUMASS(JL,JK,JMOD), JMOD==2)

               !---ammonium-nitrate and MSA do not exit in mode 2 and 4:
               ZNO3MASS(JL) = MERGE(0._JPRB, PNO3MASS(JL,JK), JMOD==2 .OR. JMOD==4)
               ZMSAMASS(JL) = MERGE(0._JPRB, PMSAMASS(JL,JK), JMOD==2 .OR. JMOD==4)

               NNA(JL) = ZSSMASS(JL) / WNACL
               NCL(JL) = NNA(JL)
               NSO4(JL) = PSO4MASS(JL,JK,JMOD) / WSO4
               NNA2SO4(JL) = MIN(NNA(JL)/2._JPRB, NSO4(JL))
               NNA(JL) = NNA(JL) - 2._JPRB*NNA2SO4(JL)
               NNACL(JL) = MIN(NCL(JL), NNA(JL))
               NCL(JL) = NNACL(JL)
               NH2SO4(JL) = NSO4(JL) - NNA2SO4(JL)

               IF (ZVOL(JL) .GE. ZEPS) THEN !eehol: total volume per mode need to be above treshold to avoid div by zero
                  !---mode kappa = volume-weighted sum of component kappa's
                  ZKAPPA(JL,JK,JMOD) = ( (Kap_ss * NNACL(JL) * WNACL / (DNACL*1.E3_JPRB)) + &
                        & (Kap_na2so4 * NNA2SO4(JL) * WNA2SO4 / (DNA2SO4*1.E3_JPRB)) + &
                        & (Kap_su * NH2SO4(JL) * WH2SO4 / (DH2SO4*1.E3_JPRB)) + &
                        & (Kap_bc * PBCMASS(JL,JK,JMOD) / (DBC*1.E3_JPRB))       + &
                        & (Kap_pom * POMMASS(JL,JK,JMOD) / (DOC*1.E3_JPRB))       + &
                        & (Kap_du * ZDUMASS(JL) / (DDUST*1.E3_JPRB)) + &
                        & (Kap_no3 * ZNO3MASS(JL) * NH4NO3_FACTOR / (DNH4NO3*1.E3_JPRB)) + &
                        & (Kap_msa * ZMSAMASS(JL) / (DMSA*1.E3_JPRB)) )   / &
                        & ZVOL(JL)

                  !---defensive step: minimum kappa to avoid divide by zero errors
                  ZKAPPA(JL,JK,JMOD) = MERGE(ZKAPPA(JL,JK,JMOD), 0.04_JPRB, ZKAPPA(JL,JK,JMOD) > 0.04_JPRB )
                  ZKAPPA(JL,JK,JMOD)=MIN(ZKAPPA(JL,JK,JMOD),1.2)
               ELSE !eehol: if total volume per mode is too small, use minimum kappa
                  ZKAPPA(JL,JK,JMOD)=0.04_JPRB
               END IF
            END IF !eehol: LCLOUD
          END DO
       END DO
    END DO

    !--- 2) Calculate maximum supersaturation and cloud droplet number concentration, 
    !       averaged over the updraft velocity PDF

    DO JK=KTDIA, KLEV
       DO JL=KIDIA,KFDIA
          IF (LCLOUD(JL,JK)) THEN

            DO JMOD=2,NSOL
            !Shift mode index
               MODEI(JMOD-1) = 1   ! Kohler mode
               TPI(JMOD-1) = PAERONUM(JL,JK,JMOD) * PRHO(JL,JK) ! Number concentration (#/m3)
               DPGI(JMOD-1) = 2._JPRB * PRDRY(JL,JK,JMOD)   ! Modal diameter (m)
               SIGI(JMOD-1) = SIGMA(JMOD)  ! Geometric dispersion (sigma_g)
               AKKI(JMOD-1) = ZKAPPA(JL,JK,JMOD)  ! Hygroscopicity parameter (kappa)
            END DO
            TPARC = PT(JL,JK) ! Temperature (K)
            PPARC = PAP(JL,JK) ! Pressure (Pa)
         
            IF ( ANY(TPI(:) .GE. ZEPS) .AND. ANY(DPGI(:) .GE. 1e-9_JPRB) .AND. TPARC.GE.(273.15_JPRB-35.0_JPRB) ) THEN !eehol: any num con, diam and temperature need to be over treshold
               ! Convert aerosol data into CCN
               CALL CCNSPEC (TPI,DPGI,SIGI,MODEI,TPARC,PPARC,NSOL-1,AKKI,A,B,SG) 

               ! xxx To be done:
               ! Save CCN spectra for supersaturations:
               ! S = 0.05, 0.1, 0.2, 0.3, 0.5, 1.0 %
               ! New routine needs as input TPI, SG, S and returns CCN(S)
               ! which needs to be put into the output as 6 3-D fields

               !eehol: give sigma_w a value depending on the input variable
               SIGW = MAX(0.1_JPRB,PSIGMA_W(JL,JK)) !treshold sigma to min value

               IF ( ZWLARGE(JL,JK).GE.ZEPS ) THEN
                  ! Calculate activated droplet number
                  IF (CHAR_VELOCITY) THEN
                     ! Use characteristic updraft velocity (m/s)
                     ! from Morales and Nenes (2010)
                     WPARC = 0.79*SIGW
                     ! Call activation for a single velocity
                     ! equal to the characteristic velocity
                     CALL PDFACTIV (WPARC,TPI,AKKI,A,B,ACCOM,SG,0.d0,TPARC,PPARC,NACT,SMAX) 
                  ELSE
                     ! Call activation for velocity PDF (SIGW is non-zero)
                     ! with WPARC set equal to the large-scale velocity (m/s)
                     WPARC = ZWLARGE(JL,JK)
                     CALL PDFACTIV (WPARC,TPI,AKKI,A,B,ACCOM,SG,SIGW,TPARC,PPARC,NACT,SMAX) 

                  ENDIF

                  ! convert CDNC to # cm-3
                  PCDNC(JL,JK) = 1.E-6_JPRB * NACT

                  ! convert Smax to %
                  PSMAX(JL,JK) = 100._JPRB * SMAX 
               END IF
            END IF
          END IF ! LCLOUD
       END DO !jl
    END DO !jk

    IF (LHOOK) CALL DR_HOOK('YOE_AER_ACTIV.AER_ACTIV_MORALES_NENES_FULL',1,ZHOOK_HANDLE)

  END SUBROUTINE AER_ACTIV_MORALES_NENES_FULL 

!  SUBROUTINE AER_ACTIV_MORALES_NENES_FULL_OLDPDF(KIDIA, KFDIA, KTDIA, KLON, KLEV, KPDF, LCLOUD, PT, PAP, PRHO, & 
!       & PW, PWPDF, PSO4MASS, PBCMASS, POMMASS, PSSMASS, &
!       & PDUMASS,PNO3MASS, PMSAMASS, PAERONUM, PRDRY, PCDNC, KFLDX, PEXTRA, PSLON, PGEMU)
!
!
!    ! *aer_activ_morales_nenes_full_oldpdf* calculates the number of activated aerosol 
!    !              particles from the aerosol size-distribution,
!    !              composition and ambient supersaturation
!    !
!    ! Author:
!    ! -------
!    ! Twan van Noije, KNMI
!    ! Thanos Nenes, EPFL/FORTH
!    !
!    ! References:
!    ! -----------
!    ! Abdul-Razzak et al., JGR, 103, D6, 6123-6131, 1998
!    ! Abdul-Razzak and Ghan, JGR, 105, D5, 6837-6844, 2000 (referred to as ARG)
!    ! Pruppbacher and Klett, Kluewer Ac. Pub., 1997
!    ! Ghan et al., JAMES 3, M10001, 2011
!    ! Nenes and Seinfeld, JGR, 108, D14, 4415, 2003 (referred to as NS)
!    ! Fountoukis and Nenes, JGR, 110, D11212, 2005 (referred to as FN)
!    ! Morales Betancourt and Nenes, GMD, 7, 2345-2357, 2014 (referred to as MN)
!    ! Kumar et al., ACP, 9, 2517-2532, 2009
!    ! Seinfeld and Pandis, Atmospheric Chemistry and Physics, Second Edition (referred to as SP)
!    ! Morales and Nenes, JGR, D18220, 2010
!
!    USE YOMCST,              ONLY: RG, RPI
!    USE YOMHOOK,             ONLY: LHOOK, DR_HOOK
!    USE YOMLUN,              ONLY: NULOUT
!    USE YOE_AERO_M7_DATA,    ONLY: NMOD, NSOL, SIGMA, SIGMALN, CMR2RAM, &
!         & DH2SO4, DBC, DOC, DNACL, DDUST, &
!         & DNA2SO4, DNH4NO3, DMSA, NH4NO3_FACTOR, &
!         & PPKAPPA_H2SO4, PPKAPPA_NACL, PPKAPPA_NA2SO4, &
!         & PPKAPPA_BC, PPKAPPA_OC, PPKAPPA_DU, &
!         & PPKAPPA_NH4NO3, PPKAPPA_MSA, &
!         & WSO4, WH2SO4, WNACL, WNA2SO4, &
!         & WH2O, WDAIR
!
!    IMPLICIT NONE
!
!    !---subroutine interface
!    INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA 
!    INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA 
!    INTEGER(KIND=JPIM), INTENT(IN) :: KTDIA 
!    INTEGER(KIND=JPIM), INTENT(IN) :: KLON 
!    INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
!    INTEGER(KIND=JPIM), INTENT(IN) :: KPDF
!    INTEGER(KIND=JPIM), INTENT(IN) :: KFLDX
!
!    LOGICAL, INTENT(IN)            :: LCLOUD(KLON,KLEV)
!
!    REAL(KIND=JPRB), INTENT(IN)    :: PT(KLON,KLEV)
!    REAL(KIND=JPRB), INTENT(IN)    :: PAP(KLON,KLEV)
!    REAL(KIND=JPRB), INTENT(IN)    :: PRHO(KLON,KLEV)
!    REAL(KIND=JPRB), INTENT(IN)    :: PW(KLON,KLEV,KPDF)
!    REAL(KIND=JPRB), INTENT(IN)    :: PWPDF(KLON,KLEV,KPDF)
!    REAL(KIND=JPRB), INTENT(IN)    :: PSO4MASS(KLON,KLEV,NSOL) ! [KG(SO4)/KG(AIR)]
!    REAL(KIND=JPRB), INTENT(IN)    :: PBCMASS(KLON,KLEV,NSOL)  ! [KG(BC)/KG(AIR)]
!    REAL(KIND=JPRB), INTENT(IN)    :: POMMASS(KLON,KLEV,NSOL)  ! [KG(OM)/KG(AIR)]
!    REAL(KIND=JPRB), INTENT(IN)    :: PSSMASS(KLON,KLEV,NSOL)  ! [KG(SS)/KG(AIR)]
!    REAL(KIND=JPRB), INTENT(IN)    :: PDUMASS(KLON,KLEV,NSOL)  ! [KG(DU)/KG(AIR)]
!    REAL(KIND=JPRB), INTENT(IN)    :: PNO3MASS(KLON,KLEV)      ! [KG(NO3)/KG(AIR)]
!    REAL(KIND=JPRB), INTENT(IN)    :: PMSAMASS(KLON,KLEV)      ! [KG(MSA)/KG(AIR)]
!    REAL(KIND=JPRB), INTENT(IN)    :: PAERONUM(KLON,KLEV,NSOL)
!    REAL(KIND=JPRB), INTENT(IN)    :: PRDRY(KLON,KLEV,NSOL)
!    REAL(KIND=JPRB), INTENT(IN)    :: PEXTRA(KLON,KLEV,KFLDX)
!    REAL(KIND=JPRB), INTENT(IN)    :: PSLON(KLON), PGEMU(KLON)
!
!    !   Output:
!    REAL(KIND=JPRB), INTENT(OUT) :: PCDNC(KLON,KLEV) ! # cm-3
!    !REAL(KIND=JPRB), INTENT(OUT) :: PSMAX(KLON,KLEV) ! maximum supersaturation in %
!
!    !---local data
!    INTEGER(KIND=JPIM) :: JL, JK, JMOD, JW
!
!    REAL(KIND=JPRB)    :: ZN(KLON,KLEV,NSOL)      ! aerosol number concentration for each mode [m-3]
!                                                  ! with the mode number median radius
!    REAL(KIND=JPRB)    :: ZVOL(KLON)              ! total dry particle volume
!    REAL(KIND=JPRB)    :: ZKAPPA(KLON,KLEV,NSOL)  ! volume-weighted kappa
!    REAL(KIND=JPRB)    :: ZWLARGE(KLON,KLEV)      ! large-scale velocity (m/s)
!    REAL(KIND=JPRB)    :: PSMAX(KLON,KLEV)        ! maximum supersaturation in %
!
!    REAL(KIND=JPRB)    :: ZNACT_WSUM(KLON,KLEV)   ! Weighted sum of activated number concentration
!    REAL(KIND=JPRB)    :: ZSMAX_WSUM(KLON,KLEV)   ! Weighted sum of Smax
!    REAL(KIND=JPRB)    :: ZPDF_NORM(KLON,KLEV)    ! Normalization factor for ZNACT_WSUM
!
!    REAL(KIND=JPRB)    :: ZSSMASS(KLON)           ! Sea salt MMR
!    REAL(KIND=JPRB)    :: ZDUMASS(KLON)           ! Dust MMR
!    REAL(KIND=JPRB)    :: ZNO3MASS(KLON)          ! Nitrate MMR
!    REAL(KIND=JPRB)    :: ZMSAMASS(KLON)          ! MSA MMR
!
!    REAL(KIND=JPRB)    :: NSO4(KLON), NH2SO4(KLON) ! Particle numbers [kmol/kg air]
!    REAL(KIND=JPRB)    :: NNACL(KLON), NNA(KLON), NCL(KLON), NNA2SO4(KLON) 
!
!    ! Miscellaneous
!    REAL(KIND=JPRB)    :: ZEPS
!    REAL(KIND=JPRB)    :: Z4PIOVER3
!    REAL(KIND=JPRB)    :: ZHOOK_HANDLE
!
!    ! Variables to interface with Nenes routines
!    ! These are declared double precision as they are in those routines.
!    ! Can be changed to IFS style.
!    DOUBLE PRECISION TPI(NSOL-1), DPGI(NSOL-1), SIGI(NSOL-1), AKKI(NSOL-1), &
!                   & TPARC, PPARC, WPARC, SG(NSOL-1), NACT, SMAX
!
!    DOUBLE PRECISION, PARAMETER :: A = 2.25_JPRB ! Default FHH adsorption parameters (in the case of FHH-AT)
!    DOUBLE PRECISION, PARAMETER :: B = 1.20_JPRB ! See Kumar et al., (2011) ACP
!    DOUBLE PRECISION, PARAMETER :: ACCOM = 1.0_JPRB ! Accommodation coefficient
!
!    INTEGER MODEI(NSOL-1)
!
!    !--- executable procedure
!    IF (LHOOK) CALL DR_HOOK('YOE_AER_ACTIV.AER_ACTIV_MORALES_NENES_FULL_OLDPDF',0,ZHOOK_HANDLE)
!
!    !--- 0) Initializations:
!
!    ZNACT_WSUM(KIDIA:KFDIA,KTDIA:KLEV) = 0._JPRB
!    ZSMAX_WSUM(KIDIA:KFDIA,KTDIA:KLEV) = 0._JPRB
!    ZPDF_NORM(KIDIA:KFDIA,KTDIA:KLEV) = 0._JPRB
!    PCDNC(KIDIA:KFDIA,KTDIA:KLEV) = 0._JPRB
!    PSMAX(KIDIA:KFDIA,KTDIA:KLEV) = 0._JPRB
!
!    ZEPS=EPSILON(1._JPRB)
!
!    !---miscellaneous
!    Z4PIOVER3 = 4._JPRB*RPI/3._JPRB
!
!    !---per-mode calculations:
!    !   soluble mode number and critical supersaturation, ignore nucleation mode (mode 1)
!    DO JMOD=2, NSOL
!       DO JK=KTDIA,KLEV
!          DO JL=KIDIA,KFDIA
!          IF (LCLOUD(JL,JK)) THEN
!          
!             !---total volume per mode [m-3 / kg(air)], used for kappa calculation
!             ZVOL(JL) = Z4PIOVER3 * PAERONUM(JL,JK,JMOD) *   &
!                      & (CMR2RAM(JMOD)*PRDRY(JL,JK,JMOD))**3 
!
!             !--- Number per unit volume [# m-3] for each mode:
!             ZN(JL,JK,JMOD) = PAERONUM(JL,JK,JMOD)*PRHO(JL,JK)
!
!             !---sea salt and dust do not exist in mode 2:
!             ZSSMASS(JL) = MERGE(0._JPRB, PSSMASS(JL,JK,JMOD), JMOD==2)
!             ZDUMASS(JL) = MERGE(0._JPRB, PDUMASS(JL,JK,JMOD), JMOD==2)
!
!             !---ammonium-nitrate and MSA do not exit in mode 2 and 4:
!             ZNO3MASS(JL) = MERGE(0._JPRB, PNO3MASS(JL,JK), JMOD==2 .OR. JMOD==4)
!             ZMSAMASS(JL) = MERGE(0._JPRB, PMSAMASS(JL,JK), JMOD==2 .OR. JMOD==4)
!
!             NNA(JL) = ZSSMASS(JL) / WNACL
!             NCL(JL) = NNA(JL)
!             NSO4(JL) = PSO4MASS(JL,JK,JMOD) / WSO4
!             NNA2SO4(JL) = MIN(NNA(JL)/2._JPRB, NSO4(JL))
!             NNA(JL) = NNA(JL) - 2._JPRB*NNA2SO4(JL)
!             NNACL(JL) = MIN(NCL(JL), NNA(JL))
!             NCL(JL) = NNACL(JL)
!             NH2SO4(JL) = NSO4(JL) - NNA2SO4(JL)
!
!             !---mode kappa = volume-weighted sum of component kappa's
!             ZKAPPA(JL,JK,JMOD) = ( (PPKAPPA_NACL * NNACL(JL) * WNACL / (DNACL*1.E3_JPRB)) + &
!                  & (PPKAPPA_NA2SO4 * NNA2SO4(JL) * WNA2SO4 / (DNA2SO4*1.E3_JPRB)) + &
!                  & (PPKAPPA_H2SO4 * NH2SO4(JL) * WH2SO4 / (DH2SO4*1.E3_JPRB)) + &
!                  & (PPKAPPA_BC * PBCMASS(JL,JK,JMOD) / (DBC*1.E3_JPRB))       + &
!                  & (PPKAPPA_OC * POMMASS(JL,JK,JMOD) / (DOC*1.E3_JPRB))       + &
!                  & (PPKAPPA_DU * ZDUMASS(JL) / (DDUST*1.E3_JPRB)) + &
!                  & (PPKAPPA_NH4NO3 * ZNO3MASS(JL) * NH4NO3_FACTOR / (DNH4NO3*1.E3_JPRB)) + &
!                  & (PPKAPPA_MSA * ZMSAMASS(JL) / (DMSA*1.E3_JPRB)) )   / &
!                  & ZVOL(JL)
!
!             !---defensive step: minimum kappa to avoid divide by zero errors
!             ZKAPPA(JL,JK,JMOD) = MERGE(ZKAPPA(JL,JK,JMOD), 0.04_JPRB, ZKAPPA(JL,JK,JMOD) > 0.04_JPRB )
!
!          END IF
!          END DO
!       END DO
!    END DO
!
!    !--- 2) Calculate maximum supersaturation and cloud droplet number concentration, 
!    !       averaged over the updraft velocity PDF
!
!    DO JW=1,KPDF
!       DO JK=KTDIA, KLEV
!          DO JL=KIDIA,KFDIA
!             IF (LCLOUD(JL,JK)) THEN
!
!                DO JMOD=2,NSOL
!                !Shift mode index
!                   MODEI(JMOD-1) = 1   ! Kohler mode
!                   TPI(JMOD-1) = PAERONUM(JL,JK,JMOD) * PRHO(JL,JK) ! Number concentration (#/m3)
!                   DPGI(JMOD-1) = 2._JPRB * PRDRY(JL,JK,JMOD)   ! Modal diameter (m)
!                   SIGI(JMOD-1) = SIGMA(JMOD)  ! Geometric dispersion (sigma_g)
!                   AKKI(JMOD-1) = ZKAPPA(JL,JK,JMOD)  ! Hygroscopicity parameter (kappa)
!                END DO
!                TPARC = PT(JL,JK) ! Temperature (K)
!                PPARC = PAP(JL,JK) ! Pressure (Pa)
!
!                ! Convert aerosol data into CCN
!                CALL CCNSPEC (TPI,DPGI,SIGI,MODEI,TPARC,PPARC,NSOL-1,AKKI,A,B,SG) 
!
!                ! xxx To be done:
!                ! Save CCN spectra for supersaturations:
!                ! S = 0.05, 0.1, 0.2, 0.3, 0.5, 1.0 %
!                ! New routine needs as input TPI, SG, S and returns CCN(S)
!                ! which needs to be put into the output as 6 3-D fields
!
!                ! Set WPARC to the updraft vertical velocity (m/s) 
!                WPARC = PW(JL,JK,JW)
!                ! Call activation for a single velocity value:
!                CALL PDFACTIV (WPARC,TPI,AKKI,A,B,ACCOM,SG,0.d0,TPARC,PPARC, NACT, SMAX)
!
!                !IF(JW==5 .AND. PSLON(JL)>0.17 .AND. PSLON(JL)<0.19 .AND. PGEMU(JL)>0.76 .AND. PGEMU(JL)<0.78) THEN
!                !   open(unit=666, file='parame', access='append', status='unknown')
!                !   write(666,*) modei, tpi, dpgi, sigi, akki, tparc, pparc, wparc, accom, sigw, sg, nact, smax
!                !   close(666)
!                !ENDIF
!
!                !IF(JW==5 .AND. PSLON(JL)>-0.35 .AND. PSLON(JL)<-0.33 .AND. PGEMU(JL)>0.76 .AND. PGEMU(JL)<0.78) THEN
!                !   open(unit=666, file='paramp', access='append', status='unknown')
!                !   write(666,*) modei, tpi, dpgi, sigi, akki, tparc, pparc, wparc, accom, sigw, sg, nact, smax
!                !   close(666)
!                !ENDIF
!
!                !---Sum up the total number of activated particles, weighted by the updraft PDF [m-3]:
!                ZNACT_WSUM(JL,JK) = ZNACT_WSUM(JL,JK) + NACT*PWPDF(JL,JK,JW)
!                ZSMAX_WSUM(JL,JK) = ZSMAX_WSUM(JL,JK) + SMAX*PWPDF(JL,JK,JW)
!
!                ZPDF_NORM(JL,JK) = ZPDF_NORM(JL,JK) + PWPDF(JL,JK,JW)
!
!             END IF ! LCLOUD
!          END DO !jl
!       END DO !jk
!    END DO ! PDF
!
!    DO JK=KTDIA,KLEV
!       DO JL=KIDIA,KFDIA
!          !---normalize over the PDF
!          IF (LCLOUD(JL,JK) .AND. ZPDF_NORM(JL,JK) > ZEPS ) THEN
!             ! convert CDNC to # cm-3
!             PCDNC(JL,JK) = 1.E-6_JPRB * ZNACT_WSUM(JL,JK) / ZPDF_NORM(JL,JK)
!             
!             ! convert Smax to %
!             PSMAX(JL,JK) = 100._JPRB * ZSMAX_WSUM(JL,JK) / ZPDF_NORM(JL,JK)
!          END IF
!       END DO
!    END DO
!
!    IF (LHOOK) CALL DR_HOOK('YOE_AER_ACTIV.AER_ACTIV_MORALES_NENES_FULL_OLDPDF',1,ZHOOK_HANDLE)
!
!  END SUBROUTINE AER_ACTIV_MORALES_NENES_FULL_OLDPDF


!  SUBROUTINE AER_ACTIV_ABDULRAZZAK_GHAN (KIDIA, KFDIA, KTDIA, KLON, KLEV, KPDF, LCLOUD, PT, PAP, PRHO, & 
!       & PQ,  PW, PWPDF, PSO4MASS, PBCMASS, POMMASS, PSSMASS, &
!       & PDUMASS, PNO3MASS, PMSAMASS, PAERONUM, PRDRY, PCDNC, KFLDX, PEXTRA)
!
!
!    ! *aer_activ_abdulrazzak_ghan* calculates the number of activated aerosol 
!    !              particles from the aerosol size-distribution,
!    !              composition and ambient supersaturation
!    !
!    ! Author:
!    ! -------
!    ! Philip Stier, MPI-MET, University of Oxford  2002-2009
!    !
!    ! Method:
!    ! -------
!    ! The calculation of the activation can be reduced to 3 tasks:
!    ! 
!    ! I)   Calculate the maximum supersaturation
!    ! II)  Calculate the corresponding radius of activation
!    !      for each mode
!    ! III) Calculate the number of particles that are larger
!    !      then the radius of activation for each mode.
!    ! 
!    ! III) Calculation of the number of activated particles:
!    !      See the routine ham_logtail below.
!    !
!    ! The calculations are now performed separately for 
!    ! stratiform and convective updraft velocities.
!    !
!    ! References:
!    ! -----------
!    ! Abdul-Razzak et al., JGR, 103, D6, 6123-6131, 1998.
!    ! Abdul-Razzak and Ghan, JGR, 105, D5, 6837-6844, 2000.
!    ! Pruppbacher and Klett, Kluewer Ac. Pub., 1997.
!    ! Ghan et al., JAMES 3, M10001, 2011
!
!    !---inherited functions, types, variables and constants
!    USE YOMCST,              ONLY: RPI, R, RCPD, RG, RMD, RV, RTT, RLVTT, RLSTT, RMD, RMV
!    USE YOECLDP,             ONLY: RTHOMO, PPRHO_WAT, LACI_DIAG
!    USE YOMHOOK,             ONLY: LHOOK, DR_HOOK
!    USE YOE_AERO_M7_DATA,     ONLY: NSOL, SIGMALN, CMR2RAM,  &
!         & PPKAPPA_H2SO4, PPKAPPA_NACL, PPKAPPA_NA2SO4, &
!         & PPKAPPA_BC, PPKAPPA_OC, PPKAPPA_DU, &
!         & PPKAPPA_NH4NO3, PPKAPPA_MSA, &
!         & DH2SO4, DNACL, DNA2SO4, &
!         & DBC, DOC, DDUST, &
!         & DNH4NO3, DMSA, NH4NO3_FACTOR, &
!         & WSO4, WH2SO4, WNACL, WNA2SO4
!    USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
!         & R4IES    ,R5LES    ,R5IES    ,R5ALVCP  ,R5ALSCP  ,&
!         & RALVDCP  ,RALSDCP  ,RTWAT    ,&
!         & RTICE    ,RTICECU  ,&
!         & RTWAT_RTICE_R      ,RTWAT_RTICECU_R,&
!         & RKOOP1   ,RKOOP2
!
!    IMPLICIT NONE
!
!    !---included functions from header files
!#include "fcttre.h"
!
!    !---subroutine interface
!    INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA 
!    INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA 
!    INTEGER(KIND=JPIM), INTENT(IN) :: KTDIA 
!    INTEGER(KIND=JPIM), INTENT(IN) :: KLON 
!    INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
!    INTEGER(KIND=JPIM), INTENT(IN) :: KPDF
!    INTEGER(KIND=JPIM), INTENT(IN) :: KFLDX
!
!    LOGICAL, INTENT(IN)            :: LCLOUD(KLON,KLEV)
!
!    REAL(KIND=JPRB), INTENT(IN)    :: PT(KLON,KLEV)
!    REAL(KIND=JPRB), INTENT(IN)    :: PAP(KLON,KLEV)
!    REAL(KIND=JPRB), INTENT(IN)    :: PRHO(KLON,KLEV)
!    REAL(KIND=JPRB), INTENT(IN)    :: PQ(KLON,KLEV)
!    REAL(KIND=JPRB), INTENT(IN)    :: PW(KLON,KLEV,KPDF)
!    REAL(KIND=JPRB), INTENT(IN)    :: PWPDF(KLON,KLEV,KPDF)
!    REAL(KIND=JPRB), INTENT(IN)    :: PSO4MASS(KLON,KLEV,NSOL) ! [KG(SO4)/KG(AIR)]
!    REAL(KIND=JPRB), INTENT(IN)    :: PBCMASS(KLON,KLEV,NSOL)  ! [KG(BC)/KG(AIR)]
!    REAL(KIND=JPRB), INTENT(IN)    :: POMMASS(KLON,KLEV,NSOL)  ! [KG(OM)/KG(AIR)]
!    REAL(KIND=JPRB), INTENT(IN)    :: PSSMASS(KLON,KLEV,NSOL)  ! [KG(SS)/KG(AIR)]
!    REAL(KIND=JPRB), INTENT(IN)    :: PDUMASS(KLON,KLEV,NSOL)  ! [KG(DU)/KG(AIR)]
!    REAL(KIND=JPRB), INTENT(IN)    :: PNO3MASS(KLON,KLEV)      ! [KG(NO3)/KG(AIR)]
!    REAL(KIND=JPRB), INTENT(IN)    :: PMSAMASS(KLON,KLEV)      ! [KG(MSA)/KG(AIR)]
!    REAL(KIND=JPRB), INTENT(IN)    :: PAERONUM(KLON,KLEV,NSOL)
!    REAL(KIND=JPRB), INTENT(IN)    :: PRDRY(KLON,KLEV,NSOL)
!    
!    !   Output:
!    REAL(KIND=JPRB), INTENT(INOUT) :: PCDNC(KLON,KLEV)
!    REAL(KIND=JPRB), INTENT(INOUT) :: PEXTRA(KLON,KLEV,KFLDX)
!
!    !---local data
!    INTEGER(KIND=JPIM) :: JL, JK, JMOD, JW
!
!    REAL(KIND=JPRB)    :: ZRC                     ! critical radius of a dry aerosol particle that becomes
!                                                  ! activated at the ambient supersaturation
!    REAL(KIND=JPRB)    :: ZN(KLON,KLEV,NSOL)      ! aerosol number concentration for each mode [m-3]
!    REAL(KIND=JPRB)    :: ZSM(KLON,KLEV,NSOL)     ! critical supersaturation for activating particles
!                                                  ! with the mode number median radius
!    REAL(KIND=JPRB)    :: ZVOL(KLON )             ! total dry particle volume per mode
!    REAL(KIND=JPRB)    :: ZKAPPA(KLON)            ! volume-weighted kappa per mode
!    REAL(KIND=JPRB)    :: ZSMAX(KLON,KLEV,KPDF)   ! maximum supersaturation 
!
!    REAL(KIND=JPRB)    :: ZKELV(KLON,KLEV)        ! Kelvin term 
!    REAL(KIND=JPRB)    :: ZESW(KLON,KLEV)         ! saturation vapour pressure
!    REAL(KIND=JPRB)    :: ZDIF(KLON,KLEV)         ! diffusivity of air
!    REAL(KIND=JPRB)    :: ZALPHA(KLON,KLEV)       ! Intermediate term in supersaturation calulation
!    REAL(KIND=JPRB)    :: ZGAMMA(KLON,KLEV)       ! Intermediate term in supersaturation calulation
!    REAL(KIND=JPRB)    :: ZSUM(KLON,KLEV,KPDF)    ! Intermediate term in supersaturation calulation
!    REAL(KIND=JPRB)    :: ZGROWTH(KLON,KLEV,NSOL) ! Growth due to vapour condensation
!    REAL(KIND=JPRB)    :: ZFRACN_TOP(KLON,KLEV,NSOL)
!    REAL(KIND=JPRB)    :: ZFRACN_BOT(KLON,KLEV,NSOL)
!
!    REAL(KIND=JPRB)    :: ZAMW                    ! molecular weight of water [kg mol-1]
!    REAL(KIND=JPRB)    :: ZAMD                    ! molecular weight of dry air [kg mol-1]
!
!    REAL(KIND=JPRB)    :: ZFRACN                  ! fraction of aerosol number activated for each mode 
!    REAL(KIND=JPRB)    :: ZF(NSOL), ZG(NSOL)      ! aerosol mode attributes, see AR&G
!    REAL(KIND=JPRB)    :: ZKA(KLON), ZKV(KLON)    ! thermal conductivity of dry air/water vapour
!    REAL(KIND=JPRB)    :: ZK(KLON,KLEV)           ! thermal conductivity of moist air
!    REAL(KIND=JPRB)    :: ZSQTERM(KLON)           ! Intermediate terms in supersaturation calulation
!    REAL(KIND=JPRB)    :: ZXV, ZXI, ZETA          ! Intermediate terms in supersaturation calulation
!    REAL(KIND=JPRB)    :: ZDIFMOD(KLON)           ! Modified diffusion coefficient
!    REAL(KIND=JPRB)    :: ZKMOD(KLON)             ! Intermediate terms in supersaturation calulation
!    REAL(KIND=JPRB)    :: ZT                      ! Intermediate term in activation calulation
!    REAL(KIND=JPRB)    :: ZTERM1(KLON)            ! Intermediate term in activation calulation
!    REAL(KIND=JPRB)    :: ZTERM2(KLON)            ! Intermediate term in activation calulation
!    REAL(KIND=JPRB)    :: ZTERM3(KLON)            ! Intermediate term in activation calulation
!    REAL(KIND=JPRB)    :: ZSSMASS(KLON)           ! Sea salt MMR
!    REAL(KIND=JPRB)    :: ZDUMASS (KLON)          ! Dust MMR
!    REAL(KIND=JPRB)    :: ZNO3MASS(KLON)          ! Nitrate MMR
!    REAL(KIND=JPRB)    :: ZMSAMASS(KLON)          ! MSA MMR
!
!    REAL(KIND=JPRB)    :: NSO4(KLON), NH2SO4(KLON)! Particle numbers [kmol/kg air]
!    REAL(KIND=JPRB)    :: NNACL(KLON), NNA(KLON), NCL(KLON), NNA2SO4(KLON) 
!
!    REAL(KIND=JPRB)    :: ZEPS
!    REAL(KIND=JPRB)    :: Z4PIOVER3, ZSQRT2
!    REAL(KIND=JPRB), PARAMETER :: PPEPSSEC = 1.E-25_JPRB  ! used to avoid division by 0
!
!    ! the following constants from Abdul-Razzak et al. 1998 Ch. 3:
!    ! thermal jump length [m]
!    REAL(KIND=JPRB), PARAMETER :: PPDELTA_T = 2.16E-7_JPRB
!
!    ! vapour jump length [m]
!    REAL(KIND=JPRB), PARAMETER :: PPDELTA_V = 1.096E-7_JPRB
!
!    ! mass accomodation coefficient 
!    REAL(KIND=JPRB), PARAMETER :: PPALPHA_C = 1.0_JPRB
!
!    ! thermal accomodation coefficient
!    REAL(KIND=JPRB), PARAMETER :: PPALPHA_T = 0.96_JPRB
!
!    REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!    !---executable procedure
!    IF (LHOOK) CALL DR_HOOK('YOE_AER_ACTIV.AER_ACTIV_ABDULRAZZAK_GHAN', 0, ZHOOK_HANDLE)
!
!    !---Initializations:
!
!    ZSUM(KIDIA:KFDIA,KTDIA:KLEV,:)= 0._JPRB
!    ZSMAX(KIDIA:KFDIA,KTDIA:KLEV,:)= 0._JPRB
!    ZSM(KIDIA:KFDIA,KTDIA:KLEV,:) = 0._JPRB
!    ZFRACN_TOP(KIDIA:KFDIA,KTDIA:KLEV,:) = 0._JPRB
!    ZFRACN_BOT(KIDIA:KFDIA,KTDIA:KLEV,:) = 0._JPRB
!    PCDNC(KIDIA:KFDIA,KTDIA:KLEV) = 0._JPRB
!
!    !--- Conversions to SI units [g mol-1 to kg mol-1]:    
!    ZAMW=RMV*1.E-3_JPRB
!    ZAMD=RMD*1.E-3_JPRB
!
!    !---miscellaneous
!    Z4PIOVER3 = 4._JPRB*RPI/3._JPRB
!    ZEPS=EPSILON(1._JPRB)
!    ZSQRT2 = SQRT(2._JPRB)
!
!    !   Abdul-Razzak and Ghan (2000):
!    !   (Equations numbers from this paper unless otherwise quoted)
!
!    DO JK=KTDIA,KLEV
!       DO JL=KIDIA,KFDIA
!          
!          IF (LCLOUD(JL,JK)) THEN
!             !---Kelvin (curvature) term (equation (5) in Abdul-Razzak et al. 1998) 
!             !   also (A1) in Ghan et al. 2011
!             ZKELV(JL,JK) = 2._JPRB * ZAMW * PPSURFTEN / (R * PPRHO_WAT * PT(JL,JK))
!
!             !--- Abdul-Razzak et al. (1998) (Eq. 11):
!             ZALPHA(JL,JK) = (RG*ZAMW*RLVTT) / (RCPD*R*PT(JL,JK)**2) - &
!                  & (RG*ZAMD) / (R*PT(JL,JK))
!
!             ! Saturation water vapour pressure:
!             ZESW(JL,JK) = FOEEWM(PT(JL,JK))
!
!             !--- Abdul-Razzak et al. (1998) (Eq. 12):
!             ZGAMMA(JL,JK) = (R*PT(JL,JK)) / (ZESW(JL,JK)*ZAMW) +  &
!                  & (ZAMW*RLVTT**2) / (RCPD*PAP(JL,JK)*ZAMD*PT(JL,JK))
!
!             !--- Diffusivity of water vapour in air (P&K, 13.3) [m2 s-1]:
!             
!             ZDIF(JL,JK)=0.211_JPRB * (PT(JL,JK)/RTT)**1.94_JPRB * (101325._JPRB/PAP(JL,JK)) *1.E-4_JPRB
!
!             !--- Thermal conductivity zk (P&K, 13.18) [cal cm-1 s-1 K-1]:
!
!             ! Mole fraction of water:
!
!             ZKA(JL)=(5.69_JPRB+0.017_JPRB*(PT(JL,JK)-273.15_JPRB))*1.E-5_JPRB
!
!             ZKV(JL)=(3.78_JPRB+0.020_JPRB*(PT(JL,JK)-273.15_JPRB))*1.E-5_JPRB
!
!             ! Moist air, convert to [J m-1 s-1 K-1]:
!
!             ZK(JL,JK) = ZKA(JL)*(1._JPRB-(1.17_JPRB-1.02_JPRB*ZKV(JL)/ZKA(JL))*PQ(JL,JK)*(ZAMD/ZAMW)) &
!                  & * 4.1868_JPRB*1.E2_JPRB
!          END IF
!
!       END DO
!    END DO
!
!    DO JMOD=2,NSOL
!       ! (7):
!
!       ZF(JMOD)=0.5_JPRB*EXP(2.5_JPRB*SIGMALN(JMOD)**2)
!
!       ! (8):
!
!       ZG(JMOD)=1._JPRB+0.25_JPRB*SIGMALN(JMOD)
!
!       ! (9):
!       !>>dod kappa
!       ! replaced by eqn. (2) from Ghan et al (2011)
!       DO JK=KTDIA,KLEV
!          DO JL=KIDIA,KFDIA
!
!             IF (LCLOUD(JL,JK)) THEN
!                !---total volume per mode [m-3 / kg(air)], used for kappa calculation
!                ZVOL(JL) = Z4PIOVER3 * PAERONUM(JL,JK,JMOD) *       &
!                     (CMR2RAM(JMOD)*PRDRY(JL,JK,JMOD))**3
!
!                !--- Number per unit volume [# m-3] for each mode:
!                ZN(JL,JK,JMOD) = PAERONUM(JL,JK,JMOD)*PRHO(JL,JK)
!
!                !---sea salt and dust do not exist in mode 2:
!                ZSSMASS(JL) = MERGE(0._JPRB, PSSMASS(JL,JK,JMOD), JMOD==2)
!                ZDUMASS(JL) = MERGE(0._JPRB, PDUMASS(JL,JK,JMOD), JMOD==2)
!                !---ammonium-nitrate and MSA do not exit in mode 2 and 4:
!                ZNO3MASS(JL) = MERGE(0._JPRB, PNO3MASS(JL,JK), JMOD==2 .OR. JMOD==4)
!                ZMSAMASS(JL) = MERGE(0._JPRB, PMSAMASS(JL,JK), JMOD==2 .OR. JMOD==4)
!
!                NNA(JL) = ZSSMASS(JL) / WNACL
!                NCL(JL) = NNA(JL)
!                NSO4(JL) = PSO4MASS(JL,JK,JMOD) / WSO4
!                NNA2SO4(JL) = MIN(NNA(JL)/2._JPRB, NSO4(JL))
!                NNA(JL) = NNA(JL) - 2._JPRB*NNA2SO4(JL)
!                NNACL(JL) = MIN(NCL(JL), NNA(JL))
!                NCL(JL) = NNACL(JL)
!                NH2SO4(JL) = NSO4(JL) - NNA2SO4(JL)
!
!                !---mode kappa = volume-weighted sum of component kappa's
!                ZKAPPA(JL) = ( (PPKAPPA_NACL * NNACL(JL) * WNACL / (DNACL*1.E3_JPRB)) + &
!                     & (PPKAPPA_NA2SO4 * NNA2SO4(JL) * WNA2SO4 / (DNA2SO4*1.E3_JPRB)) + &
!                     & (PPKAPPA_H2SO4 * NH2SO4(JL) * WH2SO4 / (DH2SO4*1.E3_JPRB)) + &
!                     & (PPKAPPA_BC * PBCMASS(JL,JK,JMOD) / (DBC*1.E3_JPRB))       + &
!                     & (PPKAPPA_OC * POMMASS(JL,JK,JMOD) / (DOC*1.E3_JPRB))       + &
!                     & (PPKAPPA_DU * ZDUMASS(JL) / (DDUST*1.E3_JPRB)) + &
!                     & (PPKAPPA_NH4NO3 * ZNO3MASS(JL) * NH4NO3_FACTOR / (DNH4NO3*1.E3_JPRB)) + &
!                     & (PPKAPPA_MSA * ZMSAMASS(JL) / (DMSA*1.E3_JPRB)) )   / &
!                     & ZVOL(JL)
!
!                !---defensive step: minimum kappa to avoid divide by zero errors
!                ZKAPPA(JL) = MERGE(ZKAPPA(JL), 0.04_JPRB, ZKAPPA(JL) > 0.04_JPRB )
!
!                ZSM(JL,JK,JMOD) = SQRT(4._JPRB*ZKELV(JL,JK)**3 / &
!                     & (27._JPRB * ZKAPPA(JL) * PRDRY(JL,JK,JMOD)**3) )
!
!                !--- modified diffusivity and thermal conductivity
!                ZSQTERM(JL) =  SQRT(2._JPRB*RPI*ZAMW/(R*PT(JL,JK)))
!
!                !--- Abdul-Razzak et al. (1998) (Eq. 17):
!                ZDIFMOD(JL) = ZDIF(JL,JK) /   &
!                     & ( (PRDRY(JL,JK,JMOD)/(PRDRY(JL,JK,JMOD)+PPDELTA_V)) + &
!                     &   (ZDIF(JL,JK)/(PRDRY(JL,JK,JMOD)*PPALPHA_C)) * ZSQTERM(JL) )
!
!                !--- Abdul-Razzak et al. (1998) (Eq. 18):
!                ZKMOD(JL) = ZK(JL,JK) /   &
!                     & ( (PRDRY(JL,JK,JMOD)/(PRDRY(JL,JK,JMOD)+PPDELTA_T)) + &
!                     &   (ZK(JL,JK)/(PRDRY(JL,JK,JMOD)*PPALPHA_T*RCPD)) * ZSQTERM(JL) )
!             
!                !--- growth coefficient due to gas kinetic effects:
!          
!                !--- Abdul-Razzak et al. (1998) (Eq. 16):
!                ZTERM1(JL) = (PPRHO_WAT*R*PT(JL,JK)) / (ZESW(JL,JK)*ZDIFMOD(JL)*ZAMW)
!          
!                ZTERM2(JL) = (RLVTT*PPRHO_WAT) / (ZKMOD(JL)*PT(JL,JK))
!
!                ZTERM3(JL) = (RLVTT*ZAMW) / (R*PT(JL,JK))-1._JPRB
!
!                ZGROWTH(JL,JK,JMOD) = 1._JPRB / &
!                     & (ZTERM1(JL) + ZTERM2(JL) * ZTERM3(JL))
!             END IF
!          END DO
!       END DO
!
!    END DO
!       !<<dod
!
!    !--- Summation for equation (6):
!
!    DO JW=1,KPDF
!       DO JMOD=2,NSOL
!          DO JK=KTDIA,KLEV
!             DO JL=KIDIA,KFDIA
!                IF (ZN(JL,JK,JMOD) > ZEPS     .AND. &
!                  & PRDRY(JL,JK,JMOD)>1.E-9_JPRB .AND. &
!                  & LCLOUD(JL,JK) .AND. PW(JL,JK,JW) > ZEPS) THEN
!
!                   ! (10):                   
!                   ZXI = (2._JPRB*ZKELV(JL,JK)/3._JPRB) * &
!                       & SQRT(ZALPHA(JL,JK)*PW(JL,JK,JW)/ZGROWTH(JL,JK,JMOD))
!
!                   ! (11):
!                   ZETA = ((ZALPHA(JL,JK)*PW(JL,JK,JW)/ZGROWTH(JL,JK,JMOD))**1.5_JPRB) / &
!                        &  (2._JPRB*RPI*PPRHO_WAT*ZGAMMA(JL,JK)*ZN(JL,JK,JMOD))
!
!                   ! (6):
!                   ZSUM(JL,JK,JW) = ZSUM(JL,JK,JW) + & 
!                                & (1._JPRB/ZSM(JL,JK,JMOD)**2 *  &
!                                & (ZF(JMOD)*(ZXI/ZETA)**1.5_JPRB + &
!                                &  ZG(JMOD)*(ZSM(JL,JK,JMOD)**2._JPRB/(ZETA+3._JPRB*ZXI))**0.75_JPRB ) )
!                   
!                END IF
!
!             END DO
!          END DO
!       END DO
!
!       DO JK=KTDIA,KLEV
!          DO JL=KIDIA,KFDIA
!             IF (LCLOUD(JL,JK)) THEN
!                IF (ZSUM(JL,JK,JW) > ZEPS) THEN
!                   ZSMAX(JL,JK,JW) = SQRT(1._JPRB / ZSUM(JL,JK,JW))
!                END IF
!             END IF
!          END DO
!       END DO
!    END DO
!
!    !---Calculate activation:
!    DO JW=1,KPDF
!       DO JMOD=2, NSOL
!          DO JK=KTDIA, KLEV
!             DO JL=KIDIA,KFDIA
!
!                IF ( LCLOUD(JL,JK)         .AND. &
!                   & ZSMAX(JL,JK,JW)>ZEPS  .AND. &
!                   & ZSM(JL,JK,JMOD)>ZEPS  .AND. &
!                   & ZN(JL,JK,JMOD)>ZEPS   .AND. &
!                   & PAERONUM(JL,JK,JMOD)>1.E-9_JPRB ) THEN
!
!                   !---critical radius (12):
!                   ZRC=PRDRY(JL,JK,JMOD)*(ZSM(JL,JK,JMOD)/ZSMAX(JL,JK,JW))**(2._JPRB/3._JPRB)
!
!                   !---fraction of particles that are *larger* than Rc:
!                   ZT=(LOG(ZRC)-LOG(PRDRY(JL,JK,JMOD))) / (ZSQRT2*SIGMALN(JMOD))
!                   ZFRACN = 0.5_JPRB * ERFC(ZT)
!
!                   !---Sum up the total number of activated particles, integrating over updraft PDF [m-3]:
!                   ZFRACN_TOP(JL,JK,JMOD) = ZFRACN_TOP(JL,JK,JMOD) &
!                                           + ZFRACN*PWPDF(JL,JK,JW)
!                END IF
!
!                ZFRACN_BOT(JL,JK,JMOD) = ZFRACN_BOT(JL,JK,JMOD) + PWPDF(JL,JK,JW)
!
!             END DO ! jl
!          END DO ! jk
!       END DO
!    END DO
!
!    DO JMOD=2,NSOL
!       DO JK=KTDIA,KLEV
!          DO JL=KIDIA,KFDIA
!             
!             !---sum the total number of activated particles [# cm-3]
!             IF (LCLOUD(JL,JK)) THEN
!                PCDNC(JL,JK) = PCDNC(JL,JK) +  1.E-6_JPRB * ZN(JL,JK,JMOD) * &
!                     & ZFRACN_TOP(JL,JK,JMOD) / ZFRACN_BOT(JL,JK,JMOD)
!             END IF
!
!          END DO
!       END DO
!    END DO
!
!    IF (LHOOK) CALL DR_HOOK('YOE_AER_ACTIV.AER_ACTIV_ABDULRAZZAK_GHAN', 1, ZHOOK_HANDLE)
!
!  END SUBROUTINE AER_ACTIV_ABDULRAZZAK_GHAN
!
!  SUBROUTINE AER_ACTIV_MENON(KIDIA, KFDIA, KLON, KTDIA, KLEV, PT, PRHO, PLSM, PSO4MASS, PSSMASS, POMMASS, &
!                          &  PCDNC)
!
!    ! Description
!    ! -----------
!    ! aer_activ_menon calculates the CDNC concentration as a function of the bulk mass
!    ! concentrations of sulfate, sea salt and organics
!
!    ! Reference: 
!    ! ----------
!    ! GCM Simulations of the Aerosol Indirect Effect: Sensitivity to Cloud Parameterization
!    ! and Aerosol Burden
!    ! Menon et al., J. Atm. Sci. 59 692-713, 2002
!
!    !---inherited functions, types, variables and constants:
!    USE YOMHOOK,             ONLY: LHOOK, DR_HOOK
!    USE YOMCST,              ONLY: RPI
!    USE YOECLDP,             ONLY: RTHOMO, RCDNCSU, RCDNCSS, RCDNCOM
!
!    !---subroutine interface
!    INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA                    ! beginning of horizontal block
!    INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA                    ! beginning of horizontal block
!    INTEGER(KIND=JPIM), INTENT(IN) :: KLON                     ! horizontal dimension
!    INTEGER(KIND=JPIM), INTENT(IN) :: KTDIA                    ! highest level with liquid cloud
!    INTEGER(KIND=JPIM), INTENT(IN) :: KLEV                     ! number of model vertical levels
!    REAL(KIND=JPRB), INTENT(IN)    :: PT(KLON,KLEV)            ! air temperature [K]
!    REAL(KIND=JPRB), INTENT(IN)    :: PRHO(KLON,KLEV)          ! air density [kg m-3]
!    REAL(KIND=JPRB), INTENT(IN)    :: PLSM(KLON)               ! land-sea mask
!    REAL(KIND=JPRB), INTENT(IN)    :: PSO4MASS(KLON,KLEV)      ! so4 mass mixing ratio [kg/kg(air)]
!    REAL(KIND=JPRB), INTENT(IN)    :: PSSMASS(KLON,KLEV)       ! sea salt mass mixing ratio [kg/kg(air)]
!    REAL(KIND=JPRB), INTENT(IN)    :: POMMASS(KLON,KLEV)       ! OM mass mixing ratio [kg/kg(air)]
!    REAL(KIND=JPRB), INTENT(INOUT) :: PCDNC(KLON,KLEV)         ! number concentration of activated cloud droplets
!                                                               ! [#/cm3]
!
!    !---local types, variables and constants:
!    REAL(KIND=JPRB)    :: ZSO4MASS                             ! sulfate mass [ug m-3]
!    REAL(KIND=JPRB)    :: ZOCMASS                              ! organic mass [ug m-3]
!    REAL(KIND=JPRB)    :: ZSSMASS                              ! sea salt mass [ug m-3]
!    REAL(KIND=JPRB)    :: ZCDNC_LAND(KLON)                     ! CDNC according to 'over land' formulation in Menon et al.
!    REAL(KIND=JPRB)    :: ZCDNC_OCEAN(KLON)                    ! CDNC according to 'over ocean' formulation in Menon et al.
!    INTEGER(KIND=JPIM) :: JK, JL                               ! loop indices
!
!    REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!
!    !---executable procedure
!    IF (LHOOK) CALL DR_HOOK('YOE_AER_ACTIV.AER_ACTIV_MENON', 0, ZHOOK_HANDLE)
!
!    PCDNC(KIDIA:KFDIA,KTDIA:KLEV) = 0.0_JPRB
!    
!    DO JK=KTDIA,KLEV
!       DO JL=KIDIA,KFDIA
!
!          ! SO4, BC, OC mass in ug m-3:
!          ZSO4MASS = PSO4MASS(JL,JK)*PRHO(JL,JK)*1.E9_JPRB
!          ZSSMASS = PSSMASS(JL,JK)*PRHO(JL,JK)*1.E9_JPRB
!          ZOCMASS = POMMASS(JL,JK)*PRHO(JL,JK)*1.E9_JPRB
!
!          !---equation 1(b) in Menon et al., reformulated without log function.
!          !   Note: this gives CDNC in #/cm3
!          ZCDNC_LAND(JL) = 257._JPRB*(ZSO4MASS**RCDNCSU)*(ZOCMASS**RCDNCOM)
!          ZCDNC_OCEAN(JL) = ZCDNC_LAND(JL)*(ZSSMASS**RCDNCSS)
!                    
!       END DO
!
!       PCDNC(KIDIA:KFDIA,JK)=MERGE(ZCDNC_LAND(KIDIA:KFDIA), &
!       &   ZCDNC_OCEAN(KIDIA:KFDIA),PLSM(KIDIA:KFDIA)>0.5)
!
!    END DO
!
!    IF (LHOOK) CALL DR_HOOK('YOE_AER_ACTIV.AER_ACTIV_MENON', 1, ZHOOK_HANDLE)
!
!  END SUBROUTINE AER_ACTIV_MENON

  SUBROUTINE GET_HAMM7_AERO_PROP(KIDIA, KFDIA, KLON, KTDIA, KLEV, KSTGLO, LMODE, LBULK, &
                             & PAPH,     PGELAM,  PGEMU, PXTM1, KTRAC,                &
                             & PDRYRSOL, PAERONUM, &
                             & PSO4MASS, PBCMASS, POMMASS, PSSMASS, PDUMASS,   &
                             & PSO4BULK, PBCBULK, POMBULK, PSSBULK, PDUBULK,   &
                             & PNO3MASS, PMSAMASS)

    !---inherited functions, types, variables and constants
    USE YOMHOOK,             ONLY: LHOOK, DR_HOOK
    USE YOMCST,              ONLY: RPI 
    USE TM5M7_DATA,          ONLY: NSOL

    !<--eehol: add tracer variables for tracer indices to reduce hardcoding
    USE MO_HAM, ONLY:     &
         sizeclass,       & ! aerosol classes in HAM
         aerocomp           ! aerosol compounds by size class in HAM
    USE mo_ham_m7ctl,     ONLY: inucs,  iaits,  iaccs,  icoas,   &
                              iaiti,  iacci,  icoai,           &
                              iso4ns, iso4ks, iso4as, iso4cs,  &
                              ibcks,  ibcas,  ibccs,  ibcki,   &
                              iocks,  iocas,  ioccs,  iocki,   &
                              issas,  isscs,                   &
                              iduas,  iducs,  iduai,  iduci
    !-->eehol
    !USE YOERAD,              ONLY: LTM5AER, LCMIP6_PI_AEROSOLS, NRADFR
    !USE YOE_AERO_M7_DATA 
    !USE CPLNG,               ONLY: CPLNG_FLD, CPLNG_IDX, LEV_IDX_TM5, &
    !                             & NLEV_TM5, NLEV_TM5_AER 
    !USE YOE_PI_AERO         
    !USE YOMCT3,              ONLY: NSTEP

    IMPLICIT NONE

    !---subroutine interface
    !   *GET_HAMM7_AERO_PROP* is called here before the activation calculations
    !   
    !   INPUT:
    INTEGER(KIND=JPIM), INTENT(IN)         :: KIDIA   ! beginning of horizontal block
    INTEGER(KIND=JPIM), INTENT(IN)         :: KFDIA   ! end of horizontal block
    INTEGER(KIND=JPIM), INTENT(IN)         :: KLON    ! horizontal dimension
    INTEGER(KIND=JPIM), INTENT(IN)         :: KTDIA   ! highest level with liquid cloud
    INTEGER(KIND=JPIM), INTENT(IN)         :: KLEV    ! number of model vertical levels
    INTEGER(KIND=JPIM), INTENT(IN)         :: KSTGLO  ! offset of horizontal block in coupling arrays
    INTEGER(KIND=JPIM), INTENT(IN)         :: KTRAC   ! number of tracers

    LOGICAL, INTENT(IN) :: LMODE                      ! Per-mode data requested
    LOGICAL, INTENT(IN) :: LBULK                      ! Bulk aerosol masses requested
    
    REAL(KIND=JPRB), INTENT(IN)            :: PAPH(KLON,KLEV+1) ! half-level pressure
    REAL(KIND=JPRB), INTENT(IN)            :: PGELAM(KLON)      ! longitude
    REAL(KIND=JPRB), INTENT(IN)            :: PGEMU(KLON)       ! sine of latitude
    REAL(KIND=JPRB), INTENT(IN)            :: PXTM1(KLON,KLEV,KTRAC)       ! sine of latitude

    !   OUTPUT:
    REAL(KIND=JPRB), INTENT(OUT)           :: PDRYRSOL(KLON,KLEV,NSOL) ! [M]
    REAL(KIND=JPRB), INTENT(OUT)           :: PAERONUM(KLON,KLEV,NSOL) ! [#/KG(AIR)]
    REAL(KIND=JPRB), INTENT(OUT)           :: PSO4MASS(KLON,KLEV,NSOL) ! [KG(SO4)/KG(AIR)]
    REAL(KIND=JPRB), INTENT(OUT)           :: PBCMASS(KLON,KLEV,NSOL)  ! [KG(BC)/KG(AIR)]
    REAL(KIND=JPRB), INTENT(OUT)           :: POMMASS(KLON,KLEV,NSOL)  ! [KG(OM)/KG(AIR)]
    REAL(KIND=JPRB), INTENT(OUT)           :: PSSMASS(KLON,KLEV,NSOL)  ! [KG(SS)/KG(AIR)]
    REAL(KIND=JPRB), INTENT(OUT)           :: PDUMASS(KLON,KLEV,NSOL)  ! [KG(DU)/KG(AIR)]
    REAL(KIND=JPRB), INTENT(OUT)           :: PSO4BULK(KLON,KLEV)      ! [KG(SO4)/KG(AIR)]
    REAL(KIND=JPRB), INTENT(OUT)           :: PBCBULK(KLON,KLEV)       ! [KG(BC)/KG(AIR)]
    REAL(KIND=JPRB), INTENT(OUT)           :: POMBULK(KLON,KLEV)       ! [KG(OM)/KG(AIR)]
    REAL(KIND=JPRB), INTENT(OUT)           :: PSSBULK(KLON,KLEV)       ! [KG(SS)/KG(AIR)]
    REAL(KIND=JPRB), INTENT(OUT)           :: PDUBULK(KLON,KLEV)       ! [KG(DU)/KG(AIR)]
    REAL(KIND=JPRB), INTENT(OUT)           :: PNO3MASS(KLON,KLEV)      ! [KG(NO3)/KG(AIR)]
    REAL(KIND=JPRB), INTENT(OUT)           :: PMSAMASS(KLON,KLEV)      ! [KG(MSA)/KG(AIR)]

    !---local data
    REAL(KIND=JPRB)    :: NSO4, NH2SO4, NNACL, NNA, NCL, NNA2SO4 ! Particle numbers [kmol/kg air]
    REAL(KIND=JPRB)    :: ZDRYVOL2, ZDRYVOL3, ZDRYVOL4      ! Volume per particle [m3/#]
    REAL(KIND=JPRB)    :: ZDUMMY(KLON,KLEV)       ! For fields we don't use in subroutine call
    REAL(KIND=JPRB)    :: Z4PIOVER3
    REAL(KIND=JPRB)    :: ZHOOK_HANDLE
    INTEGER(KIND=JPIM) :: JL, JK, IBL, IL, IK

    !---executable procedure
    IF (LHOOK) CALL DR_HOOK('YOE_AER_ACTIV.GET_HAMM7_AERO_PROP', 0, ZHOOK_HANDLE)

!    IF (LCMIP6_PI_AEROSOLS) THEN
!
!       IBL = (KSTGLO-1)/KLON+1
!
!       IF (LPI_AERO_UPDATED .AND. MOD(NSTEP,NRADFR)==0) THEN
! 
!          !---CMIP6 preindustrial aerosols: interpolate from TM5 grid to physics grid
!          CALL CMIP6_PIAER_MXR_INTERP(KIDIA, KFDIA,  KLON, KLEV, 1, 0,   &
!                                 & PAPH,  PGELAM, PGEMU,             &
!                                 & PAERONUM(:,:,2),  PAERONUM(:,:,3), PAERONUM(:,:,4), &
!                                 & PSO4MASS(:,:,2),  PSO4MASS(:,:,3), PSO4MASS(:,:,4), &
!                                 & PBCMASS(:,:,2),   PBCMASS(:,:,3),  PBCMASS(:,:,4),  ZDUMMY, &
!                                 & POMMASS(:,:,2),   POMMASS(:,:,3),  POMMASS(:,:,4),  ZDUMMY, &
!                                 & PSSMASS(:,:,3),   PSSMASS(:,:,4), &
!                                 & PDUMASS(:,:,3),   PDUMASS(:,:,4),  ZDUMMY, ZDUMMY, &
!                                 & PNO3MASS(:,:),    PMSAMASS(:,:) ) 
!
!            PI_AERO_PHY(KIDIA:KFDIA,:,IPINUM_KS,IBL) = PAERONUM(KIDIA:KFDIA,:,JP_AITS)
!            PI_AERO_PHY(KIDIA:KFDIA,:,IPINUM_AS,IBL) = PAERONUM(KIDIA:KFDIA,:,JP_ACCS)
!            PI_AERO_PHY(KIDIA:KFDIA,:,IPINUM_CS,IBL) = PAERONUM(KIDIA:KFDIA,:,JP_COAS)
!            PI_AERO_PHY(KIDIA:KFDIA,:,IPIMSO4_KS,IBL) = PSO4MASS(KIDIA:KFDIA,:,JP_AITS)
!            PI_AERO_PHY(KIDIA:KFDIA,:,IPIMSO4_AS,IBL) = PSO4MASS(KIDIA:KFDIA,:,JP_ACCS)
!            PI_AERO_PHY(KIDIA:KFDIA,:,IPIMSO4_CS,IBL) = PSO4MASS(KIDIA:KFDIA,:,JP_COAS)
!            PI_AERO_PHY(KIDIA:KFDIA,:,IPIMBC_KS,IBL) = PBCMASS(KIDIA:KFDIA,:,JP_AITS)
!            PI_AERO_PHY(KIDIA:KFDIA,:,IPIMBC_AS,IBL) = PBCMASS(KIDIA:KFDIA,:,JP_ACCS)
!            PI_AERO_PHY(KIDIA:KFDIA,:,IPIMBC_CS,IBL) = PBCMASS(KIDIA:KFDIA,:,JP_COAS)
!            ! Note: insoluble modes are not presently used on the physics grid and 
!            ! are therefore not stored. Note also the third dimension of the PxxMASS
!            ! arrays is NSOL. 
!            !  PI_AERO_PHY(KIDIA:KFDIA,:,IPIMBC_KI,IBL) = PBCMASS((KIDIA:KFDIA,:,JP_AITI)
!            PI_AERO_PHY(KIDIA:KFDIA,:,IPIMOC_KS,IBL) = POMMASS(KIDIA:KFDIA,:,JP_AITS)
!            PI_AERO_PHY(KIDIA:KFDIA,:,IPIMOC_AS,IBL) = POMMASS(KIDIA:KFDIA,:,JP_ACCS)
!            PI_AERO_PHY(KIDIA:KFDIA,:,IPIMOC_CS,IBL) = POMMASS(KIDIA:KFDIA,:,JP_COAS)
!            ! PI_AERO_PHY(KIDIA:KFDIA,:,IPIMOC_KI,IBL) = POMMASS(KIDIA:KFDIA,:,JP_AITI)
!            PI_AERO_PHY(KIDIA:KFDIA,:,IPIMSS_AS,IBL) = PSSMASS(KIDIA:KFDIA,:,JP_ACCS)
!            PI_AERO_PHY(KIDIA:KFDIA,:,IPIMSS_CS,IBL) = PSSMASS(KIDIA:KFDIA,:,JP_COAS)
!            PI_AERO_PHY(KIDIA:KFDIA,:,IPIMDU_AS,IBL) = PDUMASS(KIDIA:KFDIA,:,JP_ACCS)
!            PI_AERO_PHY(KIDIA:KFDIA,:,IPIMDU_CS,IBL) = PDUMASS(KIDIA:KFDIA,:,JP_COAS)
!            ! PI_AERO_PHY(KIDIA:KFDIA,:,IPIMDU_AI,IBL) = PDUMASS(KIDIA:KFDIA,:,JP_ACCI)
!            ! PI_AERO_PHY(KIDIA:KFDIA,:,IPIMDU_CI,IBL) = PDUMASS(KIDIA:KFDIA,:,JP_COAI)
!            PI_AERO_PHY(KIDIA:KFDIA,:,IPIMNO3,IBL) = PNO3MASS(KIDIA:KFDIA,:)
!            PI_AERO_PHY(KIDIA:KFDIA,:,IPIMMSA,IBL) = PMSAMASS(KIDIA:KFDIA,:)
!       ELSE
!            PAERONUM(KIDIA:KFDIA,:,JP_AITS) = PI_AERO_PHY(KIDIA:KFDIA,:,IPINUM_KS,IBL)
!            PAERONUM(KIDIA:KFDIA,:,JP_ACCS) = PI_AERO_PHY(KIDIA:KFDIA,:,IPINUM_AS,IBL)
!            PAERONUM(KIDIA:KFDIA,:,JP_COAS) = PI_AERO_PHY(KIDIA:KFDIA,:,IPINUM_CS,IBL)
!            PSO4MASS(KIDIA:KFDIA,:,JP_AITS) = PI_AERO_PHY(KIDIA:KFDIA,:,IPIMSO4_KS,IBL) 
!            PSO4MASS(KIDIA:KFDIA,:,JP_ACCS) = PI_AERO_PHY(KIDIA:KFDIA,:,IPIMSO4_AS,IBL) 
!            PSO4MASS(KIDIA:KFDIA,:,JP_COAS) = PI_AERO_PHY(KIDIA:KFDIA,:,IPIMSO4_CS,IBL)
!            PBCMASS(KIDIA:KFDIA,:,JP_AITS) = PI_AERO_PHY(KIDIA:KFDIA,:,IPIMBC_KS,IBL)
!            PBCMASS(KIDIA:KFDIA,:,JP_ACCS) = PI_AERO_PHY(KIDIA:KFDIA,:,IPIMBC_AS,IBL)
!            PBCMASS(KIDIA:KFDIA,:,JP_COAS) = PI_AERO_PHY(KIDIA:KFDIA,:,IPIMBC_CS,IBL) 
!            ! PBCMASS(KIDIA:KFDIA,:,JP_AITI) = PI_AERO_PHY(KIDIA:KFDIA,:,IPIMBC_KI,IBL)
!            POMMASS(KIDIA:KFDIA,:,JP_AITS) = PI_AERO_PHY(KIDIA:KFDIA,:,IPIMOC_KS,IBL)
!            POMMASS(KIDIA:KFDIA,:,JP_ACCS) = PI_AERO_PHY(KIDIA:KFDIA,:,IPIMOC_AS,IBL)
!            POMMASS(KIDIA:KFDIA,:,JP_COAS) = PI_AERO_PHY(KIDIA:KFDIA,:,IPIMOC_CS,IBL)
!            ! POMMASS(KIDIA:KFDIA,:,JP_AITI) = PI_AERO_PHY(KIDIA:KFDIA,:,IPIMOC_KI,IBL)
!            PSSMASS(KIDIA:KFDIA,:,JP_ACCS) = PI_AERO_PHY(KIDIA:KFDIA,:,IPIMSS_AS,IBL)
!            PSSMASS(KIDIA:KFDIA,:,JP_COAS) = PI_AERO_PHY(KIDIA:KFDIA,:,IPIMSS_CS,IBL)
!            PDUMASS(KIDIA:KFDIA,:,JP_ACCS) = PI_AERO_PHY(KIDIA:KFDIA,:,IPIMDU_AS,IBL)
!            PDUMASS(KIDIA:KFDIA,:,JP_COAS) = PI_AERO_PHY(KIDIA:KFDIA,:,IPIMDU_CS,IBL)
!            ! PDUMASS(KIDIA:KFDIA,:,JP_ACCI) = PI_AERO_PHY(KIDIA:KFDIA,:,IPIMDU_AI,IBL)
!            ! PDUMASS(KIDIA:KFDIA,:,JP_COAI) = PI_AERO_PHY(KIDIA:KFDIA,:,IPIMDU_CI,IBL)
!            PNO3MASS(KIDIA:KFDIA,:) = PI_AERO_PHY(KIDIA:KFDIA,:,IPIMNO3,IBL)
!            PMSAMASS(KIDIA:KFDIA,:) = PI_AERO_PHY(KIDIA:KFDIA,:,IPIMMSA,IBL)
!       END IF
!    ELSE
!
!       IF (LTM5AER) THEN
!          IF (ICPLSU2 == -1) THEN
!             ICPLSU2 = CPLNG_IDX('A_SU2')
!             ICPLSU3 = CPLNG_IDX('A_SU3')
!             ICPLSU4 = CPLNG_IDX('A_SU4')
!             ICPLBC2 = CPLNG_IDX('A_BC2')
!             ICPLBC3 = CPLNG_IDX('A_BC3')
!             ICPLBC4 = CPLNG_IDX('A_BC4')
!             ICPLOM2 = CPLNG_IDX('A_OM2')
!             ICPLOM3 = CPLNG_IDX('A_OM3')
!             ICPLOM4 = CPLNG_IDX('A_OM4')
!             ICPLSS3 = CPLNG_IDX('A_SS3')
!             ICPLSS4 = CPLNG_IDX('A_SS4')
!             ICPLDD3 = CPLNG_IDX('A_DD3')
!             ICPLDD4 = CPLNG_IDX('A_DD4')
!             ICPLNO3 = CPLNG_IDX('A_NO3')
!             ICPLMSA = CPLNG_IDX('A_MSA') 
!          END IF
!          ! NO3 and MSA are not set in RADINTG
!          IF (ICPLNO3 == -1) ICPLNO3 = CPLNG_IDX('A_NO3')
!          IF (ICPLMSA == -1) ICPLMSA = CPLNG_IDX('A_MSA')
!
!          !---aerosol masses:
!          DO JK=KTDIA,KLEV
!
!             ! Include offset in case upper levels have not been received from TM5
!             IK = LEV_IDX_TM5(JK) - NLEV_TM5 + NLEV_TM5_AER
!             IF (IK >= 1) THEN
!                DO JL=KIDIA,KFDIA
!                   IL = KSTGLO+JL-1
!                   !---nucleation mode is ignored
!
!                   !---accumulation mode: SU, OM, BC
!                   PSO4MASS(JL,JK,JP_AITS) = CPLNG_FLD(ICPLSU2)%D(IL,IK,1)                   
!                   POMMASS(JL,JK,JP_AITS)  = CPLNG_FLD(ICPLOM2)%D(IL,IK,1) 
!                   PBCMASS(JL,JK,JP_AITS)  = CPLNG_FLD(ICPLBC2)%D(IL,IK,1) 
!
!                   !---accumulation mode: SU, OM, BC, SS, DU
!                   PSO4MASS(JL,JK,JP_ACCS) = CPLNG_FLD(ICPLSU3)%D(IL,IK,1)  
!                   POMMASS(JL,JK,JP_ACCS) = CPLNG_FLD(ICPLOM3)%D(IL,IK,1) 
!                   PBCMASS(JL,JK,JP_ACCS) = CPLNG_FLD(ICPLBC3)%D(IL,IK,1) 
!                   PSSMASS(JL,JK,JP_ACCS) = CPLNG_FLD(ICPLSS3)%D(IL,IK,1) 
!                   PDUMASS(JL,JK,JP_ACCS) = CPLNG_FLD(ICPLDD3)%D(IL,IK,1) 
!
!                   !---coarse mode: SU, OM, BC, SS, DU
!                   PSO4MASS(JL,JK,JP_COAS) = CPLNG_FLD(ICPLSU4)%D(IL,IK,1) 
!                   POMMASS(JL,JK,JP_COAS) = CPLNG_FLD(ICPLOM4)%D(IL,IK,1) 
!                   PBCMASS(JL,JK,JP_COAS) = CPLNG_FLD(ICPLBC4)%D(IL,IK,1) 
!                   PSSMASS(JL,JK,JP_COAS) = CPLNG_FLD(ICPLSS4)%D(IL,IK,1) 
!                   PDUMASS(JL,JK,JP_COAS) = CPLNG_FLD(ICPLDD4)%D(IL,IK,1) 
!
!                   !---nitrate and MSA
!                   PNO3MASS(JL,JK)        = CPLNG_FLD(ICPLNO3)%D(IL,IK,1)
!                   PMSAMASS(JL,JK)        = CPLNG_FLD(ICPLMSA)%D(IL,IK,1)
!                END DO
!             ELSE
!                ! These levels have not been received from TM5
!                ! and should be set to zero 
!                PSO4MASS(KIDIA:KFDIA,JK,:) = 0._JPRB
!                POMMASS(KIDIA:KFDIA,JK,:) = 0._JPRB
!                PBCMASS(KIDIA:KFDIA,JK,:) = 0._JPRB
!                PSSMASS(KIDIA:KFDIA,JK,:) = 0._JPRB
!                PDUMASS(KIDIA:KFDIA,JK,:) = 0._JPRB
!                PNO3MASS(KIDIA:KFDIA,JK)  = 0._JPRB
!                PMSAMASS(KIDIA:KFDIA,JK)  = 0._JPRB
!             END IF
!
!          END DO
! 
!          !---aerosol number and radius (if needed)
!          IF (LMODE) THEN
!             IF (ICPLN2 == -1) THEN
!                ICPLN2=CPLNG_IDX('A_N2')
!                ICPLN3=CPLNG_IDX('A_N3')
!                ICPLN4=CPLNG_IDX('A_N4')
!             END IF
!             DO JK=KTDIA,KLEV
!                ! Include offset in case upper levels have not been received from TM5
!                IK = LEV_IDX_TM5(JK) - NLEV_TM5 + NLEV_TM5_AER
!                IF (IK >= 1) THEN
!                   DO JL=KIDIA,KFDIA
!                      IL = KSTGLO+JL-1
!                      PAERONUM(JL,JK,JP_AITS) = CPLNG_FLD(ICPLN2)%D(IL,IK,1)
!                      PAERONUM(JL,JK,JP_ACCS) = CPLNG_FLD(ICPLN3)%D(IL,IK,1)
!                      PAERONUM(JL,JK,JP_COAS) = CPLNG_FLD(ICPLN4)%D(IL,IK,1)
!                   END DO
!                ELSE
!                   ! These levels have not been received from TM5
!                   ! and should be set to zero
!                   PAERONUM(KIDIA:KFDIA,JK,:) = TINY(0._JPRB)
!                END IF
!             END DO
!          END IF
!       END IF          ! ltm5aer
!    END IF             ! lcmip6_pi_aerosols
!
!    !---aerosol radius (if needed) for both interactive TM5 and preindustrial climatology
!    IF (LMODE) THEN
!       Z4PIOVER3 = 4._JPRB*RPI/3._JPRB
!       DO JK=KTDIA,KLEV
!          DO JL=KIDIA,KFDIA
!             ! Aitken mode 
!             NSO4 = PSO4MASS(JL,JK,JP_AITS) / WSO4
!             NH2SO4 = NSO4
!             ! dry volume [m3]
!             ZDRYVOL2  = &
!                  & (NH2SO4*WH2SO4 / (DH2SO4 * 1.0E3_JPRB) + &
!                  &  POMMASS(JL,JK,JP_AITS) / (DOC * 1.0E3_JPRB) + & 
!                  &  PBCMASS(JL,JK,JP_AITS) / (DBC * 1.0E3_JPRB)) / &
!                  &  PAERONUM(JL,JK,JP_AITS) 
!             ! dry radius [m]
!             !PDRYRSOL(JL,JK,JP_AITS) = RAM2CMR(JP_AITS)* &
!             !                      &  (ZDRYVOL2 / Z4PIOVER3)**(1._JPRB/3._JPRB) 
!
!             ! Accumulation mode
!             NNA = PSSMASS(JL,JK,JP_ACCS) / WNACL
!             NCL = NNA
!             NSO4 = PSO4MASS(JL,JK,JP_ACCS) / WSO4
!             NNA2SO4 = MIN(NNA/2._JPRB, NSO4)
!             NNA = NNA - 2._JPRB*NNA2SO4
!             NNACL = MIN(NCL, NNA)
!             NCL = NNACL
!             NH2SO4 = NSO4 - NNA2SO4
!             ! dry volume [m3]
!             ZDRYVOL3  = &
!                  & (NNACL*WNACL / (DNACL * 1.0E3_JPRB) + &
!                  &  NNA2SO4*WNA2SO4 / (DNA2SO4 * 1.0E3_JPRB) + &
!                  &  NH2SO4*WH2SO4/ (DH2SO4 * 1.0E3_JPRB) + &
!                  &  POMMASS(JL,JK,JP_ACCS) / (DOC * 1.0E3_JPRB) + &
!                  &  PBCMASS(JL,JK,JP_ACCS) / (DBC * 1.0E3_JPRB) + &
!                  &  PDUMASS(JL,JK,JP_ACCS) / (DDUST * 1.0E3_JPRB) + &
!                  &  PNO3MASS(JL,JK) * NH4NO3_FACTOR / (DNH4NO3 * 1.0E3_JPRB) + &
!                  &  PMSAMASS(JL,JK) / (DMSA * 1.0E3_JPRB)) / &
!                  &  PAERONUM(JL,JK,JP_ACCS) 
!             ! dry radius [m]
!             !PDRYRSOL(JL,JK,JP_ACCS) = RAM2CMR(JP_ACCS)* &
!             !                      &  (ZDRYVOL3 / Z4PIOVER3)**(1._JPRB/3._JPRB) 
!
!             ! Coarse mode
!             NNA = PSSMASS(JL,JK,JP_COAS) / WNACL
!             NCL = NNA
!             NSO4 = PSO4MASS(JL,JK,JP_COAS) / WSO4
!             NNA2SO4 = MIN(NNA/2._JPRB, NSO4)
!             NNA = NNA - 2._JPRB*NNA2SO4
!             NNACL = MIN(NCL, NNA)
!             NCL = NNACL
!             NH2SO4 = NSO4 - NNA2SO4
!             ! dry volume [m3]
!             ZDRYVOL4  = &
!                  & (NNACL*WNACL / (DNACL * 1.0E3_JPRB) + &
!                  &  NNA2SO4*WNA2SO4 / (DNA2SO4 * 1.0E3_JPRB) + & 
!                  &  NH2SO4*WH2SO4/ (DH2SO4 * 1.0E3_JPRB) + &
!                  &  POMMASS(JL,JK,JP_COAS) / (DOC * 1.0E3_JPRB) + &
!                  &  PBCMASS(JL,JK,JP_COAS) / (DBC * 1.0E3_JPRB) + &
!                  &  PDUMASS(JL,JK,JP_COAS) / (DDUST * 1.0E3_JPRB)) / &
!                  &  PAERONUM(JL,JK,JP_COAS) 
!             ! dry radius [m]
!             !PDRYRSOL(JL,JK,JP_COAS) = RAM2CMR(JP_COAS)* &
!             !                      &  (ZDRYVOL4 / Z4PIOVER3)**(1._JPRB/3._JPRB) 
!         END DO
!       END DO
!    END IF
    
    !eehol: add separation of tracers NOTE: These PXTM1 indices need to go according to HAM indices and not OIFS!!
    IF (LMODE) THEN
      DO JK=KTDIA,KLEV
         DO JL=KIDIA,KFDIA
            
            !---nucleation mode is ignored

            !---Aitken soluble mode: SU, OM, BC
            PSO4MASS(JL,JK,2) = PXTM1(JL,JK,aerocomp(iso4ks)%idt)  !SO4 Ait sol
            POMMASS(JL,JK,2)  = PXTM1(JL,JK,aerocomp(iocks)%idt) !OC Ait sol
            PBCMASS(JL,JK,2)  = PXTM1(JL,JK,aerocomp(ibcks)%idt)  !BC Ait sol
            
            !PSO4MASS(JL,JK,2) = PXTM1(JL,JK,5)  !SO4 Ait sol
            !POMMASS(JL,JK,2)  = PXTM1(JL,JK,12) !OC Ait sol
            !PBCMASS(JL,JK,2)  = PXTM1(JL,JK,8)  !BC Ait sol

            !---accumulation soluble mode: SU, OM, BC, SS, DU
            PSO4MASS(JL,JK,3) = PXTM1(JL,JK,aerocomp(iso4as)%idt) !SO4 acc sol
            POMMASS(JL,JK,3) = PXTM1(JL,JK,aerocomp(iocas)%idt) !OC acc sol
            PBCMASS(JL,JK,3) = PXTM1(JL,JK,aerocomp(ibcas)%idt)  !BC acc sol
            PSSMASS(JL,JK,3) = PXTM1(JL,JK,aerocomp(issas)%idt) !SS acc sol
            PDUMASS(JL,JK,3) = PXTM1(JL,JK,aerocomp(iduas)%idt) !DU acc sol
            
            !PSO4MASS(JL,JK,3) = PXTM1(JL,JK,6) !SO4 acc sol
            !POMMASS(JL,JK,3) = PXTM1(JL,JK,13) !OC acc sol
            !PBCMASS(JL,JK,3) = PXTM1(JL,JK,9)  !BC acc sol
            !PSSMASS(JL,JK,3) = PXTM1(JL,JK,16) !SS acc sol
            !PDUMASS(JL,JK,3) = PXTM1(JL,JK,18) !DU acc sol

            !---coarse soluble mode: SU, OM, BC, SS, DU
            PSO4MASS(JL,JK,4) = PXTM1(JL,JK,aerocomp(iso4cs)%idt) !SO4 coa sol
            POMMASS(JL,JK,4) = PXTM1(JL,JK,aerocomp(ioccs)%idt) !OC coa sol
            PBCMASS(JL,JK,4) = PXTM1(JL,JK,aerocomp(ibccs)%idt) !BC coa sol
            PSSMASS(JL,JK,4) = PXTM1(JL,JK,aerocomp(isscs)%idt) !SS coa sol
            PDUMASS(JL,JK,4) = PXTM1(JL,JK,aerocomp(iducs)%idt) !DU coa sol

            !PSO4MASS(JL,JK,4) = PXTM1(JL,JK,7) !SO4 coa sol
            !POMMASS(JL,JK,4) = PXTM1(JL,JK,14) !OC coa sol
            !PBCMASS(JL,JK,4) = PXTM1(JL,JK,10) !BC coa sol
            !PSSMASS(JL,JK,4) = PXTM1(JL,JK,17) !SS coa sol
            !PDUMASS(JL,JK,4) = PXTM1(JL,JK,19) !DU coa sol

            !---nitrate and MSA eehol: Currently not used in HAM... need to be added later!
            PNO3MASS(JL,JK)        = 0._JPRB !CPLNG_FLD(ICPLNO3)%D(IL,IK,1)
            PMSAMASS(JL,JK)        = 0._JPRB !CPLNG_FLD(ICPLMSA)%D(IL,IK,1)

            PAERONUM(JL,JK,2) = PXTM1(JL,JK,sizeclass(iaits)%idt_no) !Ait sol
            PAERONUM(JL,JK,3) = PXTM1(JL,JK,sizeclass(iaccs)%idt_no) !acc sol
            PAERONUM(JL,JK,4) = PXTM1(JL,JK,sizeclass(icoas)%idt_no) !coa sol

            !PAERONUM(JL,JK,2) = PXTM1(JL,JK,23) !Ait sol
            !PAERONUM(JL,JK,3) = PXTM1(JL,JK,24) !acc sol
            !PAERONUM(JL,JK,4) = PXTM1(JL,JK,25) !coa sol
         END DO
      END DO
    END IF

    ! bulk masses, if requested
    IF (LBULK) THEN
       DO JK=KTDIA,KLEV
          DO JL=KIDIA,KFDIA
             PSO4BULK(JL,JK) = PSO4MASS(JL,JK,2) + &
                             & PSO4MASS(JL,JK,3) + &
                             & PSO4MASS(JL,JK,4)
             POMBULK(JL,JK) = POMMASS(JL,JK,2) + &
                            & POMMASS(JL,JK,3) + &
                            & POMMASS(JL,JK,4)
             PBCBULK(JL,JK) = PBCMASS(JL,JK,2) + &
                            & PBCMASS(JL,JK,3) + &
                            & PBCMASS(JL,JK,4)
             PSSBULK(JL,JK) = PSSMASS(JL,JK,3) + &
                            & PSSMASS(JL,JK,4) 
             PDUBULK(JL,JK) = PDUMASS(JL,JK,3) + &
                            & PDUMASS(JL,JK,4)
          END DO
       END DO
    END IF

  IF (LHOOK) CALL DR_HOOK('YOE_AER_ACTIV.GET_HAMM7_AERO_PROP', 1, ZHOOK_HANDLE)

  END SUBROUTINE GET_HAMM7_AERO_PROP


!  SUBROUTINE DIAGNOSE_AEROSOL_MASS ( & 
!       !  input
!       & KIDIA,    KFDIA,    KLON,   KLEV,&
!       & PT,       PQ,       PQSAT, &
!       & PAPH,     PAP, &
!       & PGELAM,   PGEMU,    PCLON,   PSLON,&
!       !  output
!       & PSO4MASS, PBCMASS, POMMASS, PSSMASS, PDUMASS)
!
!    !     DIAGNOSE_AEROSOL_MASS
!    !     ---------------------
!    !     Diagnose the bulk mass mixing ratios of sulfate, black carbon, organic matter,
!    !     sea salt and mineral dust from the optical properties given by the Tegen climatology
!    !          effects on clouds and convection
!
!    !     AUTHOR
!    !          A. Tompkins  E.C.M.W.F. (aer_clcld.F90, on which this is based)
!    !          D. O'Donnell, FMI
!    !     PURPOSE.
!    !     --------
!    
!    !     INTERFACE.
!    !     ----------
!    
!    !     *DIAGNOSE_AEROSOL_MASS* IS CALLED FROM *MOD_AER_ACTIV.MAIN*
!    
!    !     PARAMETER     DESCRIPTION                                   UNITS
!    !     ---------     -----------                                   -----
!    
!    ! -   INPUT ARGUMENTS.
!    !     -------------------
!    
!    ! KIDIA   : START OF HORIZONTAL LOOP
!    ! KFDIA   : END   OF HORIZONTAL LOOP
!    ! KLON    : HORIZONTAL DIMENSION
!    ! KLEV    : END OF VERTICAL LOOP AND VERTICAL DIMENSION
!    ! PGELAM     : LONGITUDE
!    ! PCLON      : COSINE OF LONGITUDE
!    ! PSLON      : SINE   OF LONGITUDE
!    ! PGEMU      : SINE OF LATITUDE
!    
!    ! -   OUTPUT ARGUMENTS.
!    !     -------------------
!    ! PSSMASS   : mass mixing ratio of seasalt aerosol 
!    ! PSO4MASS  : mass mixing ratio of sulfate aerosol 
!    ! POMMASS   : mass mixing ratio of organic aerosol 
!    ! PBCMASS   : mass mixing ratio of black carbon aerosol 
!    ! PDUMASS   : mass mixing ratio of mineral dust aerosol 
!    
!    USE PARKIND1  ,ONLY : JPIM     ,JPRB
!    USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
!    
!    USE YOMCST   , ONLY : RG       ,RD       ,RETV     ,&
!         & RLVTT    ,RTT     ,RPI, RLSTT
!    USE YOECLDP  , ONLY : RCLCRIT  ,RLCRITSNOW, RCLDMAX,&
!         & RNICE , LAERLIQAUTOLSP, LAERLIQCOLL, LAERICESED, LAERICEAUTO    
!    USE YOEAEROP , ONLY : ALF_SS, ALF_SU, ALF_BC, ALF_DD, ALF_OM      !>>dod <<
!    USE YOEAERSNK, ONLY : RRHTAB
!    USE YOERAD   , ONLY : LEPO3RA
!    USE YOE_AERO_M7_DATA, ONLY: NSOL
!    USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
!         & R4IES    ,R5LES    ,R5IES    ,R5ALVCP  ,R5ALSCP  ,&
!         & RALVDCP  ,RALSDCP  ,RTWAT    ,&
!         & RTICE    ,RTICECU  ,&
!         & RTWAT_RTICE_R      ,RTWAT_RTICECU_R,&
!         & RKOOP1   ,RKOOP2
!
!    IMPLICIT NONE
!    
!    ! input variables
!    INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
!    INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
!    INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
!    INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
!    REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1) 
!    REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV) 
!    REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV) 
!    REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ(KLON,KLEV) 
!    REAL(KIND=JPRB)   ,INTENT(IN)    :: PQSAT(KLON,KLEV) 
!    REAL(KIND=JPRB)   ,INTENT(IN)    :: PGELAM(KLON) 
!    REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLON(KLON) 
!    REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLON(KLON) 
!    REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEMU(KLON) 
!    
!    ! output
!    REAL(KIND=JPRB),   INTENT(OUT)   :: PSSMASS(KLON,KLEV)
!    REAL(KIND=JPRB),   INTENT(OUT)   :: PSO4MASS(KLON,KLEV)
!    REAL(KIND=JPRB),   INTENT(OUT)   :: POMMASS(KLON,KLEV)
!    REAL(KIND=JPRB),   INTENT(OUT)   :: PBCMASS(KLON,KLEV)
!    REAL(KIND=JPRB),   INTENT(OUT)   :: PDUMASS(KLON,KLEV)
!        
!    !----------------------------------------------------------------------
!
!    ! local arrays
!
!    ! Aerosol arrays
!    REAL(KIND=JPRB) :: ZAOD(KLON,6,KLEV)   ! aerosol optical depth
!    REAL(KIND=JPRB) :: ZOZONE(KLON,KLEV)   ! O3 concentration, not currently used
!    REAL(KIND=JPRB) :: ZMAERMN(6)          ! annual column mean mass of aerosol
!    REAL(KIND=JPRB) :: ZECPO3(KLON,KLEV)   ! dummy prognostic ozone array
!    REAL(KIND=JPRB) :: ZQS(KLON,KLEV)      ! saturation
!
!    REAL(KIND=JPRB) :: ZS0, ZSCRITHOMO, ZSVP, ZTEMPC
!    REAL(KIND=JPRB) :: ZNCRIT_GIERENS, ZNCRIT_REN
!    REAL(KIND=JPRB) :: ZNICEHOMO
!    !REAL(KIND=JPRB) :: ZLIQCLD, ZICERE
!    
!    REAL(KIND=JPRB) :: ZCLD
!    REAL(KIND=JPRB) :: ZRLIQ_CRIT, ZRICE_CRIT ! critical radii for autoconversion process
!    
!    ! for RH look up tables
!    INTEGER(KIND=JPIM) :: IRH(KLON,KLEV)
!    INTEGER(KIND=JPIM) :: JTYP, JTAB, IBIN
!    INTEGER(KIND=JPIM) :: JAERSS, JAERDU, JAEROM, JAERSU, JAERBC
!    
!    ! these are reduced compared to USE YOEAEROP, since 1 band only
!    !REAL(KIND=JPRB) :: ZALF_BC(1)     
!    !REAL(KIND=JPRB) :: ZALF_DD(3)   
!    !REAL(KIND=JPRB) :: ZALF_OM(12)  
!    !REAL(KIND=JPRB) :: ZALF_SS(12,3)
!    !REAL(KIND=JPRB) :: ZALF_SU(12)  
!    !REAL(KIND=JPRB) :: ZRHTAB(12)  
!    
!    REAL(KIND=JPRB) :: ZALF, ZRH, ZWTOT, ZALF_SU, ZALF_BC, ZALF_OM, ZALF_SS, ZALF_DU
!    
!    ! general arrays
!    REAL(KIND=JPRB) :: ZTHF(KLON,KLEV+1)   ! T on half levels
!    
!    REAL(KIND=JPRB) :: ZDPR
!    
!    ! misc variables
!    REAL(KIND=JPRB) :: ZEPSEC
!    INTEGER(KIND=JPIM) :: IWAVL
!    INTEGER(KIND=JPIM) :: JK, JL
!    
!    REAL(KIND=JPRB) :: ZHOOK_HANDLE
!    
!    !------------------------
!    ! interface include files
!    !------------------------
!#include "radact.intfb.h"
!    
!#include "fcttre.h"
!
!    !--------------------------------------------------------------------------
!    
!    IF (LHOOK) CALL DR_HOOK('YOE_AER_ACTIV.DIAGNOSE_AEROSOL_MASS',0,ZHOOK_HANDLE)
!    
!    !######################################################################
!    !                       1.0 Basic variables
!    !######################################################################
!    
!    ZEPSEC=1.E-10_JPRB
!    IWAVL=8               ! reference to 550 nm
!    
!    ! move to cldp module
!    ZRLIQ_CRIT=9.3E-6_JPRB  ! cloud to rain critical radius
!    ZWTOT=0.1_JPRB ! governs critical N
!    
!    !---initialization
!    PSSMASS(:,:) = 0.0_JPRB
!    PSO4MASS(:,:) = 0.0_JPRB
!    POMMASS(:,:) = 0.0_JPRB
!    PBCMASS(:,:) = 0.0_JPRB
!    PDUMASS(:,:) = 0.0_JPRB
!    
!    ! change to vector vmass 
!    DO JK=1,KLEV
!       DO JL=KIDIA,KFDIA
!          ZQS(JL,JK)=FOEEWM(PT(JL,JK))/PAP(JL,JK)
!          ZQS(JL,JK)=MIN(0.5_JPRB,ZQS(JL,JK))
!          ZQS(JL,JK)=ZQS(JL,JK)/(1.0_JPRB-RETV*ZQS(JL,JK))
!       ENDDO
!    ENDDO
!    DO JK=2,KLEV
!       DO JL=KIDIA,KFDIA
!          ZTHF(JL,JK)=(PT(JL,JK-1)*PAP(JL,JK-1)&
!               & *(PAP(JL,JK)-PAPH(JL,JK))&
!               & +PT(JL,JK)*PAP(JL,JK)*(PAPH(JL,JK)-PAP(JL,JK-1)))&
!               & *(1.0_JPRB/(PAPH(JL,JK)*(PAP(JL,JK)-PAP(JL,JK-1))))  
!       ENDDO
!    ENDDO
!    
!    DO JL=KIDIA,KFDIA
!       ZTHF(JL,KLEV+1)=PT(JL,KLEV) ! should be surface temperature
!       ZTHF(JL,1)=PT(JL,1)
!    ENDDO
!    
!    IF (LEPO3RA) THEN
!       DO JK=1,KLEV
!          DO JL=KIDIA,KFDIA
!             ZECPO3(JL,JK)=0.0_JPRB
!          ENDDO
!       ENDDO
!    ENDIF
!    
!    !######################################################################
!    !               2. Retrieve aerosols optical depths 
!    !######################################################################
!    ! line 2: KRINT=1, KDLON=KLON , P2=KLON, KSHIFT=0, 
!    ! line 4: ozone set to dummy variable.
!    
!    CALL RADACT ( KIDIA , KFDIA, KLON , KLEV,&
!         & 1    , KLON  , KLON , 0    , 1   ,&
!         & PAPH , &
!         & PGELAM, PGEMU, PCLON, PSLON, ZTHF,&
!         & PQ   , PQSAT , ZECPO3,&
!         & ZAOD, ZOZONE  )  
!        
!
!    !######################################################################
!    !        3. Retrieve Aerosol mass from Tau
!    !######################################################################
!    
!    !     RADACT order:
!    !     1=sulfate + organics
!    !     2=sea salt
!    !     3=black carbon
!    !     4=mineral dust
!    !     5=Volcanic
!    !     6=Background
!    ! set up indices
!    JAERSU=1
!    JAERSS=2
!    JAERDU=3
!    JAERBC=4
!    JAEROM=1
!    
!    
!    !========================================
!    ! conversion of aerosols from Tau to Mass 
!    !========================================
!    !-- define RH index from "clear-sky" (not yet!) relative humidity
!    
!    ! for now fix to using Bin 1, meaning the smallest SS particles
!    IBIN=1
!    
!    ! taken from YOEAEROP: 17- band data for IWAVL=8 corresponding to 550nm
!    
!    DO JK=1,KLEV
!       DO JL=KIDIA,KFDIA
!          ZRH=100.0_JPRB*PQ(JL,JK)/PQSAT(JL,JK)
!          ZRH=MIN(MAX(ZRH,1.0_JPRB),100.0_JPRB)
!          !>>dod
!          IRH(JL,JK) = 1
!          !<<dod
!          DO JTAB=1,12
!             IF (ZRH > RRHTAB(JTAB)) THEN
!                IRH(JL,JK)=JTAB
!             ENDIF
!          ENDDO
!       ENDDO
!    ENDDO
!    
!    
!    ZALF_BC=ALF_BC(IWAVL)
!    ZALF_DU=ALF_DD(IBIN,IWAVL)
!    
!    DO JK=1,KLEV
!       DO JL=KIDIA,KFDIA
!          
!          ZALF_SU=ALF_SU(IRH(JL,JK),IWAVL)
!          ZALF_OM=ALF_OM(IRH(JL,JK),IWAVL)
!          ZALF_SS=ALF_SS(IRH(JL,JK),IWAVL,IBIN)
!          
!          ZDPR=PAPH(JL,JK+1)-PAPH(JL,JK)
!
!          !---sulfate and organics: allocate 50% to the sulfate mass and 50% to organics
!          IF (ZALF_SU /= 0.0_JPRB .AND. ZDPR /=0.0_JPRB ) THEN 
!             ! PSO4MASS(JL,JK) = ZAOD(JL,JAERSU,JK)*RG/(ZDPR*ZALF_SU*1000._JPRB)
!             PSO4MASS(JL,JK) = 0.5_JPRB*ZAOD(JL,JAERSU,JK)*RG/(ZDPR*ZALF_SU*1000._JPRB)
!          END IF
!          
!          IF (ZALF_OM /= 0.0_JPRB .AND. ZDPR /=0.0_JPRB ) THEN 
!             ! POMMASS(JL,JK) = 0._JPRB
!             POMMASS(JL,JK) = 0.5_JPRB*ZAOD(JL,JAERSU,JK)*RG/(ZDPR*ZALF_OM*1000._JPRB)
!          END IF
!      
!          !---black carbon
!          IF (ZALF_BC /= 0.0_JPRB .AND. ZDPR /=0.0_JPRB ) THEN 
!             PBCMASS(JL,JK) = ZAOD(JL,JAERBC,JK)*RG/(ZDPR*ZALF_BC*1000._JPRB)
!          END IF
!
!          !---sea salt
!          IF (ZALF_SS /= 0.0_JPRB .AND. ZDPR /=0.0_JPRB ) THEN 
!             PSSMASS(JL,JK) = ZAOD(JL,JAERSS,JK)*RG/(ZDPR*ZALF_SS*1000._JPRB)
!          END IF
!          
!          !---mineral dust
!          IF (ZALF_DU /= 0.0_JPRB .AND. ZDPR /=0.0_JPRB ) THEN 
!             PDUMASS(JL,JK) = ZAOD(JL,JAERDU,JK)*RG/(ZDPR*ZALF_DU*1000._JPRB)
!          END IF
!
!       END DO
!    END DO
!
!    IF (LHOOK) CALL DR_HOOK('YOE_AER_ACTIV.DIAGNOSE_AEROSOL_MASS', 1, ZHOOK_HANDLE) 
!
!  END SUBROUTINE DIAGNOSE_AEROSOL_MASS

!  SUBROUTINE LIQ_CLOUD_RE(KIDIA, KFDIA, KLON, KTDIA, KLEV, LCLOUD, PL, PA, PRHO, PGFL)
!
!    ! *LIQ_CLOUD_RE* calculates the effective radius (Re) for liquid clouds
!
!    !---inherited functions, types, variables and constants
!    USE YOMHOOK,             ONLY: LHOOK, DR_HOOK
!    USE YOMCST,              ONLY: RPI
!    USE YOECLDP,             ONLY: RCLDMAX, PPRHO_WAT
!    USE YOM_YGFL,            ONLY: YGFL, YCDNC, YRE_LIQ
!    
!    IMPLICIT NONE
!
!    !---subroutine interface
!    !   Input:
!    INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
!    INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
!    INTEGER(KIND=JPIM), INTENT(IN) :: KLON
!    INTEGER(KIND=JPIM), INTENT(IN) :: KTDIA
!    INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
!
!    LOGICAL, INTENT(IN)            :: LCLOUD(KLON,KLEV)
!
!    REAL(KIND=JPRB), INTENT(IN)    :: PL(KLON,KLEV)
!    REAL(KIND=JPRB), INTENT(IN)    :: PA(KLON,KLEV)
!    REAL(KIND=JPRB), INTENT(IN)    :: PRHO(KLON,KLEV)
!    
!    !   Input and output
!    REAL(KIND=JPRB), INTENT(INOUT) :: PGFL(KLON,KLEV,YGFL%NDIM)
!
!    !---local data:
!    REAL(KIND=JPRB)    :: ZCLD, ZRE_LIQ(KLON,KLEV)
!    REAL(KIND=JPRB)    :: ZEPSEC
!    REAL(KIND=JPRB)    :: ZHOOK_HANDLE
!
!    INTEGER(KIND=JPIM) :: JL,JK
!
!    !---executable procedure
!    IF (LHOOK) CALL DR_HOOK('YOE_AER_ACTIV.LIQ_CLOUD_RE', 0, ZHOOK_HANDLE)
!
!    ZEPSEC = 1.E-10_JPRB
!
!! this is effective radius calculation
!    DO JK=KTDIA,KLEV
!       DO JL=KIDIA,KFDIA
!          ZCLD=PL(JL,JK)/MAX(PA(JL,JK),ZEPSEC)
!          ZCLD=MIN(MAX(ZCLD,0.0_JPRB),RCLDMAX)
!
!          ! 2.387e-10 is 3/(4*pi*rho_liq*10^6)  [10^6 for N in right units]
!          ! PRE_LIQ(JL,JK)=(2.387e-10_JPRB*PRHO(JL,JK)*ZCLD/PCDNC(JL,JK))**0.333_JPRB
!
!          ZRE_LIQ(JL,JK) = (0.75_JPRB*PRHO(JL,JK)*ZCLD / &
!                  & (RPI*PPRHO_WAT*1.E6_JPRB*PGFL(JL,JK,YCDNC%MP9_PH)))**0.333_JPRB
!       END DO
!    END DO
!
!    ZRE_LIQ = ZRE_LIQ*1.E6_JPRB
!
!    ! This taken from old radlswr:-
!    ! Limit effective radius to within defined range
!    ZRE_LIQ = MAX(ZRE_LIQ, 4.0_JPRB)
!    ZRE_LIQ = MIN(ZRE_LIQ,30.0_JPRB)
!
!    ! Set R_eff_liq only where there are clouds
!    PGFL(KIDIA:KFDIA,KTDIA:KLEV,YRE_LIQ%MP9_PH) = &
!    &  MERGE( ZRE_LIQ(KIDIA:KFDIA,KTDIA:KLEV),  &
!    &         PGFL(KIDIA:KFDIA,KTDIA:KLEV,YRE_LIQ%MP9_PH), &
!    &         LCLOUD(KIDIA:KFDIA,KTDIA:KLEV) )
!
!    IF (LHOOK) CALL DR_HOOK('YOE_AER_ACTIV.LIQ_CLOUD_RE', 1, ZHOOK_HANDLE)
!
!  END SUBROUTINE LIQ_CLOUD_RE

  SUBROUTINE ICE_CLOUD_PROP(KIDIA, KFDIA, KLON, KLEV, PT, PRHO, PI, PA, PAP, &
                         &  PQSAT, PSO4MASS, PBCMASS, PDUMASS, PGFL, YDMODEL, PRE_ICE, PICNC)

    ! This subroutine is mainly a copy/paste of the ice phase microphysics 
    ! implemented in the subroutine aer_cld.F90 
 
    !---inherited functions, types, variables and constants
    USE YOMHOOK,             ONLY: LHOOK, DR_HOOK
    USE YOMCST,              ONLY: RPI, RTT, RETV, RLSTT, RLVTT
    !USE YOECLDP,             ONLY: RNICE, PPRHO_ICE, RCLDMAX
    !USE YOE_AERO_M7_DATA,    ONLY: NSOL
    USE TM5M7_DATA,    ONLY: NSOL
    !USE YOM_YGFL,            ONLY: YGFL, YICNC, YRE_ICE
    USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
         & R4IES    ,R5LES    ,R5IES    ,R5ALVCP  ,R5ALSCP  ,&
         & RALVDCP  ,RALSDCP  ,RTWAT    ,&
         & RTICE    ,RTICECU  ,&
         & RTWAT_RTICE_R      ,RTWAT_RTICECU_R,&
         & RKOOP1   ,RKOOP2
    USE TYPE_MODEL,          ONLY: MODEL

    IMPLICIT NONE

    !---subroutine interface
    !   *ICE_CLD_PROP* is called from AER_ACTIV 
    !   
    !   Input:
    TYPE(MODEL),        INTENT(IN) :: YDMODEL
    INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
    INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
    INTEGER(KIND=JPIM), INTENT(IN) :: KLON
    INTEGER(KIND=JPIM), INTENT(IN) :: KLEV

    REAL(KIND=JPRB), INTENT(IN)    :: PT(KLON,KLEV)
    REAL(KIND=JPRB), INTENT(IN)    :: PRHO(KLON,KLEV)
    REAL(KIND=JPRB), INTENT(IN)    :: PI(KLON,KLEV)
    REAL(KIND=JPRB), INTENT(IN)    :: PA(KLON,KLEV)    
    REAL(KIND=JPRB), INTENT(IN)    :: PAP(KLON,KLEV)
    REAL(KIND=JPRB), INTENT(IN)    :: PQSAT(KLON,KLEV)
    REAL(KIND=JPRB), INTENT(IN)    :: PSO4MASS(KLON,KLEV)
    REAL(KIND=JPRB), INTENT(IN)    :: PBCMASS(KLON,KLEV)
    REAL(KIND=JPRB), INTENT(IN)    :: PDUMASS(KLON,KLEV)
    
    !   Input and output:
    REAL(KIND=JPRB), INTENT(INOUT) :: PGFL(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NDIM)
    REAL(KIND=JPRB), INTENT(INOUT) :: PICNC(KLON,KLEV)
    REAL(KIND=JPRB), INTENT(INOUT) :: PRE_ICE(KLON,KLEV)
    
    
    !---local data
    INTEGER(KIND=JPIM) :: JL,JK
    
    REAL(KIND=JPRB)    :: ZICENUCLEI 
    REAL(KIND=JPRB)    :: ZICNC
    REAL(KIND=JPRB)    :: ZRE_ICE
    REAL(KIND=JPRB)    :: ZQS        ! saturation
    REAL(KIND=JPRB)    :: ZS0
    REAL(KIND=JPRB)    :: ZTEMPC
    REAL(KIND=JPRB)    :: ZSO4MASS
    REAL(KIND=JPRB)    :: ZBCMASS
    REAL(KIND=JPRB)    :: ZDUMASS
    REAL(KIND=JPRB)    :: ZMAERMEAN_SO4
    REAL(KIND=JPRB)    :: ZMAERMEAN_BC
    REAL(KIND=JPRB)    :: ZMAERMEAN_DU
    REAL(KIND=JPRB)    :: ZNICEHOMO
    REAL(KIND=JPRB)    :: ZSO
    REAL(KIND=JPRB)    :: ZSVP
    REAL(KIND=JPRB)    :: ZSCRITHOMO
    REAL(KIND=JPRB)    :: ZWTOT
    REAL(KIND=JPRB)    :: ZRHO_ICE
    REAL(KIND=JPRB)    :: ZRICE_CRIT
    REAL(KIND=JPRB)    :: ZNCRIT_GIERENS, ZNCRIT_REN
    REAL(KIND=JPRB)    :: ZCLD
    REAL(KIND=JPRB)    :: ZZEPSEC
    REAL(KIND=JPRB)    :: ZHOOK_HANDLE

    
    !---executable procedure
    IF (LHOOK) CALL DR_HOOK('ICE_CLD_PROP', 0, ZHOOK_HANDLE)
    ASSOCIATE(YGFL=>YDMODEL%YRML_GCONF%YGFL,YDECLDP=>YDMODEL%YRML_PHY_EC%YRECLDP)
    ASSOCIATE(YICNC=>YGFL%YICNC, YRE_ICE=>YGFL%YRE_ICE, &
     & RNICE=>YDECLDP%RNICE, RCLDMAX=>YDECLDP%RCLDMAX)

    ! move to cldp module
    ZWTOT=0.1_JPRB ! governs critical N

    ! Table of column/annual mean mass of aerosol
    ! converted to microgram per m**3
    ! ZMAERMEAN_SO4=1.02E-09_JPRB*1.E9_JPRB 
    ZMAERMEAN_SO4=0._JPRB
    ! ZMAERMEAN_SS=2.12E-10_JPRB*1.E9_JPRB 
    ZMAERMEAN_DU=1.01E-09_JPRB*1.E9_JPRB 
    ZMAERMEAN_BC=3.05E-11_JPRB*1.E9_JPRB 

    ZZEPSEC=1.E-10_JPRB

    ZRHO_ICE=900.0_JPRB

    !---------------------------------------------------------------------
    ! Turn aerosol mass into a Ice Number concentration for ice processes
    !---------------------------------------------------------------------
    DO JK=1,KLEV
       DO JL=KIDIA,KFDIA

          !---aerosol masses in ug m-3
          ! ZSO4MASS = PSO4MASS(JL,JK)*PRHO(JL,JK)*1.E9_JPRB 
          ZSO4MASS = 0._JPRB
          ZBCMASS = PBCMASS(JL,JK)*PRHO(JL,JK)*1.E9_JPRB
          ZDUMASS = PDUMASS(JL,JK)*PRHO(JL,JK)*1.E9_JPRB

          !                       0.01_JPRB is "default" value from
          ! Demott et al. Ice SS=55% or Meyers et al. 1992 JAS, ISS=25%
          ! In a prognostic scheme this will be function of clear sky humidity

          ! By relating IN to Aerosol mass we are assuming that the mode of the 
          ! Aerosol size distribution lies in the accumulation or coarse mode

          ! The relationship will implicitly introduce the exponential height
          ! dependence that Sassen (1992) and K and Curry (1998) explicitly 
          ! introduced to their parameterizations.

          ZICENUCLEI=0.01_JPRB*&
               &  (ZSO4MASS + ZBCMASS + ZDUMASS) &
               & /(ZMAERMEAN_SO4 + ZMAERMEAN_BC +ZMAERMEAN_DU)
          ZICENUCLEI=MAX(ZICENUCLEI,0.0_JPRB)

          ! T in oC
          ZTEMPC=PT(JL,JK)-RTT

          ! Re form for ice crystals from Liou and Oort 1994
          ! used to derive Re(ice) as in Lohmann JC 2002 
          ZNICEHOMO=0.0_JPRB
          ZRE_ICE=0.5_JPRB*(326.3_JPRB+ZTEMPC* &
               & (12.42_JPRB + ZTEMPC*(0.197_JPRB + ZTEMPC*0.0012_JPRB)))
          ZRE_ICE=MAX(ZRE_ICE,0.0_JPRB)

          ! effect Re to volume mean from S Moss or Lohmann and Kaercher papers 200?
          ! ZRE_ICE on both LHS and RHS ??
          ZRE_ICE=(MAX(SQRT(5.113E6_JPRB+2.809E3_JPRB*ZRE_ICE**3.0_JPRB)-2.261E3_JPRB,0.0_JPRB))**0.333_JPRB
          ZRE_ICE=MAX(ZRE_ICE,1.0_JPRB)  ! diameter minimum 1.0 microns

          ! more default values if not applying
          ZICNC=RNICE ! place as default

          ZQS=PQSAT(JL,JK)
          ZQS=MIN(0.5_JPRB,ZQS)
          ZQS=ZQS/(1.0_JPRB-RETV*ZQS)

          IF (PT(JL,JK)<238._JPRB .AND. PI(JL,JK)>ZZEPSEC) THEN
             ZS0=1.3_JPRB
             ZSCRITHOMO=2.349_JPRB-PT(JL,JK)/259.0_JPRB !ren form of Koop 2000 
             ZSVP=MAX(ZZEPSEC,ZQS*PAP(JL,JK)/0.622_JPRB)

             ! Klaus Gierens critical ice nuclei: Gierens (2003)
             ZNCRIT_GIERENS=2.81E11_JPRB*(10.0_JPRB**(4.0_JPRB-0.02_JPRB*PT(JL,JK)))**0.75_JPRB&
                  &*(ZWTOT**1.5_JPRB)*PAP(JL,JK)**1.5_JPRB/&
                  &(PT(JL,JK)**5.415_JPRB*(1.5_JPRB*ZSVP)**0.5_JPRB*(ZSCRITHOMO-ZS0)**0.75_JPRB)
             ZNCRIT_GIERENS=ZNCRIT_GIERENS/1.E6_JPRB ! cm**-3
             
             ! Ren and Mackensie QJRMS 2005 critical ice nuclei
             ZNCRIT_REN=5.4E10_JPRB*(ZWTOT**1.5_JPRB)*PAP(JL,JK)**1.5_JPRB*&
                  & (ZSCRITHOMO/(ZSCRITHOMO-1.0_JPRB))**1.5_JPRB/ &
                  & (PT(JL,JK)**5.415_JPRB*(1.5_JPRB*ZSVP)**0.5_JPRB)
             ZNCRIT_REN=ZNCRIT_REN/1.E6_JPRB ! cm**-3
             
             ! from Re derive the number concentration - here ice density is 900 kg/m**3 
             ! Re is in microns, 1e18 factor
             ZCLD=PI(JL,JK)/MAX(PA(JL,JK),ZZEPSEC)
             ZCLD=MIN(MAX(ZCLD,0.0_JPRB),RCLDMAX)
             IF (ZCLD>ZZEPSEC) THEN
                ZNICEHOMO=0.75_JPRB*PRHO(JL,JK)*ZCLD/(RPI*ZRHO_ICE*1.0E-18_JPRB*ZRE_ICE**3.0_JPRB)
             ENDIF
             ZNICEHOMO = ZNICEHOMO/1.E6_JPRB ! cm**-3
             
             ! following Ren and Mackensie, 2005, linearly interpolate to get Ice number
             IF (ZICENUCLEI<ZNCRIT_REN) THEN
                ZICNC =  ZICENUCLEI+(1.0_JPRB-ZICENUCLEI/ZNCRIT_REN)*ZNICEHOMO
             ELSE
                ZICNC = ZICENUCLEI 
             ENDIF
             
             !---why is this recalculated here ? 
             ZRE_ICE=(0.75_JPRB*PRHO(JL,JK)*ZCLD/(RPI*ZRHO_ICE*1.E6_JPRB*ZICNC))**0.333_JPRB
             ZRE_ICE=ZRE_ICE*1.E6_JPRB

             !PGFL(JL,JK,YICNC%MP9_PH) = ZICNC
             PICNC(JL,JK) = ZICNC
             PRE_ICE(JL,JK) = 1.E-6_JPRB*ZRE_ICE
             !PGFL(JL,JK,YRE_ICE%MP9_PH) = ZRE_ICE
             !PGFL(JL,JK,YRE_ICE%MP9_PH) = 1.E-6_JPRB*ZRE_ICE !eehol add to PGFL in meters (convert from um to m)
          ENDIF
       ENDDO
    ENDDO

    END ASSOCIATE
    END ASSOCIATE

    IF (LHOOK) CALL DR_HOOK('ICE_CLD_PROP', 1, ZHOOK_HANDLE)
  END SUBROUTINE ICE_CLOUD_PROP

!  SUBROUTINE PDF_UPDRAFT(KIDIA, KFDIA, KLON, KTDIA, KLEV, KPDF, PRHO, PLSM, PVERVEL, PW, PWPDF) 
!
!    !---inherited functions, types, variables and constants
!    USE YOMHOOK, ONLY: LHOOK, DR_HOOK
!    USE YOMCST,  ONLY: RG, RPI
!    !USE YOECLDP, ONLY: NACTPDF
!    USE MO_ACTIV, ONLY: nw
!
!    IMPLICIT NONE
!
!    !---subroutine interface
!    !   Input:
!    INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
!    INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
!    INTEGER(KIND=JPIM), INTENT(IN) :: KLON
!    INTEGER(KIND=JPIM), INTENT(IN) :: KTDIA
!    INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
!    INTEGER(KIND=JPIM), INTENT(IN) :: KPDF
!    REAL(KIND=JPRB),    INTENT(IN) :: PRHO(KLON,KLEV)
!    REAL(KIND=JPRB),    INTENT(IN) :: PLSM(KLON)
!    REAL(KIND=JPRB),    INTENT(IN) :: PVERVEL(KLON,KLEV)
!    !   Output:
!    REAL(KIND=JPRB),    INTENT(OUT) :: PW(KLON,KLEV,KPDF)
!    REAL(KIND=JPRB),    INTENT(OUT) :: PWPDF(KLON,KLEV,KPDF)
!
!    !---local data
!    REAL(KIND=JPRB), PARAMETER :: RW_MIN = 0._JPRB
!    !xxx 4 sigma is too much
!    REAL(KIND=JPRB), PARAMETER :: RW_MAX = 4._JPRB*PPDF_SIGMA
!    REAL(KIND=JPRB), PARAMETER :: RINVSIGMA = 1._JPRB/PPDF_SIGMA
!    REAL(KIND=JPRB), PARAMETER :: R2SIGMA2 = 2._JPRB*PPDF_SIGMA**2
!
!    REAL(KIND=JPRB), PARAMETER :: RUPDRAFT_LAND = 1.0_JPRB
!    REAL(KIND=JPRB), PARAMETER :: RUPDRAFT_SEA = 0.5_JPRB
!
!    INTEGER(KIND=JPIM) :: JL,JK, JW
!    
!    REAL(KIND=JPRB) :: ZWLARGE(KLON,KLEV)
!    REAL(KIND=JPRB) :: ZW_PRESC(KLON)
!    REAL(KIND=JPRB) :: ZSQ2PI 
!    REAL(KIND=JPRB) :: ZW_WIDTH
!    REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!    !---executable procedure
!    IF (LHOOK) CALL DR_HOOK('YOE_AER_ACTIV.PDF_UPDRAFT',0,ZHOOK_HANDLE)
!
!    ZSQ2PI = 1._JPRB/SQRT(2._JPRB*RPI)
!
!    ZWLARGE(KIDIA:KFDIA,KTDIA:KLEV) = -1._JPRB* PVERVEL(KIDIA:KFDIA,KTDIA:KLEV) / &
!                                   &  (RG*PRHO(KIDIA:KFDIA,KTDIA:KLEV))
!
!    IF (KPDF > 1) THEN 
!       
!       ZW_WIDTH = (RW_MAX - RW_MIN) / REAL(nw,JPRB)
!
!       DO JW=1,KPDF
!          DO JK=KTDIA,KLEV
!             PW(KIDIA:KFDIA,JK,JW) = RW_MIN + (REAL(JW,JPRB) - 0.5_JPRB) * ZW_WIDTH
!             
!             PWPDF(KIDIA:KFDIA,JK,JW) = ZSQ2PI * RINVSIGMA * &
!                                  & EXP( -(PW(KIDIA:KFDIA,JK,JW) - ZWLARGE(KIDIA:KFDIA,JK))**2._JPRB &
!                                  & / R2SIGMA2)
!          END DO
!       END DO
!    ELSE
!       ZW_PRESC(KIDIA:KFDIA) = MERGE(RUPDRAFT_LAND, RUPDRAFT_SEA,      &
!                                  &  PLSM(KIDIA:KFDIA) > 0.5_JPRB)
!       DO JK=KTDIA,KLEV
!          PW(KIDIA:KFDIA,JK,1) = MERGE(ZW_PRESC(KIDIA:KFDIA), 0._JPRB, &
!                                  &    ZWLARGE(KIDIA:KFDIA,JK) > RW_MIN)
!          PWPDF(KIDIA:KFDIA,JK,1) = 1._JPRB
!       END DO
!    END IF
!
!    IF (LHOOK) CALL DR_HOOK('YOE_AER_ACTIV.PDF_UPDRAFT',1,ZHOOK_HANDLE)
!  END SUBROUTINE PDF_UPDRAFT
    
!  SUBROUTINE GET_CDNC_FACTOR(KIDIA,KFDIA,KLON,PGEMU,PGELAM,PMAC2SP_CDNC_FACTOR)
!
!   USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
!   USE YOMCST   , ONLY : RPI, RDAY
!   USE YOMCT0   , ONLY : LTWOTL, LNF
!   USE YOMCT2   , ONLY : NSTAR2
!   USE YOMDYN   , ONLY : TSTEP
!   USE YOMRIP   , ONLY : NINDAT, NSTADD
!   USE YOMLUN   , ONLY : NULOUT
!   USE YOERAD   , ONLY : MAC2SP_LAMDA, NCMIPFIXYR, LMAC2SP
!   USE AER_MACV2SP_MOD, ONLY: SP_AOP_PROFILE, MAC2SP_YEAR_MIN, MAC2SP_YEAR_MAX
!   USE DAY_NUMBER_MOD,  ONLY: NUMBER_OF_DAY
!
!   IMPLICIT NONE
!
!   INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
!   INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
!   INTEGER(KIND=JPIM), INTENT(IN) :: KLON
!   REAL(KIND=JPRB),    INTENT(IN) :: PGEMU(KLON)
!   REAL(KIND=JPRB),    INTENT(IN) :: PGELAM(KLON)
!   REAL(KIND=JPRB),    INTENT(OUT):: PMAC2SP_CDNC_FACTOR(KLON)
!
!   !---the CDNC factor doesn't vary with height, thus only 1 level is needed
!   INTEGER(KIND=JPIM), PARAMETER :: KLEV=1
!
!   REAL(KIND=JPRB) :: ZR2D, ZGLAT(KLON), ZGLON(KLON)
!   REAL(KIND=JPRB) :: YEAR_FR, KNDY
!   REAL(KIND=JPRB) :: PGEOH(KLON,KLEV+1)
!
!   INTEGER(KIND=JPIM) :: IDY0, IMN0, IYR0, IDY, IMN, IYR, IDOY
!   INTEGER(KIND=JPIM) :: IZT, ISTADD, ITIME
!   INTEGER(KIND=JPIM) :: ILMONTH(12)
!
!   INTEGER(KIND=JPIM) :: JL
!
!   !-- variables for MAC2SP anthropogenic aerosol
!   REAL(KIND=JPRB)::  AOD_MAC2SP(KLON,KLEV)
!   REAL(KIND=JPRB)::  SSA_MAC2SP(KLON,KLEV)
!   REAL(KIND=JPRB)::  ASY_MAC2SP(KLON,KLEV)
!
!   REAL(KIND=JPRB)    :: ZHOOK_HANDLE
!
!#include "updcal.intfb.h"
!#include "fcttim.h"
!
!   !---executable procedure
!   IF (LHOOK) CALL DR_HOOK('GET_CDNC_FACTOR', 0, ZHOOK_HANDLE)
!
!   ! fake geopotential, doesn't matter for the CDNC factor
!   PGEOH(:,1)=100._JPRB
!   PGEOH(:,2)=10._JPRB
!
!   ZR2D=180.0_JPRB/RPI
!   DO JL=KIDIA,KFDIA
!     ZGLAT(JL)= ASIN(PGEMU(JL)) * ZR2D
!     ZGLON(JL)= PGELAM(JL) * ZR2D
!   END DO
!
!   ! Prepare year fraction
!   IF(.NOT.LNF.AND.NSTADD == 0) THEN
!     ! IN CASE OF RESTART:
!     ITIME=NINT(TSTEP,KIND(ITIME))
!     IF (LTWOTL) THEN
!        IZT=NINT(TSTEP*(REAL(NSTAR2,JPRB)+0.5_JPRB),KIND(IZT))
!     ELSE
!        IZT=INT(ITIME,KIND(IZT))*INT(NSTAR2,KIND(IZT))
!     ENDIF
!     ISTADD=INT(IZT/NINT(RDAY,KIND(IZT)),KIND(ISTADD))
!   ELSE
!     ISTADD=NSTADD
!   ENDIF
!   IYR0=NCCAA(NINDAT)
!   IMN0=NMM(NINDAT)
!   IDY0=NDD(NINDAT)
!
!   CALL UPDCAL(IDY0,IMN0,IYR0, ISTADD, IDY,IMN,IYR, ILMONTH, NULOUT)
!
!   ! Day number and total number of days for that year
!   CALL NUMBER_OF_DAY(IDY,IMN,IYR,IDOY)
!   IF(MOD(IYR,4) == 0 .AND. MOD(IYR,400) /= 100 &
!        & .AND. MOD(IYR,400) /= 200 .AND. MOD(IYR,400) /= 300)THEN
!      KNDY=366._JPRB
!   ELSE
!      KNDY=365._JPRB
!   ENDIF
!
!   ! Replace IYR with NCMIPFIXYR
!   IF (NCMIPFIXYR>0) IYR=NCMIPFIXYR
!
!   ! Limit IYR to available dataset
!   IYR =  MIN(MAC2SP_YEAR_MAX, MAX(MAC2SP_YEAR_MIN, IYR))
!
!   YEAR_FR = IYR + (REAL(IDOY,JPRB) - 0.5_JPRB)/KNDY
!
!   CALL SP_AOP_PROFILE (KLEV ,KIDIA ,KFDIA ,KLON ,MAC2SP_LAMDA(1) , &
!   &            ZGLON ,ZGLAT ,YEAR_FR ,PGEOH ,PMAC2SP_CDNC_FACTOR , &
!   &            AOD_MAC2SP ,SSA_MAC2SP ,ASY_MAC2SP )
!
!   IF (LHOOK) CALL DR_HOOK('GET_CDNC_FACTOR',1,ZHOOK_HANDLE)
! END SUBROUTINE GET_CDNC_FACTOR
! SUBROUTINE SETUP_ACI_DIAG
!
!   USE YOMHOOK, ONLY: LHOOK, DR_HOOK
!   USE PHY_DIAG_MOD, ONLY: NEW_DIAG_SET, NEW_DIAGNOSTIC 
!   USE YOM_GRIB_CODES, ONLY: NGRBCDNC, NGRBREFF, NGRBLIQCLDTIME
!   USE YOMLUN, ONLY: NULOUT
!
!   IMPLICIT NONE
!
!   INTEGER(KIND=JPIM) :: ISET
!   REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!   IF (LHOOK) CALL DR_HOOK('YOE_AER_ACTIV.SETUP_ACI_DIAG',0,ZHOOK_HANDLE)
!
!   CALL NEW_DIAG_SET('ACI', ISET)
!   CALL NEW_DIAGNOSTIC(ISET, 'CDNC', 3, NGRBCDNC, .TRUE., .FALSE.,  D_CDNC)
!   CALL NEW_DIAGNOSTIC(ISET, 'REFF', 3, NGRBREFF, .FALSE., .TRUE.,  D_REFF)
!   CALL NEW_DIAGNOSTIC(ISET, 'LIQ_CLD_TIME', 3, NGRBLIQCLDTIME, .TRUE., .FALSE.,  & 
!                     & D_LIQCLDT)
!
!   IF (LHOOK) CALL DR_HOOK('YOE_AER_ACTIV.SETUP_ACI_DIAG',1,ZHOOK_HANDLE)
! END SUBROUTINE SETUP_ACI_DIAG

END MODULE YOE_AER_ACTIV
