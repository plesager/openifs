  MODULE YOE_AER_ACTIV

  !---inherited functions, types, variables and constants 
  USE PARKIND1,            ONLY: JPIM, JPRB
  USE YOMHOOK,             ONLY: LHOOK, DR_HOOK, JPHOOK
  !USE PHY_DIAG_MOD,        ONLY: T_DIAG

  IMPLICIT NONE

  !---public member functions
  PUBLIC AER_ACTIV 

  !---private member functions
  PRIVATE AER_ACTIV_MORALES_NENES_FULL
  PRIVATE GET_HAMM7_AERO_PROP
  PRIVATE ICE_CLOUD_PROP

  !---module types, variables and constants

  ! surface tension of pure water at 273.15K [J m-2]
  !REAL(KIND=JPRB), PARAMETER, PUBLIC :: PPSURFTEN = 0.075_JPRB   

  ! standard deviation of updraft PDF
  !REAL(KIND=JPRB), PARAMETER, PUBLIC :: PPDF_SIGMA = 0.8_JPRB

  ! Diagnostics
  !TYPE(T_DIAG), POINTER, PUBLIC :: D_CDNC
  !TYPE(T_DIAG), POINTER, PUBLIC :: D_REFF
  !TYPE(T_DIAG), POINTER, PUBLIC :: D_LIQCLDT
 
CONTAINS
  SUBROUTINE AER_ACTIV(KIDIA,   KFDIA,  KTDIA,   KLON,    KLEV,   KSTGLO, &
                     !&  KLEVX,   KFLDX,  KFLDX2,                   &
                     &  PAPH,    PAP,    PT,      PQ,      PQSAT,  &
                     &  PVERVEL, PA,     PL,      PI,              &
                     &  PLSM,    PGELAM,   PGEMU, & !PSLON,   PGEMU,  &
                     &  PGFL, YDMODEL, PCDNCACT, PICNC, PREFFL, PREFFI, PSMAX, PDRYRSOLU, PXTM1, KTRAC, PSIGMA_W, &
                     &  PFRACN, PPMINCDNC, PPDEFCDNC, PQLWC, LLIQCLD, LICECLD, PPREFFL_DEF, PPREFFI_DEF)!,    PEXTRA, PEXTR2)
   
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
   USE YOMCST,      ONLY: RD, RPI, RG
   USE TM5M7_DATA,  ONLY: NSOL
   USE MO_HAM,      ONLY: NCLASS
   USE TYPE_MODEL,  ONLY: MODEL

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
   REAL(KIND=JPRB), INTENT(IN)    :: PPMINCDNC              ! minimum CDNC [# cm-3]
   REAL(KIND=JPRB), INTENT(IN)    :: PPDEFCDNC              ! default (background) CDNC [# cm-3]
   REAL(KIND=JPRB), INTENT(IN)    :: PQLWC(KLON,KLEV)       ! LWC [kg kg-1]
   LOGICAL,         INTENT(IN)    :: LLIQCLD(KLON,KLEV)     ! logical for liquid cloud
   LOGICAL,         INTENT(IN)    :: LICECLD(KLON,KLEV)     ! logical for ice cloud
   REAL(KIND=JPRB), INTENT(IN)    :: PPREFFI_DEF            ! default (background) ice eff rad [mu m]
   REAL(KIND=JPRB), INTENT(IN)    :: PPREFFL_DEF            ! default (background) liq eff rad [mu m]

   REAL(KIND=JPRB), INTENT(INOUT) :: PGFL(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NDIM)
   REAL(KIND=JPRB), INTENT(INOUT) :: PCDNCACT(KLON,KLEV)      ! cloud droplet number concentration [# cm-3]
   REAL(KIND=JPRB), INTENT(INOUT) :: PICNC(KLON,KLEV)         ! ice crystal number concentration [# cm-3]
   REAL(KIND=JPRB), INTENT(INOUT) :: PREFFL(KLON,KLEV)        ! liquid droplet effective radius [um]
   REAL(KIND=JPRB), INTENT(INOUT) :: PREFFI(KLON,KLEV)        ! ice effective radius [um]
   REAL(KIND=JPRB), INTENT(INOUT) :: PSMAX(KLON,KLEV)         ! maximum supersaturation [%]
   REAL(KIND=JPRB), INTENT(OUT)   :: PFRACN(KLON,KLEV,NCLASS) ! fraction of activated particles per mode

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
   !REAL(KIND=JPRB) :: ZW(KLON,KLEV,nw)                      ! updraft speed [m/s]
   !REAL(KIND=JPRB) :: ZWPDF(KLON,KLEV,nw)                   ! updraft probability
   REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
   REAL(KIND=JPRB) :: ZTMP(KLON,KLEV)                       ! interim for diagnostics

   REAL(KIND=JPRB) :: ZCDNC(KLON,KLEV)                      ! cloud droplet number concentration [#/cm-3]
   REAL(KIND=JPRB) :: ZICNC(KLON,KLEV)                      ! ice crystal number concentration [#/cm-3]
   REAL(KIND=JPRB) :: ZRE_LIQ(KLON,KLEV)                    ! liquid droplet effective radius [um]
   REAL(KIND=JPRB) :: ZRE_ICE(KLON,KLEV)                    ! ice crystal effective radius [um]
   REAL(KIND=JPRB) :: ZSMAX(KLON,KLEV)                      ! maximum supersaturation [%]
   REAL(KIND=JPRB) :: ZNACT_TOT                             ! variables for modewise activated fraction calculations
   REAL(KIND=JPRB) :: ZNACT_AS(KLON,KLEV)
   REAL(KIND=JPRB) :: ZNACT_CS(KLON,KLEV)
   REAL(KIND=JPRB) :: ZNACT_KS(KLON,KLEV)                   ! variables for modewise activated fraction calculations
   REAL(KIND=JPRB) :: ZFRAC_KS,ZFRAC_AS,ZFRAC_CS            ! variables for modewise activated fraction calculations

   LOGICAL :: LLIQCLDD(KLON,KLEV)                           ! true if liquid cloud is present (for activation calculations)

   LOGICAL :: LL1
   LOGICAL :: LBULK, LMODE                                  ! fetch HAMM7 aerosols as bulk mass / per-mode mass

   LOGICAL :: LCALCINCLOUD = .TRUE.                         ! calculate activation only in-cloud (T) or everywhere (F)
   REAL(KIND=JPRB) :: ZEPS                                  ! epsilon(1.)
   INTEGER(KIND=JPIM) :: ITOP                               ! highest level for water cloud
   INTEGER(KIND=JPIM) :: JK, JL, JMOD                       ! loop indices
   !INTEGER(KIND=JPIM) :: IX                                 ! index to PEXTRA diagnostic

   REAL(KIND=JPRB), PARAMETER :: ZEPSEC  = 1.E-14_JPRB      ! taken from cloudsc.F90
   REAL(KIND=JPRB) :: ZTMPA

   !---executable procedure
   IF (LHOOK) CALL DR_HOOK('YOE_AER_ACTIV.AER_ACTIV', 0, ZHOOK_HANDLE)
   ASSOCIATE(YGFL=>YDMODEL%YRML_GCONF%YGFL,YDECLDP=>YDMODEL%YRML_PHY_EC%YRECLDP)
   ASSOCIATE(YCDNC=>YGFL%YCDNC, YICNC=>YGFL%YICNC, YRE_LIQ=>YGFL%YRE_LIQ, YRE_ICE=>YGFL%YRE_ICE, &
      & LAERICESED=>YDECLDP%LAERICESED, LAERICEAUTO=>YDECLDP%LAERICEAUTO, &
      & RLMIN=>YDECLDP%RLMIN, RAMIN=>YDECLDP%RAMIN, RTHOMO=>YDECLDP%RTHOMO, RNICE=>YDECLDP%RNICE)

   ! Init
   PFRACN(KIDIA:KFDIA,:,:) = 0._JPRB

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
   ZRE_ICE(KIDIA:KFDIA,1:KLEV) = PPREFFI_DEF

   !---get aerosols from HAMM7
   LMODE = .TRUE.
   LBULK = (LAERICESED .OR. LAERICEAUTO) ! .OR. NCLOUDACT==JP_ACT_MENON)
   CALL GET_HAMM7_AERO_PROP(KIDIA, KFDIA, KLON, KTDIA, KLEV, KSTGLO, LMODE, LBULK, &
                          & PAPH,        PGELAM,  PGEMU, PXTM1, KTRAC, &
                          & ZDRYRSOLOLD, ZAERONUM, &
                          & ZSO4MASS,    ZBCMASS, ZOMMASS, ZSSMASS, ZDUMASS, &
                          & ZSO4BULK,    ZBCBULK, ZOMBULK, ZSSBULK, ZDUBULK, &
                          & ZNO3MASS,    ZMSAMASS)

   ZDRYRSOL(KIDIA:KFDIA,1:KLEV,:) = PDRYRSOLU(KIDIA:KFDIA,1:KLEV,:)

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

     ! The new scheme either approximates the integral over the updraft velocity PDF
     ! by Gauss-Legendre quadrature or uses a PPDEFsingle characteristic velocity.
     ! PDF_UPDRAFT, NACTPDF, ZW and ZWPDF are not used.
      CALL AER_ACTIV_MORALES_NENES_FULL(KIDIA, KFDIA, ITOP, KLON, KLEV, LLIQCLDD, PT, PAP, ZRHO, &
                                      & PVERVEL, ZSO4MASS, ZBCMASS, ZOMMASS, ZSSMASS, &
                                      & ZDUMASS, ZNO3MASS, ZMSAMASS, ZAERONUM, ZDRYRSOL, &
                                      & ZCDNC, ZSMAX, PGEMU, PSIGMA_W) !PSLON, PGEMU)                                    
                                      !& PGFL(:,:,YCDNC%MP9_PH), ZSMAX, PGEMU, PSIGMA_W) !PSLON, PGEMU)
                                      !& PGFL(:,:,YCDNC%MP9_PH), KFLDX, PEXTRA, PSLON, PGEMU)

      !---limit CDNC to min PPMINCDNC, set default value for CDNC outside clouds
      DO JK=KTDIA,KLEV
         ZCDNC(KIDIA:KFDIA,JK)=MERGE( &
         & MAX(ZCDNC(KIDIA:KFDIA,JK),PPMINCDNC), &
         & PPDEFCDNC, LLIQCLD(KIDIA:KFDIA,JK) )
      END DO

      !---cloud liquid water: droplet effective radius is computed in radlswr now
      !CALL LIQ_CLOUD_RE(KIDIA, KFDIA, KLON, ITOP, KLEV, LLIQCLD, PL, PA, ZRHO, PGFL)

      ! liquid effective radius                                                                                                                         
      DO JK=1,KLEV
         DO JL=KIDIA,KFDIA                                                                                                                                        
            ! effective radius calculated similarly as in radlswr.F90                                                                                     
            ! 2.387e-10 is 3/(4*pi*rho_liq*10^6)  [10^6 for N in right units]                                                                             
            ZRE_LIQ(JL,JK) = 1.E+06_JPRB*(2.387e-10_JPRB*ZRHO(JL,JK)*PQLWC(JL,JK)/ZCDNC(JL,JK))**0.333_JPRB
         END DO
      END DO
      ZRE_LIQ(KIDIA:KFDIA,1:KLEV) = MERGE(ZRE_LIQ(KIDIA:KFDIA,1:KLEV),PPREFFL_DEF,LLIQCLD(KIDIA:KFDIA,1:KLEV))
      
   END IF

   !---cloud ice: ICNC and effective radius for ice crystals 
   IF (LAERICESED .OR. LAERICEAUTO) THEN
      CALL ICE_CLOUD_PROP(KIDIA, KFDIA, KLON, KLEV, PT, ZRHO, PI, PA, PAP, &
                       &  PQSAT, ZSO4BULK, ZBCBULK, ZDUBULK, PGFL, YDMODEL, ZRE_ICE, ZICNC)        
   END IF


   !eehol: calculate modewise fraction of activated particles
   !assume only KS, AS, CS modes to be activated
   DO JK=1,KLEV
     DO JL=KIDIA,KFDIA
       ZNACT_TOT = 1.0E6_JPRB * ZCDNC(JL,JK) / ZRHO(JL,JK) !calculate total number of activated particles in #/kg
       IF ( ZAERONUM(JL,JK,4) > 1.E-9_JPRB ) THEN
         ZFRAC_CS = MAX(ZNACT_TOT,0._JPRB) / ZAERONUM(JL,JK,4)
       ELSE
         ZFRAC_CS = 0._JPRB
       ENDIF
       PFRACN(JL,JK,4) = MAX(0._JPRB,MIN(ZFRAC_CS,1._JPRB)) !threshold between 0 and 1
       ZNACT_CS(JL,JK) = ZAERONUM(JL,JK,4) * PFRACN(JL,JK,4) !calculate activated number for CS mode

       IF ( ZAERONUM(JL,JK,3) > 1.E-9_JPRB ) THEN
         ZFRAC_AS = MAX(ZNACT_TOT - ZNACT_CS(JL,JK),0._JPRB) / ZAERONUM(JL,JK,3)
       ELSE
         ZFRAC_AS = 0._JPRB
       ENDIF
       PFRACN(JL,JK,3) = MAX(0._JPRB,MIN(ZFRAC_AS, 1._JPRB)) !threshold between 0 and 1
       ZNACT_AS(JL,JK) = ZAERONUM(JL,JK,3) * PFRACN(JL,JK,3) !calculate activated number for AS mode

       IF ( ZAERONUM(JL,JK,2) > 1.E-9_JPRB ) THEN
         ZFRAC_KS = MAX((ZNACT_TOT - ZNACT_CS(JL,JK) - ZNACT_AS(JL,JK)),0._JPRB) / ZAERONUM(JL,JK,2)
       ELSE
         ZFRAC_KS = 0._JPRB
       ENDIF
       PFRACN(JL,JK,2) = MAX(0._JPRB,MIN(ZFRAC_KS, 1._JPRB)) !threshold between 0 and 1
       ZNACT_KS(JL,JK) = ZAERONUM(JL,JK,2) * PFRACN(JL,JK,2) !calculate activated number for KS mode

       PFRACN(JL,JK,2) = MERGE(PFRACN(JL,JK,2),0._JPRB,LLIQCLD(JL,JK)) !fraction only where there is liquid cloud
       PFRACN(JL,JK,3) = MERGE(PFRACN(JL,JK,3),0._JPRB,LLIQCLD(JL,JK)) !fraction only where there is liquid cloud
       PFRACN(JL,JK,4) = MERGE(PFRACN(JL,JK,4),0._JPRB,LLIQCLD(JL,JK)) !fraction only where there is liquid cloud
     END DO
   END DO
         
            
   !eehol: diagnostics:
   !--CDNC
   PCDNCACT(KIDIA:KFDIA,1:KLEV) = 1.0E6_JPRB*ZCDNC(KIDIA:KFDIA,1:KLEV) ! CDNC [#/m3]

   !--ICNC
   PICNC(KIDIA:KFDIA,1:KLEV) = ZICNC(KIDIA:KFDIA,1:KLEV) ! ICNC [#/cm3]

   !--Liq eff rad
   PREFFL(KIDIA:KFDIA,1:KLEV) = MERGE(ZRE_LIQ(KIDIA:KFDIA,1:KLEV),PPREFFL_DEF,LLIQCLD(KIDIA:KFDIA,1:KLEV)) !eehol: output liq eff rad [um]

   !--Ice eff rad (only if there is ice cloud else minimum value)
   PREFFI(KIDIA:KFDIA,1:KLEV) = MERGE(ZRE_ICE(KIDIA:KFDIA,1:KLEV),PPREFFI_DEF,LICECLD(KIDIA:KFDIA,1:KLEV))
   
   !--Maximum supersaturation
   PSMAX(KIDIA:KFDIA,1:KLEV) = ZSMAX(KIDIA:KFDIA,1:KLEV) !eehol: output maximum supersaturation [%]
   

   END ASSOCIATE
   END ASSOCIATE

   IF (LHOOK) CALL DR_HOOK('YOE_AER_ACTIV.AER_ACTIV', 1, ZHOOK_HANDLE)

  END SUBROUTINE AER_ACTIV



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
    USE ND_PARAM, ONLY: CCNSPEC, PDFACTIV, NDPARAM 

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
    REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE

    ! Variables to interface with Nenes routines
    TYPE(NDPARAM)   :: BOX
    REAL(KIND=JPRB) :: TPI(NSOL-1), DPGI(NSOL-1), SIGI(NSOL-1), AKKI(NSOL-1)
    REAL(KIND=JPRB) :: TPARC, PPARC, WPARC, SG(NSOL-1), NACT, SMAX

    REAL(KIND=JPRB), PARAMETER :: A = 2.25 ! Default FHH adsorption parameters (in the case of FHH-AT)
    REAL(KIND=JPRB), PARAMETER :: B = 1.20 ! See Kumar et al., (2011) ACP
    REAL(KIND=JPRB), PARAMETER :: ACCOM = 1.0 ! Accommodation coefficient

    ! Standard deviation of the updraft velocity distribution (m/s)
    ! For the moment it is set to a constant value of 0.8 m/s,
    ! as we did for the other activation schemes.
    ! In reality it depends on turbulence characteristics.
    ! Several parameterizations have been proposed 
    ! (see e.g., Hoose et al., 2010; Zheng et al., GRL, 2016)
    ! A common approach is to use TKE or, alternatively, 
    ! the vertical diffusion coefficient 
    ! (see module VDIFLCZ in sinvect directory).
    !REAL(KIND=JPRB), PARAMETER :: SIGW = 0.6_JPRB
    REAL(KIND=JPRB)    ::  SIGW ! sigma_w is input now

    ! Logical switch to use a single characteristic velocity
    ! instead of Gauss-Legendre quadrature.
    LOGICAL, PARAMETER :: CHAR_VELOCITY = .FALSE.

    INTEGER(KIND=JPIM) :: MODEI(NSOL-1)

    !--- executable procedure
    IF (LHOOK) CALL DR_HOOK('YOE_AER_ACTIV.AER_ACTIV_MORALES_NENES_FULL',0,ZHOOK_HANDLE)

    !--- 0) Initializations:

    PCDNC(KIDIA:KFDIA,KTDIA:KLEV) = 0._JPRB
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
                  ZKAPPA(JL,JK,JMOD) = MIN(ZKAPPA(JL,JK,JMOD),1.2_JPRB)
               ELSE
                  ZKAPPA(JL,JK,JMOD) = 0.04_JPRB  ! if total volume per mode is too small, use minimum kappa
               END IF
            END IF
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
               MODEI(JMOD-1) = 1                                ! Kohler mode
               TPI(JMOD-1) = PAERONUM(JL,JK,JMOD) * PRHO(JL,JK) ! Number concentration (#/m3)
               DPGI(JMOD-1) = 2._JPRB * PRDRY(JL,JK,JMOD)       ! Modal diameter (m)
               SIGI(JMOD-1) = SIGMA(JMOD)                       ! Geometric dispersion (sigma_g)
               AKKI(JMOD-1) = ZKAPPA(JL,JK,JMOD)                ! Hygroscopicity parameter (kappa)
            END DO
            TPARC = PT(JL,JK) ! Temperature (K)
            PPARC = PAP(JL,JK) ! Pressure (Pa)
         
            IF ( ANY(TPI(:) .GE. ZEPS) .AND. ANY(DPGI(:) .GE. 1e-9_JPRB) .AND. TPARC.GE.(273.15_JPRB-35.0_JPRB) ) THEN !eehol: any num con, diam and temperature need to be over treshold

               ! Convert aerosol data into CCN, fill BOX object
               CALL CCNSPEC (TPI,DPGI,SIGI,MODEI,TPARC,PPARC,NSOL-1,AKKI,A,B,ACCOM,BOX) 

               ! xxx To be done:
               ! Save CCN spectra for supersaturations:
               ! S = 0.05, 0.1, 0.2, 0.3, 0.5, 1.0 %
               ! New routine needs as input TPI, SG, S and returns CCN(S)
               ! which needs to be put into the output as 6 3-D fields

               !eehol: give sigma_w a value depending on the input variable
               SIGW = MAX(0.1_JPRB,PSIGMA_W(JL,JK)) ! threshold sigma to min value

               ! Calculate activated droplet number
               IF (CHAR_VELOCITY) THEN
                 ! Use characteristic updraft velocity (m/s)
                 ! from Morales and Nenes (2010)
                 WPARC = 0.79*SIGW
                 ! Call activation for a single velocity
                 ! equal to the characteristic velocity
                 CALL PDFACTIV (WPARC,0._JPRB,NACT,SMAX,BOX) 
               ELSE
                 ! Call activation for velocity PDF (SIGW is non-zero)
                 ! with WPARC set equal to the large-scale velocity (m/s)
                 WPARC = ZWLARGE(JL,JK)
                 CALL PDFACTIV (WPARC,SIGW,NACT,SMAX,BOX)
               ENDIF

               ! convert CDNC to # cm-3
               PCDNC(JL,JK) = 1.E-6_JPRB * NACT

               ! convert Smax to %
               PSMAX(JL,JK) = 100._JPRB * SMAX 
            END IF
          END IF ! LCLOUD
       END DO !jl
    END DO !jk

    IF (LHOOK) CALL DR_HOOK('YOE_AER_ACTIV.AER_ACTIV_MORALES_NENES_FULL',1,ZHOOK_HANDLE)

  END SUBROUTINE AER_ACTIV_MORALES_NENES_FULL 



  SUBROUTINE GET_HAMM7_AERO_PROP(KIDIA, KFDIA, KLON, KTDIA, KLEV, KSTGLO, LMODE, LBULK, &
                               & PAPH,     PGELAM,  PGEMU, PXTM1, KTRAC,                &
                               & PDRYRSOL, PAERONUM, &
                               & PSO4MASS, PBCMASS, POMMASS, PSSMASS, PDUMASS,   &
                               & PSO4BULK, PBCBULK, POMBULK, PSSBULK, PDUBULK,   &
                               & PNO3MASS, PMSAMASS)

    !---inherited functions, types, variables and constants
    USE YOMCST,              ONLY: RPI 
    USE TM5M7_DATA,          ONLY: NSOL

    !---aerosols variables and indices
    USE MO_HAM, ONLY:     &
         SIZECLASS,       & ! Aerosol classes in HAM
         AEROCOMP           ! Aerosol compounds by size class in HAM
    USE MO_HAM_M7CTL,     ONLY: INUCS,  IAITS,  IACCS,  ICOAS,   &
                                IAITI,  IACCI,  ICOAI,           &
                                ISO4NS, ISO4KS, ISO4AS, ISO4CS,  &
                                IBCKS,  IBCAS,  IBCCS,  IBCKI,   &
                                IOCKS,  IOCAS,  IOCCS,  IOCKI,   &
                                ISSAS,  ISSCS,                   &
                                IDUAS,  IDUCS,  IDUAI,  IDUCI

    !USE YOERAD,              ONLY: LCMIP6_PI_AEROSOLS, NRADFR
    !USE YOE_AERO_M7_DATA 
    !USE YOE_PI_AERO         
    !USE YOMCT3,              ONLY: NSTEP

    IMPLICIT NONE

    !---subroutine interface
    !   *GET_HAMM7_AERO_PROP* is called from AER_ACTIV before the activation calculations
    !   
    !   INPUT:
    INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA   ! beginning of horizontal block
    INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA   ! end of horizontal block
    INTEGER(KIND=JPIM), INTENT(IN) :: KLON    ! horizontal dimension
    INTEGER(KIND=JPIM), INTENT(IN) :: KTDIA   ! highest level with liquid cloud
    INTEGER(KIND=JPIM), INTENT(IN) :: KLEV    ! number of model vertical levels
    INTEGER(KIND=JPIM), INTENT(IN) :: KSTGLO  ! offset of horizontal block in coupling arrays
    INTEGER(KIND=JPIM), INTENT(IN) :: KTRAC   ! number of tracers

    LOGICAL, INTENT(IN) :: LMODE              ! Per-mode data requested
    LOGICAL, INTENT(IN) :: LBULK              ! Bulk aerosol masses requested
    
    REAL(KIND=JPRB), INTENT(IN)    :: PAPH(KLON,KLEV+1)      ! half-level pressure
    REAL(KIND=JPRB), INTENT(IN)    :: PGELAM(KLON)           ! longitude
    REAL(KIND=JPRB), INTENT(IN)    :: PGEMU(KLON)            ! sine of latitude
    REAL(KIND=JPRB), INTENT(IN)    :: PXTM1(KLON,KLEV,KTRAC) ! tracer mixing ratios

    !   OUTPUT:
    REAL(KIND=JPRB), INTENT(OUT)   :: PDRYRSOL(KLON,KLEV,NSOL) ! [M]
    REAL(KIND=JPRB), INTENT(OUT)   :: PAERONUM(KLON,KLEV,NSOL) ! [#/KG(AIR)]
    REAL(KIND=JPRB), INTENT(OUT)   :: PSO4MASS(KLON,KLEV,NSOL) ! [KG(SO4)/KG(AIR)]
    REAL(KIND=JPRB), INTENT(OUT)   :: PBCMASS(KLON,KLEV,NSOL)  ! [KG(BC)/KG(AIR)]
    REAL(KIND=JPRB), INTENT(OUT)   :: POMMASS(KLON,KLEV,NSOL)  ! [KG(OM)/KG(AIR)]
    REAL(KIND=JPRB), INTENT(OUT)   :: PSSMASS(KLON,KLEV,NSOL)  ! [KG(SS)/KG(AIR)]
    REAL(KIND=JPRB), INTENT(OUT)   :: PDUMASS(KLON,KLEV,NSOL)  ! [KG(DU)/KG(AIR)]
    REAL(KIND=JPRB), INTENT(OUT)   :: PSO4BULK(KLON,KLEV)      ! [KG(SO4)/KG(AIR)]
    REAL(KIND=JPRB), INTENT(OUT)   :: PBCBULK(KLON,KLEV)       ! [KG(BC)/KG(AIR)]
    REAL(KIND=JPRB), INTENT(OUT)   :: POMBULK(KLON,KLEV)       ! [KG(OM)/KG(AIR)]
    REAL(KIND=JPRB), INTENT(OUT)   :: PSSBULK(KLON,KLEV)       ! [KG(SS)/KG(AIR)]
    REAL(KIND=JPRB), INTENT(OUT)   :: PDUBULK(KLON,KLEV)       ! [KG(DU)/KG(AIR)]
    REAL(KIND=JPRB), INTENT(OUT)   :: PNO3MASS(KLON,KLEV)      ! [KG(NO3)/KG(AIR)]
    REAL(KIND=JPRB), INTENT(OUT)   :: PMSAMASS(KLON,KLEV)      ! [KG(MSA)/KG(AIR)]

    !---local data
    REAL(KIND=JPRB)    :: NSO4, NH2SO4, NNACL, NNA, NCL, NNA2SO4 ! Particle numbers [kmol/kg air]
    REAL(KIND=JPRB)    :: ZDRYVOL2, ZDRYVOL3, ZDRYVOL4      ! Volume per particle [m3/#]
    REAL(KIND=JPRB)    :: ZDUMMY(KLON,KLEV)       ! For fields we don't use in subroutine call
    REAL(KIND=JPRB)    :: Z4PIOVER3
    REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE
    INTEGER(KIND=JPIM) :: JL, JK, IBL, IL, IK

    !---executable procedure
    IF (LHOOK) CALL DR_HOOK('YOE_AER_ACTIV.GET_HAMM7_AERO_PROP', 0, ZHOOK_HANDLE)

    !---Aerosols masses and numbers
    ! NOTE: PXTM1 indices need to go according to HAM indices and not OIFS!
    IF (LMODE) THEN
      DO JK=KTDIA,KLEV
         DO JL=KIDIA,KFDIA
            
            !---nucleation mode is ignored

            !---Aitken soluble mode: SU, OM, BC
            PSO4MASS(JL,JK,2) = PXTM1(JL,JK,AEROCOMP(ISO4KS)%IDT) !SO4 Ait sol
            POMMASS(JL,JK,2)  = PXTM1(JL,JK,AEROCOMP(IOCKS)%IDT)  !OC Ait sol
            PBCMASS(JL,JK,2)  = PXTM1(JL,JK,AEROCOMP(IBCKS)%IDT)  !BC Ait sol
            
            !---accumulation soluble mode: SU, OM, BC, SS, DU
            PSO4MASS(JL,JK,3) = PXTM1(JL,JK,AEROCOMP(ISO4AS)%IDT) !SO4 acc sol
            POMMASS(JL,JK,3)  = PXTM1(JL,JK,AEROCOMP(IOCAS)%IDT)  !OC acc sol
            PBCMASS(JL,JK,3)  = PXTM1(JL,JK,AEROCOMP(IBCAS)%IDT)  !BC acc sol
            PSSMASS(JL,JK,3)  = PXTM1(JL,JK,AEROCOMP(ISSAS)%IDT)  !SS acc sol
            PDUMASS(JL,JK,3)  = PXTM1(JL,JK,AEROCOMP(IDUAS)%IDT)  !DU acc sol
            
            !---coarse soluble mode: SU, OM, BC, SS, DU
            PSO4MASS(JL,JK,4) = PXTM1(JL,JK,AEROCOMP(ISO4CS)%IDT) !SO4 coa sol
            POMMASS(JL,JK,4)  = PXTM1(JL,JK,AEROCOMP(IOCCS)%IDT)  !OC coa sol
            PBCMASS(JL,JK,4)  = PXTM1(JL,JK,AEROCOMP(IBCCS)%IDT)  !BC coa sol
            PSSMASS(JL,JK,4)  = PXTM1(JL,JK,AEROCOMP(ISSCS)%IDT)  !SS coa sol
            PDUMASS(JL,JK,4)  = PXTM1(JL,JK,AEROCOMP(IDUCS)%IDT)  !DU coa sol

            !---Nitrate and MSA - TODO Currently not used in HAM... need to be added later!
            PNO3MASS(JL,JK)   = 0._JPRB
            PMSAMASS(JL,JK)   = 0._JPRB

            PAERONUM(JL,JK,2) = PXTM1(JL,JK,SIZECLASS(IAITS)%IDT_NO) !Ait sol
            PAERONUM(JL,JK,3) = PXTM1(JL,JK,SIZECLASS(IACCS)%IDT_NO) !acc sol
            PAERONUM(JL,JK,4) = PXTM1(JL,JK,SIZECLASS(ICOAS)%IDT_NO) !coa sol

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




  SUBROUTINE ICE_CLOUD_PROP(KIDIA, KFDIA, KLON, KLEV, PT, PRHO, PI, PA, PAP, &
                         &  PQSAT, PSO4MASS, PBCMASS, PDUMASS, PGFL, YDMODEL, PRE_ICE, PICNC)

    ! This subroutine is mainly a copy/paste of the ice phase microphysics 
    ! implemented in the subroutine aer_cld.F90 
 
    !---inherited functions, types, variables and constants
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
    REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE

    
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


END MODULE YOE_AER_ACTIV
