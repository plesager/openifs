SUBROUTINE HAMM7_INTERFACE( &
 & YDMODEL,   KIDIA,     KFDIA,     KLON,        KTDIA,    KLEV,     KTILES,  &
 & KFLDX,     KLEVX,     KTRAC,     KAERO,       KCHEM,    KSTGLO,   PGEOH,   &
 & PRS1,      PRSF1,     PAEROP,    PCAERO,      PCEN,     PAPHIF,            &
 & PFPLCL,    PFPLCN,    PFPLSL,    PFPLSN,      PGELAT,   PGELAM,            &
 & PAP,       PIP,       PLP,       PRP,         PSP,      PCOVPTOT,          &
 & PLU,       PO3P,      PQP,       PTP,         PTHP,     PTENC,    PCFLX,   &
 & PAERDDP,   PAERSDM,   PAERSRC,   PAERWS,      PAERGUST, PAERUST,  PAERMAP, &
 & PCLAERS,   PPRAERS,   PCHEM2AER,                                           &
 & PFRTI,     PLSM,      PSNS,      PWND,        PWS1,     PAERFLX,  PAERLIF, &
 & PAERODDF,  PTSPHY,    PGFL,                                                &
 & PODTO,     PAERO_WVL_DIAG,                                                 &
 & PAER_TAU,  PAER_SSA,  PAER_ASYM, PAER_TAU_LW,                              &
 & PTAUS_AER, PTAUA_AER, PPMAER,                                              &
 & PEXTRA,    PVERVEL,   PCCNL,     PCCNO,       PAHFSTI,  PCI,      PZ0M,    &
 !eehol: added here vertical velocity, CCN over land, CCN over ocean
 & PAHFLEV,   PUP,       PVP,       PCVL,        PCVH,     PSO2DD,   PGEMU)
 !, PTSO2, PTSO4, PTSO4_AQ, PFSO2,PFSO4,PFSO4_AQ&
 !  u-wind, v-wind, low veg. cover, high veg. cover, sine of latitude

! ╭────────────────────────────────────────────────────────────────────────────╮
! │                                                      (updated 30-APR-2024) │
! │ Purpose :                                                                  │
! │ -------                                                                    │
! │  *hamm7_interface* -                                                       │
! │                                                                            │
! │                                                                            │
! │ Interface :                                                                │
! │ ---------                                                                  │
! │   *HAMM7_INTERFACE* is called from AER_PHY3_LAYER                          │
! │                                                                            │
! │                                                                            │
! │ Input :                                                                    │
! │ -----                                                                      │
! │                                                                            │
! │                                                                            │
! │ Output :                                                                   │
! │ ------                                                                     │
! │                                                                            │
! │                                                                            │
! │ Externals :                                                                │
! │ ---------                                                                  │
! │                                                                            │
! │ Method :                                                                   │
! │ ------                                                                     │
! │                                                                            │
! │ Reference :                                                                │
! │ ---------                                                                  │
! │                                                                            │
! │ Author :                                                                   │
! │ -------                                                                    │
! │     Orginal version:                                                       │
! │     Vicent Huijen (KNMI), Tommi Bergman (FMI), Thomas Kuehn (FMI/UEF)      │
! │                                                                            │
! │ Modifications :                                                            │
! │ -------------                                                              │
! │     May.  2020 - V. Huijnen     : Modifications for TM5M7                  │
! │     Sep.  2020 - T. Bergman     : TM5M7 work                               │
! │     Apr.  2024 - Lianghai Wu    : revision for CY48r1                      │
! │     May.  2024 - R. Checa-Garcia: revision for CY48r1 and refactory        │
! │                                                                            │
! ╰────────────────────────────────────────────────────────────────────────────╯

! RCHG: (TODO)
!      This subroutine assumes YAEROUT() has a number of elements, we may have 
!      to introduce something that test that the number is consistent with what 
!      is needed. Also the namelist with YAEROUT specification may need to be  
!      consistent with how we fill things here (not sure about it). 
!
! RCHG: (TODO) 
!      For this subroutine the description of the steps is not very uniform 
!      For so long subroutines we may need a table of contents to any reader 
!      can understand what is doing or better per section headers or split 
!      with CONTAINS into different parts. The advantage of the last approach
!      is to indentify and isolate what are the input/outputs of each step.
!      The problem are the associate environments.


!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLEV*         NUMBER OF LEVELS
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KLEVS*        NUMBER OF SOIL LAYERS
!    *KTILES*       NUMBER OF TILES (I.E. SUBGRID AREAS WITH DIFFERENT 
!                   OF SURFACE BOUNDARY CONDITION)
!    *KVTYPES*      NUMBER OF biomes for land carbon
!    *KTRAC*        Number of tracers
!     KLEVXG  : number of levels to compute
!     KAERO                       : Number of aerosol fields 
!     KCHEM :  Number of chemistry tracers 
!     KFLDX   - number of the last field. of what??? TB
!     PRS1(KLON,0:KLEV)           : HALF-LEVEL PRESSURE           (Pa)
!     PRSF1(KLON,KLEV)             : FULL-LEVEL PRESSURE           (Pa)
!     PAEROP(KLON,KLEV,KAERO)     : Aerosol concentrations  (kg/kg) - Note that fields are only non-zero if NACTAERO > 0
!     PCAERO is the aerosol density                kg cm-3 (KWHAT=0)
!     PCAERO is the optical thickness              ND     (KWHAT=2)
!     PCEN(KLON,KLEV,KCHEM)       :  CONCENTRATION OF TRACERS           (kg/kg) TB:Chemistry
!     PAPHIF         : geopotential height "gz" at full levels.
!     PFPLCL (KLON,KLEV+1)        : CONVECTIVE PRECIPITATION AS RAIN  (kg/m2s)
!     PFPLCN   : convective precipitation as snow.
!    *PFPLSL*       LARGE SCALE RAIN FLUX                        KG/(M2*S)
!    *PFPLSN*       LARGE SCALE SNOW FLUX                        KG/(M2*S)
!     PGELAT(KLON)                : LATITUDE (RADIANS) 
!    *PGELAM*       LONGITUDE                                     RADIANS
!     PAP    : (KLON,KLEV)       ; CLOUD FRACTION
!     PIP     (KLON,KLEV)         :  ICWC                         (kg/kg)
!     PLP     (KLON,KLEV)         :  LWC                         (kg/kg)
!    *PSP*    (KLON,KLEV)         :  Snow water content           (kg/kg) 
!     PRP     (KLON,KLEV)         :  Rain water content           (kg/kg)
!     PCOVPTOT(KLON,KLEV)         :  PRECIP FRACTION               0..1  
!    *PLU*          LIQUID WATER CONTENT IN UPDRAFTS            KG/KG
!     PQP(KLON,KLEV)      : FULL-LEVEL HUMIDITY (W. DYN.TEND.) (kg kg-1)
!     PTP     (KLON,KLEV)         : TEMPERATURE                   (K)!
!     pthp ?  ZPPTHPW(KPROMA,KOPLEV)      ! theta'w???
!     PTENC  (KLON,KLEV,KTRAC)     : TENDENCY OF CONCENTRATION OF TRACERS including chemistry(kg/kg s-1)
!     PCFLX(KLON,KTRAC)      : SURFACE FLUX OF TRACERS              (xx m-2)
!     PAERDDP(KLON,NACTAERO) : aerosol dry deposition 
!     PAERSDM(KLON,NACTAERO)  : aerosol sedimentation 
!     PAERSRC(KLON,NACTAERO) : SOURCE FLUX                          (xx m-2) 
!     PAERWS                 : Wind speed (average of horizontal wind speed)   (m/s)
!     PAERGUST               : Wind gust (maximum 3 second gust in the hour)   (m/s)
!     PAERMAP(KLON,5)        : DUST MASK-RELATED QUANTITIES
! ??? PCLAERS(KLON)    ! aerosol
!     PPRAERS(KLON) ! radvis.f90
!     PALBD  : (KPROMA,NTSW)        ; DIFFUSE ALBEDO IN THE NSW SW INTERVALS
!    *PFRTI*      TILE FRACTIONS                                   (0-1)
!     PLSM   (KLON) : land-sea mask                       [0-1]
!     PSNS       : MASS OF SNOW PER UNIT SURFACE
!     PWND(KLON)                  : Surface wind
!     PWS1   : REAL     : TOP LAYER SOIL MOISTURE CONTENT
!     PAERFLX(KLON,12,9)     : DIAGNOSTIC DUST SOURCE FLUXES
!     PAERLIF(KLON,9)        : DIAGNOSTIC LIFTING THRESHOLD SPEED
!     PAERODDF    ! Diagnostic array with aerosol fluxes (sources and sinks)
!     PTSPHY                : TIMESTEP                  (s)
!     PGFL         : GFL fields
!     PVERVEL      : Vertical velocity [Pa s-1]
!     PCCNL        : CCN over land
!     PCCNO        : CCN over ocean
!    *PAHFSTI*      SURFACE SENSIBLE HEAT FLUX                    W/M2
!     PCI(KLON)           : FRACTION OF SEA-ICE
!     PZ0M                : roughness length for momentum              (m)
!     PAHFLEV             : latent heat flux                           (W/m2)
!     PUP                 : u-wind (m/s)
!     PVP                 : v-wind (m/s)
!     PCVL                : low vegetation cover
!     PCVH                : high vegetation cover
!     PGEMU               : sine of latitude
!     PGEOH               : geopotential at half levels (M2/S2)
!------------

!-----------------------------------------------------------------------

USE PARKIND1,     ONLY: JPIM, JPRB
USE YOMHOOK,      ONLY: LHOOK, DR_HOOK, JPHOOK
USE TYPE_MODEL,   ONLY: MODEL
USE YOMCST,       ONLY: RD, RG, RPI
USE YOMCT0,       ONLY: LIFSMIN, LIFSTRAJ
USE YOMCT3,       ONLY: NSTEP
USE TM5M7_DATA,   ONLY: MODAL_DATA, MODE_TRACERS, MODE_START,      &
                      & MODE_TRACERS_BY_MODS, MODE_END_SO4, NAERMOD, NMOD, NSOL,&
                      & IISVOC, IELVOC, IACS_N, ISO4, ISO4ACS,ISO4COS
USE YOMCHEM,      ONLY: IEXTR_WD, IEXTR_CH, IEXTR_NG, IEXTR_DD, IEXTR_CHTR

USE YOE_AERODIAG, ONLY: JPAERO_WVL_AOD, JPAERO_WVL_AODABS, JPAERO_WVL_AODFM,   &
                      & JPAERO_WVL_SSA, JPAERO_WVL_ASSIMETRY
USE YOMLUN,       ONLY: NULOUT

! [RCHG -> var non used ]  USE YOMCST,       ONLY: RMSO2, RMSO4, RMD, RNAVO
! [RCHG -> var non used ]  USE YOESRTCOP,    ONLY: RSASWA, RSASWB, RSFUA0, RSFUA1
! [RCHG -> var non used ]  USE YOMCOMPO , ONLY :  YRCOMPO  

! implementation of HAM-M7
USE MO_HAM,                  ONLY: nclass, naerocomp, sizeclass, nccndiag, subm_ngasspec
USE OIFS_TO_HAM,             ONLY: ind_oifs_ham
USE MO_HAM_SUBM,             ONLY: HAM_SUBM_INTERFACE !eehol: replaced HAM-M7 call with submodel interface
USE MO_ACTIV,                ONLY: activ_updraft,nw, idt_cdnc, idt_icnc !eehol: HAM-M7 activation updraft calculation, effective radii
USE MO_HAM_ACTIV,            ONLY: ham_activ_abdulrazzak_ghan, ham_activ_koehler_ab !eehol: HAM-M7 activation
USE MO_PARAM_SWITCHES,       ONLY: ncd_activ !eehol: for activation
USE MO_TRACDEF,              ONLY: ntrac, trlist
                            !eehol: number of tracer for mass/number mixing ratio conversion, trlist for wet deposition flags
USE MO_TRACER_PROCESSES,     ONLY: xt_borrow !eehol: conserving the negative tracer values from tendency
USE MO_SUBMODEL,             ONLY: lwetdep,lsedimentation !eehol: logical for wetdeposition, sedimentation
USE MO_TIME_CONTROL,         ONLY: time_step_len ! time step length for tendency
USE MO_HAMMOZ_WETDEP,        ONLY: wetdep_interface ! wet deposition interface call
USE MO_HAM_WETDEP,           ONLY: ham_conv_lfraq_so2
USE MO_HAMMOZ_SEDIMENTATION, ONLY: sedi_interface ! sedimentation interface call
USE MO_HAMMOZ_DRYDEP,        ONLY: drydep_interface ! dry deposition interface call
USE MO_HAM_RAD,              ONLY: ham_rad,ham_rad_cache_cleanup,ham_rad_cache

! [RCHG -> non used] USE MO_SPECIES,              ONLY: speclist !SO2 wetdep for simple sulfur scheme
! [RCHG -> non used] USE mo_ham_species,          ONLY: id_so2 !SO2 wetdep for simple sulfur scheme
! [RCHG -> non used] USE YOMMP0,                  ONLY : MYPROC, NPROC
 
USE TM5M7_OPTICS_DATA,       ONLY : NWDEP, NASWBAND, ASWBAND !,WDEP, AER_TAU, AER_SSA,AER_ASYM,AER_TAU_LW
USE TM5_PHOTOLYSIS,          ONLY : NBANDS_TROP, WAV_GRID, WAV_GRIDA
USE TM5M7_EMIS_DATA,         ONLY : VKARMAN ! von karman constant for dry deposition
!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(MODEL),        INTENT(IN) :: YDMODEL
INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA, KFDIA, KLON
INTEGER(KIND=JPIM), INTENT(IN) :: KTDIA, KLEV, KFLDX, KLEVX
INTEGER(KIND=JPIM), INTENT(IN) :: KTILES
INTEGER(KIND=JPIM), INTENT(IN) :: KTRAC
INTEGER(KIND=JPIM), INTENT(IN) :: KAERO(YDMODEL%YRML_GCONF%YGFL%NAERO)
INTEGER(KIND=JPIM), INTENT(IN) :: KCHEM(YDMODEL%YRML_GCONF%YGFL%NCHEM)
INTEGER(KIND=JPIM), INTENT(IN) :: KSTGLO

REAL(KIND=JPRB),INTENT(IN)    :: PGEOH(KLON,0:KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PRSF1(KLON,KLEV), PRS1(KLON,0:KLEV), PAPHIF(KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PAP(KLON,KLEV), PIP(KLON,KLEV), PLP(KLON,KLEV), PLU(KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PRP(KLON,KLEV), PSP(KLON,KLEV), PCOVPTOT(KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PAEROP(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NACTAERO) 
REAL(KIND=JPRB),INTENT(IN)    :: PCAERO(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NACTAERO)
REAL(KIND=JPRB),INTENT(IN)    :: PCEN(KLON,KLEV,KTRAC), PCFLX(KLON,KTRAC)
REAL(KIND=JPRB),INTENT(IN)    :: PO3P(KLON,KLEV), PQP(KLON,KLEV), PTP(KLON,KLEV), PTHP(KLON,0:KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PFPLCL(KLON,0:KLEV),PFPLCN(KLON,0:KLEV),PFPLSL(KLON,0:KLEV),PFPLSN(KLON,0:KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PGELAT(KLON)    , PGELAM(KLON)
REAL(KIND=JPRB),INTENT(IN)    :: PFRTI(KLON,KTILES)
REAL(KIND=JPRB),INTENT(IN)    :: PAERWS(KLON), PAERGUST(KLON), PAERUST(KLON), PAERMAP(KLON,5)
REAL(KIND=JPRB),INTENT(IN)    :: PAERSRC(KLON,YDMODEL%YRML_GCONF%YGFL%NACTAERO)
REAL(KIND=JPRB),INTENT(IN)    :: PCHEM2AER(KLON,KLEV,6)
REAL(KIND=JPRB),INTENT(IN)    :: PAERFLX(KLON,12,9), PAERLIF(KLON,9), PCLAERS(KLON)
REAL(KIND=JPRB),INTENT(IN)    :: PLSM(KLON)  , PSNS(KLON)    , PWND(KLON)   , PWS1(KLON)
REAL(KIND=JPRB),INTENT(IN)    :: PTSPHY
REAL(KIND=JPRB),INTENT(IN)    :: PVERVEL(KLON,KLEV) ! added vertical velocity as an input to TM5M7 (needed in HAM-M7 activation)
REAL(KIND=JPRB),INTENT(IN)    :: PCCNL(KLON) ! added CCN over land as an input to TM5M7 (needed in liquid effective radius calc.)
REAL(KIND=JPRB),INTENT(IN)    :: PCCNO(KLON) ! added CCN over ocean as an input to TM5M7 (needed in liquid effective radius calc.)
REAL(KIND=JPRB),INTENT(IN)    :: PAHFSTI(KLON,KTILES) ! added surface sensible heat flux for dry deposition
REAL(KIND=JPRB),INTENT(IN)    :: PCI(KLON) ! added fraction of sea-ice for dry deposition
REAL(KIND=JPRB),INTENT(IN)    :: PZ0M(KLON) ! added roughness length for momentum for dry deposition
REAL(KIND=JPRB),INTENT(IN)    :: PAHFLEV(KLON) ! added latent heat flux for dry deposition
REAL(KIND=JPRB),INTENT(IN)    :: PUP(KLON,KLEV) ! added u component of wind
REAL(KIND=JPRB),INTENT(IN)    :: PVP(KLON,KLEV) ! added v component of wind
REAL(KIND=JPRB),INTENT(IN)    :: PCVL(KLON) ! added low vegetation cover
REAL(KIND=JPRB),INTENT(IN)    :: PCVH(KLON) ! added high vegetation cover
REAL(KIND=JPRB),INTENT(IN)    :: PGEMU(KLON) ! sine of latitude

REAL(KIND=JPRB),INTENT(INOUT) :: PTENC(KLON,KLEV,KTRAC)
REAL(KIND=JPRB),INTENT(INOUT) :: PAERDDP(KLON,YDMODEL%YRML_GCONF%YGFL%NACTAERO)
REAL(KIND=JPRB),INTENT(INOUT) :: PAERSDM(KLON,YDMODEL%YRML_GCONF%YGFL%NACTAERO)
! Total optical depth at various wavelenghts. 
! NOTE!! These wavelength definitions are not necessarily consistent
! with what is used in IFS-AER
REAL(KIND=JPRB),INTENT(OUT)   :: PODTO(KLON)

!REAL(KIND=JPRB),INTENT(OUT)   :: PODTO469(KLON), PODTO670(KLON), PODTO865(KLON), PODTO1240(KLON)
!REAL(KIND=JPRB),INTENT(IN)    :: PALBD(KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NTSW), PFRTI(KLON,KTILES)

REAL(KIND=JPRB),INTENT(OUT)   :: PAERO_WVL_DIAG(KLON,YDMODEL%YRML_GCONF%YGFL%NAERO_WVL_DIAG,YDMODEL%YRML_GCONF%YGFL%NAERO_WVL_DIAG_TYPES)
REAL(KIND=JPRB),INTENT(OUT)   :: PAERODDF(KLON,YDMODEL%YRML_GCONF%YGFL%NACTAERO,8)
REAL(KIND=JPRB),INTENT(INOUT) :: PEXTRA(KLON,KLEVX,KFLDX)
REAL(KIND=JPRB),INTENT(INOUT) :: PAER_TAU(KLON,KLEV,14), PAER_SSA(KLON,KLEV,14),PAER_ASYM(KLON,KLEV,14)
REAL(KIND=JPRB),INTENT(INOUT) :: PAER_TAU_LW(KLON,KLEV,16)
REAL(KIND=JPRB),INTENT(OUT)   :: PTAUS_AER(KLON,KLEV,NBANDS_TROP,2),PTAUA_AER(KLON,KLEV,NBANDS_TROP,2)
REAL(KIND=JPRB),INTENT(OUT)   :: PPMAER(KLON,KLEV,NBANDS_TROP,2)
REAL(KIND=JPRB),INTENT(INOUT) :: PGFL(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NDIM), PPRAERS(KLON)
! Simple sulfur scheme variables:
REAL(KIND=JPRB),INTENT(INOUT)   :: PSO2DD(KLON)

!REAL(KIND=JPRB), INTENT(INOUT) :: PFSO2(KLON)  , PFSO4(KLON), PFSO4_AQ(KLON)
!REAL(KIND=JPRB), INTENT(INOUT) :: PTSO2(KLON, KLEV)  , PTSO4(KLON, KLEV), PTSO4_AQ(KLON, KLEV)

!REAL(KIND=JPRB), INTENT(INOUT) :: ZFSO2(KLON), ZFSO4(KLON), ZFSO4_AQ(KLON)
!REAL(KIND=JPRB), INTENT(INOUT) :: ZTSO2(KLON, KLEV), ZTSO4(KLON, KLEV), ZTSO4_AQ(KLON, KLEV)

!*   0.5    LOCAL VARIABLES
!           ---------------

INTEGER(KIND=JPIM) :: JAER, JK, JL, JWAVL, JT, JB, JN
INTEGER(KIND=JPIM) :: JEXT, ITRC, IKLEVTROP(KLON), IW
INTEGER(KIND=JPIM) :: JO, JH, JY                         ! inside loop index for OIFS contex, HAM context and YAEROUT 
INTEGER(KIND=JPIM) :: JCLASS, JTILE, JMASS, JGAS, JCLOUD ! local loop indice for activation and dry deposition and tracer indexing
INTEGER(KIND=JPIM) :: ISSO2, ISSO4, ISSO4_ACS
INTEGER(KIND=JPIM) :: IMODE 
INTEGER(KIND=JPIM) :: IFLAG

REAL(KIND=JPRB) :: ZAEROK(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NACTAERO)
REAL(KIND=JPRB) :: ZTAEROK(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NACTAERO)
REAL(KIND=JPRB) :: ZTAERO0(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NACTAERO)
REAL(KIND=JPRB) :: ZFAERO(KLON,ntrac)!YGFL%NACTAERO)
REAL(KIND=JPRB) :: ZAER(KLON,KLEV), ZAERNEG(KLON,KLEV)
REAL(KIND=JPRB) :: ZAP(KLON,KLEV) 
REAL(KIND=JPRB) :: ZSO2(KLON,KLEV), ZDP(KLON,KLEV), ZDZ(KLON,KLEV) 
REAL(KIND=JPRB) :: ZITSO2(KLON,KLEV)
REAL(KIND=JPRB) :: ZFSO2(KLON)  , ZFSO4(KLON), ZFSO4_AQ(KLON)
REAL(KIND=JPRB) :: ZTSO2(KLON, KLEV)  , ZTSO4(KLON, KLEV,1), ZTSO4_AQ(KLON, KLEV)
REAL(KIND=JPRB) :: ZQSAT(KLON,KLEV), ZRHO(KLON,KLEV)
REAL(KIND=JPRB) :: ZTAER(KLON,KLEV)
REAL(KIND=JPRB) :: ZTAERO(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NACTAERO)
REAL(KIND=JPRB) :: ZRH(KLON,KLEV),ZTENC0(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NCHEM)!,ZTSO4(KLON,KLEV)
!REAL(KIND=JPRB) :: ZM6RP(KLON,KLEV,NMOD)
REAL(KIND=JPRB) :: ZM6RP(KLON,KLEV,NMOD)! NMOD=7
REAL(KIND=JPRB) :: ZM6DRY(KLON,KLEV,NSOL)
REAL(KIND=JPRB) :: ZWW(KLON,KLEV,NMOD)
REAL(KIND=JPRB) :: ZRHOP(KLON,KLEV,NMOD)
REAL(KIND=JPRB) :: ZSVOC(KLON,KLEV)
REAL(KIND=JPRB) :: ZELVOC(KLON,KLEV)
REAL(KIND=JPRB) :: ZSO4G(KLON,KLEV)
REAL(KIND=JPRB) :: ZCEN(KLON,KLEV,KTRAC) ! local tracer number and mixing ratios and gas concentrations for not tendency updated values
REAL(KIND=JPRB) :: PODTO469(KLON), PODTO670(KLON), PODTO865(KLON), PODTO1240(KLON)
REAL(KIND=JPRB) :: ZAER_TAU(KLON,KLEV,14,1), ZAER_SSA(KLON,KLEV,14),ZAER_ASYM(KLON,KLEV,14),ZAER_TAU_LW(KLON,KLEV,16)

! Optics output fields (to be used and allocated by methods using the optics)
REAL(KIND=JPRB), DIMENSION(:,:,:),   ALLOCATABLE :: ZTAUS_AER, ZTAUA_AER, ZPMAER ! extinctions
REAL(KIND=JPRB), DIMENSION(:,:,:,:), ALLOCATABLE :: ZAOP_OUT_EXT ! extinctions
REAL(KIND=JPRB), DIMENSION(:,:,:),   ALLOCATABLE :: ZAOP_OUT_A   ! single scattering albedo
REAL(KIND=JPRB), DIMENSION(:,:,:),   ALLOCATABLE :: ZAOP_OUT_G   ! assymetry factor

!Defined here, but should be passed in Call really:
TYPE(MODAL_DATA), DIMENSION(NMOD), TARGET :: RW_MODE
TYPE(MODAL_DATA), DIMENSION(NSOL), TARGET :: RWD_MODE
TYPE(MODAL_DATA), DIMENSION(NSOL), TARGET :: H2O_MODE
TYPE(MODAL_DATA), DIMENSION(NMOD), TARGET :: DENS_MODE

REAL(KIND=JPRB), ALLOCATABLE ::    ZAERNGT(:,:)

REAL(KIND=JPRB) :: ZDEGRAD, ZEPSCOV, ZEPSWAT, ZRWSAT, ZRWPWP
REAL(KIND=JPRB) :: ZQLWP(KLON,KLEV)
REAL(KIND=JPRB) :: ZTMPA

LOGICAL         :: LLIQCLD(KLON,KLEV) ! logical for liquid cloud
LOGICAL         :: LICECLD(KLON,KLEV) ! logical for ice cloud

REAL(KIND=JPRB), PARAMETER :: ZEPSEC=1e-14

! [RCHG -> var non used ] INTEGER(KIND=JPIM) :: j_yaerom, JMMD, JSCAV, JSW, JSPEC
! [RCHG -> var. non used ] INTEGER(KIND=JPIM) :: IAER, IEX3D, IEX3DP
! [RCHG -> vas. non used ] INTEGER(KIND=JPIM) :: IEXTR2,ISHIFT1, IKPAER, IKP, ISTO, IWHERE
! [RCHG -> var. non used ] INTEGER(KIND=JPIM) :: NSO4SCHEME
! [RCHG -> non used ] REAL(KIND=JPRB) :: ZLAT, ZLON
! [RCHG -> non used ] REAL(KIND=JPRB) :: BETAB(KLON,KLEV), ZBETAI(KLON,KLEV)
! [RCHG -> non used ] REAL(KIND=JPRB) :: ZCLDWAT(KLON,KLEV), ZDUM(KLON,KLEV)
! [RCHG -> non used ] REAL(KIND=JPRB) :: ZBCPHI(KLON,KLEV), ZBCPHO(KLON,KLEV) 
! [RCHG -> non used ] REAL(KIND=JPRB) :: ZOMPHI(KLON,KLEV) , ZOMPHO(KLON,KLEV)
! [RCHG -> non used ] REAL(KIND=JPRB) :: ZTBCPHI(KLON,KLEV),ZTBCPHO(KLON,KLEV)
! [RCHG -> non used ] REAL(KIND=JPRB) :: ZTOMPHI(KLON,KLEV), ZTOMPHO(KLON,KLEV)
! [RCHG -> non used ] REAL(KIND=JPRB) :: ZITBCPHO(KLON,KLEV),ZITOMPHO(KLON,KLEV)
! [RCHG -> non used ] REAL(KIND=JPRB) :: ZAIRDM(KLON), ZRHCL(KLON,KLEV)   
! [RCHG -> non used ] REAL(KIND=JPRB) :: ZTAERI(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NACTAERO)
! [RCHG -> non used ] REAL(KIND=JPRB) :: ZAERWET(KLON,YDMODEL%YRML_GCONF%YGFL%NACTAERO)
! REAL(KIND=JPRB) :: ZAERNL(KLON,KLEV,NMOD)
! REAL(KIND=JPRB) :: ZAERML(KLON,KLEV,NAERMOD)
! [RCHG -> var non used ] INTEGER(KIND=JPIM) :: JMOD
! [RCHG -> var non used ] INTEGER(KIND=JPIM) :: JAERCLASS
! [RCHG -> non used ] REAL(KIND=JPRB) :: ZTENV(KLON)
! [RCHG -> non used ] REAL(KIND=JPRB) :: ZGDT
! [RCHG -> non used ] REAL(KIND=JPRB) :: ZVISICL, ZVISIPR, ZVISCAE, ZVISPAE 
!REAL(KIND=JPRB) :: PAER_TAU(KLON,KLEV,14), PAER_SSA(KLON,KLEV,14),PAER_ASYM(KLON,KLEV,14)
!REAL(KIND=JPRB), ALLOCATABLE ::    ZAERSRC(:,:),  ZAERNGT(:,:) , ZAERSCC(:,:)  
! [RCHG -> non used ] REAL(KIND=JPRB) :: ZQRWP, ZQSWP, ZQIWP, ZRANGE, ZRELRA, ZSIGAIR, ZSNOICE
! [RCHG -> non used ] REAL(KIND=JPRB) :: Z1CLD, ZCFLIRA, ZCFSNIC, ZDENSVIS, ZDESIC, ZCLWAT, ZLIQRAI, ZNS
! [RCHG -> non used ] REAL(KIND=JPRB) :: ZVISCON, ZVISRAY
! [RCHG -> non used ] REAL(KIND=JPRB) :: pmrateps(KLON,KLEV),pmrater(KLON,KLEV),pfevapr(KLON,KLEV)
! [RCHG -> non used ] REAL(KIND=JPRB) :: pfsubls(KLON,KLEV),pmsnowacl(KLON,KLEV)
! [RCHG -> non used ] REAL(KIND=JPRB) :: ZDPG(KLON,KLEV),zxtp1c(KLON,KLEV,KTRAC),zxtp10(KLON,KLEV,KTRAC)
! [RCHG -> non used ]    INTEGER(KIND=JPIM) ::KTOP
! [RCHG -> non used ]     REAL(KIND=JPRB) :: ZRHOP1(KLON,KLEV)!,zdummy,zdum2d(KLON,KLEV),zdum3d(KLON,KLEV,nmod)
! [RCHG -> non used ]     REAL(KIND=JPRB) :: zlfrac_so2(KLON,KLEV),ZFLXR,ZFLXS,ZFLXRB,ZFLXSB
! [RCHG -> non used ]     REAL(KIND=JPRB) :: ZAEROUT1(KLON,KLEV),ZAEROUT2(KLON,KLEV),ZAEROUT3(KLON,KLEV),ZAEROUT4(KLON,KLEV),ZAEROUT5(KLON,KLEV)

REAL(KIND=JPRB),PARAMETER :: INFINITY=HUGE(1._JPRB)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! variables for the M7 call
REAL(KIND=JPRB) :: ZGRVOL(KLON,KLEV) !grid box volume for diagnostics
REAL(KIND=JPRB) :: ZPBL(KLON) !planetary boundary layer top level (in vdfmain.F90 ITOP=1)
REAL(KIND=JPRB) :: ZXTM0(KLON,KLEV,ntrac) !tracer mixing ratios for HAM
REAL(KIND=JPRB) :: ZXTM1(KLON,KLEV,ntrac) !tracer mixing ratios for HAM
REAL(KIND=JPRB) :: ZXTTE(KLON,KLEV,ntrac) !tracer tendency for HAM
REAL(KIND=JPRB) :: ZXTTEM1(KLON,KLEV,ntrac) !tracer tendency for HAM
! added here variables for HAM-M7 activation
REAL(KIND=JPRB), ALLOCATABLE :: ZW(:,:,:) !mean or bins of updraft velocity [m s-1]
REAL(KIND=JPRB), ALLOCATABLE :: ZWPDF(:,:,:) !updraft velocity PDF over bins
REAL(KIND=JPRB), ALLOCATABLE :: ZRC(:,:,:,:) !critical radius of activation per mode [m]
REAL(KIND=JPRB), ALLOCATABLE :: ZSMAX(:,:,:) !maximum supersaturation
REAL(KIND=JPRB) :: ZTKEM1(KLON,KLEV) !turbulent kinetic energy  as zero as it is not used for now
REAL(KIND=JPRB) :: ZWCAPE(KLON) !CAPE  as zero as it is not used
REAL(KIND=JPRB) :: ZRDRY(KLON,KLEV,nclass) !dry radius for each class
REAL(KIND=JPRB) :: ZESW(KLON,KLEV) !saturation water vapor pressure
REAL(KIND=JPRB) :: ZA(KLON,KLEV,nclass) !curvature parameter A of the Koehler equation
REAL(KIND=JPRB) :: ZB(KLON,KLEV,nclass) !hygroscopicity parameter B of the Koehler equation
REAL(KIND=JPRB) :: ZSC(KLON,KLEV,nclass) !critical supersaturation [% 0-1]
REAL(KIND=JPRB) :: ZNACT(KLON,KLEV,nclass) !number of activated particles per mode [m-3]
REAL(KIND=JPRB) :: ZFRACN(KLON,KLEV,nclass) !fraction of activated particles per mode
REAL(KIND=JPRB) :: ZCDNCACT(KLON,KLEV)     !number of activated particles
REAL(KIND=JPRB) :: ZRE_LIQ(KLON,KLEV)! liquid effective radius
! variables for HAM-M7 wet deposition
REAL(KIND=JPRB) :: ZXTP1(KLON,KLEV,ntrac) !updated tracer mass/number mixing ratio
REAL(KIND=JPRB) :: ZXTP1C(KLON,KLEV,ntrac) !in-cloud tracer mass/number mixing ratio
REAL(KIND=JPRB) :: ZXTP10(KLON,KLEV,ntrac) !ambient tracer mass/number mixing ratio
REAL(KIND=JPRB) :: ZDUMMY(KLON,ntrac) !placeholder for pxtbound which is only necessary in the conv. case (conv. massfix boundary condition)
REAL(KIND=JPRB) :: ZDUM3D(KLON,KLEV,ntrac) !updraft mass flux (for conv case only)
REAL(KIND=JPRB) :: ZDUM2D(KLON,KLEV) !convective flux needed only for conv. case (see cuflx)
REAL(KIND=JPRB) :: ZLFRAC_SO2(KLON,KLEV) !liquid tracer fraction (SO2) -ham specific-
REAL(KIND=JPRB) :: ZDPG(KLON,KLEV) !dp/g
REAL(KIND=JPRB) :: ZQP(KLON,KLEV) !full level humidity with treshold
LOGICAL :: LSTRAT = .TRUE. !logical switch for stratiform or convective case (TRUE for strat., FALSE for conv.)
REAL(KIND=JPRB) :: ZFEVAPR(KLON,KLEV) !evaporation of rain [kg/m2/s]
REAL(KIND=JPRB) :: ZFSUBLS(KLON,KLEV) !sublimation of snow [kg/m2/s]
REAL(KIND=JPRB) :: ZMRATEPR(KLON,KLEV) !rain formation rate in cloudy part
REAL(KIND=JPRB) :: ZMRATEPS(KLON,KLEV) !ice formation rate in cloudy part
REAL(KIND=JPRB) :: ZMSNOWACL(KLON,KLEV) !accretion rate of snow with cloud droplets in cloudy part
REAL(KIND=JPRB) :: ZFLXR, ZFLXS, ZFLXRB, ZFLXSB !variables to calculate rain and snow evap/formation
REAL(KIND=JPRB) :: ZLP(KLON,KLEV) !temporary variable for cloud water content
REAL(KIND=JPRB) :: ZIP(KLON,KLEV) !temporary variable for cloud ice water content
! variables for ICNC calculations
REAL(KIND=JPRB) :: ZICNC(KLON,KLEV) ! ice crystal number concentration [#/cm3]
! added here variables for dry deposition and sedimentation
REAL(KIND=JPRB) :: ZTENCIH(KLON,KLEV,ntrac) !for HAM tendencies
REAL(KIND=JPRB) :: ZAHFSM(KLON)
REAL(KIND=JPRB) :: ZWND(KLON)
! variables not needed for HAM aerosol dry deposition (if gas deposition and different surfaces are taken into account these need to be revised!)
REAL(KIND=JPRB) :: ZCFML(KLON), ZCFMW(KLON), ZCFMI(KLON), ZCFNCL(KLON), ZCFNCW(KLON), ZCFNCI(KLON), ZEPDU2, ZKAP
REAL(KIND=JPRB) :: ZGEOM1(KLON,KLEV), ZRIL(KLON), ZRIW(KLON), ZRII(KLON), ZTVIR1(KLON,KLEV), ZTVL(KLON), ZTVW(KLON)
REAL(KIND=JPRB) :: ZTVI(KLON), ZAZ0(KLON), ZFRL(KLON), ZSRFL(KLON), ZFOREST(KLON), ZTSI(KLON), ZAZ0L(KLON)
REAL(KIND=JPRB) :: ZAZ0I(KLON), ZCDNI(KLON)
! variables needed for HAM aerosol dry deposition
LOGICAL         :: ZLOLAND(KLON) !logical land mask
REAL(KIND=JPRB) :: ZXTEMS(KLON,ntrac) !surface emissions modified by dry deposition
REAL(KIND=JPRB) :: ZAZ0W(KLON), ZFRW(KLON), ZCVS(KLON), ZCVW(KLON), ZVGRAT(KLON) !rough. len. wat., wat. frac., snow cov. frac., wet skin frac., veg. ratio
REAL(KIND=JPRB) :: ZCDNL(KLON), ZCDNW(KLON) !ustar (in not used variable), aerodynamic resis. on surface (in not used variable)
REAL(KIND=JPRB) :: ZXTMD1(KLON,KLEV,ntrac) !tracer mixing ratios for HAM drydep (updated with tend)
! output diagnostics
 REAL(KIND=JPRB) :: ZOUT(KLON,ntrac),ZOUT2(KLON,14),zout3(KLON,KLEV,2*(naerocomp+nclass))
!REAL(KIND=JPRB) :: ZOUT4(KLON,KLEV),ZOUT5(KLON,klev),ZOUT6(KLON,klev),ZOUT7(KLON,klev),ZOUT8(KLON,klev),ZOUT9(KLON,klev)
!REAL(KIND=JPRB) :: ZZOUT1(KLON,KLEV),ZZOUT2(KLON,klev),ZZOUT3(KLON,KLEV),ZZOUT4(KLON,KLEV),ZZOUT5(KLON,klev),ZZOUT6(KLON,KLEV)
!REAL(KIND=JPRB) :: ZZOUT7(KLON,KLEV),ZZOUT8(KLON,klev),ZZOUT9(KLON,KLEV),ZZOUT10(KLON,KLEV),ZZOUT11(KLON,klev),ZZOUT12(KLON,KLEV)
 
REAL(KIND=JPRB) :: sedout(KLON,KLEV,ktrac)   ! changed ntrack to ktrac (RCHG)
REAL(KIND=JPRB) :: ddepout(KLON,KLEV,KTRAC)
REAL(KIND=JPRB) :: Wdepout(KLON,KLEV,ktrac)
REAL(KIND=JPRB) :: sedout_2D(KLON,ktrac)
REAL(KIND=JPRB) :: ddepout_2D(KLON,KTRAC)
REAL(KIND=JPRB) :: Wdepout_2D(KLON,ktrac)
REAL(KIND=JPRB) :: M7TEND_OUT(KLON,KLEV,KTRAC)
REAL(KIND=JPRB) :: M7TEND_IN(KLON,KLEV,KTRAC)
REAL(KIND=JPRB) :: zaveragep(KLON,klev,(nclass+naerocomp))
REAL(KIND=JPRB) :: zm7kappa(KLON,klev,(nclass+naerocomp))
REAL(KIND=JPRB) :: zh2so4cs(KLON,klev,(nclass+naerocomp))
REAL(KIND=JPRB) :: zm7prodcond(KLON,klev,(nclass+naerocomp))
REAL(KIND=JPRB) :: ZVDA(KLON,YDMODEL%YRML_GCONF%YGFL%NACTAERO)
REAL(KIND=JPRB) :: ZSEDIFLUX(KLON,KLEV,ntrac) 
REAL(KIND=JPRB) :: ZDDEPFLUX(KLON,ntrac)
REAL(KIND=JPRB) :: ZDDEPFLUX_SO2(KLON)
REAL(KIND=JPRB) :: ZVDEP(KLON,ntrac) !ddep velocity for diagnostics from ham
INTEGER(kind=JPIM), parameter::ZKROW=1 ! KROW only used in ECHAM but needed inside HAM-codes so set as 1.
INTEGER(kind=JPIM) :: IBLK
REAL(KIND=JPRB)    :: reffi(KLON,KLEV,ZKROW), reffl(KLON,KLEV,ZKROW)
INTEGER(kind=JPIM) :: LWBANDS !laakso: number of LW bands
INTEGER(KIND=JPIM) :: INWAVL, ITWAVL(20)
REAL(KIND=JPRB)    :: PRS1D(KLON,KLEV)
INTEGER(kind=JPIM) :: ISO4_C, ISSO4_C ! temporary tracer index of gas-phase SO4 (retrieved from chemistry module)

! [RCHG -> non used ] REAL(KIND=JPRB) :: ZTENCI(KLON,KLEV,KTRAC) !for OIFS tendencies
! [RCHG -> non used ]  INTEGER(kind=JPIM)::ZISO4 ! temporary tracer index of gas-phase SO4 (retrieved from chemistry module)
! [RCHG -> non used ] REAL(KIND=JPRB):: zza(KLON,klev)
! [RCHG -> non used ]  REAL(KIND=JPRB)::zhenry_so2(2),zheneff(KLON,klev),ze3(KLON,klev), zqtp1(KLON,KLEV)
! [RCHG -> non used ] INTEGER(KIND=JPIM) :: IWHAT 
! [RCHG -> non used]  REAL(KIND=JPRB) :: TEST1(KLON,KLEV), TEST2(KLON,KLEV)
! [RCHG -> non used]  REAL(KIND=JPRB) :: TEST3(KLON,KLEV), TEST4(KLON,KLEV)
! [RCHG -> non used]  REAL(KIND=JPRB) :: TEST5(KLON,KLEV), TEST6(KLON,KLEV)
! [RCHG -> non used]  REAL(KIND=JPRB) :: TEST7(KLON,KLEV), TEST8(KLON,KLEV)

REAL(KIND=JPRB) :: PAOD(KLON,NASWBAND), PSSA(KLON,NASWBAND), PABS(KLON,NASWBAND), PASY(KLON,NASWBAND),PFAOD(KLON,NASWBAND)


!-----------------------------------------------------------------------

#include "surf_inq.h"
#include "aer_so2so4_v2.intfb.h"
#include "satur.intfb.h"
#include "aer_negat.intfb.h"
#include "tm5m7_optics_aop_get.intfb.h"
#include "troplev.intfb.h"
#include "chem_inext.intfb.h"
#include "m7_simple_sulfur_drydep.intfb.h"
#include "ice_effective_radius.intfb.h"
!#include "m7.intfb.h"
! include calculations for ice effective radius

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('HAMM7_INTERFACE',0,ZHOOK_HANDLE)
ASSOCIATE( &
         & YGFL         => YDMODEL%YRML_GCONF%YGFL,        & 
         & YDPHYRAD     => YDMODEL%YRML_PHY_RAD,           &
         & YDPHYAER     => YDMODEL%YRML_PHY_AER,           &
         & YREAERSRC    => YDMODEL%YRML_PHY_AER%YREAERSRC, &
         & YREAERATM    => YDMODEL%YRML_PHY_RAD%YREAERATM, &
         & YREAERLID    => YDMODEL%YRML_PHY_AER%YREAERLID, &
         & YREAERSNK    => YDMODEL%YRML_PHY_AER%YREAERSNK, &
         & YRERAD       => YDMODEL%YRML_PHY_RAD%YRERAD,    &
         & YDRIP        => YDMODEL%YRML_GCONF%YRRIP,       &
         & YDCHEM       => YDMODEL%YRML_CHEM%YRCHEM,       &
         & YRCOMPO      => YDMODEL%YRML_CHEM%YRCOMPO,      &
         & YREPHY       => YDMODEL%YRML_PHY_EC%YREPHY,     &
         & YRECLDP      => YDMODEL%YRML_PHY_EC%YRECLDP,    &
         & YDSPP_CONFIG => YDMODEL%YRML_GCONF%YRSPP_CONFIG)

 ASSOCIATE( &
         ! --- YGFL ------------------------------------------------
         & NACTAERO => YGFL%NACTAERO, NAERO     => YGFL%NAERO,     &
         & NGHG     => YGFL%NGHG,     NCHEM     => YGFL%NCHEM,     &
         & ZTRAC    => YGFL%NTRAC,    NDIM      => YGFL%NDIM,      &
         & YAEROUT  => YGFL%YAEROUT,  YR        => YGFL%YR,        &
         & YS       => YGFL%YS,       YCHEM     => YGFL%YCHEM,     &
         & YAERO    => YGFL%YAERO,    YGHG      => YGFL%YGHG,      &
         & YTRAC    => YGFL%YTRAC,    YAEROCLIM => YGFL%YAEROCLIM, &
         & YCDNC    => YGFL%YCDNC,    YICNC     => YGFL%YICNC,     &
         & YRE_LIQ  => YGFL%YRE_LIQ,  YRE_ICE   => YGFL%YRE_ICE,   &
         !included CDNC, ICNC, liq and ice eff rad
         & NAERO_WVL_DIAG => YGFL%NAERO_WVL_DIAG,                  &
         & LAERCHEM       => YGFL%LAERCHEM,                        &
         ! --- YRCOMPO ---------------------------------------------
         & LCHEM_DIA      => YRCOMPO%LCHEM_DIA,                    &
         & AERO_SCHEME    => YRCOMPO%AERO_SCHEME,                  & !CHEM_SCHEME=>YDCHEM%CHEM_SCHEME,&
         & LAERNITRATE    => YRCOMPO%LAERNITRATE,                  &
         & NSO4SCHEME     => YREAERSRC%NSO4SCHEME,                 &
         ! --- YREAERATM -------------------------------------------
         & LAERDRYDP      => YREAERATM%LAERDRYDP,                  &
         & LAERSEDIM      => YREAERATM%LAERSEDIM,                  &
         & LAERSURF       => YREAERATM%LAERSURF,                   & !add logicals for dry dep and sedi
         & LAER6SDIA      => YREAERATM%LAER6SDIA,                  &
         & LAERCLIMG      => YREAERATM%LAERCLIMG,                  &
         & LAERCLIMZ      => YREAERATM%LAERCLIMZ,                  &
         & LAERGTOP       => YREAERATM%LAERGTOP,                   &
         & LAERHYGRO      => YREAERATM%LAERHYGRO,                  &
         & LAERLISI       => YREAERATM%LAERLISI,                   &
         & LAERNGAT       => YREAERATM%LAERNGAT,                   &
         & LAERSCAV       => YREAERATM%LAERSCAV,                   &
         & LAERSCAV_CHEM  => YREAERATM%LAERSCAV_CHEM,              &
         & LAERVOL        => YREAERATM%LAERVOL,                    &
         & NXT3DAER       => YREAERATM%NXT3DAER,                   & 
         & LAERRRTM       => YREAERATM%LAERRRTM,                   &
         ! --- YRERAD ----------------------------------------------
         & LAERVISI       => YRERAD%LAERVISI,                      &
         & NTSW           => YRERAD%NTSW,                          &
         & RNS            => YRERAD%RNS,                           &
         & RSIGAIR        => YRERAD%RSIGAIR,                       & !!! YAERCLIM is now become YAEROCLIM
         & NRADFR         => YRERAD%NRADFR,                        & !FREQUENCY OF FULL RADIATION COMPUTATIONS
         & NAEROOPT       => YRERAD%NAEROOPT,                      &
         ! --- OTHERS ----------------------------------------------
         & YDAERM7        => YDPHYAER%YREAEROPT,                   & ! use this to transfer AOD, SSA and ASY to rad scheme
         & NWLID          => YREAERLID%NWLID,                      &
         & YSURF          => YREPHY%YSURF,                         &
         & RRHTAB         => YREAERSNK%RRHTAB,                     &
         & RNICE          => YRECLDP%RNICE,                        & !default for ICNC
         & RCLDMAX        => YRECLDP%RCLDMAX,                      & !max cloud value
         & NSTART         => YDRIP%NSTART                          ) 

! & NINDSCAV=>YREAERATM%NINDSCAV, NTSCAV=>YREAERATM%NTSCAV, &
! & NDDUST=>YREAERSRC%NDDUST, NTYPAER=>YREAERSRC%NTYPAER, &
!     ------------------------------------------------------------------

!*         0.     PROGNOSTIC AEROSOLS - FINAL COMPUTATIONS
!                 ----------------------------------------

!write(*,*) "PGEOH", minval(PGEOH)
!write(*,*) "PRS1",  minval(PRS1)   , minval(PRSF1)  , minval(PAEROP) 
!write(*,*) "PCAERO ",  minval(PCAERO) , minval(PCEN)  , minval(PAPHIF)
!write(*,*) "PFPLSN ",  minval(PFPLSN) , minval(PGELAT), minval(PGELAM)
!write(*,*) "PFPLCL ",  minval(PFPLCL) , minval(PFPLCN) , minval(PFPLSL)
!write(*,*) "PRP  ",  minval(PRP)  ,    minval(PSP)  
!write(*,*) "PAP",  minval(PAP)    , minval(PIP)    , minval(PLP)
!write(*,*) "1",  minval(PCOVPTOT), minval(PLU)    , minval(PO3P)  ,minval(PQP)   , minval(PTP)  
!write(*,*) "2",  minval(PTHP)   , minval(PTENC)  ,minval( PCFLX)
!write(*,*) "3",  minval(PAERDDP), minval(PAERSDM), minval(PAERSRC), minval(PAERWS)
!write(*,*) "4",   minval(PAERGUST), minval(PAERUST), minval(PAERMAP)
!write(*,*) "5",  minval(PCLAERS), minval(PPRAERS), minval(PCHEM2AER)
!write(*,*) "6.2",  minval(PFRTI)
!write(*,*) "6.3",  minval(PLSM)
!write(*,*) "6.4",  minval(PSNS) 
!write(*,*) "7",  minval(PWND)  , minval(PWS1)   , minval(PAERFLX), minval(PAERLIF)
!!write(*,*) "8",  minval(PAERODDF),PTSPHY 
!!write(*,*) "9",  minval(PODTO)  , minval(PAERO_WVL_DIAG)
!write(*,*) "10",  minval(PAER_TAU), minval(PAER_SSA),minval(PAER_ASYM),minval(PAER_TAU_LW)
!!write(*,*) "11",  minval(PTAUS_AER),minval(PTAUA_AER),minval(PPMAER)
!write(*,*) "12",  minval(PEXTRA), minval(PVERVEL), minval(PCCNL), minval(PCCNO)
!write(*,*) "13",  minval(PAHFSTI), minval(PCI), minval(PZ0M), minval(PAHFLEV)
!write(*,*) "14",  minval(PUP), minval(PVP), minval(PCVL), minval(PCVH), minval(PGEMU)!
!write(*,*) "15",  minval(PSO2DD)


ZAHFSM  = 0._JPRB
ZEPSCOV = 1.E-03_JPRB
ZEPSWAT = 1.E-18_JPRB
ZDEGRAD = 180._JPRB/RPI

CALL SURF_INQ(YSURF,PRWSAT=ZRWSAT)
CALL SURF_INQ(YSURF,PRWPWP=ZRWPWP)

!*         0.1    SWITCHING ON POSSIBLE DEBUG PRINTS AND ALLOCATING MEMORY
!                 --------------------------------------------------------

ALLOCATE( ZAERNGT(KLON,NACTAERO) )

! RCHG -> there is a recurrent error: initialize without KIDIA:KFDIA 
ZAERNGT(KIDIA:KFDIA,1:NACTAERO) = 0._JPRB

!ALLOCATE( ZAERSCC(KLON,NACTAERO) )
!ALLOCATE( ZAERSRC(KLON,NACTAERO) )
!ZAERSRC(KIDIA:KFDIA,1:NACTAERO)       =0._JPRB
!ZAERSCC(KIDIA:KFDIA,1:NACTAERO)       =0._JPRB
!ZAEROUT1(:,:) =0._JPRB
!ZAEROUT2(:,:) =0._JPRB
!ZAEROUT3(:,:) =0._JPRB
!ZAEROUT4(:,:) =0._JPRB
!ZAEROUT5(:,:) =0._JPRB

!
! RCHG -> FIXME be careful. This subroutine has KIDIA, KFDIA in the 
!         arguments so probably it is designed for parallel grid. 
!         this initialization may overwrite zeros the proper way might
!         be ARRAY(KIDIA:KFDIA,:,:) = 0._JPRB 
!
 ZOUT(:,:)    = 0._JPRB
 ZOUT2(:,:)   = 0._JPRB
 ZOUT3(:,:,:) = 0._JPRB
!!! initialize those arrays are importnat for PGFL under GNU, Lianghai Wu
 ZVDEP(:,:) = 0._JPRB ! ddep velocity as zero
 ZXTEMS(:,:) = 0._JPRB ! surface emissions as zero for input
 ZXTMD1(:,:,:) = 0._JPRB 
!ZOUT4(:,:)   = 0._JPRB
!ZOUT5(:,:)   = 0._JPRB
!ZOUT6(:,:)   = 0._JPRB
!ZOUT7(:,:)   = 0._JPRB
!ZOUT8(:,:)   = 0._JPRB
!ZOUT9(:,:)   = 0._JPRB
!ZZOUT1(:,:)  = 0._JPRB
!ZZOUT2(:,:)  = 0._JPRB
!ZZOUT3(:,:)  = 0._JPRB
!ZZOUT4(:,:)  = 0._JPRB
!ZZOUT5(:,:)  = 0._JPRB
!ZZOUT6(:,:)  = 0._JPRB
!ZZOUT7(:,:)  = 0._JPRB
!ZZOUT8(:,:)  = 0._JPRB
!ZZOUT9(:,:)  = 0._JPRB
!ZZOUT10(:,:) = 0._JPRB
!ZZOUT11(:,:) = 0._JPRB
!ZZOUT12(:,:) = 0._JPRB

M7TEND_IN(:,:,:)   = 0._JPRB
M7TEND_OUT(:,:,:)  = 0._JPRB
SEDOUT(:,:,:)      = 0._JPRB
DDEPOUT(:,:,:)     = 0._JPRB
WDEPOUT(:,:,:)     = 0._JPRB
ZSEDIFLUX(:,:,:)   = 0._JPRB
ZDDEPFLUX(:,:)     = 0._JPRB
ZDDEPFLUX_SO2(:)   = 0._JPRB
ZVDA(:,:)          = 0._JPRB
SEDOUT_2D(:,:)     = 0._JPRB
DDEPOUT_2D(:,:)    = 0._JPRB
WDEPOUT_2D(:,:)    = 0._JPRB

zaveragep(:,:,:)   = 0.0_JPRB
zm7kappa(:,:,:)    = 0.0_JPRB
zh2so4cs(:,:,:)    = 0.0_JPRB
zm7prodcond(:,:,:) = 0.0_JPRB

ZTAERO(:,:,:)      = 0._JPRB

ZCEN(:,:,:) = 0._JPRB

!ZAERSRC(KIDIA:KFDIA,1:NACTAERO)=PAERSRC(KIDIA:KFDIA,1:NACTAERO) 

! computation of tropopause level 
CALL TROPLEV(KLON,KIDIA,KFDIA,KLEV,.FALSE.,PTP,PQP,PRSF1,IKLEVTROP)

! Initializing tracer number and mixing ratios and gas concentrations to not be tendency updated values
ZCEN(KIDIA:KFDIA,1:KLEV,:) = 0._JPRB
ITRC=0
DO JEXT=1,NGHG
  ITRC=ITRC+1
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZCEN(JL,JK,ITRC) = PGFL(JL,JK,YGHG(JEXT)%MP9_PH)
    ENDDO
  ENDDO
ENDDO

DO JEXT=1,ZTRAC
  ITRC=ITRC+1
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
       ZCEN(JL,JK,ITRC) = PGFL(JL,JK,YTRAC(JEXT)%MP9_PH)
    ENDDO
  ENDDO
ENDDO

DO JEXT=1,NAERO
  ITRC=ITRC+1
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZCEN(JL,JK,ITRC) = PGFL(JL,JK,YAERO(JEXT)%MP9_PH)
    ENDDO
  ENDDO
ENDDO

IF(LAERCHEM) then 
  DO JEXT=1,NCHEM
    ITRC=ITRC+1
    DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
        ZCEN(JL,JK,ITRC) = PGFL(JL,JK,YCHEM(JEXT)%MP9_PH)
      ENDDO
    ENDDO
  ENDDO
ENDIF

!DO JAER=1,NACTAERO
!  DO JK=1,KLEV
!    DO JL=KIDIA,KFDIA
!      ! check on max-values here to prevent excessive values in OPTICS_AOP_GET.. probably to be removed again!!
!      ZAEROK(JL,JK,JAER) =MIN(ZCEN (JL,JK,KAERO(JAER)), 1e10)
!      ZTAEROK(JL,JK,JAER)=PTENC(JL,JK,KAERO(JAER))
!    ENDDO
!  ENDDO
!ENDDO

! Initialize output aerosol tendencies
DO JAER=1,NACTAERO
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZAEROK(JL,JK,JAER) =ZCEN (JL,JK,KAERO(JAER))
      ZTAEROK(JL,JK,JAER)=PTENC(JL,JK,KAERO(JAER))
    ENDDO
  ENDDO
ENDDO

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZRHO(JL,JK) = PRSF1(JL,JK)/(RD*PTP(JL,JK))
    ZDP(JL,JK)  = PRS1(JL,JK) - PRS1(JL,JK-1)
    ZDPG(JL,JK) = (PRS1(JL,JK) - PRS1(JL,JK-1))/RG ! calculate dp/g
    ZDZ(JL,JK)  = ZDP(JL,JK) / (ZRHO(JL,JK)*RG)
  ENDDO
ENDDO

ZM6RP(KIDIA:KFDIA,:,:)  = 0.0_JPRB
ZM6DRY(KIDIA:KFDIA,:,:) = 0.0_JPRB
ZRHOP(KIDIA:KFDIA,:,:)  = 0.0_JPRB
ZWW(KIDIA:KFDIA,:,:)    = 0.0_JPRB
!ZAERML(KIDIA:KFDIA,:,:)=0.0_JPRB
!ZAERNL(KIDIA:KFDIA,:,:)=0.0_JPRB

IF (LCHEM_DIA) THEN
  ZTAERO0(KIDIA:KFDIA,1:KLEV,1:NACTAERO) =  ZTAEROK(KIDIA:KFDIA,1:KLEV,1:NACTAERO)
  ZTENC0(KIDIA:KFDIA,1:KLEV, :) = 0._JPRB
ENDIF 

! RCHG -> I would move this block to a subroutine (here after a "contains") 
!*         0.2    Preliminary computation
!*         GAS-TO-PARTICLE CONVERSION (SO2 -> SO4)
!*         (SO2 -> SO4)

IF(.NOT. LAERCHEM) THEN

  DO JGAS=1,2
    IF (TRIM(YAERO(ind_oifs_ham%ind_gas_OIFS(JGAS))%CNAME)=='SO2') THEN
      ISSO2=ind_oifs_ham%ind_gas_OIFS(JGAS)
    ELSE IF (TRIM(YAERO(ind_oifs_ham%ind_gas_OIFS(JGAS))%CNAME)=='SO4_gas')THEN
      ISSO4=ind_oifs_ham%ind_gas_OIFS(JGAS)
    ELSE
      CALL ABOR1('ABORT: IN AER_PHY3_layer, SO2 not defined. Wrong table in use')
    END IF
  END DO
  
  DO JAER=1,NACTAERO
    IF (TRIM(YAERO(JAER)%CNAME)=='SO4_AS') THEN
      ISSO4_ACS=JAER
      EXIT
    END IF
  END DO

  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
 !!!ZDP(JL,JK)= PAUX%PRS1(JL,JK) - PAUX%PRS1(JL,JK-1)
      ZTSO2(JL,JK)    = PTENC(JL,JK,ISSO2)
      ZTSO4(JL,JK, 1) = PTENC(JL,JK,ISSO4)
      ZTSO4_AQ(JL,JK) = 0.0_JPRB
      ZSO2(JL,JK)     = PAEROP(JL,JK,ISSO2)
      ZITSO2(JL,JK)   = PTENC(JL,JK,ISSO2)
    END DO
  END DO

  !CALL AER_SO2SO4_V2 ( &
  !     KIDIA, KFDIA, KLON, KLEV, &
  !     ! TSPHY, STATE%T, PAUX%PRSF1 , PRAD%PNEB, PRAD%PQLI, PAUX%PGELAM,&
  !     ! take the liquidwater and cloud fraction from state variables i
  !     ! instead of radiation data structures, in 43r3 PRAD% are zero
  !     ! between radiation time steps
  !     TSPHY, STATE%T, PAUX%PRSF1 , STATE%A, STATE%CLD(:,:,NCLDQL), PAUX%PGELAM,&
  !     ZSO2, ZITSO2, PGFL(:,:,YAERCLIM(1)%MP),&
  !     PGFL(:,:,YAERCLIM(2)%MP), PGFL(:,:,YAERCLIM(3)%MP) ,&
  !     ZTSO2 , ZTSO4, ZTSO4_AQ, ZFSO2, ZFSO4, ZFSO4_AQ, ZDP )
  
  CALL AER_SO2SO4_V2( YDRIP,                                                &
                    & KIDIA , KFDIA , KLON  , KLEV  ,                       &
                    & PTSPHY, PTP   , PRSF1 , PAP  , PLP  , PGELAT, PGELAM, &
                    & ZSO2  , ZITSO2,                                       &
                    & PGFL(:,:, YGFL%YAEROCLIM(1)%MP),                      &
                    & PGFL(:,:,YGFL%YAEROCLIM(2)%MP),                       &
                    & PGFL(:,:,YGFL%YAEROCLIM(3)%MP) ,                      &
                    & ZTSO2 , ZTSO4(:,:,1), ZTSO4_AQ, ZFSO2, ZFSO4, ZFSO4_AQ, ZDP)

  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ! ZTSO2
      PTENC(JL,JK,ISSO2)=PTENC(JL,JK,ISSO2)+ZTSO2(JL,JK)
      ! ZTSO4 
      PTENC(JL,JK,ISSO4)=PTENC(JL,JK,ISSO4)+ZTSO4(JL,JK,1)
      ! SO4 formed in clouds is applied to Accumulati 
      !write(7700,*)ISSO4_ACS,shape(GEMSL%ZTENC(JL,JK,:))
      PTENC(JL,JK,ISSO4_ACS)=PTENC(JL,JK,ISSO4_ACS)+ZTSO4_AQ(JL,JK)
      !This is done insed aer_phy3.F90:
      !       ELSE IF(AERO_SCHEME == "aer") THEN
      !          GEMSL%ZTENC(JL,JK,ISSO4)=GEMSL%ZTENC(JL,JK,ISSO4)+ZTSO4(JL,JK)
      
    END DO
  END DO

ENDIF
!end of *         (SO2 -> SO4)


!RCHG -> I would move this to a second subroutine (again after "CONTAINS") 
!
!*         1.1    COMPUTE RELATIVE HUMIDITY WITHOUT VERTICAL SMOOTING
!                 ---------------------------------------------------
! Q at saturation for RH calculation
IFLAG=2
CALL SATUR(KIDIA , KFDIA , KLON  , KTDIA , KLEV,  YDMODEL%YRML_PHY_SLIN%YREPHLI%LPHYLIN, &
  & PRSF1, PTP    , ZQSAT , IFLAG)  

! RH calculation
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZQP(JL,JK)=MAX(0.0_JPRB,PQP(JL,JK)) !add treshold for full level humid
    ZRH(JL,JK)=ZQP(JL,JK)/(MAX(1.E-30_JPRB,ZQSAT(JL,JK)))
    ZRH(JL,JK)=MIN(1.0_JPRB,MAX(0.0_JPRB,ZRH(JL,JK)))
    ZAP(JL,JK)=MIN(1.0_JPRB,MAX(0.0_JPRB,PAP(JL,JK))) !add threshold for cloud cover
  ENDDO
ENDDO

! TB apparently unnecessary in current implementation, but ISSO4_C still needed for chem_inext in the code.
! needs to be reviewed if it can be removed. 
IF(LAERCHEM .AND. NCHEM>0) THEN
  DO JT=1,NCHEM
    IF(TRIM(YCHEM(JT)%CNAME)== 'SO4' ) THEN
      ISSO4_C=KCHEM(JT)
      ISO4_C=JT
    ENDIF
  ENDDO
  ZSO4G(KIDIA:KFDIA,1:KLEV)=ZCEN(KIDIA:KFDIA,1:KLEV,ISSO4_C)
ELSE
  ZSO4G(KIDIA:KFDIA,1:KLEV)=0._JPRB
ENDIF


! Convert from kg kg-1 to molec cm-3

IF (LAERCHEM) THEN
   ZELVOC(KIDIA:KFDIA,1:KLEV)=ZCEN(KIDIA:KFDIA,1:KLEV,KAERO(ielvoc))
   ZSVOC(KIDIA:KFDIA,1:KLEV)=ZCEN(KIDIA:KFDIA,1:KLEV,KAERO(iisvoc))
ELSE
   ZELVOC(KIDIA:KFDIA,1:KLEV)=0.0_JPRB
   ZSVOC(KIDIA:KFDIA,1:KLEV)=0.0_JPRB
END IF

! RCHG => before it has ZICNC(:,:)=0._JPRB which is a semantic error.
!calculate ICNC
ZICNC(KIDIA:KFDIA,1:KLEV) = 0._JPRB
ZICNC(KIDIA:KFDIA,1:KLEV) = RNICE
!Put tracer mixing ratios and tendencies to from OIFS to HAM
!initialize to zero
ZXTM1(KIDIA:KFDIA,1:KLEV,:) = 0._JPRB
ZXTTE(KIDIA:KFDIA,1:KLEV,:) = 0._JPRB
ZXTTEM1(KIDIA:KFDIA,1:KLEV,:) = 0._JPRB

!number
DO JCLASS=1,nclass
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZXTM1(JL,JK,ind_oifs_ham%ind_class_HAM(JCLASS)) = ZCEN(JL,JK,KAERO(ind_oifs_ham%ind_class_OIFS(JCLASS)))
      ZXTTE(JL,JK,ind_oifs_ham%ind_class_HAM(JCLASS)) = PTENC(JL,JK,KAERO(ind_oifs_ham%ind_class_OIFS(JCLASS)))
    END DO
  END DO
END DO
!mass
DO JMASS=1,naerocomp
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZXTM1(JL,JK,ind_oifs_ham%ind_mass_HAM(JMASS)) = ZCEN(JL,JK,KAERO(ind_oifs_ham%ind_mass_OIFS(JMASS)))
      ZXTTE(JL,JK,ind_oifs_ham%ind_mass_HAM(JMASS)) = PTENC(JL,JK,KAERO(ind_oifs_ham%ind_mass_OIFS(JMASS)))
      ! in case of simple sulfur scheme add SO4_AQ part into SO4_ACS
      ! both original tendency and m7tendency

    END DO
  END DO
END DO
!gas
DO JGAS=1,subm_ngasspec
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      IF (LAERCHEM) THEN
        ZXTM1(JL,JK,ind_oifs_ham%ind_gas_HAM(JGAS)) = ZCEN(JL,JK,KCHEM(ind_oifs_ham%ind_gas_OIFS(JGAS)))
        ZXTTE(JL,JK,ind_oifs_ham%ind_gas_HAM(JGAS)) = PTENC(JL,JK,KCHEM(ind_oifs_ham%ind_gas_OIFS(JGAS)))
        !ADD SO4 from wet chemistry to tendencies
        IF (ind_oifs_ham%ind_gas_OIFS(JGAS) == ISO4ACS )THEN
          ZXTTE(KIDIA:KFDIA,1:KLEV,ind_oifs_ham%ind_gas_OIFS(JGAS))=ZXTTE(KIDIA:KFDIA,1:KLEV,ind_oifs_ham%ind_gas_OIFS(JGAS))+PCHEM2AER(KIDIA:KFDIA,1:KLEV,2)
        ELSEIF( ind_oifs_ham%ind_gas_OIFS(JGAS) == ISO4COS)THEN
          ZXTTE(KIDIA:KFDIA,1:KLEV,ind_oifs_ham%ind_gas_OIFS(JGAS))=ZXTTE(KIDIA:KFDIA,1:KLEV,ind_oifs_ham%ind_gas_OIFS(JGAS))+PCHEM2AER(KIDIA:KFDIA,1:KLEV,3)
        END IF
      ELSE !simple sulfur scheme  
        IF(TRIM(YAERO(ind_oifs_ham%ind_gas_OIFS(JGAS))%CNAME)=='SO2')THEN
          ZXTM1(JL,JK,ind_oifs_ham%ind_gas_HAM(JGAS)) = ZCEN(JL,JK,KAERO(ind_oifs_ham%ind_gas_OIFS(JGAS)))
          ZXTTE(JL,JK,ind_oifs_ham%ind_gas_HAM(JGAS)) = PTENC(JL,JK,KAERO(ind_oifs_ham%ind_gas_OIFS(JGAS)))
          ZXTTEM1(JL,JK,ind_oifs_ham%ind_gas_HAM(JGAS)) = PTENC(JL,JK,KAERO(ind_oifs_ham%ind_gas_OIFS(JGAS)))
          
        ELSE IF (TRIM(YAERO(ind_oifs_ham%ind_gas_OIFS(JGAS))%CNAME)=='SO4_gas')THEN
          ZXTM1(JL,JK,ind_oifs_ham%ind_gas_HAM(JGAS)) = ZCEN(JL,JK,KAERO(ind_oifs_ham%ind_gas_OIFS(JGAS)))
          ZXTTE(JL,JK,ind_oifs_ham%ind_gas_HAM(JGAS)) = PTENC(JL,JK,KAERO(ind_oifs_ham%ind_gas_OIFS(JGAS))) 
          ZXTTEM1(JL,JK,ind_oifs_ham%ind_gas_HAM(JGAS)) = PTENC(JL,JK,KAERO(ind_oifs_ham%ind_gas_OIFS(JGAS))) 
        END IF

      END IF
    END DO
  END DO
END DO


! RCHG -> This will produce segmentation fault if CDNC are not in the namelist 
!         we need to test these things and do a CALL ABORT1() 
!         for these cases. For this we need to find the best flags or combinations 
!         of flags. FIXME

! IF (LAERRRTM) THEN
!cloud
DO JCLOUD=1,2 !CDNC and ICNC
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZXTM1(JL,JK,ind_oifs_ham%ind_cloud_HAM(JCLOUD)) = ZCEN(JL,JK,KAERO(ind_oifs_ham%ind_cloud_OIFS(JCLOUD)))
      ZXTTE(JL,JK,ind_oifs_ham%ind_cloud_HAM(JCLOUD)) = PTENC(JL,JK,KAERO(ind_oifs_ham%ind_cloud_OIFS(JCLOUD)))
    END DO
  END DO
END DO
!ENDIF

! implementation of HAM-M7
ZWND(:) = 0._JPRB
DO JL=KIDIA,KFDIA
  ZWND(JL)=MAX(1.0E-10_JPRB,PAERWS(JL)) ! make threshold for wind speed
  ! Compute average fluxes over tiles
  DO JTILE=1,KTILES
    ZAHFSM(JL)=ZAHFSM(JL)+PFRTI(JL,JTILE)*PAHFSTI(JL,JTILE)
  ENDDO
ENDDO
 
! --> calling the correct microphysics scheme
 
SELECT CASE (TRIM(AERO_SCHEME))
  CASE ("tm5m7")
  ! ** OBSOLETE ** (PLS)
  
  !alaak: commented out because model crashes
  !kg(air)-1 to cm-3
!!$
  !do JK =1,KLEV
  !   do JL =KIDIA,KFDIA
  !      ZAIRDM(JL) = (7.24291E16_JPRB*PRSF1(JL,JK)/PTP(JL,JK)) * RMD
  !      !end do
  !      DO JMOD=1,NMOD
  !         ZAERNL(JL,JK,JMOD)=MAX(PAEROP(JL,JK,MODE_START(JMOD))/ZAIRDM(JL),0.0_JPRB)
  !      END DO
  !
  !   end do
  !end do
  
!!$     DO JMOD=1,NMOD
!!$        DO JAERCLASS=1,NAERMOD
!!$           IF (MODE_TRACERS_BY_MODS(JAERCLASS,JMOD)>0) THEN
!!$              ZAERML(KIDIA:KFDIA,1:KLEV,JAERCLASS)=ZCEN(KIDIA:KFDIA,1:KLEV,MODE_TRACERS_BY_MODS(JAERCLASS,JMOD))
!!$           END IF
!!$        END DO
!!$     END DO
!!$
!!$     call M7(KIDIA  , KFDIA  , KLON   , KLEV , &
!!$          ! &  KFLDX, KLEVX,  KTRAC  , KAERO , KCHEM, &
!!$          PRSF1,ZRH,PTP,&
!!$          ZSO4G,ZELVOC, ZSVOC,ZAERML,ZAERNL,&
!!$          ZRHOP,ZWW,ZM6RP,ZM6DRY,&
!!$          PTSPHY)

  CASE("hamm7")
    ! Initializations for submodel interface
    ZGRVOL(KIDIA:KFDIA,1:KLEV) = 1.79e12_JPRB ! ZGRVOL is only used for diagnostics (only when HAMMOZ is on)
    ZPBL = 1 ! boundary layer top = 1 (ITOP=1)
    ! Allocate variables for aerosol processes
    ALLOCATE(ZW(KLON,KLEV,nw))
    ALLOCATE(ZWPDF(KLON,KLEV,nw))
    !IF (.NOT. ALLOCATED(w_large)) ALLOCATE(w_large(KLON,KLEV,ZKROW))
    !IF (.NOT. ALLOCATED(w_turb)) ALLOCATE(w_turb(KLON,KLEV,ZKROW))
    ALLOCATE(ZRC(KLON,KLEV,nclass,nw))
    ALLOCATE(ZSMAX(KLON,KLEV,nw))
    !IF (.NOT. ALLOCATED(reffi)) ALLOCATE(reffi(KLON,KLEV,ZKROW))
    !IF (.NOT. ALLOCATED(reffl)) ALLOCATE(reffl(KLON,KLEV,ZKROW))
    DO IMODE=1,NMOD
      ALLOCATE(RW_MODE (IMODE)%d2(KLON,KLEV))
      ALLOCATE(DENS_MODE(IMODE)%d2(KLON,KLEV))
      IF (sizeclass(IMODE)%lsoluble) THEN
        ALLOCATE(RWD_MODE(IMODE)%d2(KLON,KLEV))
        ALLOCATE(H2O_MODE(IMODE)%d2(KLON,KLEV))
      END IF
    ENDDO
    ! End allocate variables for aerosol processes
    
    !-----------------------------------------------------------------
    ! Submodel interface call (HAM aerosol microphysics)
    
    ! changes to particles:
    ! Update SO4ACS and SO4OCS tendencies to include the contribution from wet chemistry
    IF (LAERCHEM) THEN
      ZTAEROK(KIDIA:KFDIA,1:KLEV,ISO4ACS)=ZTAEROK(KIDIA:KFDIA,1:KLEV,ISO4ACS)+PCHEM2AER(KIDIA:KFDIA,1:KLEV,2)
      ZTAEROK(KIDIA:KFDIA,1:KLEV,ISO4COS)=ZTAEROK(KIDIA:KFDIA,1:KLEV,ISO4COS)+PCHEM2AER(KIDIA:KFDIA,1:KLEV,3)
    END IF

    CALL GSTATS(2501,0)

    CALL HAM_SUBM_INTERFACE(&
         & KFDIA, KLON,  KLEV, ZKROW, & !dimension indices
         & ntrac, PRSF1, PRS1,        & !number of tracers, pressure full levels, pressure half levels
         & PTP,   ZQP,   ZQSAT,       & !temperature, specific humidity, saturation spec. hum.
         & ZXTM1, ZXTTE,              & !tracer mass/number mr, tracer tendencies
         & ZM6RP, ZM6DRY,             & !mean mode actual radius [m], dry radius for soluble modes [m] 
         & ZRHOP, ZWW,                & !mean mode particle density [kg m-3], aerosol water content for each mode [kg(water)-3(air)]
         & ZAP,   ZGRVOL,             & !cloud fraction, grid box volume (only for diagnostics)
         & ZPBL,  zout3)                !boundary layer top level, outputs
    
    CALL GSTATS(2501,1)

    ! updating wet and dry radii, aerosol density and water content
    DO IMODE=1,NMOD
      RW_MODE(IMODE)%d2(KIDIA:KFDIA,1:KLEV)=ZM6RP(KIDIA:KFDIA,1:KLEV,IMODE) ! m ( , KLEV, NMOD)
      DENS_MODE(IMODE)%d2(KIDIA:KFDIA,1:KLEV)=ZRHOP(KIDIA:KFDIA,1:KLEV,IMODE) ! kg/m3
    ENDDO
    DO IMODE=1,NMOD
      IF (sizeclass(IMODE)%lsoluble) THEN
        RWD_MODE(IMODE)%d2(KIDIA:KFDIA,1:KLEV) = ZM6DRY(KIDIA:KFDIA,1:KLEV,IMODE) ! m
        H2O_MODE(IMODE)%d2(KIDIA:KFDIA,1:KLEV) = ZWW(KIDIA:KFDIA,1:KLEV,IMODE) ! ?
      END IF
    ENDDO
    !-----------------------------------------------------------------

    CALL GSTATS(2502,0)

    !---calculate updraft velocity
    ZTKEM1(KIDIA:KFDIA,1:KLEV) = 0._JPRB ! turbulent kinetic energy as zero for now as it is not used (YET!)
    ZWCAPE(KIDIA:KFDIA) = 0._JPRB !  CAPE as zero as it is not used

    CALL activ_updraft(KFDIA, KLON, KLEV, ZKROW, & ! krow = 1
         ZTKEM1, ZWCAPE, PVERVEL, ZRHO, & ! turbulent kinetic energy, CAPE contr. to conv. vert. veloc. [m s-1], large scale vert. veloc.
         ZW, ZWPDF)

    IF (ncd_activ == 2 .OR. nccndiag > 0) THEN
      CALL ham_activ_koehler_ab(KFDIA, KLON, KLEV, ZKROW, KTDIA, & ! krow=1 ktdia=1
           ZXTM1, PTP, ZA, ZB)
    END IF

    DO JCLASS = 1,nclass! nclass=7
      IF (sizeclass(jclass)%lsoluble) THEN 
        ZRDRY(KIDIA:KFDIA,1:KLEV,JCLASS) = ZM6DRY(KIDIA:KFDIA,1:KLEV,JCLASS) !soluble modes rdry from rdry_m7
      ELSE
        ZRDRY(KIDIA:KFDIA,1:KLEV,JCLASS) = ZM6RP(KIDIA:KFDIA,1:KLEV,JCLASS)  !insoluble modes rdry from rwet_m7
      END IF
    END DO

    !calculate saturation water vapor pressure from saturation specific humidity
    DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
        ZESW(JL,JK)=(ZQSAT(JL,JK)*PRSF1(JL,JK))/(0.62198_JPRB)
      ENDDO
    ENDDO

    CALL HAM_ACTIV_ABDULRAZZAK_GHAN( KFDIA, KLON, KLEV, ZKROW, KTDIA,  & ! in original 1 = ktdia... for diagnostics so krow=1 and ktdia=1
                                   & ZCDNCACT, ZESW, ZRHO,             & ! number of activated particles, saturation vapor pressure, air density
                                   & ZXTM1, PTP, PRSF1, ZQP,           & ! tracer mix rat, temperature, air pressure, spec. humid.
                                   & ZW, ZWPDF, ZA, ZB, ZRDRY,         & ! mean udr veloc, pdf of udr. veloc, Koehler A, Koehler B, dry radius
                                   & ZNACT, ZFRACN, ZSC, ZRC, ZSMAX)     ! num. act. part., frac ", crit. ssat., crit. radius, max ssat
    CALL GSTATS(2502,1)
    
    !<-- End activation for HAM-M7
    !-----------------------------------------------------------------
    
    
    !<-- Store CDNC (number of activated particles) and ICNC as a number mixing ratio to tracer values and to PGFL fields
    ZXTM1(KIDIA:KFDIA,1:KLEV,idt_cdnc) = ZCDNCACT(KIDIA:KFDIA,1:KLEV)/ZRHO(KIDIA:KFDIA,1:KLEV) ! [#/kg]
    ZXTM1(KIDIA:KFDIA,1:KLEV,idt_icnc) = (1.0E6_JPRB)*ZICNC(KIDIA:KFDIA,1:KLEV)/ZRHO(KIDIA:KFDIA,1:KLEV) !ice crystal number conc = #/cm3 --> number mix rat [#/kg]
    PGFL(KIDIA:KFDIA,1:KLEV,YCDNC%MP9_PH) = (1.0E-6_JPRB)*ZCDNCACT(KIDIA:KFDIA,1:KLEV) ! add CDNC to PGFL field (convert from #/m3 to #/cm3)
    PGFL(KIDIA:KFDIA,1:KLEV,YICNC%MP9_PH) = ZICNC(KIDIA:KFDIA,1:KLEV) ! add ICNC to PGFL field (does not need convert)
    !--> End store CDNC and ICNC


    !-----------------------------------------------------------------
    !--> Calculation for effective radii and put to PGFL fields

    ! put default values for effective radii
    reffl(KIDIA:KFDIA,1:KLEV,ZKROW) = 4._JPRB ! comes from liquid effective radius routine (PP_MIN_RE_UM)
    reffi(KIDIA:KFDIA,1:KLEV,ZKROW) = 80._JPRB*0.64952_JPRB ! comes from ice effective radius routine (ZDEFAULT_RE_UM)

    ! liquid effective radius
    DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
        ZTMPA = 1.0_JPRB/MAX(ZAP(JL,JK),ZEPSEC)
        LLIQCLD(JL,JK) = ( PLP(JL,JK)*ZTMPA  ) > ZEPSEC ! logical for liquid cloud
        LICECLD(JL,JK) = ( PIP(JL,JK)*ZTMPA  ) > ZEPSEC ! logical for ice cloud
        ZQLWP(JL,JK) = PLP(JL,JK)/MAX(ZAP(JL,JK),1.E-10_JPRB) ! calculate lwp
        ZQLWP(JL,JK) = MIN(MAX(ZQLWP(JL,JK),0.0_JPRB),RCLDMAX) ! treshold lwp
        ! effective radius calculated similarly as in radlswr.F90
        ! 2.387e-10 is 3/(4*pi*rho_liq*10^6)  [10^6 for N in right units]
        ZRE_LIQ(JL,JK) = 1.E+06_JPRB*(2.387e-10_JPRB*ZRHO(JL,JK)*ZQLWP(JL,JK)/(MAX(PGFL(JL,JK,YCDNC%MP9_PH),10._JPRB)))**0.333_JPRB ! calculate effective radius in um (use minimum value for CDNC if CDNC is small)
      END DO
    END DO

    ! Add liq. eff. rad. to HAM variables (only if there is liquid cloud else minimum value)
    REFFL(KIDIA:KFDIA,1:KLEV,ZKROW) = MERGE(ZRE_LIQ(KIDIA:KFDIA,1:KLEV),4._JPRB,LLIQCLD(KIDIA:KFDIA,1:KLEV))
    CALL ICE_EFFECTIVE_RADIUS(YRERAD, YDSPP_CONFIG, KIDIA, KFDIA, KLON, KLEV, &
         &  PRSF1, PTP, ZAP, PIP, PSP, PGEMU, & ! pressure, temp, cloud fr., IWC, SWC, sine of latitude
         &  reffi(1:KLON,1:KLEV,ZKROW)) ! ice effective radius (updated to mo_activ variable 'reffi' which used in mo_ham_wetdep)

    ! only if there is ice cloud else minimum value
    REFFI(KIDIA:KFDIA,1:KLEV,ZKROW) = MERGE(REFFI(KIDIA:KFDIA,1:KLEV,ZKROW), 20._JPRB, LICECLD(KIDIA:KFDIA,1:KLEV))

    ! add effective radii to PGFL fields
    PGFL(KIDIA:KFDIA,1:KLEV,YRE_LIQ%MP9_PH) = 1.0E-06_JPRB * reffl(KIDIA:KFDIA,1:KLEV,ZKROW) ! convert um to meters and save to PGFL fields
    PGFL(KIDIA:KFDIA,1:KLEV,YRE_ICE%MP9_PH) = 1.0E-06_JPRB * reffi(KIDIA:KFDIA,1:KLEV,ZKROW) ! convert um to meters and save to PGFL fields

    !<-- End calculation for effective radii
    !-----------------------------------------------------------------

    !--- Mass conserving correction of negative tracer values:
    CALL xt_borrow(KFDIA, KLON,  KLEV, KLEV+1, ntrac, &
         PRSF1, PRS1, &
         ZXTM1, ZXTTE)

    ! RCHG -> This is wetdep interface it can be a subroutine afer "CONTAINS" 
    !
    !--> Wet deposition for HAM-M7
    CALL GSTATS(2503,0)

    IF ( LAERSCAV ) THEN

      !--> initialize mixing ratios for wet deposition
      !-- initialise in-cloud and interstitial mixing ratios
      !   set both equal to tracer mixing ratio as starting point
      !   ham_wet_chemistry will re-compute these values if lham=true
      DO JT = 1,NTRAC
        ZXTP1(KIDIA:KFDIA,1:KLEV,JT)  = ZXTM1(KIDIA:KFDIA,1:KLEV,JT) + ZXTTE(KIDIA:KFDIA,1:KLEV,JT) * time_step_len
        ZXTP1C(KIDIA:KFDIA,1:KLEV,JT) = ZXTP1(KIDIA:KFDIA,1:KLEV,JT)
        ZXTP10(KIDIA:KFDIA,1:KLEV,JT) = ZXTP1(KIDIA:KFDIA,1:KLEV,JT)
      END DO

      !<-- call wetdep interface for wet deposition
      !-- interface to wet deposition routine (also from cuflx_subm)
      IF ( lwetdep .AND. ANY(trlist%ti(:)%nwetdep > 0) ) THEN

        ! for calculating the rain and snow evaporation/formation variables used in wet deposition
        DO JK=1,KLEV
          DO JL=KIDIA,KFDIA
            ZFLXR=PFPLCL(JL,JK-1)
            ZFLXS=PFPLCN(JL,JK-1)
            ZFLXRB=PFPLCL(JL,JK)
            ZFLXSB=PFPLCN(JL,JK)
            IF (PCOVPTOT(JL,JK) > 1e-40_JPRB) THEN
              ZMRATEPR(JL,JK) = ( ZFLXRB-ZFLXR ) / PCOVPTOT(JL,JK)
              ZMRATEPS(JL,JK) = ( ZFLXSB-ZFLXS ) / PCOVPTOT(JL,JK)
              !same formula negatives/positives for evap or formation
              ZFEVAPR(JL,JK) = ( ZFLXRB-ZFLXR ) / PCOVPTOT(JL,JK)
              ZFSUBLS(JL,JK) = ( ZFLXSB-ZFLXS ) / PCOVPTOT(JL,JK)
            ELSE
              ZMRATEPR(JL,JK) = 0.0_JPRB
              ZMRATEPS(JL,JK) = 0.0_JPRB
              !same formula negatives/positives for evap or formation
              ZFEVAPR(JL,JK) = 0.0_JPRB
              ZFSUBLS(JL,JK) = 0.0_JPRB
            END IF
          END DO
        END DO
        ZMSNOWACL(KIDIA:KFDIA,1:KLEV) = ZMRATEPS(KIDIA:KFDIA,1:KLEV) !?

        ZDUMMY(:,:) = 0._JPRB ! initialize dummy variables for conv. case only
        ZDUM2D(:,:) = 0._JPRB ! initialize dummy variables for conv. case only
        ZDUM3D(:,:,:) = 0._JPRB ! initialize dummy variables for conv. case only
        ZLFRAC_SO2(:,:) = 0._JPRB ! zlfrac_so2 only needed in gas scavenging and this is off for now (put this zero)

        ZLP(KIDIA:KFDIA,1:KLEV) = PLP(KIDIA:KFDIA,1:KLEV) ! temporary variable for cloud water content (modified in wetdep)
        ZIP(KIDIA:KFDIA,1:KLEV) = PIP(KIDIA:KFDIA,1:KLEV) ! temporary variable for cloud ice water content (modified in wetdep)
        ZTENCIH(:,:,:) = 0._JPRB
        ZTENCIH(KIDIA:KFDIA,1:KLEV,1:NTRAC)=ZXTTE(KIDIA:KFDIA,1:KLEV,1:NTRAC)
        IF (.NOT.LAERCHEM)THEN
          CALL HAM_CONV_LFRAQ_SO2(KFDIA,KLON,KLEV,PTP,zxtm1,zrho,ZLP,zlfrac_so2)
        END IF

        CALL WETDEP_INTERFACE(KFDIA, KLON, KLEV, 1, ZKROW, LSTRAT, & !eehol: ktop = 1 (top level index), lstrat = TRUE for strat. case
                ZDPG,  ZMRATEPR, ZMRATEPS,   ZMSNOWACL,          & !eehol: dp/g, evap. of rain, subl. of snow, accr. rate of snow with cl. drop in-cl.
                ZLP,  ZIP,                                       & !eehol: cloud water content, cloud ice water content
                ZM6RP,  ZM6DRY,                                  & ! m7 aerosol: to replace rwet_m7, dry radius for soluble modes [cm]
                reffi,  reffl,                                   & ! m7 aerosol
                ZNACT, ZFRACN,                                   &
                PTP, ZXTM1, ZLFRAC_SO2, ZXTTE, ZXTP10, ZXTP1C,   & !eehol: temperature, zlfrac_so2 only needed in gas scavenging (0 for now)
                PFPLSL, PFPLSN, ZFEVAPR, ZFSUBLS,                & !eehol: rain flux, snow flux, 
                ZDUM2D, ZDUM3D,                                  & !eehol: as zero as these are not needed in strat. case
                ZAP,  PCOVPTOT, ZRHO, ZDUMMY)  !eehol: cloud frac., precip. frac., air dens.

        DO JMASS=1,NAEROCOMP
          JY=KAERO(ind_oifs_ham%ind_mass_OIFS(JMASS))
          WDEPOUT_2D(KIDIA:KFDIA, JY) = WDEPOUT_2D(KIDIA:KFDIA,JY) + ZDUMMY(KIDIA:KFDIA,ind_oifs_ham%ind_mass_ham(JMASS))
        END DO

        DO JCLASS=1,NCLASS
          JY=KAERO(ind_oifs_ham%ind_class_OIFS(JCLASS))
          WDEPOUT_2D(KIDIA:KFDIA,JY) = WDEPOUT_2D(KIDIA:KFDIA,JY) + ZDUMMY(KIDIA:KFDIA, ind_oifs_ham%ind_class_ham(JCLASS))
        END DO

      END IF
    END IF
    CALL GSTATS(2503,1)

    !<-- End wet deposition for HAM-M7
    !-----------------------------------------------------------------

    !--- Mass conserving correction of negative tracer values:
    CALL XT_BORROW(KFDIA, KLON,  KLEV, KLEV+1, ntrac, &
         PRSF1, PRS1, &
         ZXTM1, ZXTTE)

    !-----------------------------------------------------------------
    !--> Sedimentation for HAM-M7
    CALL GSTATS(2504,0)

    IF (LAERSEDIM) THEN
      IF ( lsedimentation .AND. ANY(trlist%ti(:)%nsedi > 0) ) THEN

        ZTENCIH(:,:,:) = 0._JPRB
        ZTENCIH(KIDIA:KFDIA,1:KLEV,1:ntrac)=ZXTTE(KIDIA:KFDIA,1:KLEV,1:ntrac)     

        CALL SEDI_INTERFACE(KLON, KFDIA, KLEV, ZKROW,   &
             PTP, ZQP, PRSF1, PRS1, & !eehol: temperature, specific humidity, pressure at full level, pressure at half level
             ZM6RP, ZRHOP, & !mean mode actual radius [m], mean mode particle density [kg m-3]
             ZXTM1, ZXTTE, ZSEDIFLUX) !eehol: tracer mixing ratios and tendency (sediflux for diagnostics)

        SEDOUT(KIDIA:KFDIA, 1:KLEV,:)=(ZTENCIH(KIDIA:KFDIA, 1:KLEV,:)-ZXTTE(KIDIA:KFDIA, 1:KLEV,:))
        DO JK=1,KLEV
          DO JCLASS=1,NCLASS
            SEDOUT_2D(KIDIA:KFDIA,KAERO(ind_oifs_ham%ind_class_OIFS(JCLASS)))=SEDOUT_2D(KIDIA:KFDIA,KAERO(ind_oifs_ham%ind_class_OIFS(JCLASS))) + ZSEDIFLUX(KIDIA:KFDIA, JK,ind_oifs_ham%ind_class_HAM(JCLASS))
          END DO
          DO JMASS=1,NAEROCOMP
            SEDOUT_2D(KIDIA:KFDIA,KAERO(ind_oifs_ham%ind_mass_OIFS(JMASS)))=SEDOUT_2D(KIDIA:KFDIA,KAERO(ind_oifs_ham%ind_mass_OIFS(JMASS))) + ZSEDIFLUX(KIDIA:KFDIA,JK,ind_oifs_ham%ind_mass_HAM(JMASS))
          END DO
        END DO
      END IF
    ENDIF
    CALL GSTATS(2504,1)

    !<-- End sedimentation for HAM-M7
    !-----------------------------------------------------------------
    
    !--- Mass conserving correction of negative tracer values:
    CALL xt_borrow(KFDIA, KLON,  KLEV, KLEV+1, ntrac, &
         PRSF1, PRS1, &
         ZXTM1, ZXTTE)

    !-----------------------------------------------------------------
    !--> Dry deposition for HAM-M7
    IF (.NOT.LAERSURF) THEN
      PAERDDP(:,:)  =0._JPRB
    ELSEIF (LAERSURF) THEN

      !*             DRY DEPOSITION INCLUDED AS MODIFICATION TO SURFACE FLUXES
      !              ---------------------------------------------------------
      CALL GSTATS(2505,0)

      IF (LAERDRYDP) THEN

        !--> variables not needed for aerosol dry deposition
        ZCFML(:) = 0._JPRB
        ZCFMW(:) = 0._JPRB
        ZCFMI(:) = 0._JPRB
        ZCFNCL(:) = 0._JPRB
        ZCFNCW(:) = 0._JPRB
        ZCFNCI(:) = 0._JPRB
        ZEPDU2 = 0._JPRB
        ZKAP = 0._JPRB
        ZGEOM1(:,:) = 0._JPRB
        ZRIL(:) = 0._JPRB
        ZRIW(:) = 0._JPRB
        ZRII(:) = 0._JPRB
        ZTVIR1(:,:) = 0._JPRB
        ZTVL(:) = 0._JPRB
        ZTVW(:) = 0._JPRB
        ZTVI(:) = 0._JPRB
        ZAZ0(:) = 0._JPRB
        ZFRL(:) = 0._JPRB
        ZSRFL(:) = 0._JPRB
        ZFOREST(:) = 0._JPRB
        ZTSI(:) = 0._JPRB
        ZAZ0L(:) = 0._JPRB
        ZAZ0I(:) = 0._JPRB
        ZCDNI(:) = 0._JPRB

        !--> variables calculated for dry deposition
        DO JL = KIDIA,KFDIA
          IF ( PLSM(JL) < 0.99_JPRB ) THEN
            ZLOLAND(JL) = .FALSE.
            ZAZ0W(JL) = PZ0M(JL)
          ELSE
            ZLOLAND(JL) = .TRUE.
            ZAZ0W(JL) = 0._JPRB
          END IF
          ZAZ0W(JL)  = MAX(1.0E-5_JPRB,ZAZ0W(JL))  ! treshold roughness length to min value
          ZFRW(JL)   = MAX(0.,1.-PLSM(JL)-PCI(JL)) ! water fraction = 1 - land mask - sea ice fraction
          ZCVS(JL)   = PFRTI(JL,5)+PFRTI(JL,7)     ! snow cover fraction = Snow on low-veg + snow on bare-soil + snow under high-veg
          ZCVW(JL)   = PFRTI(JL,3)                 ! wet skin fraction
          ZVGRAT(JL) = PCVL(JL)+PCVH(JL)           ! vegetation ratio = low veg. cover + high veg. cover
          ZCDNL(JL)  = PAERUST(JL)                 ! adding ustar to not used variable
          ZCDNW(JL)  = LOG(ZDZ(JL,KLEV)/PZ0M(JL))/(VKARMAN*PAERUST(JL)) ! calculate aerodyn. resistance on surface to not used variable
        END DO
        
        !--> init values
        ZTENCIH(:,:,:) = 0._JPRB
        ZTENCIH(KIDIA:KFDIA,1:KLEV,:) = ZXTTE(KIDIA:KFDIA,1:KLEV,:) ! init tendency before drydep
        ZXTEMS(KIDIA:KFDIA,:)         = 0._JPRB                     ! surface emissions as zero for input
        ZXTMD1(KIDIA:KFDIA,1:KLEV,:)  = 0._JPRB
        ZXTMD1(KIDIA:KFDIA,1:KLEV,:)  = ZXTM1(KIDIA:KFDIA,1:KLEV,:) + (ZTENCIH(KIDIA:KFDIA,1:KLEV,:) * time_step_len) ! update mixrat with tendency
        ZVDEP(KIDIA:KFDIA,:)          = 0._JPRB                     ! ddep velocity as zero

        ! RCHG: Recommendation, those subroutines specific of m7 should have 
        !       m7 in the name not sure if this is specific or general/common 
        !       but adapted to m7. Like m7_simple_sulfur_drydep below. 
        CALL DRYDEP_INTERFACE(KLON, KFDIA,  KLEV, ZKROW,                    &
            & ZQP(:,KLEV), ZQSAT(:,KLEV), PTP(:,KLEV), ZCFML, ZCFMW, ZCFMI, &
            & ZCFNCL, ZCFNCW, ZCFNCI,                                       &
            & ZEPDU2, ZKAP, PUP, PVP, ZGEOM1, ZRIL, ZRIW,                   &
            & ZRII,                                                         &
            & ZTVIR1, ZTVL, ZTVW, ZTVI, ZAZ0,                               &
            & PTP(:,KLEV), ZLOLAND,                                         &
            & ZM6RP, ZRHOP,                                                 & ! M7
            & ZFRL,   ZFRW,  PCI,     ZCVS,   ZCVW,     ZVGRAT,             &
            & ZSRFL,  PUP(:,KLEV),     PVP(:,KLEV),                         & !eehol: FIXME 10m u and v wind from lowest level.. needs to be revised in future!!
            & ZXTEMS, ZXTMD1, ZRHO(:,KLEV), PRS1, ZFOREST, ZTSI,            & !air dens lowest, air press at int.
            & ZAZ0L, ZAZ0W, ZAZ0I, ZCDNL, ZCDNW, ZCDNI, ZDDEPFLUX, ZVDEP)     !ZCDNL and ZCDNW used for ustar and aerodyn. resist.
           
        IF (.NOT. LAERCHEM) THEN
          CALL M7_SIMPLE_SULFUR_DRYDEP(YDMODEL, KIDIA,KFDIA, KLON, KLEV, &
               Zxtm1, PCFLX(:,KAERO(1):KAERO(NACTAERO)),  &
               ZDP, PGEOH, ZRHO, ZXTTE, PTSPHY,&
               PSO2DD, PGELAM, &
               ZFAERO, Zxtp1, ZDDEPFLUX_SO2)
          ZDDEPFLUX(KIDIA:KFDIA,2)=ZDDEPFLUX_SO2(KIDIA:KFDIA)
        END IF
        
        !--> modify tendency at surface according to changes in surface emissions           
        DO JT = 1,NTRAC
          DO JL = KIDIA,KFDIA
            ZXTTE(JL,KLEV,JT) = ZTENCIH(JL,KLEV,JT) + ((ZXTEMS(JL,JT)*RG)/(ZDP(JL,KLEV)))
          END DO
        END DO

      ENDIF ! LAERDRYDP
    END IF
    CALL GSTATS(2505,1)

    !<-- End dry deposition for HAM-M7
    !-----------------------------------------------------------------

    !--- Mass conserving correction of negative tracer values:
    CALL xt_borrow(KFDIA, KLON,  KLEV, KLEV+1, ntrac, &
         PRSF1, PRS1, &
         ZXTM1, ZXTTE)

    !-----------------------------------------------------------------
    !--> Add HAM modified tendency back to PTENC (OIFS values)

    !number
    DO JCLASS=1,NCLASS
      PTENC(KIDIA:KFDIA,1:KLEV,KAERO(ind_oifs_ham%ind_class_OIFS(JCLASS))) = ZXTTE(KIDIA:KFDIA,1:KLEV,ind_oifs_ham%ind_class_HAM(JCLASS))
    END DO
    !mass
    DO JMASS=1,NAEROCOMP
      PTENC(KIDIA:KFDIA,1:KLEV,KAERO(ind_oifs_ham%ind_mass_OIFS(JMASS))) = ZXTTE(KIDIA:KFDIA,1:KLEV,ind_oifs_ham%ind_mass_HAM(JMASS))
    END DO
    !gas
    IF(LAERCHEM) THEN
      DO JGAS=1,SUBM_NGASSPEC
        PTENC(KIDIA:KFDIA,1:KLEV,KCHEM(ind_oifs_ham%ind_gas_OIFS(JGAS))) = ZXTTE(KIDIA:KFDIA,1:KLEV,ind_oifs_ham%ind_gas_HAM(JGAS))
      END DO
    ELSE
      DO JGAS=1,SUBM_NGASSPEC
        PTENC(KIDIA:KFDIA,1:KLEV,KAERO(ind_oifs_ham%ind_gas_OIFS(JGAS))) = ZXTTE(KIDIA:KFDIA,1:KLEV,ind_oifs_ham%ind_gas_HAM(JGAS))
      END DO
    END IF

    ! RCHG -> not sure best way to solve here. I commented to avoid segmentation fault 
    !         but it may be avoided with other more specific flag. Anyway something was 
    !         needed to avoid core-dump. 
    !         The problem with these arrays in loop below is that the indices:
    !          ind_oifs_ham%ind_cloud_HAM(JCLOUD)) 
    !          KAERO(ind_oifs_ham%ind_cloud_OIFS(JCLOUD)))
    !         are not propoperly set up by hamm7_init.F90
    !         -- this need to be solved probably in hamm7_init.F90 which detect the 
    !         problems with these tracers about CCN. 
    !cloud variables
    DO JCLOUD=1,2 !CDNC and ICNC
      PTENC(KIDIA:KFDIA,1:KLEV,KAERO(ind_oifs_ham%ind_cloud_OIFS(JCLOUD))) = ZXTTE(KIDIA:KFDIA,1:KLEV,ind_oifs_ham%ind_cloud_HAM(JCLOUD))
    END DO
    !<-- End adding HAM modified tendency back to PTENC
    !-----------------------------------------------------------------
    
  CASE DEFAULT
    ! this case should never occur, as it is handled in the calling subroutine
    CALL ABOR1(" AEROSOL SCHEME "//TRIM(AERO_SCHEME)//" IS NOT HANDLED IN HAMM7" )

END SELECT

! write flux to extra fields for diagnostic of aerosol 'chemical' conversion  
IF (LCHEM_DIA) THEN
  CALL CHEM_INEXT( KIDIA , KFDIA  , KLON , KLEV , 1, 1,  &
     &    ZDP, PTSPHY, ZTSO4, ZTENC0,PEXTRA(:,ISO4_C,IEXTR_CH))
END IF


!*         4.    ELIMINATION OF NEGATIVE PROGNOSTIC AEROSOL CONCENTRATIONS
!                 ---------------------------------------------------------

IF (LAERNGAT) THEN

  IF (LCHEM_DIA) THEN
    ZTAERO0(KIDIA:KFDIA,1:KLEV,1:NACTAERO) =  ZTAEROK(KIDIA:KFDIA,1:KLEV,1:NACTAERO)
  ENDIF

  DO JAER=1,NACTAERO

      DO JK=1,KLEV
        DO JL=KIDIA,KFDIA
          ZAER(JL,JK) = ZCEN(JL,JK,KAERO(JAER))
          ZTAER(JL,JK)= PTENC(JL,JK,KAERO(JAER))
        ENDDO
      ENDDO

    CALL AER_NEGAT &
         & ( YREAERATM, KIDIA   , KFDIA, KLON , KLEV, &
         &   PTSPHY  , &
         &   ZAER    , ZTAER, PRS1, &
         &   ZAERNEG )  

    DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
        ZTAERO(JL,JK,JAER)=ZTAER(JL,JK)
      ENDDO
    ENDDO
    DO JL=KIDIA,KFDIA
      ZAERNGT(JL,JAER)=ZAERNEG(JL,KLEV)
    ENDDO

  ENDDO

  ! collect neg fix tendencies
  IF (LCHEM_DIA) THEN
    CALL CHEM_INEXT( KIDIA , KFDIA  , KLON , KLEV , NACTAERO, NACTAERO ,  &
         &    ZDP, PTSPHY, ZTAERO ,ZTAERO0, PEXTRA(:,NCHEM+1:NCHEM+NACTAERO,IEXTR_NG))
  ENDIF


! do not fix the tendencies for now, number concentration fixes will break the 
! correlation between mass and number
  PTENC(KIDIA:KFDIA,1:KLEV,KAERO(1):KAERO(NACTAERO)) = ZTAERO(KIDIA:KFDIA,1:KLEV,1:NACTAERO)

!ELSE
!  ZAERNGT(:,:)=0._JPRB  
ENDIF

!------------------------------------------------------------------------------
!*         7.       STORE ALL AEROSOL VERTICALLY INTEGRATED FLUXES
!                   ----------------------------------------------

DO JAER=1,NACTAERO
  DO JL=KIDIA,KFDIA
    PAERODDF(JL,JAER,1)=PAERSRC(JL,JAER) !aerosol so4 source term 
    PAERODDF(JL,JAER,2)=PAERDDP(JL,JAER) ! aerosol dry deposition
    PAERODDF(JL,JAER,3)=PAERSDM(JL,JAER) ! aerosol sedimentation 
    PAERODDF(JL,JAER,4)=0.0!ZAERSCL(JL,JAER) ! so2 sink added to scavenging
    PAERODDF(JL,JAER,5)=0.0!ZAERSCC(JL,JAER) ! scavenging (in-cloud & below cloud) so wet deposition
    PAERODDF(JL,JAER,6)=ZAERNGT(JL,JAER)
    PAERODDF(JL,JAER,7)=0.0!ZAERTAUT(JL,JAER,1) !total AOD?
  ENDDO
ENDDO

!-----------------------------------------------------------------------


!*         5.      OPTICAL DEPTH
!                  -------------------------------------------------- 
!calculate optical properties only when radiation is calculated 
!radiation is calculated before microphysics -> nstep+1

INWAVL = 20
ITWAVL( 1)= 9   ! 550 nm
ITWAVL( 2)= 1   ! 340 nm
ITWAVL( 3)= 2   ! 355 nm 
ITWAVL( 4)= 3   ! 380 nm
ITWAVL( 5)= 4   ! 400 nm
ITWAVL( 6)= 5   ! 440 nm
ITWAVL( 7)= 6   ! 469 nm
ITWAVL( 8)= 7   ! 500 nm
ITWAVL( 9)= 8   ! 532 nm 
ITWAVL(10)=10   ! 645 nm
ITWAVL(11)=11   ! 670 nm
ITWAVL(12)=12   ! 800 nm
ITWAVL(13)=13   ! 858 nm                  
ITWAVL(14)=14   ! 865 nm
ITWAVL(15)=15   ! 1020 nm
ITWAVL(16)=16   ! 1064 nm
ITWAVL(17)=17   ! 1240 nm
ITWAVL(18)=18   ! 1640 nm
ITWAVL(19)=19   ! 2130 nm
ITWAVL(20)=20   ! 10 microns

IBLK=(KSTGLO-1)/KLON + 1

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    DO JAER=1,14
      PAER_TAU(JL,JK,JAER)   = YDAERM7%M7AOD(JL,JK,JAER,IBLK)
      PAER_SSA(JL,JK,JAER)   = YDAERM7%M7SSA(JL,JK,JAER,IBLK)
      PAER_ASYM(JL,JK,JAER ) = YDAERM7%M7ASYM(JL,JK,JAER,IBLK)
    ENDDO
    DO JAER=1,16    
      PAER_TAU_LW(JL,JK,JAER)= YDAERM7%M7AODLW(JL,JK,JAER,IBLK)
    ENDDO
  ENDDO
ENDDO

IF(MOD(NSTEP+1,NRADFR) == 0) THEN
CALL GSTATS(2506,0)
ZAER_TAU(KIDIA:KFDIA,:,:,:)  = 0.0_JPRB
ZAER_SSA(KIDIA:KFDIA,:,:)    = 0.0_JPRB
ZAER_ASYM(KIDIA:KFDIA,:,:)   = 0.0_JPRB
ZAER_TAU_LW(KIDIA:KFDIA,:,:) = 0.0_JPRB

SELECT CASE (NAEROOPT) 

CASE (0)

   ! Optical properties is not calculated

CASE (1)
   ! Use TM5 codes to calculate optical properties (optical properties for LW = 0)
      !--> Add HAM updated tendency to ZTAERO and use that in optics
   ZTAERO(:,:,:) = 0._JPRB
   DO JAER=1,NACTAERO
      DO JK=1,KLEV
         DO JL=KIDIA,KFDIA
            ZTAERO(JL,JK,JAER)= PTENC(JL,JK,KAERO(JAER))
         ENDDO
      ENDDO
   ENDDO

   ZAEROK = ZAEROK+ZTAERO*time_step_len

   CALL TM5M7_OPTICS_AOP_GET( YGFL, YREAERSRC, KIDIA,KFDIA, KLON, KLEV,NACTAERO, &
   &                          NASWBAND, ASWBAND, 1, .false., &
   &                          ZRHO, ZAEROK,RW_MODE,RWD_MODE,H2O_MODE,&
   &                          ZAER_TAU, ZAER_SSA, ZAER_ASYM)

   PAOD(KIDIA:KFDIA,:)=0._JPRB
   DO JK = 1, KLEV
     DO JL = KIDIA,KFDIA
       DO IW=1,NASWBAND         
         PAER_TAU(JL,JK,IW)=ZAER_TAU(JL,JK,IW,1)*(PGEOH(JL,JK-1) - PGEOH(JL,JK))
         PAER_SSA(JL,JK,IW)=ZAER_SSA(JL,JK,IW)
         PAER_ASYM(JL,JK,IW)=ZAER_ASYM(JL,JK,IW)
         !PAOD(JL,IW)=ZAER_TAU(JL,JK,IW,1)*(PGEOH(JL,JK-1) - PGEOH(JL,JK))+PAOD(JL,IW)
       ENDDO
       DO IW=1,16
         PAER_TAU_LW(JL,JK,IW)=0.0_JPRB
       END DO
     ENDDO
   ENDDO

   ! "I am not sure if the rest is needed" [who is I? FIXME]
   
   ALLOCATE( ZAOP_OUT_EXT( KLON, KLEV, NWDEP, 1)) ; ZAOP_OUT_EXT = 0.0_JPRB
   ALLOCATE( ZAOP_OUT_A  ( KLON, KLEV, NWDEP)   ) ; ZAOP_OUT_A   = 0.0_JPRB
   ALLOCATE( ZAOP_OUT_G  ( KLON, KLEV, NWDEP)   ) ; ZAOP_OUT_G   = 0.0_JPRB


   ALLOCATE(ZTAUS_AER (KLON, KLEV,NWDEP)); ZTAUS_AER = 0.0
   ALLOCATE(ZTAUA_AER (KLON, KLEV,NWDEP)); ZTAUA_AER = 0.0  
   ALLOCATE(ZPMAER    (KLON, KLEV,NWDEP)); ZPMAER    = 0.0

   DO JB=1, NBANDS_TROP
     PTAUS_AER(:,:,JB,1) = ZTAUS_AER(:,:,WAV_GRID(JB))
     PTAUA_AER(:,:,JB,1) = ZTAUA_AER(:,:,WAV_GRID(JB))
     PPMAER   (:,:,JB,1) = ZPMAER   (:,:,WAV_GRID(JB))

     PTAUS_AER(:,:,JB,2) = ZTAUS_AER(:,:,WAV_GRIDA(JB))
     PTAUA_AER(:,:,JB,2) = ZTAUA_AER(:,:,WAV_GRIDA(JB))
     PPMAER   (:,:,JB,2) = ZPMAER   (:,:,WAV_GRIDA(JB))
   ENDDO
   DEALLOCATE(ZTAUS_AER)
   DEALLOCATE(ZTAUA_AER)
   DEALLOCATE(ZPMAER)

   DEALLOCATE(ZAOP_OUT_EXT)
   DEALLOCATE(ZAOP_OUT_A  )
   DEALLOCATE(ZAOP_OUT_G  )

 CASE (2)
   ! Use HAM codes to calculate optical properties
   LWBANDS=16
   PRS1D(KIDIA:KFDIA,:) = PRS1(KIDIA:KFDIA,1:KLEV)-PRS1(KIDIA:KFDIA,0:KLEV-1)
   !CALL ham_rad_cache(KLON,KLEV)
   ZXTM0(KIDIA:KFDIA,1:KLEV,:) = ZXTM1(KIDIA:KFDIA,1:KLEV,:) + ZXTTE(KIDIA:KFDIA,1:KLEV,:)*time_step_len
   CALL HAM_RAD(KFDIA, KLON, KLEV, ZKROW, LWBANDS, NASWBAND, ZXTM0, PRS1D, &
        & ZAER_TAU(:,:,:,1), ZAER_SSA, ZAER_ASYM, ZAER_TAU_LW, ZM6RP)
   !CALL ham_rad_cache_cleanup

   DO JK = 1, KLEV
     DO JL = KIDIA,KFDIA
       DO IW=1,NASWBAND         
         PAER_TAU(JL,JK,IW)=ZAER_TAU(JL,JK,IW,1)!*(PGEOH(JL,JK-1) - PGEOH(JL,JK))
         PAER_SSA(JL,JK,IW)=ZAER_SSA(JL,JK,IW)
         PAER_ASYM(JL,JK,IW)=ZAER_ASYM(JL,JK,IW)
         !PAOD(JL,IW)=ZAER_TAU(JL,JK,IW,1)+PAOD(JL,IW)
       ENDDO
       DO IW=1,16
         PAER_TAU_LW(JL,JK,IW)=ZAER_TAU_LW(JL,JK,IW)
       END DO
     ENDDO
   ENDDO

 END SELECT


CALL GSTATS(2506,1)
ENDIF  ! (MOD(NSTEP+1,NRADFR) == 0)

PAOD(:,:) =0._JPRB
PABS(:,:) =0._JPRB
PFAOD(:,:)=0._JPRB
PSSA(:,:) =0._JPRB
PASY(:,:) =0._JPRB
DO JK = 1, KLEV
  DO JL = KIDIA,KFDIA
    DO IW=1,NASWBAND  
      PAOD(JL,IW)=PAER_TAU(JL,JK,IW)+PAOD(JL,IW)         
      PSSA(JL,IW)=PAER_SSA(JL,JK,IW)*PAER_TAU(JL,JK,IW)+PSSA(JL,IW)         
      PASY(JL,IW)=PAER_ASYM(JL,JK,IW)*PAER_TAU(JL,JK,IW)+PASY(JL,IW)         
    END DO
    !IF(PAOD(JL,IW)>0._JPRB) THEN
    !    PSSA(JL,:) =PSSA(JL,:)/PAOD(JL,IW)! AOD AVERAGE     
    !    PASY(JL,:) =PASY(JL,:)/PAOD(JL,IW)! AOD AVERAGE
    !ENDIF                     
  END DO
END DO

DO JL = KIDIA,KFDIA
  DO IW=1,NASWBAND
    IF(PAOD(JL,IW)>0._JPRB) THEN
      PSSA(JL,IW) =PSSA(JL,IW)/PAOD(JL,IW)! AOD AVERAGE     
      PASY(JL,IW) =PASY(JL,IW)/PAOD(JL,IW)! AOD AVERAGE
    ENDIF
  END DO
END DO

!------------------------------------------------------------------------------
!*         6.0     Fill selective aerosol OD fields in structure as available in IFS-AER

DO JWAVL=1,MIN(INWAVL,NAERO_WVL_DIAG)
  DO JL=KIDIA,KFDIA
    IF (YDMODEL%YRML_GCONF%YGFL%NAERO_WVL_DIAG_TYPES >= JPAERO_WVL_AOD) THEN
      PAERO_WVL_DIAG(JL,JWAVL,JPAERO_WVL_AOD) = PAOD(JL,JWAVL)
    ENDIF
    IF (YDMODEL%YRML_GCONF%YGFL%NAERO_WVL_DIAG_TYPES >= JPAERO_WVL_AODABS) THEN
      PAERO_WVL_DIAG(JL,JWAVL,JPAERO_WVL_AODABS) = PABS(JL,JWAVL)! 0.0
    ENDIF
    IF (YDMODEL%YRML_GCONF%YGFL%NAERO_WVL_DIAG_TYPES >= JPAERO_WVL_AODFM) THEN
      PAERO_WVL_DIAG(JL,JWAVL,JPAERO_WVL_AODFM) = PFAOD(JL,JWAVL)! 0.0
    ENDIF
    IF (YDMODEL%YRML_GCONF%YGFL%NAERO_WVL_DIAG_TYPES >= JPAERO_WVL_SSA) THEN
      PAERO_WVL_DIAG(JL,JWAVL,JPAERO_WVL_SSA) = PSSA(JL,JWAVL)
    ENDIF
    IF (YDMODEL%YRML_GCONF%YGFL%NAERO_WVL_DIAG_TYPES >= JPAERO_WVL_ASSIMETRY) THEN
      PAERO_WVL_DIAG(JL,JWAVL,JPAERO_WVL_ASSIMETRY) = PASY(JL,JWAVL)
    ENDIF
  ENDDO
ENDDO

!*
!* VH - Requires full checking by expert if this makes sense - FIXME
!*
DO JL=KIDIA,KFDIA
  PODTO(JL)    =PTAUS_AER(JL,KLEV,1,1)
  PODTO469(JL) =PTAUS_AER(JL,KLEV,2,1)
  PODTO670(JL) =PTAUS_AER(JL,KLEV,3,1)
  PODTO865(JL) =PTAUS_AER(JL,KLEV,4,1)
  PODTO1240(JL)=PTAUS_AER(JL,KLEV,5,1)
ENDDO

!*         6.1     STORE IN AEROUT-1 
!                  ------------------------------

!-- total instantaneous optical depth
!  DO JB=1,NBANDS_TROP
!    DO JL=KIDIA,KFDIA
!      ZAEROUT1(JL,JB)= PTAUS_AER(JL,KLEV,JB,1)
!      ZAEROUT1(JL,NBANDS_TROP+JB)= PTAUA_AER(JL,KLEV,JB,1)
!    ENDDO
!  ENDDO 

!-- 

!*         6.1     STORE IN AEROUT2-AEROUT4
!                  ------------------------------

!-- the total extinction coefficient at wavelengths ?? nm is archived in GFL%AEROUT
!    DO JK=1,KLEV
!      DO JL=KIDIA,KFDIA
!        ZAEROUT2(JL,JK)=PTAUS_AER(JL,JK, 1,1)
!        ZAEROUT3(JL,JK)=PTAUA_AER(JL,JK, 1,1)
!        ZAEROUT4(JL,JK)=Ppmaer(JL,JK,1,1)
!      ENDDO
!    ENDDO


! RCHG -> TODO explain here LIFSMIN and LIFSTRAJ ---
IF(.NOT.LIFSMIN  .AND. .NOT.LIFSTRAJ) THEN
  ! input for HAM-M7
  PGFL(KIDIA:KFDIA,1,YAEROUT(1)%MP)=PAOD(KIDIA:KFDIA,10) !PAER_TAU(KIDIA:KFDIA,1:KLEV,10) !533nm
  DO JK=1,KLEV 
    PGFL(KIDIA:KFDIA,2,YAEROUT(1)%MP)= PGFL(KIDIA:KFDIA,2,YAEROUT(1)%MP) + PAER_SSA(KIDIA:KFDIA,JK,10)   !PAER_TAU(KIDIA:KFDIA,1:KLEV,10) !533nm
    PGFL(KIDIA:KFDIA,3,YAEROUT(1)%MP)= PGFL(KIDIA:KFDIA,3,YAEROUT(1)%MP) + PAER_ASYM(KIDIA:KFDIA,JK,10)  !PAER_TAU(KIDIA:KFDIA,1:KLEV,10) !533nm
  END DO

  ! RCHG: It seems that here 1:14 refers to 14 species. PLS: AOD of 14 short wavelengths stored in YAEROUT(12)%MP of first 14 tracers (juggling!!) 
  PGFL(KIDIA:KFDIA,1:14,YAEROUT(12)%MP)=PAOD(KIDIA:KFDIA,1:14)

  DO JN=1,NAEROCOMP !ntrac!NACTAERO
    PGFL(KIDIA:KFDIA,ind_oifs_ham%ind_mass_OIFS(JN),YAEROUT(2)%MP)=ZDDEPFLUX(KIDIA:KFDIA,ind_oifs_ham%IND_mass_HAM(JN))
  END DO
  DO JN=1,NCLASS
    PGFL(KIDIA:KFDIA,ind_oifs_ham%ind_class_OIFS(JN),YAEROUT(2)%MP)=ZDDEPFLUX(KIDIA:KFDIA,ind_oifs_ham%IND_class_HAM(JN))
  END DO
  !do JN=1,NACTAERO   !ktrac
  ! if (NSTEP==10)then
  !    do JL=KIDIA,KFDIA
  !       write(3345,*)JN,WDEPOUT_2D(JL,KAERO(JN))
  !    end do
  ! end if
  !end do

  DO JN=1,NACTAERO   !ktrac
    PGFL(KIDIA:KFDIA,KAERO(JN),YAEROUT(3)%MP)=WDEPOUT_2D(KIDIA:KFDIA,KAERO(JN))
    PGFL(KIDIA:KFDIA,KAERO(JN),YAEROUT(4)%MP)=SEDOUT_2D(KIDIA:KFDIA,KAERO(JN))
    ! PGFL(KIDIA:KFDIA,JN,YAEROUT(5)%MP)=PAERSRC(KIDIA:KFDIA,KAERO(JN))
  END DO

  PGFL(KIDIA:KFDIA,1,YAEROUT(6)%MP)=PAERFLX(KIDIA:KFDIA,3,9)
  PGFL(KIDIA:KFDIA,:,YAEROUT(7)%MP)=0.0_JPRB
  ! 1=TOtal SO4 not in use
  ! 2 Ammonium not in use
  ! 3 Nitrate not in use
  DO JK=1,KLEV
    ! save load for each N/M as one level
    ! kg/kg -> kg/m2 N/kg-> N/m2
    DO JN=1,NAEROCOMP    !ntrac!NACTAERO   
      JO=ind_oifs_ham%ind_mass_OIFS(JN)  ! JO -> index context OIFS 
      JH=ind_oifs_ham%IND_mass_HAM(JN)   ! JH -> index context HAM 
      JY=YAEROUT(7)%MP
      PGFL(KIDIA:KFDIA,JO,JY) = PGFL(KIDIA:KFDIA, JO, JY) + (ZXTM1(KIDIA:KFDIA,JK,JH)+(ZXTTE(KIDIA:KFDIA,JK,JH)*time_step_len)) * ZDPG(KIDIA:KFDIA,JK)
    END DO
    DO JN=1,NCLASS
      JO=ind_oifs_ham%ind_mass_OIFS(JN)  ! JO -> index context OIFS 
      JH=ind_oifs_ham%IND_mass_HAM(JN)   ! JH -> index context HAM 
      JY=YAEROUT(7)%MP
      PGFL(KIDIA:KFDIA,JO,JY) = PGFL(KIDIA:KFDIA, JO, JY) + (ZXTM1(KIDIA:KFDIA,JK,JH)+(ZXTTE(KIDIA:KFDIA,JK,JH)*time_step_len)) * ZDPG(KIDIA:KFDIA,JK)
    END DO
  END DO

  IF (.NOT. LAERCHEM) THEN
    PGFL(KIDIA:KFDIA,1,YAEROUT(8)%MP)       = ZDDEPFLUX_SO2(KIDIA:KFDIA) ! aergn7 SO4 gas
    PGFL(KIDIA:KFDIA,1:KLEV,YAEROUT(9)%MP)  = PGFL(KIDIA:KFDIA,1:KLEV,YAEROCLIM(1)%MP)
    PGFL(KIDIA:KFDIA,1:KLEV,YAEROUT(10)%MP) = PGFL(KIDIA:KFDIA,1:KLEV,YAEROCLIM(2)%MP)
    PGFL(KIDIA:KFDIA,1:KLEV,YAEROUT(11)%MP) = PGFL(KIDIA:KFDIA,1:KLEV,YAEROCLIM(3)%MP)
  END IF
  PGFL(KIDIA:KFDIA,1:KLEV,YAEROUT(12)%MP) = PGFL(KIDIA:KFDIA,1:KLEV,YRE_LIQ%MP9_PH)
  PGFL(KIDIA:KFDIA,1:KLEV,YAEROUT(13)%MP) = PGFL(KIDIA:KFDIA,1:KLEV,YRE_ICE%MP9_PH)
  PGFL(KIDIA:KFDIA,1:KLEV,YAEROUT(14)%MP) = PGFL(KIDIA:KFDIA,1:KLEV,YCDNC%MP9_PH)
  PGFL(KIDIA:KFDIA,1:KLEV,YAEROUT(15)%MP) = PGFL(KIDIA:KFDIA,1:KLEV,YICNC%MP9_PH)
  PGFL(KIDIA:KFDIA,KLEV,YAEROUT(16)%MP)   = ZDUMMY(KIDIA:KFDIA,17)         ! SS CS wdepflux
  PGFL(KIDIA:KFDIA,KLEV-1,YAEROUT(16)%MP) = ZDUMMY(KIDIA:KFDIA,7)          ! SO4 CS wdepflux
  PGFL(KIDIA:KFDIA,KLEV-2,YAEROUT(16)%MP) = ZDUMMY(KIDIA:KFDIA,25)         ! NUM CS wdepflux

  PGFL(KIDIA:KFDIA,KLEV,YAEROUT(17)%MP)   = ZDDEPFLUX(KIDIA:KFDIA,25)      ! drydep flux NUM CS ham
  PGFL(KIDIA:KFDIA,KLEV-1,YAEROUT(17)%MP) = ZDDEPFLUX(KIDIA:KFDIA,17)      ! drydep flux SS CS ham
  PGFL(KIDIA:KFDIA,KLEV-2,YAEROUT(17)%MP) = ZDDEPFLUX(KIDIA:KFDIA,7)       ! drydep flux SO4 CS ham
  PGFL(KIDIA:KFDIA,KLEV,YAEROUT(18)%MP)   = ZVDEP(KIDIA:KFDIA,25)          ! ddepveloc NUM CS ham
  PGFL(KIDIA:KFDIA,KLEV-1,YAEROUT(18)%MP) = ZVDEP(KIDIA:KFDIA,17)          ! ddepveloc SS CS ham
  PGFL(KIDIA:KFDIA,KLEV-2,YAEROUT(18)%MP) = ZVDEP(KIDIA:KFDIA,7)           ! ddepveloc SO4 CS ham
  PGFL(KIDIA:KFDIA,KLEV,YAEROUT(19)%MP)   = ZXTTE(KIDIA:KFDIA,KLEV,3)      ! tendency SS CS ham after update surface
  PGFL(KIDIA:KFDIA,KLEV-1,YAEROUT(20)%MP) = ZTENCIH(KIDIA:KFDIA,KLEV,17)   ! tendency SS CS ham before update surface
  PGFL(KIDIA:KFDIA,1:KLEV,YAEROUT(21)%MP) = ZXTMD1(KIDIA:KFDIA,1:KLEV,17)  ! mix rat SS CS ham before update


  ! It is not clear when NACTERO and when NTRAC 
  DO JN=1,NACTAERO
    PGFL(KIDIA:KFDIA,JN,YAEROUT(22)%MP)=DDEPOUT(KIDIA:KFDIA,KLEV,KAERO(JN))
    PGFL(KIDIA:KFDIA, JN, YAEROUT(28)%MP) = PAERSRC(KIDIA:KFDIA,KAERO(JN)) ! Emissions per specie
  END DO
  DO JN=1,NTRAC
    PGFL(KIDIA:KFDIA,JN,YAEROUT(23)%MP)=WDEPOUT(KIDIA:KFDIA,KLEV,JN)
    PGFL(KIDIA:KFDIA,JN,YAEROUT(24)%MP)=SEDOUT(KIDIA:KFDIA,KLEV,JN)
    PGFL(KIDIA:KFDIA,JN,YAEROUT(39)%MP)=ZXTEMS(KIDIA:KFDIA,JN)
  END DO
  DO JN=1,SUBM_NGASSPEC
    PGFL(KIDIA:KFDIA,JN,YAEROUT(25)%MP)=zxtm1(KIDIA:KFDIA,KLEV,ind_oifs_ham%ind_gas_HAM(JN))
  END DO

  PGFL(KIDIA:KFDIA,1:KLEV,YAEROUT(26)%MP)=ZSEDIFLUX(KIDIA:KFDIA,1:KLEV,25)

  IF (LAERCHEM)THEN
    DO JGAS=1,SUBM_NGASSPEC
      !ZXTM1(JL,JK,ind_gas_HAM(JGAS)) = MAX(0._JPRB,ZCEN(JL,JK,KCHEM(ind_gas_OIFS(JGAS)))) !eehol: remove negative values
      !PGFL(KIDIA:KFDIA,JN,YAEROUT(28+JGAS)%MP)=ZCEN(JL,JK,KCHEM(ind_gas_OIFS(JGAS)))
      PGFL(KIDIA:KFDIA,1:KLEV,YAEROUT(28+JGAS)%MP)= PTENC(KIDIA:KFDIA,1:KLEV,KCHEM(ind_oifs_ham%ind_gas_OIFS(JGAS)))
    END DO
    !ELSE
    !   DO JGAS=1,subm_ngasspec
    !       PGFL(KIDIA:KFDIA,1:KLEV,YAEROUT(28+JGAS)%MP)= PTENC(KIDIA:KFDIA,1:KLEV,KAERO(ind_oifs_ham%ind_gas_OIFS(JGAS)))
    !   ENDDO  
  ENDIF

  DO IMODE=1,NMOD
    PGFL(KIDIA:KFDIA,1:KLEV,YAEROUT(30+IMODE)%MP)=RW_MODE(IMODE)%d2(KIDIA:KFDIA,1:KLEV) ! m
  ENDDO

  IF (.NOT.LAERCHEM) THEN
    PGFL(KIDIA:KFDIA,KLEV,YAEROUT(40)%MP)   = ZFSO2(KIDIA:KFDIA)           ! tendency SS CS ham after update surface
    PGFL(KIDIA:KFDIA,KLEV,YAEROUT(41)%MP)   = ZFSO4(KIDIA:KFDIA)           ! tendency SS CS ham after update surfac
    PGFL(KIDIA:KFDIA,KLEV,YAEROUT(42)%MP)   = ZFSO4_AQ(KIDIA:KFDIA)        ! tendency SS CS ham after update surface
    PGFL(KIDIA:KFDIA,1:KLEV,YAEROUT(43)%MP) = ZTSO4(KIDIA:KFDIA,1:KLEV,1)  ! tendency SS CS ham after update surface
    PGFL(KIDIA:KFDIA,1:KLEV,YAEROUT(44)%MP) = ZTSO4_AQ(KIDIA:KFDIA,1:KLEV) ! tendency SS CS ham after update surface
    PGFL(KIDIA:KFDIA,1:KLEV,YAEROUT(45)%MP) = ZTSO2(KIDIA:KFDIA,1:KLEV)    ! tendency SS CS ham after update surface
  END IF

ENDIF

!------------------------------------------------------------------------------
!*         9.       RELEASE LOCAL MEMORY
!                   --------------------
!DEALLOCATE( ZAERSRC )
IF (ALLOCATED(ZAERNGT ) ) DEALLOCATE( ZAERNGT )
!DEALLOCATE( ZAERSCC )

IF (ALLOCATED(ZWPDF ) ) DEALLOCATE( ZWPDF )
IF (ALLOCATED(ZW ) )    DEALLOCATE( ZW )
IF (ALLOCATED(ZRC) )    DEALLOCATE( ZRC )
IF (ALLOCATED(ZSMAX))   DEALLOCATE( ZSMAX )

DO IMODE=1,NMOD
  IF (ASSOCIATED(RW_MODE(IMODE)%D2)) DEALLOCATE(RW_MODE(IMODE)%D2)
  IF (ASSOCIATED(DENS_MODE(IMODE)%D2)) DEALLOCATE(DENS_MODE(IMODE)%D2)
ENDDO

DO IMODE=1,NMOD
  IF (SIZECLASS(IMODE)%LSOLUBLE) THEN
    IF (ASSOCIATED(RWD_MODE(IMODE)%D2)) DEALLOCATE(RWD_MODE(IMODE)%D2)
    IF (ASSOCIATED(H2O_MODE(IMODE)%D2)) DEALLOCATE(H2O_MODE(IMODE)%D2)
  END IF
ENDDO

END ASSOCIATE
END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('HAMM7_INTERFACE',1,ZHOOK_HANDLE)

END SUBROUTINE HAMM7_INTERFACE
