MODULE TM5M7_DATA

  ! *TM5M7_DATA*  Contains parameters, switches and initialization
  !               routines for the tm5m7 aerosol scheme.
  !
  ! Authors:
  ! --------
  ! E. Vignati (JRC/IES)						   2005
  ! P. Stier (MPI)                                 2001/2002 
  ! E. Vignati and J. Wilson (JRC / EI)            2000
  ! V. Huijnen (KNMI): copied from TM5 to IFS code 2020




USE PARKIND1, ONLY : JPIM, JPRB

IMPLICIT NONE

SAVE
  
! Location of lookup table data. Default location is at the directory of the executable
!CHARACTER(LEN=256) :: TM5M7_DATADIR='/perm/ms/nl/nk9/oifs43r3/data/cifs_input/tm5m7/optics/'
CHARACTER(LEN=256) :: TM5M7_DATADIR='./'
  
TYPE MODAL_DATA
  REAL(KIND=JPRB), DIMENSION(:,:), POINTER :: d2   ! KLON, KLEV
  REAL(KIND=JPRB), DIMENSION(:)  , POINTER :: surf ! KLON
END TYPE MODAL_DATA

  !--- 1) Define and pre-set switches for the processes of M7: -----------------------

  !--- Physical:

  LOGICAL            :: lsnucl     = .TRUE., & ! nucleation
     &                  lscoag     = .TRUE., & ! coagulation
     &                  lscond     = .TRUE.    ! condensation of H2SO4
     
  INTEGER(KIND=JPIM) :: nnucl      = 1         ! Choice of the nucleation scheme:
                                               !    nnucl = 1  Vehkamaeki (2002)
                                               !          = 2  Kulmala (1998) -NOT RECOMMENDED-
  !--- Technical:

  LOGICAL :: lmass_diag = .FALSE.   ! mass balance check in m7_interface

  !--- 2) Numbers of compounds and modes of m7: --------------------------------------

  INTEGER(KIND=JPIM), PARAMETER :: naermod=23,     & !number of all compounds
     &                             nmod=7,         & !number of modes
     &                             nss=2,          & !number of sea salt compounds 
     &                             nsol=4,         & !number of soluble  compounds 
     &                             ngas=3,         & !number of gaseous  compounds 
     &                             nsulf=4,        & !number of sulfate  compounds 
     &                             ncomp=5           !number of compounds

  !--- 3) List of indexes corresponding to the compound masses and mode numbers:------

  !--- 3.1) Mass index (in array aerml and ttn): 
  !
  !         Attention:
  !         The mass of sulfate compounds is always given in [molec. cm-3] 
  !         whilst the mass of other compounds is given in [ug cm-3].
  !
  !         Compounds:
  !
  !           so4 = sulphate
  !           bc  = black carbon
  !           oc  = organic carbon, 
  !           ss  = sea salt
  !           du  = dust 
  !
  !         Modes:
  !
  !           n   = nucleation mode
  !           k   = aitken mode 
  !           a   = accumulation mode
  !           c   = coarse mode
  !
  !         Type:
  !
  !           s   = soluble mode
  !           i   = insoluble mode
                                                                                  !  COMPOUND:
  INTEGER(KIND=JPIM), PARAMETER ::                                               &
     &      iso4ns=1, iso4ks=2, iso4as=3, iso4cs=4,                              & !- Sulfate
     &                ibcks =5, ibcas =6, ibccs =7, ibcki =8,                    & !- Black Carbon
     &                iocks =9, iocas=10, ioccs=11, iocki=12,                    & !- Organic Carbon
     &                          issas=13, isscs=14,                              & !- Sea Salt
     &                          iduas=15, iducs=16,           iduai=17, iduci=18,& !- Dust  
     &      isoans=19,isoaks=20,isoaas=21,isoacs=22,isoaki=23                      !- SOA
  ! MODE:            |         |         |         |         |
  !         nucl .   | aitk.   | acc.    | coar.   | aitk.   | acc.    | coar.   |
  !         solub le | soluble | soluble | soluble | insol.  | insol.  | insol.  |

   
  !--- 3.2) Number index (in array aernl):
  !

  INTEGER(KIND=JPIM), PARAMETER ::                                  &
     &      inucs=1,  iaits=2,  iaccs=3,  icoas=4,  iaiti=5,  iacci=6,  icoai=7    
  ! MODE:           |         |         |         |         |
  !         nucl.   | aitk.   | acc.    | coar.   | aitk.   | acc.    | coar.   |
  !         soluble | soluble | soluble | soluble | insol.  | insol.  | insol.  |


  !--- 4) Definition of the modes of M7: ------------------------------------------------------

  !--- 4.1) Threshold radii between the different modes [cm]:
  !         Used for the repartititioning in m7_dconc.
  !         crdiv(jmod) is the lower bound and crdiv(jmod+1) is 
  !         the upper bound of the respective geometric mode
  !         Default value for nucleation mode is modified by the 
  !         choice of the nuclation scheme.

  REAL(KIND=JPRB) :: crdiv(4)=(/ 0.0005E-4, 0.005E-4, 0.05E-4, 0.5E-4 /)    
  !                                       |         |        |      
  !                                       |         |        |
  !                           nucleation -- aitken  -  accum -- coarse mode

  !--- 4.2) Standard deviation for the modes:

  REAL(KIND=JPRB), PARAMETER :: sigma(nmod)=(/ 1.59, 1.59, 1.59, 2.00, 1.59, 1.59, 2.00 /)

  !--- Natural logarithm of the standard deviation of each mode:
  !    Calulated in m7_initialize. 

  REAL(KIND=JPRB)            :: sigmaln(nmod)

  !--- 5) Conversion factors for lognormal particle size distributions: -------------
  !       Calulated in m7_initialize. 

  REAL(KIND=JPRB)            :: cmr2ras(nmod) ! Conversion factor: count median radius to radius of average surface

  REAL(KIND=JPRB)            :: cmr2mmr(nmod) ! Conversion factor: count median radius to mass mean radius

  REAL(KIND=JPRB)            :: cmedr2mmedr(nmod) ! Conversion factor: count median radius to mass median radius

  REAL(KIND=JPRB)            :: cmr2ram(nmod) ! Conversion factor: count median radius to radius of average mass

  REAL(KIND=JPRB)            :: ram2cmr(nmod) ! Conversion factor: radius of average mass to count median radius


  !--- 6) Assumed thresholds for occurence of specific quantities: -------------
  !@@@    To be done!
  
  !  REAL, PARAMETER :: cmin_aerml     = 1.E-15 , ! Aerosol mass
  !                     cmin_aernl     = 1.E-10 , ! Aerosol number
  !                     
  
  !--- 7) Chemical constants: ----------------------------------------------------
  !
  !--- Accomodation coefficient of H2SO4 on aerosols:
  !    (reduced for insoluble modes)

  REAL(KIND=JPRB), PARAMETER :: caccso4(nmod) = (/ 1.0, 1.0, 1.0, 1.0, 0.3, 0.3, 0.3 /)

  !--- Critical relative humidity:

  REAL(KIND=JPRB), PARAMETER :: crh    = 0.45 ! Assumed relative humidity for the 
                                              ! Na2SO4 / NaCl system below which 
                                              ! crystalization occurs.
                                              ! (estimated from Tang, I.N.; JGR 102, D2 1883-1893)

  !--- 8) Physical constants: ----------------------------------------------------
  !
  !--- 8.1) General physical constants: 

  REAL(KIND=JPRB), PARAMETER :: &
                  &   bk      = 1.38e-16,   & ! Bolzman constant []
                  &   avo     = 6.02217E+23,& ! Avogadro number [mol-1]
                  &   rerg    = 8.314E+7,   & ! Ideal gas constant [erg.K-1.mole-1]
                  &   r_kcal  = 1.986E-3      ! Ideal gas constant [kcal K-1.mole-1]
  
  !--- 8.2) Type specific physical constants:
  !
  REAL(KIND=JPRB), PARAMETER :: &
                 &   dh2so4  = 1.841,      & ! Density          H2SO4  [g cm-3]
                 &   ddust   = 2.650,      & ! Density          du     [g cm-3]
                     !>>> TvN 
                     ! The density of BC is in the range 1.7 to 1.9 g/cm3. 
                     ! (Bond and Bergstrom, Aerosol Sci. Technol., 2006).
                     ! We therefore adopt a value of 1.8 g/cm3.
                     ! Details can be found in Bond et al. (JGR, 2013),
                     ! and references therein:
                     ! Park et al. (J. Nanoparticle Research, 2004) measured
                     ! 1.77 +- 0.07 g/cm3 for the non-volatile components of diesel soot, 
                     ! and give a range 1.7-1.8 g/cm3 in their conclusions.
                     ! Kondo et al. (Aerosol Sci. Techn., 2011) measured 
                     ! 1.718 +- 0.004 g/cm3 for fullerene soot.
                     ! Schmid et al. (Environ. Sci. Technol., 2009)
                     ! derive a value 1.8 +- 0.2 g/cm3
                     ! for elemental carbon from biomass burning.
                     ! For comparison, in GLOMAP a value of 1.5 g/cm3 is used.
                     ! Note that these density estimates measure the mass per volume 
                     ! occupied by the spherules, as should be the case:
                     ! "If the radiative forcing of BC particles is to be
                     ! calculated from their mass concentrations, 
                     ! as it is usually the case, density should represent
                     ! the material density of the spherules, and not that
                     ! of their ramiform (branched) or aciniform (packed) aggregates."
                     ! (A. Gelencser, Carbonaceous Aerosol, Springer, 2004, p. 228).
                     ! This is explained in more detail by Bond and Bergstrom 
                     ! (Aerosol Sci. Technol., 2005).
                     ! 
                     !dbc     = 2.,         & ! Density          bc     [g cm-3]
                 &   dbc     = 1.8,        & ! Density          bc     [g cm-3]
                     ! The density of OA is highly variable,
                     ! but in any case substantially lower than 2 g/cm3.
                     ! We adopt an average value of 1.3 g/cm3,
                     ! based on a number of studies:
                     ! Turpin and Lim (Aerosol Sci. Technol., 2001)
                     ! suggest that 1.2 g/cm3 is a reasonable estimate.
                     ! Lee et al. (ACP, 2010) measured 1.26 +- 0.24 g/cm3 
                     ! at during FAME-2009 at Finokalia after evaporation
                     ! Cross et al. (Aerosol Sci. Technol., 2007) assume an
                     ! average bulk density of 1.27 g/cm3
                     ! Nakao et al. (Atmos. Environ., 2013) measure average densities 
                     ! for SOA between 1.22 and 1.42 g/cm3, depending of species.
                     ! This is in line with predictions from the OA density model
                     ! by Kuwata et al. (Environ. Sci. Technol., 2012),
                     ! who find a range between 1.23 and 1.46 g/cm3 for SOA.
                     ! This model can be used to estimate OA density 
                     ! as function of O:C and H:C elemental ratios,
                     ! with an accuracy of 12% or more.
                     ! As a further simplication, it is often assumed that
                     ! H:C = 2 - O:C (e.g. Murphy et al., ACP, 2011).
                     ! The model, however, is restricted to particle components
                     ! having negligible quantities of additional elements, 
                     ! most notably nitrogen.
                     ! Schmid et al. (Environ. Sci. Technol., 2009)
                     ! derive a value of 1.39 +- 0.13 for OA from biomass burning.
                     !
                     !doc     = 2.,         & ! Density          oc     [g cm-3]
                 &   doc     = 1.3,        & ! Density          POM     [g cm-3]
                     !<<< TvN
                 &   dnacl   = 2.165,      & ! Density          NaCl   [g cm-3]
                 &   dna2so4 = 2.68,       & ! Density          Na2SO4 [g cm-3]
                 &   dnahso4 = 2.435,      & ! Density          NaHSO4 [g cm-3]
                 &   dh2o    = 1.0,        & ! Density          H2O    [g cm-3]

                 &   wh2so4  = 98.0734,    & ! Molecular weight H2SO4  [g mol-1]
                 &   wh2o    = 18.0,       & ! Molecular weight H2O    [g mol-1]
                 &   wso4    = 96.0576,    & ! Molecular weight SO4    [g mol-1]
                 &   wso2    = 64.0,       & ! Molecular weight SO2    [g mol-1]
                 &   wna     = 22.99,      & ! Atomic    weight Na     [g mol-1]
                 &   wcl     = 35.453,     & ! Atomic    weight Cl     [g mol-1]
                 &   wnacl   = 58.443,     & ! Molecular weight NaCl   [g mol-1]
                 &   wna2so4 = 142.0376,   & ! Molecular weight Na2SO4 [g mol-1]
                 &   wnahso4 = 120.0555      ! Molecular weight NaHSO4 [g mol-1]


  !--- 9) Assumed parameters: ------------------------------------------------------

  REAL(KIND=JPRB), PARAMETER :: critn=100.,& ! Assumed mass of an nucleated sulfate 
                                             ! particle for the Kulmala scheme [molecules]
                &    fmax=0.95,            & ! Factor that limits the condensation 
                                             ! of sulfate to fmax times the available
                                             ! sulfate in the gas phase [1].
                                             ! (m7_dgas)
                &    cLayerThickness = 1.0   ! Assumed required layer thickness of
                                             ! sulfate to transfer an insoluble 
                                             ! particle to a soluble mode. It is 
                                             ! given in units of layers of 
                                             ! monomolecular sulfate. Determines the
                                             ! transfer rate from insoluble to 
                                             ! soluble modes. 

  !--- 10) Computational constants: ------------------------------------------------

  REAL(KIND=JPRB), PARAMETER :: sqrt2=1.4142136,  pi=3.141592654

  !--- 11) Data used for the calculation of the aerosol properties -----------------
  !       under ambient conditions:
  !       (Included the conversion from Pa to hPa in the first parameter.)

  REAL(KIND=JPRB), PARAMETER :: wvb(17)=                                        &
             &       (/   95.80188,     -28.5257,     -1.082153,     0.1466501, &
             &           -20627.51,    0.0461242,     -0.003935,      -3.36115, &
             &         -0.00024137,  0.067938345, 0.00000649899,   8616124.373, &
             &         1.168155578, -0.021317481,   0.000270358, -1353332314.0, &
             &        -0.002403805                                              /)

  REAL(KIND=JPRB), PARAMETER :: gmb(9)=                                      &
             &       (/ 1.036391467, 0.00728531, -0.011013887, -0.068887407, &
             &          0.001047842, 0.001049607, 0.000740534, -1.081202685, &
             &         -0.0000029113                                         /)

  !--- 4) Logical mask for coagulation kernel: -------------------------------------
  !       (The coagulation kernel mask is symmetric and not all 
  !        values are used for physical considerations. As its 
  !        calculation is very expensive, a mask is used to 
  !        calculate only the necessarey elements.)

  LOGICAL :: locoagmask(nmod,nmod)

  DATA locoagmask(1:nmod,1) / .TRUE.,  .TRUE.,  .TRUE.,  .TRUE.,  .TRUE.,  .TRUE.,  .TRUE.  /

  DATA locoagmask(1:nmod,2) / .FALSE., .TRUE.,  .TRUE.,  .TRUE.,  .TRUE.,  .TRUE.,  .TRUE.  /

  DATA locoagmask(1:nmod,3) / .FALSE., .FALSE., .TRUE.,  .FALSE., .TRUE.,  .FALSE., .FALSE. /

  DATA locoagmask(1:nmod,4) / .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE. /

  DATA locoagmask(1:nmod,5) / .FALSE., .FALSE., .FALSE., .FALSE., .TRUE.,  .FALSE., .FALSE. /

  DATA locoagmask(1:nmod,6) / .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE. /

  DATA locoagmask(1:nmod,7) / .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE. /
  






  ! TM5M7-aerosol-specific elements, taken from chem_param.F90

  !
  ! component indices and families, to identify tracer
  !
  ! RCHG -> Should be each of these identifers numbers equal to
  !         the "i" of YAERO_NL(i) in the namelist? 
  INTEGER(KIND=JPIM), PARAMETER :: ISO4 = 1
  INTEGER(KIND=JPIM), PARAMETER :: INH4  = 2  ! Note: these are tracer fields in the scope of aerosol module. 
  INTEGER(KIND=JPIM), PARAMETER :: INO3_A = 3 ! Check how to treat when running with chemistry.  
  INTEGER(KIND=JPIM), PARAMETER :: IACS_N = 4
  INTEGER(KIND=JPIM), PARAMETER :: ISO4ACS = 5
  INTEGER(KIND=JPIM), PARAMETER :: IBCACS = 6
  INTEGER(KIND=JPIM), PARAMETER :: IPOMACS = 7
  INTEGER(KIND=JPIM), PARAMETER :: ISSACS = 8
  INTEGER(KIND=JPIM), PARAMETER :: IDUACS = 9
  INTEGER(KIND=JPIM), PARAMETER :: ISOANUS= 10
  INTEGER(KIND=JPIM), PARAMETER :: ISOAAIS = 11
  INTEGER(KIND=JPIM), PARAMETER :: ISOAACS = 12
  INTEGER(KIND=JPIM), PARAMETER :: ISOACOS = 13
  INTEGER(KIND=JPIM), PARAMETER :: ISOAAII = 14
  INTEGER(KIND=JPIM), PARAMETER :: IH2OPART = 15
  INTEGER(KIND=JPIM), PARAMETER :: IAII_N = 16
  INTEGER(KIND=JPIM), PARAMETER :: IBCAII = 17
  INTEGER(KIND=JPIM), PARAMETER :: IPOMAII = 18
  INTEGER(KIND=JPIM), PARAMETER :: IACI_N = 19
  INTEGER(KIND=JPIM), PARAMETER :: IDUACI = 20
  INTEGER(KIND=JPIM), PARAMETER :: IAIS_N = 21
  INTEGER(KIND=JPIM), PARAMETER :: ISO4AIS = 22
  INTEGER(KIND=JPIM), PARAMETER :: IBCAIS = 23
  INTEGER(KIND=JPIM), PARAMETER :: IPOMAIS = 24
  INTEGER(KIND=JPIM), PARAMETER :: ICOI_N = 25
  INTEGER(KIND=JPIM), PARAMETER :: IDUCOI = 26
  INTEGER(KIND=JPIM), PARAMETER :: ICOS_N= 27
  INTEGER(KIND=JPIM), PARAMETER :: ISO4COS = 28
  INTEGER(KIND=JPIM), PARAMETER :: IBCCOS= 29
  INTEGER(KIND=JPIM), PARAMETER :: IPOMCOS = 30
  INTEGER(KIND=JPIM), PARAMETER :: ISSCOS = 31
  INTEGER(KIND=JPIM), PARAMETER :: IDUCOS = 32
  INTEGER(KIND=JPIM), PARAMETER :: INUS_N = 33
  INTEGER(KIND=JPIM), PARAMETER :: ISO4NUS = 34
  INTEGER(KIND=JPIM), PARAMETER :: IELVOC = 35
  INTEGER(KIND=JPIM), PARAMETER :: IISVOC = 36
  INTEGER(KIND=JPIM), PARAMETER :: IMSA = 37


  !
  ! molar weights of selected components.
  !
  REAL(KIND=JPRB), PARAMETER :: xmair=28.94 ! mass of air, g/mol
  REAL(KIND=JPRB), PARAMETER :: xmh=1.0079
  REAL(KIND=JPRB), PARAMETER :: xmn=14.0067
  REAL(KIND=JPRB), PARAMETER :: xmc=12.01115
  REAL(KIND=JPRB), PARAMETER :: xms=32.064
  REAL(KIND=JPRB), PARAMETER :: xmo=15.9994
  REAL(KIND=JPRB), PARAMETER :: xmna=22.990
  REAL(KIND=JPRB), PARAMETER :: xmcl=35.453


  REAL(KIND=JPRB), PARAMETER :: xmno3=xmn+3.*xmo
  REAL(KIND=JPRB), PARAMETER :: xmh2so4=2.*xmh+xms+4.*xmo
  REAL(KIND=JPRB), PARAMETER :: xmdust=xmair
  REAL(KIND=JPRB), PARAMETER :: xmnumb=xmair
  ! attention xmso2: conversion emissions done when added...
  REAL(KIND=JPRB), PARAMETER :: xmso2=xms+2.*xmo
  REAL(KIND=JPRB), PARAMETER :: xmdms=xms+2*xmc+6*xmh
  ! attention xmnh3: conversion emissions when added...
  REAL(KIND=JPRB), PARAMETER :: xmnh3=xmn+3.*xmh
  ! attention: conversion emissions when added...
  REAL(KIND=JPRB), PARAMETER :: xmnh4=xmn+4.*xmh
  REAL(KIND=JPRB), PARAMETER :: xmmsa=xms+xmc+3*xmo+4*xmh
  REAL(KIND=JPRB), PARAMETER :: xmnh2=xmn+xmh*2.
  REAL(KIND=JPRB), PARAMETER :: xmnh2o2=xmnh2+2.*xmo
  REAL(KIND=JPRB), PARAMETER :: xmso4=xms+4.*xmo



  ! densities (kg/m3) used in emission and/or optics routines
  REAL(KIND=JPRB), PARAMETER         :: density_ref = 1800.0   ! for 'reference' density calculations
  REAL(KIND=JPRB), PARAMETER         :: ss_density   = dnacl * 1.e3
  REAL(KIND=JPRB), PARAMETER         :: dust_density = ddust * 1.e3
  REAL(KIND=JPRB), PARAMETER         :: carbon_density = dbc * 1.e3
  REAL(KIND=JPRB), PARAMETER         :: pom_density = doc * 1.e3 ! Note that doc actually is the density of POM not OC
  REAL(KIND=JPRB), PARAMETER         :: soa_density = pom_density ! TB first order approx. same as pom
  ! H2-SO4 particle density:
  REAL(KIND=JPRB), PARAMETER         :: so4_density = dh2so4 * 1.e3  
  REAL(KIND=JPRB), PARAMETER         :: h2so4_factor = xmh2so4 / xmso4
  ! Ammonium-nitrate particle density  used in the optics routine:
  ! Value based on Lowenthal et al. (Atmos. Environ., 2000) (see also De Meij et al., ACP, 2006).
  !REAL(KIND=JPRB), PARAMETER         :: nh4no3_density = 1700.  
  REAL(KIND=JPRB), PARAMETER         :: nh4no3_density = 1.73 * 1.e3
  REAL(KIND=JPRB), PARAMETER         :: nh4no3_factor = (xmnh4+xmno3)/xmno3
  REAL(KIND=JPRB), PARAMETER         :: msa_density = 1.48 * 1.e3

  ! Kappa values
  REAL(KIND=JPRB), PARAMETER         :: Kap_su = 0.6
  REAL(KIND=JPRB), PARAMETER         :: Kap_pom = 0.1
  REAL(KIND=JPRB), PARAMETER         :: Kap_soa = 0.1
  REAL(KIND=JPRB), PARAMETER         :: Kap_bc = 0.
  REAL(KIND=JPRB), PARAMETER         :: Kap_ss = 1.0
  REAL(KIND=JPRB), PARAMETER         :: Kap_du = 0.
  REAL(KIND=JPRB), PARAMETER         :: Kap_na2so4 = 0.95
  REAL(KIND=JPRB), PARAMETER         :: Kap_msa = 0.6
  REAL(KIND=JPRB), PARAMETER         :: Kap_no3 = 0.6

  REAL(KIND=JPRB), DIMENSION(NMOD), PARAMETER    ::  sigma_lognormal = (/ 1.59, 1.59, 1.59, 2.00, 1.59, 1.59, 2.00 /)
  !
  ! mode numbers
  !
  INTEGER(KIND=JPIM), PARAMETER :: mode_nuc = 1
  INTEGER(KIND=JPIM), PARAMETER :: mode_ais = 2
  INTEGER(KIND=JPIM), PARAMETER :: mode_acs = 3
  INTEGER(KIND=JPIM), PARAMETER :: mode_cos = 4
  INTEGER(KIND=JPIM), PARAMETER :: mode_aii = 5
  INTEGER(KIND=JPIM), PARAMETER :: mode_aci = 6
  INTEGER(KIND=JPIM), PARAMETER :: mode_coi = 7

  !  mode_number => mode_start
  INTEGER(KIND=JPIM), PARAMETER :: mode_start   (nmod) = (/  inus_n,  iais_n,  iacs_n,  icos_n,  iaii_n, iaci_n, icoi_n /)  ! first tracer in mode
  INTEGER(KIND=JPIM), PARAMETER :: mode_end_so4 (nmod) = (/ iso4nus, iso4ais, iso4acs, iso4cos,       0,      0,      0 /)
  INTEGER(KIND=JPIM), PARAMETER :: mode_end_bc  (nmod) = (/       0,  ibcais,  ibcacs,  ibccos,  ibcaii,      0,      0 /)
  INTEGER(KIND=JPIM), PARAMETER :: mode_end_pom (nmod) = (/       0, ipomais, ipomacs, ipomcos, ipomaii,      0,      0 /)
  INTEGER(KIND=JPIM), PARAMETER :: mode_end_ss  (nmod) = (/       0,       0,  issacs,  isscos,       0,      0,      0 /)
  INTEGER(KIND=JPIM), PARAMETER :: mode_end_dust(nmod) = (/       0,       0,  iduacs,  iducos,       0, iduaci, iducoi /)
  INTEGER(KIND=JPIM), PARAMETER :: mode_end_soa (nmod) = (/ isoanus, isoaais, isoaacs, isoacos, isoaaii,      0,      0 /)
  INTEGER(KIND=JPIM), PARAMETER :: mode_nm      (nmod) = (/       2,       4,       6,       6,       3,      1,      1 /)  ! # tracers in mode
  INTEGER(KIND=JPIM), PARAMETER :: mode_nm_sed  (nmod) = (/       2,       4,       9,       6,       3,      1,      1 /)  ! # tracers in mode
  INTEGER(KIND=JPIM), PARAMETER :: mode_end     (nmod) = mode_start + mode_nm                                               ! last tracer in mode
  INTEGER(KIND=JPIM), PARAMETER :: mode_tracers(0:6,nmod) = &
                    &   RESHAPE( (/ inus_n, iso4nus, isoanus, 0, 0, 0, 0, &
                    &               iais_n, iso4ais, ibcais, ipomais,isoaais, 0, 0, &
                    &               iacs_n, iso4acs, ibcacs, ipomacs, issacs, iduacs, isoaacs, &
                    &               icos_n, iso4cos, ibccos, ipomcos, isscos, iducos, isoacos, &
                    &               iaii_n, ibcaii,  ipomaii,isoaaii, 0, 0, 0, &
                    &               iaci_n, iduaci, 0, 0, 0, 0, 0, &
                    &               icoi_n, iducoi, 0, 0, 0, 0, 0 /), (/ 7, nmod/) )
  INTEGER(KIND=JPIM), PARAMETER :: mode_tracers_by_mods(0:6,nmod) = &
                    &   RESHAPE( (/ inus_n, iso4nus, 0     , 0      , isoanus, 0      , 0      , &
                    &               iais_n, iso4ais, ibcais, ipomais, isoaais, 0      , 0      , &
                    &               iacs_n, iso4acs, ibcacs, ipomacs, isoaacs, issacs , iduacs, &
                    &               icos_n, iso4cos, ibccos, ipomcos, isoacos, isscos , iducos, &
                    &               iaii_n, 0      , ibcaii, ipomaii, isoaaii, 0      , 0      , &
                    &               iaci_n, 0      , 0     , 0      , 0      , 0      , iduaci ,  &
                    &               icoi_n, 0      , 0     , 0      , 0      , 0      , iducoi  /), (/ 7, nmod/) )
  INTEGER(KIND=JPIM), PARAMETER :: mode_tracers_sed(0:9,nmod) = RESHAPE( (/&
  &     inus_n, iso4nus, isoanus, 0,       0,       0,      0,       0,    0,      0,    &
  &     iais_n, iso4ais, ibcais,  ipomais, isoaais, 0,      0,       0,    0,      0,    &
  &     iacs_n, iso4acs, ibcacs,  ipomacs, issacs,  iduacs, isoaacs, inh4, ino3_a, imsa, &
  &     icos_n, iso4cos, ibccos,  ipomcos, isscos,  iducos, isoacos, 0,    0,      0,    &
  &     iaii_n, ibcaii,  ipomaii, isoaaii, 0,       0,      0,       0,    0,      0,    &
  &     iaci_n, iduaci,  0,       0,       0,       0,      0,       0,    0,      0,    &
  &     icoi_n, iducoi,  0,       0,       0,       0,      0,       0,    0,      0/),  (/ 10, nmod/) )



  ! number of aerosol bins used for deposition:
  integer, parameter  ::  nrdep = 23
  
  ! aerosol radii used for each bin:
  real, parameter     ::  lur(nrdep) = &
                 (/  0.001,  0.01,   0.05,  0.1,  0.3,  &
                     0.5  ,  0.7 ,   0.8 ,  0.9,  1.0,  &
                     1.2  ,  1.5 ,   2.0 ,  3.0,  4.0,  &
                     5.0  ,  6.0 ,   8.0 , 10.0, 15.0,  &
                    20.0  , 50.0 , 100.0                /)



  ! ********************************************************************
  ! wet deposition
  ! ********************************************************************

  ! nscav       : selected species for scavenging
  ! nscav_index : index for scavenging:
  ! nscav_type  : type of scavenging:
  !               0 no scavenging
  !               1 scavenging 100 % solubility assumed
  !               2 scavenging henry solubility assumed
  !               3 scavenging, aerosol removal assumed
  !               4 scavenging, special case for SO2 with aq phase diss.
  !
  integer,parameter                    :: nscav=30

  integer,parameter,dimension(nscav)   :: nscav_index  = (/ &
                        inus_n,  iais_n,  iacs_n,  icos_n, iaii_n, iaci_n, icoi_n, &
                        iso4nus, iso4ais, iso4acs, iso4cos, &
                        ibcais,  ibcacs,  ibccos,  ibcaii,  &
                        ipomais, ipomacs, ipomcos, ipomaii, &
                        issacs,  isscos, &
                        iduacs,  iducos,  iduaci,  iducoi,  &
                        isoanus, isoaais, isoaacs, isoacos, isoaaii  &
       /)

  ! nscav_type = 5  : nu mode soluble aerosol
  ! nscav_type = 6  : ai mode soluble aerosol
  ! nscav_type = 7  : ac mode soluble aerosol
  ! nscav_type = 8  : co mode soluble aerosol
  ! nscav_type = 9  : ai mode insoluble aerosol
  ! nscav_type = 10 : ac mode insoluble aerosol
  ! nscav_type = 11 : co mode insoluble aerosol

! in m7-version so4 is treated as gas-phase sulphuric acid (H2SO4) (scav-type 2)
! Since ammonium-nitrate and MSA are assumed to be in the soluble accumulation mode,
! their scavenging efficiency has been changed to the value for that mode.
! For consistency, the same value is used for Pb210.
  integer, dimension(nscav),parameter :: &
       nscav_type = (/    &
                          5,    6,    7,    8,    9,    10,   11,  &     ! particle number
                          5,    6,    7,    8,                     &     ! sulphate mass
                          6,    7,    8,    9,                     &     ! BC mass
                          6,    7,    8,    9,                     &     ! POM mass
                          7,    8,                                 &     ! SS mass
                          7,    8,    10,   11,                    &     ! DUST mass
                          5,    6,    7,    8,    9     /)               ! SOA mass

END MODULE TM5M7_DATA
