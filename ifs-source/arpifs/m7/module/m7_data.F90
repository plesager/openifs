MODULE M7_DATA

USE PARKIND1, ONLY : JPIM  , JPRB

  ! *M7_DATA* contains phyiscal switches and parameters 
  !           for the ECHAM/HAM aerosol model.
  !
  ! Author:
  ! -------
  ! Philip Stier, MPI-MET                    12/2002
  !

USE TM5M7_DATA,       ONLY: nmod

IMPLICIT NONE

  !--- 0) Submodel ID:

  INTEGER(KIND=JPIM) :: id_ham

  !--- 1) Switches:

  !--- 1.1) Physical:

  !--- Define control variables and pre-set with default values: 

  LOGICAL :: lm7         = .TRUE.      ! Aerosol dynamics and thermodynamics scheme M7


  INTEGER(KIND=JPIM) :: ncdnc       = 0,       &  ! CDNC activation scheme:
                                       !
                                       !    ncdnc = 0  OFF => standard ECHAM5
                                       !
                                       !          = 1  Lohmann et al. (1999) + Lin and Leaitch (1997)
                                       !          = 2  Lohmann et al. (1999) + Abdul-Razzak and Ghan (2000)
                                       !          = 3  Lohmann et al. (1999) + ( Nenes et al. (2003) ) 
                                       !
             nicnc       = 0,       &  ! ICNC scheme:
                                       !
                                       !    ncdnc = 0  OFF
                                       !          = 1  Kaercher and Lohmann (2002)
                                       !
             nauto       = 1,       &  ! Autoconversion scheme:
                                       !
                                       !    nauto = 1  Beheng (1994) - ECHAM5 Standard
                                       !          = 2  Khairoutdinov and Kogan (2000)
                                       !
             ndust       = 2,       &  ! Dust emission scheme:
                                       ! 
                                       !    ndust = 1  Balkanski et al. (2002)
                                       !          = 2  Tegen et al. (2002)
                                       !
!mo_ham.f90 is used             nseasalt    = 2,       &  ! Sea Salt emission scheme:
!mo_ham.f90 is used                                       ! 
!mo_ham.f90 is used                                       !    nseasalt = 1  Monahan (1986)
!mo_ham.f90 is used                                       !             = 2  Schulz et al. (2002)
!mo_ham.f90 is used                                       !
             npist       = 3,       &  ! DMS emission scheme:
                                       !
                                       !    npist = 1 Liss & Merlivat (1986) 
                                       !          = 2 Wanninkhof (1992)
                                       !          = 3 Nightingale (2000)
                                       !
             nemiss      = 1,       &  ! Emission inventory
                                       !    
                                       !    nemiss =1 old version
                                       !    nemiss =2 AEROCOM emissions 2000
                                       !
             nsoa        = 2           ! SOA formation scheme:
                                       !
                                       !    nsoa = 0 POM mass emission into both Aitken modes (standard TM5)
                                       !           1 POM mass emission + distribution according to volatility assumptions to 5 modes
                                       !           2 atmospheric formation from precursors + distribution according to volatility assumptions to 5 modes

  LOGICAL :: lodiag      = .FALSE.     ! Extended diagnostics

  LOGICAL :: laero_rad   = .FALSE.     ! Radiation calculation

  LOGICAL :: lorad(nmod) = .FALSE.     !    switch for each mode

  LOGICAL :: lodiagrad   = .FALSE.     ! Extended radiation diagnostics

  INTEGER(KIND=JPIM) :: nwv         = 0           !    nwv: number of additional wavelengths
                                       !         for the radiation calculations
                                       !         (max currently set to 10)

  REAL(KIND=JPRB)    :: cwv(10)     = 0.          !    cwv: array of additional wavelengths
                                       !         for the radiation calculations [m]

  LOGICAL :: lomassfix   = .TRUE.      ! Mass fixer in convective scheme

  !--- 1.2) Technical:

!!$  INTEGER(KIND=JPIM) :: NFILETYPE   = GRIB        ! Output stream filetypes


  !--- 2) Parameters:

  !-- 2.1) Number of aerosol compounds: (needs to be harmonized with nmode in mo_aero_m7)

  INTEGER(KIND=JPIM), PARAMETER :: ntype=6

  !--- 2.2) Mode names:

  CHARACTER(LEN=2), PARAMETER :: cmode(nmod)=(/'NS','KS','AS','CS','KI','AI','CI'/)

  !--- 2.3) Compound names:

  CHARACTER(LEN=3), PARAMETER :: ctype(ntype)=(/'SO4','BC ','OC ','SS ','DU ','WAT'/)

  !--- 2.4) Index field of tracer indices for the aerosol numbers in each mode:

  INTEGER(KIND=JPIM) :: nindex(nmod)

  !--- 2.5) Emissions:

  !--- Carbon Emissions


!  REAL(KIND=JPRB), PARAMETER         :: zbb_wsoc_perc  = 0.65,      & ! Biom. Burn. Percentage of Water Soluble OC (WSOC) [1]
                                                           ! (M.O. Andreae; Talk: Smoke and Climate)

  REAL(KIND=JPRB), PARAMETER         :: zbb_wsoc_perc  = 0.95,      & ! TB:
                                                           ! To reduce the AOD over china and outflow region of 
                                                           ! Africa the water soluble fraction was increasde to 95%
                                                           ! in preparation for CMIP6.
                                                           !
                                                           ! Some basis for the choice can be found here:
                                                           ! e (e.g. Janhall et al., 2010; 
                                                           ! https://doi.org/10.5194/acp-10-1427-2010 ; Winijkul et al., 2015;
                                                           ! https://doi.org/10.1016/j.atmosenv.2015.02.037; Li et al., 2009;
                                                           ! https://pubs.acs.org/doi/abs/10.1021/es803330j).


                             zbge_wsoc_perc = 0.65,      & ! Assume same Percentage of WSOC for biogenic OC
                             !>>> TvN
                             ! The value of 1.4 for the POM to OC mass ratio is an outdated estimate.
                             ! In the current code we can apply different ratios
                             ! for emissions from different sources.
                             ! For further details, see comment in emission_data.F90.     
                             ! The use of a single constant value, on the other hand,
                             ! would have the advantage that the simulated POM concentrations
                             ! can easily be converted to OC.
                             ! An average value of 1.8 seems reasonable.
                             ! Assuming that there are no substantial contributions from
                             ! elements other than H and O, a value of 1.8 can be obtained
                             ! with an H:C atomic ratio of 1.6 and and O:C ratio of 0.5,
                             ! which are well within the range of oxidation states 
                             ! presented by Heald et al. (GRL, 2010).
                             ! According to the model of Kuwata et al. (Environ. Sci. Technol., 2012),
                             ! the resulting particle density would be close to the value 
                             ! assumed in the model (doc = 1.3 g/cm3 in mo_aero_m7.F90).
                             !zom2oc         = 1.4,       & ! Mass ratio organic species to organic carbon
                                                           ! (Seinfeld and Pandis, 1998, p709;
                                                           !  Ferek et al., JGR, 1998) 
                             !
                             ! The emission radii for carbonaceous aerosols of the original code below
                             ! correspond to the values recommended by AeroCom (Dentener et al., ACP, 2006),
                             ! but adapted to sigma = 1.59 as used in M7 (Stier et al., ACP, 2005).
                             cmr_ff         = 0.03E-6,   & ! Fossil fuel emissions:
                                                           ! assumed number median radius of the emitted
                                                           ! particles with the standard deviation given
                                                           ! in mo_aero_m7 [m]. Has to lie within the 
                                                           ! Aitken mode for the current setup!
                             cmr_bb         = 0.075E-6,  & ! Biomass burning emissions:
                                                           ! Assumed number median radius of the emitted
                                                           ! particles with the standard deviation given
                                                           ! in mo_aero_m7 [m]. Has to lie within the 
                                                           ! Accumulation mode for the current setup!
                             cmr_bg         = 0.03E-6,    &! Biogenic secondary particle formation:
                                                           ! Assumed number median radius of the emitted
                                                           ! particles with the standard deviation given
                                                           ! in mo_aero_m7 [m]. Has to lie within the 
                                                           ! Aitken mode for the current setup!
                             cmr_sk         = 0.03E-6,    &! SO4 primary emission  ---> aitken mode
                                                           ! Assumed number median radius of the emitted
                                                           ! particles with the standard deviation given
                                                           ! in mo_aero_m7 [m]. Has to lie within the 
                                                           ! Aitken mode for the current setup!
                             cmr_sa         = 0.075E-6,   &! SO4 primary emission  ---> accumulation mode
                                                           ! Assumed number median radius of the emitted
                                                           ! particles with the standard deviation given
                                                           ! in mo_aero_m7 [m]. Has to lie within the 
                                                           ! Accumulation mode for the current setup!
                             cmr_sc         = 0.75E-6,    &! SO4 primary emission  ---> coarse mode
                                                           ! Assumed number median radius of the emitted
                                                           ! particles with the standard deviation given
                                                           ! in mo_aero_m7 [m]. Has to lie within the 
                                                           ! Coarse mode for the current setup!
                             facso2         = 0.975,      &! factor to scale primary SO4 emissions 
                                                           ! AEROCOM assumption 2.5 % of the SO2 emissions 
                                                           ! in the from of SO4
                             so2ts          = 1./1.998     ! conversion factor SO2 to S

  REAL(KIND=JPRB), PUBLIC :: zm2n_bcki_ff, zm2n_bcki_bb, &
                             zm2n_bcks_bb, zm2n_ocki_ff, &
                             zm2n_ocki_bb, zm2n_ocki_bg, &
                             zm2n_ocks_bb, zm2n_ocks_bg, &
                             zm2n_s4ks_sk, zm2n_s4as_sa, &
                             zm2n_s4cs_sc


END MODULE M7_DATA
