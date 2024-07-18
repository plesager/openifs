MODULE TM5M7_EMIS_DATA

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE M7_DATA, ONLY : zbb_wsoc_perc, zbge_wsoc_perc, cmr_ff,cmr_bb, facso2, &
    & cmr_sk, cmr_sa, cmr_sc
USE TM5M7_DATA, ONLY : XMAIR    

IMPLICIT NONE

SAVE
  
! Array to collect emissions. Note that this is a 3D array, and not 4D as in TM5.
TYPE MODAL_EMISSIONS
  REAL(KIND=JPRB), DIMENSION(:,:,:), POINTER :: d3 ! KLON, BB_LM, mode_nm(mode)
END TYPE MODAL_EMISSIONS


  ! Count median / geometric mean radii of 
  ! primary carbonaceous and sulfate emissions
  ! based on AeroCom-I recommendations (Dentener et al., ACP, 2006).
  ! The corresponding values for the M7 modes are given by Stier et al. (ACP, 2005). 
  ! The same values have also been adopted in GLOMAP (Mann et al., 2010).
  !
  ! Count median radii for carbonaceous aerosol emissions from Dentener et al.,
  ! corresponding to sigma = 1.8:
  !REAL(KIND=JPRB), PARAMETER         :: rad_emi_ff = 0.015e-6
  !REAL(KIND=JPRB), PARAMETER         :: rad_emi_bb = 0.04e-6
  ! The corresponding values for sigma = 1.59 used in M7
  ! are cmr_ff and cmr_bb, specified in m7_data.F90.

  ! For comparison, Bond et al. (JGR, 2013) give number median radii
  ! between 25 and 40 nm for fresh BC in the urban areas 
  ! of Tokyo, Nagoya, and Seoul, 
  ! of 60 nm in plumes associated with wildfires,
  ! and about 15 nm from aircraft jet engines.
  ! These values are volume-equivalent radii (see their Fig. 4).
  !
  ! According to the original paper by Schwarz et al. (GRL, 2008) 
  ! the corresponding geometric standard deviation
  ! is sigma = 1.71 for the urban BC 
  ! and 1.43 for the biomass burning aerosol.
  !
  ! For BC in biomass burning plumes,
  ! Kondo et al (JGR, 2011) estimated
  ! number median radii in the range 68-70.5 nm (+- 6-8 nm)
  ! and geometric standard deviation between 1.32 and 1.36 (+- 0.01-0.04),
  ! for particles thickly coated by organics.
  !
  ! Janhaell et al. (ACP, 2010) have compiled measurements of
  ! particle size in fresh biomass burning smoke from vegetation fires.
  ! They mention that particles from biomass burning are dominated
  ! by an accumulation mode.
  ! They also present a relation between the geometric mean diameter Dg 
  ! and geometric standard deviation sigma for fresh smoke:
  ! Dg (um) = (584 +- 5) - (269 +-1) sigma
  ! This gives a geometric mean radius of 78 um for sigma = 1.59,
  ! in close agreement with the value used by Stier et al.
  !
  ! Another way to account for differences in sigma
  ! is to modify the number median radius rg
  ! such that the number of particles (N) 
  ! emitted for a certain mass (M) is the same.
  ! N is proportional to M/rv^3,
  ! where rv is the volume mean radius, 
  ! which is related to rg by
  ! rv = rg * exp(1.5*(ln(sigma))^2).
  ! According to Janhall, fresh smoke has an average
  ! Dg of 117 +- 13 nm and sigma of 1.7 +- 0.1.
  ! At sigma=1.59, this would translate into Dg of about 129 nm,
  ! or rg of about 65 nm.
  ! For the estimates from Kondo et al. and Schwarz et al.,
  ! this would give somewhat smaller rg values 
  ! of about 58 and 53 nm, resp.
  !
  ! Particles emitted by grass and savannah fires are generally
  ! somewhat smaller than those from wood burning.
  ! Janhaell et al. estimate that the mean emission radii
  ! for grass and savannah fires, resp., are 12.5 and 10 nm smaller.
  ! These differences is not accounted for in the model.
  !
  ! In a later version of ECHAM-HAM particles the emission radius
  ! for biomass burning was reduced to the value for fossil fuel
  ! (Zhang et al., ACP, 2012).  
  ! However, such as a small value seems inconsistent with measurements.
  !
  ! In the CMIP6 emission data set,
  ! the contributions from solid biofuel combustion
  ! are included in the 2-D anthropogenic sectors,
  ! and provided separately in a supplementary data set.
  ! In CMIP6, biofuel is only non-zero
  ! for the energy, industry and residential and
  ! transportation sectors.
  ! Most of the residential emissions are due
  ! to solid biofuel combustion.
  ! The contribution to the industrial sector
  ! is only substantial in developing countries,
  ! while the contributions to the energy 
  ! and transportation sectors are generally very small.
  ! 
  ! For inventories other than CMIP6,
  ! no distinction between fossil and biofuel emissions
  ! is made in the model.
  ! 
  ! Winijkul et al. (Atm. Env., 2015) have measured
  ! size distributions from energy-related combustion
  ! sources for the residential, industrial, power and 
  ! transportation sectors.
  ! They give regional and global estimates of
  ! mass median diameters.
  ! In their supplementary material they give
  ! an overview of results from other studies.
  ! These results indicate that the size distribution 
  ! for both biofuel and fossil fuel combustion sources
  ! are strongly dependent on the technique.
  ! Count median radii for the residential sector,
  ! presented in their Table S1, vary between about
  ! 15 nm for modern (improved) wood-fueled stoves,
  ! and 25-30 for fireplaces,
  ! to 160 nm for (regular) wood-fueled heating stoves,
  ! to 270 nm for traditional cookstoves.

  ! General fossil fuel (2-d sectors)
  ! emitted in the Aitken mode
  REAL(KIND=JPRB), PARAMETER         :: rad_emi_ff_sol    = cmr_ff  
  REAL(KIND=JPRB), PARAMETER         :: rad_emi_ff_insol  = cmr_ff
  
  ! Energy sector
  ! emitted in the Aitken mode.
  REAL(KIND=JPRB), PARAMETER         :: rad_emi_ene_sol   = cmr_ff
  REAL(KIND=JPRB), PARAMETER         :: rad_emi_ene_insol = cmr_ff

  ! Industry sector
  ! emitted in the Aitken mode.
  REAL(KIND=JPRB), PARAMETER         :: rad_emi_ind_sol   = cmr_ff
  REAL(KIND=JPRB), PARAMETER         :: rad_emi_ind_insol = cmr_ff

  ! Transportation sector
  ! emitted in the Aitken mode.
  REAL(KIND=JPRB), PARAMETER         :: rad_emi_tra_sol   = cmr_ff
  REAL(KIND=JPRB), PARAMETER         :: rad_emi_tra_insol = cmr_ff

  ! Shipping sector
  ! emitted in the Aitken mode.
  ! Currently set to cmr_ff,
  ! but could be changed.
  REAL(KIND=JPRB), PARAMETER         :: rad_emi_shp_sol   = cmr_ff
  REAL(KIND=JPRB), PARAMETER         :: rad_emi_shp_insol = cmr_ff

  ! Aircraft sector
  ! emitted in the Aitken mode.
  ! Currently set to cmr_ff,
  ! but a smaller value could be used,
  ! e.g. 10 or 15 nm.
  REAL(KIND=JPRB), PARAMETER         :: rad_emi_air_sol   = cmr_ff
  REAL(KIND=JPRB), PARAMETER         :: rad_emi_air_insol = cmr_ff

  ! Open biomass burning 
  ! Soluble part emitted in the accumulation mode,
  ! insoluble part in the Aitken mode.
  ! Stier et al. apply cmr_ff to the insoluble part,
  ! but a slightly larger value seems more realistic,
  ! e.g. 40 nm.
  REAL(KIND=JPRB), PARAMETER         :: rad_emi_bb_sol   = cmr_bb
  REAL(KIND=JPRB), PARAMETER         :: rad_emi_bb_insol = cmr_ff

  ! Solid biofuel combustion
  ! Currently treated as biomass burning,
  ! but could be changed (see comments above).
  ! The value cmr_bb = 75 nm corresponds
  ! with the value of 50 nm at sigma=2
  ! assumed by Kondros et al. (ACP, 2015)
  ! for emissions from biofuel combusion
  ! in their baseline run.
  ! In one of their sensitivity runs,
  ! they increase it by a factor of 2.
  REAL(KIND=JPRB), PARAMETER         :: rad_emi_bf_sol   = cmr_bb 
  REAL(KIND=JPRB), PARAMETER         :: rad_emi_bf_insol = cmr_ff
 
  ! not used anymore:
  !REAL(KIND=JPRB), PARAMETER         :: rad_soa = 0.01e-6      ! soa average radius 
                                                    ! assuming 3nm particle formation and growth to
                                                    ! that size in half an hour 

  ! Count median radii for sulfate aerosol emissions adapted to the M7 modes:
  REAL(KIND=JPRB), PARAMETER         :: rad_so4_ait = cmr_sk  ! aitken mode radius
  REAL(KIND=JPRB), PARAMETER         :: rad_so4_acc = cmr_sa  ! accumulation mode radius
  REAL(KIND=JPRB), PARAMETER         :: rad_so4_coa = cmr_sc  ! coarse mode radius

  ! Count median dry radii for sea salt emissions.
  ! These values have been updated 
  ! following Vignati et al. (Atmos. Environ., 2010)
  ! For further explanations, see emission_ss.F90.
  !REAL(KIND=JPRB), PARAMETER         :: radius_ssa = 0.0794e-6
  !REAL(KIND=JPRB), PARAMETER         :: radius_ssc = 0.63e-6
  REAL(KIND=JPRB), PARAMETER         :: radius_ssa = 0.09e-6  ! accumulation mode
  REAL(KIND=JPRB), PARAMETER         :: radius_ssc = 0.794e-6 ! coarse mode

  ! Soluble fraction of POM mass
  ! According to Janhall et al. (ACP, 2010)
  ! 40 to 80% of the organic matter from
  ! vegetation fires is water soluble.
  ! Kondros et al. (ACP, 2015) use 80%
  ! for POM from biofuel combustion
  ! in their base run.
  ! For the moment, emisions from
  ! biofuel combustion and open biomass burning
  ! are treated in the same way.
  !
  ! open biomass burning
  REAL(KIND=JPRB), PARAMETER         :: frac_pom_sol_bb = zbb_wsoc_perc
  ! solid biofuel combustion
  REAL(KIND=JPRB), PARAMETER         :: frac_pom_sol_bf = zbb_wsoc_perc 
  !REAL(KIND=JPRB), PARAMETER         :: frac_pom_sol_bf = 0.8           ! alternative value
  ! fossil fuel combustion
  !REAL(KIND=JPRB), PARAMETER         :: frac_pom_sol_ff = 0.65          ! original value
  REAL(KIND=JPRB), PARAMETER         :: frac_pom_sol_ff = 0.0            ! value since 2015 revision

  ! Soluble fraction of BC mass
  ! In the original code, 
  ! fresh BC was assumed 100% insoluble,
  ! and therefore emitted into the Aitken mode.
  ! The new code allows to use a non-zero fraction,
  ! to account for non-resolved ageing close to the source.
  ! See e.g. Kondros et al.,
  ! who use 50% for biofuel combustion 
  ! in their base run.
  ! In the standard EMEP MSC-W model,
  ! 20% of the elemental carbon (EC) from
  ! anthropogenic sources 
  ! (representative of fossil fuel combustion)
  ! and all EC from open biomass fires
  ! is assumed hygroscopic.
  !
  ! open biomass burning
  !REAL(KIND=JPRB), PARAMETER         :: frac_bc_sol_bb = 0.0            ! original value
  !REAL(KIND=JPRB), PARAMETER         :: frac_bc_sol_bb = 0.5
  REAL(KIND=JPRB), PARAMETER         :: frac_bc_sol_bb = 0.95            ! To reduce the AOD over china and outflow region of 
                                                              ! Africa the water soluble fraction was increasde to 95%
                                                              ! in preparation for CMIP6.
  ! solid biofuel combustion
  !REAL(KIND=JPRB), PARAMETER         :: frac_bc_sol_bf = 0.0            ! original value
  !REAL(KIND=JPRB), PARAMETER         :: frac_bc_sol_bf = 0.5            ! Aerocom
  REAL(KIND=JPRB), PARAMETER         :: frac_bc_sol_bf = 0.95            ! TB:
                                                              ! To reduce the AOD over china and outflow region of 
                                                              ! Africa the water soluble fraction was increasde to 95%
                                                              ! in preparation for CMIP6.
                                                              !
                                                              ! Some basiss for the choice can be found here:
                                                              !  (e.g. Janhall et al., 2010; 
                                                              ! https://doi.org/10.5194/acp-10-1427-2010 ; Winijkul et al., 2015;
                                                              ! https://doi.org/10.1016/j.atmosenv.2015.02.037; Li et al., 2009;
                                                              ! https://pubs.acs.org/doi/abs/10.1021/es803330j).

  ! fossil fuel combustion
  REAL(KIND=JPRB), PARAMETER         :: frac_bc_sol_ff = 0.0

  ! Soluble fraction of surrogate SOA emissions
  ! POM from SOA is considered 65% soluble,
  ! as recommended by AeroCom.
  ! The paper by Kanakidou et al. (ACP, 2004)
  ! (mentioned in emission_pom.F90)
  ! does not seem to support 100% solubility,
  ! as was previously assumed.
  !REAL(KIND=JPRB), PARAMETER         :: frac_soa_sol    = 1.0           ! original value
  REAL(KIND=JPRB), PARAMETER         :: frac_soa_sol    = zbge_wsoc_perc ! value since 2015 revision

  ! Fraction of SOx mass emitted directly as sulfate
  REAL(KIND=JPRB), PARAMETER         :: frac_so4=1.-facso2


  ! The value of 1.4 for the POM to OC mass ratio,
  ! set in m7_data.F90, is an outdated estimate,
  ! see e.g. Turpin and Lim (Aerosol Sci. Technol., 2001) and
  ! Aiken et al. (Environ. Sci. Technol., 2008).
  ! Turpin and Lim estimate a ratio of 1.6 +- 0.2 for urban aerosol,
  ! and 2.1 +- 0.2 for aged (nonurban) aerosol;
  ! They also note that aerosols heavily impacted by woodsmoke can
  ! have an even higher ratio (2.2 to 2.6).
  ! According to Reid et al. (ACP, 2005),
  ! the POM to OC ratio in fresh biomass burning smoke is very uncertain,
  ! somewhere in between 1.4 and ~2.
  ! Aiken et al. measure ambient aerosol values 
  ! between 1.25 for O/C = 0 to 2.44 for O/C = 1.0.
  ! For OA from biomass burning, they measure 1.56-1.70,
  ! lower than the estimates from Turpin and Lim (~2.0).
  ! They find the highest ratios for 
  ! aged and freshly formed SOA (~2.4 and ~1.9, respectively)
  ! and lowest values for primary OA from urban combustion.
  ! Based on these studies, we apply different values
  ! for emissions from biomass burning versus other emissions,
  ! as is also done in some other models
  ! (see Table 1 in Tsigaridis et al., ACP, 2014),
  ! For SOA (seem emission_pom.F90)
  ! we apply a relatively high value valid for aged SOA.
  ! This compensates for the lack of SOA formation from isoprene,
  ! and improves the agreement with aerosol optical depth (AOD)
  ! derived from satellite observations (MODIS).
  !
  ! We could go even further and apply different values for the 
  ! water soluble and insoluble fractions (Turpin and Lim).
  ! 
  ! It should be acknowledged that the representation of OA 
  ! with a single tracer is very simplistic.
  ! In particular, increase in OA mass due to ageing
  ! is not properly accounted for.
  ! 
  !REAL(KIND=JPRB), PARAMETER  ::  oc2pom = zom2oc    !factor for conversion of OC mass to POM
  REAL(KIND=JPRB), PARAMETER  ::  oc2pom_ff = 1.6  ! fossil fuel
  REAL(KIND=JPRB), PARAMETER  ::  oc2pom_bf = 1.6  ! solid biofuel combustion
  REAL(KIND=JPRB), PARAMETER  ::  oc2pom_bb = 1.6  ! open biomass burning
  REAL(KIND=JPRB), PARAMETER  ::  oc2pom_soa = 2.4 ! SOA
  
  
  ! -----------
  ! Dust data
  ! -----------
  
  
    ! parameters for online dust calculations
    INTEGER(KIND=JPIM), PARAMETER              :: ntraced=8                     ! number of coarse-grained bins 
                                                                     ! in the original emission model
    INTEGER(KIND=JPIM), PARAMETER              :: nbin=24                       ! number of discretization points per bin
    INTEGER(KIND=JPIM), PARAMETER              :: nclass=ntraced*nbin           ! total number of discretization points
    INTEGER(KIND=JPIM), PARAMETER              :: nats=12                       ! number of soil types
    INTEGER(KIND=JPIM), PARAMETER              :: nmode=4                       ! number of particle size distributions in soils,
                                                                     ! which distinguishes between clay, silt, 
                                                                     ! medium/fine sand, and coarse sand
    INTEGER(KIND=JPIM), PARAMETER              :: nspe=nmode*3+2                ! for explanation, see below


  ! parameters for online emission input file ("onlinedust.nc")
  ! fields on 1x1 deg grid
  INTEGER(KIND=JPIM), PARAMETER :: nsoilph =  5, &
                        nfpar   = 12, &
                        nz0     = 13  ! number of {soilph, par, z0} fields
                                      ! entry nz0 indicates the annual mean.

    ! von Karman constant
    REAL(KIND=JPRB), PARAMETER                :: VKARMAN=0.4

    ! Constants used in the parameterization of the efficient friction velocity ratio,
    ! see Eqs. (17-20) in MB95: 
    REAL(KIND=JPRB), PARAMETER                 :: aeff=0.35                     
    REAL(KIND=JPRB), PARAMETER                 :: xeff=10.                     
    !
    ! -- scaling factor for threshold friction velocity
    ! u1fac is a tuning parameter necessary to obtain a reasonable global annual
    ! emission amount.  u1fac < 1 is used to reduce the threshold friction
    ! velocity.  In ECHAM-HAM simulations at T63 values of 0.86 and 0.56 were
    ! used by Cheng et al. (ACP, 2008).  The lower value was introduced to
    ! increase emissions when surface roughness lengths were increased from a
    ! constant value of 0.001 cm to values based on satellite measurements from
    ! Prigent et al. (JGR, 2005).  It is unclear where the value 0.66 specified
    ! below is based on.  In ECHAM-HAM2 (Zhang et al., ACP, 2012) the satellite
    ! based surface roughness values were abandoned again.
    REAL(KIND=JPRB), PARAMETER                 :: u1fac=0.6    ! 0.7 in EC-Earth 3.2.3
    
    REAL(KIND=JPRB), PARAMETER                 :: cd=1.2507E-06                 ! flux dimensioning parameter [g s^2/cm^4]
    
    !<<< TvN                                                                 ! (=roa/(grav*1.e2))
    ! ustar_min is not used:
    !REAL(KIND=JPRB), PARAMETER                 :: ustar_min=5.                  ! min. fricton velocity (cm/s)
    ! minimum surface roughness length z0 (cm)
    ! The minimum value in the data set 
    ! from Prigent et al. is 1e-3 cm.
    ! but that seems very low.
    ! For instance, the minimum value in the 
    ! measurements used in the regression
    ! in that study is 2.3e-3 cm.
    ! Also, at very low z0, volume scattering
    ! of the microwave radiation will take place 
    ! that can significantly decrease the radar 
    ! backscatter coefficient (p. 8).
    ! Furthermore, using 1e-3 cm leads to 
    ! an overestimation of AOD (compared to MODIS)
    ! in the areas concerned,
    ! in particular around the dust hot spots
    ! of the Sahara (using current u1fac value).
    ! For these reasons the minimum value
    ! has been increased. 
    !REAL(KIND=JPRB), PARAMETER           :: z0_min=1.e-3
    !REAL(KIND=JPRB), PARAMETER           :: z0_min=5.e-3
    REAL(KIND=JPRB), PARAMETER           :: z0_min=1.e-2
    !REAL(KIND=JPRB), PARAMETER           :: z0_min=2.e-2
    !<<< TvN

    REAL(KIND=JPRB), PARAMETER           :: lai_lim=0.25
    REAL(KIND=JPRB), PARAMETER           :: lai_lim2=0.5

    ! d_thrsld [cm^2.5] = 0.006/(ddust * grav*1.e2) with ddust = 2.65 g/cm^3,
    ! see Eq. (4) in MB95:
    REAL(KIND=JPRB), PARAMETER           :: d_thrsld=2.31e-6           ! threshold value
    !>>> TvN
    ! There are eight coarse-grained size bins,
    ! of which only the first four are used here.
    ! According to Tegen et al., Heinold et al.,
    ! the radius boundaries of the first seven bins are 
    ! at 0.1, 0.3, 0.9, 2.6, 8.0, 24, 72, and 220 um.
    ! However, these number don't seem to be exact.
    ! Since there is a constant ratio between the right
    ! and low boundaries, it seems this ratio is 3.0.
    ! Indeed, in Laurent et al. (JGR, 2010), 
    ! 2.6 is corrected to 2.7, which would be consistent
    ! with 8.0/3.0 = 2.67.
    ! This would imply that the radius boundaries are at
    ! 0.0987654 = 72./(3.^6), 0.296296, 0.889, 2.67, 8.0, 24, 72, 216, 
    ! and 648 um.

    ! Next, each bin is discretized with 24 size points,
    ! where d(n+1) = d(n) * exp(Dstep).
    ! Thus, Dstep = ln(3.)/24 = 0.04577551202.
    ! Dmin is the diameter of the first size point,
    ! given by 2* 72./(3.^6)) * exp(0.5*Dstep) = 0.20210403762 um.
    ! Similarly, the last size point is at a diameter
    ! 2* 648. * exp(-0.5*Dstep) = 1266.67434757 um.
    ! 
    ! With the original bin settings,
    ! the number of size points is 191 not 192 (=8*24).
    !
    !REAL(KIND=JPRB), PARAMETER           :: Dmin=0.00002                  ! minimum partic. diameter (cm)
    !REAL(KIND=JPRB), PARAMETER           :: Dmax=0.130                    ! maximum partic. diameter (cm)
    !REAL(KIND=JPRB), PARAMETER           :: Dstep=0.0460517018598807      ! diameter increment
    REAL(KIND=JPRB), PARAMETER           :: Dmin=2.0210403762e-5          ! diameter (cm) at first discretization point
    REAL(KIND=JPRB), PARAMETER           :: Dmax=0.126667434757           ! diameter (cm) at last discretization point
    REAL(KIND=JPRB), PARAMETER           :: Dstep=0.04577551202           ! diameter increment in log-space
    !<<< TvN

    ! Constants in the parameterization of the Reynolds number,
    ! see Eq. (5) in MB95:
    REAL(KIND=JPRB), PARAMETER           :: a_rnolds=1331.647             ! Reynolds constant
    REAL(KIND=JPRB), PARAMETER           :: b_rnolds=0.38194              ! Reynolds constant
    REAL(KIND=JPRB), PARAMETER           :: x_rnolds=1.561228             ! Reynolds constant
    !
    ! Air density has been made variable,
    ! to account for orographic effects.
    ! Previously, a global value for the 
    ! threshold friction velocity Uth was calculated.
    ! To keep its unit the same,
    ! roa is kept as a reference value,
    ! but its exact value is not important anymore.
    REAL(KIND=JPRB), PARAMETER           :: roa=0.001227                  ! reference air density (g/cm^3)
    REAL(KIND=JPRB), PARAMETER           :: airfac=1./8.3144*xmair*1.e-6    ! factor for rho_air
    !<<< TvN
    REAL(KIND=JPRB), PARAMETER           :: umin=13.75                    ! minimum threshold friction velocity (cm/s)
    REAL(KIND=JPRB), PARAMETER           :: ZZ=1000.                      ! wind measurement height (cm)

    ! parameters for the grouping in 2 modes
    ! The code follows the ECHAM-HAM implementation
    ! of Stier et al. (JGR, 2005),
    ! where the emission distribution is
    ! fitted onto three log-normal modes
    ! corresponding to the accumulation, coarse and super-coarse mode.
    ! (see presentation E. Vignati, TM meeting, 6 June 2008).
    !
    ! According to Heinold et al., 
    ! the three largest dust bins
    ! are less important for long-range transport,
    ! so particles with radius larger than 24 um
    ! can safely be neglected.
    ! However, a substantial part of the emitted mass
    ! is carried by particles with a radius larger than 10 um
    ! (see Tegen et al., Table 5).
    !
    ! The amounts of mass emitted in the accumulation and coarse modes
    ! are calculated from the masses emitted in the bin model,
    ! using two size ranges:
    ! r1 from 0.0987654 to 0.296296 um, and 
    ! r2 from 0.296296 to 8.0 um. 
    !
    ! Boundaries for Acc. mode
    INTEGER(KIND=JPIM), PARAMETER        :: min_ai=1
    INTEGER(KIND=JPIM), PARAMETER        :: max_ai=1
    ! Boundaries for Coa. mode
    INTEGER(KIND=JPIM), PARAMETER        :: min_ci=2
    INTEGER(KIND=JPIM), PARAMETER        :: max_ci=4
    !
    ! These size ranges include only part of
    ! the mass in the accumulation and coarse modes.
    ! The corresponding mass fractions are given by 
    ! mf(rmin,rmax) = 0.5*(
    !                erf(ln(rmax/mmr)/(sqrt(2)*ln(sigma)))-
    !                erf(ln(rmin/mmr)/(sqrt(2)*ln(sigma))) ),
    ! where mmr is the mass median radius.
    ! Applying this formula, 
    ! we find the following numbers:
    ! mf_acc(0,0.0987654)=0.00219913
    ! mf_acc_r1=mf_acc(0.0987654,0.296296)=0.313758
    ! mf_acc_r2=mf_acc(0.296296,8.0)=0.684043
    ! mf_acc(0.296296,inf)=0.684043
    !
    ! mf_coa(0,0.296296)=0.00519991
    ! mf_coa_r1=mf_coa(0.0987654,0.296296)=0.00518309 
    ! mf_coa_r2=mf_coa(0.296296,8.0)=0.980634
    ! mf_coa(8.0,inf)=0.0141665
    !
    REAL(KIND=JPRB), PARAMETER :: mf_acc_r1 = 0.313758
    REAL(KIND=JPRB), PARAMETER :: mf_acc_r2 = 0.684043
    REAL(KIND=JPRB), PARAMETER :: mf_coa_r1 = 0.00518309
    REAL(KIND=JPRB), PARAMETER :: mf_coa_r2 = 0.980634
    !
    ! Most importantly, r1 contains only about 31.4% 
    ! of the mass in the accumulation mode!
    ! This implies that we cannot just put the emissions
    ! from r1 to the accumulation mode, 
    ! and those from r2 to the coarse mode!
    !
    ! Instead, the modal emissions are determined
    ! by the following system of linear equations:
    ! mf_acc_r1 * flux_ai + mf_coa_r1 * flux_ci = flux_r1
    ! mf_acc_r2 * flux_ai + mf_coa_r2 * flux_ci = flux_r2,
    ! which relates the mass emitted in the ranges r1 and r2
    ! to the mass emitted in the accumulation and coarse modes.
    ! The solution is expressed using 
    ! the following parameters:
    !
    REAL(KIND=JPRB), PARAMETER :: ratio_coa = mf_coa_r1/mf_coa_r2
    REAL(KIND=JPRB), PARAMETER :: ratio_acc = mf_acc_r2/mf_acc_r1
    REAL(KIND=JPRB), PARAMETER :: denom_acc_inv = 1./(mf_acc_r1-ratio_coa*mf_acc_r2)
    REAL(KIND=JPRB), PARAMETER :: denom_coa_inv = 1./(mf_coa_r2-ratio_acc*mf_coa_r1)
    REAL(KIND=JPRB), PARAMETER :: mf_acc_r12_inv = 1./(mf_acc_r1+mf_acc_r2)
    REAL(KIND=JPRB), PARAMETER :: mf_coa_r12_inv = 1./(mf_coa_r1+mf_coa_r2) 
    ! 
    ! Source mass median radius (cm)
    ! Stier et al. (2005) uses very similar numbers 
    ! for mass median radii, 
    ! but uses 0.37 um for the accumulation mode.
    ! Thus, it seems these numbers are not mass mean,
    ! but mass median radii.
    !
    ! The super-coarse mode has 
    ! a mass median radius of 15.0 and sigma=2.0,
    ! but is not included.
    !
    ! The AeroCom recommendation of Dentener et al. (ACP, 2006)
    ! is to use a number median radius
    ! of 0.65 um for the coarse mode,
    ! which corresponds to mass median radius of 2.75 um
    ! (the conversion factor is exp(3.0*ln(sigma)^2),
    ! see Zender, Particle Size Distributions:
    ! Theory and Application to Aerosols, Clouds, and Soils, 2002).
    ! 
    !REAL(KIND=JPRB), PARAMETER           :: mmr_ai=0.35E-4
    REAL(KIND=JPRB), PARAMETER           :: mmr_ai=0.37E-4
    REAL(KIND=JPRB), PARAMETER           :: mmr_ci=1.75E-4
    !<<< TvN

    !----------------------------------------------------------------
    ! SOIL CARACTERISTICS:
    ! ZOBLER texture classes 
    !----------------------------------------------------------------

  ! solspe includes for each soil type (first dimension)
  ! the mass median diameter (cm) and standard deviation (see Table 1, MB95)
  ! and the relative contribution (Table 2, MP95) for the four size populations.
  ! The two additional entries describe the saltation efficiency alpha (cm^-1),
  ! and the residual moisture, which is currently not used.
  ! Efficiencies are calculated as averages over the four populations
  ! (as in Eq. (8) in Marticorena et al. (JGR, 1997),
  ! where 1e-7, 1e-6 and 1e-5 cm^-1 is used for coarse sand, 
  ! medium/fine sand and silt, respectively,
  ! and 1e-6 for clay for soils with clay fractions below 45%
  ! and 1e-7 for clay for soils with clay fractions above 45%.
  ! (Tegen et al.).
  REAL(KIND=JPRB), PARAMETER, DIMENSION(nats,nspe) :: solspe=RESHAPE( (/ &
    !--     soil type 1 : Coarse
         0.0707, 2.,  0.43 ,      &
         0.0158, 2.,  0.4 ,       &
         0.0015, 2.,  0.17 ,      &
         0.0002 ,2.,  0. ,        &
         2.1E-06,   0.2, &
    !--     soil type 2 : Medium
         0.0707, 2.,  0. ,            &
         0.0158, 2.,  0.37 ,          &
         0.0015, 2.,  0.33 ,          &
         0.0002, 2.,  0.3 ,           &
         4.0e-6,    0.25, &
    !--     soil type 3 : Fine
         0.0707, 2.,  0. ,            &
         0.0158, 2.,  0. ,            &
         0.0015, 2.,  0.33 ,          &
         0.0002, 2.,  0.67 ,          &
         !>>> TvN
         ! 33% x 1e-5 + 67% x 1e-7 = 3.367e-6 cm^-1
         !1.E-07,   0.5, &
         3.4e-6,   0.5, &
         !<<< TvN
    !--     soil type 4 : Coarse Medium
         0.0707, 2.,  0.1 ,           &
         0.0158, 2.,  0.5 ,           &
         0.0015, 2.,  0.2 ,           &
         0.0002, 2.,  0.2 ,           &
         2.7E-06,   0.23, &
    !--     soil type 5 : Coarse Fine
         0.0707, 2.,  0. ,            &
         0.0158, 2.,  0.5 ,           &
         0.0015, 2.,  0.12 ,          &
         0.0002, 2.,  0.38 ,          &
         !>>> TvN
         ! 50% x 1e-6 + 12% x 1e-5 + 38% x 1e-6 = 2.08e-6 cm^-1
         !2.8E-06,   0.25, &
         2.1e-6,   0.25, &
         !<<< TvN
    !--     soil type 6 : Medium Fine
         0.0707, 2.,  0.   ,          &
         0.0158, 2.,  0.27 ,          &
         0.0015, 2.,  0.25 ,          &
         0.0002, 2.,  0.48 ,          &
         !>>> TvN
         ! 27% x 1e-6 + 25% x 1e-5 + 48% x 1e-7 = 2.818e-6 cm^-1
         !1e-07,   0.36, &
         2.8e-6,   0.36, &
         !<<< TvN
    !--     soil type 7 : Coarse, Medium, Fine
         0.0707, 2.,  0.23 ,          &
         0.0158, 2.,  0.23 ,          &
         0.0015, 2.,  0.19 ,          &
         0.0002, 2.,  0.35 ,          &
         2.5E-06,  0.25, &
    !--     soil type 8 : Organic
         0.0707, 2.,  0.25 ,          &
         0.0158, 2.,  0.25 ,          &
         0.0015, 2.,  0.25 ,          &
         0.0002, 2.,  0.25 ,          &
         0.,   0.5, &
    !--     soil type 9 : Ice
         0.0707,  2.,  0.25 ,         &
         0.0158,  2.,  0.25 ,         &
         0.0015,  2.,  0.25 ,         &
         0.0002,  2.,  0.25 ,         &
         0.,       0.5, &
    !--     soil type 10 : Potential Lakes (additional)
    !       GENERAL CASE
         0.0707,  2.,  0. ,            &
         0.0158,  2.,  0. ,            &
         0.0015,  2.,  1. ,            &
         0.0002,  2.,  0. ,            &
         1.E-05,  0.25, &
    !--     soil type 11 : Potential Lakes (clay)
    !       GENERAL CASE
         0.0707,  2.,  0. ,            &
         0.0158,  2.,  0. ,            &
         0.0015,  2.,  0. ,            &
         0.0002,  2.,  1. ,            &
         1.E-05,  0.25, &
    !--     soil type 12 : Potential Lakes Australia
         0.0707,  2.,  0. ,            &
         0.0158,  2.,  0. ,            &
         0.0027,  2.,  1. ,            &
         0.0002,  2.,  0. ,            &
         1.E-05,  0.25 /),(/nats,nspe/),order=(/2,1/) )   

   REAL(KIND=JPRB) ::    UTH  (     NCLASS)
   REAL(KIND=JPRB) ::    SREL (NATS,NCLASS)
   REAL(KIND=JPRB) ::    SRELV(NATS,NCLASS)
   REAL(KIND=JPRB) :: SU_SRELV(NATS,NCLASS)
   
    

  
END MODULE TM5M7_EMIS_DATA
