SUBROUTINE TM5M7_SCAV &
 & ( KIDIA , KFDIA  , KLON , KLEV , KAER , KSTEP, PTSPHY, &
 &   PRSF1, PDP , PTP, PFLXR, PFLXS, PCLCOV, PCLWAT, PCLICE, PRAIN, PSNOW, PCEN, PTENC0, &
 &   PTENC1, PFAERO, PPRCOV  )

!*** * TM5M7_SCAV* - IN-CLOUD AND BELOW CLOUD SCAVENGING OF TRACERS
!      CALLED SEPARATELY FOR CONVECTIVE AND LARGE-SCALE PRECIP. 
! COPY of CHEM_SCAV 

! INPUTS:
! -------
! KSTEP :  Time step number
! KIDIA :  Start of Array
! KFDIA :  End  of Array
! KLON  :  Length of Arrays
! KLEV  :  Number of Levels
! KCHEM :  Number of aerosol tracers 


! PTSPHY:  Time step length in seconds
! PDP(KLON,KLEV)              :  PRESSURE DELTA in PRESSURE UNITES      (Pa)
! PRSF1(KLON,KLEV)            :  Mid-level pressure           (Pa)
! PTP(KLON,KLEV)              :  Temperature in            (T)
! PCLWAT  (KLON,KLEV)         :  Cloud water content    (kg/kg) - stratiform and convective 
! PCLICE  (KLON,KLEV)         :  Cloud ice water content    (kg/kg) - stratiform and convective
! PRAIN (KLON,KLEV)           :  Rain water content    (kg/kg) for stratiform precip
! PSNOW (KLON,KLEV)           :  Snow water content    (kg/kg) for stratiform precip
! PFLXR   (KLON,KLEV+1)       :  Precip flux rain       (kg/m2s) - either stratiform or convective (depending on argument)
! PFLXS   (KLON,KLEV+1)       :  Precip flux snow       (kg/m2s) - either stratoform or convective (depending on argument)
! PCLCOV  (KLON,KLEV)         :  cloud fraction   0..1
! PPRCOV  (KLON,KLEV)         :  precipitation fraction   0..1
! PCEN(KLON,KLEV,KCHEM)       :  CONCENTRATION OF TRACERS           (kg/kg)
! PTENC0(KLON,KLEV,KCHEM)     :  TOTAL TENDENCY OF CONCENTRATION OF TRACERS BEFORE(kg/kg s-1)
!
! NB: PCLWAT is the in-cloud water mixing ratio
! OUTPUTS:
! -------
! PTENC1 (KLON,KLEV,KCHEM)     : TENDENCY OF CONCENTRATION OF TRACERS after (kg/kg s-1)
!
!**   INTERFACE.
!     ----------
!          *TM5M7_SCAV* IS CALLED FROM *TM5M7*.
!
!     AUTHOR.
!     -------
!        Johannes Flemming (original CHEM_SCAV)
!        
!
!     MODIFICATIONS.
!     --------------
!        ORIGINAL (J Flemming) : 2009-11-09
!        V. Huijnen : Modification following routine wet_deposition.F90 from TM5 2021-08-26
!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE YOMCST    ,ONLY : RG, RD, RLVTT ,RLSTT ,RTT 
USE YOMLUN   , ONLY : NULOUT,NULERR
USE YOEDBUG   ,ONLY : YREDBUG
USE YOM_YGFL , ONLY : YGFL
USE YOECLDP  , ONLY : YRECLDP 
USE YOEAERSRC ,ONLY : YREAERSRC
USE YOEAERATM ,ONLY : YREAERATM
USE YOEAERSNK ,ONLY : YREAERSNK
USE YOMCHEM  , ONLY : YRCHEM 
USE YOETHF   , ONLY :  R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
 & R4IES    ,R5LES    ,R5IES    ,R5ALVCP  ,R5ALSCP  ,&
 & RALVDCP  ,RALSDCP  ,RTWAT    ,&
 & RTICE    ,RTICECU  ,&
 & RTWAT_RTICE_R      ,RTWAT_RTICECU_R
USE TM5M7_DATA, ONLY : NSCAV, NSCAV_INDEX, NSCAV_TYPE, &
 & INUS_N,IAIS_N,IACS_N,ICOS_N,IAII_N,IACI_N,ICOI_N


IMPLICIT NONE

!-----------------------------------------------------------------------

!*       0.1  ARGUMENTS
!             ---------

INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA, KFDIA, KLON, KLEV, KAER,  KSTEP

REAL(KIND=JPRB),INTENT(IN)    :: PDP(KLON,KLEV) , PRSF1(KLON,KLEV) , PTP(KLON,KLEV)   
REAL(KIND=JPRB),INTENT(IN)    :: PCLCOV(KLON,KLEV) , PCLWAT(KLON,KLEV),PCLICE(KLON,KLEV) , PRAIN(KLON,KLEV), PSNOW(KLON,KLEV)  
REAL(KIND=JPRB),OPTIONAL, INTENT(IN)    ::  PPRCOV(KLON,KLEV)  
REAL(KIND=JPRB),INTENT(IN)    :: PFLXR(KLON,KLEV+1), PFLXS(KLON,KLEV+1)
REAL(KIND=JPRB),INTENT(IN)    :: PTENC0(KLON,KLEV,KAER), PCEN(KLON,KLEV,KAER)
REAL(KIND=JPRB),INTENT(IN)    :: PTSPHY

REAL(KIND=JPRB),INTENT(OUT)   :: PTENC1(KLON,KLEV,KAER)
REAL(KIND=JPRB),INTENT(OUT)   :: PFAERO(KLON,KAER)

!*       0.5   LOCAL VARIABLES
!              ---------------

INTEGER(KIND=JPIM) :: JK, JL, JT, JSCAV, IWETDEP

REAL(KIND=JPRB) :: ZSCAV
REAL(KIND=JPRB) :: ZBETA,  ZBETAR, ZBETAR_M7
REAL(KIND=JPRB) :: INCLOUD, BELOWCLOUD
REAL(KIND=JPRB) :: ZDP, ZDX, ZDZ, ZFRAC, ZFUNC, ZRHO,ZTEMP
REAL(KIND=JPRB) :: ZEPSQLIQ, ZEPSFLX, ZEPSQLIQ2
REAL(KIND=JPRB) :: ZRD, ZLMMR2VMR, ZH2R, ZRET, ZHNRYEFT, ZHNRYEF    
REAL(KIND=JPRB) :: ZLIQ2TOT, ZICE2TOT, ZLIQ2GAS, ZICE2GAS, ZRAINW, ZFALLSP
REAL(KIND=JPRB) :: ZPRCOV, ZCLCOV, ZFLXR, ZFLXS,ZFLXRB, ZFLXSB,  ZCLWAT, ZCLICE,ZCLTOT, ZMINCLCOV, ZMINPRCOV 
 ! Interstitial Fraction: 30% of aerosol remains in atmosphere
REAL(KIND=JPRB)  ::  ZINTERST_FR  
REAL(KIND=JPRB),PARAMETER  :: ZDGHNO3 = 0.136      ! viscosity of HNO3 in [cm2/s] 
REAL(KIND=JPRB),PARAMETER  :: ZDGAIR  = 0.133      ! viscosity of air in [cm2/s] 
REAL(KIND=JPRB),PARAMETER  :: ZXMHNO3    =1.008_JPRB + 14.007_JPRB + 3*16.0_JPRB ! HNO3 tracer mass
REAL(KIND=JPRB),PARAMETER :: ZRDRAD2 = (1E-5)**2  ! square of raindroplet radius (20 microns)
REAL(KIND=JPRB),PARAMETER :: ZMAX_LWC   =2.E-3    ! kg/m3 
REAL(KIND=JPRB),PARAMETER :: ZHPLUS =3.16227E-6_JPRB ! is rain water ph=5.5 H+ 

! how much less efficient is tracer scavenged from ice
! cloud droplet compared to water cloud droplet. 
! This should be tracer dependent. 
REAL(KIND=JPRB),PARAMETER :: ZICE_EFF=0.2

REAL(KIND=JPRB)  :: ZRLWC,ZRDRAD, ZRU, ZNRE,ZNSC,ZNSH,  ZKG, ZTR, ZKSO2, ZKHSO3, ZFACTSO2
REAL(KIND=JPRB)            :: ZRL      ! composite factor of Rgas and liquid water content of raining cloud
                                       ! rgas (8.314 J/mol/K) ---> 0.08314 atm/(mol/l)/K
                                       ! 1e-6 corresponds to 1 g/m3 dimensionless 
REAL(KIND=JPRB)   :: ZWLOSS(KLON,KAER)

LOGICAL :: LLWDAER 
LOGICAL :: LLPRINT, LLCHEM_WDFR, LLCONV

REAL(KIND=JPRB) :: ZHOOK_HANDLE

#include "fcttre.func.h"
!#include "fccld.h"


!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TM5M7_SCAV',0,ZHOOK_HANDLE)
ASSOCIATE(YCHEM=>YGFL%YCHEM,YAERO=>YGFL%YAERO, &
 & LCHEM_WDFR=>YRCHEM%LCHEM_WDFR, &
 & RCOVPMIN=>YRECLDP%RCOVPMIN, &
 & NINDAER=>YREAERSRC%NINDAER, &
 & RFRAER=>YREAERSNK%RFRAER, RFRBC=>YREAERSNK%RFRBC, RFRDD=>YREAERSNK%RFRDD, &
 & RFRIF=>YREAERSNK%RFRIF, RFROM=>YREAERSNK%RFROM, RFRSO4=>YREAERSNK%RFRSO4, &
 & RFRSS=>YREAERSNK%RFRSS)



LLPRINT=.FALSE.
LLCHEM_WDFR=LCHEM_WDFR
LLCONV=.FALSE.

ZINTERST_FR = 0.3_JPRB
!VH IF ( LLCHEM_WDFR ) ZINTERST_FR = 0.0_JPRB 
  
! deposition for convective Precip
IF (.NOT. PRESENT(PPRCOV) )  LLCONV=.TRUE.

!* set flux to zero 
ZWLOSS(:,:)=0._JPRB

ZEPSFLX =1.E-18_JPRB
ZEPSQLIQ=2.E-6_JPRB
ZEPSQLIQ2=2.E-7_JPRB

!* mini cloud cover and precip cover
ZMINCLCOV=0.001_JPRB
ZMINPRCOV=0.001_JPRB

!* precip fall speed 
ZFALLSP = 5.0_JPRB

!* ideal gas constant in atm M-1
ZRD=1000.0_JPRB * RD * 9.8692_JPRB / 1000000.0_JPRB
ZRD=0.082_JPRB
ZRL=RG/1E2_JPRB*1E-6_JPRB 
PFAERO(:,:)=0._JPRB

!* initialisation    
!* re-evaporation fraction to account for drop shrinking without releasing species
! Jacob says 
!ZFRAC=0.5_JPRB
ZFRAC=0.2_JPRB

!* update tendecies  

!* if not LLCHEM_WDFR 
ZPRCOV =  1.0_JPRB
ZCLCOV =  1.0_JPRB

PTENC1(KIDIA:KFDIA,1:KLEV,1:KAER) =PTENC0(KIDIA:KFDIA,1:KLEV,1:KAER)

!* LOOP OVER LAYERS  
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA

!- temperature
    ZTEMP = PTP(JL,JK)
!- precip flux at box top 
    ZFLXR =  PFLXR(JL,JK)
    ZFLXS =  PFLXS(JL,JK)
    ZFLXRB =  PFLXR(JL,JK+1)
    ZFLXSB =  PFLXS(JL,JK+1)

!--  Precip flux change grid box average - indicator for rain/snow formation or evaporation
    ZBETA = ( ZFLXRB - ZFLXR ) + ( ZFLXSB - ZFLXS  )  ! in kg m-2 s-1 
    ZCLWAT =  PCLWAT(JL,JK) 
    ZCLICE =  PCLICE(JL,JK) 

    IF (  ZBETA ==  0.0_JPRB .AND. ZFLXR < ZEPSFLX  )  CYCLE 
    IF ( LLCHEM_WDFR ) THEN 
!* use effective flux and cloud cover for liquid water content
      IF (.NOT. LLCONV ) THEN
        ZPRCOV = MIN(MAX(RCOVPMIN,PPRCOV(JL,JK)),1.0_JPRB)  ! grid scale precip
        ZCLCOV = MIN(1.0_JPRB,MAX(0.0_JPRB,PCLCOV(JL,JK)))
      ELSE
        ZPRCOV = 0.05    ! convective precip  
        ZCLCOV = 0.05    ! convective precip  
       ! ZPRCOV = 0.1    ! convective precip  
      ENDIF  
      IF (ZPRCOV >= ZMINPRCOV) THEN
        ZFLXR =  ZFLXR / ZPRCOV 
        ZFLXS =  ZFLXS / ZPRCOV  
        ZFLXRB =  ZFLXRB / ZPRCOV 
        ZFLXSB =  ZFLXSB / ZPRCOV  
        ZBETA = ZBETA / ZPRCOV
      ELSE
        ZFLXR =  0.0_JPRB 
        ZFLXS =  0.0_JPRB  
        ZFLXRB =  0.0_JPRB 
        ZFLXSB =  0.0_JPRB  
        ZBETA =  0.0_JPRB
      ENDIF
      IF (ZCLCOV >= ZMINCLCOV ) THEN 
        ZCLWAT =  MAX(0.0_JPRB,PCLWAT(JL,JK) /  ZCLCOV)
        ZCLICE =  MAX(0.0_JPRB,PCLICE(JL,JK) /  ZCLCOV)
      ELSE
        ZCLWAT = 0.0_JPRB  
        ZCLICE = 0.0_JPRB  
      ENDIF 
    ENDIF

!* calculate air density, layer depth 
     ZRHO=PRSF1(JL,JK)/(RD*ZTEMP)
     ZDZ= PDP(JL,JK) / (ZRHO*RG)
     ZDP =  PDP(JL,JK)
! clwc mass mixing ratio to volume mixing ratio : rho air / rho cloud water - note it is liquid water not vapor
    ZLMMR2VMR=ZRHO/1000.0_JPRB 
    IWETDEP=0

!* LOOP over species
    DO JSCAV=1,NSCAV
       
       ! Identify aerosol tracer
       JT=NSCAV_INDEX(JSCAV)
       LLWDAER=.TRUE. 
       IWETDEP=IWETDEP+1

       !VH This whole 'selection mechanism is ideally coded in an alternative way,
       !   using a lookup table approach. That should be comuptationally much more efficient.
       SELECT CASE (NSCAV_TYPE(JSCAV))
       case(0)
          incloud    = 0.0
          belowcloud = 0.0
       case(1)   ! 100% solubility
          !incloud    = 1.
          !belowcloud = corr_diff
	  CALL ABOR1(" option not supported")
       case(2)   ! henry solubility assumed
          ! rtl = rtl*henry(n,itemp) / ( 1.0 + rtl*henry(n,itemp) )
          ! incloud    = rloss1(region)%d3(i,j,k)*rtl
          ! belowcloud = rloss2(region)%d3(i,j,k)*rtl*corr_diff
          CALL ABOR1("Henry solubility currently not supported")
       case(3)   ! bulk aerosol 
          incloud    = (1.0 - ZINTERST_FR)
          !>>>TvN
          ! Alternative would be to make the interstitial fraction for bulk aerosols
          ! consistent the values used for the M7 modes, 
          ! which are taken from Bourgeois and Bey (JGR, 2011)
          ! and distinguish between warm, mixed and ice clouds
          !<<<TvN
          belowcloud = 1.
       case(4)   ! SO2
          ! ztr=(1./ZTEMP-1./298)
          ! dkso2 =1.7e-2*exp(2090.*ztr)	!so2<=>hso3m+hplus
          ! dkhso3  = 6.6e-8*exp(1510.*ztr)	!hso3m<=>so3-- + hplus
          ! factor = 1.0 + dkso2/hplus + (dkso2*dkhso3)/(hplus**2)
          ! heff = factor*henry(n,itemp)
          ! rtl = rtl*heff/ ( 1.0 + rtl*heff )
          ! incloud    = rloss1(region)%d3(i,j,k)*rtl !
          ! belowcloud = rloss2(region)%d3(i,j,k)*rtl*corr_diff
          CALL ABOR1("SO2 solubility currently not supported")


!>>>TvN
! The in-cloud scavenging coefficients are defined as the fraction of the tracer 
! in the cloudy part of the grid box that is embedded in the cloud liquid or ice water,
! i.e. the non-interstitial part.
! We distinguish between liquid, mixed and ice stratiform clouds (Stier et al., 2005),
! depending on the local temperature in the grid cell (Croft et al., ACP, 2010).
! The in-cloud scavenging coefficients depend on size and composition;
! revised values for the M7 modes were provided by Bourgeois and Bey (JGR, 2011).
! For mixed clouds, an alternative method was presented by Zhang et al. (ACP, 2012), 
! which uses a continuous temperature dependency.
! Note that these in-cloud scavenging coefficients account for both nucleation scavenging
! and impaction scavenging (Croft et al., ACP, 2009; 2010).
! Thus, the below-cloud scavenging rates should only account for
! the impaction scavenging by precipitation coming from clouds above the current level.
!
! Estimates for below-cloud scavenging coefficients can be derived 
! from Fig. 2 of Dana and Hales (AE, 1975).
! For estimating these values from the figure, I used aerodynamic radii of 
! 0.007, 0.07, and 0.7 micron as the boundaries of the M7 modes
! (corresponding to a particle density of about 1800 g/cm^3).
! As in Stier et al. (2005), we do not distinguish between soluble and insoluble modes.
! Thus, dry particle radii can be used for estimating the scavenging coefficients from the figure
! (see also the mode boundaries applied in Fig. 2 in Croft et al., 2009).
! I thus arrive at the following rough estimates for below-cloud mass scavenging coefficients
! for the nucleation, aitken, accumulation, and coarse modes: ~0.01, 0.002, 0.01, and 1 mm^-1.
! These numbers are close to the estimates derived earlier from the same figure
! by Elisabetta Vignati, which were previously used: 0.005, 0.002, 0.008, and 1 mm^-1.
!
! However, both sets of estimates based on Dana and Hales are substantially higher
! than the values presented by Croft et al. (2009).
! From the curves presented in their Fig. 2 for the standard Marshall-Palmer rain distribution,
! rough estimates of the mass scavenging coefficients for the four size modes can be obtained.
! My estimates are 0.002, 0.0002, 0.03, and 0.7 mm^-1. 
! Note that especially the value for the accumulation mode is very sensitive to the
! actual mean particle size, and hard to estimate from the figure.
! Since the mean droplet size of the Marshall-Palmer distribution depends on the rain intensity,
! these estimates are only valid for a rain rate of 1 mm/hr.
! For simplicity, we assume that the scavenging coefficients derived from the figure at 1 mm/hr
! can also be applied at other rain intensities.
! 
! In the new implementation particle masses and numbers are scavenged at different rates.
! Rough estimates of the number scavenging coefficients for the four size modes
! can be obtained from the same figure in Croft et al.
! My estimates are 0.02, 0.001, 0.0003, and 0.3 mm^-1.

! Ideally, the below-cloud mass/number scavenging coefficients should be calculated
! using look-up tables to describe the dependence on median radius and precipitation rate, 
! e.g. following the formulation/curves presented by Croft et al.
!


         case(5)   ! soluble nu
            if (JT /= inus_n) then
              belowcloud=0.5 ! 0.5*0.004 = 0.002 mm^-2
            else
               belowcloud=5. ! 5.*0.004 = 0.02 mm^-2
           endif
            incloud=0.06

         case(6)   ! soluble ai
            if (jt /= iais_n) then
              belowcloud=0.05 ! 0.05*0.004 = 0.0002 mm^-2
            else
               belowcloud=0.25 ! 0.25*0.004 = 0.001 mm^-2
           endif
            if (ZTEMP.gt.273.15) then 
               incloud=0.25
            else
               incloud=0.06
            endif

         case(7)   ! soluble ac
            if (jt /= iacs_n) then
              belowcloud=7.5  ! 7.5*0.004 = 0.03 mm^-1
            else  
               belowcloud=0.075  ! 0.075*0.004 = 0.0003 mm^-1
           endif

            if (ZTEMP.gt.273.15) then 
               incloud=0.85
            else
               incloud=0.06
            endif 

         case(8)   ! soluble co
            if (jt /= icos_n) then
              belowcloud=175. ! 175*0.004 = 0.7 mm^-1
            else
               belowcloud=75. ! 75*0.004 = 0.3 mm^-1
           endif
            if (ZTEMP.gt.273.15) then
               incloud=0.99
            else if (ZTEMP.gt.238.15) then
               incloud=0.75
            else
               incloud=0.06
            endif

         case(9)   ! insoluble ai
            if (jt /= iaii_n) then
              belowcloud=0.05
            else
               belowcloud=0.25
            endif
            !incloud=0.0
            if (ZTEMP.gt.273.15) then
               incloud=0.2
            else
               incloud=0.06
            endif
            
         case(10)   ! insoluble ac
            if (jt /= iaci_n) then
              belowcloud=7.5
            else
               belowcloud=0.075
           endif
            if (ZTEMP.gt.273.15) then
               incloud=0.4
            else
               incloud=0.06
            endif
            
         case(11)   ! insoluble co
            if (jt /= icoi_n) then
              belowcloud=175.
            else
               belowcloud=75.
           endif
            if (ZTEMP.gt.238.15) then
               incloud=0.4
            else
               incloud=0.06
            endif 

         case default
            incloud    = 0.0
            belowcloud = 0.0
         end select





         !* total cloud water and ice
         ZCLTOT=ZCLWAT + ZCLICE 
         !* Rain-out in Cloud
         IF (ZCLTOT > ZEPSQLIQ .AND. ZBETA > 0.0_JPRB ) THEN
            !* water loss due to precip formation relative to total (water and ice) cloud
            ZBETAR_M7=ZBETA / ( ZCLTOT*ZDP / RG ) ! in s-1
            !* safety check
            ZBETAR_M7=MIN( 200._JPRB, MAX( ZBETAR_M7, 0._JPRB ))

            !* water loss due to precip formation relative to cloud, with reduced weight for ice cloud.
	    !* Relevant to trace gas scavenging, but not treated here.
            ! ZINC_RDF=ZCLICE/ZCLTOT*ICE_EFF + ZCLWAT/ZCLTOT
            ! ZBETAR=ZBETA *ZINC_RDF / ( ZCLTOT*ZDP / RG ) ! in s-1
            !* safety check
            ! ZBETAR=MIN( 200._JPRB, MAX( ZBETAR, 0._JPRB ))


            !* scavenging coefficient  in s-1 
            !In-cloud scavenging is different for aerosol than for gas-phase 
            ZSCAV=ZBETAR_M7 *INCLOUD

            ZFUNC=EXP(-ZSCAV*PTSPHY)                     ! N.D.       (N.D.)
            IF (.NOT. LLCONV) THEN
              ZFUNC=MAX(0.95_JPRB, ZFUNC)
            ENDIF
            ZDX = PCEN(JL,JK,JT)*(ZFUNC - 1._JPRB) * ZPRCOV      ! in kg kg-1
            PTENC1 (JL,JK,JT) = PTENC1(JL,JK,JT) + ZDX/PTSPHY                ! in kg kg-1 s-1
            ZWLOSS(JL,IWETDEP) = ZWLOSS(JL,IWETDEP) - ZDX/PTSPHY * ZDP/RG             ! in kg m-2 s-1
         ENDIF





!* reevaporation  ! bottom flux smaller and top flux above limit
      IF (ZBETA < 0.0_JPRB .AND. ZFLXR+ZFLXS > ZEPSFLX ) THEN
        ZBETAR=ZBETA/(ZFLXR+ZFLXS)
        IF ( ZFLXRB+ZFLXSB  <= ZEPSFLX) THEN
!-- total reevaporation, bottom flux <= ZEPSFLX
         ZBETAR=MIN(MAX(0._JPRB, -ZBETAR), 1._JPRB)              ! N.D.
       ELSE
!-- partial reevaporation 
         ZBETAR=MIN(MAX(0._JPRB, -ZBETAR)*ZFRAC, 1._JPRB)        ! N.D.
       ENDIF
       ZDX = (ZBETAR*ZWLOSS(JL,IWETDEP) * (RG/ZDP) *PTSPHY) * ZPRCOV 
       PTENC1(JL,JK,JT) = PTENC1(JL,JK,JT) + ZDX/PTSPHY
       ZWLOSS(JL,IWETDEP) = ZWLOSS(JL,IWETDEP) - ZDX/PTSPHY * ZDP/RG             ! in kg m-2 s-1
       IF (LLPRINT .AND. ZBETAR > 0._JPRB) THEN
         WRITE(UNIT=NULOUT,FMT='(1x,''ZSCAV3'',5I3,(11E10.3))') KSTEP,JL,JK,&
&        PTENC0(JL,JK,JT),PCEN(JL,JK,JT),PTENC1(JL,JK,JT),ZWLOSS(JL,IWETDEP),ZBETAR,ZDX
       ENDIF
     ENDIF




!!*   wash-out with rain  when flux is positive and no rain formation or evaporation (below-cloud)
!    IF (ZFLXR > ZEPSFLX .AND. ABS(ZBETA) <= ZEPSFLX ) THEN    
!*   wash-out with rain  when flux is positive  
    IF (ZFLXR > ZEPSFLX ) THEN    

         ! calculate fraction in rainwater according to Henry 
         !* rain water in box in kg/kg
         IF (LLCONV ) THEN 
           ZRAINW=(ZFLXR/(ZFALLSP*ZRHO)) / ZPRCOV  
         ELSE        
           ZRAINW = (PRAIN(JL,JK)+PSNOW(JL,JK)) / ZPRCOV ! + PSNOW(JL,JK) ! use prognostic variable from input 
         ENDIF

           ! Rain
           ZSCAV=ZFLXR * BELOWCLOUD  
           ZFUNC=EXP(-ZSCAV*PTSPHY)                          ! N.D.       (N.D.)
           ZDX = PCEN(JL,JK,JT)*(ZFUNC -1._JPRB)*ZPRCOV      ! in kg kg-1
           ZFUNC=MAX(0.97_JPRB, ZFUNC)

           ! snow
           ZSCAV=ZFLXS*5.E-3                                      ! first assumption - requires review/update!
           ZFUNC=EXP(-ZSCAV*PTSPHY)                               ! N.D.       (N.D.)
           ZDX = ZDX + PCEN(JL,JK,JT)*(ZFUNC - 1._JPRB)*ZPRCOV    ! in kg kg-1
           ZFUNC=MAX(0.97_JPRB, ZFUNC)

        PTENC1 (JL,JK,JT) = PTENC1(JL,JK,JT) + ZDX/PTSPHY                ! in kg kg-1 s-1
        ZWLOSS(JL,IWETDEP) = ZWLOSS(JL,IWETDEP) - ZDX/PTSPHY * ZDP/RG             ! in kg m-2 s-1

     ENDIF  ! rainout
    ENDDO
  ENDDO
ENDDO



!* LOOP OVER LAYERS  
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    PFAERO(JL,:) = PFAERO(JL,:) -(PTENC1(JL,JK,:)-PTENC0(JL,JK,:))*(PDP(JL,JK))/RG
  ENDDO
ENDDO


!! debug undo anything  
! PTENC1(KIDIA:KFDIA,1:KLEV,1:KCHEM) =PTENC0(KIDIA:KFDIA,1:KLEV,1:KCHEM)
!-----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('TM5M7_SCAV',1,ZHOOK_HANDLE)
END SUBROUTINE TM5M7_SCAV 
