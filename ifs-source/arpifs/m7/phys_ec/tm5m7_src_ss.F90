SUBROUTINE TM5M7_SRC_SS &
  &( KIDIA, KFDIA, KLON, KLEV, &
  &  PCI  , PCLAKE, PLSM , PSST, PWIND, &
  &  emis_mass, emis_number & 
  &)

!*** * TM5M7_SRC_SS* - SOURCE TERMS FOR SEA SALT AEROSOLS

!**   INTERFACE.
!     ----------
!          *TM5M7_SRC_SS* IS CALLED FROM *TM5M7_SRC*.

!     AUTHOR.
!     -------
!
!     T. van Noije et al.
!
!     SOURCE.
!     -------
! !DESCRIPTION: Calculate emitted numbers and mass as function of the ten-meter
!                wind speed as described in Vignati et al. (2010) and Gong
!                (2003). Sea salt is emitted in both the soluble accumulation
!                mode and the soluble coarse mode.
!              
!          Ref: Vignati et al. (2010) : Global scale emission and
!                distribution of sea-spray aerosol: Sea-salt and organic
!                enrichment, Atmos. Environ., 44, 670-677,
!                doi:10.1016/j.atmosenv.2009.11.013
!
!               Gong (2003) : A parameterization of sea-salt aerosol source
!                function for sub- and super-micron particles, Global
!                Biogeochem. Cy., 17, 1097, doi:10.1029/2003GB002079
  
!     MODIFICATIONS.
!     --------------
!
!    ? ??? 2006 - EV and MK - changed for introducing the sea salt 
!                 source function developed from Gong (2003) in two modes 
!    ? ??? ???? - AJS - switch from default aerocom emissions to Gong 
!                 parameterisation if 'seasalt_emis_gong' is defined.
!    1 Sep 2010 - Achim Strunk - deleted procedures
!                              - removed with_seasalt-switch
!                              - introduced vertical splitting 
!   25 Jun 2012 - Ph. Le Sager - adapted for lon-lat MPI domain decomposition
!    April 2015 - T. van Noije - modified mode prefactors
!      Sep 2020 - V. Huijnen: first (partial) introduction into OpenIFS
!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST, ONLY : RPI
USE TM5M7_DATA, ONLY: NMOD, MODE_ACS, MODE_COS, sigma_lognormal, SS_DENSITY 
USE TM5M7_EMIS_DATA, ONLY : MODAL_EMISSIONS, radius_ssa, radius_ssc

USE MO_HAM_M7_EMI_SEASALT, ONLY: seasalt_emissions_gong_SST
USE MO_HAM,          ONLY: nseasalt

IMPLICIT NONE

!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------
INTEGER(KIND=JPIM),    INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),    INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),    INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),    INTENT(IN)    :: KLEV
REAL(KIND=JPRB),       INTENT(IN)    :: PCI(KLON), PCLAKE(KLON), PLSM(KLON), PSST(KLON)
REAL(KIND=JPRB),       INTENT(IN)    :: PWIND(KLON)  ! 10m wind speed, see tm5m7_src.F90
TYPE(MODAL_EMISSIONS), INTENT(INOUT) :: emis_mass(NMOD)
TYPE(MODAL_EMISSIONS), INTENT(INOUT) :: emis_number(NMOD)


!*       0.5   LOCAL VARIABLES
!              ---------------

INTEGER(KIND=JPIM) :: JL
REAL(KIND=JPRB)    :: NORM, XSEA, AREA_FRAC, TT, T_SCALE, DENS, RG1, RG2 
REAL(KIND=JPRB)    :: EMIS_FAC(KLON)
REAL(KIND=JPRB)    :: NUMBER(KLON), MASS(KLON)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TM5M7_SRC_SS',0,ZHOOK_HANDLE)

!ASSOCIATE(RSSFLX=>YDEAERSRC%RSSFLX)

IF (NSEASALT==0) THEN
    !>>> TvN
    ! The parameterization of Gong (2003)
    ! gives the particle number flux as a function 
    ! of the radius and the 10-m wind speed u_10:
    ! dF/dr80 = f(u_10) x g(r80).
    ! The dependence on wind speed is given by 
    ! the power law f(u_10) = u_10^3.41,
    ! and does not affect the size distribution.
    ! r80 is the radius at 80% humidity, 
    ! which is almost exactly 2.0 times the dry radius
    ! (Lewis and Schwartz, Sea Salt Aerosol Production, 2004).
    ! Note also that in Eq. (2) of Gong 
    ! B is defined in terms of log(10,r) not ln(r).
    !
    ! The number mean radii defined in chem_param.F90,
    ! i.e. 0.09 um for the accumulation mode,
    ! and 0.794 for the coarse mode,
    ! were obtained by fitting an accumulation
    ! and coarse mode to the dF/dln(r),
    ! with r the dry particle radius
    ! (see presentation E. Vignati, 21 December 2005).
    ! It was verified that the size distribution of Gong 
    ! can be reasonable well described using these mode radii.
    !
    ! It is not totally clear how the corresponding prefactors 
    ! for the two modes were obtained.
    ! One way to calculate these prefactors,
    ! is to define two size ranges
    ! and require that the numbers of particles emitted
    ! in both ranges correspond to Gong
    ! This procedure is similar to that used in emission_dust.F90,
    ! but for particle number instead of mass.
    ! Using ranges r1 and r2 for the dry particle radii
    ! with r1 from 0.05 to 0.5 um and r2 from 0.5 to 5 um,
    ! the resulting prefactors are 9.62013e-3 and 9.05658e-4
    ! for the accumulation and coarse mode, respectively.
    ! These numbers are very close to the values
    ! of 9e-3 and 9e-4, obtained by E. Vignati.
    ! They are also insensitive to the
    ! value used for the upper bound of r2.
    ! (Using a value of 10 instead of 5 um,
    ! the prefactors would be 9.62536e-3 and 8.91184e-3.)
    !

    norm=1E4 * 1._JPRB * 1._JPRB ! Conversion from #/cm2 to #/m2

    !TB added zero init  in case no sea
    emis_fac(:) = 0.0_JPRB

    DO JL=KIDIA,KFDIA

!VH     norm = 1.e4 * dxy11(JL) * sec_month ! Used in TM5.. 
!VH   but dxy11 is not available.. use unity area and per sec. instead..

      ! sea fraction
      xsea=1.-PLSM(JL)

      ! sea salt is emitted only over sea without ice cover
      area_frac = xsea * (1.-PCI(JL))
      !write(6566,*)area_frac,xsea,PCI(JL),PCLAKE(JL)
      if (area_frac .LT. 1.e-10_JPRB .or. PCLAKE(JL)>0.0_JPRB) CYCLE

      emis_fac(JL) = norm * area_frac

      ! Wind speed dependence following W10 parameterization
      ! of Salisbury et al. (JGR, 2013; GRL, 2014),
      ! which replaces the dependence of Gong.
      ! Salisbury et al. suggest that "at this stage ...
      ! use of W10 is preferable to W37".
      ! In TM5, W10 gives better results than W37,
      ! which leads to high AOD values compared to MODIS C6,
      ! at mid- and high latitudes.
      ! The same is true for the wind dependence
      ! proposed by Albert et al. (ACPD, 2015).
      ! xwind =  SQRT(u10m_dat(iglbsfc)%data(JL,1)**2+v10m_dat(iglbsfc)%data (JL,1)**2)

      ! Revert to original wind speed dependence of Monahan et al. (1986)
      ! used by Gong (2003):
      !emis_fac(JL) = emis_fac(JL) * W10_fac * xwind**W10_exp
      emis_fac(JL) = emis_fac(JL) * PWIND(JL)**3.41
    ENDDO
    !<<< TvN

    !===================
    ! Accumulation mode
    !===================

    ! emitted numbers
    ! ---------------
    DO JL=KIDIA,KFDIA
       !>>> TvN
       ! sea fraction
       !xsea=1.-PLSM(JL)
 
       !xwind=SQRT(u10m_dat(iglbsfc)%data(JL,1)**2+v10m_dat(iglbsfc)%data (JL,1)**2) 

       !number(JL) = 0.009*xwind**(3.4269)*1e4*dxy11(j)*xsea*(1.-ci_dat(iglbsfc)%data(JL,1))*sec_month   ! #/gridbox/month
       number(JL) = emis_fac(JL)*9.62013e-3 ! #/m2/sec
       !TB emis_fac includes all the rest but first multiplier  
       !although dxy11 is not there, assumed 1 earlier, is that correct?
       ! multiplier needs an explanation!
       ! cm2->m2*land*ice*wind*prefacor (see above)
       !  1e4* (1-PLSM)*(1-pci)* pwind**3.41*9.602e-3

       !Include temperature dependence following Salter et al. (ACP, 2015).
       !We assume that the emissions in our accumulation mode 
       !follow the temperature dependence of the smallest mode described by Santer et al.
       !In reality our mode is in between the smallest and middle mode of Santer et al.,
       !so the temperature dependence might also be.
       !The corresponding third-order polynomial for the smallest mode 
       !decreases with temperature between -1 and about 15 degC,
       !remains almost constant between 15 degC and 25 degC,
       !and shows a very small decrease between 25 degC and 30 degC.
       !The latter dependence seems an artefact of the fitting procedure,
       !and is neglected here.
       !As for the coarse mode,
       !we rescale the polynomial of Salter et al.
       !so that it goes through 1 at some reference temperature,
       !which can currently be set to either 15 or 20 degC.
       !Because of the small difference between 15 and 20 degC,
       !it really doesn't matter which reference value 
       !is chosen for the accumulation mode.
       !This temperature dependency is used to produce 
       !the pre-industrial aerosol climatology v4.0,
       !which is used in all CMIP6 versions of EC-Earth3 that use MACv2-SP,
       !and will therefore also be used in EC-Earth3-AerChem.
       ! PSST in Kelvin
       tt=max(-1.,PSST(JL)-273.15_JPRB)

       !For a reference temperature of 20 degC, use:
       !if (tt .lt. 20.) then
       !  t_scale = -8.88108055e-5*tt*tt*tt + 5.64728654e-3*tt*tt &
       !            -0.118363619*tt + 1.81884421
       !  number(JL) = number(JL) * t_scale
       !endif
       !For a reference temperature of 15 degC, use:
       if (tt .lt. 15.) then
         t_scale = -8.75593266e-5*tt*tt*tt + 5.56770771e-3*tt*tt &
                   -0.116695696*tt + 1.79321393
         number(JL) = number(JL) * t_scale
       endif

       !Using the above relation we still underestimate CDNC by a factor 2 to 4 
       !over the Soutern Ocean.
       !The enhancement factor from Salter et al. (2015) is only 1.9 at -1 degC,
       !while earlier studies measured a stronger enhancement 
       !in the number of submicron particles as SST decreases
       !from about 9 to -1 (or -1.3) degC
       !(Salter et al., JGR, 2014, Bowyer et al., JGR, 1990; 
       !Zabori et al. (ACP, 2012a; 2012b)
       !In these studies the enhancement at about -1 degC varies
       !from ~3 (Salter et al.; Bowyer et al., 1990), 
       !~7 (Zabori et al., 2012b), to ~10 (Zabori et al., 2012a).
       !As the studies don't provide an expression for the enhancement factor
       !as a function of temperature,
       !I have approximated it by a quadratic function 
       !which reaches 1 at t0=9 degC with a zero slope:
       ![(F-1)/(t0-tmin)^2]*(t-t0)^2 + 1 
       !where F is the enhancement factor at tmin=-1 degC.
       !Using t0=9 and tmin=-1, the quadratic coefficient
       !becomes (F-1)/100.
       !I have implemented this expression for F=3, F=5, and F=7.
       !With F=5 or F=7, the zonal mean AOD over the Southern Ocean
       !is overestimated at high latitudes
       !compared to observational estimates from AATSR.
       !For v5.0 of the pre-industrial aerosol climatology
       !we have therefore used F=3.
       !tt=max(-1.,temperature(JL))
       !if (tt .lt. 9.) then
         !F=3:
         !t_scale = 0.02*(tt-9.)*(tt-9.) + 1.0
         !F=5:
         !t_scale = 0.04*(tt-9.)*(tt-9.) + 1.0
         !F=7:
         !t_scale = 0.06*(tt-9.)*(tt-9.) + 1.0
         !number(JL) = number(JL) * t_scale
       !endif
       !<<< TvN
    END DO

    ! vertically distribute according to sector
    ! CALL emission_vdist_by_sector( splittype, 'SS', region, emis_temp(region), emis3d, status )
    
    ! For now fill in surface layer
    emis_number(mode_acs)%d3(KIDIA:KFDIA,KLEV,4)   = number(KIDIA:KFDIA)  !#/m2/sec
    !write(6565,*) number(KIDIA:KFDIA)
    
    ! emitted mass
    ! ------------
    dens = ss_density*0.001  ! in g/cm3
    rg1  = radius_ssa*1e6    ! in micron
    mass(KIDIA:KFDIA) = number(KIDIA:KFDIA)*RPI*4./3. *(rg1*1e-4)**3 &
                    & * EXP(9./2.*log(sigma_lognormal(3))**2)*dens*1e-3  ! kg/m2/sec


    ! vertically distribute according to sector
    ! CALL emission_vdist_by_sector( splittype, 'SS', region, emis_temp(region), emis3d, status )
    emis_mass(mode_acs)%d3(KIDIA:KFDIA,KLEV,4)   = mass(KIDIA:KFDIA)   !kg/m2/sec


    !===================
    ! Coarse mode
    !===================
        
    ! emitted numbers
    ! ---------------
    DO JL=KIDIA,KFDIA
       ! >>> TvN
       ! sea fraction
       !xsea=1.-lsmask_dat(iglbsfc)%data(JL,1)/100.

       !xwind=SQRT(u10m_dat(iglbsfc)%data(JL,1)**2+v10m_dat(iglbsfc)%data (JL,1)**2) 

       !number(JL) = 0.0009*xwind**(3.4195)*dxy11(j)*1e4*xsea*(1.-ci_dat(iglbsfc)%data(JL,1))*sec_month  !#/cm2/s-->#/gridbox/month 
       number(JL) = emis_fac(JL)*9.05658e-4  !#/cm2/s-->#/m2/sec
       !TB
       ! multiplier needs an explanation!

       !Include temperature dependence following Salter et al. (ACP, 2015).
       !We assume that the emissions in our coarse mode 
       !follow the temperature dependence of the largest mode described by Santer et al.
       !Indeed, the mode radii of these modes are close (0.794 vs. 0.75 micron).
       !The corresponding second-order polynomial (valid between -1 and 30 degC)
       !shows almost linear increase with temperature.
       !The higher order terms are therefore irrelevant.
       !We rescale the resulting linear function
       !so that it goes through 1 at some reference temperature,
       !which can currently be set to either 15 or 20 degC. 
       !The resulting expression is then used as
       !an additional factor multiplying the expression from Gong.
       !Gong uses the relation between whitecap coverage and wind speed
       !from Monahan and Muircheartaigh (1980).
       !It was derived from data in the Atlantic and Pacific Ocean,
       !where most measurements were made at sea water temperatures 
       !between 20 and 30 degC: 46 out of 54 points in the Atlantic data set
       !and all 36 points in the Pacific data set. 
       !Based on this, one would argue that a reference temperature
       !in the range 20-25 degC would be most reasonable.
       !As the emissions in the coarse mode increase with temperature,
       !the smaller the reference temperature the higher the emissions.
       !For example, the emissions will increase by 18.76%
       !when the reference temperature is reduced from 20 to 15 degC.
       !The resulting temperature dependence is somewhat weaker than 
       !the quasi-linear dependence obtained by Ovadnevaite et al. 
       ![ACP, 2014; see also Gantt et al., GMD, 2015, Eq. (2)].
       !The function proposed by Jaegle et al. (JGR, 2011)
       !oscillates around or mostly below 
       !(reference temperature at 20 degC or 15 degC, resp.) 
       !our linear function up to about 20 degC, 
       !but increases more strongly at higher temperatures.
       !For CMIP6 a reference temperature of 15 degC is used.
       ! PSST in Kelvin
       tt=min(30.,max(-1.0,PSST(JL)-273.15_JPRB))
       !For a reference temperature of 20 degC, use:
       !t_scale = 0.03159982*tt + 0.36800362
       !For a reference temperature to 15 degC, use:
       t_scale = 0.03752944*tt + 0.43705846
       number(JL) = number(JL) * t_scale
       ! <<< TvN
    ENDDO

    !VH: For now emitted at surface, to be fixed
    emis_number(mode_cos)%d3(KIDIA:KFDIA,KLEV,4)   = number(KIDIA:KFDIA)   !#/grid/month


    ! emitted mass
    ! ------------
    dens = ss_density*0.001  ! im g/cm3
    rg2  = radius_ssc*1e6    ! in micron
    mass(KIDIA:KFDIA) = number(KIDIA:KFDIA)*RPI*4./3. *(rg2*1e-4)**3 * EXP(9./2.*log(sigma_lognormal(4))**2)*dens*1e-3  ! kg/m2/sec

    !For now introduce emissions in surface layer. Should be fixed.
    emis_mass(mode_cos)%d3(KIDIA:KFDIA,KLEV,4)   = mass(KIDIA:KFDIA)    !kg/m2/sec

  ELSEIF (NSEASALT==8) THEN  !HAM gong_SST
     
    CALL SEASALT_EMISSIONS_GONG_SST(KFDIA, KLON, 1,&
         & PSST,PWIND, ss_density, PLSM, PCLAKE, PCI, &
         & emis_mass(mode_acs)%d3(KIDIA:KFDIA,KLEV,4)  , emis_mass(mode_cos)%d3(KIDIA:KFDIA,KLEV,4),&
         & emis_number(mode_acs)%d3(KIDIA:KFDIA,KLEV,4), emis_number(mode_cos)%d3(KIDIA:KFDIA,KLEV,4))
  ELSE

    CALL ABOR1('ABORT: IN TM5_SRC_SS, NSEASALT is NOT 0 or 8!')

  END IF

  ! RCHG -> In AER scheme there is a flag named LVDFTRAC that might be related with not vertical diffusion. In that case,
  !        the tendencies seems to be re-scaled in vertical layers "manually". The flux themselves are not changed.

IF (LHOOK) CALL DR_HOOK('TM5M7_SRC_SS',1,ZHOOK_HANDLE)
END SUBROUTINE TM5M7_SRC_SS

