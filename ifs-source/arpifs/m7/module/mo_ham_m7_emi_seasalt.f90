!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! mo_ham_m7_emi_seasalt.f90
!!
!! \brief
!! This module contains several sea salt emission schemes.
!!
!! \author Kai Zhang (MPI-Met) (see also the individual subroutines)
!!
!! \responsible_coder
!! Kai Zhang, kai.zhang@pnnl.gov
!!
!! \revision_history
!!   -# Kai Zhang (MPI-Met) - original code (2009)
!!
!! \limitations
!! [ Start an optional warning here ]
!!
!! \details
!! Use nseasalt in namelist hamctl to select the scheme you want:
!!   - seasalt_emissions_monahan:   nseasalt=1
!!   - seasalt_emissions_lsce:      nseasalt=2
!!   - seasalt_emissions_mh:        nseasalt=4
!!   - seasalt_emissions_guelle:    nseasalt=5
!!   - seasalt_emissions_gong:      nseasalt=6
!!   - seasalt_emissions_long       nseasalt=7
!!   - seasalt_emissions_gong_SST   nseasalt=8
!!
!! \bibliographic_references
!!   - see individual seasalt scheme routines
!!
!! \belongs_to
!!  HAMMOZ
!!
!! \copyright
!! Copyright and licencing conditions are defined in the ECHAM-HAMMOZ
!! licencing agreement to be found at:
!! https://redmine.hammoz.ethz.ch/projects/hammoz/wiki/1_Licencing_conditions
!! The ECHAM-HAMMOZ software is provided "as is" and without warranty of any kind.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Anton Laakso: Currently 
!----

MODULE mo_ham_m7_emi_seasalt

    !---inherited types, data and functions
    USE mo_kind,           ONLY: dp
  
    IMPLICIT NONE
  
    !---public member functions
    PUBLIC :: start_emi_seasalt
    !PUBLIC :: seasalt_emissions_monahan     ! nseasalt=1
    !PUBLIC :: seasalt_emissions_lsce        ! nseasalt=2
    !!PUBLIC :: seasalt_emissions_martensson  !nseasalt=3
    !PUBLIC :: seasalt_emissions_mh          !nseasalt=4
    !PUBLIC :: seasalt_emissions_guelle      !nseasalt=5
    !PUBLIC :: seasalt_emissions_gong        !nseasalt=6 
    !PUBLIC :: seasalt_emissions_long        !nseasalt=7
    PUBLIC :: seasalt_emissions_gong_SST    !nseasalt=8
  
    !---module data
    REAL(dp), PARAMETER, PRIVATE :: ppww = 3.41_dp      ! exponent of wind speed |u| (|u|**ppww) 
  
    INTEGER,  PARAMETER, PRIVATE :: nbin = 300          ! number of bins for the bin schemes
                                                        ! (Monahan (nseasalt=4), Guelle or Gong, Long)
    REAL(dp), PARAMETER, PRIVATE :: dmta = 0.100E-06_dp ! lower diameter [m], bin schemes 
    REAL(dp), PARAMETER, PRIVATE :: dmtd = 1.000E-05_dp ! upper diameter [m], bin schemes 
  
    REAL(dp), PARAMETER, PRIVATE :: dmtb_gong = 0.221E-06_dp ! diameter limit, [m], corresponding to 
                                                             ! 0.2um wet radius at 80% RH (G03), 0.2*2/1.814 
    REAL(dp), PARAMETER, PRIVATE :: dmtb_guelle = 8.000E-06_dp ! dry diameter limit, [m], corresonding 
                                                               ! to 4um radius in G01 
  
    REAL(dp), PARAMETER, PRIVATE :: dbeg(3) = (/0.050E-6_dp, 0.100E-6_dp, 1.000E-6_dp/)  ! ait, acc, coa
    REAL(dp), PARAMETER, PRIVATE :: dend(3) = (/0.100E-6_dp, 1.000E-6_dp, 1.000E-5_dp/)  ! ait, acc, coa
  
    REAL(dp), PRIVATE :: dmt(nbin)   
  
    REAL(dp), PRIVATE :: rm(nbin)
    REAL(dp), PRIVATE :: rd(nbin)
    REAL(dp), PRIVATE :: bmn(nbin) 
  
    REAL(dp), PRIVATE  :: ss1_mon                   ! sea salt flux factor 1, Monahan (nseasalt=1) scheme    
    REAL(dp), PRIVATE  :: ss2_mon                   ! sea salt flux factor 2, Monahan (nseasalt=1) scheme    
  
  
  CONTAINS
  
    SUBROUTINE start_emi_seasalt
  
      USE mo_kind,         ONLY: dp
      USE mo_math_constants, ONLY: pi
      USE mo_ham,          ONLY: nseasalt
      !USE mo_exception,    ONLY: message, message_text, em_param
      !USE mo_util_string,  ONLY: separator
  
      IMPLICIT NONE
  
      !---local variables
      !   intermediate variables in calculating flux factors ss1 and ss2 for the monahan (nseasalt=1) scheme
      REAL(dp) :: zr1, zr2, zb1, zb2, zx1, zx2, zdr1, zdr2
      REAL(dp) :: zfact
  
      !   variables for constructing bins in the bin schemes
      REAL(dp) :: zdx, zdd
      INTEGER :: m
  
      !---executable procedure
  
      !---initialize the monahan (nseasalt=1) scheme (copy/paste ECHAM5/HAM2 code)
      zr1=0.416_dp
      zr2=3.49_dp
      zb1=0.58_dp-1.54_dp*LOG10(zr1)
      zb2=0.58_dp-1.54_dp*LOG10(zr2)
      zx1=10._dp**(1.19_dp*EXP(-zb1*zb1))
      zx2=10._dp**(1.19_dp*EXP(-zb2*zb2))
      zdr1=0.5_dp
      zdr2=4.5_dp
      zfact=1.373_dp*4._dp/3._dp*pi*1.15e3_dp
      ss1_mon=zfact*(1._dp+0.057_dp*zr1**1.05_dp)*zx1*zdr1*1.e-18_dp
      ss2_mon=zfact*(1._dp+0.057_dp*zr2**1.05_dp)*zx2*zdr2*1.e-18_dp
  
      !IF (nseasalt==1) THEN
      !   CALL message('', separator)
      !   CALL message('', 'Monahan seasalt emissions:', level=em_param)
      !   WRITE(message_text, '(a,e25.15,a,e25.15)') 'Factors for sea salt fluxes: 1st = ',ss1_mon, '2nd = ', ss2_mon
      !   CALL message('', message_text, level=em_param)
      !   CALL message('', separator)
      !END IF
  
      !---construction of bins
      ! dDp, take LOG scale 
  
      zdx = (LOG(dmtd) - LOG(dmta) ) / REAL(nbin,dp)
  
      zdd = 0._dp
  
      DO m = 1, nbin
         dmt(m) = EXP(LOG(dmta) + zdd)
         zdd = zdd + zdx
      END DO
  
      ! dry radius (m)  
    
      rd(:) = dmt(:) * 0.5_dp 
  
      ! wet radius (um) at RH=80%
  
      rm(:) = 2.0_dp*rd(:)*1.E+06_dp
#ifdef HAMMOZ    
      ! B: monahan and guelle schemes, also for larger particle in Gong scheme
      bmn(:) = ( 0.380_dp - LOG10( rm(:) ) ) / 0.650_dp
  
      ! B: overwrite for smaller particles in the Gong scheme
      IF (nseasalt == 6 .or. nseasalt == 8 ) THEN
         DO m = 2,nbin
            IF (dmt(m).GT.dmta .and. dmt(m).le.dmtb_gong) THEN 
               bmn(m) = ( 0.433_dp - LOG10( rm(m) ) ) / 0.433_dp
            END IF
         END DO
      END IF
#else  
   
  
   ! alaak: openIFS - Gong used for all bins
   IF (nseasalt == 6 .or. nseasalt == 8 ) THEN
      bmn(:) = ( 0.433_dp - LOG10( rm(:) ) ) / 0.433_dp
   END IF


#endif   
    END SUBROUTINE start_emi_seasalt
      

    SUBROUTINE seasalt_emissions_gong_SST(kproma, kbdim, krow , sst, wind10m, ss_density, slf, alake, seaice, pmassf_as, pmassf_cs,pnumf_as, pnumf_cs, SSCAL)

        !  
        ! Description:
        ! ------------
        ! Calculates the emitted sea salt flux from the 10m wind speed following
        ! Gong, 2003 
        !
        !   method: M86/lab
        !   size: 0.07um < r80 < 20um
        !   wind speed: N.A.
        !   SST: Explicit SST dependence, according to Sofiev et al 2011
        ! 
        ! currently Aitken mode particles are negelected. 
        !   
        ! Authors:
        ! ------------ 
        ! Kai Zhang, MPI-Met, 2009, modified by I. Tegen, 2016
        !   
        ! References:
        ! ------------
        ! 1. Monahan, E. C., D. E. Spiel, and K. L. Davidson, 
        !    A model of marine aerosol generation via whitecaps and wave disruption, 
        !    in Oceanic Whitecaps, edited by E. C. Monahan and G. MacNiochaill, 
        !    pp. 167–193, D. Reidel, Norwell, Mass., 1986. (M86)  
        ! 
        ! 2. S.L. Gong, 
        !    A parameterization of sea-salt aerosol source function for sub- and super-micron particles, 
        !    Global Biogeochemical Cycles 17 (2003) (4), p. 1097.   (G03) 
        !
        ! 3. Sofiev, M., Soares,J.,Prank, M., deLeeuw, G., Kukkonen, J., A
        ! regional-to-global model of emission and transport of sea salt particles
        ! in the atmosphere, JGR (doi:0148‐0227/11/2010JD014713)
      
        !
      
        USE mo_kind,         ONLY: dp
        USE mo_math_constants, ONLY: pi
        USE mo_species,      ONLY: speclist
        !USE mo_ham_species,  ONLY: id_ss
        !USE mo_memory_g3b,   ONLY: slf, alake, seaice
        !USE mo_vphysc,       ONLY: vphysc
      
        IMPLICIT NONE
      
        !--- Parameters:
        !    -
      
      
        !--- I/O:
      
        INTEGER, INTENT(in)    :: kproma               !kproma
        INTEGER, INTENT(in)    :: kbdim                !column  
        INTEGER, INTENT(in)    :: krow                 !chunk 
        REAL(dp), INTENT(in)    :: sst(kbdim),wind10m(kbdim),ss_density,slf(kbdim),alake(kbdim),seaice(kbdim)
        REAL(dp),INTENT(out)   :: pmassf_as(kbdim)    ! mass flux of ss acc particles
        REAL(dp),INTENT(out)   :: pmassf_cs(kbdim)    ! mass flux of ss coa particles
        REAL(dp),INTENT(out)   :: pnumf_as(kbdim)    ! number flux of ss acc particles
        REAL(dp),INTENT(out)   :: pnumf_cs(kbdim)    ! number flux of ss coa particles
        REAL(dp), INTENT(in)   ::SSCAL          !SEASALT deactivation !Mch
        !--- Local:
      
        REAL(dp):: zseafrac(kbdim)         ! fraction of the gridcell covered by
                                           ! non-iced sea water [0.-1.]
        REAL(dp):: zmassf_ks(kbdim)        ! mass   flux of ss ait particles  (currently not supported)
        REAL(dp):: znumf_ks(kbdim)         ! number flux of ss ait particles  (currently not supported)
      
      
        
      
      
      
        REAL(dp):: fi(kbdim,nbin) 
        REAL(dp):: p0,p1,p2,p3,dr
        REAL(dp):: zav !particle volumn  
        REAL(dp):: SST_corr(1:kproma),SST_corr_1,SST_corr_2,dmtum(nbin)
        REAL(dp):: SST_mask1(1:kproma),SST_mask2(1:kproma),SST_mask3(1:kproma),SST_corr_all(1:kproma)
      
      
      
      
        INTEGER :: m
      
        ! initialize number flux for each bin (#/m2/s) 
      
        fi = 0._dp
      
       ! calculate fraction of the gridcell of non ice-covered water
#ifdef HAMMOZ      
        zseafrac(1:kproma) = (1._dp-slf(1:kproma)-alake(1:kproma))*(1._dp-seaice(1:kproma))
#else
        !alaak: seaice fraction of gridbox in oifs:
        zseafrac(1:kproma) = 1._dp-slf(1:kproma)-alake(1:kproma)-seaice(1:kproma)
#endif
        zseafrac(1:kproma) = MAX(0._dp,MIN(zseafrac(1:kproma),1._dp))
      
        !>>SF #458 (replacing where statements)
        zseafrac(1:kproma) = MERGE( &
                                   0._dp, &
                                   zseafrac(1:kproma), &
                                   (slf(1:kproma) > 0.5_dp))
        !<<SF #458 (replacing where statements)
      
        ! initialization 
      
        zmassf_ks(1:kproma) = 0._dp 
        pmassf_as(1:kproma) = 0._dp 
        pmassf_cs(1:kproma) = 0._dp 
      
        znumf_ks(1:kproma) = 0._dp
        pnumf_as(1:kproma) = 0._dp
        pnumf_cs(1:kproma) = 0._dp
      
        ! loop over bins 
      
        DO m = 2,nbin
      
           !dLOGdrm 
      
           dr = rm(m) - rm(m-1) 
#ifdef HAMMOZ       
           IF (dmt(m).GT.dmta .and. dmt(m).le.dmtb_gong) THEN 
      
             p0 = 4.7_dp*(1._dp+30._dp*rm(m))**(-0.017_dp*rm(m)**(-1.44_dp)) 
             p1 = 1.373_dp*rm(m)**(-p0) 
             p2 = 1._dp + 0.057_dp*rm(m)**3.45_dp
             p3 = 10**(1.607_dp*EXP(-bmn(m)**2)) 
      
             fi(1:kproma,m) = p1*p2*p3*dr*wind10m(1:kproma)**ppww
      
           ELSEIF (dmt(m).GT.dmtb_gong .and. dmt(m).le.dmtd) THEN 
      
             p1 = 1.373_dp*rm(m)**(-3) 
             p2 = 1._dp + 0.057_dp*rm(m)**1.05_dp
             p3 = 10**(1.19_dp*EXP(-bmn(m)**2))  
      
             fi(1:kproma,m) = p1*p2*p3*dr*wind10m(1:kproma)**ppww
      
           END IF
#else
   ! alaak: openIFS - Gong used for all bins
   IF (dmt(m).GT.dmta .and. dmt(m).le.dmtd) THEN    
      p0 = 4.7_dp*(1._dp+30._dp*rm(m))**(-0.017_dp*rm(m)**(-1.44_dp)) 
      p1 = 1.373_dp*rm(m)**(-p0) 
      p2 = 1._dp + 0.057_dp*rm(m)**3.45_dp
      p3 = 10**(1.607_dp*EXP(-bmn(m)**2)) 

      fi(1:kproma,m) = p1*p2*p3*dr*wind10m(1:kproma)**ppww
   END IF
#endif
           ! SST correction according to Sofiev et al 2011
       
           dmtum(m)=dmt(m)*1.e06  ! dmt -> micrometers
       
            SST_corr_1 = 0.092e0*dmtum(m)**(-0.96e0) ! valid for Tw=271.15K ; -2°C      ! T base (Long)  25 deg
            SST_corr_2 = 0.15e0*dmtum(m)**(-0.88e0) ! valid for Tw=278.15K ; 5°C
       
           ! SST_corr_1 = 0.19e0*dmtum(m)**(-0.60e0) ! valid for Tw=271.15K ; -2°CC      T base  (Long) 15 deg
           ! SST_corr_2 = 0.31e0*dmtum(m)**(-0.56e0) ! valid for Tw=278.15K ; 5°C
       
           !SST_corr_1 = 0.13e0*dmtum(m)**(-0.78e0) ! valid for Tw=271.15K ; -2°CC        T base (Long) 20 deg
           !SST_corr_2 = 0.22e0*dmtum(m)**(-0.70e0) ! valid for Tw=278.15K ; 5°C
       
       
       
           SST_corr(1:kproma) = (SST_corr_1*(278.15e0-sst(1:kproma)) &
                                     +SST_corr_2*(sst(1:kproma)-271.15e0))/7.e0
       
           SST_mask1(1:kproma) = MERGE (SST_corr(1:kproma),0._dp,(sst(1:kproma) .LE. 278.15 ))
       
            SST_corr_1 = 0.15e0*dmtum(m)**(-0.88e0) ! valid for Tw=278.15K ; 5°C
            SST_corr_2 = 0.48e0*dmtum(m)**(-0.36e0) ! valid for Tw=288.15K ; 15°C
       
           ! SST_corr_1 = 0.31e0*dmtum(m)**(-0.56e0) ! valid for Tw=278.15K ; 5°C
           ! SST_corr_2 = 1.e0 ! valid for Tw=288.15K ; 15°
       
           !SST_corr_1 = 0.22e0*dmtum(m)**(-0.70e0) ! valid for Tw=278.15K ; 5°C
           !SST_corr_2 = 0.70e0*dmtum(m)**(-0.18e0)  ! valid for Tw=288.15K ; 15°C
       
       
           SST_corr(1:kproma) = (SST_corr_1*(288.15e0-sst(1:kproma)) &
                                   +  SST_corr_2*(sst(1:kproma)-278.15e0))/1.e1
       
           SST_mask2(1:kproma) = MERGE (SST_corr(1:kproma), 0._dp,        &
                                             (sst(1:kproma) .GT. 278.15 &
                                        .AND. sst(1:kproma) .LE. 288.15))
       
       
            SST_corr_1 = 0.48e0*dmtum(m)**(-0.36e0) ! valid for Tw=288.15K ; 15°C
            SST_corr_2 = 1.e0 ! valid for Tw=298.15K ; 25°C
       
           ! SST_corr_1 = 1.e0! valid for Tw=288.15K ; 15°C
           ! SST_corr_2 = 2.08e0*dmtum(m)**(0.36e0) ! valid for Tw=298.15K ; 25°C
       
       
           !SST_corr_1 = 0.70e0*dmtum(m)**(-0.18e0)  ! valid for Tw=288.15K ; 15°C
           !SST_corr_2 = 1.45e0*dmtum(m)**(0.18e0) ! valid for Tw=298.15K ; 25°C
       
           SST_corr(1:kproma) = (SST_corr_1*(298.15e0-sst(1:kproma)) &
                                     +SST_corr_2*(sst(1:kproma)-288.15e0))/1.e1
       
           SST_mask3(1:kproma) = MERGE (SST_corr(1:kproma), 0._dp, &
                                             (sst(1:kproma) .GT. 288.15))
       
           SST_corr_all(1:kproma)=SST_mask1(1:kproma)+SST_mask2(1:kproma)+SST_mask3(1:kproma)
           !SST_corr_all(1:kproma) = MERGE (1._dp, SST_corr_all(1:kproma), &
                !(vphysc%tsw(1:kproma) .GT. 298.15))  ! optional, limit T dependence to <25 deg C
       
           fi(1:kproma,m)=fi(1:kproma,m)*SST_corr_all(1:kproma)
      
      
           zav = ss_density*4._dp/3._dp*pi*rd(m)**3
           ! mass flux (kg/m2/s) and number flux (#/m2/s) 
      
           IF (dmt(m).GT.dbeg(1) .AND. dmt(m).LE.dend(1) ) THEN
              znumf_ks(1:kproma) = znumf_ks(1:kproma) + fi(1:kproma,m)*zseafrac(1:kproma)
              zmassf_ks(1:kproma) = zmassf_ks(1:kproma) + fi(1:kproma,m)*zseafrac(1:kproma)*zav 
           END IF
      
           IF (dmt(m).GT.dbeg(2) .AND. dmt(m).LE.dend(2) ) THEN
              pnumf_as(1:kproma) = pnumf_as(1:kproma) + fi(1:kproma,m)*zseafrac(1:kproma)    * SSCAL
              pmassf_as(1:kproma) = pmassf_as(1:kproma) + fi(1:kproma,m)*zseafrac(1:kproma)*zav * SSCAL
           END IF
      
           IF (dmt(m).GT.dbeg(3) .AND. dmt(m).LE.dend(3) ) THEN
              pnumf_cs(1:kproma) = pnumf_cs(1:kproma) + fi(1:kproma,m)*zseafrac(1:kproma) * SSCAL
              pmassf_cs(1:kproma) = pmassf_cs(1:kproma) + fi(1:kproma,m)*zseafrac(1:kproma)*zav * SSCAL 
           END IF
      
        END DO
      
        END SUBROUTINE seasalt_emissions_gong_SST
      
  
  
  END MODULE mo_ham_m7_emi_seasalt
  
