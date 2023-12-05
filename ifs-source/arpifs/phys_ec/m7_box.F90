!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! m7_box.f90
!!
!! \brief
!! Boxmodel of the M7 aerosol microphysics scheme currently in ECHAM6-HAMMOZ. 
!! To find out what is exactly the corresponding version number of ECHAM6-HAMMOZ,
!! please read:
!! echam6-hammoz_base_version.txt
!!
!! \author Grazia Frontoso, Nov 2012 (C2SM-ETHZ)
!!
!! \responsible_coder
!! Grazia Frontoso, grazia.frontoso@env.ethz.ch
!!
!! \revision_history
!!
!! \limitations
!! None
!!
!! \bibliographic_references
!!   - Vignati E., Wilson J. and Stier P., M7: a size resolved aerosol mixture module 
!!     for the use in global aerosol models, JGR 109, D22 202, doi:10.1029/2003JD004 485, 2004.
!!   - Stier P. et al, The aerosol-climate model ECHAM5-HAM, ACP, 5, 1125-1156, 2005
!!
!!
!! \copyright
!! Copyright and licencing conditions are defined in the ECHAM-HAMMOZ
!! licencing agreement to be found at:
!! https://redmine.hammoz.ethz.ch/projects/hammoz/wiki/1_Licencing_conditions
!! The ECHAM-HAMMOZ software is provided "as is" and without warranty of any kind.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM m7_box
  
  USE mo_kind,           ONLY: dp
  USE mo_ham_m7,         ONLY: m7
  USE mo_ham_init,       ONLY: ham_init_memory, ham_free_memory, &
                               ham_initialize, start_ham
  USE mo_ham_m7ctl,      ONLY: kproma, kbdim, klev,              &
                               nmod, nsol, naerocompmax,         &
                               iso4ns, iso4ks, iso4as, iso4cs,   &
                               ibcks,  ibcas,  ibccs,  ibcki,    &
                               iocks,  iocas,  ioccs,  iocki,    &
                               issas,  isscs,                    &
                               iduas,  iducs,  iduai,  iduci,    &
                               inucs,  iaits,  iaccs,  icoas,    &
                               iaiti,  iacci,  icoai
  USE mo_namelist,       ONLY: open_nml, position_nml, POSITIONED
  
  IMPLICIT NONE 
  !--- Initial conditions for m7.box :
  REAL(dp):: pap(kbdim,klev)   = 1013.25_dp ! Pressure [Pa]
  REAL(dp):: pt(kbdim,klev)    =   265.0_dp ! Temperature at [K]
  REAL(dp):: prh(kbdim,klev)   =     0.3_dp ! Relative humidity at [0,1]
  REAL(dp):: pipr(kbdim,klev)  =     5.0_dp ! Ionization rate [cm-3 s-1]
  REAL(dp):: pdgso4(kbdim,klev)=   5.0E4_dp ! d[H2SO4(g)]/dt [cm-3 s-1]
  REAL(dp):: pgso4(kbdim,klev) =   1.0E4_dp ! [H2SO4(g)] [cm-3]
  REAL(dp):: paclc(kbdim,klev) =     0.0_dp ! Cloud cover [0,1]
  REAL(dp):: pforest(kbdim)    =     1.0_dp ! Forest fraction
  REAL(dp):: paerml(kbdim,klev,naerocompmax) 
             ! aerosol mass concentrations for individual compounds
             ! [molec. cm-3 for sulfate and ug m-3 for bc, oc, ss, and dust]
  REAL(dp):: paernl(kbdim,klev,nmod)
             ! aerosol number concentrations for each mode [cm-3]
  REAL(dp):: pm6rp(kbdim,klev,nmod) 
             ! mean mode actual radius for each mode (wet for soluble and dry for insoluble modes) [cm]
  REAL(dp):: pm6dry(kbdim,klev,nsol) 
             ! dry radius for soluble modes [cm]
  REAL(dp):: prhop(kbdim,klev,nmod) 
             ! mean mode particle density for each mode [g cm-3]
  REAL(dp):: pww(kbdim,klev,nmod)           
             ! aerosol water content for each mode [kg(water) m-3(air)]
  INTEGER :: ierr, inml, iunit
  INCLUDE 'ini_m7ctl.inc'  
  ! Initialization:
  paerml(1:kproma,:,:) = 0.0_dp
  paernl(1:kproma,:,:) = 0.0_dp
  pm6rp(1:kproma,:,:)  = 0.0_dp
  pm6dry(1:kproma,:,:) = 0.0_dp
  prhop(1:kproma,:,:)  = 0.0_dp
  pww(1:kproma,:,:)    = 0.0_dp
  CALL start_ham
  CALL ham_initialize  
  CALL ham_init_memory
 ! Read initial values from namelist: ini_m7ctl.inc
      inml = open_nml('namelist.echam') 
      iunit = position_nml ('INI_M7CTL', inml, status=ierr)
      SELECT CASE (ierr)
      CASE (POSITIONED)
      READ (iunit, ini_m7ctl)
      END SELECT
  CALL m7(pap,     pt,      prh,           &  ! Pressure, temperature, RH 
          pgso4,   pdgso4,                 &  ! [H2SO4(g)], derivative (production rate)
          paerml,  paernl,                 &  ! Aerosol mass and number
          pm6rp,   pm6dry, prhop,   pww,   &  ! Aerosol properties
          pipr,    paclc,  pforest         )  ! Ionization rate, cloud cover, forest fraction
  !
  ! Print the output:
  !
  
  WRITE(*,'(A)') '-------------------'
  WRITE(*,'(A)') 'M7 boxmodel outputs'
  WRITE(*,'(A)') '-------------------'
      
      ! H2SO4 gas phase concentration:
      
      WRITE(*,'(A,E14.7)') '[H2SO4] (cm-3)                      : ',pgso4(1:kproma,:)
      
      ! Aerosol mass concentrations for individual compounds, (molec. cm-3 for sulfate and
      ! ug m-3 for bc, oc, ss, and dust):
      
      ! Sulfate:
      
      WRITE(*,'(A)') ''
      WRITE(*,'(A,E14.7)') 'NS SO4-- concents (molec. cm-3)     : ',paerml(1:kproma,:,iso4ns)
      WRITE(*,'(A,E14.7)') 'KS SO4-- concents (molec. cm-3)     : ',paerml(1:kproma,:,iso4ks)
      WRITE(*,'(A,E14.7)') 'AS SO4-- concents (molec. cm-3)     : ',paerml(1:kproma,:,iso4as)
      WRITE(*,'(A,E14.7)') 'CS SO4-- concents (molec. cm-3)     : ',paerml(1:kproma,:,iso4cs)
      
      ! Black Carbon:
      
      WRITE(*,'(A)') ''
      WRITE(*,'(A,E14.7)') 'BC KS (ug cm-3)                     : ',paerml(1:kproma,:,ibcks)
      WRITE(*,'(A,E14.7)') 'BC KI (ug cm-3)                     : ',paerml(1:kproma,:,ibcki)
      WRITE(*,'(A,E14.7)') 'BC AS (ug cm-3)                     : ',paerml(1:kproma,:,ibcas)
      WRITE(*,'(A,E14.7)') 'BC CS (ug cm-3)                     : ',paerml(1:kproma,:,ibccs)
      
      ! Organic Carbon:
      
      WRITE(*,'(A)') ''
      WRITE(*,'(A,E14.7)') 'OC KS (ug cm-3)                     : ',paerml(1:kproma,:,iocks)
      WRITE(*,'(A,E14.7)') 'OC KI (ug cm-3)                     : ',paerml(1:kproma,:,iocki)
      WRITE(*,'(A,E14.7)') 'OC AS (ug cm-3)                     : ',paerml(1:kproma,:,iocas)
      WRITE(*,'(A,E14.7)') 'OC CS (ug cm-3)                     : ',paerml(1:kproma,:,ioccs)
      
      ! Sea Salt:
      
      WRITE(*,'(A)') ''
      WRITE(*,'(A,E14.7)') 'SS AS (ug cm-3)                     : ',paerml(1:kproma,:,issas)
      WRITE(*,'(A,E14.7)') 'SS CS (ug cm-3)                     : ',paerml(1:kproma,:,isscs)
      
      ! Dust:
      
      WRITE(*,'(A)') ''
      WRITE(*,'(A,E14.7)') 'DU AS (ug cm-3)                     : ',paerml(1:kproma,:,iduas)
      WRITE(*,'(A,E14.7)') 'DU AI (ug cm-3)                     : ',paerml(1:kproma,:,iduai)
      WRITE(*,'(A,E14.7)') 'DU CS (ug cm-3)                     : ',paerml(1:kproma,:,iducs)
      WRITE(*,'(A,E14.7)') 'DU CI (ug cm-3)                     : ',paerml(1:kproma,:,iduci)
      
      ! Particle number concentrations for each mode (cm-3):
      
      WRITE(*,'(A)') ''
      WRITE(*,'(A,E14.7)') '[NS] (cm-3)                         : ',paernl(1:kproma,:,inucs)
      WRITE(*,'(A,E14.7)') '[KI] (cm-3)                         : ',paernl(1:kproma,:,iaiti)
      WRITE(*,'(A,E14.7)') '[KS] (cm-3)                         : ',paernl(1:kproma,:,iaits)
      WRITE(*,'(A,E14.7)') '[AI] (cm-3)                         : ',paernl(1:kproma,:,iacci)
      WRITE(*,'(A,E14.7)') '[AS] (cm-3)                         : ',paernl(1:kproma,:,iaccs)
      WRITE(*,'(A,E14.7)') '[CI] (cm-3)                         : ',paernl(1:kproma,:,icoai)
      WRITE(*,'(A,E14.7)') '[CS] (cm-3)                         : ',paernl(1:kproma,:,icoas)
      
      ! Particle radii for each mode (nm):
 
      WRITE(*,'(A)') ''
      WRITE(*,'(A,E14.7)') 'NS geom. mean rad. (nm)             : ',1.0E9_dp*pm6rp(1:kproma,:,inucs)
      WRITE(*,'(A,E14.7)') 'KI geom. mean rad. (nm)             : ',1.0E9_dp*pm6rp(1:kproma,:,iaiti)
      WRITE(*,'(A,E14.7)') 'KS geom. mean rad. (nm)             : ',1.0E9_dp*pm6rp(1:kproma,:,iaits)
      WRITE(*,'(A,E14.7)') 'AI geom. mean rad. (nm)             : ',1.0E9_dp*pm6rp(1:kproma,:,iacci)
      WRITE(*,'(A,E14.7)') 'AS geom. mean rad. (nm)             : ',1.0E9_dp*pm6rp(1:kproma,:,iaccs)
      WRITE(*,'(A,E14.7)') 'CI geom. mean rad. (nm)             : ',1.0E9_dp*pm6rp(1:kproma,:,icoai)
      WRITE(*,'(A,E14.7)') 'CS geom. mean rad. (nm)             : ',1.0E9_dp*pm6rp(1:kproma,:,icoas)
      
      WRITE(*,'(A)') ''
      WRITE(*,'(A,E14.7)') 'NS dry rad. (nm)                    : ',1.0E9_dp*pm6dry(1:kproma,:,inucs)
      WRITE(*,'(A,E14.7)') 'KS dry rad. (nm)                    : ',1.0E9_dp*pm6dry(1:kproma,:,iaits)
      WRITE(*,'(A,E14.7)') 'AS dry rad. (nm)                    : ',1.0E9_dp*pm6dry(1:kproma,:,iaccs)
      WRITE(*,'(A,E14.7)') 'CS dry rad. (nm)                    : ',1.0E9_dp*pm6dry(1:kproma,:,icoas)
      ! Mean mode particle density for each mode [g cm-3]
      WRITE(*,'(A)') ''
      WRITE(*,'(A,E14.7)') 'NS particle density (g cm-3)        : ',prhop(1:kproma,:,inucs)
      WRITE(*,'(A,E14.7)') 'KI particle density (g cm-3)        : ',prhop(1:kproma,:,iaiti)
      WRITE(*,'(A,E14.7)') 'KS particle density (g cm-3)        : ',prhop(1:kproma,:,iaits)
      WRITE(*,'(A,E14.7)') 'AI particle density (g cm-3)        : ',prhop(1:kproma,:,iacci)
      WRITE(*,'(A,E14.7)') 'AS particle density (g cm-3)        : ',prhop(1:kproma,:,iaccs)
      WRITE(*,'(A,E14.7)') 'CI particle density (g cm-3)        : ',prhop(1:kproma,:,icoai)
      WRITE(*,'(A,E14.7)') 'CS particle density (g cm-3)        : ',prhop(1:kproma,:,icoas)
      
      ! Water content for each mode [kg(water) m-3(air)]
  
      WRITE(*,'(A)') ''
      WRITE(*,'(A,E14.7)') 'NS water content (kg(H2O) m-3(air)) : ',pww(1:kproma,:,inucs)
      WRITE(*,'(A,E14.7)') 'KI water content (kg(H2O) m-3(air)) : ',pww(1:kproma,:,iaiti)
      WRITE(*,'(A,E14.7)') 'KS water content (kg(H2O) m-3(air)) : ',pww(1:kproma,:,iaits)
      WRITE(*,'(A,E14.7)') 'AI water content (kg(H2O) m-3(air)) : ',pww(1:kproma,:,iacci)
      WRITE(*,'(A,E14.7)') 'AS water content (kg(H2O) m-3(air)) : ',pww(1:kproma,:,iaccs)
      WRITE(*,'(A,E14.7)') 'CI water content (kg(H2O) m-3(air)) : ',pww(1:kproma,:,icoai)
      WRITE(*,'(A,E14.7)') 'CS water content (kg(H2O) m-3(air)) : ',pww(1:kproma,:,icoas)
  
      WRITE(*,'(A)') ''
  
  ! Clean up:
  
    CALL ham_free_memory
END PROGRAM m7_box
