!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename
!! OIFS_to_HAM.f90
!!
!! \brief
!! Contains all the variables for needed to make HAM compatible with OIFS.
!!
!! \author Eemeli Holopainen (FMI)
!!
!! \responsible_coder
!! Eemeli Holopainen, eemeli.holopainen@fmi.fi
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE OIFS_to_HAM

  USE mo_ham, ONLY: nclass, naerocomp, subm_ngasspec

  IMPLICIT NONE

  PUBLIC :: init_ind_oifs_ham
  
  !-->eehol: allocatable integer list for HAM and OIFS
  type  ind_oifs_ham_type
      INTEGER, PUBLIC, ALLOCATABLE :: ind_class_OIFS(:)    !eehol: index list for sizeclass OIFS tracers
      INTEGER, PUBLIC, ALLOCATABLE :: ind_mass_OIFS(:)    !eehol: index list for mass OIFS tracers
      INTEGER, PUBLIC, ALLOCATABLE :: ind_gas_OIFS(:)     !eehol: index list for gas OIFS tracers
      INTEGER, PUBLIC, ALLOCATABLE :: ind_cloud_OIFS(:)     !eehol: index list for cloud OIFS tracers
      INTEGER, PUBLIC, ALLOCATABLE :: ind_class_HAM(:)     !eehol: index list for sizeclass HAM tracers
      INTEGER, PUBLIC, ALLOCATABLE :: ind_mass_HAM(:)     !eehol: index list for mass HAM tracers
      INTEGER, PUBLIC, ALLOCATABLE :: ind_gas_HAM(:)      !eehol: index list for gas HAM tracers
      INTEGER, PUBLIC, ALLOCATABLE :: ind_cloud_HAM(:)     !eehol: index list for cloud HAM tracers
  end type ind_oifs_ham_type
  !<--eehol

  TYPE(ind_oifs_ham_type) :: ind_oifs_ham
  !!$OMP THREADPRIVATE(ind_oifs_ham)

CONTAINS

  SUBROUTINE init_ind_oifs_ham(knclass,knaerocomp,ksubm_ngasspec,kcloudind)
    
    ! *ind_oif_ham* allocates and initializes the index list for OIFS tracers and HAM tracers
    ! Authors:
    ! -------
    ! Eemeli Holopainen, FMI                4/2022

    INTEGER, INTENT(IN) :: knclass,knaerocomp,ksubm_ngasspec,kcloudind
    
    IF (ALLOCATED(ind_oifs_ham%ind_class_OIFS)) DEALLOCATE(ind_oifs_ham%ind_class_OIFS)
    IF (ALLOCATED(ind_oifs_ham%ind_class_HAM))  DEALLOCATE(ind_oifs_ham%ind_class_HAM)
    IF (ALLOCATED(ind_oifs_ham%ind_mass_OIFS))  DEALLOCATE(ind_oifs_ham%ind_mass_OIFS)
    IF (ALLOCATED(ind_oifs_ham%ind_mass_HAM))   DEALLOCATE(ind_oifs_ham%ind_mass_HAM)
    IF (ALLOCATED(ind_oifs_ham%ind_gas_OIFS))   DEALLOCATE(ind_oifs_ham%ind_gas_OIFS)
    IF (ALLOCATED(ind_oifs_ham%ind_gas_HAM))    DEALLOCATE(ind_oifs_ham%ind_gas_HAM)
    IF (ALLOCATED(ind_oifs_ham%ind_cloud_OIFS)) DEALLOCATE(ind_oifs_ham%ind_cloud_OIFS)
    IF (ALLOCATED(ind_oifs_ham%ind_cloud_HAM))  DEALLOCATE(ind_oifs_ham%ind_cloud_HAM)

    ALLOCATE(ind_oifs_ham%ind_class_OIFS(knclass))
    ALLOCATE(ind_oifs_ham%ind_class_HAM(knclass))
    ALLOCATE(ind_oifs_ham%ind_mass_OIFS(knaerocomp))
    ALLOCATE(ind_oifs_ham%ind_mass_HAM(knaerocomp))
    ALLOCATE(ind_oifs_ham%ind_gas_OIFS(ksubm_ngasspec))
    ALLOCATE(ind_oifs_ham%ind_gas_HAM(ksubm_ngasspec))
    ALLOCATE(ind_oifs_ham%ind_cloud_OIFS(kcloudind))
    ALLOCATE(ind_oifs_ham%ind_cloud_HAM(kcloudind))

    ind_oifs_ham%ind_class_OIFS(:) = 0
    ind_oifs_ham%ind_mass_OIFS(:) = 0
    ind_oifs_ham%ind_gas_OIFS(:) = 0
    ind_oifs_ham%ind_cloud_OIFS(:) = 0
    ind_oifs_ham%ind_class_HAM(:) = 0
    ind_oifs_ham%ind_mass_HAM(:) = 0
    ind_oifs_ham%ind_gas_HAM(:) = 0
    ind_oifs_ham%ind_cloud_HAM(:) = 0
    
  END SUBROUTINE init_ind_oifs_ham
  
END MODULE OIFS_to_HAM
