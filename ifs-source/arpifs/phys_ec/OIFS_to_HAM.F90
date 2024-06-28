MODULE OIFS_to_HAM

  USE mo_ham,   ONLY: nclass, naerocomp, subm_ngasspec

  IMPLICIT NONE

  PUBLIC :: init_ind_oifs_ham
  
! ╒════════════════════════════════════════════════════════════════════════════╕
! │ MODULE OFS_to_HAM                                    (updated 22-MAY-2024) │
! │                                                                            │
! │  Contains all the variables for needed to make HAM compatible with OIFS.   │
! │                                                                            │
! │ TYPES:                                                                     │
! │ - TYPE_GFL_COMP   -> ind_oifs_ham_type                                     │
! │                                                                            │
! │                                                                            │
! │ Author : Eemeli Holopainen (FMI)   eemeli.holopainen@fmi.fi                │
! │ -------                                                                    │
! │                                                                            │
! │ Modifications :                                                            │
! │ -------------                                                              │
! │                                                                            │
! ╘════════════════════════════════════════════════════════════════════════════╛

  ! Index list for HAM and OIFS tracers
  type  ind_oifs_ham_type
      INTEGER, PUBLIC, ALLOCATABLE :: ind_class_OIFS(:)    ! indices of sizeclass OIFS tracers
      INTEGER, PUBLIC, ALLOCATABLE :: ind_mass_OIFS(:)     ! indices of mass OIFS tracers
      INTEGER, PUBLIC, ALLOCATABLE :: ind_gas_OIFS(:)      ! indices of gas OIFS tracers
      INTEGER, PUBLIC, ALLOCATABLE :: ind_cloud_OIFS(:)    ! indices of cloud OIFS tracers
      INTEGER, PUBLIC, ALLOCATABLE :: ind_class_HAM(:)     ! indices of sizeclass HAM tracers
      INTEGER, PUBLIC, ALLOCATABLE :: ind_mass_HAM(:)      ! indices of mass HAM tracers
      INTEGER, PUBLIC, ALLOCATABLE :: ind_gas_HAM(:)       ! indices of gas HAM tracers
      INTEGER, PUBLIC, ALLOCATABLE :: ind_cloud_HAM(:)     ! indices of cloud HAM tracers
  end type ind_oifs_ham_type

  TYPE(ind_oifs_ham_type) :: ind_oifs_ham
  !!$OMP THREADPRIVATE(ind_oifs_ham)

CONTAINS

  SUBROUTINE init_ind_oifs_ham(knclass, knaerocomp, ksubm_ngasspec, kcloudind)
    
    ! *ind_oif_ham* allocates and initializes the index list for OIFS tracers and HAM tracers
    ! Authors:
    ! -------
    ! Eemeli Holopainen, FMI                4/2022

    INTEGER, INTENT(IN) :: knclass, knaerocomp, ksubm_ngasspec, kcloudind
    
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
    ind_oifs_ham%ind_mass_OIFS(:)  = 0
    ind_oifs_ham%ind_gas_OIFS(:)   = 0
    ind_oifs_ham%ind_cloud_OIFS(:) = 0
    ind_oifs_ham%ind_class_HAM(:)  = 0
    ind_oifs_ham%ind_mass_HAM(:)   = 0
    ind_oifs_ham%ind_gas_HAM(:)    = 0
    ind_oifs_ham%ind_cloud_HAM(:)  = 0
    
  END SUBROUTINE init_ind_oifs_ham
  
END MODULE OIFS_to_HAM
