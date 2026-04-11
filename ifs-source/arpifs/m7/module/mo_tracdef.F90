MODULE mo_tracdef


  USE mo_kind,           ONLY: dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: trlist                         ! tracer info list variable                       
  PUBLIC :: t_trlist, t_trinfo, t_flag     ! "" type definitions
  PUBLIC :: ln, ll, lf, nf                 ! length of names in trlist

  PUBLIC :: jptrac                         ! max. number of tracers
  PUBLIC :: ntrac                          ! number of tracers actually defined (*)
                                           ! Note (*): ntrac always equal trlist% ntrac. 
                                           ! It's just more convenient.

  !-- flag values
  !   general flag values
  PUBLIC :: OK                             
  PUBLIC :: ON                             
  PUBLIC :: OFF
  !   initialization mode 'ninit'
  PUBLIC :: INITIAL
  PUBLIC :: RESTART
  PUBLIC :: CONSTANT
  PUBLIC :: PERIODIC
  !   'nsoluble'
  PUBLIC :: SOLUBLE
  PUBLIC :: INSOLUBLE
  !   'nphase'
  PUBLIC :: GAS
  PUBLIC :: AEROSOL
  PUBLIC :: GAS_OR_AEROSOL
  PUBLIC :: AEROSOLMASS
  PUBLIC :: AEROSOLNUMBER
  PUBLIC :: UNDEFINED
  !   choice of advection scheme 'ntran'
  !   tracer type 'itrtype'
  PUBLIC :: ITRNONE
  PUBLIC :: ITRPRESC
  PUBLIC :: ITRDIAG
  PUBLIC :: ITRPROG
  ! Limits
  !
  INTEGER, PARAMETER :: jptrac = 500  ! maximum number of tracers allowed
  !
  ! Individual settings for each tracer
  !
  INTEGER, PARAMETER :: ln =  24  ! length of name (char) components  
  INTEGER, PARAMETER :: ll = 256  ! length of longname and standardname 
  INTEGER, PARAMETER :: lf =   8  ! length of flag character string
  INTEGER, PARAMETER :: nf =  10  ! number of user defined flags
  INTEGER, PARAMETER :: ns =  20  ! max number of submodels
 
  !
  ! Constants (argument values for new_tracer routine in mo_tracer)
  !
  !
  ! error return values
  !
  INTEGER, PARAMETER :: OK         = 0
  INTEGER, PARAMETER :: NAME_USED  = 2
  INTEGER, PARAMETER :: NAME_MISS  = 3
  INTEGER, PARAMETER :: TABLE_FULL = 4
  !
  ! general flags
  !
  INTEGER, PARAMETER :: OFF        = 0
  INTEGER, PARAMETER :: ON         = 1
  !
  ! initialisation flag  (ninit)
  !
  INTEGER, PARAMETER :: CONSTANT = 1
  INTEGER, PARAMETER :: RESTART  = 2
  INTEGER, PARAMETER :: INITIAL  = 4
  INTEGER, PARAMETER :: PERIODIC = 8
  ! 
  ! Tracer type
  ! 
  INTEGER, PARAMETER :: ITRNONE    = 0   ! No ECHAM tracer
  INTEGER, PARAMETER :: ITRPRESC   = 1   ! Prescribed tracer, read from file
  INTEGER, PARAMETER :: ITRDIAG    = 2   ! Diagnostic tracer, no transport
  INTEGER, PARAMETER :: ITRPROG    = 3   ! Prognostic tracer, transport by ECHAM
  !
  ! soluble flag  (nsoluble)
  !
  INTEGER, PARAMETER :: INSOLUBLE = 0 ! insoluble
  INTEGER, PARAMETER :: SOLUBLE   = 1 ! soluble
  !
  ! phase indicator  (nphase)
  !
  INTEGER, PARAMETER :: GAS            = 1  ! gas
  INTEGER, PARAMETER :: AEROSOL        = 2  ! aerosol (for species definition)
  INTEGER, PARAMETER :: GAS_OR_AEROSOL = 3  ! gas or aerosol (for species definition)
  INTEGER, PARAMETER :: AEROSOLMASS    = 2  ! aerosol mass
  INTEGER, PARAMETER :: AEROSOLNUMBER  = -2 ! particle number concentration
  INTEGER, PARAMETER :: UNDEFINED      = 0  ! other tracers (e.g. CDNC)

  !
  ! General purpose flag (additional tracer property)
  !
  TYPE t_flag
    CHARACTER(len=lf) :: c      ! character string
    REAL(dp)          :: v      ! value
  END TYPE t_flag
 
  !
  ! Individual settings for each tracer
  ! Default values are marked with * *
  !
 !=============!
  TYPE t_trinfo
 !=============!
    !
    ! identification of transported quantity
    !
    CHARACTER(len=ln) :: basename   ! name (instead of xt..)
    CHARACTER(len=ln) :: subname    ! optional for 'colored' tracer
    CHARACTER(len=ln) :: fullname   ! name_subname
    CHARACTER(len=ln) :: modulename ! name of requesting sub-model
    CHARACTER(len=ln) :: units      ! units
    CHARACTER(len=ll) :: longname   ! long name 
    CHARACTER(len=ll) :: standardname   ! CF standard name
    INTEGER           :: trtype     ! type of tracer: 0=undef., 1=prescribed, 2=diagnostic (no transport),
                                    !                 *3*=prognostic (transported)
    INTEGER           :: spid       ! species id (index in speclist) where physical/chemical 
                                    ! properties are defined 
    INTEGER           :: nphase     ! phase (1=GAS, 2=AEROSOLMASS, 3=AEROSOLNUMBER, ...??)
                                    ! [add liquid or ice phase ??]
    INTEGER           :: mode       ! aerosol mode or bin number  (default 0)
    REAL(dp)          :: moleweight ! molecular mass (copied from species upon initialisation)
!   INTEGER           :: tag        ! tag of requesting routine 

    !
    ! Requested resources ...
    ! 
    INTEGER           :: burdenid   ! index in burden diagnostics
    INTEGER           :: nbudget    ! calculate budgets        (default 0)
    INTEGER           :: ntran      ! perform transport        (default 1)
    INTEGER           :: nfixtyp    ! type of mass fixer for semi lagrangian adv. (default 1)
    INTEGER           :: nconvmassfix ! use xt_conv_massfix in cumastr
    INTEGER           :: nvdiff     ! vertical diffusion flag  (default 1)
    INTEGER           :: nconv      ! convection flag          (default 1) 
    INTEGER           :: ndrydep    ! dry deposition flag: *0*=no drydep, 1=prescribed vd,
                                    !                       2=Ganzeveld 
    INTEGER           :: nwetdep    ! wet deposition flag      (default 0)
    INTEGER           :: nsedi      ! sedimentation flag       (default 0) 
    REAL(dp)          :: tdecay     ! decay time (exponential) (default 0.sec)
    INTEGER           :: nemis      ! surface emission flag    (default 0) 
                                    ! emission flag:  *0*=no emissions, 
                                    !                  1=surface flux cond.,
                                    !                  2=tendency (2D emis.) [additive] 
!   INTEGER           :: nint       ! integration flag         (default 1)  
    
    !
    ! initialization and restart
    !
    INTEGER           :: ninit      ! initialization request flag
    INTEGER           :: nrerun     ! rerun flag
    REAL(dp)          :: vini       ! initialisation value     (default 0.) 
    INTEGER           :: init       ! initialisation method actually used 
    !
    ! Flags used for postprocessing
    !
    INTEGER           :: nwrite     ! write flag            (default 1)
    INTEGER           :: code       ! tracer code,          (default 235...)
    INTEGER           :: table      ! tracer code table     (default 0)
    INTEGER           :: gribbits   ! bits for encoding     (default 16) 
    INTEGER           :: nint       ! integration (accumulation) flag  (default 1)  
    
    !
    ! Flags to be used by chemistry or tracer modules
    ! 
!### henry constant will be removed from trlist -- available as species property instead
!!mgs!!    REAL(dp)          :: henry      ! Henry coefficient [?] (default 1.e-10)
!!mgs!!    REAL(dp)          :: dryreac    ! dry reactivity coeff. (default 0.)
    INTEGER           :: nsoluble   ! soluble flag          (default 0) 
    
    TYPE(t_flag)      :: myflag (nf)! user defined flag
    !
    ! Indicate actions actually performed by ECHAM
    ! 
    
  END TYPE t_trinfo

  
  !
  ! Reference to memory buffer information for each tracer
  ! used to access the 'restart' flags
  !

  !
  ! Basic data type definition for tracer info list
  !
 !=============!
  TYPE t_trlist
 !=============!
    !
    ! global tracer list information
    !
    INTEGER         :: ntrac        ! number of tracers specified
    INTEGER         :: anyfixtyp    ! mass fixer types used
    INTEGER         :: anywetdep    ! wet deposition requested for any tracer
    INTEGER         :: anydrydep    ! wet deposition requested for any tracer
    INTEGER         :: anysedi      ! sedimentation  requested for any tracer
    INTEGER         :: anysemis     ! surface emission flag for any tracer
    INTEGER         :: anyconv      ! convection flag
    INTEGER         :: anyvdiff     ! vertical diffusion flag
    INTEGER         :: anyconvmassfix  ! 
    INTEGER         :: nadvec       ! number of advected tracers
    LOGICAL         :: oldrestart   ! true to read old restart format
    !
    ! individual information for each tracer
    !
    TYPE (t_trinfo) :: ti  (jptrac) ! Individual settings for each tracer
    !
    ! reference to memory buffer info
    !
  END TYPE t_trlist
  !
  ! module variables
  !

  TYPE(t_trlist)   ,SAVE ,TARGET :: trlist        ! tracer list 
  INTEGER          ,SAVE         :: ntrac = 0     ! number of tracers actually defined
 
END MODULE mo_tracdef

