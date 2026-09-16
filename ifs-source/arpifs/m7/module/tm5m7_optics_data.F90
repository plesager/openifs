MODULE TM5M7_OPTICS_DATA

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE TM5M7_DATA, ONLY: NMOD

IMPLICIT NONE

SAVE

  ! wavelength type (to be used by methods using the optics)
  type, public :: WAVELENDEP 
     real  :: wl             ! user requested wavelength    unit = um (e.g. 0.550)
     real, dimension(7) :: n ! SO4, BC, OC, SOA, SS, DU, WATER
     real, dimension(7) :: k ! SO4, BC, OC, SOA, SS, DU, WATER
     logical :: split, insitu
  end type WAVELENDEP

  ! AOP input type and field
  type aopi
     REAL(KIND=JPRB), dimension(nmod) :: SO4, NO3, BC, OC, SOA, SS, DU, H2O, numdens, rg, rgd
  end type aopi


 
 ! Characteristics of the lookup-table
 INTEGER(KIND=JPIM), PARAMETER  :: N_RII=15
 INTEGER(KIND=JPIM), PARAMETER  :: N_RIR=40
 INTEGER(KIND=JPIM), PARAMETER  :: N_X=100
 
 REAL(KIND=JPRB),DIMENSION(N_RII)   :: lkval        ! -log img part refr. index
 REAL(KIND=JPRB),DIMENSION(N_RII)   :: kval         ! img part refr. index, 10^(-lkval)
 REAL(KIND=JPRB),DIMENSION(N_RIR)   :: n1r          ! real part refr. index
 
 REAL(KIND=JPRB),DIMENSION(N_X)   :: XS
 REAL(KIND=JPRB),DIMENSION(N_X,N_RIR,N_RII),TARGET :: CEXT_159,A_159,G_159
 REAL(KIND=JPRB),DIMENSION(N_X,N_RIR,N_RII),TARGET :: CEXT_200,A_200,G_200
 
 INTEGER(KIND=JPIM), PARAMETER :: OPACDIM = 61
 INTEGER(KIND=JPIM), PARAMETER :: ECHAMHAMDIM=49
 INTEGER(KIND=JPIM), PARAMETER :: SEGELSTEINDIM=1261
 REAL(KIND=JPRB),DIMENSION(:,:), ALLOCATABLE:: OPAC,ECHAMHAM,SEGELSTEIN

 INTEGER(KIND=JPIM) :: NWDEP
 TYPE(WAVELENDEP),DIMENSION(:),ALLOCATABLE :: WDEP
 

 INTEGER(KIND=JPIM) :: NASWBAND
 TYPE(WAVELENDEP),DIMENSION(:),ALLOCATABLE :: ASWBAND

 INTEGER(KIND=JPIM), PARAMETER  :: NALWBAND =16
 REAL(KIND=JPRB),DIMENSION(NALWBAND)   ::  ALWWN1
 REAL(KIND=JPRB),DIMENSION(NALWBAND)   ::  ALWWN2
 ! dimension for array AOP_OUT_ADD 
 INTEGER(KIND=JPIM), PARAMETER :: NADD = 5
 
END MODULE TM5M7_OPTICS_DATA
