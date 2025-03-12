 SUBROUTINE TM5M7_INIT(YDGEOMETRY, YRCOMPO, YGFL, YDERAD)

!**   DESCRIPTION 
!     ----------
!
!   init routine for IFS tm5m7 aerosol 
!
!
!
!**   INTERFACE.
!     ----------
!          *TM5M7_INIT* IS CALLED FROM *CNT4*.
!

!     Externals.  
!     ---------                  

!
!     AUTHOR.
!     -------
!       Vincent Huijnen *KNMI*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2020-08-24

USE GEOMETRY_MOD,   ONLY : GEOMETRY
USE PARKIND1,       ONLY : JPRB, JPIM
USE YOMHOOK,        ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN,         ONLY : NULOUT
!USE YOMCOMPO , ONLY : YRCOMPO
USE YOMCOMPO,       ONLY : TCOMPO
USE YOM_YGFL,       ONLY : TYPE_GFLD!YGFL
USE TM5_PHOTOLYSIS, ONLY : PHOTOLYSIS_INI, NBANDS_TROP,LMID,LMID_GRIDA,WAVE,WAV_GRID,WAV_GRIDA,LL_TM5_PHOTO_INI
USE TM5M7_DATA,     ONLY : ISO4 ,  INH4 ,  INO3_A ,  IACS_N ,  ISO4ACS ,  IBCACS , IPOMACS ,  ISSACS ,  IDUACS , &
                         & ISOANUS ,  ISOAAIS ,  ISOAACS ,  ISOACOS ,  ISOAAII ,  IH2OPART ,IAII_N ,  IBCAII , &
                         & IPOMAII ,  IACI_N ,   IDUACI ,  IAIS_N ,  ISO4AIS ,  IBCAIS ,  IPOMAIS , ICOI_N , &
                         & IDUCOI  ,  ICOS_N ,  ISO4COS ,  IBCCOS ,  IPOMCOS ,  ISSCOS ,  IDUCOS ,  INUS_N , &
                         & ISO4NUS ,  IELVOC ,  IISVOC ,  IMSA, &   
                         ! Needed for the emissions declaration..
                         & sigma_lognormal, &
                         & xmc, sigma_lognormal, pom_density, &
                         & mode_aii, mode_ais, mode_acs
USE TM5M7_OPTICS_DATA, ONLY :WAVELENDEP,NWDEP,WDEP, & 
  & ASWBAND, NASWBAND,ALWWN1,ALWWN2
USE YOERAD   , ONLY : TERAD!YRERAD
!USE YOMPRAD  , ONLY : RADGRID

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TCOMPO)      ,INTENT(IN)    :: YRCOMPO
TYPE(TYPE_GFLD)   ,INTENT(IN)    :: YGFL
TYPE(TERAD),INTENT(IN) :: YDERAD
!-----------------------------------------------------------------------
!*       0.5   LOCAL VARIABLES
!              ---------------
INTEGER(KIND=JPIM) :: JI,JL,JK
LOGICAL            :: LLFOUND
REAL(KIND=JPRB),DIMENSION(:),ALLOCATABLE :: PHOTO_WAVELENGTHS

REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
#include "abor1.intfb.h"
#include "tm5m7_src_dust_init.intfb.h"
#include "tm5m7_optics_init.intfb.h"


IF (LHOOK) CALL DR_HOOK('TM5M7_INIT',0,ZHOOK_HANDLE)
ASSOCIATE(&
    & NACTAERO    => YGFL%NACTAERO, YAERO    => YGFL%YAERO,    &
    & NAERO       => YGFL%NAERO,    LAERCHEM => YGFL%LAERCHEM, &
    & AERO_SCHEME => YRCOMPO%AERO_SCHEME )

 
!*             Init aerosol scheme 
!              ---------------

SELECT CASE (TRIM(AERO_SCHEME))

  CASE ("aer")         
     
    ! Setup of 'aer' configuration is done in su_aerw.F90
    IF (LHOOK) CALL DR_HOOK('TM5M7_INIT',1,ZHOOK_HANDLE)    
    RETURN
             
  CASE ("hamm7")
    
    ! Initialization of tm5m7 aerosol tracer indices is handled in hamm7_init 

    ! Initialize various dust properties
    CALL TM5M7_SRC_DUST_INIT

    !IF(.not.LAERCHEM)THEN
      ! Initialize optics:
      ! Make sure that 'WAVE' is already initialized (in tm5_init.F90)
      !IF (.NOT. LL_TM5_PHOTO_INI) THEN
      !   if (.NOT. LAERCHEM)THEN
      !      call PHOTOLYSIS_INI
      !   ELSE
      !      CALL ABOR1('tm5-based photolysis not yet initialized!!')
      !   END if
     
      !ENDIF
    !END IF

    ! Define wavelengths for optics calculations
    nwdep = nbands_trop + count(lmid.ne.lmid_gridA)
    wav_grid  = 0
    wav_gridA = 0
    allocate(photo_wavelengths(nwdep))

    JL=1
    do JI=1,nbands_trop
      if (lmid(JI)==lmid_gridA(JI)) then
        photo_wavelengths(JL) = wave(lmid(JI))*1.e4 ! cm to um
        wav_grid(JI) = JL
        wav_gridA(JI) = JL
        JL=JL+1   
      else
        photo_wavelengths(JL) = wave(lmid(JI))*1.e4 ! cm to um
        photo_wavelengths(JL+1) = wave(lmid_gridA(JI))*1.e4 ! cm to um
        wav_grid(JI) = JL
        wav_gridA(JI) = JL+1
        JL=JL+2
      endif
    enddo
    allocate(wdep(nwdep))
    wdep(:)%wl = photo_wavelengths
    wdep(:)%split = .false.
    wdep(:)%insitu = .false.

    CALL TM5M7_OPTICS_INIT(NWDEP,WDEP)

    if (allocated(photo_wavelengths)) deallocate(photo_wavelengths)
    if (allocated(wdep)) deallocate(wdep)

    ! A.Laakso: Taken from ecearth_optics (TM5-ECEARTH3) 
    ! HAM aerosol optics are using these too
    NASWBAND=YDERAD%NTSW
    if (allocated(ASWBAND)) deallocate(ASWBAND)
    allocate(ASWBAND(YDERAD%NTSW))
    ASWBAND(13)%wl = 0.257_JPRB
    ASWBAND(12)%wl = 0.313_JPRB
    ASWBAND(11)%wl = 0.398_JPRB
    ASWBAND(10)%wl = 0.530_JPRB
    ASWBAND( 9)%wl = 0.697_JPRB
    ASWBAND( 8)%wl = 0.973_JPRB
    ASWBAND( 7)%wl = 1.269_JPRB
    ASWBAND( 6)%wl = 1.447_JPRB
    ASWBAND( 5)%wl = 1.767_JPRB
    ASWBAND( 4)%wl = 2.040_JPRB
    ASWBAND( 3)%wl = 2.308_JPRB
    ASWBAND( 2)%wl = 2.752_JPRB
    ASWBAND( 1)%wl = 3.407_JPRB
    ASWBAND(14)%wl = 5.254_JPRB

    ASWBAND(:)%split = .false. 
    ASWBAND(:)%insitu = .false. 

    CALL TM5M7_OPTICS_INIT(NASWBAND,ASWBAND)
    
    !LW wavenumbers for ham optics
    ALWWN1 = (/ & !< Spectral band lower boundary in wavenumbers
    &   10._JPRB, 350._JPRB, 500._JPRB, 630._JPRB, 700._JPRB, 820._JPRB, &
    &  980._JPRB,1080._JPRB,1180._JPRB,1390._JPRB,1480._JPRB,1800._JPRB, &
    & 2080._JPRB,2250._JPRB,2380._JPRB,2600._JPRB/)
    ALWWN2 = (/ & !< Spectral band upper boundary in wavenumbers
    &  350._JPRB, 500._JPRB, 630._JPRB, 700._JPRB, 820._JPRB, 980._JPRB, &
    & 1080._JPRB,1180._JPRB,1390._JPRB,1480._JPRB,1800._JPRB,2080._JPRB, &
    & 2250._JPRB,2380._JPRB,2600._JPRB,3250._JPRB/)

  CASE DEFAULT
  
    ! Option not implemented
    CALL ABOR1(" NO AEROSOL SCHEME "//TRIM(AERO_SCHEME))

END SELECT 

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('TM5M7_INIT',1,ZHOOK_HANDLE)
END SUBROUTINE TM5M7_INIT 
