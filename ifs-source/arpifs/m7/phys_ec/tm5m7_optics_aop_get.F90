SUBROUTINE TM5M7_OPTICS_AOP_GET(YGFL, YREAERSRC, KIDIA,KFDIA,KLON, KLEV, NACTAERO, &
          &  nwav, wdep, ncontr, ecearth_units, &
          &  PRHO, PAERO, RW_MODE,RWD_MODE,H2O_MODE, &
          &  Paop_out_ext, Paop_out_a, Paop_out_g, aop_out_add )

!*** * TM5M7_OPTICS_AOP_GET* 
!
!
!-----------------------------------------------------------------------------
!                    TM5                                                     !
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE:	OPTICS_AOP_GET
!
! !DESCRIPTION: Initialise the fields "aop_in" and then calculate the
!		optical properties through a call to optics_calculate_aop. 
!\\
!\\
  
!     SOURCE.
!     -------
!
!     Taken from TM5-code
!
!     MODIFICATIONS.
!     --------------
!
!      Sep 2021 - V. Huijnen: first introduction into OpenIFS
!
!
!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMLUN    ,ONLY : NULOUT
USE YOEAERSRC, ONLY : TEAERSRC !YREAERSRC
!USE YOMCST, ONLY : 
USE TM5M7_DATA, ONLY : H2SO4_FACTOR,NH4NO3_FACTOR, SO4_DENSITY,MSA_DENSITY, &
  & NH4NO3_DENSITY, &
  & MODAL_DATA,NMOD,NSOL
USE TM5M7_OPTICS_DATA, ONLY: WAVELENDEP, AOPI,NADD
USE YOM_YGFL , ONLY : TYPE_GFLD!YGFL
USE TM5M7_DATA, ONLY : &
  & INO3_A  ,  IACS_N ,  ISO4ACS ,  IBCACS ,  IPOMACS ,  ISSACS ,  IDUACS , &
  & ISOANUS ,  ISOAAIS ,  ISOAACS ,  ISOACOS ,  ISOAAII ,  IH2OPART ,IAII_N ,  IBCAII , &
  & IPOMAII ,  IACI_N ,   IDUACI ,  IAIS_N ,  ISO4AIS ,  IBCAIS ,  IPOMAIS , ICOI_N , &
  & IDUCOI  ,  ICOS_N ,  ISO4COS ,  IBCCOS ,  IPOMCOS ,  ISSCOS ,  IDUCOS ,  INUS_N , &
  & ISO4NUS ,  IELVOC ,  IISVOC ,  IMSA 


IMPLICIT NONE

!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

TYPE(TEAERSRC), INTENT(IN) :: YREAERSRC
TYPE(TYPE_GFLD)   ,INTENT(IN) :: YGFL
INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA,KFDIA,KLON,KLEV,NACTAERO, NWAV,NCONTR
TYPE(WAVELENDEP), DIMENSION(NWAV), INTENT(IN):: WDEP
LOGICAL,         INTENT(IN) :: ECEARTH_UNITS
REAL(KIND=JPRB), INTENT(IN) :: PRHO(KLON,KLEV)     ! air Density [kg/m3]
REAL(KIND=JPRB), INTENT(IN) :: PAERO(KLON,KLEV,NACTAERO)


TYPE(MODAL_DATA), INTENT(IN) :: RW_MODE(NMOD)
TYPE(MODAL_DATA), INTENT(IN) :: RWD_MODE(NSOL)
TYPE(MODAL_DATA), INTENT(IN) :: H2O_MODE(NSOL)

REAL(KIND=JPRB), INTENT(OUT):: PAOP_OUT_EXT(KLON,KLEV,NWAV,NCONTR)
REAL(KIND=JPRB), INTENT(OUT):: PAOP_OUT_A(KLON,KLEV,NWAV)
REAL(KIND=JPRB), INTENT(OUT):: PAOP_OUT_G(KLON,KLEV,NWAV)
REAL(KIND=JPRB), INTENT(OUT), OPTIONAL :: AOP_OUT_ADD(KLON,KLEV,NWAV,NADD) ! additional parameters


!*       0.2   LOCAL VARIABLES
!              ---------------

TYPE(AOPI), dimension(:,:), allocatable :: aop_in 
INTEGER(KIND=JPRB) :: IMODE

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
INTEGER(kind=JPIM)::NCHEM
    
!-----------------------------------------------------------------------

#include "tm5m7_optics_calculate_aop.intfb.h"
associate(NCHEM => YGFL%NCHEM)
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TM5M7_OPTICS_AOP_GET',0,ZHOOK_HANDLE)
  


    ! --- start ------------------------------


    allocate( aop_in(KLON,KLEV) )

    ! Initialize full array to zero.
    DO IMODE = 1, nmod
       aop_in(1:KLON,1:KLEV)%so4    (IMODE) = 0.0_JPRB ; aop_in(1:KLON,1:KLEV)%bc   (IMODE) = 0.0_JPRB
       aop_in(1:KLON,1:KLEV)%oc     (IMODE) = 0.0_JPRB ; aop_in(1:KLON,1:KLEV)%soa  (IMODE) = 0.0_JPRB
       aop_in(1:KLON,1:KLEV)%ss     (IMODE) = 0.0_JPRB
       aop_in(1:KLON,1:KLEV)%du     (IMODE) = 0.0_JPRB ; aop_in(1:KLON,1:KLEV)%h2o  (IMODE) = 0.0_JPRB
       aop_in(1:KLON,1:KLEV)%numdens(IMODE) = 0.0_JPRB ; aop_in(1:KLON,1:KLEV)%rg   (IMODE) = 0.0_JPRB
       aop_in(1:KLON,1:KLEV)%rgd    (IMODE) = 0.0_JPRB ; aop_in(1:KLON,1:KLEV)%no3  (IMODE) = 0.0_JPRB
    END DO
    
!>>> TvN
! In M7 sulphate is assumed to be H2-SO4 with corresponding particle density so4_density
! The sulphate mass should therefore also include the small contribution from the H atoms
    ! NUS
    aop_in(KIDIA:KFDIA,1:KLEV)%so4(1) = 1.E9_JPRB * h2so4_factor * PAERO(KIDIA:KFDIA,1:KLEV,iso4nus) 
    aop_in(KIDIA:KFDIA,1:KLEV)%soa(1) = 1.E9_JPRB * PAERO(KIDIA:KFDIA,1:KLEV,isoanus) 

    ! AIS 
    aop_in(KIDIA:KFDIA,1:KLEV)%so4(2) = 1.E9_JPRB * h2so4_factor * PAERO(KIDIA:KFDIA,1:KLEV,iso4ais) 
    aop_in(KIDIA:KFDIA,1:KLEV)%bc (2) = 1.E9_JPRB * PAERO(KIDIA:KFDIA,1:KLEV,ibcais ) 
    aop_in(KIDIA:KFDIA,1:KLEV)%oc (2) = 1.E9_JPRB * PAERO(KIDIA:KFDIA,1:KLEV,ipomais) 
    aop_in(KIDIA:KFDIA,1:KLEV)%soa(2) = 1.E9_JPRB * PAERO(KIDIA:KFDIA,1:KLEV,isoaais) 
    ! ACS (additional: NO3)
! The contribution from methane sulfonate (MSA-) aerosol is added
! to that for sulfate.
! As the addition is done by volume,
! we need to account for the difference in densities
! (as done below for ammonium nitrate).
    if(NCHEM<1) then
       aop_in(KIDIA:KFDIA,1:KLEV)%so4(3) = 1.E9_JPRB * ( h2so4_factor * PAERO(KIDIA:KFDIA,1:KLEV,iso4acs))
    ELSE
       aop_in(KIDIA:KFDIA,1:KLEV)%so4(3) = 1.E9_JPRB * ( h2so4_factor * PAERO(KIDIA:KFDIA,1:KLEV,iso4acs) + &
            (so4_density / msa_density) * PAERO(KIDIA:KFDIA,1:KLEV,imsa) ) 
    END if
! Since nh4no3_density is the density of NH4NO3, the contribution from NH4 should be included.
! Moreover, assuming the same refractive index for NH4NO3 as for H2-SO4,
! the contributions from both components can be added by volume;
! thus we need to account for the difference in densities.
! Estimates of the refractive index of NH4NO3 are available from literature
! (e.g. Lowenthal et al., Atmos. Environ., 2000).
! For practical purposes, it can be set equal to the value used for sulfate,
! i.e. the value for a solution containing 75% H2SO4 (Fenn et al., 1985).
    if(NCHEM>0) then
       aop_in(KIDIA:KFDIA,1:KLEV)%no3(3) = 1.E9_JPRB * nh4no3_factor * (so4_density / nh4no3_density) * &
                                       PAERO(KIDIA:KFDIA,1:KLEV,ino3_a ) 
    END if
       aop_in(KIDIA:KFDIA,1:KLEV)%bc (3) =  1.E9_JPRB * PAERO(KIDIA:KFDIA,1:KLEV,ibcacs ) 
    aop_in(KIDIA:KFDIA,1:KLEV)%oc (3) =  1.E9_JPRB * PAERO(KIDIA:KFDIA,1:KLEV,ipomacs) 
    aop_in(KIDIA:KFDIA,1:KLEV)%soa(3) =  1.E9_JPRB * PAERO(KIDIA:KFDIA,1:KLEV,isoaacs) 
    aop_in(KIDIA:KFDIA,1:KLEV)%ss (3) =  1.E9_JPRB * PAERO(KIDIA:KFDIA,1:KLEV,issacs ) 
    aop_in(KIDIA:KFDIA,1:KLEV)%du (3) =  1.E9_JPRB * PAERO(KIDIA:KFDIA,1:KLEV,iduacs ) 
    ! COS
    aop_in(KIDIA:KFDIA,1:KLEV)%so4(4) =  1.E9_JPRB * h2so4_factor * PAERO(KIDIA:KFDIA,1:KLEV,iso4cos) 
!<<< TvN
    aop_in(KIDIA:KFDIA,1:KLEV)%bc (4) =  1.E9_JPRB * PAERO(KIDIA:KFDIA,1:KLEV,ibccos )
    aop_in(KIDIA:KFDIA,1:KLEV)%oc (4) =  1.E9_JPRB * PAERO(KIDIA:KFDIA,1:KLEV,ipomcos)
    aop_in(KIDIA:KFDIA,1:KLEV)%soa(4) =  1.E9_JPRB * PAERO(KIDIA:KFDIA,1:KLEV,isoacos)
    aop_in(KIDIA:KFDIA,1:KLEV)%ss (4) =  1.E9_JPRB * PAERO(KIDIA:KFDIA,1:KLEV,isscos )
    aop_in(KIDIA:KFDIA,1:KLEV)%du (4) =  1.E9_JPRB * PAERO(KIDIA:KFDIA,1:KLEV,iducos )
    ! AII
    aop_in(KIDIA:KFDIA,1:KLEV)%bc (5) =  1.E9_JPRB * PAERO(KIDIA:KFDIA,1:KLEV,ibcaii ) 
    aop_in(KIDIA:KFDIA,1:KLEV)%oc (5) =  1.E9_JPRB * PAERO(KIDIA:KFDIA,1:KLEV,ipomaii) 
    aop_in(KIDIA:KFDIA,1:KLEV)%soa(5) =  1.E9_JPRB * PAERO(KIDIA:KFDIA,1:KLEV,isoaaii) 
    ! ACI
    aop_in(KIDIA:KFDIA,1:KLEV)%du (6) = 1.E9_JPRB * PAERO(KIDIA:KFDIA,1:KLEV,iduaci ) 
    ! COI
    aop_in(KIDIA:KFDIA,1:KLEV)%du (7) = 1.E9_JPRB * PAERO(KIDIA:KFDIA,1:KLEV,iducoi ) 
    ! Water in (hydrophillic) modes
    aop_in(KIDIA:KFDIA,1:KLEV)%h2o(1) = 1.E9_JPRB * h2o_mode(1)%d2(KIDIA:KFDIA,1:KLEV)/PRHO(KIDIA:KFDIA,1:KLEV)
    aop_in(KIDIA:KFDIA,1:KLEV)%h2o(2) = 1.E9_JPRB * h2o_mode(2)%d2(KIDIA:KFDIA,1:KLEV)/PRHO(KIDIA:KFDIA,1:KLEV) 
    aop_in(KIDIA:KFDIA,1:KLEV)%h2o(3) = 1.E9_JPRB * h2o_mode(3)%d2(KIDIA:KFDIA,1:KLEV)/PRHO(KIDIA:KFDIA,1:KLEV) 
    aop_in(KIDIA:KFDIA,1:KLEV)%h2o(4) = 1.E9_JPRB * h2o_mode(4)%d2(KIDIA:KFDIA,1:KLEV)/PRHO(KIDIA:KFDIA,1:KLEV) 

    aop_in(KIDIA:KFDIA,1:KLEV)%rg (1) = 1.E6_JPRB * rw_mode (1)%d2(KIDIA:KFDIA,1:KLEV)
    aop_in(KIDIA:KFDIA,1:KLEV)%rg (2) = 1.E6_JPRB * rw_mode (2)%d2(KIDIA:KFDIA,1:KLEV)
    aop_in(KIDIA:KFDIA,1:KLEV)%rg (3) = 1.E6_JPRB * rw_mode (3)%d2(KIDIA:KFDIA,1:KLEV)
    aop_in(KIDIA:KFDIA,1:KLEV)%rg (4) = 1.E6_JPRB * rw_mode (4)%d2(KIDIA:KFDIA,1:KLEV)
    aop_in(KIDIA:KFDIA,1:KLEV)%rg (5) = 1.E6_JPRB * rw_mode (5)%d2(KIDIA:KFDIA,1:KLEV)
    aop_in(KIDIA:KFDIA,1:KLEV)%rg (6) = 1.E6_JPRB * rw_mode (6)%d2(KIDIA:KFDIA,1:KLEV)
    aop_in(KIDIA:KFDIA,1:KLEV)%rg (7) = 1.E6_JPRB * rw_mode (7)%d2(KIDIA:KFDIA,1:KLEV)

    ! dry radius for soluble modes / rest equals the usual radii
    aop_in(KIDIA:KFDIA,1:KLEV)%rgd(1) = 1.E6_JPRB * rwd_mode(1)%d2(KIDIA:KFDIA,1:KLEV)
    aop_in(KIDIA:KFDIA,1:KLEV)%rgd(2) = 1.E6_JPRB * rwd_mode(2)%d2(KIDIA:KFDIA,1:KLEV)
    aop_in(KIDIA:KFDIA,1:KLEV)%rgd(3) = 1.E6_JPRB * rwd_mode(3)%d2(KIDIA:KFDIA,1:KLEV)
    aop_in(KIDIA:KFDIA,1:KLEV)%rgd(4) = 1.E6_JPRB * rwd_mode(4)%d2(KIDIA:KFDIA,1:KLEV)
    aop_in(KIDIA:KFDIA,1:KLEV)%rgd(5) = aop_in(KIDIA:KFDIA,1:KLEV)%rg (5)
    aop_in(KIDIA:KFDIA,1:KLEV)%rgd(6) = aop_in(KIDIA:KFDIA,1:KLEV)%rg (6)
    aop_in(KIDIA:KFDIA,1:KLEV)%rgd(7) = aop_in(KIDIA:KFDIA,1:KLEV)%rg (7)

       
    !VH original:       rm/vol ->units  [N/gridbox] / [m3/gridbox] -> [Number/m3]  
    !VH new   PAERO/ DENS-> units [N/kg] * [kg/m3] = [N/m3]
    aop_in(KIDIA:KFDIA,1:KLEV)%numdens(1) = PAERO(KIDIA:KFDIA,1:KLEV,inus_n) * PRHO(KIDIA:KFDIA,1:KLEV)
    aop_in(KIDIA:KFDIA,1:KLEV)%numdens(2) = PAERO(KIDIA:KFDIA,1:KLEV,iais_n) * PRHO(KIDIA:KFDIA,1:KLEV)
    aop_in(KIDIA:KFDIA,1:KLEV)%numdens(3) = PAERO(KIDIA:KFDIA,1:KLEV,iacs_n) * PRHO(KIDIA:KFDIA,1:KLEV)
    aop_in(KIDIA:KFDIA,1:KLEV)%numdens(4) = PAERO(KIDIA:KFDIA,1:KLEV,icos_n) * PRHO(KIDIA:KFDIA,1:KLEV)
    aop_in(KIDIA:KFDIA,1:KLEV)%numdens(5) = PAERO(KIDIA:KFDIA,1:KLEV,iaii_n) * PRHO(KIDIA:KFDIA,1:KLEV)
    aop_in(KIDIA:KFDIA,1:KLEV)%numdens(6) = PAERO(KIDIA:KFDIA,1:KLEV,iaci_n) * PRHO(KIDIA:KFDIA,1:KLEV)
    aop_in(KIDIA:KFDIA,1:KLEV)%numdens(7) = PAERO(KIDIA:KFDIA,1:KLEV,icoi_n) * PRHO(KIDIA:KFDIA,1:KLEV)
    !TB safeguard for early timesteps when very small negatives may appear
    where(aop_in(KIDIA:KFDIA,1:KLEV)%numdens(1) .lt. 1.E-15_JPRB )aop_in(KIDIA:KFDIA,1:KLEV)%numdens(1) = 1.E-15_JPRB
    where(aop_in(KIDIA:KFDIA,1:KLEV)%numdens(2).lt. 1.E-15_JPRB ) aop_in(KIDIA:KFDIA,1:KLEV)%numdens(2) = 1.E-15_JPRB
    where(aop_in(KIDIA:KFDIA,1:KLEV)%numdens(3) .lt. 1.E-15_JPRB )aop_in(KIDIA:KFDIA,1:KLEV)%numdens(3) = 1.E-15_JPRB
    where(aop_in(KIDIA:KFDIA,1:KLEV)%numdens(4).lt. 1.E-15_JPRB )aop_in(KIDIA:KFDIA,1:KLEV)%numdens(4) = 1.E-15_JPRB 
    where(aop_in(KIDIA:KFDIA,1:KLEV)%numdens(5).lt. 1.E-15_JPRB ) aop_in(KIDIA:KFDIA,1:KLEV)%numdens(5) = 1.E-15_JPRB
    where(aop_in(KIDIA:KFDIA,1:KLEV)%numdens(6) .lt. 1.E-15_JPRB ) aop_in(KIDIA:KFDIA,1:KLEV)%numdens(6) = 1.E-15_JPRB
    where(aop_in(KIDIA:KFDIA,1:KLEV)%numdens(7).lt. 1.E-15_JPRB ) aop_in(KIDIA:KFDIA,1:KLEV)%numdens(7) = 1.E-15_JPRB
    
    ! check valid ranges in particle sizes (might be zero)
    where( aop_in(KIDIA:KFDIA,1:KLEV)%rg (1) .lt. 1.E-15_JPRB ) aop_in(KIDIA:KFDIA,1:KLEV)%rg (1) = 1.E-15_JPRB
    where( aop_in(KIDIA:KFDIA,1:KLEV)%rgd(1) .lt. 1.E-15_JPRB ) aop_in(KIDIA:KFDIA,1:KLEV)%rgd(1) = 1.E-15_JPRB
    where( aop_in(KIDIA:KFDIA,1:KLEV)%rg (2) .lt. 1.E-15_JPRB ) aop_in(KIDIA:KFDIA,1:KLEV)%rg (2) = 1.E-15_JPRB
    where( aop_in(KIDIA:KFDIA,1:KLEV)%rgd(2) .lt. 1.E-15_JPRB ) aop_in(KIDIA:KFDIA,1:KLEV)%rgd(2) = 1.E-15_JPRB
    where( aop_in(KIDIA:KFDIA,1:KLEV)%rg (3) .lt. 1.E-15_JPRB ) aop_in(KIDIA:KFDIA,1:KLEV)%rg (3) = 1.E-15_JPRB
    where( aop_in(KIDIA:KFDIA,1:KLEV)%rgd(3) .lt. 1.E-15_JPRB ) aop_in(KIDIA:KFDIA,1:KLEV)%rgd(3) = 1.E-15_JPRB
    where( aop_in(KIDIA:KFDIA,1:KLEV)%rg (4) .lt. 1.E-15_JPRB ) aop_in(KIDIA:KFDIA,1:KLEV)%rg (4) = 1.E-15_JPRB
    where( aop_in(KIDIA:KFDIA,1:KLEV)%rgd(4) .lt. 1.E-15_JPRB ) aop_in(KIDIA:KFDIA,1:KLEV)%rgd(4) = 1.E-15_JPRB
    where( aop_in(KIDIA:KFDIA,1:KLEV)%rg (5) .lt. 1.E-15_JPRB ) aop_in(KIDIA:KFDIA,1:KLEV)%rg (5) = 1.E-15_JPRB
    where( aop_in(KIDIA:KFDIA,1:KLEV)%rgd(5) .lt. 1.E-15_JPRB ) aop_in(KIDIA:KFDIA,1:KLEV)%rgd(5) = 1.E-15_JPRB
    where( aop_in(KIDIA:KFDIA,1:KLEV)%rg (6) .lt. 1.E-15_JPRB ) aop_in(KIDIA:KFDIA,1:KLEV)%rg (6) = 1.E-15_JPRB
    where( aop_in(KIDIA:KFDIA,1:KLEV)%rgd(6) .lt. 1.E-15_JPRB ) aop_in(KIDIA:KFDIA,1:KLEV)%rgd(6) = 1.E-15_JPRB
    where( aop_in(KIDIA:KFDIA,1:KLEV)%rg (7) .lt. 1.E-15_JPRB ) aop_in(KIDIA:KFDIA,1:KLEV)%rg (7) = 1.E-15_JPRB
    where( aop_in(KIDIA:KFDIA,1:KLEV)%rgd(7) .lt. 1.E-15_JPRB ) aop_in(KIDIA:KFDIA,1:KLEV)%rgd(7) = 1.E-15_JPRB


    ! Initialization to zero is done in tm5m7_optics_calculate_aop 
    Paop_out_ext=0.0_JPRB
    Paop_out_a=0.0_JPRB
    Paop_out_g=0.0_JPRB
    if (present(aop_out_add)) then
       call tm5m7_optics_calculate_aop(KIDIA,KFDIA, KLON,KLEV, nwav,NCONTR, wdep,  ecearth_units, &
       &                           AOP_IN, &       
       &                           Paop_out_ext, Paop_out_a, Paop_out_g, aop_out_add )
    else
       call tm5m7_optics_calculate_aop( KIDIA,KFDIA,KLON,KLEV, nwav, NCONTR, wdep, ecearth_units, &
       &                           AOP_IN, &       
       &                           Paop_out_ext, Paop_out_a, Paop_out_g )
    endif

    ! OK
    Deallocate( aop_in )






IF(LHOOK) CALL DR_HOOK('TM5M7_OPTICS_AOP_GET',1,ZHOOK_HANDLE)
end associate
END SUBROUTINE TM5M7_OPTICS_AOP_GET
