SUBROUTINE TM5M7_GET_REFR_IDX(wdep, SO4, BC, OC, SOA, SS, DU, water, mode, m_eff)

!*** * TM5M7_GET_REFR_IDX* 
!
!
!-----------------------------------------------------------------------------
!                    TM5                                                     !
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE:	TM5M7_GET_REFR_IDX
!
! !DESCRIPTION: Compute refractive index of internally mixed aerosols by use
!               of effective medium theory for the size-dependent aerosol
!               mixtures assumed in M7.
!
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
!
! !REVISION HISTORY: 
!
!   12 Aug 2008 - Michael Kahnert, SMHI
!    6 Feb 2011 - Achim Strunk -
!
!      Sep 2021 - V. Huijnen: first introduction into OpenIFS
!
!
!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMLUN    ,ONLY : NULERR

USE YOMCST, ONLY : RPI
USE TM5M7_DATA, ONLY : ss_density, dust_density, carbon_density, &
  & pom_density, so4_density, soa_density 
USE TM5M7_OPTICS_DATA, ONLY: WAVELENDEP

IMPLICIT NONE

!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

type(wavelendep), INTENT(IN)    :: wdep   ! wavelength properties (wavelength, re/img part of refractive index)
REAL(KIND=JPRB),             INTENT(IN)    :: SO4, BC, OC,SOA    ! mass mixing ratios or concentrations of sulphate, black carbon, organic carbon
REAL(KIND=JPRB),             INTENT(IN)    :: SS,  DU, water ! sea salt, dust, and water
INTEGER(KIND=JPIM),          INTENT(IN)    :: mode   ! mode number (M7)
!
! !OUTPUT PARAMETERS:
!
COMPLEX,                     INTENT(OUT)  :: m_eff  ! effective refractive index of mixture
!INTEGER(KIND=JPIM),          INTENT(OUT)  :: status



!*       0.2   LOCAL VARIABLES
!              ---------------

! refractive indices
COMPLEX            :: m_SO4, m_BC, m_OC, m_SOA, m_SS, m_DU, m_water

! volume fractions
REAL(KIND=JPRB)   :: v_SO4, v_BC, v_OC, v_SOA, v_SS, v_DU, v_water, water_iv

REAL(KIND=JPRB)    :: vtot, v2
COMPLEX            :: m00, m0, m1, m2

!VH REAL(KIND=JPRB)    :: rpls, ipls
REAL(KIND=JPRB),PARAMETER :: ROL = 1000. ! kg/m^3


REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TM5M7_GET_REFR_IDX',0,ZHOOK_HANDLE)
  
    ! Get the refractive indices from the lookup-tables and put them into COMPLEX numbers.
    m_so4    = cmplx( wdep%n(1),wdep%k(1) ) ! H2-SO4 + NH4NO3
    m_bc     = cmplx( wdep%n(2),wdep%k(2) ) ! BC
    m_oc     = cmplx( wdep%n(3),wdep%k(3) ) ! POM
    m_soa    = cmplx( wdep%n(4),wdep%k(4) ) ! SOA
    m_ss     = cmplx( wdep%n(5),wdep%k(5) ) ! SS
    m_du     = cmplx( wdep%n(6),wdep%k(6) ) ! DU
    m_water  = cmplx( wdep%n(7),wdep%k(7) ) ! Water

    !status = 0
    !     no mixing for mode=6,7:
    if(mode.ge.6)then
       m_eff=m_DU
       IF(LHOOK) CALL DR_HOOK('TM5M7_GET_REFR_IDX',1,ZHOOK_HANDLE)
       RETURN
    endif

    !     compute volume fractions:
    v_SO4=0._JPRB
    v_BC=0._JPRB
    v_OC=0._JPRB
    v_SOA=0._JPRB
    v_SS=0._JPRB
    v_DU=0._JPRB
    v_water=0._JPRB
    vtot=0._JPRB

    ! Added sanity check (15-7-2010 - P. Le Sager) : Avoid negative water
    ! mixing ratio!
    ! The bruggeman logically assumes that v_water is between 0 and 1, but
    ! this is never checked in the call chain : 
    ! ECEarth_Optics_Step -> calculate_aop -> get_refr_idx [here] ->
    ! Bruggeman
    ! We do it here, with a warning since it reflects a problem upstream:
    if(water.lt.0.0_JPRB)then
       !write (gol,*)" WARNING - [Get_refr_idx] : negative relative humidity..." 
       !write (gol,*)" WARNING - [Get_refr_idx] : .....set to 0" 
       water_iv=0.0_JPRB
    else
       water_iv=water
    endif

    if(mode.le.4)then
       v_SO4=SO4/so4_density
       v_water=water_iv/rol
       vtot=vtot+v_SO4+v_water
    endif
    if(mode.le.5)then
       !v_OC=OC/pom_density
       v_SOA=SOA/soa_density
       vtot=vtot+v_SOA
    end if
    if(mode.ge.2.and.mode.le.5)then
       v_BC=BC/carbon_density
       v_OC=OC/pom_density
       vtot=vtot+v_BC+v_OC
    endif
    if(mode.ge.3.and.mode.le.4)then
       v_SS=SS/ss_density
       vtot=vtot+v_SS
    endif
    if(mode.ge.3.and.mode.ne.5)then
       v_DU=DU/dust_density
       vtot=vtot+v_DU
    endif
    ! If vtot is zero, we will get 0.0/0.0's causing NaNs. In that case, the 
    ! refractive index does not matter and will be set to (1.0,1.0e-9). The 
    ! reason not to take (1.0,0.0) is that someone with humour might take 
    ! the logarithm of the imaginary part. Dust particles get their usual 
    ! refractive index, because they already returned m_DU. But that does 
    ! not matter, because there are zero aerosols in this case.
    if (vtot .le. 1E-20_JPRB) then
       m_eff = Cmplx(1.0,1.0e-9)
    else
       v_SO4=MAX(1E-20,v_SO4/vtot)
       v_BC=MAX(1E-20,v_BC/vtot)
       v_OC=MAX(1E-20,v_OC/vtot)
       v_SOA=MAX(1E-20,v_SOA/vtot)
       v_SS=MAX(1E-20,v_SS/vtot)
       v_DU=MAX(1E-20,v_DU/vtot)
       v_water=MAX(1E-20,v_water/vtot)

       !-----------------------------------------------------------------------
       !     effective medium computations for each mode
       !-----------------------------------------------------------------------
       if(mode.eq.1)then
          !        Bruggeman mixing rule for SO4 OC and water:
          m1=m_SO4
          m2=m_SOA
          vtot=v_SO4+v_SOA
          v2=v_SOA/vtot
          call Bruggeman(m1,m2,v2,m0)
          m1=m0
          m2=m_water
          v2=v_water
          call Bruggeman(m1,m2,v2,m0)
       elseif(mode.eq.2)then
          !        iterative Bruggeman mixing rule for SO4, OC, and water:
          m1=m_SO4
          m2=m_OC
          vtot=v_SO4+v_OC
          v2=v_OC/vtot
          call Bruggeman(m1,m2,v2,m00)
          m1=m00
          m2=m_SOA
          vtot=vtot+v_SOA
          v2=v_SOA/vtot
          call Bruggeman(m1,m2,v2,m00)
          m1=m00
          m2=m_water
          vtot=vtot+v_water
          v2=v_water/vtot
          call Bruggeman(m1,m2,v2,m00)
          !        Maxwell-Garnett mixing rule for BC inclusions:
          m1=m00
          m2=m_BC
          v2=v_BC
          call Maxwell_Garnett(m1,m2,v2,m0)
       elseif(mode.eq.3.or.mode.eq.4)then
          !        iterative Bruggeman mixing rule for SO4, OC, SS, and water:
          m1=m_SO4
          m2=m_OC
          vtot=v_SO4+v_OC
          if ( vtot < TINY( vtot ) ) then
             v2=0.0_JPRB
          else
             v2=v_OC/vtot
          end if
          call Bruggeman(m1,m2,v2,m00)
          m1=m00
          m2=m_SOA
          vtot=vtot+v_SOA
          if ( vtot < TINY( vtot ) ) then
             v2=0.0_JPRB
          else
             v2=v_SOA/vtot
          end if
          call Bruggeman(m1,m2,v2,m00)
          m1=m00
          m2=m_SS
          vtot=vtot+v_SS
          if ( vtot < TINY( vtot ) ) then
             v2=0.0_JPRB
          else
             v2=v_SS/vtot
          end if
          call Bruggeman(m1,m2,v2,m00)
          m1=m00
          m2=m_water
          vtot=vtot+v_water
          v2=v_water/vtot
          call Bruggeman(m1,m2,v2,m00)
          !        iterative Maxwell-Garnett mixing rule for BC and dust 
          !        inclusions:
          m1=m00
          m2=m_BC
          vtot=vtot+v_BC
          if ( vtot < TINY( vtot ) ) then
             v2=0.0_JPRB
          else
             v2=v_BC/vtot
          end if
          call Maxwell_Garnett(m1,m2,v2,m00)
          m1=m00
          m2=m_DU
          v2=v_DU
          call Maxwell_Garnett(m1,m2,v2,m0)
       elseif(mode.eq.5)then

          m1=m_SOA
          m2=m_OC
          vtot=v_SOA+v_OC
          v2=v_OC/vtot
          call Bruggeman(m1,m2,v2,m00)
          !        Maxwell-Garnett mixing rule for BC inclusions:
          m1=m00
          m2=m_BC
          v2=v_BC
          call Maxwell_Garnett(m1,m2,v2,m0)          

       endif
       m_eff = m0
    End If

    ! Debug : trap for a NAN (13-7-2010 - P. Le Sager)
    ! rpls=real(m_eff)
    ! ipls=imag(m_eff)
    ! IF ((rpls.NE.rpls).or.(ipls.NE.ipls)) then
    !   status = 1
    !   write (NULERR,'(" GET_REFR_IDX-NAN:  ", 3(E16.4,2x),i4,2x,7(E16.4,2x))') rpls, ipls, vtot, mode,&
    !        &    SO4,BC,OC,SOA,SS,DU,water
    ! endif

IF(LHOOK) CALL DR_HOOK('TM5M7_GET_REFR_IDX',1,ZHOOK_HANDLE)
END SUBROUTINE TM5M7_GET_REFR_IDX




!--------------------------------------------------------------------------
!                    TM5                                                  !
!--------------------------------------------------------------------------
!BOP
!
! !IROUTINE:    BRUGGEMAN
!
! !DESCRIPTION: Compute effective refractive index of a mixture of 2 components
!               by use of the Bruggeman mixing rule
!\\
!\\
! !INTERFACE:
!
PURE SUBROUTINE BRUGGEMAN(m1,m2,v2,m0)
  !
  ! !INPUT PARAMETERS:
  !
  USE PARKIND1  ,ONLY : JPIM     ,JPRB
  COMPLEX, INTENT(IN)             :: m1,m2
  REAL(KIND=JPRB),    INTENT(IN)  :: v2
  !
  ! !OUTPUT PARAMETERS:
  !
  COMPLEX, INTENT(OUT):: m0
  !
  ! !REVISION HISTORY: 
  !   12 Aug 2008 - Michael Kahnert, SMHI
  !    6 Feb 2011 - Achim Strunk -
  !
  ! !REMARKS:
  !
  !EOP
  !------------------------------------------------------------------------
  !BOC
  
  !local:
  COMPLEX         ::  m1s,m2s,mt
  REAL(KIND=JPRB) :: fac1,fac2
  
  !Begin
  
  if(v2.eq.1.0_JPRB)then
     m0=m2
     ! IF(LHOOK) CALL DR_HOOK('routine_name',1,ZHOOK_HANDLE)
     RETURN
  elseif(v2.eq.0.0_JPRB)then
     m0=m1
     ! IF(LHOOK) CALL DR_HOOK('routine_name',1,ZHOOK_HANDLE)
     RETURN
  endif
  
  fac1=2._JPRB-3._JPRB*v2
  fac2=3._JPRB*v2-1._JPRB
  m1s=m1**2_JPIM
  m2s=m2**2_JPIM
  mt=m1s*fac1+m2s*fac2
  m0=1._JPRB/16._JPRB*mt**2_JPIM+0.5_JPRB*m1s*m2s
  m0=csqrt(m0)
  m0=m0+0.25_JPRB*mt
  m0=csqrt(m0)

 END SUBROUTINE BRUGGEMAN
!EOC





!--------------------------------------------------------------------------
!                    TM5                                                  !
!--------------------------------------------------------------------------
!BOP
!
! !IROUTINE:    MAXWELL_GARNETT
!
! !DESCRIPTION: Compute effective refractive index for a medium consisting of
!               a matrix with refractive index m1 and inclusions with refractive
!               index m2 and volume fraction v2 by use of the Maxwell-Garnett 
!               mixing rule.
!\\
!\\
! !INTERFACE:
!
PURE SUBROUTINE MAXWELL_GARNETT( m1, m2, v2, m0)
  !
  ! !INPUT PARAMETERS:
  !
  USE PARKIND1  ,ONLY : JPIM     ,JPRB
  COMPLEX, INTENT(IN)             :: m1, m2
  REAL(KIND=JPRB),    INTENT(IN)  :: v2
  !
  ! !OUTPUT PARAMETERS:
  !
  COMPLEX, INTENT(OUT):: m0
  !
  ! !REVISION HISTORY: 
  !   12 Aug 2008 - Michael Kahnert, SMHI
  !    6 Feb 2011 - Achim Strunk -
  !
  ! !REMARKS:
  !
  !EOP
  !------------------------------------------------------------------------
  !BOC

  !local:
  COMPLEX         :: m1s,m2s
  REAL(KIND=JPRB) :: fac1,fac2
  
  ! Begin:

  fac1=3.0_JPRB-2.0_JPRB*v2
  fac2=3.0_JPRB-v2
  m1s=m1**2_JPIM
  m2s=m2**2_JPIM

  m0=m2s*(fac1*m1s+2.0*v2*m2s)/(v2*m1s+fac2*m2s)
  m0=csqrt(m0)

END SUBROUTINE MAXWELL_GARNETT
!EOC


