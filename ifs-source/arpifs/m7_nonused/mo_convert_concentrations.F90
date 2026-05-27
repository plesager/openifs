!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! mo_convert_concentrations.f90
!!
!! \brief
!! Contains subroutines and functions to convert initial concentrations to pxtm1 and then pxtm1 back
!! to concentrations and mixing ratios. In addition, the gas conversion to mixing ratio
!! is done in a different routine as it is calculated for different time steps.
!!         
!!
!! \author Eemeli Holopainen (FMI)
!!
!! \responsible_coder
!! Eemeli Holopainen, eemeli.holopainen@fmi.fi
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE mo_convert_concentrations

  PUBLIC :: conc2mmr
  PUBLIC :: mmr2conc
  
CONTAINS

  SUBROUTINE conc2mmr(kproma,kbdim,klev, &
       ktrac,pxtm1,paerml,paernl,prhoa)

    USE mo_kind, ONLY : dp
    
    USE mo_species, ONLY: &
       speclist

    USE mo_ham_species, ONLY: &
       id_so4

    USE mo_ham,          ONLY: naerocomp, aerocomp, subm_aerounitconv, subm_aero_idx, &
         sizeclass, mw_so4, nclass, immr2molec
    
    USE mo_physical_constants, ONLY: avo

    USE mo_exception,          ONLY: finish, message, message_text

    IMPLICIT NONE

    INTEGER, PARAMETER :: nmod = 7
    
    !-- input output variables --------
    INTEGER, INTENT(IN) :: &
         kproma, &
         kbdim,  &
         klev,   &    ! number of vertical levels
         ktrac

    REAL(dp), INTENT(IN) :: prhoa(kbdim,klev), &  ! air mass density [kg m-3]
                  paerml(kbdim,klev,naerocomp), &
                  paernl(kbdim,klev,nclass)
    
    REAL(dp), INTENT(INOUT) :: &
         pxtm1(kbdim,klev,ktrac) ! tracer mass/number mixing ratio
       
    
    !-- local variables --------

    INTEGER :: ii, jj, kk
    INTEGER :: jt, jn ,jl, jk, jspec                    ! for indexing
         
    REAL(dp):: zfac,          &
             zfacm,           &
             zfacn,           &
             zfac1, &
             zqunitfac, &
             zmvsu
    
    
    zmvsu = (speclist(id_so4)%moleweight/1000.)/avo/speclist(id_so4)%density
    
    !--- Factor to transform mass SO4 in kg into molecules per kg:
    zfacm  = 6.022e+20_dp/mw_so4
    
    !--- Factor to transform kg into micro gram:
    
    zfac   = 1.e09_dp
    
    !--- Factor to transform N/m**3 into N/cm**3:
    
    zfacn  = 1.0e-06_dp
    
    !--- Calculate mixing ratios from paerml
    DO jn=1, naerocomp
       jt    = aerocomp(jn)%idt       ! get tracer id
       jspec = aerocomp(jn)%spid      ! get species id
       jl    = subm_aero_idx(jspec)     ! get index to subm_aerospec list
       !!mgs=old code!!     IF (aerocomp(jn)%species%m7unitconv == immr2molec) THEN
       IF (jl <= 0) THEN
          WRITE(message_text,*) 'SUBM_AERO_IDX Mapping error !! No index for jspec=',jspec
          CALL finish('ham_subm_interface', message_text)
       END IF
       IF (subm_aerounitconv(jl) == immr2molec) THEN
          zfac1 = zfacm
       ELSE
          zfac1 = zfac
       END IF
       
       pxtm1(1:kproma,:,jt) = paerml(1:kproma,:,jn)/ (zfac1*prhoa(1:kproma,:))  
       
    END DO
    
    !--- Calculate particle numbers from paernl
    
    DO jn=1, nclass
       jt = sizeclass(jn)%idt_no
       !paernl(1:kproma,:,jn) = paernl(1:kproma,:,jn)*1.e-6_dp !eehol: convert 1/m3 to 1/cm3
       pxtm1(1:kproma,:,jt) = paernl(1:kproma,:,jn) / (zfacn*prhoa(1:kproma,:))
    END DO
    
  END SUBROUTINE conc2mmr

  SUBROUTINE mmr2conc(kproma,kbdim,klev, &
       ktrac,pxtm1,paerml,paernl,prhoa)

    USE mo_kind, ONLY : dp

    !USE mo_ham_salsactl,  ONLY: fn2b
    
    USE mo_ham,          ONLY: naerocomp, aerocomp, subm_aerounitconv, subm_aero_idx, &
         sizeclass, mw_so4, nclass, immr2molec
    
    USE mo_exception,          ONLY: finish, message, message_text

    !-- input output variables --------

    INTEGER, INTENT(in) :: &
         kproma, &
         kbdim,  &
         klev,   &    ! number of vertical levels
         ktrac
    
    REAL(dp), INTENT(INOUT) :: &
         paerml(kbdim,klev,naerocomp), &
         paernl(kbdim,klev,nclass)

    REAL(dp), INTENT(IN) ::  &
         pxtm1(kbdim,klev,ktrac), & ! tracer mass/number mixing ratio
         prhoa(kbdim,klev)   ! air mass density [kg m-3]
    
    !-- local variables --------

    REAL(dp):: zfac,            &
             zfacm,           &
             zfacn,           &
             zfac1, &
             zqunitfac, &
             zmvsu

    INTEGER :: jt, jn ,jl, jk, jspec                    ! for indexing


    !--- Factor to transform mass SO4 in kg into molecules per kg:
    zfacm  = 6.022e+20_dp/mw_so4
    
    !--- Factor to transform kg into micro gram:
    
    zfac   = 1.e09_dp
    
    !--- Factor to transform N/m**3 into N/cm**3:
    
    zfacn  = 1.0e-06_dp
    
    !--- Calculate mixing ratios from pxtm1
    DO jn=1, naerocomp
       jt    = aerocomp(jn)%idt       ! get tracer id
       jspec = aerocomp(jn)%spid      ! get species id
       jl    = subm_aero_idx(jspec)     ! get index to subm_aerospec list
       !!mgs=old code!!     IF (aerocomp(jn)%species%m7unitconv == immr2molec) THEN
       IF (jl <= 0) THEN
          WRITE(message_text,*) 'SUBM_AERO_IDX Mapping error !! No index for jspec=',jspec
          CALL finish('ham_subm_interface', message_text)
       END IF
       IF (subm_aerounitconv(jl) == immr2molec) THEN
          zfac1 = zfacm
       ELSE
          zfac1 = zfac
       END IF
       
       paerml(1:kproma,:,jn) = pxtm1(1:kproma,:,jt)*(zfac1*prhoa(1:kproma,:))
       
    END DO
    
    !--- Calculate particle numbers from pxtm1
    
    DO jn=1, nclass
       jt = sizeclass(jn)%idt_no
       paernl(1:kproma,:,jn) = zfacn*prhoa(1:kproma,:)*pxtm1(1:kproma,:,jt)!*1.e6_dp !eehol: convert 1/cm3 to 1/m3
    END DO

  END SUBROUTINE mmr2conc
END MODULE mo_convert_concentrations
