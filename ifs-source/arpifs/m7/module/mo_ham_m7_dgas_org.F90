!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! mo_ham_m7_dgas_org.F90
!!
!! \brief
!! mo_ham_m7_dgas_org provides routines for calculating aerosol organic condensation.
!!
!! \authors
!! J. Wilson, E. Vignati, JRC/EI (original source)                05/2000
!! R. Makkonen, UHEL                                              2015
!! Hermanni Halonen (adapted from ECE3)
!!
!! \responsible_coder
!! Hermanni Halonen, hermanni.halonen@uef.fi
!!
!! \limitations
!! None
!!
!! \details
!! This module provides routines for calculating organic nucleation for SVOC and ELVOC
!!
!! \belongs_to
!!  HAMMOZ
!!
!!  SPDX-License-Identifier: BSD-3-Clause
!! Copyright (c) 2021 hammoz
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE mo_ham_m7_dgas_org
  USE mo_kind, ONLY: dp
  USE mo_time_control,  ONLY: time_step_len
  USE mo_math_constants, ONLY: pi
  USE mo_physical_constants, ONLY: avo
  USE mo_ham_m7ctl,     ONLY: caccso4, iocks, iocas, wh2so4, rerg, &
      ioccs, iocki, isoans, isoaks, isoaas, isoacs, isoaki
  USE mo_ham_m7_cond_sink, ONLY: condensation_sink
  USE mo_ham, ONLY: naerocomp, nclass
  IMPLICIT NONE
  !
  !--- Parameter list:
  !
  ! pm6rp     = mean mode actual radius (wet radius for soluble modes 
  !             and dry radius for insoluble modes) [cm]
  ! 
  !--- Local Variables:
  !
  ! zde       = molecular diffusion []
  ! zvelb     = velocity []
  ! zfcond    = total mass of condensing sulfate for one timestep []

CONTAINS

  SUBROUTINE m7_dgas_org(kproma, kbdim, klev, pelvoc, psvoc, paerml, paernl, &
                         ptp1, papp1, pm6rp)
    ! Calculates the transfer of mass due to organic condensation
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: kproma, kbdim, klev
    REAL(dp), DIMENSION(kbdim,klev), INTENT(IN) :: ptp1, papp1, pelvoc, psvoc
    REAL(dp), DIMENSION(kbdim,klev,nclass), INTENT(INOUT) :: paernl, pm6rp
    REAL(dp), DIMENSION(kbdim,klev,naerocomp), INTENT(INOUT) :: paerml

    INTEGER :: jl, jk, jmod
    REAL(dp) :: zfcond, zftot, zpbyone, zde2, zvelb, zxibc, zm6rp, zf1, zqtmst
    REAL(dp) :: lambda, cs(7), cs_sum, modeweight(7), modeweight_sum
    
    ! Initialize variables
    zqtmst = 1.0_dp / time_step_len

    !--- 1) Calculate condensation rate for cm diameter sulphate aerosols: ---
    !
    DO jmod = 1, 7
      DO jk = 1, klev
        DO jl = 1, kproma
          IF (pm6rp(jl,jk,jmod) > 0.0_dp) THEN
            !--- Diffusion coefficient (Reference???):
	        zpbyone = 1000.0_dp / (papp1(jl,jk) / 100.0_dp)
	        
            zde2 = 0.073_dp * zpbyone * (ptp1(jl,jk) / 298.15_dp)**1.5_dp
            
	        !--- Mean molecule velocity (Moore, 1962 (S+P equ. 8.2)):
            zvelb = SQRT(8.0_dp * rerg * ptp1(jl,jk) / pi / wh2so4)
            
            !--- ???Fuchs???
	        zxibc = 8.0_dp * zde2 / pi / zvelb
	        
            ! Use count median radius:
	        zm6rp = pm6rp(jl,jk,jmod)
	        
	        !--- Distance from particle up to which the kinetic regime applies:
            zf1 = ((zm6rp + zxibc)**3.0_dp - (zm6rp**2.0_dp + zxibc**2.0_dp)**1.5_dp) / &
                  (3.0_dp * zm6rp * zxibc) - zm6rp

          END IF
        END DO
      END DO
    END DO

    !--- 2) Calculation of the new organic aerosol masses and of the ---------
    !       mass of organic condensing on the respective modes:
    DO jk = 1, klev
      DO jl = 1, kproma
        lambda = 6.6E-8_dp
        !
        !--- Calculate condensation sink and OA-mass of each mode
        !
        CALL condensation_sink(7, paernl(jl,jk,:), 0.02_dp * pm6rp(jl,jk,:), lambda, cs)
        modeweight = 0.0_dp
        modeweight(2) = MIN(MAX(paerml(jl,jk,iocks) + paerml(jl,jk,isoaks), 0.0_dp), 1.0E5_dp)
        modeweight(3) = MIN(MAX(paerml(jl,jk,iocas) + paerml(jl,jk,isoaas), 0.0_dp), 1.0E5_dp)
        modeweight(4) = MIN(MAX(paerml(jl,jk,ioccs) + paerml(jl,jk,isoacs), 0.0_dp), 1.0E5_dp)
        modeweight(5) = MIN(MAX(paerml(jl,jk,iocki) + paerml(jl,jk,isoaki), 0.0_dp), 1.0E5_dp)
        cs_sum = SUM(cs(1:5))
        modeweight_sum = SUM(modeweight)
        
        ! Condense S/LVOC according to mode OA mass
        IF (modeweight_sum > 1.0E-13_dp .AND. psvoc(jl,jk) > 1.0E-16_dp .AND. psvoc(jl,jk) < 1.0E10_dp) THEN
          paerml(jl,jk,isoaks) = paerml(jl,jk,isoaks) + psvoc(jl,jk) * modeweight(2) / modeweight_sum
          paerml(jl,jk,isoaas) = paerml(jl,jk,isoaas) + psvoc(jl,jk) * modeweight(3) / modeweight_sum
          paerml(jl,jk,isoacs) = paerml(jl,jk,isoacs) + psvoc(jl,jk) * modeweight(4) / modeweight_sum
          paerml(jl,jk,isoaki) = paerml(jl,jk,isoaki) + psvoc(jl,jk) * modeweight(5) / modeweight_sum
        
	    !--- If only small amount of existing OA, condensing to accumulation mode
        ELSE IF (psvoc(jl,jk) > 1.0E-16_dp .AND. psvoc(jl,jk) < 1.0E10_dp) THEN
          paerml(jl,jk,isoaas) = paerml(jl,jk,isoaas) + psvoc(jl,jk)
        END IF
        
        ! Condense ELVOC according to condensation sink
        IF (cs_sum > 1.0E-15_dp .AND. cs_sum < 1.0E5_dp .AND. pelvoc(jl,jk) > 1.0E-16_dp .AND. pelvoc(jl,jk) < 1.0E10_dp) THEN
          paerml(jl,jk,isoans) = paerml(jl,jk,isoans) + pelvoc(jl,jk) * cs(1) / cs_sum
          paerml(jl,jk,isoaks) = paerml(jl,jk,isoaks) + pelvoc(jl,jk) * cs(2) / cs_sum
          paerml(jl,jk,isoaas) = paerml(jl,jk,isoaas) + pelvoc(jl,jk) * cs(3) / cs_sum
          paerml(jl,jk,isoacs) = paerml(jl,jk,isoacs) + pelvoc(jl,jk) * cs(4) / cs_sum
          paerml(jl,jk,isoaki) = paerml(jl,jk,isoaki) + pelvoc(jl,jk) * cs(5) / cs_sum
       !--- If only small CS, condensing to accumulation mode 
       ELSE IF (pelvoc(jl,jk) > 1.0E-16_dp .AND. pelvoc(jl,jk) < 1.0E10_dp) THEN
          paerml(jl,jk,isoaas) = paerml(jl,jk,isoaas) + pelvoc(jl,jk)
        END IF
      END DO
    END DO

  END SUBROUTINE m7_dgas_org

END MODULE mo_ham_m7_dgas_org