!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! mo_param_switches.f90
!!
!! \belongs_to
!!  HAMMOZ
!!
!!  SPDX-License-Identifier: BSD-3-Clause
!! Copyright (c) 2021 hammoz
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE mo_param_switches

  IMPLICIT NONE

  LOGICAL :: lcdnc_progn  !   true for prognostic cloud droplet activation !HK
  INTEGER :: ncd_activ    !   type of cloud droplet activation scheme (see physctl.inc) !HK
  INTEGER :: nactivpdf    !   sub-grid scale PDF of updraft velocities (see physctl.inc) !ZK
  INTEGER ::  cdnc_min_fixed ! minimum number of cloud droplet number concentration cm-3
END MODULE mo_param_switches
