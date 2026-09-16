!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! mo_kind.f90
!!
!! \belongs_to
!!  HAMMOZ
!!
!!  SPDX-License-Identifier: BSD-3-Clause
!! Copyright (c) 2021 hammoz
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE mo_kind
  ! --> thk: bug fix:
  ! changed the order of the two lines
  ! removed commented-out code
  USE parkind1, ONLY: dp => JPRB !TeMi 
  USE parkind1, ONLY: wp => JPRB !eehol
  USE parkind1, ONLY: JPRD
  IMPLICIT NONE
  ! <-- thk

  ! Threshold for branch cutoffs. Historically the double precision
  ! epsilon was used throughout the M7 code. That value is now also
  ! used to define the **physical cutoffs** in single precision.
  REAL(KIND=dp), PARAMETER :: THRESHOLD = REAL(EPSILON(1.0_JPRD), KIND=dp)
  
END MODULE mo_kind
