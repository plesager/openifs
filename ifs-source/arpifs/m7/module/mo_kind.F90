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

  ! Threshold for branch cutoffs **and** to check for safe divisions
  !
  ! Historically (when the model was run only in double precision) the
  ! double precision epsilon was used throughout the M7 code. That DP
  ! value is now also used to define the physical cutoffs in single
  ! precision runs. It is also better than (too large) SP epsilon to
  ! check on safe divisions. Consequently we are still using the same
  ! threshold for both goals and regardless of the precision.
  ! 
  REAL(KIND=dp), PARAMETER :: THRESHOLD = REAL(EPSILON(1.0_JPRD), KIND=dp)
  
END MODULE mo_kind
