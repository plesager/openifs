!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! mo_radiation_parmeters.f90
!!
!! \belongs_to
!!  HAMMOZ
!!
!!  SPDX-License-Identifier: BSD-3-Clause
!! Copyright (c) 2021 hammoz
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE mo_radiation_parameters

  USE mo_kind,            ONLY: wp
  USE mo_math_constants,  ONLY: pi

IMPLICIT NONE

  PRIVATE


  PUBLIC :: iaero, decl_sun_cur

  INTEGER :: iaero = 1  !< aerosol model
  REAL(wp) :: decl_sun_cur                  !< solar declination at current time step

END MODULE mo_radiation_parameters
