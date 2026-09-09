!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! mo_submodel.f90
!!
!! \belongs_to
!!  HAMMOZ
!!
!!  SPDX-License-Identifier: BSD-3-Clause
!! Copyright (c) 2021 hammoz
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE mo_submodel

  IMPLICIT NONE

  PRIVATE
  
  PUBLIC :: id_ham, lham, lmoz
  LOGICAL :: lham              = .TRUE. ! .true. for aerosol module HAM
  LOGICAL :: lmoz              = .FALSE. ! .true. for gas-phase chemistry module MOZ
  INTEGER :: id_ham 

END MODULE mo_submodel
