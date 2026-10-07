!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!! \filename 
!! mo_exception.F90
!!
!! \belongs_to
!!  HAMMOZ
!!
!!  SPDX-License-Identifier: BSD-3-Clause
!! Copyright (c) 2021 hammoz


MODULE mo_exception

  USE mo_io_units, ONLY: nerr, nlog
  
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: message_text
  PUBLIC :: message, finish
  PUBLIC :: em_none, em_info, em_warn, em_error, em_param, em_debug

  INTEGER, PARAMETER :: em_none  = 0   ! normal message
  INTEGER, PARAMETER :: em_info  = 1   ! informational message
  INTEGER, PARAMETER :: em_warn  = 2   ! warning message: number of warnings counted
  INTEGER, PARAMETER :: em_error = 3   ! error message: number of errors counted
  INTEGER, PARAMETER :: em_param = 4   ! report parameter value
  INTEGER, PARAMETER :: em_debug = 5   ! debugging message

  CHARACTER(len=256) :: message_text = ''         !++mgs

  LOGICAL :: l_log   = .FALSE.

  INTEGER :: number_of_warnings  = 0
  INTEGER :: number_of_errors    = 0

#include "abor1.intfb.h"

CONTAINS
  
  SUBROUTINE finish (name, text, exit_no)

    CHARACTER(len=*), INTENT(in)           :: name
    CHARACTER(len=*), INTENT(in), OPTIONAL :: text
    INTEGER,          INTENT(in), OPTIONAL :: exit_no

    INTEGER           :: iexit

    CALL ABOR1(text)
       
  END SUBROUTINE finish

  SUBROUTINE message (name, text, out, level, all_print, adjust_right)

    CHARACTER (len=*), INTENT(in) :: name
    CHARACTER (len=*), INTENT(in) :: text
    INTEGER,           INTENT(in), OPTIONAL :: out
    INTEGER,           INTENT(in), OPTIONAL :: level
    LOGICAL,           INTENT(in), OPTIONAL :: all_print
    LOGICAL,           INTENT(in), OPTIONAL :: adjust_right

    INTEGER :: iout
    INTEGER :: ilevel
    LOGICAL :: lprint
    LOGICAL :: ladjustr     !++mgs renamed from ladjust to ladjustr

    CHARACTER(len=32) :: prefix

    CHARACTER(len=LEN(message_text)) :: write_text

    WRITE(nlog,*) name, text
    
  END SUBROUTINE message

END MODULE mo_exception
