MODULE mo_filename
  ! -----------------------------------------------------------------
  !
  ! module *mo_filename* - quantities needed for file names etc.
  !
  ! -----------------------------------------------------------------
  !
  ! L. Kornblueh, MPI, April 1998, added NWP forecast mode
  ! I. Kirchner,  MPI, April 2001, revision
  ! A. Rhodin,    MPI, June  2001, adapted to stream interface
  IMPLICIT NONE
  PRIVATE                           
  PUBLIC :: find_next_free_unit
CONTAINS
!------------------------------------------------------------------------------
  FUNCTION find_next_free_unit(istart,istop) RESULT(unit)
    INTEGER :: istart, istop, unit
    LOGICAL :: found, opened
    INTEGER :: i
    found = .FALSE.
    DO i=istart,istop
       INQUIRE(unit=i,opened=opened)
       IF (.NOT.opened) THEN
          unit = i
          found = .TRUE.
          EXIT
       END IF
    END DO
    IF (.NOT. found) THEN
!gf
!       WRITE(message_text,'(a,i2.2,a,i2.2,a)') &
!         'No unit in range <',istart,':',istop,'> free.'
!       CALL finish('find_next_free_unit',message_text)
       WRITE(*,'(a,i2.2,a,i2.2,a)') 'No unit in range <',istart,':',istop,'> free.'
       STOP
!gf 
    END IF
  END FUNCTION find_next_free_unit
!------------------------------------------------------------------------------
END MODULE mo_filename
