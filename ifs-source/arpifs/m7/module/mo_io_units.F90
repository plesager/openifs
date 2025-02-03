MODULE mo_io_units

  IMPLICIT NONE
  
  PUBLIC
  
  INTEGER, PARAMETER :: filename_max = 1024
  
  ! Standard I/O-units
  
  INTEGER, PARAMETER :: nerr  = 0     ! error output
  INTEGER, PARAMETER :: nlog  = 1     ! standard log file unit
  INTEGER, PARAMETER :: nin   = 5     ! standard input
  INTEGER, PARAMETER :: nout  = 6     ! standard output  
  
  INTEGER, PARAMETER, PRIVATE :: none = -1  ! unit given back, when nothing 
                                            ! in the allowed range is available
  
END MODULE mo_io_units











