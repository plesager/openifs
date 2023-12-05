Module mo_time_control
  USE mo_kind, only : dp
  USE YOMRIP   , ONLY : YRRIP
  IMPLICIT NONE
  ! <-- thk: bug fix
  !REAL(dp), PUBLIC   :: time_step_len    = 720.0_dp
  REAL(dp), PUBLIC :: time_step_len
  !REAL(dp), PUBLIC :: delta_time !eehol: not needed
  ! this does not work:
  !ASSOCIATE(time_step_len=>YRPHY2%TSPHY) !TeMi

  CONTAINS
    SUBROUTINE init_mo_time_control
      ! copy the time step length from IFS control structure
      time_step_len = YRRIP%TSTEP
    END SUBROUTINE init_mo_time_control
  ! --> thk
END MODULE mo_time_control

