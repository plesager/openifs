MODULE mo_radiation_parameters

  USE mo_kind,            ONLY: wp
  USE mo_math_constants,  ONLY: pi

IMPLICIT NONE

  PRIVATE


  PUBLIC :: iaero, decl_sun_cur

  INTEGER :: iaero = 2  !< aerosol model
  REAL(wp) :: decl_sun_cur                  !< solar declination at current time step

END MODULE mo_radiation_parameters
