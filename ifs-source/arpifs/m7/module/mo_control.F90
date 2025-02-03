MODULE mo_control

  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  SAVE
  
  INTEGER :: nn             !   max meridional wave number for m=0.
  INTEGER :: ngl            !   number of gaussian latitudes. 
  INTEGER :: nlon           !   max number of points on each latitude line. 
  INTEGER :: nlev           !   number of vertical levels. 
  INTEGER :: nmp1           !   max zonal wave number + 1. 
  INTEGER :: nlevp1         !   *nlev+1.
  INTEGER :: nsp            !   number of spectral coefficients. 
  INTEGER :: nhgl           !   (number of gaussian latitudes)/2. 
  INTEGER :: nvclev         !   number of levels with vertical coefficients. 
  LOGICAL :: lrce      = .FALSE. 
  LOGICAL :: ltimer            = .FALSE. ! to use timer
  LOGICAL :: ldebugio          = .FALSE. ! to debug IO 

END MODULE mo_control
