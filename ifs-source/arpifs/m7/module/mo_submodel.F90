  MODULE mo_submodel

  IMPLICIT NONE

  PRIVATE
  
  PUBLIC :: id_ham, lham, lmoz
  LOGICAL :: lham              = .TRUE. ! .true. for aerosol module HAM
  LOGICAL :: lmoz              = .FALSE. ! .true. for gas-phase chemistry module MOZ
  INTEGER :: id_ham 

END MODULE mo_submodel
