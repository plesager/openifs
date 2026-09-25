MODULE mo_ham_m7_cond_sink
  USE mo_kind, ONLY: dp
  IMPLICIT NONE

CONTAINS

  SUBROUTINE condensation_sink(NMODE, N, MEAN_IN, LAMBDA_IN, OUTPUT)
    ! This is a very approximative tool to calculate the division factors for organic condensation
    ! Assumes sigma = 1.5 (will underestimate condensation to coarse mode by a factor of ~1-5)
    ! AA (2006)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NMODE       ! Number of modes considered
    REAL(dp), INTENT(IN) :: LAMBDA_IN  ! Input lambda value
    REAL(dp), DIMENSION(NMODE), INTENT(IN) :: N       ! Number concentrations [arbitrary]
    REAL(dp), DIMENSION(NMODE), INTENT(IN) :: MEAN_IN ! Current geometric median diameters of the modes [m]
    REAL(dp), DIMENSION(NMODE), INTENT(OUT) :: OUTPUT ! Output coefficients

    REAL(dp) :: LAMBDA, LLAM
    REAL(dp), DIMENSION(NMODE) :: MEAN, CS, LMEAN
    REAL(dp), DIMENSION(3) :: COEF
    ! Fit of parameters in respect to log10(lambda)
    REAL(dp), DIMENSION(3), PARAMETER :: C1 = (/ 0.0292_dp, 0.4376_dp, 1.5011_dp /) ! 2nd order
    REAL(dp), DIMENSION(3), PARAMETER :: C2 = (/ 0.4289_dp, 6.6937_dp, 25.4629_dp /) ! 1st order
    REAL(dp), DIMENSION(3), PARAMETER :: C3 = (/ 1.4346_dp, 23.0498_dp, 87.5852_dp /) ! 0th order
    
    REAL(dp), DIMENSION(NMODE) :: NCS
    REAL(dp) :: summa
    INTEGER :: indeksi

    ! Assign input values to local variables
    LAMBDA = LAMBDA_IN
    MEAN = MEAN_IN

    ! Check parameter validity
    IF (LAMBDA > 7.1220e-07_dp) LAMBDA = 7.1220e-07_dp
    IF (LAMBDA < 6.5100e-08_dp) LAMBDA = 6.5100e-08_dp
    WHERE (MEAN < 1e-10_dp) MEAN = 1e-10_dp
    WHERE (MEAN > 10e-6_dp) MEAN = 10e-6_dp

    ! Calculate logarithms
    LMEAN = LOG10(MEAN)
    LLAM = LOG10(LAMBDA)

    ! Calculate coefficients
    COEF(1) = LLAM * (LLAM * C1(1) + C1(2)) + C1(3)
    COEF(2) = LLAM * (LLAM * C2(1) + C2(2)) + C2(3)
    COEF(3) = LLAM * (LLAM * C3(1) + C3(2)) + C3(3)

    ! Calculate condensation sink
    CS = 10.0_dp ** (LMEAN * (LMEAN * COEF(1) + COEF(2)) + COEF(3))
    
    IF (SUM(N * CS) > 1e-20_dp) THEN
      OUTPUT = (N * CS) / SUM(N * CS)
    ELSE
      OUTPUT = 1.0_dp / REAL(NMODE, dp)
    END IF

  END SUBROUTINE condensation_sink

END MODULE mo_ham_m7_cond_sink