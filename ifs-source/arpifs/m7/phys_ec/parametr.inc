C=======================================================================
C
C *** INCLUDE FILE 'PARAMETR.INC'
C *** THIS FILE CONTAINS THE DECLARATIONS OF THE GLOBAL CONSTANTS
C     AND VARIABLES. 
C
C *** WRITTEN BY ATHANASIOS NENES 
C *** MODIFIED BY ATHANASIOS NENES AND PRASHANT KUMAR
C
C=======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      ! Three points is enough for PDF integration
      ! using Gauss-Legendre quadrature
      PARAMETER (NpGauss=3)
      !PARAMETER (NpGauss=10)
      ! Maximum number of lognormal modes,
      ! here set to three for M7.
      PARAMETER (NSMX = 3)
      REAL      XGS, WGS
      LOGICAL   CRIT2, CCNSPST
      LOGICAL   FIRST_GAULEG

      COMMON /INPUTS/ TEMP, PRES
C
     &       /CCNSPC/ DPG(NSMX),  SIG(NSMX), Dpc(NSMX),
     &                MODE(NSMX), NMD 
C
     &       /ACTVPR/ AKOH, SSPLT, ALFA, BET1, CRIT2, CCNSPST
C
     &       /THERMO/ AMW, AMA, GRAV, DENW, CPAIR, DHV, RGAS, AKA,
     &                PSAT, DAIR, SURT, Dw
C
     &       /SLNPAR/ EPS, MAXIT, NITER
C
     &       /GAUSSL/ XGS(NpGauss), WGS(NpGauss), FIRST_GAULEG
C
     &       /OTHER/  PI, ZERO, GREAT, SQ2PI,
     &                D11,D12,D13,D14,D15,
     &                D21,D22,D23,D24,D25,
     &                D31,D32,D33,D34,D35,
     &                D41,D42,D43,D44,D45
