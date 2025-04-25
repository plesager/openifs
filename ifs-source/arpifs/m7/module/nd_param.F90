MODULE ND_PARAM
  !=======================================================================
  ! This is a F90 version of a F77 code
  ! *** WRITTEN BY ATHANASIOS NENES
  ! *** MODIFIED BY PRASHANT KUMAR AND ATHANASIOS NENES
  ! *** MODIFIED FOR EC-EARTH3 BY TWAN VAN NOIJE AND ATHANASIOS NENES
  ! *** MODIFIED FOR OpenIFS 48R1 BY PHILIPPE LE SAGER (KNMI)
  !=======================================================================

  !---Inherited functions, types, variables and constants 
  USE PARKIND1, ONLY: JPIM, JPRB

  IMPLICIT NONE
  PRIVATE
  
  !---public member functions
  PUBLIC ND_PARAM_SETUP  
  PUBLIC CCNSPEC
  PUBLIC PDFACTIV
  
  ! Three points is enough for PDF integration using Gauss-Legendre quadrature
  INTEGER(KIND=JPIM), PARAMETER :: NPGAUSS=3
  REAL(KIND=JPRB) :: XGS(NPGAUSS), WGS(NPGAUSS)      ! Points and weights for GAUSS integration

  ! Maximum number of lognormal modes, set to three for M7 (but basically should be .GE. "NSOL-1 = 3" from TM5M7_DATA)
  INTEGER(KIND=JPIM), PARAMETER :: NSMX = 3

  TYPE, PUBLIC :: NDPARAM
     REAL(KIND=JPRB) :: A     ! Default FHH adsorption parameters (in the case of FHH-AT) [READ]
     REAL(KIND=JPRB) :: B     ! See Kumar et al., (2011) ACP                              [READ]
     REAL(KIND=JPRB) :: ACCOM ! Accommodation coefficient                                 [READ]
     ! Parcel props
     REAL(KIND=JPRB) :: TEMP  ! Temperture (K)                                            [READ]
     REAL(KIND=JPRB) :: PRES  ! Pressure (Pa)                                             [READ]
     REAL(KIND=JPRB) :: ALFA  !                                                           [COMPUTED]
     REAL(KIND=JPRB) :: AKOH  ! Kelvin parameter                                          [COMPUTED]
     REAL(KIND=JPRB) :: SURT  ! Surface Tension for water (J m-2)                         [COMPUTED]
     ! Aerosol params
     INTEGER(KIND=JPIM) :: NMD        ! Number of lognormal modes effectively used in calculations [READ]. Ideally should be used to allocate DPG/SIG/DPC/TP, for now just be sure that NMD .LE. NSMX
     REAL(KIND=JPRB)    :: DPG(NSMX)  ! Modal diameter (m)               [READ for 1:NMD]
     REAL(KIND=JPRB)    :: SIG(NSMX)  ! Geometric dispersion (sigma_g)   [READ for 1:NMD]
     REAL(KIND=JPRB)    :: TP(NSMX)   ! Number concentration (#/m3)      [READ for 1:NMD]
     INTEGER(KIND=JPIM) :: MODE(NSMX) ! Kohler mode                      [READ for 1:NMD]
     REAL(KIND=JPRB)    :: DPC(NSMX)  ! Critical particle diameter       [COMPUTED for 1:NMD]
     REAL(KIND=JPRB)    :: SG(NSMX)   !                                  [COMPUTED for 1:NMD]
  END TYPE NDPARAM

  ! Physical constants
  REAL(KIND=JPRB), PARAMETER :: AMA  = 29.E-3_JPRB    ! Air molecular weight
  REAL(KIND=JPRB), PARAMETER :: GRAV = 9.81_JPRB      ! g constant
  REAL(KIND=JPRB), PARAMETER :: RGAS = 8.31_JPRB      ! Universal gas constant
  REAL(KIND=JPRB), PARAMETER :: Dw   = 2.75E-10_JPRB  ! Water Molecule Diameter
  REAL(KIND=JPRB), PARAMETER :: AMW  = 18.E-3_JPRB    ! Water molecular weight
  REAL(KIND=JPRB), PARAMETER :: DENW = 1E+3_JPRB      ! Water density
  REAL(KIND=JPRB), PARAMETER :: DHV  = 2.25E+6_JPRB   ! Water enthalpy of vaporization
  REAL(KIND=JPRB), PARAMETER :: CPAIR= 1.0061E+3_JPRB ! Air Cp

  ! Data for FHH exponent calculation
  ! for C1
  REAL(KIND=JPRB), PARAMETER :: D11 = -0.1907_JPRB
  REAL(KIND=JPRB), PARAMETER :: D12 = -1.6929_JPRB
  REAL(KIND=JPRB), PARAMETER :: D13 = 1.4963_JPRB
  REAL(KIND=JPRB), PARAMETER :: D14 = -0.5644_JPRB 
  REAL(KIND=JPRB), PARAMETER :: D15 = 0.0711_JPRB
  ! for C2
  REAL(KIND=JPRB), PARAMETER :: D21 = -3.9310_JPRB
  REAL(KIND=JPRB), PARAMETER :: D22 = 7.0906_JPRB
  REAL(KIND=JPRB), PARAMETER :: D23 = -5.3436_JPRB
  REAL(KIND=JPRB), PARAMETER :: D24 = 1.8025_JPRB 
  REAL(KIND=JPRB), PARAMETER :: D25 = -0.2131_JPRB
  ! for C3
  REAL(KIND=JPRB), PARAMETER :: D31 = 8.4825_JPRB
  REAL(KIND=JPRB), PARAMETER :: D32 = -14.9297_JPRB
  REAL(KIND=JPRB), PARAMETER :: D33 = 11.4552_JPRB
  REAL(KIND=JPRB), PARAMETER :: D34 = -3.9115_JPRB 
  REAL(KIND=JPRB), PARAMETER :: D35 = 0.4647_JPRB
  ! for C4
  REAL(KIND=JPRB), PARAMETER :: D41 = -5.1774_JPRB
  REAL(KIND=JPRB), PARAMETER :: D42 = 8.8725_JPRB
  REAL(KIND=JPRB), PARAMETER :: D43 = -6.8527_JPRB
  REAL(KIND=JPRB), PARAMETER :: D44 = 2.3514_JPRB 
  REAL(KIND=JPRB), PARAMETER :: D45 = -0.2799_JPRB
  !
  INTEGER(KIND=JPIM), PARAMETER :: MAXIT = 30        ! Max iterations for solution
  REAL(KIND=JPRB),    PARAMETER :: EPS = 1.E-5_JPRB  ! Convergence criterion
  !
  REAL(KIND=JPRB), PARAMETER :: PI    = 3.1415927_JPRB         ! Some constants
  REAL(KIND=JPRB), PARAMETER :: ZERO  = 0.0_JPRB
  REAL(KIND=JPRB), PARAMETER :: GREAT = 1.E+30_JPRB
  REAL(KIND=JPRB), PARAMETER :: SQ2PI = 2.5066282746_JPRB
    
CONTAINS

  SUBROUTINE ND_PARAM_SETUP
    !
    ! ** CALCULATE GAUSS QUADRATURE POINTS
    !
    CALL GAULEG (XGS, WGS, NPGAUSS)
    
  END SUBROUTINE ND_PARAM_SETUP
  
!=======================================================================
!
! *** SUBROUTINE CCNSPEC
! *** THIS SUBROUTINE CALCULATES THE CCN SPECTRUM OF THE AEROSOL USING
!     THE APPROPRIATE FORM OF KOHLER THEORY
!
! *** ORIGINALLY WRITTEN BY ATHANASIOS NENES FOR ONLY KOHLER PARTICLES
! *** MODIFIED BY PRASHANT KUMAR AND ATHANSIOS NENES TO INCLUDE 
! *** ACTIVATION BY FHH PARTICLES
!
!=======================================================================

  SUBROUTINE CCNSPEC (TPI,DPGI,SIGI,MODEI,TPARC,PPARC,NMODES, &
       AKKI,A,B,ACCOM,BOX)

    REAL(KIND=JPRB),    INTENT(IN) :: TPI(NMODES)
    REAL(KIND=JPRB),    INTENT(IN) :: DPGI(NMODES)
    REAL(KIND=JPRB),    INTENT(IN) :: SIGI(NMODES)
    INTEGER(KIND=JPIM), INTENT(IN) :: MODEI(NMODES)
    REAL(KIND=JPRB),    INTENT(IN) :: TPARC
    REAL(KIND=JPRB),    INTENT(IN) :: PPARC
    INTEGER(KIND=JPIM), INTENT(IN) :: NMODES
    REAL(KIND=JPRB),    INTENT(IN) :: AKKI(NMODES), A, B, ACCOM

    TYPE(NDPARAM),     INTENT(OUT) :: BOX
    
    REAL(KIND=JPRB), PARAMETER :: ZEPS = 1E-9_JPRB  !eehol: small value for diameter (1e-8 = 0.01 um)
    REAL(KIND=JPRB)    :: Dpcm, PAR1, PAR2
    INTEGER(KIND=JPIM) :: I,K
    
#include "abor1.intfb.h"
    
    ! Check
    IF (NMODES > NSMX ) THEN
      CALL ABOR1('CCNSPEC: NMODES must be .LE. NSMX. Increase NSMX')
    ENDIF

    ! Store aerosol params, parcel properties, etc.
    !
    BOX%A = A         ! Default FHH adsorption parameters (in the case of FHH-AT)
    BOX%B = B         ! See Kumar et al., (2011) ACP
    BOX%ACCOM = ACCOM ! Accommodation coefficient
    BOX%NMD = NMODES
    DO I=1,BOX%NMD
      BOX%MODE(I)= MODEI(I)
      BOX%DPG(I) = DPGI(I)
      BOX%SIG(I) = SIGI(I)
      BOX%TP(I)  = TPI(I)
    ENDDO

    BOX%TEMP = TPARC           ! Temperature (K)
    BOX%PRES = PPARC           ! Pressure (Pa)

    ! Thermophysical properties
    BOX%SURT = SFT(TPARC)                          ! Surface Tension for water (J m-2)
    BOX%AKOH = 4._JPRB*AMW*BOX%SURT/RGAS/BOX%TEMP/DENW ! Kelvin parameter
    BOX%ALFA = GRAV*AMW*DHV/CPAIR/RGAS/TPARC/TPARC - GRAV*AMA/RGAS/TPARC ! Need to store this intermediate variable, because it is used in SINTEGRAL
      
    DO K=1,BOX%NMD
      IF (MODEI(K).EQ.1) THEN                    ! Kohler modes
        IF (DPGI(K).GE.ZEPS) THEN !eehol: add treshold if median diam is low
          PAR1   = 4._JPRB/27._JPRB/AKKI(K)/DPGI(K)**3         
          PAR2   = SQRT(PAR1*BOX%AKOH**3)
          BOX%SG(K)  = EXP(PAR2) - 1._JPRB
        ELSE
          BOX%SG(K) = ZERO
        END IF
      ELSEIF (MODEI(K).EQ.2) THEN                ! FHH modes
        CALL DpcFHH(DPGI(K),BOX,Dpcm)
        BOX%DPC(K) = Dpcm
        BOX%SG(K)  = (BOX%AKOH/Dpcm)+(-A*(((Dpcm-DPGI(K))/(2*Dw))**(-B)))
      ENDIF
    ENDDO

    !open(unit=667, file='stuffxxx', access='append', status='unknown')
    !write(667,*) TEMP, PRES, AKOH, AMW, SURT, RGAS, DENW, SG
    !close(667)

  END SUBROUTINE CCNSPEC
  
!C=======================================================================
!C
!C *** SUBROUTINE DpcFHH
!C *** THIS SUBROUTINE CALCULATES THE CRITICAL PARTICLE DIAMETER
!C     ACCORDING TO THE FHH ADSOSPRTION ISOTHERM THEORY.
!C
!C *** WRITTEN BY PRASHANT KUMAR AND ATHANASIOS NENES
!C
!C=======================================================================

  SUBROUTINE DPCFHH(DDRY,BOX,DC)

    REAL(KIND=JPRB), INTENT(IN)  :: DDRY
    TYPE(NDPARAM),   INTENT(IN)  :: BOX
    REAL(KIND=JPRB), INTENT(OUT) :: DC

    REAL(KIND=JPRB) :: mu,mu1,mu2,mu3,X1,X2l,Dpcm,Dpcl,Dpcu,TEMP,SURT, A, B
    REAL(KIND=JPRB) :: X3l,X2u,X3u,FDpcl,FDpcu,FDpcm,X2m,X3m

    ! inputs
    SURT = BOX%SURT
    TEMP = BOX%TEMP
    A = BOX%A
    B = BOX%B
                 
    mu=(4._JPRB*SURT*AMW)/(RGAS*TEMP*DENW)
    mu1=(mu*2._JPRB*Dw)/((A*B)*((2._JPRB*Dw)**(B+1._JPRB)))
    mu2=1._JPRB/mu1
    mu3=1._JPRB-(mu2**(1._JPRB/(1._JPRB+B)))
    
    Dpcl = 0._JPRB         !Lower Limit
    Dpcu = 10.E-4_JPRB     !Upper Limit
    
    DO
          
      X1 = mu2**(1._JPRB/(1._JPRB+B))
      X2l = Dpcl**(2._JPRB/(1._JPRB+B))
      X3l = X1*X2l
      FDpcl=((Dpcl-X3l)/Ddry)-1._JPRB

      X1 = mu2**(1._JPRB/(1._JPRB+B))
      X2u = Dpcu**(2._JPRB/(1._JPRB+B))
      X3u = X1*X2u
      FDpcu=((Dpcu-X3u)/Ddry)-1._JPRB

      Dpcm = (Dpcu+Dpcl)/2._JPRB

      X1= mu2**(1._JPRB/(1._JPRB+B))
      X2m= Dpcm**(2._JPRB/(1._JPRB+B))
      X3m= X1*X2m
      FDpcm=((Dpcm-X3m)/Ddry)-1._JPRB

      IF ((FDPCL*FDPCM).LE.0._JPRB) THEN

        IF (ABS(FDPCM).LE.10E-8_JPRB) THEN
          EXIT
        ELSE
          DPCL = DPCL
          DPCU = DPCM
        END IF

      ELSE IF ((FDPCL*FDPCM).GE.0._JPRB) THEN

        IF (ABS(FDPCM).LE.10E-8_JPRB) THEN
          EXIT
        ELSE
          DPCL = DPCM
          DPCU = DPCU
        END IF

!     ELSE IF ((FDPCL*FDPCM).EQ.0) THEN
      ELSE
        EXIT
      END IF
      
    END DO

    DC = DPCM

  END SUBROUTINE DPCFHH

!C=======================================================================
!C
!C *** SUBROUTINE PDFACTIV
!C *** THIS SUBROUTINE CALCULATES THE CCN ACTIVATION FRACTION ACCORDING
!C     TO THE Nenes and Seinfeld (2003) PARAMETERIZATION, WITH
!C     MODIFICATION FOR NON-CONTUNUUM EFFECTS AS PROPOSED BY Fountoukis
!C     and Nenes (2004). THIS ROUTINE CALCULATES FOR A PDF OF
!C     UPDRAFT VELOCITIES.
!C
!C *** WRITTEN BY ATHANASIOS NENES
!C
!C=======================================================================

  SUBROUTINE PDFACTIV (WPARC,SIGW,NACT,SMAX,BOX)

    REAL(KIND=JPRB), INTENT(IN)  :: WPARC, SIGW
    REAL(KIND=JPRB), INTENT(OUT) :: NACT, SMAX
    TYPE(NDPARAM),   INTENT(IN)  :: BOX

    REAL(KIND=JPRB) :: NACTI, SMAXI, DENOM, WPI
    REAL(KIND=JPRB) :: A, B, ACCOM, PDF, WHI, WLO, SCAL, PROBI
    INTEGER(KIND=JPIM) :: I

    REAL(KIND=JPRB), PARAMETER :: PLIMT = 1.0E-3_JPRB
    !
    ! *** Single updraft case
    !
    IF (SIGW.LT.1.e-10_JPRB) THEN
      !
      ! *** Case where updraft is very small
      !
      IF (WPARC.LE.1E-6_JPRB) THEN
        SMAX  = 0.0_JPRB
        NACT  = 0.0_JPRB
        RETURN
      ENDIF

      CALL ACTIVATE (WPARC,NACT,SMAX,BOX)
      !
      ! *** PDF of updrafts
      !
    ELSE
      NACT  = ZERO
      SMAX  = ZERO
      DENOM = ZERO
      PROBI = SQRT(-2.0_JPRB*LOG(PLIMT*SIGW*SQ2PI))                  ! Probability of High Updraft limit
      WHI   = WPARC + SIGW*PROBI                                     ! Upper updrft limit
      ! No need to cut off the PDF at 0.05 m/s. 
      ! Using a lower value will change the normalization.
      ! WLO   = 0.05                                                 ! Low updrft limit
      WLO = 0.0_JPRB
      SCAL  = 0.5_JPRB*(WHI-WLO)                                     ! Scaling for updrafts
      !open(unit=667,file='pgaussxx',access='append',status='unknown')
      DO I=1,NPGAUSS
        ! Points are symmetric around zero, 
        ! so the sign of the XGS term is irrelevant.
        ! As the convention is to use a plus sign,
        ! we change the minus from the original code into a plus:
        !WPI  = WLO + SCAL*(1.0-XGS(i))                              ! Updraft
        WPI  = WLO + SCAL*(1.0_JPRB+XGS(i))                          ! Updraft

        ! Catch very small velocities using the same cutoff as above
        IF (WPI.LE.1E-6_JPRB) THEN
          SMAXI = 0.0_JPRB
          NACTI = 0.0_JPRB
        ELSE
          CALL ACTIVATE (WPI,NACTI,SMAXI,BOX)                         ! # of drops
        ENDIF
        PDF  = (1.0_JPRB/SQ2PI/SIGW)*EXP(-0.5_JPRB*((WPI-WPARC)/SIGW)**2)     ! Prob. of updrafts
        NACT = NACT + WGS(i)*(PDF*NACTI)                              ! Integral for drops
        SMAX = SMAX + WGS(i)*(PDF*SMAXI)                              ! Integral for Smax
        DENOM = DENOM + WGS(i)*PDF
        IF (PDF.LT.PLIMT) EXIT
        !write(667,*) NpGauss, i, nacti, smaxi
      ENDDO
      NACT = NACT/DENOM
      SMAX = SMAX/DENOM
      !close(667)
    ENDIF

  END SUBROUTINE PDFACTIV

!C=======================================================================
!C
!C *** SUBROUTINE ACTIVATE
!C *** THIS SUBROUTINE CALCULATES THE CCN ACTIVATION FRACTION ACCORDING
!C     TO THE Nenes and Seinfeld (2003) PARAMETERIZATION, WITH
!C     MODIFICATION FOR NON-CONTUNUUM EFFECTS AS PROPOSED BY Fountoukis
!C     and Nenes (in preparation).
!C
!C *** WRITTEN BY ATHANASIOS NENES FOR KOHLER PARTICLES
!C *** MODIFIED BY PRASHANT KUMAR AND ATHANASIOS NENES TO INCLUDE FHH 
!C     PARTICLES 
!C=======================================================================


  SUBROUTINE ACTIVATE (WPARC,NDRPL,SMAX,BOX)

    REAL(KIND=JPRB), INTENT(IN)  :: WPARC
    REAL(KIND=JPRB), INTENT(OUT) :: NDRPL, SMAX
    TYPE(NDPARAM),   INTENT(IN)  :: BOX
    
    REAL(KIND=JPRB) :: TEMP, PRES, AKOH, A, B, ACCOM
    REAL(KIND=JPRB) :: PRESA, AKA, DAIR, PSAT, DV, DBIG, DLOW, COEF, ALFA, BET1
    REAL(KIND=JPRB) :: WPARCEL, BET2, BETA, CF1, CF2
    REAL(KIND=JPRB) :: C1, C2, C3, C4, X_FHH
    REAL(KIND=JPRB) :: X1, X2, X3, Y1, Y2, Y3, SINTEG1, SINTEG2, SINTEG3
    INTEGER(KIND=JPIM) :: I
    
    ! Inputs
    TEMP = BOX%TEMP
    PRES = BOX%PRES
    AKOH = BOX%AKOH
    ALFA = BOX%ALFA
    A = BOX%A
    B = BOX%B
    ACCOM = BOX%ACCOM
    !
    ! *** Setup common block variables
    !
    PRESA = PRES/1.013E+5_JPRB                         ! Pressure (Pa)
    AKA   = (4.39_JPRB+0.071_JPRB*TEMP)*1.E-3_JPRB     ! Air thermal conductivity
    DAIR  = PRES*AMA/RGAS/TEMP                         ! Air density
    PSAT  = VPRES(TEMP)*(1.E+5_JPRB/1.0E+3_JPRB)       ! Saturation vapor pressure
    DV    = (0.211_JPRB/PRESA)*(TEMP/273._JPRB)**1.94
    DV    = DV*1.E-4_JPRB                              ! Water vapor diffusivity in air
    DBIG  = 5.0E-6_JPRB
    DLOW  = 0.207683_JPRB*((ACCOM)**(-0.33048_JPRB))
    DLOW  = DLOW*1.E-6_JPRB
    !
    ! *** Compute an average diffusivity Dv as a function of ACCOM
    !
    COEF  = ((2._JPRB*PI*AMW/(RGAS*TEMP))**0.5_JPRB)
    DV    = (DV/(DBIG-DLOW))*((DBIG-DLOW)-(2._JPRB*DV/ACCOM)*COEF* &
         (LOG((DBIG+(2._JPRB*DV/ACCOM)*COEF)/(DLOW+(2._JPRB*DV/ACCOM)* &
         COEF))))                      ! Non-continuum effects

    WPARCEL = WPARC
    !
    ! *** Setup constants
    !
    BET1 = PRES*AMA/PSAT/AMW + AMW*DHV*DHV/CPAIR/RGAS/TEMP/TEMP
    BET2 = RGAS*TEMP*DENW/PSAT/DV/AMW/4._JPRB + &
         DHV*DENW/4._JPRB/AKA/TEMP*(DHV*AMW/RGAS/TEMP - 1._JPRB)
    BETA = 0.5_JPRB*PI*BET1*DENW/BET2/ALFA/WPARC/DAIR
    CF1  = 0.5_JPRB*(((1._JPRB/BET2)/(ALFA*WPARC))**0.5_JPRB)
    CF2  = AKOH/3._JPRB
    !
    !     DETERMINATION OF EXPONENT FOR FHH PARTICLES
    !
    C1     = (D11)+(D12/A)+(D13/(A*A))+(D14/(A*A*A))+(D15/(A*A*A*A))
    C2     = (D21)+(D22/A)+(D23/(A*A))+(D24/(A*A*A))+(D25/(A*A*A*A))
    C3     = (D31)+(D32/A)+(D33/(A*A))+(D34/(A*A*A))+(D35/(A*A*A*A))
    C4     = (D41)+(D42/A)+(D43/(A*A))+(D44/(A*A*A))+(D45/(A*A*A*A))
    X_FHH  = (C1) + (C2/B) + (C3/(B*B)) + (C4/(B*B*B))
    !
    ! *** INITIAL VALUES FOR BISECTION *************************************
    !     
    X1   = 1.0E-5_JPRB   ! Min cloud supersaturation -> 0
    CALL SINTEGRAL (X1,NDRPL,WPARCEL,X_FHH,BET2, &
         SINTEG1,SINTEG2,SINTEG3,BOX)
    Y1   = (SINTEG1*CF1+SINTEG2*CF2+SINTEG3*CF1)*BETA*X1 - 1._JPRB
    !     
    X2   = 0.1_JPRB      ! MAX cloud supersaturation = 10%
    CALL SINTEGRAL (X2,NDRPL,WPARCEL,X_FHH,BET2, &
         SINTEG1,SINTEG2,SINTEG3,BOX)
    Y2   = (SINTEG1*CF1+SINTEG2*CF2+SINTEG3*CF1)*BETA*X2 - 1._JPRB
    !
    ! *** PERFORM BISECTION ************************************************
    !
    DO I=1,MAXIT
      X3   = 0.5_JPRB*(X1+X2)
      CALL SINTEGRAL (X3,NDRPL,WPARCEL,X_FHH,BET2, &
           SINTEG1,SINTEG2,SINTEG3,BOX)
      Y3 = (SINTEG1*CF1+SINTEG2*CF2+SINTEG3*CF1)*BETA*X3 - 1._JPRB
      !
      IF (SIGN(1._JPRB,Y1)*SIGN(1._JPRB,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
        Y2    = Y3
        X2    = X3
      ELSE
        Y1    = Y3
        X1    = X3
      ENDIF
      !
      IF (ABS(X2-X1) .LE. EPS*X1) EXIT

    ENDDO

    ! *** CONVERGED ; RETURN ***********************************************
    X3   = 0.5_JPRB*(X1+X2)

    CALL SINTEGRAL (X3,NDRPL,WPARCEL,X_FHH,BET2, &
         SINTEG1,SINTEG2,SINTEG3,BOX)
    Y3   = (SINTEG1*CF1+SINTEG2*CF2+SINTEG3*CF1)*BETA*X3 - 1._JPRB

    SMAX = X3

  END SUBROUTINE ACTIVATE

!C=======================================================================
!C
!C *** SUBROUTINE SINTEGRAL
!C *** THIS SUBROUTINE CALCULATES THE CONDENSATION INTEGRALS, ACCORDING
!C     TO THE POPULATION SPLITTING ALGORITHM AND THE SUBSEQUENT VERSIONS:
!C
!C       - Nenes and Seinfeld (2003)       Population Splitting
!C       - Fountoukis and Nenes (2004)     Modal formulation
!C       - Barahona and Nenes (2010)       Approach for large CCN
!C       - Morales and Nenes (2014)        Population Splitting revised
!C
!C *** WRITTEN BY ATHANASIOS NENES for Kohler Particles
!C *** MODFIFIED BY PRASHANT KUMAR AND ATHANASIOS NENES TO INCLUDE FHH
!C     PARTICLES
!C=======================================================================

  SUBROUTINE SINTEGRAL (SPAR, SUMMA, WPARCEL, XFHH, BET2, &
       SUM, SUMMAT, SUMFHH, BOX)

    REAL(KIND=JPRB), INTENT(IN) :: SPAR
    REAL(KIND=JPRB), INTENT(IN) :: WPARCEL
    REAL(KIND=JPRB), INTENT(IN) :: XFHH
    REAL(KIND=JPRB), INTENT(IN) :: BET2

    REAL(KIND=JPRB), INTENT(OUT) :: SUMMA, SUM, SUMMAT, SUMFHH

    TYPE(NDPARAM),   INTENT(IN) :: BOX
    
    REAL(KIND=JPRB) :: ND(NSMX), NDF(NSMX)
    REAL(KIND=JPRB) :: INTEG1(NSMX), INTEG2(NSMX), INTEG1F(NSMX)
    REAL(KIND=JPRB) :: ALFA, AKOH, TP(NSMX), SG(NSMX), SIG(NSMX)
    
    REAL(KIND=JPRB)  :: ERF1C,ERF2,ERF3,ERF4C,ERF5C,ERF6C,ERF4FC,ERF5F
    REAL(KIND=JPRB)  :: ORISM1, ORISM2, ORISM3, ORISM4, ORISM5,ORISM6
    REAL(KIND=JPRB)  :: INTAUX1P1, INTAUX1P2, DLGSP, DLGSP1, DLGSP2, DLGSPF, DLGSGF, DLGSG
    REAL(KIND=JPRB)  :: RATIO, SCRIT, DW3, SSPLT1, SSPLT2, DESCR, DEQ, EKTH, SQTWO

    REAL(KIND=JPRB) :: ORISM1F, ORISM2F, ORISM3F, ORISM4F, ORISM5F
    REAL(KIND=JPRB) :: ORISM6F, ORISM7F, ORISM8F, ORISM9F
    LOGICAL :: CRIT2
    INTEGER(KIND=JPIM) :: I,J

    SQTWO  = SQRT(2._JPRB)

    ! input
    ALFA = BOX%ALFA
    AKOH = BOX%AKOH
    DO I=1,BOX%NMD
      TP(I) = BOX%TP(I)
      SG(I) = BOX%SG(I)
      SIG(I) = BOX%SIG(I)
    ENDDO

    ! ** Population Splitting -- Modified by Ricardo Morales 2014

    DESCR  = 1._JPRB - (16._JPRB/9._JPRB)*ALFA*WPARCEL*BET2*(AKOH/SPAR**2)**2
    IF (DESCR.LE.0.0_JPRB) THEN
      CRIT2  = .TRUE.             
      scrit  = ((16._JPRB/9._JPRB)*ALFA*WPARCEL*BET2*(AKOH**2))**(0.25_JPRB)            ! Scrit - (only for DELTA < 0 )
      RATIO  = (2.0E+7_JPRB/3.0_JPRB)*AKOH*(SPAR**(-0.3824_JPRB)-scrit**(-0.3824_JPRB)) ! Computing sp1 and sp2 (sp1 = sp2)
      RATIO  = 1._JPRB/SQTWO + RATIO
      IF (RATIO.GT.1.0_JPRB) RATIO = 1.0_JPRB
      SSPLT2 = SPAR*RATIO
    ELSE
      CRIT2  = .FALSE.
      SSPLT1 = 0.5_JPRB*(1._JPRB-SQRT(DESCR))     ! min root --> sp1
      SSPLT2 = 0.5_JPRB*(1._JPRB+SQRT(DESCR))     ! max root --> sp2
      SSPLT1 = SQRT(SSPLT1)*SPAR                  ! Multiply ratios with Smax
      SSPLT2 = SQRT(SSPLT2)*SPAR
    ENDIF

    !
    ! *** Computing the condensation integrals I1 and I2
    !
    SUM       = 0.0_JPRB   !Contribution of integral 1 for Kohler 
    SUMMAT    = 0.0_JPRB   !Contribution of integral 2 for kohler
    SUMMA     = 0.0_JPRB   !Variable that stores all droplets
    SUMFHH    = 0.0_JPRB   !Contribution of FHH integral

    DO J = 1, BOX%NMD

      IF (SG(J).GT.0.0_JPRB) THEN !eehol: do not calculate if SG=0 or less
        IF (BOX%MODE(J).EQ.1) THEN          ! Kohler modes

          DLGSG  = LOG(SIG(J))                            !ln(sigmai)
          DLGSP  = LOG(SG(J)/SPAR)                        !ln(sg/smax)
          DLGSP2 = LOG(SG(J)/SSPLT2)                      !ln(sg/sp2)

          ORISM1 = 2._JPRB*DLGSP2/(3._JPRB*SQTWO*DLGSG)          ! u(sp2)
          ORISM2 = ORISM1 - 3._JPRB*DLGSG/(2._JPRB*SQTWO)        ! u(sp2)-3ln(sigmai)/(2sqrt(2)
          ORISM5 = 2._JPRB*DLGSP/(3._JPRB*SQTWO*DLGSG)           ! u(smax)
          ORISM3 = ORISM5 - 3._JPRB*DLGSG/(2._JPRB*SQTWO)        ! u(smax)-3ln(sigmai)/(2sqrt(2)
          DEQ    = AKOH*2._JPRB/SG(j)/3._JPRB/SQRT(3._JPRB)      ! Dp0 = Dpc/sqrt(3) - Equilibrium diameter

          ERF2   = erf(ORISM2)
          ERF3   = erf(ORISM3)

          INTEG2(J) = (EXP(9._JPRB/8._JPRB*DLGSG*DLGSG)*TP(J)/SG(J))* &
               (ERF2 - ERF3)                          ! I2(sp2,smax)

          IF (CRIT2) THEN     

            ORISM6 = (SQTWO*DLGSP2/3._JPRB/DLGSG)-(1.5_JPRB*DLGSG/SQTWO)
            ERF6C  = erfc(ORISM6)

            INTEG1(J) = 0.0_JPRB
            DW3       = TP(j)*DEQ*EXP(9._JPRB/8._JPRB*DLGSG*DLGSG)* &  ! 'inertially' limited particles
                 ERF6C*((BET2*ALFA*WPARCEL)**0.5_JPRB)

          ELSE

            EKTH    = EXP(9._JPRB/2._JPRB*DLGSG*DLGSG)
            DLGSP1  = LOG(SG(J)/SSPLT1)                      ! ln(sg/sp1)
            ORISM4  = ORISM1 + 3._JPRB*DLGSG/SQTWO           ! u(sp2) + 3ln(sigmai)/sqrt(2)
            ERF1C   = erfc(ORISM1)
            ERF4C   = erfc(ORISM4)

            intaux1p2 =  TP(J)*SPAR*(ERF1C - &
                 0.5_JPRB*((SG(J)/SPAR)**2)*EKTH*ERF4C)  ! I1(0,sp2)

            ORISM1  = 2._JPRB*DLGSP1/(3._JPRB*SQTWO*DLGSG)       ! u(sp1)
            ORISM4  = ORISM1 + 3._JPRB*DLGSG/SQTWO               ! u(sp1) + 3ln(sigmai)/sqrt(2)
            ORISM6  = (SQTWO*DLGSP1/3._JPRB/DLGSG)-(1.5_JPRB*DLGSG/SQTWO)

            ERF1C = erfc(ORISM1)
            ERF4C = erfc(ORISM4)
            ERF6C = erfc(ORISM6)

            intaux1p1 = TP(J)*SPAR*(ERF1C - &
                 0.5_JPRB*((SG(J)/SPAR)**2)*EKTH*ERF4C)    ! I1(0,sp1)

            INTEG1(J) = (intaux1p2-intaux1p1)                   ! I1(sp1,sp2) = I1(0,sp2) - I1(0,sp1)
            !
            DW3 = TP(j)*DEQ*EXP(9._JPRB/8._JPRB*DLGSG*DLGSG)* &           ! 'inertially' limited particles.
                 ERF6C*((BET2*ALFA*WPARCEL)**0.5_JPRB)

          ENDIF

          ! *** Calculate number of Drops

          ERF5C    = erfc(ORISM5)
          ! 
          Nd(J)    = (TP(J)/2.0_JPRB)*ERF5C
          SUM      = SUM    + INTEG1(J) + DW3           !SUM OF INTEGRAL 1 FOR KOHLER
          SUMMAT   = SUMMAT + INTEG2(J)                 !SUM OF INTEGRAL 2 FOR KOHLER
          SUMMA    = SUMMA  + Nd(J)                     !SUM OF ACTIVATED KOHLER PARTICLES

        ELSEIF (BOX%MODE(J).EQ.2) THEN                      ! FHH modes

          DLGSGF  = LOG(SIG(J))                          ! ln(sigma,i)
          DLGSPF  = LOG(SG(J)/SPAR)                      ! ln(sg/smax)
          ORISM1F = (SG(J)*SG(J))/(SPAR*SPAR)            ! (sg/smax)^2
          ORISM2F = EXP(2._JPRB*XFHH*XFHH*DLGSGF*DLGSGF) ! exp(term)
          ORISM3F = SQTWO*XFHH*DLGSGF                    ! sqrt(2).x.ln(sigma,i)
          ORISM4F = DLGSPF/(-1._JPRB*ORISM3F)            ! Umax
          ORISM5F = ORISM3F - ORISM4F
          ERF5F   = erf(ORISM5F)
          ORISM6F = ERF5F
          ORISM7F = ORISM6F + 1._JPRB
          ORISM8F = 0.5_JPRB*ORISM1F*ORISM2F*ORISM7F
          ERF4FC  = erfc(ORISM4F)
          ORISM9F = ORISM8F - ERF4FC

          INTEG1F(J) =-1._JPRB*TP(J)*SPAR*ORISM9F

          ! *** Calculate number of drops activated by FHH theory
          ERF4FC  = erfc(ORISM4F)

          NdF(J)  = (TP(J)/2.0_JPRB)*ERF4FC
          SUMFHH  = SUMFHH + INTEG1F(J)         !Sum of Integral 1 for FHH
          SUMMA   = SUMMA + NdF(J)              !Sum of ACTIVATED Kohler + FHH particles

        ENDIF
      ENDIF !eehol: if SG.GT.ZERO
    ENDDO
    
  END SUBROUTINE SINTEGRAL
  
  !C=======================================================================
  !C
  !C *** SUBROUTINE PROPS
  !C *** THIS SUBROUTINE CALCULATES THE THERMOPHYSICAL PROPERTIES
  !C
  !C *** WRITTEN BY ATHANASIOS NENES
  !C
  !C=======================================================================

!  SUBROUTINE PROPS
!    
!    !REAL(KIND=JPRB) :: VPRES, SFT
!
!    !PRESA = PRES/1.013d5                  ! Pressure (Pa)
!    !DAIR  = PRES*AMA/RGAS/TEMP            ! Air density
!    !AKA   = (4.39+0.071*TEMP)*1d-3        ! Air thermal conductivity
!    !PSAT  = VPRES(SNGL(TEMP))*(1e5/1.0d3) ! Saturation vapor pressure
!    !SURT  = SFT(SNGL(TEMP))               ! Surface Tension for water (J m-2)
!
!  END SUBROUTINE PROPS

  !C=======================================================================
  !C
  !C *** FUNCTION VPRES
  !C *** THIS FUNCTION CALCULATES SATURATED WATER VAPOUR PRESSURE AS A
  !C     FUNCTION OF TEMPERATURE. VALID FOR TEMPERATURES BETWEEN -50 AND
  !C     50 C.
  !C
  !C========================= ARGUMENTS / USAGE ===========================
  !C
  !C  INPUT:
  !C     [T]
  !C     REAL variable.
  !C     Ambient temperature expressed in Kelvin.
  !C  OUTPUT:
  !C     [VPRES]
  !C     REAL variable.
  !C     Saturated vapor pressure expressed in mbar.
  !C
  !C=======================================================================

  REAL(KIND=JPRB) PURE FUNCTION VPRES (T)

    REAL(KIND=JPRB), INTENT(IN) :: T
    REAL(KIND=JPRB)    :: A(0:6), TTEMP
    INTEGER(KIND=JPIM) :: I
    
    A=(/ 6.107799610E+0_JPRB, 4.436518521E-1_JPRB, 1.428945805E-2_JPRB, &
         2.650648471E-4_JPRB, 3.031240396E-6_JPRB, 2.034080948E-8_JPRB, &
         6.136820929E-11_JPRB /)

    TTEMP = T-273._JPRB
    VPRES = A(6)*TTEMP
    DO I=5,1,-1
      VPRES = (VPRES + A(I))*TTEMP
    ENDDO
    VPRES = VPRES + A(0)
    RETURN
  END FUNCTION VPRES

!C=======================================================================
!C
!C *** FUNCTION SFT
!C *** THIS FUNCTION CALCULATES WATER SURFACE TENSION AS A
!C     FUNCTION OF TEMPERATURE. VALID FOR TEMPERATURES BETWEEN -40 AND
!C     40 C.
!C
!C ======================== ARGUMENTS / USAGE ===========================
!C
!C  INPUT:
!C     [T]
!C     REAL variable.
!C     Ambient temperature expressed in Kelvin.
!C
!C  OUTPUT:
!C     [SFT]
!C     REAL variable.
!C     Surface Tension expressed in J m-2.
!C
!C=======================================================================

  REAL(KIND=JPRB) PURE FUNCTION SFT (T)
    REAL(KIND=JPRB), INTENT(IN) :: T
    REAL(KIND=JPRB) :: TPARS
    
    TPARS = T-273._JPRB
    SFT   = 0.0761_JPRB - 1.55E-4_JPRB*TPARS

    RETURN
  END FUNCTION SFT

  !C ***********************************************************************
  !C Calculation of points and weights for N point GAUSS integration
  !C ***********************************************************************
  SUBROUTINE GAULEG (X,W,N)

    INTEGER(KIND=JPIM), INTENT(IN)  :: N
    REAL(KIND=JPRB),    INTENT(OUT) :: X(N), W(N)
    
    REAL(KIND=JPRB),    PARAMETER  :: EPS=1.E-6_JPRB
    REAL(KIND=JPRB),    PARAMETER  :: X1=-1.0_JPRB, X2=1.0_JPRB

    REAL(KIND=JPRB)    :: XM, XL, Z, Z1, P1, P2, P3, PP
    INTEGER(KIND=JPIM) :: I,J,M
    !
    ! Calculation
    !
    M=(N+1)/2
    XM=0.5_JPRB*(X2+X1)
    XL=0.5_JPRB*(X2-X1)
    DO I=1,M
      Z=COS(3.141592654_JPRB*(I-.25_JPRB)/(N+.5_JPRB))
      DO
        P1=1._JPRB
        P2=0._JPRB
        DO J=1,N
          P3=P2
          P2=P1
          P1=((2._JPRB*J-1._JPRB)*Z*P2-(J-1._JPRB)*P3)/J
        ENDDO
        PP=N*(Z*P1-P2)/(Z*Z-1._JPRB)
        Z1=Z
        Z=Z1-P1/PP
        IF (ABS(Z-Z1).LE.EPS) EXIT
      END DO

      X(I)=XM-XL*Z
      X(N+1-I)=XM+XL*Z
      W(I)=2._JPRB*XL/((1._JPRB-Z*Z)*PP*PP)
      W(N+1-I)=W(I)
    ENDDO
    RETURN
  END SUBROUTINE GAULEG

!C=======================================================================
!C
!C *** REAL FUNCTION erfp
!C *** THIS SUBROUTINE CALCULATES THE ERROR FUNCTION USING A
!C *** POLYNOMIAL APPROXIMATION
!C
!C=======================================================================

!USE-ERF/ERFC-BUILTINS  REAL(KIND=JPRB) PURE FUNCTION ERFP(X)
!USE-ERF/ERFC-BUILTINS    REAL(KIND=JPRB), INTENT(IN) :: X
!USE-ERF/ERFC-BUILTINS    REAL(KIND=JPRB) :: AXX, Y
!USE-ERF/ERFC-BUILTINS    REAL(KIND=JPRB), DIMENSION(4), PARAMETER :: AA = (/0.278393_JPRB, 0.230389_JPRB, 0.000972_JPRB, 0.078108_JPRB/)
!USE-ERF/ERFC-BUILTINS
!USE-ERF/ERFC-BUILTINS    Y = ABS(X)
!USE-ERF/ERFC-BUILTINS    AXX = 1._JPRB + Y*(AA(1)+Y*(AA(2)+Y*(AA(3)+Y*AA(4))))
!USE-ERF/ERFC-BUILTINS    AXX = AXX*AXX
!USE-ERF/ERFC-BUILTINS    AXX = AXX*AXX
!USE-ERF/ERFC-BUILTINS    AXX = 1._JPRB - (1._JPRB/AXX)
!USE-ERF/ERFC-BUILTINS    IF(X.LE.0._JPRB) THEN
!USE-ERF/ERFC-BUILTINS      ERFP = -AXX
!USE-ERF/ERFC-BUILTINS    ELSE
!USE-ERF/ERFC-BUILTINS      ERFP = AXX
!USE-ERF/ERFC-BUILTINS    ENDIF
!USE-ERF/ERFC-BUILTINS    RETURN
!USE-ERF/ERFC-BUILTINS  END FUNCTION ERFP

END MODULE ND_PARAM
