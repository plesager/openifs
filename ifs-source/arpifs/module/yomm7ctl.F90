!****  YOMM7CTL
!
!     PURPOSE.
!     --------
!       MODULE CONTAINING VARIABLES FOR HAM-M7's M7CTL 
!
!     PARAMETER        DESCRIPTION                                  
!     ---------         -----------                                  !      
!
!     REFERENCE.
!     ----------
!                   
!     AUTHOR.
!     -------
!     2020-11-17   Tero Mielonen (TeMi)

MODULE YOMM7CTL

USE PARKIND1  ,ONLY : JPIM

IMPLICIT NONE
SAVE

TYPE :: TM7CTL

!Aerosol water uptake scheme: 
!NWATER = 0 Jacobson et al., JGR 1996
!       = 1 Kappa-Koehler theory based approach (Petters and Kreidenweis, ACP 2007)
     INTEGER(KIND=JPIM) :: NWATER

!Choice of the sulfate aerosol nucleation scheme:
!  NSNUCL = 0 off
!         = 1 Vehkamaeki et al., JGR 2002
!         = 2 Kazil and Lovejoy, ACP 2007
     INTEGER(KIND=JPIM) :: NSNUCL

!Choice of the organic aerosol nucleation scheme:
!NONUCL = 0 off
!         = 1 Activation nucleation, Kulmala et al., ACP 2006
!         = 2 Kinetic nucleation, Laakso et al., ACP 2004
     INTEGER(KIND=JPIM) :: NONUCL

END TYPE TM7CTL 

TYPE(TM7CTL), POINTER :: YRM7CTL => NULL()

END MODULE YOMM7CTL 
