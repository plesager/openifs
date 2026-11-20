! (C) Copyright 2026- KNMI.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE HAMM7_CLD ( & 
!  input
 & YDMODEL,KIDIA , KFDIA, KLON,  KLEV, &
 & PAP, PT, PGFL, &
 & PL    , PI   , PA   , &
!  output
 & PLCRIT_AER, PICRIT_AER, &
 & PRE_LIQ   , PRE_ICE, &
 & PCCN      , PNICE &
 & )    

!     HAMM7_CLD 
!     --------
!          Sets up all information for subsequent calculation 
!          of the effect of prognostic aerosols on clouds and convection
!          (2nd indirect effect)

!     AUTHOR
!          Codes built upon aer_cld.F90

!     PURPOSE.
!     --------

!     INTERFACE.
!     ----------

!          *HAMM7_CLD* IS CALLED FROM *AER_CLOUD_LAYER*

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----

! -   INPUT ARGUMENTS.
!     -------------------

! KIDIA    : START OF HORIZONTAL LOOP
! KFDIA    : END   OF HORIZONTAL LOOP
! KLON     : HORIZONTAL DIMENSION
! KLEV     : END OF VERTICAL LOOP AND VERTICAL DIMENSION

! -   OUTPUT ARGUMENTS.
!     -------------------
! PLCRIT_AER : critical liquid mmr for autoconversion process
! PICRIT_AER : critical liquid mmr for autoconversion process
! PRE_LIQ : liq Re
! PRE_ICE : ice Re
! PCCN    : liquid cloud condensation nuclei
! PNICE   : ice number concentration (cf. CCN)


USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK, JPHOOK

USE YOMCST   , ONLY : RD       ,RPI
USE TYPE_MODEL, ONLY : MODEL

IMPLICIT NONE

! input variables
TYPE(MODEL)       ,INTENT(IN)    :: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PL(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PI(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PA(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGFL(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NDIM)
! output
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLCRIT_AER(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PICRIT_AER(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRE_LIQ(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRE_ICE(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCCN(KLON,KLEV)     ! liquid cloud condensation nuclei
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PNICE(KLON,KLEV)    ! ice number concentration (cf. CCN)
!----------------------------------------------------------------------

! local arrays

! Aerosol arrays
REAL(KIND=JPRB) :: ZRHO_ICE, ZRHO_LIQ ! density of pristine ice crystals and cloud
REAL(KIND=JPRB) :: ZRLIQ_CRIT, ZRICE_CRIT ! critical radii for autoconversion process

! general arrays
REAL(KIND=JPRB) :: ZRHO(KLON,KLEV)            ! density

! misc variables
REAL(KIND=JPRB) :: ZEPSEC
INTEGER(KIND=JPIM) :: ISRCCN
INTEGER(KIND=JPIM) :: JK, JL

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!------------------------

!#include "fcttre.func.h"
!#include "fccld.func.h"

!--------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('HAMM7_CLD',0,ZHOOK_HANDLE)
ASSOCIATE(YGFL=>YDMODEL%YRML_GCONF%YGFL, YDECLDP=>YDMODEL%YRML_PHY_EC%YRECLDP)
ASSOCIATE(LAERICEAUTO=>YDECLDP%LAERICEAUTO, LAERICESED=>YDECLDP%LAERICESED, &
 & LAERLIQAUTOLSP=>YDECLDP%LAERLIQAUTOLSP, LAERLIQCOLL=>YDECLDP%LAERLIQCOLL, &
 & RCLDMAX=>YDECLDP%RCLDMAX)  

!###############################################################
!                1. Basic variables
!###############################################################

ZEPSEC =1.E-10_JPRB

! move to cldp module
ZRHO_ICE=900.0_JPRB
ZRHO_LIQ=1000.0_JPRB
ZRICE_CRIT=60.0E-6_JPRB ! ice to snow critical radius
ZRLIQ_CRIT=9.3E-6_JPRB  ! cloud to rain critical radius

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZRHO(JL,JK) = PAP(JL,JK)/(RD*PT(JL,JK))

    PCCN(JL,JK)  = MAX(PGFL(JL,JK,YGFL%YCDNC%MP9_PH), 0.0_JPRB)
    PRE_LIQ(JL,JK) = MAX(PGFL(JL,JK,YGFL%YRE_LIQ%MP9_PH),0.0_JPRB)

    PNICE(JL,JK) = MAX(PGFL(JL,JK,YGFL%YICNC%MP9_PH), 0.0_JPRB)
    PRE_ICE(JL,JK) = MAX(1.0E+06_JPRB*PGFL(JL,JK,YGFL%YRE_ICE%MP9_PH),1.0_JPRB)
  ENDDO
ENDDO

!###############################################################
!                2. WARM PHASE MICROPHYSICS
!###############################################################

IF (LAERLIQAUTOLSP.OR.LAERLIQCOLL) THEN
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ISRCCN=2
      ! liquid critical mmr (Rotstayn & Penner 2001), same formula as AER_CLD
      IF (PCCN(JL,JK) > 0.0_JPRB) THEN
        PLCRIT_AER(JL,JK) = 1.333E6_JPRB*RPI*ZRHO_LIQ*PCCN(JL,JK)*ZRLIQ_CRIT**3/ZRHO(JL,JK)
        IF (ISRCCN == 2) THEN
          PLCRIT_AER(JL,JK) = MAX(PLCRIT_AER(JL,JK), 0.1_JPRB*YDECLDP%RCLCRIT)
          PLCRIT_AER(JL,JK) = MIN(PLCRIT_AER(JL,JK), 10.0_JPRB*YDECLDP%RCLCRIT)
        ENDIF
      ELSE
        PLCRIT_AER(JL,JK) = YDECLDP%RCLCRIT
      ENDIF
    ENDDO
  ENDDO
ENDIF

!###############################################################
!                3. ICE PHASE MICROPHYSICS
!###############################################################
IF (LAERICESED.OR.LAERICEAUTO) THEN
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ! ice critical mmr, same condition structure as AER_CLD
      IF (PT(JL,JK) < 238.0_JPRB .AND. PI(JL,JK) > ZEPSEC) THEN ! or PNICE? xw
        PICRIT_AER(JL,JK) = 1.333E6_JPRB*RPI*ZRHO_ICE*PNICE(JL,JK)*ZRICE_CRIT**3/ZRHO(JL,JK)
        PICRIT_AER(JL,JK) = MAX(PICRIT_AER(JL,JK), 0.1_JPRB*YDECLDP%RLCRITSNOW)
        PICRIT_AER(JL,JK) = MIN(PICRIT_AER(JL,JK), 10.0_JPRB*YDECLDP%RLCRITSNOW)
      ELSE
        PICRIT_AER(JL,JK) = YDECLDP%RLCRITSNOW
      ENDIF
      ! CLOUDSC computes (RNICE/PNICE)**0.333 in the LAERICEAUTO branch: avoid PNICE=0
      IF (PNICE(JL,JK) <= 0.0_JPRB) PNICE(JL,JK) = YDECLDP%RNICE
    ENDDO
  ENDDO
ENDIF

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('HAMM7_CLD',1,ZHOOK_HANDLE)
END SUBROUTINE HAMM7_CLD
