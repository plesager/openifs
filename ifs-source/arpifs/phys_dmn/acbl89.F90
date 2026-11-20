SUBROUTINE ACBL89  ( YDPHY,YDPHY0,KIDIA,  KFDIA,    KLON,   KTDIAN, KLEV,& !eehol: YDCST taken out (first entry)
 !-----------------------------------------------------------------------
 ! - INPUT  2D .
 & PAPHI,  PAPHIF,   PAPRS,  PAPRSF, PT, PU, PV, PA, &
 & PECT,   PQV,      PQICE,  PQLI, PRI, PNLAB, PNLABCVP ,&
 ! - INPUT  1D .
 & PGZ0,   PTS, PEIS, PLSM, PKHFL, PKQFL,&
 ! - OUTPUT 2D .
 & PUSLE,  PLMECT,   PPHI3, PCHI3, PMRIFPP)

!**** *ACBL89 - CALCUL DES LONGUEURS DE MELANGE ET DE DISSIPATION,
!               D'APRES LA METHODE DE BOUGEAULT-LACARRERRE (1989),
!               MODIFIEE EN 2000 COMME DANS MESO-NH (LES EXPOSANTS
!               -2/3 ET -3/2), AVEC DES SECURITES POUR LE CLIMAT :    
!               L > MIN(RG*ALMAVE, VKARMN*ZPHIH). 

!     Sujet.
!     ------
!     - ROUTINE DE CALCUL ACTIF .
!       CALCUL DE LA LONGUEUR DE MELANGE "PLMECT=g*L_mix" ET DE 
!       L'INVERSE DE LA LONGUEUR DE DISSIPATION "PUSLE=1/(ALD*g*L_diss)"

!     - COMPUTATION OF THE MIXING LENGTH "PLMECT=g*L_mix"
!       AND THE DISSIPATIVE LENGTH "PUSLE=1/(ALD*g*L_diss)".

!**   Interface.
!     ----------
!        *CALL* *ACBL89*

!-----------------------------------------------------------------------
! WARNING: THE ENGLISH VERSION OF VARIABLES' NAMES IS TO BE READ IN THE
!          "APLPAR" CODE, EXCEPT FOR KTDIAN.
!-----------------------------------------------------------------------

! -   ARGUMENTS D'ENTREE.
!     -------------------

! - NOM DES PARAMETRES DE DIMENSIONNEMENT DE LA PHYSIQUE.

! KIDIA      : INDICE DE DEPART DES BOUCLES VECTORISEES SUR L'HORIZONT.
! KFDIA      : INDICE DE FIN DES BOUCLES VECTORISEES SUR L'HORIZONTALE.
! KLON       : DIMENSION HORIZONTALE DES TABLEAUX.
! KTDIAN     : INDICE DE DEPART DES BOUCLES VERTICALES (1 EN GENERAL)
!              POUR LES CALCULS DE TURBULENCE + NEBULOSITE.
! KLEV       : DIMENSION VERTICALE DES TABLEAUX "FULL LEVEL".

! - NOM DES VARIABLES DE LA PHYSIQUE (PAR ORDRE ALPHABETIQUE DANS CHAQUE
!   CATEGORIE).

! - 2D (0:KLEV) .

! PAPHI      : GEOPOTENTIEL AUX DEMI-NIVEAUX.
! PAPRS      : PRESSION AUX DEMI-NIVEAUX.

! - 2D (1:KLEV) .

! PAPHIF     : GEOPOTENTIEL AUX NIVEAUX DES COUCHES.
! PAPRSF     : PRESSION AUX NIVEAUX DES COUCHES.
! PT         : TEMPERATURE (APRES AJUSTEMENT CONVECTIF).
! PU         : COMPOSANTE EN X DU VENT.
! PV         : COMPOSANTE EN Y DU VENT.
! PA         : CLOUD FRACTION  
! PECT       : ENERGIE CINETIQUE TURBULENTE.
! PQV        : HUMIDITE SPECIFIQUE DE LA VAPEUR D'EAU.
! PQICE      : HUMIDITE SPECIFIQUE SOLIDE.
! PQLI       : HUMIDITE SPECIFIQUE LIQUIDE.
! PRI        : GRADIENT RICHARDSON NUMBER  
! PNLAB      : Si 1 Presence d'un nuage Shallow
! PNLABCVP   : Si 1 Presence d'un nuage Deep

! - 1D (DIAGNOSTIQUE) .

! PGZ0       : G FOIS LA LONGUEUR DE RUGOSITE COURANTE.
! PTS        : TEMPERATURE DE SURFACE
! PEIS       : EIS INDEX
! PLSM       : Land sea mask  

!-----------------------------------------------------------------------

! -   ARGUMENTS DE SORTIE.
!     --------------------

! - NOM DES VARIABLES DE LA PHYSIQUE (PAR ORDRE ALPHABETIQUE DANS CHAQUE
!   CATEGORIE).

! - 2D (1:KLEV) .

! PUSLE      : INVERSE DE LA LONGUEUR DE DISSIPATION MULTIPLIEE PAR G
! PLMECT     : UNE LONGUEUR DE MELANGE (FOIS G) POUR ACNEBR
! PPHI3      : LA FONCTION DE REDESLPERGER POUR K_T - ST. FUNCTION FOR SCALARS
! PCHI3      : ST. FUNCTION FOR MOMENTUM
! PMRIFPP    : FLUX RICHARDSON NUMBER  

!-----------------------------------------------------------------------

! -   ARGUMENTS IMPLICITES.
!     ---------------------

! COMMON/YOMCST /
! COMMON/YOMPHY0/

!-----------------------------------------------------------------------

!     Externes.
!     ---------

!     Methode.
!     --------

!     Auteur.
!     -------
!      2003-03, P. Marquet.
!              - - - - - - - - - - - - - - - - - - - - - - - - - - -
!               From the FIRST part of the old ACCOEFKE code,
!               written in 1993 in ARPEGE format by P. Lacarrere
!               from the old PERIDOT code, then tested by C. Bossuet
!               in 1997/98 using an Eulerian T42 Arpege GCM, then 
!               improved and tested by P. Marquet with the next 
!               versions of Arpege GCM (semi-lagrangian, SCM,
!               EUROCS Shallow convection case, with the use of
!               the new ideas coming from Meso-NH developments).
!              - - - - - - - - - - - - - - - - - - - - - - - - - - -

!     Modifications.
!     --------------
!      2004-06-03, P. Marquet : a limitation of ZG2L2SLD2 according to
!                  the use of ALMAVE>0, for which the minimal value of
!                  2.*ARSC1 must be granted in stable cases.
!                  If ALMAVE>0, then (ZGLMCBR)*2 is no longer equal to 
!                  2*e*Theta/[dTheta/dphi] and the MIN(_TWO_, ...) in
!                  factor of ARSC1 ensure the limit value of 2.*ARSC1,
!                  even if ZGLMCBR=ALMAVE>0.
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!      2007-10-26, E. Bazile and S. Malardel : bugs correction for Lm_down
!                  ZDLUP1 replaced by ZDLDN1 and add deltaZ to Lm_down   
!      2008-07-18  Y. Bouteloup : 1) Modification of mixing length in case of 
!                                    shallow or deep convection cloud
!                                 2) : lup(j) = max(lup(j),lup(j+1)-delta_phi)
!                                 3) : ldn(j) = max(ldn(j),ldn(j-1)-delta_phi)
!      2008-10-03 E. Bazile : optimisation du calcul L=Lup**-2/3+Ldw**-2/3
!      K. Yessad (Jul 2009): remove CDLOCK + some cleanings
!      R. El Khatib 22-Jul-2014 Vectorizations
!     R. El Khatib 22-Jun-2022 A contribution to simplify phasing after the refactoring of YOMCLI/YOMCST/YOETHF.
!-----------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST   , ONLY : RG       ,RETV     ,RKAPPA   ,RATM ,&
       & RLVTT, RLSTT, RCPD, &
       & RCPV, RD, RV, RCW, RCS !TCST  !eehol: use params straight from YOMCST
USE YOMPHY   , ONLY : TPHY
USE YOMPHY0  , ONLY : TPHY0

!-----------------------------------------------------------------------

IMPLICIT NONE

!TYPE (TCST), INTENT (IN) :: YDCST !eehol: not needed
TYPE(TPHY)        ,INTENT(IN)    :: YDPHY
TYPE(TPHY0)       ,INTENT(IN)    :: YDPHY0
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTDIAN 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PNLAB(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PNLABCVP(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHI(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHIF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPRS(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPRSF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PA(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PV(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PECT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQV(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQICE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQLI(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRI(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGZ0(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTS(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEIS(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSM(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKHFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKQFL(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUSLE(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLMECT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PPHI3(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCHI3(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PMRIFPP(KLON,KLEV)

!-----------------------------------------------------------------------

REAL(KIND=JPRB) :: ZGDZF   (KLON,KLEV),ZGDZH  (KLON,KLEV)
REAL(KIND=JPRB) :: ZTHETAH (KLON,0:KLEV),ZTHETA (KLON,KLEV),ZDTHETA (KLON,KLEV)
REAL(KIND=JPRB) :: ZTHETAP (KLON,KLEV)

REAL(KIND=JPRB) :: ZEN(KLON),ZGZBOT(KLON),ZGZTOP(KLON),ZGZBOTCVP(KLON),ZGZTOPCVP(KLON)
REAL(KIND=JPRB) :: ZGLMUP (KLON,KLEV),ZGLMDN (KLON,KLEV),ZPHIH(KLON,KLEV)
REAL(KIND=JPRB) :: ZTESTSAVE(KLON)

! Tableaux pour les sorties sur listing (1D seulement)

INTEGER(KIND=JPIM) :: JJLEV, JLEV, JLON

REAL(KIND=JPRB) ::  ZDLDN,  ZDLDN1,  ZDLDN2,   ZQV, &
 & ZDLUP,  ZDLUP1,  ZDLUP2,   ZDPHI, &
 & ZEPSX,  ZINCR,   ZG2L2SLD2, &
 & ZGLDIS, ZGLKARMN,ZGLMCBR, &
 & ZGLMINF,ZGLMIX,ZGLMIX2,  ZPHI3MAX, &
 & ZPREF,   ZQC,      ZTEST,   ZTEST0, &
 & ZTESTM, ZUSX,    ZX,       ZZDTHVL, ZZDTHVLP, &
 & ZZTHVL, ZZTHVLP, ZQCS, ZGZLCVPUP, ZGZLCVPDN, &
 & Z2SQRT2, ZLUP(KLON), ZLDN,ZLWK0,ZLWK1(KLON)
REAL(KIND=JPRB) ::  ZEPS_QT, ZEPS_PBLH, ZMHPBLH, ZEPSL, &
 & ZCPBLH, ZBL89EXP, ZUSRBL89, ZLSCALET, ZLAS_MAX, ZGLMAX, &
 & ZEPS_RIF, ZINCR_SH, ZEPS_SHEAR, ZCPDOCP, ZRDOR, ZRIG, ZDC,&
 & ZLL1, ZLL2, ZALMAVEP
REAL(KIND=JPRB) :: ZSHEAR (KLON,0:KLEV), ZTKEH(KLON,0:KLEV), &
 & ZPBLH(KLON), ZFC(KLON,KLEV), ZMC(KLON,KLEV), &
 & ZTHETA_S1(KLON,KLEV), ZDTHETA_S1 (KLON,KLEV),ZQT(KLON,KLEV),&
 & ZDQT(KLON,KLEV), ZTHETA_S1H (KLON,0:KLEV), ZQTH (KLON,0:KLEV),&
 & ZMCH(KLON,0:KLEV), ZFCH (KLON,0:KLEV), ZRVH (KLON,0:KLEV),&
 & ZALMAVE(KLON), ZTKE_LUPW2(KLON),ZTKE_LSC(KLON),ZSCOEF(KLON,KLEV),&
 & ZCOEF_STCU(KLON),ZTKE_LUPW(KLON),ZTHETA_XH (KLON,0:KLEV),& 
 & ZTHETA_X (KLON,KLEV),ZTKE_LAS_M(KLON),ZTKE_LSC_ST,ZMFUBN_SC(KLON),&
 & ZMFUBN_XX
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ACBL89',0,ZHOOK_HANDLE)
ASSOCIATE(ALMAVX=>YDPHY0%ALMAVX, ACBRPHIM=>YDPHY0%ACBRPHIM, &
 & VKARMN=>YDPHY0%VKARMN, ALD=>YDPHY0%ALD, ECTMIN=>YDPHY0%ECTMIN, &
 & ALMAVE=>YDPHY0%ALMAVE, ARSC1=>YDPHY0%ARSC1, AKN=>YDPHY0%AKN, &
 !& RATM=>YDCST%RATM, RETV=>YDCST%RETV, RG=>YDCST%RG, RKAPPA=>YDCST%RKAPPA, & !eehol: used straight from YOMCST
 & ALPHAT=>YDPHY0%ALPHAT, USHEARM=>YDPHY0%USHEARM,&
 & TKE_LSC=>YDPHY0%TKE_LSC, TKE_LSC_MST=>YDPHY0%TKE_LSC_MST,&
 & TKE_LAS_STCU=>YDPHY0%TKE_LAS_STCU, TKE_LAS_M=>YDPHY0%TKE_LAS_M,&
 & TKE_LAS_CU=>YDPHY0%TKE_LAS_CU, TKE_LAS_ST=>YDPHY0%TKE_LAS_ST,&
 & TKE_EIS_STCUL=>YDPHY0%TKE_EIS_STCUL, TKE_EIS_STCUH=>YDPHY0%TKE_EIS_STCUH,&
 & TKE_DEIS=>YDPHY0%TKE_DEIS, TKE_LCR17=>YDPHY0%TKE_LCR17,&
 & TKE_LSL_MH=>YDPHY0%TKE_LSL_MH, TKE_LSL_RPBL=>YDPHY0%TKE_LSL_RPBL,&
 & TKE_ST_R=>YDPHY0%TKE_ST_R, TKE_ST_P=>YDPHY0%TKE_ST_P,&
 & TKE_ST_C3=>YDPHY0%TKE_ST_C3,TKE_ST_RIFM=>YDPHY0%TKE_ST_RIFM,&
 & TKE_LUPW=>YDPHY0%TKE_LUPW,TKE_LUPW2S=>YDPHY0%TKE_LUPW2S,&
 & TKE_LUPWS=>YDPHY0%TKE_LUPWS,&
 & ENT_LAMBDA=>YDPHY0%ENT_LAMBDA,TKEMULT=>YDPHY0%TKEMULT,&
 & LBL_SCALE=>YDPHY%LBL_SCALE,LBL_ASP=>YDPHY%LBL_ASP,&
 & LBL_R17=>YDPHY%LBL_R17,LBL_SL=>YDPHY%LBL_SL,&
 & LBL_S2PEAK=>YDPHY%LBL_S2PEAK,LBL_NLOC=>YDPHY%LBL_NLOC,&
 & LTKE_SFRIF=>YDPHY%LTKE_SFRIF,LTKE_THS1=>YDPHY%LTKE_THS1,&
 & LBL_SQ1=>YDPHY%LBL_SQ1,&
 & LECSHAL=>YDPHY%LECSHAL, LECDEEP=>YDPHY%LECDEEP, &
 & LTKE_NOCU_STCU=>YDPHY%LTKE_NOCU_STCU)
!-----------------------------------------------------------------------

!*
!     ------------------------------------------------------------------
!     1 - CALCULS PRELIMINAIRES
!     ------------------------------------------------------------------
ZEPS_QT    = 1e-6_JPRB    !minimal QT
ZEPS_PBLH  = 100.0_JPRB   !minimal PBLH
ZMHPBLH    = 1e4_JPRB     !PBLH computation  limit
ZEPSL      = 1e-10_JPRB   !minimal Lup/Ldw
ZCPBLH     = 1.75_JPRB    !coefficient for PBLH computation
!ZBL89EXP   = 1.5_JPRB
ZBL89EXP   = 2.0_JPRB/3.0_JPRB
ZUSRBL89   = 1._JPRB/ZBL89EXP
ZLSCALET   = SQRT(SQRT((AKN/ALPHAT)))/AKN*VKARMN !coef. for L in surface layer
!ZLSCALET   = TKE_LSC !coef. for L in surface layer
ZLAS_MAX   = 1e4_JPRB     !height at witch as. L reaches its maximum
!ZLAS_MAX   = 3e3_JPRB     !height at witch as. L reaches its maximum
ZGLMAX     = 1e+4_JPRB*RG !maximla value of RG*L
ZEPS_RIF   = 1E-6_JPRB    !minimal delta in Rif computation
ZEPS_SHEAR = 1e-2_JPRB    !minimal shear * dz^2
ZTKE_LSC_ST=TKE_LSC*TKE_LSC_MST
  
ZEPSX    = ECTMIN

!   UNE CONSTANTE POUR LA FONCTION "PHI3" (Cuxart Bougeault Redels.)
!   A CONSTANT FOR "PHI3" FUNCTION        (Quart. J. Roy. Met. Soc.)

ZPHI3MAX= (1.0_JPRB-ACBRPHIM)/ACBRPHIM

!   TABLEAUX DE TRAVAIL
!   WORK ARRAYS ONCE FOR ALL

DO JLEV=KTDIAN,KLEV
  DO JLON=KIDIA,KFDIA
    ZGDZH(JLON,JLEV) = PAPHIF(JLON,JLEV)  -PAPHI (JLON,JLEV)
  ENDDO ! JLON
ENDDO   ! JLEV
DO JLEV=KTDIAN+1,KLEV
  DO JLON=KIDIA,KFDIA
    ZGDZF(JLON,JLEV) = PAPHIF(JLON,JLEV-1)-PAPHIF(JLON,JLEV)
  ENDDO ! JLON
ENDDO   ! JLEV


!!  Calcul de ZGZTOP et ZGZBOT de la cvpp et de la cvp

ZGZTOP(:) = 0.0_JPRB
ZGZBOT(:) = 100000.0_JPRB

ZGZTOPCVP(:) = 0.0_JPRB
ZGZBOTCVP(:) = 100000.0_JPRB

DO JLEV=KTDIAN,KLEV
  DO JLON=KIDIA,KFDIA
     ZGZTOP(JLON) = MAX(ZGZTOP(JLON),PAPHIF(JLON,JLEV)*PNLAB(JLON,JLEV))
     ZGZBOT(JLON) = (1.0_JPRB-PNLAB(JLON,JLEV))*ZGZBOT(JLON)&
    &             + PNLAB(JLON,JLEV)*MIN(ZGZBOT(JLON),PAPHIF(JLON,JLEV))
  ENDDO
ENDDO    

DO JLEV=KTDIAN,KLEV
  DO JLON=KIDIA,KFDIA
     ZGZTOPCVP(JLON) = MAX(ZGZTOPCVP(JLON),PAPHIF(JLON,JLEV)*PNLABCVP(JLON,JLEV))
     ZGZBOTCVP(JLON) = (1.0_JPRB-PNLABCVP(JLON,JLEV))*ZGZBOTCVP(JLON)&
    &             + PNLABCVP(JLON,JLEV)*MIN(ZGZBOTCVP(JLON),PAPHIF(JLON,JLEV))
  ENDDO
ENDDO    

!*
!     ------------------------------------------------------------------
!     IV - DEFINITION DE TABLEAUX DE TRAVAIL POUR LES CALCULS A SUIVRE.
!          (CALCUL DES TEMPERATURES POTENTIELLES SECHES ET HUMIDES)
!     ------------------------------------------------------------------

! - - - - - - - - - - -
! CALCULS DE THETA (sec)
! - - - - - - - - - - -
DO JLEV=KTDIAN,KLEV
  DO JLON=KIDIA,KFDIA
    ZPREF              = PAPRSF(JLON,JLEV)
    ZTHETA(JLON,JLEV)  = PT(JLON,JLEV)*(RATM/ZPREF)**(RKAPPA)
  ENDDO
ENDDO

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! CALCUL DE (THETA)vl = THETA * ( 1 + RETV*Qv - (Ql+Qi) )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
DO JLEV=KTDIAN,KLEV
  DO JLON=KIDIA,KFDIA
    ZQV               =   PQV (JLON,JLEV)
    ZQC               =   PQLI(JLON,JLEV)+PQICE(JLON,JLEV)
    ZTHETA_S1(JLON,JLEV)  = MAX(ZTHETA(JLON,JLEV)*EXP(ENT_LAMBDA*(ZQC+ZQV))*&
            &EXP(-(RLVTT*PQLI(JLON,JLEV)+RLSTT*PQICE(JLON,JLEV))/(RCPD*PT(JLON,JLEV))),1.0_JPRB)
    ZTHETA_X(JLON,JLEV)  = MAX(ZTHETA(JLON,JLEV)*&
            &EXP(-(RLVTT*PQLI(JLON,JLEV)+RLSTT*PQICE(JLON,JLEV))/(RCPD*PT(JLON,JLEV))),1.0_JPRB)
    ZQT(JLON,JLEV)  = MAX(ZQC+ZQV,ZEPS_QT)
    ZTHETA(JLON,JLEV) = ZTHETA(JLON,JLEV)*(1.0_JPRB+RETV*ZQV-ZQC)
    ZDC = (1.0_JPRB+RV/RD*ZQV)*(RLVTT*ZQV/&
            &(RD*PT(JLON, JLEV)*((1.0_JPRB+RETV*ZQV-ZQC))))
    !ibd IF(PQLI(JLON,JLEV).GT.1e-5_JPRB) THEN !liquid clouds only
    ZFC(JLON,JLEV) = 1.0_JPRB+PA(JLON,JLEV)*(RLVTT*(RD+(RV-RD)*ZQV)/&
            &((RCPD+(RCPV-RCPD)*ZQV)*RV*PT(JLON, JLEV)))
    !ibd ELSE
    !ibd  ZFC(JLON,JLEV) = 1.0_JPRB+0.0_JPRB*(RLVTT*(RD+(RV-RD)*ZQV)/&
    !ibd         &((RCPD+(RCPV-RCPD)*ZQV)*RV*PT(JLON, JLEV)))
    !ibd ENDIF
    ZMC(JLON,JLEV) = (1.0_JPRB+ZDC)/(1.0_JPRB+ZDC*ZFC(JLON,JLEV))
  ENDDO
ENDDO

! - - - - - - - - - - - - - - - - -
! CALCULS DE (THETA)vl/Half-level
! - - - - - - - - - - - - - - - - -
DO JLEV=KTDIAN,KLEV-1
  DO JLON=KIDIA,KFDIA
     ZTHETAH(JLON,JLEV) = 0.5_JPRB*(ZTHETA(JLON,JLEV)+ZTHETA(JLON,JLEV+1))
     ZTHETA_S1H(JLON,JLEV) = 0.5_JPRB*(ZTHETA_S1(JLON,JLEV)+ZTHETA_S1(JLON,JLEV+1))
     ZTHETA_XH(JLON,JLEV) = 0.5_JPRB*(ZTHETA_X(JLON,JLEV)+ZTHETA_X(JLON,JLEV+1))
     ZQTH(JLON,JLEV) = 0.5_JPRB*(ZQT(JLON,JLEV)+ZQT(JLON,JLEV+1))
     ZRVH(JLON,JLEV) = 0.5_JPRB*(PQV(JLON,JLEV)+PQV(JLON,JLEV+1))/ZQTH(JLON,JLEV)
     ZMCH(JLON,JLEV) = 0.5_JPRB*(ZMC(JLON,JLEV)+ZMC(JLON,JLEV+1))
     ZFCH(JLON,JLEV) = 0.5_JPRB*(ZFC(JLON,JLEV)+ZFC(JLON,JLEV+1))
     ZTKEH(JLON,JLEV) = MAX(0.5_JPRB*(PECT(JLON,JLEV)+PECT(JLON,JLEV+1)),ECTMIN)
     ZSHEAR(JLON,JLEV)=SQRT(MAX(ZEPS_SHEAR,(PU(JLON,JLEV)-PU(JLON,JLEV-1))**2&
          & +(PV(JLON,JLEV)-PV(JLON,JLEV-1))**2))
  ENDDO
ENDDO
DO JLON=KIDIA,KFDIA
   ZTHETAH(JLON,KLEV)   = PTS(JLON)*(RATM/PAPRS(JLON,KLEV))**(RKAPPA)
   ZQTH(JLON,KLEV)   = 0.0_JPRB
   ZRVH(JLON,KLEV)   = 0.0_JPRB
   ZMCH(JLON,KLEV)   = 0.0_JPRB
   ZFCH(JLON,KLEV)   = 0.0_JPRB
   ZSHEAR(JLON,KLEV)=SQRT(MAX(ZEPS_SHEAR,(PU(JLON,KLEV))**2&
        & +(PV(JLON,KLEV))**2))
   ZTKEH(JLON,KLEV) = MAX(PECT(JLON,KLEV),ECTMIN)
ENDDO

! - - - - - - - - - - - -
! CALCULS DE d(THETA)vl
! - - - - - - - - - - - -
!DO JLON=KIDIA,KFDIA
!  ZDTHETA(JLON,KTDIAN) = ZTHETA(JLON,MAX(KTDIAN-1,1))&
!                    &   -ZTHETA(JLON,MAX(KTDIAN  ,1))
!ENDDO
DO JLEV=KTDIAN+1,KLEV
  DO JLON=KIDIA,KFDIA
     ZDTHETA(JLON,JLEV) = ZTHETA(JLON,JLEV-1)-ZTHETA(JLON,JLEV)
     ZDTHETA_S1(JLON,JLEV) = ZTHETA_S1(JLON,JLEV-1)-ZTHETA_S1(JLON,JLEV)
     ZDQT(JLON,JLEV) = ZQT(JLON,JLEV-1)-ZQT(JLON,JLEV)
  ENDDO
ENDDO

!*
!     ------------------------------------------------------------------
!     V - CALCUL DES LONGUEURS DE MELANGE ET DE DISSIPATION.
!         ON ATTRIBUT A UNE PARTICULE L'ENERGIE CINETIQUE MOYENNE
!         DONT ELLE EST ISSUE, ON REGARDE SES DEPLACEMENTS MAXIMA
!         VERS LE HAUT "ZGLMUP" ET VERS LE BAS "ZGLMDN". ON FABRIQUERA
!         "PLMECT" ET "PUSLE" EN MOYENNANT CES DEUX QUANTITES.
!     ------------------------------------------------------------------

!         Pour ne pas ruiner les vectorisations on calcule les deux cas
!         ZEN(JLON)>ZINCR ==> ZDLUP1 (on atteint le niveau cible)
!         ZEN(JLON)<ZINCR ==> ZDLUP2 (on s'arrete entre 2 niveaux)
!         Et on utilise : ZTEST=1 pour garder ZDLUP1 et ZTEST=0
!         pour garder ZDLUP2.

!         De plus, il ne faut incrementer ZGLMUP(JLON) ou ZGLMDN(JLON)
!         que si ZEN(JLON) est positif au depart : ZTEST0=1., car
!         on est amene a continuer les boucles vers le haut et
!         vers le bas jusqu'a ce que TOUS les points ont atteint
!         leur seuil de flottabilite. ZEN(JLON) peut donc etre negatif
!         pour les premiers points a atteindre ce seuil...

!         Le calcul de ZDLDN2 est general sauf si ZUSX=1/ZX est proche
!         de 0 (si ZINCR -> 0). Dans ce cas, comme la fonction est
!         continue et egale a SQRT(PZEN(JLON)/ZINCR)*ZDLDN1 avec
!         0<PZEN(JLON)/ZINCR<1, il suffit de prendre une valeur
!         seuil ZEPSX pour la variable ZUSX afin d'assurer la
!         continuite de ZDLDN2(ZUSX,ZEN(JLON)/ZINCR).

!         Les divers ABS() qui interviennent dans les SQRT() ne
!         devraient pas etre necessaires, sinon pour remedier au
!         probleme de faire tous les calculs ZDLUP1/ZDLUP2 et
!         ZDLDN1/ZDLDN2 pour garantir la vectorisation, y compris
!         dans des cas ou les formules degenerent mal (mais
!         seules les formules pertinantes sont retenues ensuite
!         en fonction de ZTEST).

!     ------------------------------------------------------------------

ZQCS = 1.0E-6_JPRB

!     --------------------------------------
DO JLEV=KTDIAN,KLEV-1 ! BOUCLE GENERALE
!     --------------------------------------

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Pour les modifications "wet-BL89" de l'INM
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Calcul de (Theta)vl de la particule apres deplacement
  ! du demi-niveaux (JK) vers autres niveaux pleins, en
  ! conservant (Theta)l et Qtot au cours du deplacement.
  ! On a (THETA)vl = [ (THETA)l + Lv*THETA/Cp/T*(Ql+Qi) ]
  !                 *[ 1 + (1+RETV)*Qv - Qtot ]
  ! et ici : JLEV+1 varie entre KTDIAN+1 et KLEV (donc OK)
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  DO JJLEV=KTDIAN,KLEV
    DO JLON=KIDIA,KFDIA
       !ZTHETAP(JLON,JJLEV) = ZTHETAH(JLON, JLEV)
       ZTHETAP(JLON,JJLEV) = ZTHETA_XH(JLON, JLEV)
    ENDDO ! JLON
  ENDDO ! JJLEV

  ! - - - - - - - - - - - - - - - - - - - - -
  ! BOUCLE VERS LE HAUT QUI CALCULE ZGLMUP :
  ! - - - - - - - - - - - - - - - - - - - - -

  ! L'energie (ZEN) et valeur initiale de "Lup"

  DO JLON=KIDIA,KFDIA
    ZEN   (JLON)=PECT(JLON,JLEV)*TKEMULT
    ZGLMUP(JLON,JLEV)=0.0_JPRB
  ENDDO ! JLON

  ! Le passage du demi-niveau JLEV au niveau
  ! JLEV juste au dessus :

  DO JLON=KIDIA,KFDIA
    ZDLUP1 = ZGDZH(JLON,JLEV)
    !ZINCR  =(ZTHETA (JLON,JLEV)-ZTHETAP(JLON,JLEV))&
    ! & /ZTHETAH(JLON,JLEV)/2.0_JPRB*ZDLUP1
    ZINCR  =(ZTHETA_X (JLON,JLEV)-ZTHETAP(JLON,JLEV))&
         & /ZTHETA_XH(JLON,JLEV)/2.0_JPRB*ZDLUP1  
    ZINCR_SH = 0.5_JPRB*TKE_LCR17*ZSHEAR(JLON,JLEV+1)*SQRT(ZTKEH(JLON,JLEV+1))
    IF(LBL_R17) THEN
     !de Roddier 2017 wind shear extension of LBL89
       ZINCR    = ZINCR+ZINCR_SH
    ENDIF
    ZINCR=SIGN(1._JPRB,ZINCR)*MAX(1.E-10_JPRB,ABS(ZINCR))
    ZTEST0 =0.5_JPRB+SIGN(0.5_JPRB,ZEN(JLON))
    ZTEST  =0.5_JPRB+SIGN(0.5_JPRB,ZEN(JLON)-ZINCR)
    ZTESTSAVE(JLON) = ZTEST
    !ZDLUP2 =SQRT(ABS(ZEN(JLON)/ZINCR))*ZDLUP1
    IF(LBL_SQ1) THEN
       ZDLUP2 =SQRT(ABS(ZEN(JLON)/ZINCR))*ZDLUP1
    ELSE
       ZDLUP2 =(ABS(ZEN(JLON)/ZINCR))*ZDLUP1
    ENDIF
    ZDLUP  =ZTEST*ZDLUP1+(1.0_JPRB-ZTEST)*ZDLUP2
    ZGLMUP(JLON,JLEV)=ZGLMUP(JLON,JLEV)+ZDLUP*ZTEST0
    ZEN   (JLON)=ZEN  (JLON)-ZINCR*ZTEST0
  ENDDO
  ZTESTM=SUM(ZTESTSAVE(KIDIA:KFDIA))

  ! On boucle ensuite en passant du niveau JJLEV
  ! au niveau JJLEV-1 :

  DO JJLEV=JLEV,KTDIAN+1,-1
    IF (ZTESTM > 0.0_JPRB) THEN
      DO JLON=KIDIA,KFDIA
        ZDLUP1   = ZGDZF(JLON,JJLEV)
        !ZZTHVL   =(ZTHETA (JLON,JJLEV)+ZTHETA (JLON,JJLEV-1))/2.0_JPRB
        ZZTHVL   =(ZTHETA_X (JLON,JJLEV)+ZTHETA_X (JLON,JJLEV-1))/2.0_JPRB
        ZZTHVLP  =(ZTHETAP(JLON,JJLEV)+ZTHETAP(JLON,JJLEV-1))/2.0_JPRB
        !ZINCR    =(ZZTHVL-ZZTHVLP)/ZTHETAH(JLON,JLEV)*ZDLUP1
        ZINCR    =(ZZTHVL-ZZTHVLP)/ZTHETA_XH(JLON,JLEV)*ZDLUP1
        IF(LBL_R17) THEN
           ZINCR    = ZINCR+ZINCR_SH
        ENDIF
        ZINCR=SIGN(1._JPRB,ZINCR)*MAX(1.E-10_JPRB,ABS(ZINCR))
        ZTEST0   =0.5_JPRB+SIGN(0.5_JPRB,ZEN(JLON))
        ZTEST    =0.5_JPRB+SIGN(0.5_JPRB,ZEN(JLON)-ZINCR)
        ZTESTSAVE(JLON) = ZTEST
        !ZZDTHVL  = ZTHETA (JLON,JJLEV-1)-ZTHETA (JLON,JJLEV)
        ZZDTHVL  = ZTHETA_X (JLON,JJLEV-1)-ZTHETA_X (JLON,JJLEV)
        ZZDTHVLP = ZTHETAP(JLON,JJLEV-1)-ZTHETAP(JLON,JJLEV)
        !ZX    =(ZZDTHVL-ZZDTHVLP)/ZTHETAH(JLON,JLEV)/ZINCR*ZDLUP1
        IF(LBL_R17) THEN
           !ibd ZX    =((ZZDTHVL-ZZDTHVLP)/ZTHETAH(JLON,JLEV)+&
           !ibd  & ZINCR_SH/ZDLUP1)/ZINCR*ZDLUP1
           ZX    =((ZZDTHVL-ZZDTHVLP)/ZTHETA_XH(JLON,JLEV)+&
                & ZINCR_SH/ZDLUP1)/ZINCR*ZDLUP1
        ELSE
           !ibd ZX    =(ZZDTHVL-ZZDTHVLP)/ZTHETAH(JLON,JLEV)/ZINCR*ZDLUP1
           ZX    =(ZZDTHVL-ZZDTHVLP)/ZTHETA_XH(JLON,JLEV)/ZINCR*ZDLUP1
        ENDIF
        ZUSX  =1.0_JPRB/(SIGN(1.0_JPRB,ZX)*MAX(ZEPSX,ABS(ZX)))
        ZDLUP2=( -(ZUSX-0.5_JPRB) +SIGN(1.0_JPRB,ZX)*&
         & SQRT(ABS((ZUSX-0.5_JPRB)*(ZUSX-0.5_JPRB)&
         & +2.0_JPRB*ZUSX*ABS(ZEN(JLON)/ZINCR))) )*ZDLUP1  
        ZDLUP =ZTEST*ZDLUP1+(1.0_JPRB-ZTEST)*ZDLUP2
        ZGLMUP(JLON,JLEV)=ZGLMUP(JLON,JLEV)+ZDLUP*ZTEST0
        ZEN   (JLON)=ZEN  (JLON)-ZINCR*ZTEST0
      ENDDO
      ZTESTM=SUM(ZTESTSAVE(KIDIA:KFDIA))
    ENDIF
  ENDDO

  ! - - - - - - - - - - - - - - - - - - - -
  ! BOUCLE VERS LE BAS QUI CALCULE ZGLMDN :
  ! - - - - - - - - - - - - - - - - - - - -

  ! L'energie (ZEN) et valeur initiale de "Ldown"

  DO JLON=KIDIA,KFDIA
    ZEN   (JLON)=PECT(JLON,JLEV)*TKEMULT
    ZGLMDN(JLON,JLEV)=0.0_JPRB
  ENDDO ! JLON

  ! Le passage du demi-niveau JLEV au niveau
  ! JLEV+1 juste au dessous :

  DO JLON=KIDIA,KFDIA
    ZDLDN1 = ZGDZF  (JLON,JLEV+1)-ZGDZH (JLON,JLEV)
    !ZINCR  =(ZTHETAP(JLON,JLEV+1)-ZTHETA(JLON,JLEV+1))&
    ! & /ZTHETAH(JLON,JLEV)/2.0_JPRB*ZDLDN1
    ZINCR  =(ZTHETAP(JLON,JLEV+1)-ZTHETA_X(JLON,JLEV+1))&
         & /ZTHETA_XH(JLON,JLEV)/2.0_JPRB*ZDLDN1  
    ZINCR_SH = 0.5_JPRB*TKE_LCR17*ZSHEAR(JLON,JLEV+1)*SQRT(ZTKEH(JLON,JLEV+1))
    IF(LBL_R17) THEN
       ZINCR    = ZINCR+ZINCR_SH
    ENDIF
    ZINCR=SIGN(1._JPRB,ZINCR)*MAX(1.E-10_JPRB,ABS(ZINCR))
    ZTEST0 =0.5_JPRB+SIGN(0.5_JPRB,ZEN(JLON))
    ZTEST  =0.5_JPRB+SIGN(0.5_JPRB,ZEN(JLON)-ZINCR)
    ZTESTSAVE(JLON) = ZTEST
    !ZDLDN2 =SQRT(ABS(ZEN(JLON)/ZINCR))*ZDLDN1
    IF(LBL_SQ1) THEN
       ZDLDN2 =SQRT(ABS(ZEN(JLON)/ZINCR))*ZDLDN1
    ELSE
       ZDLDN2 =(ABS(ZEN(JLON)/ZINCR))*ZDLDN1
    ENDIF
    ZDLDN  =ZTEST*ZDLDN1+(1.0_JPRB-ZTEST)*ZDLDN2
    ZGLMDN(JLON,JLEV)=ZGLMDN(JLON,JLEV)+ZDLDN*ZTEST0
    ZEN   (JLON)=ZEN  (JLON)-ZINCR*ZTEST0
  ENDDO
  ZTESTM=SUM(ZTESTSAVE(KIDIA:KFDIA))

  ! On boucle ensuite en passant du niveau JJLEV
  ! au niveau JJLEV+1 :
  DO JJLEV=JLEV+1,KLEV-1
    IF (ZTESTM > 0.0_JPRB) THEN
      ZTESTM=0.0_JPRB
      DO JLON=KIDIA,KFDIA
        ZDLDN1   = ZGDZF(JLON,JJLEV+1)
        !ZZTHVL   =(ZTHETA (JLON,JJLEV)+ZTHETA (JLON,JJLEV+1))/2.0_JPRB
        ZZTHVL   =(ZTHETA_X (JLON,JJLEV)+ZTHETA_X (JLON,JJLEV+1))/2.0_JPRB
        ZZTHVLP  =(ZTHETAP(JLON,JJLEV)+ZTHETAP(JLON,JJLEV+1))/2.0_JPRB
        ZINCR    =(ZZTHVLP-ZZTHVL)/ZTHETA_XH(JLON,JLEV)*ZDLDN1
        IF(LBL_R17) THEN
           ZINCR    = ZINCR+ZINCR_SH
        ENDIF
        ZINCR=SIGN(1._JPRB,ZINCR)*MAX(1.E-10_JPRB,ABS(ZINCR))
        ZTEST0   =0.5_JPRB+SIGN(0.5_JPRB,ZEN(JLON))
        ZTEST    =0.5_JPRB+SIGN(0.5_JPRB,ZEN(JLON)-ZINCR)
        ZTESTSAVE(JLON) = ZTEST
        !ZZDTHVL  = ZTHETA (JLON,JJLEV)-ZTHETA (JLON,JJLEV+1)
        ZZDTHVL  = ZTHETA_X (JLON,JJLEV)-ZTHETA_X (JLON,JJLEV+1)
        ZZDTHVLP = ZTHETAP(JLON,JJLEV)-ZTHETAP(JLON,JJLEV+1)
        !ZX    =(ZZDTHVL-ZZDTHVLP)/ZTHETAH(JLON,JLEV)/ZINCR*ZDLDN1
        IF(LBL_R17) THEN
           !ibd ZX    =((ZZDTHVL-ZZDTHVLP)/ZTHETAH(JLON,JLEV)+&
           !ibd & ZINCR_SH/ZDLDN1)/ZINCR*ZDLDN1
           ZX    =((ZZDTHVL-ZZDTHVLP)/ZTHETA_XH(JLON,JLEV)+&
                & ZINCR_SH/ZDLDN1)/ZINCR*ZDLDN1
        ELSE
           !ibd ZX    =(ZZDTHVL-ZZDTHVLP)/ZTHETAH(JLON,JLEV)/ZINCR*ZDLDN1
           ZX    =(ZZDTHVL-ZZDTHVLP)/ZTHETA_XH(JLON,JLEV)/ZINCR*ZDLDN1
        ENDIF
        ZUSX  =1.0_JPRB/(SIGN(1.0_JPRB,ZX)*MAX(ZEPSX,ABS(ZX)))
        ZDLDN2=( -(ZUSX-0.5_JPRB) +SIGN(1.0_JPRB,ZX)*&
         & SQRT(ABS((ZUSX-0.5_JPRB)*(ZUSX-0.5_JPRB)&
         & +2.0_JPRB*ZUSX*ABS(ZEN(JLON)/ZINCR))) )*ZDLDN1  
        ZDLDN=ZTEST*ZDLDN1+(1.0_JPRB-ZTEST)*ZDLDN2
        ZGLMDN(JLON,JLEV)=ZGLMDN(JLON,JLEV)+ZDLDN*ZTEST0
        ZEN   (JLON)=ZEN  (JLON)-ZINCR*ZTEST0
      ENDDO
      ZTESTM=SUM(ZTESTSAVE(KIDIA:KFDIA))
    ENDIF
  ENDDO

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! CALCUL EFFECTIF DES LONGUEURS DE MELANGE : "Lmel" et "Ldiss"
  ! ET CALCUL DES COEFFICIENTS DE MELANGE "PKUROV" ET "PKTROV".
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! ON PREVOIT UNE LIMITATION INFERIEURE (PAS SUPERIEURE CAR LA
  ! LOI POUR FAIBLES "z" EST PLUTOT EN "2.8*z" QU'EN "0.4*z") :
  !    -> ALMAVE       :  DANS L'ATMOSPHERE LIBRE
  !    -> KARMANN*"z"  :  DANS LA CLP
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  
  DO JLON=KIDIA,KFDIA

    ! On limite Ldown par la hauteur :

    ZPHIH(JLON,JLEV)   =PAPHI(JLON,JLEV)-PAPHI(JLON,KLEV)+PGZ0(JLON)
    ZTEST0   =0.5_JPRB+SIGN(0.5_JPRB,ZEN(JLON))
    ZGLMDN(JLON,JLEV)=ZGLMDN(JLON,JLEV)+ZTEST0*(PAPHIF(JLON,KLEV)-PAPHI(JLON,KLEV))

    ZGLMDN(JLON,JLEV) =MIN(ZPHIH(JLON,JLEV), ZGLMDN(JLON,JLEV))
    
    ! On majore par la hauteur de cvpp et de cvp
    
    IF (LECSHAL) THEN
      ZGLMUP(JLON,JLEV) = MAX(ZGLMUP(JLON,JLEV),PNLAB(JLON,JLEV)*(ZGZTOP(JLON)-PAPHI(JLON,JLEV)))
      ZGLMDN(JLON,JLEV) = MAX(ZGLMDN(JLON,JLEV),PNLAB(JLON,JLEV)*(PAPHI(JLON,JLEV)-ZGZBOT(JLON)))
    ENDIF
    IF (LECDEEP) THEN  
      ZGZLCVPUP = MIN(RG*ALMAVX,ZGZTOPCVP(JLON)-PAPHI(JLON,JLEV))
      ZGZLCVPDN = MIN(RG*ALMAVX,PAPHI(JLON,JLEV)-ZGZBOTCVP(JLON))
      ZGLMUP(JLON,JLEV) = MAX(ZGLMUP(JLON,JLEV),PNLABCVP(JLON,JLEV)*ZGZLCVPUP)
      ZGLMDN(JLON,JLEV) = MAX(ZGLMDN(JLON,JLEV),PNLABCVP(JLON,JLEV)*ZGZLCVPDN)
    ENDIF  

  ENDDO  ! jlon  
ENDDO ! Jlev    

IF(LBL_NLOC) THEN
!     --------------------------------------
DO JLEV=KLEV-2,KTDIAN,-1 ! BOUCLE GENERALE Nr 2
!     --------------------------------------
  DO JLON=KIDIA,KFDIA

    !  On veut que ça monte et descende au moins à la même hauteur !!
    
    ZGLMUP(JLON,JLEV) = MAX(ZGLMUP(JLON,JLEV),ZGLMUP(JLON,JLEV+1)&
    &  + PAPHI(JLON,JLEV+1) - PAPHI(JLON,JLEV))

  ENDDO  ! jlon  
ENDDO ! Jlev    

!     --------------------------------------
DO JLEV=KTDIAN+1,KLEV-1 ! BOUCLE GENERALE Nr 3
!     --------------------------------------
  DO JLON=KIDIA,KFDIA

    ZGLMDN(JLON,JLEV) = MAX(ZGLMDN(JLON,JLEV),ZGLMDN(JLON,JLEV-1)&
    &  + PAPHI(JLON,JLEV)   - PAPHI(JLON,JLEV-1))

  ENDDO  ! jlon  
ENDDO ! Jlev    
END IF

Z2SQRT2=2._JPRB*SQRT(2._JPRB)

ZPBLH(:)=ZEPS_PBLH
!     --------------------------------------
DO JLEV=KTDIAN+1,KLEV-1 ! PBL Height - Bastak Duran et al 2022
!     --------------------------------------
  DO JLON=KIDIA,KFDIA
    
     IF (PAPHIF(JLON,JLEV)/RG.LT.ZMHPBLH) THEN !only under ZMHPBLH
       ZPBLH(JLON)=ZPBLH(JLON)+TKE_LSC*MAX(ZGLMUP(JLON,JLEV),ZEPSL)*&
             &(PAPHIF(JLON,JLEV-1)-PAPHIF(JLON,JLEV))/RG**2
     ENDIF

  ENDDO
ENDDO


DO JLON=KIDIA,KFDIA
  ZPBLH(JLON)=MIN(MAX(ZCPBLH*SQRT(MAX(ZPBLH(JLON),&
  & 1.0_JPRB)),1.0_JPRB),ZMHPBLH)
ENDDO

!EIS to identify StCu regime
DO JLON=KIDIA,KFDIA
  !increase mixing when shallow convection param is decreased
  IF (LTKE_NOCU_STCU) THEN
    ZMFUBN_XX=(PEIS(JLON)-TKE_EIS_STCUL)/(TKE_EIS_STCUH-TKE_EIS_STCUL)
    ZMFUBN_SC(JLON)=MIN(1.0_JPRB,MAX(0.0_JPRB,1.0_JPRB/(1.0_JPRB+20_JPRB*ZMFUBN_XX**2+50.0_JPRB*ZMFUBN_XX**3)))
    !ZMFUBN_SC(JLON)=MAX(0.0_JPRB,1.0_JPRB-MAX(0.0_JPRB,&
    !      &1.0_JPRB-5.0_JPRB*(PEIS(JLON)-TKE_EIS_STCUL)/(TKE_EIS_STCUH-TKE_EIS_STCUL)))
    ZTKE_LSC_ST=TKE_LSC+ZMFUBN_SC(JLON)*(TKE_LSC_MST-1.0_JPRB)*TKE_LSC
  ENDIF
  IF(PEIS(JLON).LT.TKE_EIS_STCUL) THEN
   ZALMAVE(JLON)=TKE_LAS_CU
   ZTKE_LAS_M(JLON)=TKE_LAS_CU
   !ibd ZTKE_LAS_M(JLON)=TKE_LAS_M
   ZTKE_LUPW2(JLON)=TKE_LUPW
   ZTKE_LUPW(JLON)=TKE_LUPW
   ZTKE_LSC(JLON)=TKE_LSC
   ZCOEF_STCU(JLON)=0.0_JPRB
  !ELSE IF (PEIS(JLON).LT.(TKE_DEIS+TKE_EIS_STCUL).AND.PLSM(JLON).LT.0.9_JPRB) THEN
  ELSE IF (PEIS(JLON).LT.(TKE_DEIS+TKE_EIS_STCUL).AND.&
          & (RCPD*PKHFL(JLON)<-0.1_JPRB.AND.-RCPD*PKHFL(JLON) < -0.25_JPRB*RLVTT*PKQFL(JLON)) ) THEN
   ZTKE_LAS_M(JLON)=TKE_LAS_M
   ZALMAVE(JLON)=TKE_LAS_CU-(TKE_LAS_CU-TKE_LAS_STCU)/&
   &       TKE_DEIS*(TKE_DEIS+TKE_EIS_STCUL-PEIS(JLON))
   ZTKE_LUPW2(JLON)=MAX(0.01_JPRB,MIN(0.99_JPRB,&
   &          TKE_LUPW-(TKE_LUPW-TKE_LUPW2S)&
   &          /TKE_DEIS*(TKE_DEIS+TKE_EIS_STCUL-PEIS(JLON))))
   ZTKE_LUPW(JLON)=MAX(0.01_JPRB,MIN(0.99_JPRB,&
   &          TKE_LUPW-(TKE_LUPW-TKE_LUPWS)&
   &          /TKE_DEIS*(TKE_DEIS+TKE_EIS_STCUL-PEIS(JLON))))
   ZTKE_LSC(JLON)=TKE_LSC-(TKE_LSC-ZTKE_LSC_ST)/&
   &       TKE_DEIS*(TKE_DEIS+TKE_EIS_STCUL-PEIS(JLON))
   ZCOEF_STCU(JLON)=0.0_JPRB-(0.0_JPRB-1.0_JPRB)/&
   &       TKE_DEIS*(TKE_DEIS+TKE_EIS_STCUL-PEIS(JLON))
  !ELSE IF (PEIS(JLON).LT.TKE_EIS_STCUH.AND.PLSM(JLON).LT.0.9_JPRB) THEN
  ELSE IF (PEIS(JLON).LT.TKE_EIS_STCUH.AND.&
          & (RCPD*PKHFL(JLON)<-0.1_JPRB.AND.-RCPD*PKHFL(JLON) < -0.25_JPRB*RLVTT*PKQFL(JLON)) ) THEN
   ZTKE_LAS_M(JLON)=TKE_LAS_M
   ZALMAVE(JLON)=TKE_LAS_STCU
   ZTKE_LUPW2(JLON)=TKE_LUPW2S
   ZTKE_LUPW(JLON)=TKE_LUPWS
   ZTKE_LSC(JLON)=ZTKE_LSC_ST
   ZCOEF_STCU(JLON)=1.0_JPRB
  !ELSE IF (PEIS(JLON).LT.(TKE_EIS_STCUH+TKE_DEIS).AND.PLSM(JLON).LT.0.9_JPRB) THEN
  ELSE IF (PEIS(JLON).LT.(TKE_EIS_STCUH+TKE_DEIS).AND.&
          & (RCPD*PKHFL(JLON)<-0.1_JPRB.AND.-RCPD*PKHFL(JLON) < -0.25_JPRB*RLVTT*PKQFL(JLON)) ) THEN
   ZTKE_LAS_M(JLON)=TKE_LAS_M
   ZALMAVE(JLON)=TKE_LAS_STCU-(TKE_LAS_STCU-TKE_LAS_ST)/&
   &       TKE_DEIS*(TKE_DEIS+TKE_EIS_STCUH-PEIS(JLON))
   ZTKE_LUPW2(JLON)=MAX(0.01_JPRB,MIN(0.99_JPRB,&
   &       TKE_LUPW2S-(TKE_LUPW2S-TKE_LUPW)/&
   &       TKE_DEIS*(TKE_DEIS+TKE_EIS_STCUH-PEIS(JLON))))
   ZTKE_LUPW(JLON)=MAX(0.01_JPRB,MIN(0.99_JPRB,&
   &       TKE_LUPWS-(TKE_LUPWS-TKE_LUPW)/&
   &       TKE_DEIS*(TKE_DEIS+TKE_EIS_STCUH-PEIS(JLON))))
   ZTKE_LSC(JLON)=TKE_LSC-(TKE_LSC-ZTKE_LSC_ST)/&
   &       TKE_DEIS*(TKE_DEIS+TKE_EIS_STCUH-PEIS(JLON))
   ZCOEF_STCU(JLON)=0.0_JPRB-(0.0_JPRB-1.0_JPRB)/&
   &       TKE_DEIS*(TKE_DEIS+TKE_EIS_STCUH-PEIS(JLON))
  ELSE
   ZTKE_LAS_M(JLON)=TKE_LAS_M
   ZALMAVE(JLON)=TKE_LAS_ST
   ZTKE_LUPW2(JLON)=TKE_LUPW
   ZTKE_LUPW(JLON)=TKE_LUPW
   ZTKE_LSC(JLON)=TKE_LSC
   ZCOEF_STCU(JLON)=0.0_JPRB
  ENDIF
ENDDO

!     --------------------------------------
DO JLEV=KTDIAN,KLEV-1 ! BOUCLE GENERALE Nr 4
!     --------------------------------------
  DO JLON=KIDIA,KFDIA
    
    ! La "nouvelle" longueur de melange de meso-NH 
    ! (devant approcher "2.8*g*z" pres de la surface) :

     ZLUP(JLON)=MAX(ZGLMUP(JLON,JLEV),1.E-10_JPRB)
     ZLDN=MAX(ZGLMDN(JLON,JLEV),1.E-10_JPRB)
     ZLWK0=ZLUP(JLON)/ZLDN
     !ZLWK1(JLON)=1._JPRB+ZLWK0**(2.0_JPRB/3._JPRB)
     ZLWK1(JLON)=1._JPRB+ZLWK0**ZBL89EXP

  ENDDO ! JLON

!DEC$ IVDEP
  DO JLON=KIDIA,KFDIA
     IF (LBL_S2PEAK) THEN
        !doubple peak on L profile for StCu
        ZLL1 = 1._JPRB/&
             & (         ZTKE_LUPW(JLON)/MAX(ZGLMUP(JLON,JLEV),ZEPSL)**ZBL89EXP+&
             &  (1._JPRB-ZTKE_LUPW(JLON))/MAX(ZGLMDN(JLON,JLEV),ZEPSL)**ZBL89EXP&
             & )**ZUSRBL89
        ZLL2 = 1._JPRB/&
             & (         ZTKE_LUPW2(JLON)/MAX(ZGLMUP(JLON,JLEV),ZEPSL)**ZBL89EXP+&
             &  (1._JPRB-ZTKE_LUPW2(JLON))/MAX(ZGLMDN(JLON,JLEV),ZEPSL)**ZBL89EXP&
             & )**ZUSRBL89
        ZGLMCBR= MAX(ZLL1,ZLL2)
     ELSE
        !ibd ZGLMCBR = Z2SQRT2*ZLUP(JLON)/(ZLWK1(JLON)*SQRT(ZLWK1(JLON)))
        ZGLMCBR=ZLUP(JLON)*(2._JPRB/ZLWK1(JLON))**ZUSRBL89
     ENDIF
     
     !ZGLMCBR=Z2SQRT2*ZLUP(JLON)/(ZLWK1(JLON)*SQRT(ZLWK1(JLON)))
     !scale LBL89
    IF(LBL_SCALE) THEN
      !ZGLMCBR = TKE_LSC*ZGLMCBR
      ZGLMCBR = ZTKE_LSC(JLON)*ZGLMCBR
    ENDIF

    IF(LBL_SL) THEN
     !Surface layer LBL89 modification
     IF ((ZPHIH(JLON,JLEV) .LT. RG*TKE_LSL_RPBL*ZPBLH(JLON)).OR.&
         (ZPHIH(JLON,JLEV) .LT. RG*TKE_LSL_MH)) THEN
      !only in surface layer
      ZGLKARMN= ZLSCALET*ZPHIH(JLON,JLEV)
     ELSE
      ZGLKARMN=ZEPSL
     ENDIF
    ELSE
     ZGLKARMN= VKARMN*ZPHIH(JLON,JLEV)
    ENDIF
  
    ! La limitation : L > LINF = MIN( 0.4*(G*z) , ALMAVE )
    !ZGLKARMN= VKARMN*ZPHIH(JLON,JLEV)
    ZGLMINF = MIN(RG*ALMAVE, ZGLKARMN)

    !ZGLMIX   = MAX(ZGLMINF, ZGLMCBR)
    !ZGLDIS   = MAX(ZGLMINF, ZGLMCBR)

    IF(LBL_ASP) THEN
      ZALMAVEP=ALMAVE-MAX(0.001_JPRB,ZLAS_MAX-ZPHIH(JLON,JLEV)/RG)*&
      & MAX(ZEPSL,ALMAVE-ZALMAVE(JLON))/ZLAS_MAX
      !ZLAS_MAX=3.5_JPRB*ZPBLH(JLON)
      !ZLAS_MIN=2.5_JPRB*ZPBLH(JLON)
      !IF(ZPHIH(JLON,JLEV)/RG.LT.ZLAS_MIN) THEN
      !  ZALMAVEP=ZALMAVE(JLON)
      !ELSE IF(ZPHIH(JLON,JLEV)/RG.LT.ZLAS_MAX) THEN
      !  ZALMAVEP=MAX(ALMAVE-(ZLAS_MAX-ZPHIH(JLON,JLEV)/RG)*&
      !  & (ALMAVE-ZALMAVE(JLON))/ZLAS_MAX,0.001_JPRB)
      !ELSE
      !  ZALMAVEP=ALMAVE
      !ENDIF
      ZGLMIX   = MIN(ZGLMAX,&
      & MAX(ZGLMINF, ZGLMCBR, RG*ZALMAVEP))
      !ZALMAVEP=ALMAVE-MAX(0.001_JPRB,ZLAS_MAX-ZPHIH(JLON,JLEV)/RG)*&
      !& MAX(ZEPSL,ALMAVE-TKE_LAS_ST)/ZLAS_MAX
      !ZALMAVEP=ALMAVE-MAX(0.001_JPRB,ZLAS_MAX-ZPHIH(JLON,JLEV)/RG)*&
      !& MAX(ZEPSL,ALMAVE-TKE_LAS_CU)/ZLAS_MAX
      ZALMAVEP=ALMAVE-MAX(0.001_JPRB,ZLAS_MAX-ZPHIH(JLON,JLEV)/RG)*&
      & MAX(ZEPSL,ALMAVE-ZTKE_LAS_M(JLON))/ZLAS_MAX
      !IF(ZPHIH(JLON,JLEV)/RG.LT.ZLAS_MIN) THEN
      !  ZALMAVEP=ZTKE_LAS_M(JLON)
      !ELSE IF(ZPHIH(JLON,JLEV)/RG.LT.ZLAS_MAX) THEN
      ! ZALMAVEP=MAX(ALMAVE-(ZLAS_MAX-ZPHIH(JLON,JLEV)/RG)*&
      ! & (ALMAVE-ZTKE_LAS_M(JLON))/(ZLAS_MAX-ZLAS_MIN),0.001_JPRB)
      !ELSE
      !  ZALMAVEP=ALMAVE
      !ENDIF
      !ibd Lasim for momentum not changed for convective regime
      !ZALMAVEP=ZTKE_LAS_M(JLON)
      ZGLMIX2   = MIN(ZGLMAX,&
      & MAX(ZGLMINF, ZGLMCBR, RG*ZALMAVEP))
      ZSCOEF(JLON,JLEV)=ZGLMIX2/ZGLMIX
      ZGLDIS   = ZGLMIX
    ELSE
      ZGLMIX   = MAX(ZGLMINF, ZGLMCBR)
      ZGLDIS   = MAX(ZGLMINF, ZGLMCBR)
    ENDIF
    
    ! Les deux longueurs : (1) de melange ; (2) de dissipation

    PLMECT(JLON,JLEV)=ZGLMIX
    PUSLE (JLON,JLEV)=1.0_JPRB/(ALD*ZGLDIS)

    ZDPHI    =ZGDZF(JLON,JLEV+1)
    
    !ZG2L2SLD2=MAX( ZPHI3MAX, ARSC1* MIN( 2.0_JPRB, ZGLMCBR*ZGLMCBR&
    ! & *ZDTHETA(JLON,JLEV+1)/ZDPHI&
    ! & /PECT(JLON,JLEV)/ZTHETAH(JLON,JLEV) )&
    ! & )  
    !PPHI3 (JLON,JLEV)=1.0_JPRB/(1.0_JPRB+ZG2L2SLD2)

    IF(LTKE_SFRIF) THEN
      !st. functions according to Bastak Duran et al. 2014
      ZCPDOCP=RCPD/(RCPD+(RCPD-RCPV)*ZQTH(JLON,JLEV+1))
      ZRDOR= RV/(RD+(RD-RV)*ZQTH(JLON,JLEV+1))
      !IF(LTKE_THS1) THEN
      IF(LTKE_THS1.AND.ZPHIH(JLON,JLEV)/RG.LT.1.5_JPRB*ZPBLH(JLON)) THEN
      !IF(LTKE_THS1.AND.ZPHIH(JLON,JLEV)/RG.LT.1.1_JPRB*ZPBLH(JLON)) THEN
        !Rif computed from entropy potential temperature formulation based on Marquet et al. 2013
        ZRIG  = (ZMCH(JLON,JLEV)*ZCPDOCP*ZDTHETA_S1(JLON,JLEV+1)/ZTHETA_S1H(JLON,JLEV)+&
                & (ZMCH(JLON,JLEV)*(ZFCH(JLON,JLEV)*(1.0_JPRB+ZRVH(JLON,JLEV))*&
                & ZRDOR-ZCPDOCP*ENT_LAMBDA)-1.0_JPRB/(1.0_JPRB-ZQTH(JLON,JLEV)))*&
                & ZDQT(JLON,JLEV+1))/&
                & ZSHEAR(JLON,JLEV+1)**2*ZDPHI
      ELSE
        ZRIG  = RG*ZDTHETA(JLON,JLEV+1)/ZTHETAH(JLON,JLEV)/ZSHEAR(JLON,JLEV+1)**2*ZDPHI/RG
      ENDIF
      PMRIFPP(JLON,JLEV)=MAX(TKE_ST_RIFM,&
      & MIN(-(SQRT((TKE_ST_R*(ZRIG*TKE_ST_C3+TKE_ST_P))**2 &
      & -4.0_JPRB*ZRIG*TKE_ST_C3*TKE_ST_R*TKE_ST_P**2) &
      & -(ZRIG*TKE_ST_C3+TKE_ST_P)*TKE_ST_R)/(2.0_JPRB*TKE_ST_P),&
      & TKE_ST_P-ZEPS_RIF))
      PPHI3 (JLON,JLEV)=(1.0_JPRB-PMRIFPP(JLON,JLEV)/TKE_ST_P)/&
      & (1.0_JPRB-PMRIFPP(JLON,JLEV))
      PCHI3 (JLON,JLEV)=(1.0_JPRB-PMRIFPP(JLON,JLEV)/TKE_ST_R)/&
      & (1.0_JPRB-PMRIFPP(JLON,JLEV))
      !ibd chi3 mod
      PCHI3 (JLON,JLEV)=PCHI3 (JLON,JLEV)*ZSCOEF(JLON,JLEV)
    ELSE
      ! La fonction "PHI3" de Redeslperger ; neutre=1./(1.+2.*ARSC1)
      ZG2L2SLD2=MAX( ZPHI3MAX, ARSC1* MIN( 2.0_JPRB, ZGLMCBR*ZGLMCBR&
       & *ZDTHETA(JLON,JLEV+1)/ZDPHI&
       & /PECT(JLON,JLEV)/ZTHETAH(JLON,JLEV) )&
       & )  
      PPHI3 (JLON,JLEV)=1.0_JPRB/(1.0_JPRB+ZG2L2SLD2)
      PCHI3 (JLON,JLEV)=1.0_JPRB
    ENDIF

  ENDDO ! JLON

!     ------------------------------
ENDDO ! JLEV (BOUCLE GENERALE)
!     ------------------------------

!*
!     --------------------------------------------------------------------
!     VI  - AU DERNIER NIVEAU : CALCULS DE LA LONGEUR DE MELANGE "PLMECT"
!           DE L'INVERSE DE LA LONGUEUR DE DISSIPATION "PUSLE" ET DE LA
!           FONCTION "PPHI3" DE REDELSPERGER.
!     --------------------------------------------------------------------

!DEC$ IVDEP
DO JLON=KIDIA,KFDIA
  PLMECT(JLON,KLEV)= 0.5_JPRB*PLMECT(JLON,KLEV-1)
  PUSLE (JLON,KLEV)= 1.0_JPRB/(ALD*0.5_JPRB*PLMECT(JLON,KLEV-1))
  PPHI3 (JLON,KLEV)= PPHI3(JLON,KLEV-1)
  PCHI3 (JLON,KLEV)= PCHI3(JLON,KLEV-1)
ENDDO ! JLON
!*
!     ------------------------------------------------------------------
!Rif in surface layer - for consistency
DO JLON=KIDIA,KFDIA
    ZRIG=PRI(JLON,KLEV)!full or half level?
    PMRIFPP(JLON,KLEV)=MAX(TKE_ST_RIFM,&
    & MIN(-(SQRT((TKE_ST_R*(ZRIG*TKE_ST_C3+TKE_ST_P))**2 &
    & -4.0_JPRB*ZRIG*TKE_ST_C3*TKE_ST_R*TKE_ST_P**2) &
    & -(ZRIG*TKE_ST_C3+TKE_ST_P)*TKE_ST_R)/(2.0_JPRB*TKE_ST_P),&
    & TKE_ST_P-ZEPS_RIF))
ENDDO
 
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ACBL89',1,ZHOOK_HANDLE)
END SUBROUTINE ACBL89
