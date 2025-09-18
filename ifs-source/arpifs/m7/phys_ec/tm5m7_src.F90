SUBROUTINE TM5M7_SRC( &
 & YDGEOMETRY, YDMODEL, KIDIA, KFDIA, KLON , KTDIA, KLEV, KTILES, KSTART, KSTEP ,KSTGLO,  &
 & KSW  , KTRAC, KAERO,                                                                   &
 & PALB , PALBD, PAPHI ,                                                                  &
 & PAERDEP, PAERLTS, PAERSCC, PAERGUST, PALTH ,                                           &
 & PAPH , PAP  , PCI , PCLAKE, PINJF, PBLH, PDELP, PGELAM, PGELAT, PGEMU, PFRTI, PHSDFOR, &
 & PLSM , PSST , PQ  , PRHO  , PSNS , PT  , PTL  , PTSPHY, PZ0M, KCHEM,                   &
 & PWIND, PWS1 ,PSOIL_TYPE,                                                               &
 & PCVL, PCVH, KTVL, KTVH,                                                                &
 & PLDAY,  PAERFLX, PCFLX , PCEN  , PTENC, PEMIDIAG, PSO2SRC,PSO4SRC )

! ╭────────────────────────────────────────────────────────────────────────────╮
! │                                                      (updated 03-Jun-2024) │
! │ Purpose :                                                                  │
! │ -------                                                                    │
! │  *tm5m7_src* -                                                             │
! │                                                                            │
! │                                                                            │
! │ Interface :                                                                │
! │ ---------                                                                  │
! │   *TM5M7_SRC* IS CALLED FROM *TM5M7_PHY2"                                  │
! │                                                                            │
! │                                                                            │
! │ Input :                                                                    │
! │ -----                                                                      │
! │                                                                            │
! │                                                                            │
! │ Output :                                                                   │
! │ ------                                                                     │
! │                                                                            │
! │                                                                            │
! │ Externals :                                                                │
! │ ---------                                                                  │
! │                                                                            │
! │ Method :                                                                   │
! │ ------                                                                     │
! │                                                                            │
! │ Reference :                                                                │
! │ ---------                                                                  │
! │                                                                            │
! │ Author :                                                                   │
! │ -------                                                                    │
! │     Orginal version:                                                       │
! │     Vicent Huijen (KNMI) - 2020-08-25                                      │
! │                                                                            │
! │ Modifications :                                                            │
! │ -------------                                                              │
! │     May.  2024 - R. Checa-Garcia: revision for CY48r1 and refactory        │
! │                                                                            │
! ╰────────────────────────────────────────────────────────────────────────────╯



! IFS model --------------------------------------------------------------------

USE GEOMETRY_MOD, ONLY : GEOMETRY
USE TYPE_MODEL,   ONLY : MODEL
USE PARKIND1,     ONLY : JPIM,   JPRB,    JPRD
USE YOMHOOK,      ONLY : LHOOK,  DR_HOOK, JPHOOK
USE YOMLUN,       ONLY : NULOUT, NULERR
USE YOMCST,       ONLY : RA,     RPI,     RDAY, RG
USE YOMRIP0,      ONLY : NINDAT, NSSSSS

! M7 modules -------------------------------------------------------------------

USE TM5M7_DATA ,     ONLY :  NMOD, MODE_NM, MODE_NM_SED, MODE_TRACERS_SED,     &
  &                    xmc, sigma_lognormal, pom_density, carbon_density,      &
  &                    mode_aii, mode_ais, mode_acs, mode_aci,iduai,iaii_n,    &
  &                    INO3_A, INH4,IMSA, issacs, isscos, iduaci, iducoi,      &
  &                    iaci_n, iacs_n, icoi_n, icos_n, mode_cos, mode_coi
USE TM5M7_EMIS_DATA, ONLY : MODAL_EMISSIONS, &
  &                    rad_emi_ff_insol,  rad_emi_ene_insol,rad_emi_ind_insol, &
  &                    rad_emi_tra_insol, rad_emi_shp_insol,rad_emi_air_insol, &
  &                    rad_emi_bf_insol,  rad_emi_bb_insol,                    &
  &                    rad_emi_ff_sol,    rad_emi_ene_sol,rad_emi_ind_sol,     &
  &                    rad_emi_tra_sol,   rad_emi_shp_sol,rad_emi_air_sol,     &
  &                    rad_emi_bf_sol,    rad_emi_bb_sol,                      &
  &                    frac_pom_sol_bf,   frac_pom_sol_bb, frac_pom_sol_ff,    &
  &                    frac_bc_sol_bf,    frac_bc_sol_bb,  frac_bc_sol_ff
USE OIFS_TO_HAM, ONLY: ind_oifs_ham !% ind_gas_OIFS

IMPLICIT NONE

!-----------------------------------------------------------------------
!*   0.1   ARGUMENTS
!          ---------

TYPE(GEOMETRY),     INTENT(IN)    :: YDGEOMETRY
TYPE(MODEL),        INTENT(INOUT) :: YDMODEL
INTEGER(KIND=JPIM), INTENT(IN)    :: KLON, KIDIA, KFDIA
INTEGER(KIND=JPIM), INTENT(IN)    :: KLEV, KTDIA, KSTGLO
INTEGER(KIND=JPIM), INTENT(IN)    :: KTILES
INTEGER(KIND=JPIM), INTENT(IN)    :: KSTEP, KSTART
INTEGER(KIND=JPIM), INTENT(IN)    :: KSW
INTEGER(KIND=JPIM), INTENT(IN)    :: KTRAC
INTEGER(KIND=JPIM), INTENT(IN)    :: KAERO(YDMODEL%YRML_GCONF%YGFL%NAERO)

REAL(KIND=JPRB),    INTENT(IN)    :: PALB(KLON), PALBD(KLON,KSW)
REAL(KIND=JPRB),    INTENT(IN)    :: PAPHI(KLON,0:KLEV), PALTH(KLON,0:KLEV)
REAL(KIND=JPRB),    INTENT(IN)    :: PAERDEP(KLON), PAERLTS(KLON), PAERSCC(KLON)
REAL(KIND=JPRB),    INTENT(IN)    :: PAERGUST(KLON), PHSDFOR(KLON)
REAL(KIND=JPRB),    INTENT(IN)    :: PAP(KLON,KLEV), PAPH(KLON,0:KLEV)
REAL(KIND=JPRB),    INTENT(IN)    :: PGELAM(KLON), PGELAT(KLON), PGEMU(KLON)
REAL(KIND=JPRB),    INTENT(IN)    :: PFRTI(KLON,KTILES)
REAL(KIND=JPRB),    INTENT(IN)    :: PCI(KLON), PCLAKE(KLON), PLSM(KLON), PSST(KLON)
REAL(KIND=JPRB),    INTENT(IN)    :: PINJF(KLON)
REAL(KIND=JPRB),    INTENT(IN)    :: PBLH(KLON)
REAL(KIND=JPRB),    INTENT(IN)    :: PDELP(KLON,KLEV)
REAL(KIND=JPRB),    INTENT(IN)    :: PQ(KLON,KLEV), PRHO(KLON,KLEV), PSNS(KLON)
REAL(KIND=JPRB),    INTENT(IN)    :: PT(KLON,KLEV)
REAL(KIND=JPRB),    INTENT(IN)    :: PTL(KLON)
REAL(KIND=JPRB),    INTENT(IN)    :: PWIND(KLON)
REAL(KIND=JPRB),    INTENT(IN)    :: PWS1(KLON)
REAL(KIND=JPRB),    INTENT(IN)    :: PSOIL_TYPE(KLON)
REAL(KIND=JPRB),    INTENT(IN)    :: PTSPHY
REAL(KIND=JPRB),    INTENT(IN)    :: PZ0M(KLON)
INTEGER(KIND=JPIM), INTENT(IN)    :: KCHEM(YDMODEL%YRML_GCONF%YGFL%NCHEM)

! RCHG -> try to understand what is 12 and 9 here.
REAL(KIND=JPRB),    INTENT(INOUT) :: PAERFLX(KLON,12,9)
REAL(KIND=JPRB),    INTENT(INOUT) :: PCFLX(KLON,KTRAC)

REAL(KIND=JPRB),    INTENT(IN)    :: PCVL(KLON), PCVH(KLON) ! Low/High vegetation cover
INTEGER(KIND=JPIM), INTENT(IN)    :: KTVL(KLON), KTVH(KLON) ! Low/High vegetation type

REAL(KIND=JPRB),    INTENT(INOUT) :: PCEN(KLON,KLEV,KTRAC)
REAL(KIND=JPRB),    INTENT(INOUT) :: PTENC(KLON,KLEV,KTRAC)
REAL(KIND=JPRB),    INTENT(INOUT) :: PLDAY(KLON)
REAL(KIND=JPRB),    INTENT(INOUT) :: PEMIDIAG(KLON,YDMODEL%YRML_GCONF%YGFL%NACTAERO)
REAL(KIND=JPRB),    INTENT(INOUT) :: PSO4SRC(KLON,KLEV),PSO2SRC(KLON,KLEV)

!*    0.5   LOCAL VARIABLES
!           ---------------

REAL(KIND=JPRB)    :: ZBCBF(KLON), ZBCFF(KLON), ZBCGF(KLON)     ! BC related 
REAL(KIND=JPRB)    :: ZOMBF(KLON), ZOMFF(KLON), ZOMGF(KLON)     ! OM related
INTEGER(KIND=JPIM) :: JAER, JK, JL, IMODE, INMODE, JN, II, JGAS
INTEGER(KIND=JPIM) :: IGLGLO, IHTST

! TM5-M7 data

! Arrays to collect emissions
TYPE(MODAL_EMISSIONS), DIMENSION(NMOD), TARGET :: EMIS_MASS
TYPE(MODAL_EMISSIONS), DIMENSION(NMOD), TARGET :: EMIS_NUMBER

REAL(KIND=JPRB) :: ZAEROCLIS(KLON,KLEV,2) 
REAL(KIND=JPRB) :: ZCFLX(KLON,KTRAC)
REAL(KIND=JPRB) :: ZFAERO(KLON,YDMODEL%YRML_GCONF%YGFL%NACTAERO)    
REAL(KIND=JPRB) :: ZAEROK(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NACTAERO)
REAL(KIND=JPRB) :: ZTAEROK(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NACTAERO)

REAL(KIND=JPRB) :: ZGLAT(KLON), ZGLON(KLON)
REAL(KIND=JPRB) :: ZHDD, ZHSS
REAL(KIND=JPRB) :: ZDETAH(KLON,KLEV), ZETA(KLON,KLEV) , ZETAH(KLON,0:KLEV)

!-- various sources
REAL(KIND=JPRB) :: ZLOCALTIM   , ZDIURN(KLON)
REAL(KIND=JPRB) :: ZWNDDU(KLON), ZWNDSS(KLON) 
REAL(KIND=JPRB) :: ZBCSOURC, ZOMSOURC
REAL(KIND=JPRB) :: ZDEGRAD  
REAL(KIND=JPRB) :: ZRWPWP , ZRWSAT 

!-- volcano-related variables
INTEGER(KIND=JPIM) :: IYY, IMM, IDD, IMDATE 
INTEGER(KIND=JPIM) :: IY0, IM0, ID0, INC, IMON(12)

REAL(KIND=JPRB)    :: ZGRDLON(KLON) , ZGELAT(KLON), ZGDLAT(KLON), ZGDLON(KLON)
REAL(KIND=JPRB)    :: ZDLAT, ZDLON
REAL(KIND=JPRB)    :: ZGRDLAT, ZGRDLAT2, ZGRDLON2, ZINCLAT, Z1GP

!-- map 
REAL(KIND=JPRB)    :: ZLAT, ZLON

!-- QnD oceanic DMS
REAL(KIND=JPRB)    :: ZCOS0, ZSIN0, ZRAD2DEG
REAL(KIND=JPRB)    :: ZGEMU(KLON), ZLATK(KLON)


REAL(KIND=JPRB) ::  numbscale_exp,       mass2numb_fact,        &
                 &  mass2numb_ff_sol,    mass2numb_ff_insol,    &
                 &  mass2numb_ene_sol,   mass2numb_ene_insol,   &
                 &  mass2numb_ind_sol,   mass2numb_ind_insol,   &
                 &  mass2numb_tra_sol,   mass2numb_tra_insol,   &
                 &  mass2numb_shp_sol,   mass2numb_shp_insol,   &
                 &  mass2numb_air_sol,   mass2numb_air_insol,   &
                 &  mass2numb_bf_sol,    mass2numb_bf_insol,    &
                 &  mass2numb_bb_sol,    mass2numb_bb_insol,    &
                 &  mass2numb_nonbf_sol, mass2numb_nonbf_insol, &
                 &  oc2pom
REAL(KIND=JPRB)    :: ZSOA(KLON)

REAL(KIND=JPRB) :: FRAC_BF(KLON), EMIT(KLON,KLEV) 

INTEGER(KIND=JPIM) :: ISSO2, ISSO4
!-- Injection height for biomass burning emissions
INTEGER(KIND=JPIM) :: ILINJ1, ILINJ2, IX(1)
REAL(KIND=JPRB)    :: ZDELP

!RCHG -> try to understand what is 5 here
REAL(KIND=JPRB)    :: ZAERMAP(KLON,5)
#ifdef __PGI
REAL(KIND=JPRB) :: ERF
#else
INTRINSIC ERF
#endif

REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "updcal.intfb.h"
#include "fcttim.func.h"

#include "surf_inq.h"

#include "tm5m7_src_ss.intfb.h"
#include "tm5m7_src_dust.intfb.h"
!#include "satur.intfb.h"
!#include "aer_volce.intfb.h"
!#include "aer_stratcl.intfb.h"

IF (LHOOK) CALL DR_HOOK('TM5M7_SRC',0,ZHOOK_HANDLE)

!-----------------------------------------------------------------------
ASSOCIATE(&
  & YDCSGLEG  => YDGEOMETRY%YRCSGLEG,            &
  & YDEPHY    => YDMODEL%YRML_PHY_EC%YREPHY,     &
  & YDEAERMAP => YDMODEL%YRML_PHY_AER%YREAERMAP, &
  & YGFL      => YDMODEL%YRML_GCONF%YGFL,        &
  & YDCOMPO   => YDMODEL%YRML_CHEM%YRCOMPO,      &
  & YDEAERSRC => YDMODEL%YRML_PHY_AER%YREAERSRC, &
  & YDRIP     => YDMODEL%YRML_GCONF%YRRIP)

ASSOCIATE(&
  & YAERO   => YGFL%YAERO,              NACTAERO  => YGFL%NACTAERO,            &
  & NAERO   => YGFL%NAERO,              NDGLG     => YDGEOMETRY%YRDIM%NDGLG,   &
  & RHGMT   => YDRIP%RHGMT,             RSTATI    => YDRIP%RSTATI,             &
  & RSIDECA => YDEAERSRC%RSIDECA,       NAERWND   => YDEAERSRC%NAERWND,        &
  & RSIVSRA => YDEAERSRC%RSIVSRA,       RCODECA   => YDEAERSRC%RCODECA,        &
  & RCOVSRA => YDEAERSRC%RCOVSRA,                                              &
  & NLOENG  => YDGEOMETRY%YRGEM%NLOENG, NGLOBALAT => YDGEOMETRY%YRMP%NGLOBALAT,&
  & YSURF   => YDEPHY%YSURF,            LAERCHEM  => YGFL%LAERCHEM)

!VH maybe 43r3, only??  
!& LAERODIU=>YDCOMPO%LAERODIU, YAERO=>YGFL%YAERO, LFIRE=>YDCOMPO%LFIRE, LINJ=>YDCOMPO%LINJ, &

! N.B.: In ECMWF model conventions, flux going upward from the surface 
! are negative
! All surface fluxes PCFLUX in kg m-2 s-1

!-----------------------------------------------------------------------

!*       0.1   TIME AND DATE OF THE MODEL
!              --------------------------
IY0=NCCAA(NINDAT)
IM0=NMM(NINDAT)
ID0=NDD(NINDAT)
INC=(NSSSSS + NINT(RSTATI)/NINT(RDAY))
CALL UPDCAL(ID0, IM0, IY0, INC,  IDD, IMM, IYY, IMON, -1)
IMDATE=IYY*10000+IMM*100+IDD
ZRAD2DEG = 180._JPRB/RPI 
!
!*       0.2   A LENGTH OF DAY INDEX
!              ---------------------
DO JL=KIDIA,KFDIA
  IGLGLO=NGLOBALAT(KSTGLO+JL-1)
  ZGEMU(JL)=YDCSGLEG%RMU(IGLGLO)                      ! sine of latitude
  ZLAT=ASIN(YDCSGLEG%RMU(IGLGLO))*ZRAD2DEG
  ZLATK(JL)=ZLAT
  ZCOS0=1._JPRB
  ZSIN0=0._JPRB
  PLDAY(JL)=MAX( RSIDECA*ZGEMU(JL)&
   & -RCODECA*RCOVSRA*SQRT(1.0_JPRB-ZGEMU(JL)**2)* ZCOS0&
   & +RCODECA*RSIVSRA*SQRT(1.0_JPRB-ZGEMU(JL)**2)* ZSIN0&
   & ,0.0_JPRB) ! PLDAY should be positive 
ENDDO

!-----------------------------------------------------------------------

!*       0.3   CONSTANTS AND ACCESS TO DATA (remove?)
!              ----------------------------

CALL SURF_INQ(YSURF,PRWPWP=ZRWPWP)
CALL SURF_INQ(YSURF,PRWSAT=ZRWSAT)



!-----------------------------------------------------------------------

!*       0.4   Time data
!              ----------------------------

ZETAH(KIDIA:KFDIA,0)=0._JPRB
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZETA(JL,JK) =PAP(JL,JK) /PAPH(JL,KLEV)
    ZETAH(JL,JK)=PAPH(JL,JK)/PAPH(JL,KLEV)
  ENDDO
ENDDO

ZDEGRAD= 180._JPRB/RPI
ZDLAT  = 180._JPRB / NDGLG      ! distance in degrees between latitude lines
ZGRDLAT= RPI / NDGLG            ! distance in radians between latitude lines
ZGRDLAT2=ZGRDLAT*0.55_JPRB

DO JL=KIDIA,KFDIA
  IGLGLO=NGLOBALAT(KSTGLO+JL-1)
  Z1GP=1.0_JPRB/REAL(NLOENG(IGLGLO),JPRB)
  ZDLON=Z1GP*2.0_JPRB*RPI      ! distance in radians between longitude points on a given latitude line
  ZGRDLON(JL)=ZDLON
  ZLAT=ASIN(YDCSGLEG%RMU(IGLGLO))             ! latitude in radians

  ZGLON(JL)=PGELAM(JL)*ZDEGRAD
  ZGLAT(JL)=ZLAT*ZDEGRAD
  ZGDLAT(JL)=ZDLAT
  ZGDLON(JL)=360._JPRB*Z1GP 

  ZLOCALTIM =RHGMT + ZGLON(JL)/360._JPRB*RDAY
  ZDIURN(JL)=COS( ((ZLOCALTIM-54000._JPRB)/RDAY) * 2._JPRB*RPI)+1._JPRB
ENDDO
!VH IF (.NOT.LAERODIU) THEN
!VH  ZDIURN(KIDIA:KFDIA)=1.0_JPRB
!VH ENDIF



!-----------------------------------------------------------------------
!*       0.5   Array initializations
!              ----------------------

IHTST=20

DO IMODE=1,NMOD
  ALLOCATE(EMIS_NUMBER(IMODE)%d3(KIDIA:KFDIA,KLEV,MODE_NM(IMODE)))
  ALLOCATE(EMIS_MASS(IMODE)%d3(KIDIA:KFDIA,KLEV,MODE_NM(IMODE)))

  EMIS_NUMBER(IMODE)%d3(KIDIA:KFDIA,1:KLEV,1:MODE_NM(IMODE))=0.0_JPRB
  EMIS_MASS(IMODE)%d3(KIDIA:KFDIA,1:KLEV,1:MODE_NM(IMODE))=0.0_JPRB
ENDDO

!DO JAER=1,NACTAERO
!   DO JL=KIDIA,KFDIA
!      PCFLX(JL,KAERO(JAER))=0._JPRB
!   END DO
!END DO

ZFAERO (KIDIA:KFDIA,         1:NACTAERO) = 0.0_JPRB
ZAEROK (KIDIA:KFDIA, 1:KLEV, 1:NACTAERO) = PCEN (KIDIA:KFDIA, 1:KLEV, KAERO(1):KAERO(NACTAERO)) 
ZTAEROK(KIDIA:KFDIA, 1:KLEV, 1:NACTAERO) = PTENC(KIDIA:KFDIA, 1:KLEV, KAERO(1):KAERO(NACTAERO))
PEMIDIAG(KIDIA:KFDIA,        1:NACTAERO) = 0.0_JPRB

! RCHG: FIXME -> there were are recurrent sematic error ARRAY(:) = 0.0_JPRB is dangerous.
ZOMBF(KIDIA:KFDIA) = 0.0_JPRB
ZOMFF(KIDIA:KFDIA) = 0.0_JPRB
ZOMGF(KIDIA:KFDIA) = 0.0_JPRB
ZBCFF(KIDIA:KFDIA) = 0.0_JPRB
ZBCBF(KIDIA:KFDIA) = 0.0_JPRB
ZBCGF(KIDIA:KFDIA) = 0.0_JPRB
ZCFLX(KIDIA:KFDIA,:) = 0.0_JPRB

!-----------------------------------------------------------------------

!*       0.6   SURFACE WIND VARIABLE RELEVANT FOR SS AND DU EMISSIONS
!              ------------------------------------------------------

IF (NAERWND == 0) THEN
!-- no gust accounted for
  ZWNDDU(KIDIA:KFDIA) = PWIND(KIDIA:KFDIA)
  ZWNDSS(KIDIA:KFDIA) = PWIND(KIDIA:KFDIA)
ELSEIF (NAERWND == 1) THEN
!-- gust only for SS, 10-m wind for DU
  ZWNDDU(KIDIA:KFDIA) = PWIND(KIDIA:KFDIA)
  ZWNDSS(KIDIA:KFDIA) = PAERGUST(KIDIA:KFDIA)
ELSEIF (NAERWND == 2) THEN
!-- gust only for DU, 10-m wind for SS
  ZWNDDU(KIDIA:KFDIA) = PAERGUST(KIDIA:KFDIA)
  ZWNDSS(KIDIA:KFDIA) = PWIND(KIDIA:KFDIA)
ELSEIF (NAERWND == 3) THEN
!-- gust for both SS and DU
  ZWNDDU(KIDIA:KFDIA) = PAERGUST(KIDIA:KFDIA)
  ZWNDSS(KIDIA:KFDIA) = PAERGUST(KIDIA:KFDIA)
ENDIF

! correction to account for the decrease of mean wind and gusts with decreasing
! time step
IF (PTSPHY < 1000) THEN
  ZWNDDU(KIDIA:KFDIA)=1.06_JPRB*ZWNDDU(KIDIA:KFDIA)
  ZWNDSS(KIDIA:KFDIA)=1.08_JPRB*ZWNDSS(KIDIA:KFDIA)
ENDIF

!-----------------------------------------------------------------------
!*   1.0   SEA SALT
!          --------

!- INFO: Simplistic lifting from surface based on 10-m wind and land-sea mask
!        (currently not used!)
!        ZHSS=8434._JPRB/1000._JPRB 

! RCHG: PCI, PCLAKE, PLSM, PSST 
CALL TM5M7_SRC_SS( KIDIA, KFDIA,  KLON, KLEV,         &
                 & PCI,   PCLAKE, PLSM, PSST, ZWNDSS, &
                 & emis_mass, emis_number )



!-----------------------------------------------------------------------
!*  2.0   DESERT DUST
!         -----------

!INFO: Simplistic lifting from surface based on 10-m wind and surface albedo
!      ZHDD=MAX(1.0_JPRB,8434._JPRB/1000._JPRB)  -> RCHG non-used and non-sense

! RCHG: define what is 1:12, 1:9 and 1:5 with a meaninful name 
!
PAERFLX(KIDIA:KFDIA,1:12,1:9)=0._JPRB
ZAERMAP(KIDIA:KFDIA,1:5)=0._JPRB
CALL TM5M7_SRC_DUST( YDEPHY, YDEAERMAP, YDEAERSRC, KIDIA, KFDIA, KLON, KLEV, KTILES, KSW,&
                   & PLSM, ZWNDDU, PSNS, PZ0M, &
                   & PAP(:,KLEV), PTL,  PSOIL_TYPE, &
                   & PFRTI, PCVL, PCVH, KTVL, KTVH, &
                   & emis_mass, emis_number, PAERFLX, ZGLON, ZGLAT, &
                   & ZRWPWP, ZRWSAT, ZAERMAP, PALB, PALBD, PWS1, PHSDFOR)

!-----------------------------------------------------------------------
!*       3.0   PARTICULATE ORGANIC MATTER
!              ---------------------------------------------------------
! CALL POM_SRC_43R3() 
!
!-----------------------------------------------------------------------
!*       4.0   BLACK CARBON
!              ------------
! CALL BC_SRC_43R3() 
!
!----------------------------------------------------------------------
!*       5.0   Convert emissions into tendencies: loop over tracers in mode
!              ------------
! CALL ADD_TOFLUX() ! -> only used for CY43R3
!
!
! RCHG -> in the case of CY48R1 the emissions non-interactive (all except SS and DUST)
!         are directly added as PCFLX and PTENC per tracer (remember there are 37 tracers)
!         Now we meed to add SS and DUST into the correct tracer identifier.
!         
!         After PCFLUX and PTENC we also fix PEMIDIAG. Note that PEMIDIAG is only 
!         used to transfer to PAERSRC array which is transfered to PGFL object to 
!         store emissions. Probably we can directly store in PAERSRC array, but 
!         I keep current implementation. 

DO JL=KIDIA,KFDIA
  PCFLX(JL,KAERO(iacs_n)) = 0.0
  PCFLX(JL,KAERO(icos_n)) = 0.0
  PCFLX(JL,KAERO(issacs)) = 0.0
  PCFLX(JL,KAERO(isscos)) = 0.0
  PCFLX(JL,KAERO(iaci_n)) = 0.0
  PCFLX(JL,KAERO(icoi_n)) = 0.0
  PCFLX(JL,KAERO(iduaci)) = 0.0
  PCFLX(JL,KAERO(iducoi)) = 0.0

  ZCFLX(JL,KAERO(iacs_n)) = emis_number(mode_acs)%d3(JL,KLEV,4)*(-1._JPRB)
  ZCFLX(JL,KAERO(icos_n)) = emis_number(mode_cos)%d3(JL,KLEV,4)*(-1._JPRB)
  ZCFLX(JL,KAERO(issacs)) = emis_mass(mode_acs)%d3(JL,KLEV,4)*(-1._JPRB)
  ZCFLX(JL,KAERO(isscos)) = emis_mass(mode_cos)%d3(JL,KLEV,4)*(-1._JPRB)
  ZCFLX(JL,KAERO(iaci_n)) = emis_number(mode_aci)%d3(JL,KLEV,1)*(-1._JPRB)
  ZCFLX(JL,KAERO(icoi_n)) = emis_number(mode_coi)%d3(JL,KLEV,1)*(-1._JPRB)
  ZCFLX(JL,KAERO(iduaci)) = emis_mass(mode_aci)%d3(JL,KLEV,1)*(-1._JPRB)
  ZCFLX(JL,KAERO(iducoi)) = emis_mass(mode_coi)%d3(JL,KLEV,1)*(-1._JPRB)

  PTENC(JL,KLEV, KAERO(iacs_n)) = PTENC(JL,KLEV, KAERO(iacs_n)) + emis_number(mode_acs)%d3(JL,KLEV,4) * RG / PDELP(JL,KLEV) 
  PTENC(JL,KLEV, KAERO(icos_n)) = PTENC(JL,KLEV, KAERO(icos_n)) + emis_number(mode_cos)%d3(JL,KLEV,4) * RG / PDELP(JL,KLEV) 
  PTENC(JL,KLEV, KAERO(issacs)) = PTENC(JL,KLEV, KAERO(issacs)) + emis_mass(mode_acs)%d3(JL,KLEV,4) * RG / PDELP(JL,KLEV) 
  PTENC(JL,KLEV, KAERO(isscos)) = PTENC(JL,KLEV, KAERO(isscos)) + emis_mass(mode_cos)%d3(JL,KLEV,4) * RG / PDELP(JL,KLEV) 

  PTENC(JL,KLEV, KAERO(iaci_n)) = PTENC(JL,KLEV, KAERO(iaci_n)) + emis_number(mode_aci)%d3(JL,KLEV,1) * RG / PDELP(JL,KLEV) 
  PTENC(JL,KLEV, KAERO(icoi_n)) = PTENC(JL,KLEV, KAERO(icoi_n)) + emis_number(mode_coi)%d3(JL,KLEV,1) * RG / PDELP(JL,KLEV) 
  PTENC(JL,KLEV, KAERO(iduaci)) = PTENC(JL,KLEV, KAERO(iduaci)) + emis_mass(mode_aci)%d3(JL,KLEV,1) * RG / PDELP(JL,KLEV) 
  PTENC(JL,KLEV, KAERO(iducoi)) = PTENC(JL,KLEV, KAERO(iducoi)) + emis_mass(mode_coi)%d3(JL,KLEV,1) * RG / PDELP(JL,KLEV) 

ENDDO

DO JL=KIDIA,KFDIA
  DO IMODE=1,NMOD                                 ! loop in each mode 
    DO INMODE=0,MODE_NM_SED(IMODE)                ! loop in aerosols species per mode 
       JN = MODE_TRACERS_SED(INMODE,IMODE)        ! retrieve indentifier of each specie
       PEMIDIAG(JL,KAERO(JN))= PEMIDIAG(JL,KAERO(JN)) + ZCFLX(JL,KAERO(JN))*(-1._JPRB) ! assign ZCFLX to emissions (we still not added dep. to PCFLX) 
    ENDDO
  ENDDO
ENDDO


!-----------------------------------------------------------------------
!*       6.0   De-allocate arrays
!              ------------

DO IMODE=1,NMOD
  IF(associated(EMIS_NUMBER(IMODE)%d3)) DEALLOCATE(EMIS_NUMBER(IMODE)%d3)
  IF(associated(EMIS_MASS(IMODE)%d3))   DEALLOCATE(EMIS_MASS(IMODE)%d3)
ENDDO

END ASSOCIATE
END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('TM5M7_SRC',1,ZHOOK_HANDLE)

CONTAINS 
SUBROUTINE POM_SRC_43R3()

    ! mass to number conversion factors for the relevant modes
    numbscale_exp  = EXP(1.5*(LOG(sigma_lognormal(mode_aii)))**2)
    mass2numb_fact = 3./(4.*RPI*(numbscale_exp**3)*pom_density)
    mass2numb_ff_insol  = mass2numb_fact/(rad_emi_ff_insol**3)
    mass2numb_ene_insol = mass2numb_fact/(rad_emi_ene_insol**3)
    mass2numb_ind_insol = mass2numb_fact/(rad_emi_ind_insol**3)
    mass2numb_tra_insol = mass2numb_fact/(rad_emi_tra_insol**3)
    mass2numb_shp_insol = mass2numb_fact/(rad_emi_shp_insol**3)
    mass2numb_air_insol = mass2numb_fact/(rad_emi_air_insol**3)
    mass2numb_bf_insol  = mass2numb_fact/(rad_emi_bf_insol**3)
    mass2numb_bb_insol  = mass2numb_fact/(rad_emi_bb_insol**3)

    numbscale_exp  = EXP(1.5*(LOG(sigma_lognormal(mode_ais)))**2)
    mass2numb_fact = 3./(4.*RPI*(numbscale_exp**3)*pom_density)
    mass2numb_ff_sol  = mass2numb_fact/(rad_emi_ff_sol**3)
    mass2numb_ene_sol = mass2numb_fact/(rad_emi_ene_sol**3)
    mass2numb_ind_sol = mass2numb_fact/(rad_emi_ind_sol**3)
    mass2numb_tra_sol = mass2numb_fact/(rad_emi_tra_sol**3)
    mass2numb_shp_sol = mass2numb_fact/(rad_emi_shp_sol**3)
    mass2numb_air_sol = mass2numb_fact/(rad_emi_air_sol**3)

    numbscale_exp  = EXP(1.5*(LOG(sigma_lognormal(mode_acs)))**2)
    mass2numb_fact = 3./(4.*RPI*(numbscale_exp**3)*pom_density)
    mass2numb_bf_sol = mass2numb_fact/(rad_emi_bf_sol**3)
    !mass2numb_bb_sol = mass2numb_fact/(rad_emi_bb_sol**3) 
    mass2numb_nonbf_sol = mass2numb_ff_sol
    mass2numb_nonbf_insol = mass2numb_ff_insol

    frac_bf(KIDIA:KFDIA)=1.0_JPRB
    ! calculate mass fraction related to solid biofuel
    where ( ZOMFF(KIDIA:KFDIA) > 1E-30_JPRB )
       frac_bf(KIDIA:KFDIA) = ZOMBF(KIDIA:KFDIA) / ZOMFF(KIDIA:KFDIA)
    elsewhere
       frac_bf(KIDIA:KFDIA) = 0.0_JPRB
    endwhere

    ! for safety, prevent fractions larger than unity.
    where (frac_bf(KIDIA:KFDIA) > 1.0_JPRB )
       frac_bf(KIDIA:KFDIA) = 1.0_JPRB
    endwhere

   ! add to emis target arrays. 
   ! For now treat all sectors identical, and put all emissions in lowest model layer (KLEV)

   ! Fossil fuel categories..
   DO JL=KIDIA,KFDIA
        emis_mass  (mode_aii)%d3(JL,KLEV,2) = &
     &  emis_mass  (mode_aii)%d3(JL,KLEV,2) + ZOMFF(JL) * &
     &    ( (1.-frac_bf(JL)) * (1.-frac_pom_sol_ff) + &
     &          frac_bf(JL)  * (1.-frac_pom_sol_bf) ) 

        emis_number(mode_aii)%d3(JL,KLEV,2) = &
     &  emis_number(mode_aii)%d3(JL,KLEV,2) + ZOMFF(JL) * &
     &    ( (1.-frac_bf(JL)) * (1.-frac_pom_sol_ff) * mass2numb_nonbf_insol + &
     &          frac_bf(JL)  * (1.-frac_pom_sol_bf) * mass2numb_bf_insol )

        emis_mass  (mode_ais)%d3(JL,KLEV,3) = &
     &  emis_mass  (mode_ais)%d3(JL,KLEV,3) + ZOMFF(JL) * &
     &    ( (1.-frac_bf(JL)) * frac_pom_sol_ff )

        emis_number(mode_ais)%d3(JL,KLEV,3) = &
     &  emis_number(mode_ais)%d3(JL,KLEV,3) + ZOMFF(JL) * &
     &    ( (1.-frac_bf(JL)) * frac_pom_sol_ff * mass2numb_nonbf_sol )

        emis_mass  (mode_acs)%d3(JL,KLEV,3) = &
     &  emis_mass  (mode_acs)%d3(JL,KLEV,3) + ZOMFF(JL) * &
     &    (        frac_bf(JL)  * frac_pom_sol_bf )

        emis_number(mode_acs)%d3(JL,KLEV,3) = &
     &  emis_number(mode_acs)%d3(JL,KLEV,3) + ZOMFF(JL) * &
     &    (        frac_bf(JL)  * frac_pom_sol_bf * mass2numb_bf_sol )
    ENDDO
    
    ! Biofuel categories ? (ZOMBF emissions)
    
!!$    IF (.not. LAERCHEM)THEN
!!$       ! SOA from CO
!VH     DO JL=KIDIA,KFDIA
!VH        ZSOA(JL)=0._JPRB
!VH        IF (LAERSOA_CHEM) THEN
!VH           ZSOA(JL)=MAX(PSOACO(JL),PSOA(JL))
!VH        ELSE
!VH           ZSOA(JL)=PSOACO(JL)
!VH        ENDIF
!VH        ZOMSOURC=ZOMSOURC+ZSOA(JL)
!VH     END DO
    ! These do not apply for M7
!!$       PCFLX(JL,KAERO(INBAER+1))= -ZOMSOURC * ROMPHIL
!!$       PCFLX(JL,KAERO(INBAER+2))= -ZOMSOURC * ROMPHOB 
!!$    END IF
    ! biomass burning

!VH   IF (LFIRE) THEN
!VH     IF (LINJ) THEN
!VH       DO JL=KIDIA,KFDIA
!VH       ! Height of injection for biomass burning emissions : update emis_mass
!VH         IF (PINJF(JL) > 200._JPRB .AND. PBLH(JL) > 1500._JPRB) THEN
!VH           IX=MINLOC( ABS( (PAPHI(JL,1:KLEV)-PAPHI(JL,KLEV))/RG - PINJF(JL)))
!VH           ILINJ1=IX(1)
!VH           ILINJ2=ILINJ1
!VH           ! calculate total deltap over injected levels
!VH           ZDELP=0.0_JPRB
!VH           DO JK = ILINJ1, ILINJ2
!VH              ZDELP = ZDELP + PDELP(JL,JK)
!VH           ENDDO
!VH           DO JK = ILINJ1, ILINJ2
!VH        
!VH            ! add to emis target arrays
!VH               emis_mass  (mode_aii)%d3(JL,JK,2) = &
!VH            &  emis_mass  (mode_aii)%d3(JL,JK,2) + ZOMGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
!VH            &    (1.-frac_pom_sol_bb)
!VH       
!VH               emis_number(mode_aii)%d3(JL,JK,2) = &
!VH            &  emis_number(mode_aii)%d3(JL,JK,2) + ZOMGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
!VH            &    (1.-frac_pom_sol_bb) * mass2numb_bb_insol 
!VH 
!VH               emis_mass  (mode_acs)%d3(JL,JK,3) = &
!VH            &  emis_mass  (mode_acs)%d3(JL,JK,3) + ZOMGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
!VH            &     frac_pom_sol_bb 
!VH 
!VH               emis_number(mode_acs)%d3(JL,JK,3) = &
!VH            &  emis_number(mode_acs)%d3(JL,JK,3) + ZOMGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
!VH            &    frac_pom_sol_bb * mass2numb_bb_sol       
!VH           ENDDO
!VH         ELSE
!VH           ZDELP=0.0_JPRB
!VH           DO JK = KLEV-3, KLEV-2
!VH                ZDELP = ZDELP + PDELP(JL,JK)
!VH           ENDDO
!VH           DO JK = KLEV-3, KLEV-2
!VH               ! add to emis target arrays
!VH               emis_mass  (mode_aii)%d3(JL,JK,2) = &
!VH            &  emis_mass  (mode_aii)%d3(JL,JK,2) + ZOMGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
!VH            &    (1.-frac_pom_sol_bb)
!VH 
!VH               emis_number(mode_aii)%d3(JL,JK,2) = &
!VH            &  emis_number(mode_aii)%d3(JL,JK,2) + ZOMGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
!VH            &    (1.-frac_pom_sol_bb) * mass2numb_bb_insol 
!VH 
!VH               emis_mass  (mode_acs)%d3(JL,JK,3) = &
!VH            &  emis_mass  (mode_acs)%d3(JL,JK,3) + ZOMGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
!VH            &    frac_pom_sol_bb 
!VH 
!VH               emis_number(mode_acs)%d3(JL,JK,3) = &
!VH            &  emis_number(mode_acs)%d3(JL,JK,3) + ZOMGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
!VH            &    frac_pom_sol_bb * mass2numb_bb_sol
!VH           ENDDO
!VH         ENDIF
!VH       ENDDO
!VH     ELSE ! LINJ=false: always injection at lowest levels
!VH       DO JL=KIDIA,KFDIA
!VH         ZDELP=0.0_JPRB
!VH         DO JK = KLEV-2, KLEV-1
!VH           ZDELP = ZDELP + PDELP(JL,JK)
!VH         ENDDO
!VH         DO JK = KLEV-2, KLEV
!VH             ! add to emis target arrays
!VH             emis_mass  (mode_aii)%d3(JL,JK,2) = &
!VH          &  emis_mass  (mode_aii)%d3(JL,JK,2) + ZOMGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
!VH          &        (1.-frac_pom_sol_bb)
!VH 
!VH             emis_number(mode_aii)%d3(JL,JK,2) = &
!VH          &  emis_number(mode_aii)%d3(JL,JK,2) + ZOMGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
!VH          &        (1.-frac_pom_sol_bb) * mass2numb_bb_insol 
!VH 
!VH             emis_mass  (mode_acs)%d3(JL,JK,3) = &
!VH          &  emis_mass  (mode_acs)%d3(JL,JK,3) + ZOMGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
!VH          &        frac_pom_sol_bb 
!VH 
!VH             emis_number(mode_acs)%d3(JL,JK,3) = &
!VH          &  emis_number(mode_acs)%d3(JL,JK,3) + ZOMGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
!VH          &        frac_pom_sol_bb * mass2numb_bb_sol
!VH         ENDDO
!VH       ENDDO
!VH     ENDIF ! LINJ
!VH   ENDIF ! LFIRE

END SUBROUTINE


SUBROUTINE BC_SRC_43R3() 

    ! mass to number conversion factors for the relevant modes
    numbscale_exp  = EXP(1.5*(LOG(sigma_lognormal(mode_aii)))**2)
    mass2numb_fact = 3./(4.*RPI*(numbscale_exp**3)*carbon_density)
    mass2numb_ff_insol  = mass2numb_fact/(rad_emi_ff_insol**3)
    mass2numb_ene_insol = mass2numb_fact/(rad_emi_ene_insol**3)
    mass2numb_ind_insol = mass2numb_fact/(rad_emi_ind_insol**3)
    mass2numb_tra_insol = mass2numb_fact/(rad_emi_tra_insol**3)
    mass2numb_shp_insol = mass2numb_fact/(rad_emi_shp_insol**3)
    mass2numb_air_insol = mass2numb_fact/(rad_emi_air_insol**3)
    mass2numb_bf_insol  = mass2numb_fact/(rad_emi_bf_insol**3)
    mass2numb_bb_insol  = mass2numb_fact/(rad_emi_bb_insol**3)

    numbscale_exp  = EXP(1.5*(LOG(sigma_lognormal(mode_ais)))**2)
    mass2numb_fact = 3./(4.*RPI*(numbscale_exp**3)*carbon_density)
    mass2numb_ff_sol  = mass2numb_fact/(rad_emi_ff_sol**3)
    mass2numb_ene_sol = mass2numb_fact/(rad_emi_ene_sol**3)
    mass2numb_ind_sol = mass2numb_fact/(rad_emi_ind_sol**3)
    mass2numb_tra_sol = mass2numb_fact/(rad_emi_tra_sol**3)
    mass2numb_shp_sol = mass2numb_fact/(rad_emi_shp_sol**3)
    mass2numb_air_sol = mass2numb_fact/(rad_emi_air_sol**3)

    numbscale_exp  = EXP(1.5*(LOG(sigma_lognormal(mode_acs)))**2)
    mass2numb_fact = 3./(4.*RPI*(numbscale_exp**3)*carbon_density)
    mass2numb_bf_sol = mass2numb_fact/(rad_emi_bf_sol**3)
    !mass2numb_bb_sol = mass2numb_fact/(rad_emi_bb_sol**3)
    mass2numb_nonbf_sol = mass2numb_ff_sol
    mass2numb_nonbf_insol = mass2numb_ff_insol
    
!    frac_bf(KIDIA:KFDIA)=1.0_JPRB
    ! calculate mass fraction related to solid biofuel
    where ( ZBCFF(KIDIA:KFDIA) > 1E-30_JPRB )
       frac_bf(KIDIA:KFDIA) = ZBCBF(KIDIA:KFDIA) / &
                                       ZBCFF(KIDIA:KFDIA)
    elsewhere
       frac_bf(KIDIA:KFDIA) = 0.0_JPRB
    endwhere

    ! for safety, prevent fractions larger than unity.
    where (frac_bf(KIDIA:KFDIA) > 1.0_JPRB )
       frac_bf(KIDIA:KFDIA) = 1.0_JPRB
    endwhere



   ! add to emis target arrays. 
   ! For now treat all sectors identical, and put all emissions in lowest model layer (KLEV)

   ! Fossil fuel categories..
   DO JL=KIDIA,KFDIA
        emis_mass  (mode_aii)%d3(JL,KLEV,1) = &
     &  emis_mass  (mode_aii)%d3(JL,KLEV,1) + ZBCFF(JL) * &
     &    ( (1.-frac_bf(JL)) * (1.-frac_bc_sol_ff) + &
     &          frac_bf(JL)  * (1.-frac_bc_sol_bf) ) 

        emis_number(mode_aii)%d3(JL,KLEV,1) = &
     &  emis_number(mode_aii)%d3(JL,KLEV,1) + ZBCFF(JL) * &
     &    ( (1.-frac_bf(JL)) * (1.-frac_bc_sol_ff) * mass2numb_nonbf_insol + &
     &          frac_bf(JL)  * (1.-frac_bc_sol_bf) * mass2numb_bf_insol )

        emis_mass  (mode_ais)%d3(JL,KLEV,2) = &
     &  emis_mass  (mode_ais)%d3(JL,KLEV,2) + ZBCFF(JL) * &
     &    ( (1.-frac_bf(JL)) * frac_bc_sol_ff )

        emis_number(mode_ais)%d3(JL,KLEV,2) = &
     &  emis_number(mode_ais)%d3(JL,KLEV,2) + ZBCFF(JL) * &
     &    ( (1.-frac_bf(JL)) * frac_bc_sol_ff * mass2numb_nonbf_sol )

        emis_mass  (mode_acs)%d3(JL,KLEV,2) = &
     &  emis_mass  (mode_acs)%d3(JL,KLEV,2) + ZBCFF(JL) * &
     &    (        frac_bf(JL)  * frac_bc_sol_bf )

        emis_number(mode_acs)%d3(JL,KLEV,2) = &
     &  emis_number(mode_acs)%d3(JL,KLEV,2) + ZBCFF(JL) * &
     &    (        frac_bf(JL)  * frac_bc_sol_bf * mass2numb_bf_sol )
    ENDDO
    
    ! Biofuel categories ? (ZBCBF emissions - currently not treated..)
    

    ! biomass burning

!VH   IF (LFIRE) THEN
!VH     IF (LINJ) THEN
!VH       DO JL=KIDIA,KFDIA
!VH       ! Height of injection for biomass burning emissions : update emis_mass
!VH 
!VH         IF (PINJF(JL) > 200._JPRB .AND. PBLH(JL) > 1500._JPRB) THEN
!VH           IX=MINLOC( ABS( (PAPHI(JL,1:KLEV)-PAPHI(JL,KLEV))/RG - PINJF(JL)))
!VH           ILINJ1=IX(1)
!VH           ILINJ2=ILINJ1
!VH           ! calculate total deltap over injected levels
!VH           ZDELP=0.0_JPRB
!VH           DO JK = ILINJ1, ILINJ2
!VH              ZDELP = ZDELP + PDELP(JL,JK)
!VH           ENDDO
!VH           DO JK = ILINJ1, ILINJ2
!VH        
!VH            ! add to emis target arrays
!VH               emis_mass  (mode_aii)%d3(JL,JK,1) = &
!VH            &  emis_mass  (mode_aii)%d3(JL,JK,1) + ZBCGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
!VH            &    (1.-frac_bc_sol_bb)
!VH       
!VH               emis_number(mode_aii)%d3(JL,JK,1) = &
!VH            &  emis_number(mode_aii)%d3(JL,JK,1) + ZBCGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
!VH            &    (1.-frac_bc_sol_bb) * mass2numb_bb_insol 
!VH 
!VH               emis_mass  (mode_acs)%d3(JL,JK,2) = &
!VH            &  emis_mass  (mode_acs)%d3(JL,JK,2) + ZBCGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
!VH            &     frac_bc_sol_bb 
!VH 
!VH               emis_number(mode_acs)%d3(JL,JK,2) = &
!VH            &  emis_number(mode_acs)%d3(JL,JK,2) + ZBCGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
!VH            &    frac_bc_sol_bb * mass2numb_bb_sol       
!VH           ENDDO
!VH         ELSE
!VH           ZDELP=0.0_JPRB
!VH           DO JK = KLEV-3, KLEV-2
!VH                ZDELP = ZDELP + PDELP(JL,JK)
!VH           ENDDO
!VH           DO JK = KLEV-3, KLEV-2
!VH               ! add to emis target arrays
!VH               emis_mass  (mode_aii)%d3(JL,JK,1) = &
!VH            &  emis_mass  (mode_aii)%d3(JL,JK,1) + ZBCGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
!VH            &    (1.-frac_bc_sol_bb)
!VH 
!VH               emis_number(mode_aii)%d3(JL,JK,1) = &
!VH            &  emis_number(mode_aii)%d3(JL,JK,1) + ZBCGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
!VH            &    (1.-frac_bc_sol_bb) * mass2numb_bb_insol 
!VH 
!VH               emis_mass  (mode_acs)%d3(JL,JK,2) = &
!VH            &  emis_mass  (mode_acs)%d3(JL,JK,2) + ZBCGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
!VH            &    frac_bc_sol_bb 
!VH 
!VH               emis_number(mode_acs)%d3(JL,JK,2) = &
!VH            &  emis_number(mode_acs)%d3(JL,JK,2) + ZBCGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
!VH            &    frac_bc_sol_bb * mass2numb_bb_sol
!VH           ENDDO
!VH         ENDIF
!VH       ENDDO
!VH     ELSE ! LINJ=false: always injection at lowest levels
!VH       DO JL=KIDIA,KFDIA
!VH         ZDELP=0.0_JPRB
!VH         DO JK = KLEV-2, KLEV-1
!VH           ZDELP = ZDELP + PDELP(JL,JK)
!VH         ENDDO
!VH         DO JK = KLEV-2, KLEV
!VH             ! add to emis target arrays
!VH             emis_mass  (mode_aii)%d3(JL,JK,1) = &
!VH          &  emis_mass  (mode_aii)%d3(JL,JK,1) + ZBCGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
!VH          &        (1.-frac_bc_sol_bb)
!VH 
!VH             emis_number(mode_aii)%d3(JL,JK,1) = &
!VH          &  emis_number(mode_aii)%d3(JL,JK,1) + ZBCGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
!VH          &        (1.-frac_bc_sol_bb) * mass2numb_bb_insol 
!VH 
!VH             emis_mass  (mode_acs)%d3(JL,JK,2) = &
!VH          &  emis_mass  (mode_acs)%d3(JL,JK,2) + ZBCGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
!VH          &        frac_bc_sol_bb 
!VH 
!VH             emis_number(mode_acs)%d3(JL,JK,2) = &
!VH          &  emis_number(mode_acs)%d3(JL,JK,2) + ZBCGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
!VH          &        frac_bc_sol_bb * mass2numb_bb_sol
!VH         ENDDO
!VH       ENDDO
!VH     ENDIF ! LINJ
!VH   ENDIF ! LFIRE

END SUBROUTINE

! SUBROUTINE ADD_TOFLUX()
!   emit(KIDIA:KFDIA,:) = 0.0
!   DO IMODE=1,NMOD                           ! loop in each mode 
!     DO INMODE=0,MODE_NM_SED(IMODE)          ! loop in aerosols species per mode 
!        JN = MODE_TRACERS_SED(INMODE,IMODE)  ! retrieve indentifier of each specie
!        if (JN==ino3_a.or.JN==inh4.or.JN==imsa) then
!           emit(KIDIA:KFDIA,:) = 0.0
!        else if(inmode == 0) then            ! 1st specie of mode is number part. in mode 
!           emit(KIDIA:KFDIA,:) = 0.0
!           do ii=1,mode_nm(IMODE)            ! add up all number emissions in the mode 'imode'...
!              emit(KIDIA:KFDIA,:) = emit(KIDIA:KFDIA,:) + emis_number(IMODE)%d3(KIDIA:KFDIA,:,ii)
!           enddo
!        else                                 ! this is a 'mass' emission with index nmode
!           emit(KIDIA:KFDIA,:) = emis_mass(IMODE)%d3(KIDIA:KFDIA,:,inmode)
!        endif
!        ! Change units from kg/m2/sec to kg/kg/sec and update tendency..
!        DO JL=KIDIA,KFDIA
!          ! Should limit to troposphere?! (for now sfc only)
!           !JK=91
!           !if (JN==iduai)then
!           !   write(2020,*)jk,emit(jl,jk)
!           !end if
!           !if (JN==iaii_n)then
!           !   write(2020,*)jk,emit(jl,jk)
!           !end if
!           !write(2929,*)JN,KAERO(JN)
!           PEMIDIAG(JL,KAERO(JN))=PEMIDIAG(JL,KAERO(JN))+sum(emit(JL,:))
!          DO JK=1,KLEV
!             !PCFLX(JL,KAERO(JN))=PCFLX(JL,KAERO(JN))+emit(JL,JK)
            
!             PTENC(JL,JK,KAERO(JN))=PTENC(JL,JK,KAERO(JN))+emit(JL,JK) * RG /PDELP(JL,JK)
!          ENDDO
         
!        ENDDO
!     ENDDO
!  ENDDO
! END SUBROUTINE ADD_TOFLUX

END SUBROUTINE TM5M7_SRC

