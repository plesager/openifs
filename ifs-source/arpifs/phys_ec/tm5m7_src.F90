SUBROUTINE TM5M7_SRC &
 &( YDGEOMETRY, YDMODEL,  KIDIA, KFDIA, KLON , KTDIA, KLEV, KTILES, KSTART, KSTEP ,KSTGLO,&
 &  KSW  , KTRAC, KAERO,&
 &  PALB , PALBD, PAPHI ,&
 &  PAERDEP, PAERLTS, PAERSCC, PAERGUST, PALTH ,&
 &  PBCBF, PBCFF, PBCGF, POMBF, POMFF, POMGF,&
 &  PAPH , PAP  , PCI  , PCLAKE, PINJF, PBLH, PDELP, PGELAM, PGELAT, PGEMU , PFRTI , PHSDFOR,&
 &  PLSM , PSST , PQ   , PRHO , PSNS  , PT    , PTL   , PTSPHY, PZ0M, KCHEM,&
 &  PWIND, PWS1 ,PSOIL_TYPE, &
 &  PCVL, PCVH, KTVL, KTVH, &
 &  PLDAY,  PAERFLX, PCFLX , PCEN  , PTENC, PEMIDIAG, PSO2SRC,PSO4SRC,PSOA,PSOACO)

!*** * TM5M7_SRC* - SOURCE TERMS FOR TM5M7 AEROSOL SCHEME

!**   INTERFACE.
!     ----------
!          *TM5M7_SRC* IS CALLED FROM *TM5M7_PHY2*


!     AUTHOR.
!     -------
!        Vincent Huijnen  *KNMI*
!        ORIGINAL : 2020-08-25

!     MODIFICATIONS.
!     --------------
!
!
!-----------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE TYPE_MODEL   , ONLY : MODEL
USE PARKIND1  ,ONLY : JPIM, JPRB
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK
USE YOMLUN    ,ONLY : NULOUT, NULERR
!USE YOM_YGFL  ,ONLY : YGFL
USE YOMCST    ,ONLY : RA, RPI, RDAY, RG
USE YOMRIP0   ,ONLY : NINDAT, NSSSSS
!USE YOMRIP    ,ONLY : YRRIP
!USE YOEPHY    ,ONLY : YREPHY
!USE YOEAERATM ,ONLY : YREAERATM
!USE YOEAERMAP ,ONLY : YREAERMAP
!USE YOEAERSRC ,ONLY : YREAERSRC
!USE YOEAERVOL ,ONLY : YREAERVOL
!USE YOEAERSNK, ONLY : YREAERSNK
!USE YOMCOMPO,  ONLY : YRCOMPO  
USE TM5M7_DATA, ONLY :  NMOD, MODE_NM, MODE_NM_SED, MODE_TRACERS_SED, &
  & xmc, sigma_lognormal, pom_density, carbon_density, &
  & mode_aii, mode_ais, mode_acs, mode_aci,iduai,iaii_n,&
  & INO3_A, INH4,IMSA
USE TM5M7_EMIS_DATA, ONLY : MODAL_EMISSIONS, &
  & rad_emi_ff_insol,  rad_emi_ene_insol,rad_emi_ind_insol, &
  & rad_emi_tra_insol, rad_emi_shp_insol,rad_emi_air_insol, &
  & rad_emi_bf_insol,  rad_emi_bb_insol,&
  & rad_emi_ff_sol,    rad_emi_ene_sol,rad_emi_ind_sol, &
  & rad_emi_tra_sol,   rad_emi_shp_sol,rad_emi_air_sol, &
  & rad_emi_bf_sol,    rad_emi_bb_sol, &
  & frac_pom_sol_bf,   frac_pom_sol_bb, frac_pom_sol_ff, &
  & frac_bc_sol_bf,    frac_bc_sol_bb,  frac_bc_sol_ff
  USE OIFS_TO_HAM, ONLY: ind_oifs_ham!% ind_gas_OIFS

IMPLICIT NONE

!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(MODEL)       ,INTENT(IN)    :: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON, KIDIA, KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV, KTDIA, KSTGLO
INTEGER(KIND=JPIM),INTENT(IN)    :: KTILES
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTEP, KSTART
INTEGER(KIND=JPIM),INTENT(IN)    :: KSW
INTEGER(KIND=JPIM),INTENT(IN)    :: KTRAC
INTEGER(KIND=JPIM),INTENT(IN)    :: KAERO(YDMODEL%YRML_GCONF%YGFL%NAERO)

REAL(KIND=JPRB)   ,INTENT(IN)    :: PALB(KLON), PALBD(KLON,KSW)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHI(KLON,0:KLEV), PALTH(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAERDEP(KLON), PAERLTS(KLON), PAERSCC(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBCBF(KLON), PBCFF(KLON), PBCGF(KLON), POMBF(KLON), POMFF(KLON), POMGF(KLON)
REAL(KIND=JPRB),   INTENT(IN)    :: PAERGUST(KLON), PHSDFOR(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV), PAPH(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGELAM(KLON), PGELAT(KLON), PGEMU(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCI(KLON), PCLAKE(KLON), PLSM(KLON), PSST(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PINJF(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBLH(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDELP(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ(KLON,KLEV), PRHO(KLON,KLEV), PSNS(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTL(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWIND(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWS1(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSOIL_TYPE(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZ0M(KLON)
INTEGER(KIND=JPIM),INTENT(IN)    :: KCHEM(YDMODEL%YRML_GCONF%YGFL%NCHEM)

REAL(KIND=JPRB)   ,INTENT(INOUT) :: PAERFLX(KLON,12,9)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCFLX(KLON,KTRAC)

REAL(KIND=JPRB)  , INTENT(IN)     :: PCVL(KLON), PCVH(KLON) ! Low/High vegetation cover
INTEGER(KIND=JPIM), INTENT(IN)    :: KTVL(KLON), KTVH(KLON) ! Low/High vegetation type

REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCEN(KLON,KLEV,KTRAC)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENC(KLON,KLEV,KTRAC)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLDAY(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEMIDIAG(KLON,YDMODEL%YRML_GCONF%YGFL%NACTAERO)
REAL(KIND=JPRB),INTENT(INOUT)    :: PSO4SRC(KLON,KLEV),PSO2SRC(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSOACO(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSOA(KLON)


!*       0.5   LOCAL VARIABLES
!              ---------------

INTEGER(KIND=JPIM) :: JAER, JK, JL, IMODE, INMODE, JN, II, JGAS
INTEGER(KIND=JPIM) :: IGLGLO, IHTST

! TM5-M7 data

! Arrays to collect emissions
TYPE(MODAL_EMISSIONS), DIMENSION(NMOD), TARGET :: EMIS_MASS
TYPE(MODAL_EMISSIONS), DIMENSION(NMOD), TARGET :: EMIS_NUMBER


REAL(KIND=JPRB) :: ZFAERO(KLON,YDMODEL%YRML_GCONF%YGFL%NACTAERO)     , ZAEROCLIS(KLON,KLEV,2) 
REAL(KIND=JPRB) :: ZAEROK(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NACTAERO), ZTAEROK(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NACTAERO)

REAL(KIND=JPRB) :: ZGLAT(KLON), ZGLON(KLON)
REAL(KIND=JPRB) :: ZHDD, ZHSS
REAL(KIND=JPRB) :: ZDETAH(KLON,KLEV), ZETA(KLON,KLEV) , ZETAH(KLON,0:KLEV)

!-- various sources
REAL(KIND=JPRB) :: ZLOCALTIM   , ZDIURN(KLON)
REAL(KIND=JPRB) :: ZBCGF(KLON)  
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


REAL(KIND=JPRB) ::  numbscale_exp, mass2numb_fact, &
 &  mass2numb_ff_sol,     mass2numb_ff_insol,  mass2numb_ene_sol,     mass2numb_ene_insol, &
 &  mass2numb_ind_sol,     mass2numb_ind_insol, mass2numb_tra_sol,     mass2numb_tra_insol, &
 &  mass2numb_shp_sol,     mass2numb_shp_insol, &
 &  mass2numb_air_sol,     mass2numb_air_insol, mass2numb_bf_sol,     mass2numb_bf_insol, &
 &  mass2numb_bb_sol,     mass2numb_bb_insol,  mass2numb_nonbf_sol, mass2numb_nonbf_insol, &
 &  oc2pom

REAL(KIND=JPRB)    :: ZSOA(KLON)

REAL(KIND=JPRB) :: FRAC_BF(KLON), EMIT(KLON,KLEV) 

INTEGER(KIND=JPIM) :: ISSO2, ISSO4
!-- Injection height for biomass burning emissions
INTEGER(KIND=JPIM) :: ILINJ1, ILINJ2, IX(1)
REAL(KIND=JPRB)    :: ZDELP

REAL(KIND=JPRB)    :: ZAERMAP(KLON,5)
#ifdef __PGI
REAL(KIND=JPRB) :: ERF
#else
INTRINSIC ERF
#endif

REAL(KIND=JPRB)    :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "updcal.intfb.h"
!#include "fcttim.func.h"

#include "surf_inq.h"

#include "tm5m7_src_ss.intfb.h"
#include "tm5m7_src_dust.intfb.h"
!#include "satur.intfb.h"
!#include "aer_volce.intfb.h"
!#include "aer_stratcl.intfb.h"

IF (LHOOK) CALL DR_HOOK('TM5M7_SRC',0,ZHOOK_HANDLE)


!-----------------------------------------------------------------------
ASSOCIATE(YDVAB=>YDGEOMETRY%YRVAB,YDVETA=>YDGEOMETRY%YRVETA,YDVFE=>YDGEOMETRY%YRVFE,&
  & YDCSGLEG=>YDGEOMETRY%YRCSGLEG,YDVSPLIP=>YDGEOMETRY%YRVSPLIP,YDVSLETA=>YDGEOMETRY%YRVSLETA,&
  & YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY, &
  & YDEAERMAP=>YDMODEL%YRML_PHY_AER%YREAERMAP, &
  & YGFL=>YDMODEL%YRML_GCONF%YGFL, &
  & YDCOMPO=>YDMODEL%YRML_CHEM%YRCOMPO, &
  & YDEAERSRC=>YDMODEL%YRML_PHY_AER%YREAERSRC, &
  & YDEAERATM=>YDMODEL%YRML_PHY_RAD%YREAERATM, &
  & YDRIP=>YDMODEL%YRML_GCONF%YRRIP)



ASSOCIATE(LAERODIU=>YDCOMPO%LAERODIU, YAERO=>YGFL%YAERO, LFIRE=>YDCOMPO%LFIRE, LINJ=>YDCOMPO%LINJ, &
 & NACTAERO=>YGFL%NACTAERO, NAERO=>YGFL%NAERO, &
 & NDGLG=>YDGEOMETRY%YRDIM%NDGLG, RHGMT=>YDRIP%RHGMT, &
 & RSTATI => YDRIP%RSTATI, RSIDECA=>YDEAERSRC%RSIDECA, &
 & NAERWND=>YDEAERSRC%NAERWND, &
 & RSIVSRA=>YDEAERSRC%RSIVSRA, &
 & RCODECA=>YDEAERSRC%RCODECA, RCOVSRA=>YDEAERSRC%RCOVSRA, &
 & NLOENG=>YDGEOMETRY%YRGEM%NLOENG, &
 & NGLOBALAT=>YDGEOMETRY%YRMP%NGLOBALAT, &
 & YSURF=>YDEPHY%YSURF,LVDFTRAC=>YDEPHY%LVDFTRAC, &
 & NCHEM=>YGFL%NCHEM, YCHEM=>YGFL%YCHEM, LAERSOA_CHEM=>YDEAERATM%LAERSOA_CHEM, &
 & LAERCHEM=>YGFL%LAERCHEM)

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
CALL UPDCAL (ID0, IM0, IY0, INC,  IDD, IMM, IYY, IMON, -1)
IMDATE=IYY*10000+IMM*100+IDD

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
   & ,0.0_JPRB)
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
ZGRDLAT= RPI / NDGLG                          ! distance in radians between latitude lines
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
IF (.NOT.LAERODIU) THEN
  ZDIURN(KIDIA:KFDIA)=1.0_JPRB
ENDIF



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

PEMIDIAG(KIDIA:KFDIA,         1:NACTAERO)= 0.0_JPRB
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

!*       1.0   SEA SALT
!              --------

!- Simplistic lifting from surface based on 10-m wind and land-sea mask
! (currently not used!)
ZHSS=8434._JPRB/1000._JPRB 


  CALL TM5M7_SRC_SS(YDEAERSRC, KIDIA, KFDIA, KLON, KLEV, &
      & PCI, PCLAKE, PLSM, PSST, ZWNDSS, &
      & emis_mass, emis_number )



!-----------------------------------------------------------------------

!*       2.0   DESERT DUST
!              -----------

!- Simplistic lifting from surface based on 10-m wind and surface albedo
ZHDD=MAX(1.0_JPRB,8434._JPRB/1000._JPRB)

PAERFLX(KIDIA:KFDIA,1:12,1:9)=0._JPRB
ZAERMAP(KIDIA:KFDIA,1:5)=0._JPRB
CALL TM5M7_SRC_DUST(YDMODEL, KIDIA, KFDIA, KLON, KLEV, KTILES, KSW,&
      & PLSM, ZWNDDU, PSNS, PZ0M, &
      & PAP(:,KLEV), PTL,  PSOIL_TYPE, &
      & PFRTI, PCVL, PCVH, KTVL, KTVH, &
      & emis_mass, emis_number ,PAERFLX,ZGLON,ZGLAT&
      ,ZRWPWP,ZRWSAT,ZAERMAP,PALB,PALBD,PWS1,PHSDFOR)
!write(2345,*) 'test',ptsphy, emis_mass(mode_aci)%d3(KIDIA:KFDIA,91,1),emis_number(mode_aci)%d3(KIDIA:KFDIA,91,1),PAERFLX(:,1,:)
!write(2346,*) 'test',ptsphy, emis_mass(mode_aii)%d3(KIDIA:KFDIA,91,1),emis_number(mode_aii)%d3(KIDIA:KFDIA,91,1),PAERFLX(:,1,:)
!DO JK=1,KLEV
!   DO JL=KIDIA,KFDIA
!      ZAEROUT5(JL,JK)=emis_mass(mode_aci)%d3(Jl,JK,1)
!   END DO
!END DO
!-----------------------------------------------------------------------

!*       3.0   PARTICULATE ORGANIC MATTER
!              --------------------------------------------------------------


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
    mass2numb_bb_sol = mass2numb_fact/(rad_emi_bb_sol**3) 
    mass2numb_nonbf_sol = mass2numb_ff_sol
    mass2numb_nonbf_insol = mass2numb_ff_insol

    frac_bf(KIDIA:KFDIA)=1.0_JPRB
    ! calculate mass fraction related to solid biofuel
    where ( POMFF(KIDIA:KFDIA) > 1E-30_JPRB )
       frac_bf(KIDIA:KFDIA) = POMBF(KIDIA:KFDIA) / &
                                       POMFF(KIDIA:KFDIA)
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
     &  emis_mass  (mode_aii)%d3(JL,KLEV,2) + POMFF(JL) * &
     &    ( (1.-frac_bf(JL)) * (1.-frac_pom_sol_ff) + &
     &          frac_bf(JL)  * (1.-frac_pom_sol_bf) ) 

        emis_number(mode_aii)%d3(JL,KLEV,2) = &
     &  emis_number(mode_aii)%d3(JL,KLEV,2) + POMFF(JL) * &
     &    ( (1.-frac_bf(JL)) * (1.-frac_pom_sol_ff) * mass2numb_nonbf_insol + &
     &          frac_bf(JL)  * (1.-frac_pom_sol_bf) * mass2numb_bf_insol )

        emis_mass  (mode_ais)%d3(JL,KLEV,3) = &
     &  emis_mass  (mode_ais)%d3(JL,KLEV,3) + POMFF(JL) * &
     &    ( (1.-frac_bf(JL)) * frac_pom_sol_ff )

        emis_number(mode_ais)%d3(JL,KLEV,3) = &
     &  emis_number(mode_ais)%d3(JL,KLEV,3) + POMFF(JL) * &
     &    ( (1.-frac_bf(JL)) * frac_pom_sol_ff * mass2numb_nonbf_sol )

        emis_mass  (mode_acs)%d3(JL,KLEV,3) = &
     &  emis_mass  (mode_acs)%d3(JL,KLEV,3) + POMFF(JL) * &
     &    (        frac_bf(JL)  * frac_pom_sol_bf )

        emis_number(mode_acs)%d3(JL,KLEV,3) = &
     &  emis_number(mode_acs)%d3(JL,KLEV,3) + POMFF(JL) * &
     &    (        frac_bf(JL)  * frac_pom_sol_bf * mass2numb_bf_sol )
    ENDDO
    
    ! Biofuel categories ? (POMBF emissions)
    
!!$    IF (.not. LAERCHEM)THEN
!!$       ! SOA from CO
    DO JL=KIDIA,KFDIA
       ZSOA(JL)=0._JPRB
       IF (LAERSOA_CHEM) THEN
          ZSOA(JL)=MAX(PSOACO(JL),PSOA(JL))
       ELSE
          ZSOA(JL)=PSOACO(JL)
       ENDIF
       ZOMSOURC=ZOMSOURC+ZSOA(JL)
    END DO
    ! These do not apply for M7
!!$       PCFLX(JL,KAERO(INBAER+1))= -ZOMSOURC * ROMPHIL
!!$       PCFLX(JL,KAERO(INBAER+2))= -ZOMSOURC * ROMPHOB 
!!$    END IF
    ! biomass burning

  IF (LFIRE) THEN
    IF (LINJ) THEN
      DO JL=KIDIA,KFDIA
      ! Height of injection for biomass burning emissions : update emis_mass
        IF (PINJF(JL) > 200._JPRB .AND. PBLH(JL) > 1500._JPRB) THEN
          IX=MINLOC( ABS( (PAPHI(JL,1:KLEV)-PAPHI(JL,KLEV))/RG - PINJF(JL)))
          ILINJ1=IX(1)
          ILINJ2=ILINJ1
          ! calculate total deltap over injected levels
          ZDELP=0.0_JPRB
          DO JK = ILINJ1, ILINJ2
             ZDELP = ZDELP + PDELP(JL,JK)
          ENDDO
          DO JK = ILINJ1, ILINJ2
       
           ! add to emis target arrays
              emis_mass  (mode_aii)%d3(JL,JK,2) = &
           &  emis_mass  (mode_aii)%d3(JL,JK,2) + POMGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
           &    (1.-frac_pom_sol_bb)
      
              emis_number(mode_aii)%d3(JL,JK,2) = &
           &  emis_number(mode_aii)%d3(JL,JK,2) + POMGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
           &    (1.-frac_pom_sol_bb) * mass2numb_bb_insol 

              emis_mass  (mode_acs)%d3(JL,JK,3) = &
           &  emis_mass  (mode_acs)%d3(JL,JK,3) + POMGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
           &     frac_pom_sol_bb 

              emis_number(mode_acs)%d3(JL,JK,3) = &
           &  emis_number(mode_acs)%d3(JL,JK,3) + POMGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
           &    frac_pom_sol_bb * mass2numb_bb_sol       
          ENDDO
        ELSE
          ZDELP=0.0_JPRB
          DO JK = KLEV-3, KLEV-2
               ZDELP = ZDELP + PDELP(JL,JK)
          ENDDO
          DO JK = KLEV-3, KLEV-2
              ! add to emis target arrays
              emis_mass  (mode_aii)%d3(JL,JK,2) = &
           &  emis_mass  (mode_aii)%d3(JL,JK,2) + POMGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
           &    (1.-frac_pom_sol_bb)

              emis_number(mode_aii)%d3(JL,JK,2) = &
           &  emis_number(mode_aii)%d3(JL,JK,2) + POMGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
           &    (1.-frac_pom_sol_bb) * mass2numb_bb_insol 

              emis_mass  (mode_acs)%d3(JL,JK,3) = &
           &  emis_mass  (mode_acs)%d3(JL,JK,3) + POMGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
           &    frac_pom_sol_bb 

              emis_number(mode_acs)%d3(JL,JK,3) = &
           &  emis_number(mode_acs)%d3(JL,JK,3) + POMGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
           &    frac_pom_sol_bb * mass2numb_bb_sol
          ENDDO
        ENDIF
      ENDDO
    ELSE ! LINJ=false: always injection at lowest levels
      DO JL=KIDIA,KFDIA
        ZDELP=0.0_JPRB
        DO JK = KLEV-2, KLEV-1
          ZDELP = ZDELP + PDELP(JL,JK)
        ENDDO
        DO JK = KLEV-2, KLEV
            ! add to emis target arrays
            emis_mass  (mode_aii)%d3(JL,JK,2) = &
         &  emis_mass  (mode_aii)%d3(JL,JK,2) + POMGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
         &        (1.-frac_pom_sol_bb)

            emis_number(mode_aii)%d3(JL,JK,2) = &
         &  emis_number(mode_aii)%d3(JL,JK,2) + POMGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
         &        (1.-frac_pom_sol_bb) * mass2numb_bb_insol 

            emis_mass  (mode_acs)%d3(JL,JK,3) = &
         &  emis_mass  (mode_acs)%d3(JL,JK,3) + POMGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
         &        frac_pom_sol_bb 

            emis_number(mode_acs)%d3(JL,JK,3) = &
         &  emis_number(mode_acs)%d3(JL,JK,3) + POMGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
         &        frac_pom_sol_bb * mass2numb_bb_sol
        ENDDO
      ENDDO
    ENDIF ! LINJ
  ENDIF ! LFIRE


!-----------------------------------------------------------------------

!*       4.0   BLACK CARBON
!              ------------



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
    mass2numb_bb_sol = mass2numb_fact/(rad_emi_bb_sol**3)
    mass2numb_nonbf_sol = mass2numb_ff_sol
    mass2numb_nonbf_insol = mass2numb_ff_insol
    
!    frac_bf(KIDIA:KFDIA)=1.0_JPRB
    ! calculate mass fraction related to solid biofuel
    where ( PBCFF(KIDIA:KFDIA) > 1E-30_JPRB )
       frac_bf(KIDIA:KFDIA) = PBCBF(KIDIA:KFDIA) / &
                                       PBCFF(KIDIA:KFDIA)
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
     &  emis_mass  (mode_aii)%d3(JL,KLEV,1) + PBCFF(JL) * &
     &    ( (1.-frac_bf(JL)) * (1.-frac_bc_sol_ff) + &
     &          frac_bf(JL)  * (1.-frac_bc_sol_bf) ) 

        emis_number(mode_aii)%d3(JL,KLEV,1) = &
     &  emis_number(mode_aii)%d3(JL,KLEV,1) + PBCFF(JL) * &
     &    ( (1.-frac_bf(JL)) * (1.-frac_bc_sol_ff) * mass2numb_nonbf_insol + &
     &          frac_bf(JL)  * (1.-frac_bc_sol_bf) * mass2numb_bf_insol )

        emis_mass  (mode_ais)%d3(JL,KLEV,2) = &
     &  emis_mass  (mode_ais)%d3(JL,KLEV,2) + PBCFF(JL) * &
     &    ( (1.-frac_bf(JL)) * frac_bc_sol_ff )

        emis_number(mode_ais)%d3(JL,KLEV,2) = &
     &  emis_number(mode_ais)%d3(JL,KLEV,2) + PBCFF(JL) * &
     &    ( (1.-frac_bf(JL)) * frac_bc_sol_ff * mass2numb_nonbf_sol )

        emis_mass  (mode_acs)%d3(JL,KLEV,2) = &
     &  emis_mass  (mode_acs)%d3(JL,KLEV,2) + PBCFF(JL) * &
     &    (        frac_bf(JL)  * frac_bc_sol_bf )

        emis_number(mode_acs)%d3(JL,KLEV,2) = &
     &  emis_number(mode_acs)%d3(JL,KLEV,2) + PBCFF(JL) * &
     &    (        frac_bf(JL)  * frac_bc_sol_bf * mass2numb_bf_sol )
    ENDDO
    
    ! Biofuel categories ? (PBCBF emissions - currently not treated..)
    

    ! biomass burning

  IF (LFIRE) THEN
    IF (LINJ) THEN
      DO JL=KIDIA,KFDIA
      ! Height of injection for biomass burning emissions : update emis_mass

        IF (PINJF(JL) > 200._JPRB .AND. PBLH(JL) > 1500._JPRB) THEN
          IX=MINLOC( ABS( (PAPHI(JL,1:KLEV)-PAPHI(JL,KLEV))/RG - PINJF(JL)))
          ILINJ1=IX(1)
          ILINJ2=ILINJ1
          ! calculate total deltap over injected levels
          ZDELP=0.0_JPRB
          DO JK = ILINJ1, ILINJ2
             ZDELP = ZDELP + PDELP(JL,JK)
          ENDDO
          DO JK = ILINJ1, ILINJ2
       
           ! add to emis target arrays
              emis_mass  (mode_aii)%d3(JL,JK,1) = &
           &  emis_mass  (mode_aii)%d3(JL,JK,1) + PBCGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
           &    (1.-frac_bc_sol_bb)
      
              emis_number(mode_aii)%d3(JL,JK,1) = &
           &  emis_number(mode_aii)%d3(JL,JK,1) + PBCGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
           &    (1.-frac_bc_sol_bb) * mass2numb_bb_insol 

              emis_mass  (mode_acs)%d3(JL,JK,2) = &
           &  emis_mass  (mode_acs)%d3(JL,JK,2) + PBCGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
           &     frac_bc_sol_bb 

              emis_number(mode_acs)%d3(JL,JK,2) = &
           &  emis_number(mode_acs)%d3(JL,JK,2) + PBCGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
           &    frac_bc_sol_bb * mass2numb_bb_sol       
          ENDDO
        ELSE
          ZDELP=0.0_JPRB
          DO JK = KLEV-3, KLEV-2
               ZDELP = ZDELP + PDELP(JL,JK)
          ENDDO
          DO JK = KLEV-3, KLEV-2
              ! add to emis target arrays
              emis_mass  (mode_aii)%d3(JL,JK,1) = &
           &  emis_mass  (mode_aii)%d3(JL,JK,1) + PBCGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
           &    (1.-frac_bc_sol_bb)

              emis_number(mode_aii)%d3(JL,JK,1) = &
           &  emis_number(mode_aii)%d3(JL,JK,1) + PBCGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
           &    (1.-frac_bc_sol_bb) * mass2numb_bb_insol 

              emis_mass  (mode_acs)%d3(JL,JK,2) = &
           &  emis_mass  (mode_acs)%d3(JL,JK,2) + PBCGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
           &    frac_bc_sol_bb 

              emis_number(mode_acs)%d3(JL,JK,2) = &
           &  emis_number(mode_acs)%d3(JL,JK,2) + PBCGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
           &    frac_bc_sol_bb * mass2numb_bb_sol
          ENDDO
        ENDIF
      ENDDO
    ELSE ! LINJ=false: always injection at lowest levels
      DO JL=KIDIA,KFDIA
        ZDELP=0.0_JPRB
        DO JK = KLEV-2, KLEV-1
          ZDELP = ZDELP + PDELP(JL,JK)
        ENDDO
        DO JK = KLEV-2, KLEV
            ! add to emis target arrays
            emis_mass  (mode_aii)%d3(JL,JK,1) = &
         &  emis_mass  (mode_aii)%d3(JL,JK,1) + PBCGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
         &        (1.-frac_bc_sol_bb)

            emis_number(mode_aii)%d3(JL,JK,1) = &
         &  emis_number(mode_aii)%d3(JL,JK,1) + PBCGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
         &        (1.-frac_bc_sol_bb) * mass2numb_bb_insol 

            emis_mass  (mode_acs)%d3(JL,JK,2) = &
         &  emis_mass  (mode_acs)%d3(JL,JK,2) + PBCGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
         &        frac_bc_sol_bb 

            emis_number(mode_acs)%d3(JL,JK,2) = &
         &  emis_number(mode_acs)%d3(JL,JK,2) + PBCGF(JL)*ZDIURN(JL)*PDELP(JL,JK)/ZDELP * &
         &        frac_bc_sol_bb * mass2numb_bb_sol
        ENDDO
      ENDDO
    ENDIF ! LINJ
  ENDIF ! LFIRE




!-----------------------------------------------------------------------

!*       5.0   Convert emissions into tendencies: loop over tracers in mode
!              ------------

  DO IMODE=1,NMOD
    DO INMODE=0,MODE_NM_SED(IMODE)
       JN = MODE_TRACERS_SED(INMODE,IMODE)
       if (JN==ino3_a.or.JN==inh4.or.JN==imsa) then
          emit(KIDIA:KFDIA,:) = 0.0
       else if(inmode == 0) then
          emit(KIDIA:KFDIA,:) = 0.0
          do ii=1,mode_nm(IMODE)   ! add up all number emissions in the mode 'imode'...
             emit(KIDIA:KFDIA,:) = emit(KIDIA:KFDIA,:) + emis_number(IMODE)%d3(KIDIA:KFDIA,:,ii)
          enddo
       else ! this is a 'mass' emission with index nmode
          emit(KIDIA:KFDIA,:) = emis_mass(IMODE)%d3(KIDIA:KFDIA,:,inmode)
       endif
       ! Change units from kg/m2/sec to kg/kg/sec and update tendency..
       DO JL=KIDIA,KFDIA
         ! Should limit to troposphere?! (for now sfc only)
          !JK=91
          !if (JN==iduai)then
          !   write(2020,*)jk,emit(jl,jk)
          !end if
          !if (JN==iaii_n)then
          !   write(2020,*)jk,emit(jl,jk)
          !end if
          !write(2929,*)JN,KAERO(JN)
          PEMIDIAG(JL,KAERO(JN))=PEMIDIAG(JL,JN)+sum(emit(JL,:))
         DO JK=1,KLEV
            !PCFLX(JL,KAERO(JN))=PCFLX(JL,KAERO(JN))+emit(JL,JK)
            
            PTENC(JL,JK,KAERO(JN))=PTENC(JL,JK,KAERO(JN))+emit(JL,JK) * RG /PDELP(JL,JK)
         ENDDO
         
       ENDDO
    ENDDO
 ENDDO
! write(9292,*)PSO2SRC(JL,JK),PSO4SRC(JL,JK)
 IF (.not. LAERCHEM) THEN
    DO JL=KIDIA,KFDIA
       DO JK=1,KLEV
          DO JGAS=1,2
             IF (TRIM(YAERO(ind_oifs_ham%ind_gas_OIFS(JGAS))%CNAME)=='SO2') THEN
                ISSO2=ind_oifs_ham%ind_gas_OIFS(JGAS)
                                !write(9292,*)ISSO2,PSO2SRC(JL,JK)
                PTENC(JL,JK,KAERO(ISSO2))=PTENC(JL,JK,KAERO(ISSO2))+ PSO2SRC(JL,JK)
                !PCFLX(JL,KAERO(ISSO2))=PCFLX(JL,KAERO(ISSO2)) + PSO2SRC(JL,JK)
                PEMIDIAG(JL,KAERO(ISSO2))=PEMIDIAG(JL,KAERO(ISSO2))+ PSO2SRC(JL,JK)
             ELSE IF (TRIM(YAERO(ind_oifs_ham%ind_gas_OIFS(JGAS))%CNAME)=='SO4_gas') THEN
                ISSO4=ind_oifs_ham%ind_gas_OIFS(JGAS)
                !write(9292,*)ISSO4, PSO4SRC(JL,JK)
                PTENC(JL,JK,KAERO(ISSO4))=PTENC(JL,JK,KAERO(ISSO4))+ PSO4SRC(JL,JK)
                !PCFLX(JL,KAERO(ISSO4))=PCFLX(JL,KAERO(ISSO4)) + PSO4SRC(JL,JK)
                PEMIDIAG(JL,KAERO(ISSO4))=PEMIDIAG(JL,KAERO(ISSO4)) + PSO4SRC(JL,JK)
             END IF
          END DO
       END DO


       ! For add SOA from CO into ISVOC tracer
!!$       DO JGAS=1,NACTAERO
!!$          IF (TRIM(YAERO(JGAS)%CNAME)=='ISVOC') THEN
!!$             
!!$             PTENC(JL,JK,KAERO(JGAS))=PTENC(JL,JK,KAERO(JGAS))+ PSOACO(JL)
!!$             PEMIDIAG(JL,KAERO(JGAS))=PEMIDIAG(JL,KAERO(JGAS)) + PSOACO(JL)
!!$          END IF
!!$       END DO

    END DO
 END IF
!  PGFL(KIDIA:KFDIA,:,YAEROUT(5)%MP)=ZAEROUT5(KIDIA:KFDIA,:)

!-----------------------------------------------------------------------

!*       6.0   De-allocate arrays
!              ------------


DO IMODE=1,NMOD
  DEALLOCATE(EMIS_NUMBER(IMODE)%d3)
  DEALLOCATE(EMIS_MASS(IMODE)%d3)
ENDDO

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('TM5M7_SRC',1,ZHOOK_HANDLE)
END SUBROUTINE TM5M7_SRC

