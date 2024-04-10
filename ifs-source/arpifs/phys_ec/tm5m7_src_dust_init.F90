SUBROUTINE TM5M7_SRC_DUST_INIT

!*** * TM5M7_SRC_DUST_INIT* - SOURCE TERMS FOR DUST AEROSOLS
!
! Online dust emissions based on Tegen/Vignati/Strunk
! ---------------------------------------------------
!
! Please read the section above for background information about the underlying 
! approach. An improved and modified online implementation has been accomplished 
! from which. It can be activated by setting 
!
!    input.emis.dust : ONLINE
!
! in the rc-file. An additional netcdf file is needed for some input parameters. 
! The path to which needs to be defined in the key
!
!    input.emis.dust.dir : /ms_perm/TM/TM5/emissions/other/Dust_online/onlinedust.nc
! 
! For every time step there will be particles emitted, scaled to monthly 
! amounts (both mass and numbers) in order to keep compliance with assumptions
! about the aerosol emissions in sedimentation.F90. 
!
!
!**   INTERFACE.
!     ----------
!          *TM5M7_SRC_DUST_INIT* IS CALLED FROM *TM5M7_INIT*.

!     AUTHOR.
!     -------
!
!     T. van Noije et al. (?)
!
!     SOURCE.
!     -------
!
!     MODIFICATIONS.
!     --------------
!
!      Nov 2011 - Achim Strunk - v0
!
!          xxxx - ??
!
!      Sep 2020 - V. Huijnen: first (partial) introduction into OpenIFS
!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST, ONLY : RPI, RG
USE TM5M7_DATA, ONLY : DDUST
USE TM5M7_EMIS_DATA, ONLY : nmode, nats,nclass, solspe, UTH, SREL, SRELV, SU_SRELV, &
    & DMIN, DMAX,DSTEP, A_RNOLDS,B_RNOLDS,X_RNOLDS, ROA, D_THRSLD

IMPLICIT NONE

!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

!       ... None ...


!*       0.2   LOCAL VARIABLES
!              ---------------

INTEGER(KIND=JPIM) :: NN, ND, NS, KK, NM, NSI, NP
REAL(KIND=JPRB)    :: BB, CCC, DDD, EE, FF, DP, XK, STOTAL,STOTALV
REAL(KIND=JPRB)    :: SU, SUV, SU_LOC, SU_LOCV, XL, XM, XN, XNV

!REAL, DIMENSION(:), ALLOCATABLE :: su_class, su_classv, utest
REAL(KIND=JPRB) :: su_class(nclass), su_classv(nclass), utest(nats)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TM5M7_SRC_DUST_INIT',0,ZHOOK_HANDLE)


! global fields, which have to be available throughout the whole


! only needed within "initial"
!ALLOCATE( su_class ( nclass ) )
!ALLOCATE( su_classv( nclass ) )
!ALLOCATE( utest    ( nats   ) )


! ---------------------------------------------
! read input file ... CHECK HOW TO HANDLE THIS ... 
! ---------------------------------------------



!---------------------------------------------------------------------------------------
!	 initializations
!---------------------------------------------------------------------------------------
uth	 = 0. 
srel	 = 0.	       ! fraction of the grid area correspondent to each soil population
srelV	 = 0.	       ! fraction of volume
su_srelV = 0.
utest	 = 0.

!---------------------------------------------------------------------------------------
!	Uth calculation 
!	Threshold friction velocity dependent on the particle diameter
!	following Eqs. (3-5) in MB95.
!---------------------------------------------------------------------------------------
nn = 0
dp = Dmin
DO WHILE( dp <= Dmax + 1E-05 )
   nn  = nn + 1
   BB  = a_rnolds * (dp ** x_rnolds) + b_rnolds
   XK  = SQRT(ddust * RG *100. * dp / roa)	   ! grav should be in cm s-2
   CCC = SQRT(1. + d_thrsld /(dp ** 2.5))
   IF( BB < 10. ) THEN
      DDD = SQRT(1.928 * (BB ** 0.092) - 1.)
      Uth(nn) = 0.129 * XK * CCC / DDD
   ELSE
      EE = -0.0617 * (BB - 10.) 
      FF = 1. -0.0858 * EXP(EE)
      Uth(nn) = 0.12 * XK * CCC * FF	   
   END IF
   dp = dp * EXP(Dstep)  
END DO


!---------------------------------------------------------------------------------------
!	surface calculation - calculation of the soil size distribution
!	Through all soil particle diameter the calculation of the relative contribution
!	in surface and volume of the soil population independently of the grid 
!---------------------------------------------------------------------------------------
DO ns = 1, nats ! soil types

   Stotal    = 0.
   StotalV   = 0.
   su_class  = 0. 
   su_classV = 0. 

   kk = 0
   dp = Dmin
   DO WHILE( dp <= Dmax + 1.0E-5 )	  
      kk  = kk + 1
      su  = 0.
      suV = 0.
      DO nm = 1, Nmode  	  ! particle size populations in soils
	 nd  = ((nm - 1) *3 ) + 1 ! index to mass median diameter
	 nsi = nd + 1		  ! index to standard deviation
	 np  = nd + 2		  ! index to relative contribution
	 !
	 !   based on soil type and contribution of population of the soil type the soil size
	 !   distribution population is calculated
	 !

	 !>>> TvN
	 ! Bug in the original code: nd should be np
	 ! Since solspe(ns,nd) is never zero
	 ! and the final result is proportional to solspe(ns,np),
	 ! the bug has no impact on the results. 
	 !IF (solspe(ns,nd).EQ.0.) THEN 	   
         IF (solspe(ns,np).EQ.0. .or. solspe(ns,nsi).EQ.0..or. solspe(ns,nd).EQ.0.) THEN
	 !<<< TvN
	    su_loc = 0.
	    su_locV=0.
	 ELSE
	    xk      = solspe(ns,np)/(SQRT(2.* RPI)*LOG(solspe(ns,nsi)))
            xl      = ( (LOG(dp) - LOG( solspe(ns,nd ) ))**2 ) / &
		 (2.*(LOG( solspe(ns,nsi) ))**2 )
	    xm      =  xk * EXP(-xl)	     ! value of the lognormal mass size distribution
					     ! dM/dln(dp) in Eq. (29) in MB95
					     ! (Aerosol Sci. Technol., 1994)
	    xn      =  ddust*(2./3.)*(dp/2.) ! surface
					     ! cf. the denominator in Eq. (30) in MB95
					     ! The factor 2 difference is irrelevant,
					     ! since only relative contributions are used.
	    xnV     =  1. !volume
	    su_loc  = (xm*Dstep/xn)	     ! Eq. (30) in MB95
	    su_locV = (xm*Dstep/xnV)	    
	 END IF !
	 su  = su  + su_loc
	 suV = suV + su_locV
      END DO !Nmode

      su_class(kk)   = su
      su_classV(kk)  = suV
      Stotal	     = Stotal + su
      StotalV	     = StotalV + suV
      dp	     = dp * EXP(Dstep)
   END DO !dp

   DO nn = 1,Nclass
      IF (Stotal.EQ.0.)THEN
	 srel (ns,nn) = 0.
	 srelV(ns,nn) = 0.
      ELSE
	 srel	 (ns,nn) = su_class(nn)/Stotal
	 srelV   (ns,nn) = su_classV(nn)/StotalV
	 utest   (ns   ) = utest(ns)+srelV(ns,nn)
	 su_srelV(ns,nn) = utest(ns) 
      END IF
   END DO !j=1,nclass
END DO !ns (soil type)


!DEALLOCATE( su_class, su_classV, utest )


IF (LHOOK) CALL DR_HOOK('TM5M7_SRC_DUST_INIT',1,ZHOOK_HANDLE)
END SUBROUTINE TM5M7_SRC_DUST_INIT

