SUBROUTINE TM5M7_OPTICS_INIT(NWAV,WDEP)

!*** * TM5M7_OPTICS_INIT* 
!
!
!-----------------------------------------------------------------------------
!                    TM5                                                     !
!-----------------------------------------------------------------------------
!BOP
!
! !MODULE:      OPTICS
!
! !DESCRIPTION: Optics module to calculate optical depth from m7 output,
!               based on the AOP_Package op Michael Kahnert.
!               
!               The optics code should be used in the following way (see
!               examples in photolysis.F90, ecearth_optics.F90, and
!               user_output_aerocom.F90): 
!               
!               1) prepare the optics calculation by calling 
!                  'OPTICS_INIT' --> lookuptables etc. 
!                  this routine reads the wavelengths, lookupable and calculates
!                  refr. indices at these wavelengths.
!               
!               2) allocate AOP arrays aop_out_ext/a/g with 
!                  (n_gridcells, n_wavelengths, n_split), with:
!                  'n_split' = 1 for (split == .false.) or 
!                  'n_split' = 11 for (split == .true. ), ie. 
!                    partial contributions by 
!                     Total, SO4, BC, OC, SS, DU, NO3, Water, Fine, Fine Dust, Fine SS
!                  additional fields are also provided for split==.true. in 
!                   aop_out_add, which has to have size 
!                   (n_gridcells, n_wavelengths, n_add) with nadd = 2 for 
!                    surface PM10 dry extinction and surface PM10 dry absorption
!               3) Compute AOP for current conditions by calling 'OPTICS_AOP_GET'
!               4) done: 'OPTICS_DONE'
!
!  IMPORTANT:   *) Skipped the ZOOM options! (temporary) 
!               *) OC is actually POM. Remember converting OC to POM 
!                  when sending it to optics_calculate_aop
!\\
!
!
!**   INTERFACE.
!     ----------
!          *TM5M7_OPTICS_INIT* IS CALLED FROM *TM5M7_INIT*.

!     AUTHOR.
!     -------
!
!    6 Feb 2011 - Achim Strunk - worked on DECOUPLING optics routines 
!				 from (optics)_output, in order to 
!				 (re)establish a flexible way of using it
!			       - remaining parts have been moved to
!				 (the new routine) user_output_optics
!
!   24 Jun 2011 - Achim Strunk - added NO3 explicitly in order to allow
!				 for a slightly better split of the
!				 optical properties of (SO4+NO3)
  
!     SOURCE.
!     -------
!
!     Taken from TM5-code
!
!     MODIFICATIONS.
!     --------------
!
!      Sep 2021 - V. Huijnen: first introduction into OpenIFS
!
!
!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMLUN    ,ONLY : NULOUT

USE YOMCST, ONLY : 
USE TM5M7_DATA, ONLY : TM5M7_DATADIR
USE TM5M7_EMIS_DATA, ONLY : 
USE TM5M7_OPTICS_DATA, ONLY: WAVELENDEP, N_RIR,N_RII,N_X, &
  & lkval, &     ! -log img part refr. index
  & kval, &      ! img part refr. index, 10^(-lkval)
  & n1r, &       ! real part refr. index
  & xs, cext_159, a_159, g_159, cext_200, a_200, g_200, &
  & opacdim,echamhamdim,segelsteindim,&
  & opac,echamham,segelstein


IMPLICIT NONE

!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

INTEGER(KIND=JPIM), INTENT(IN) :: NWAV
TYPE(WAVELENDEP), DIMENSION(NWAV), INTENT(INOUT) :: WDEP


!*       0.2   LOCAL VARIABLES
!              ---------------

  character(len=256) :: lookuptable, refractive_indices
  character(len=256) :: CL_FILE_NAME
  INTEGER(KIND=JPIM) :: JL


REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TM5M7_OPTICS_INIT',0,ZHOOK_HANDLE)

  allocate( opac      (7, opacdim      ) )
  allocate( echamham  (5, echamhamdim  ) )
  allocate( segelstein(3, segelsteindim) )

  ! identify filenames from input
  lookuptable = "lookup_table.nc"
  !CL_FILE_NAME="/perm/ms/nl/nm6/TM5_INPUT/photolysis/lookup_table.nc"
  CL_FILE_NAME = TRIM(TM5M7_DATADIR)//lookuptable


  CALL LOAD_TM5M7_OPTICS_DATA_1D(CL_FILE_NAME, 'mr', n_rir,n1r)
  CALL LOAD_TM5M7_OPTICS_DATA_1D(CL_FILE_NAME, 'mi', n_rii,kval)
  CALL LOAD_TM5M7_OPTICS_DATA_1D(CL_FILE_NAME, 'x' , n_x  ,xs)

  lkval = -log10(kval)

  CALL LOAD_TM5M7_OPTICS_DATA_3D(CL_FILE_NAME, 'ext_159', n_x,n_rir,n_rii,cext_159)
  CALL LOAD_TM5M7_OPTICS_DATA_3D(CL_FILE_NAME, 'ext_200', n_x,n_rir,n_rii,cext_200)
  CALL LOAD_TM5M7_OPTICS_DATA_3D(CL_FILE_NAME, 'a_159', n_x,n_rir,n_rii,a_159)
  CALL LOAD_TM5M7_OPTICS_DATA_3D(CL_FILE_NAME, 'a_200', n_x,n_rir,n_rii,a_200)
  CALL LOAD_TM5M7_OPTICS_DATA_3D(CL_FILE_NAME, 'g_159', n_x,n_rir,n_rii,g_159)
  CALL LOAD_TM5M7_OPTICS_DATA_3D(CL_FILE_NAME, 'g_200', n_x,n_rir,n_rii,g_200)


  refractive_indices = "refractive_indices_hdfstyle.nc"
  CL_FILE_NAME = TRIM(TM5M7_DATADIR)//refractive_indices

  CALL LOAD_TM5M7_OPTICS_DATA_2D(CL_FILE_NAME, 'Opac', 7,opacdim,opac)
  CALL LOAD_TM5M7_OPTICS_DATA_2D(CL_FILE_NAME, 'ECHAM-HAM', 5,echamhamdim,echamham)
  CALL LOAD_TM5M7_OPTICS_DATA_2D(CL_FILE_NAME, 'Segelstein', 3,segelsteindim,segelstein)

  call TM5M7_optics_wavelen_init( wdep,NWAV )

  write(NULOUT,*) 'Optical parameters:'
  do JL = 1, NWAV
      write(NULOUT,*) 'Wavelength :', wdep(JL)%wl 
      write(NULOUT,*) ' SO4 (real/img) :', wdep(JL)%n(1), wdep(JL)%k(1) 
      write(NULOUT,*) ' BC  (real/img) :', wdep(JL)%n(2), wdep(JL)%k(2) 
      write(NULOUT,*) ' OC  (real/img) :', wdep(JL)%n(3), wdep(JL)%k(3) 
      write(NULOUT,*) ' SOA (real/img) :', wdep(JL)%n(4), wdep(JL)%k(4) 
      write(NULOUT,*) ' SS  (real/img) :', wdep(JL)%n(5), wdep(JL)%k(5) 
      write(NULOUT,*) ' DU  (real/img) :', wdep(JL)%n(6), wdep(JL)%k(6) 
      write(NULOUT,*) ' H2O (real/img) :', wdep(JL)%n(7), wdep(JL)%k(7) 
  enddo
     

  if (allocated(opac)) deallocate( opac)
  if (allocated(echamham)) deallocate(echamham )
  if (allocated(segelstein )) deallocate(  segelstein )

  !--
  
  IF(LHOOK) CALL DR_HOOK('TM5_OPTICS_INIT',1,ZHOOK_HANDLE)
CONTAINS

  ! Load an individual field  (or receive it from the first
  ! MPI task)
  SUBROUTINE LOAD_TM5M7_OPTICS_DATA_3D(CD_FILE_NAME, CD_VAR_NAME, ILEN,JLEN,KLEN,PVAR)
    USE YOMMP0    , ONLY : MYPROC, NPROC
    USE MPL_MODULE, ONLY : MPL_BROADCAST
    USE YOMTAG    , ONLY : MTAGRAD
    USE EASY_NETCDF,ONLY : NETCDF_FILE
    USE YOMLUN    , ONLY : NULERR
    USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK

    IMPLICIT NONE

    ! Full name of the aerosol climatology file
    CHARACTER(LEN=*), INTENT(IN) :: CD_FILE_NAME
    ! Name of the variable in the file
    CHARACTER(LEN=*), INTENT(IN) :: CD_VAR_NAME
    ! Expected dimensions for input file
    INTEGER(KIND=JPIM), INTENT(IN) :: ILEN,JLEN,KLEN 
    ! Variable to be filled
    REAL(KIND=JPRB), INTENT(OUT) :: PVAR(ILEN,JLEN,KLEN)

    ! Local variable ..
    REAL(KIND=JPRB), ALLOCATABLE :: ZVAR(:,:,:)

    
    ! The NetCDF file containing the input climatology
    TYPE(NETCDF_FILE)  :: FILE

    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

    IF (LHOOK) CALL DR_HOOK('TM5M7_OPTICS_INIT:LOAD_TM5M7_OPTICS_DATA_3D',0,ZHOOK_HANDLE)

    ! First process actually loads the file - expected 3-dim array. 
    IF (MYPROC == 1) THEN
      CALL FILE%OPEN(TRIM(CD_FILE_NAME), IVERBOSE=4)
      CALL FILE%GET(TRIM(CD_VAR_NAME), ZVAR)
      CALL FILE%CLOSE()

      ! Check dimensions
      IF (SIZE(ZVAR,1) /= ILEN .OR. SIZE(ZVAR,2) /= JLEN &
           &                   .OR. SIZE(ZVAR,3) /= KLEN) THEN
        WRITE(NULERR,*) CD_VAR_NAME, 'in', CD_FILE_NAME, 'must be dimensioned (', &
             &  ILEN, ",", JLEN, ",", KLEN, ")"
        CALL ABOR1("Error reading TM5-M7 optics NetCDF file")
      ELSE
        PVAR(:,:,:)=ZVAR(:,:,:)
      ENDIF 
    ENDIF

    ! If there are more than one processes then broadcast the data
    ! (which means receive the data if MYPROC==1)
    IF (NPROC > 1) THEN
      CALL MPL_BROADCAST(PVAR, MTAGRAD, 1, &
           &   CDSTRING='TM5M7_OPTICS_INIT:LOAD_TM5M7_OPTICS_DATA_3D')
    ENDIF

    IF (LHOOK) CALL DR_HOOK('TM5M7_OPTICS_INIT:TM5M7_OPTICS_INIT:LOAD_TM5M7_OPTICS_DATA_3D',1,ZHOOK_HANDLE)

  END SUBROUTINE LOAD_TM5M7_OPTICS_DATA_3D



  ! Load an individual field  (or receive it from the first
  ! MPI task)
  SUBROUTINE LOAD_TM5M7_OPTICS_DATA_2D(CD_FILE_NAME, CD_VAR_NAME, ILEN,JLEN,PVAR)
    USE YOMMP0    , ONLY : MYPROC, NPROC
    USE MPL_MODULE, ONLY : MPL_BROADCAST
    USE YOMTAG    , ONLY : MTAGRAD
    USE EASY_NETCDF,ONLY : NETCDF_FILE
    USE YOMLUN    , ONLY : NULERR
    USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK

    IMPLICIT NONE

    ! Full name of the aerosol climatology file
    CHARACTER(LEN=*), INTENT(IN) :: CD_FILE_NAME
    ! Name of the variable in the file
    CHARACTER(LEN=*), INTENT(IN) :: CD_VAR_NAME
    ! Expected dimensions for input file
    INTEGER(KIND=JPIM), INTENT(IN) :: ILEN,JLEN 
    ! Variable to be filled
    REAL(KIND=JPRB), INTENT(OUT) :: PVAR(ILEN,JLEN)

    ! Local variable ..
    REAL(KIND=JPRB), ALLOCATABLE :: ZVAR(:,:)

    
    ! The NetCDF file containing the input climatology
    TYPE(NETCDF_FILE)  :: FILE

    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

    IF (LHOOK) CALL DR_HOOK('TM5M7_OPTICS_INIT:LOAD_TM5M7_OPTICS_DATA_2D',0,ZHOOK_HANDLE)

    ! First process actually loads the file - expected 2-dim array. 
    IF (MYPROC == 1) THEN
      CALL FILE%OPEN(TRIM(CD_FILE_NAME), IVERBOSE=4)
      CALL FILE%GET(TRIM(CD_VAR_NAME), ZVAR)
      CALL FILE%CLOSE()

      ! Check dimensions
      IF (SIZE(ZVAR,1) /= ILEN .OR. SIZE(ZVAR,2) /= JLEN) THEN
        WRITE(NULERR,*) CD_VAR_NAME, 'in', CD_FILE_NAME, 'must be dimensioned (', &
             &  ILEN, ",", JLEN, ")"
        CALL ABOR1("Error reading TM5-M7 optics NetCDF file")
      ELSE
        PVAR(:,:)=ZVAR(:,:)
      ENDIF 
    ENDIF

    ! If there are more than one processes then broadcast the data
    ! (which means receive the data if MYPROC==1)
    IF (NPROC > 1) THEN
      CALL MPL_BROADCAST(PVAR, MTAGRAD, 1, &
           &   CDSTRING='TM5M7_OPTICS_INIT:LOAD_TM5M7_OPTICS_DATA_2D')
    ENDIF

    IF (LHOOK) CALL DR_HOOK('TM5M7_OPTICS_INIT:TM5M7_OPTICS_INIT:LOAD_TM5M7_OPTICS_DATA_2D',1,ZHOOK_HANDLE)

  END SUBROUTINE LOAD_TM5M7_OPTICS_DATA_2D



  ! Load an individual field  (or receive it from the first
  ! MPI task)
  SUBROUTINE LOAD_TM5M7_OPTICS_DATA_1D(CD_FILE_NAME, CD_VAR_NAME, KLEN,PVAR)
    USE YOMMP0    , ONLY : MYPROC, NPROC
    USE MPL_MODULE, ONLY : MPL_BROADCAST
    USE YOMTAG    , ONLY : MTAGRAD
    USE EASY_NETCDF,ONLY : NETCDF_FILE
    USE YOMLUN    , ONLY : NULERR
    USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK

    IMPLICIT NONE

    ! Full name of the aerosol climatology file
    CHARACTER(LEN=*), INTENT(IN) :: CD_FILE_NAME
    ! Name of the variable in the file
    CHARACTER(LEN=*), INTENT(IN) :: CD_VAR_NAME
    ! Expected dimension lengths for input file
    INTEGER(KIND=JPIM), INTENT(IN) :: KLEN 
    ! Variable to be filled
    REAL(KIND=JPRB), INTENT(OUT) :: PVAR(KLEN)
    ! Variable to be filled
    REAL(KIND=JPRB), ALLOCATABLE :: ZVAR(:)

    
    ! The NetCDF file containing the input climatology
    TYPE(NETCDF_FILE)  :: FILE

    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

    IF (LHOOK) CALL DR_HOOK('TM5M7_OPTICS_INIT:LOAD_TM5M7_OPTICS_DATA_1D',0,ZHOOK_HANDLE)

    ! First process actually loads the file - expected 1-dim array.
    IF (MYPROC == 1) THEN
      CALL FILE%OPEN(TRIM(CD_FILE_NAME), IVERBOSE=4)
      CALL FILE%GET(TRIM(CD_VAR_NAME), ZVAR)
      CALL FILE%CLOSE()
      ! Check dimensions
      IF (SIZE(ZVAR,1) /= KLEN ) THEN
        WRITE(NULERR,*) CD_VAR_NAME, 'in', CD_FILE_NAME, 'must be dimensioned (', &
             &  KLEN, ")"
        CALL ABOR1("Error reading TM5-M7 optics NetCDF file")
      ELSE
        PVAR(:)=ZVAR(:)
      ENDIF 
    ENDIF

    ! If there are more than one processes then broadcast the data
    ! (which means receive the data if MYPROC==1)
    IF (NPROC > 1) THEN
      CALL MPL_BROADCAST(PVAR, MTAGRAD, 1, &
           &   CDSTRING='TM5M7_OPTICS_INIT:LOAD_TM5M7_OPTICS_DATA_1D')
    ENDIF

    IF (LHOOK) CALL DR_HOOK('TM5M7_OPTICS_INIT:TM5M7_OPTICS_INIT:LOAD_TM5M7_OPTICS_DATA_1D',1,ZHOOK_HANDLE)

  END SUBROUTINE LOAD_TM5M7_OPTICS_DATA_1D

END SUBROUTINE TM5M7_OPTICS_INIT




  !--------------------------------------------------------------------------
  !                    TM5                                                  !
  !--------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE:    OPTICS_WAVELEN_INIT
  !
  ! !DESCRIPTION: Initialise parameters which are depending on given 
  !               wavelengths.
  !\\
  !\\
  ! !INTERFACE:
  !
  SUBROUTINE TM5M7_OPTICS_WAVELEN_INIT( wdep,NWL )

    USE PARKIND1, ONLY : JPIM,JPRB
    USE TM5M7_OPTICS_DATA, ONLY : WAVELENDEP, OPAC,ECHAMHAM,SEGELSTEIN, &
       & opacdim, echamhamdim, segelsteindim
    USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK

    IMPLICIT NONE
    !
    ! !INPUT/OUTPUT PARAMETERS:
    !
    ! wavelength properties (wavelength itself and real/img part of refractive index)
    INTEGER(KIND=JPIM), INTENT(IN) :: NWL
    type(wavelendep), intent(inout), dimension(NWL) :: wdep 
    !
    !
    !EOP
    !------------------------------------------------------------------------
    !BOC
    
    INTEGER(KIND=JPIM) :: i, idx
    REAL(KIND=JPRB)    :: wl, h
    REAL(KIND=JPRB)    :: nscale, kscale

    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
    !-----------------------------------------------------------------------
    IF (LHOOK) CALL DR_HOOK('TM5M7_OPTICS_WAVELEN_INIT',0,ZHOOK_HANDLE)



    ! Real part n and imaginary part k of refractive index 
    ! for each wavelength:

    !>>> TvN
    ! The refractive index for 'sulfate' is taken from 
    ! the OPAC database (Hess et al., 1998; 
    ! Koepke et al., MPI report no. 243, 1997).
    ! The value used is that of 'sulfate solution',
    ! i.e. particles consisting of 75% of sulfuric acid (H2SO4),
    ! based on Fenn et al. (chapter 18 in Handbook of Geophysics
    ! and Space Environment, 1985).
    ! This is in line with Kinne et al. (JGR, 2003),
    ! who write that refractive indices for sulfate
    ! are usually based on 75% sulfuric acid solution.
    ! Actually, the OPAC refractive index agrees very well
    ! with the expression given by Boucher and Anderson (JGR, 1995)
    ! for H2SO4 at visible wavelengths.
    ! Thus, the OPAC data can be considered to apply to pure H2SO4,
    ! which is consistent of our application of mixing rules.
    
    ! For BC, the OPAC refractive index is scaled to one of
    ! the values proposed by Bond and Bergstrom (2006),
    ! valid at 550 nm (see their Table 5).
    ! Selected values for high, medium and low absorption are
    ! 1.95 + 0.79 i, 1.85 + 0.71 i, and 1.75 + 0.63 i, respectively.
    ! We select the low-absorption value,
    ! because it produces results in best agreement
    ! with AAOD from MODIS Collection 6 Deep Blue.
    ! Note that the medium-absorption value 1.85 + 0.71 i 
    ! was used in ECHAM simulations by Stier et al. (2007),
    ! and found to give reasonable results.
    ! Lowenthal et al. (Atmos. Environ., 2000)
    ! give a value of 1.90 + 0.6 i.
    ! The reference value used in the scaling
    ! is the OPAC value at 550 nm: 1.75 + 0.44 i.
    ! It would be better to include the scaling 
    ! in the input file.
    !nscale = 1.0 ! for using OPAC
    !kscale = 1.0 ! for using OPAC
    !nscale = 1.95/1.75
    !kscale = 0.79/0.44
    nscale = 1.85/1.75
    kscale = 0.71/0.44
    !nscale = 1.75/1.75
    !kscale = 0.63/0.44
    !<<< TvN

    !VH take from input: nwl = size(wdep)
    do i = 1, nwl
       wl = wdep(i)%wl

       ! Interpolate Opac data
       findOpac: Do idx = 1,opacdim
          If(opac(1,idx) .gt. wl) exit findOpac
       End Do findOpac
       If(idx .gt. opacdim) then
          idx = opacdim
          h = 1.0
       Else If(idx .eq. 1) then
          idx = 2
          h = 0.0
       Else
          h = (wl-opac(1,idx-1))/(opac(1,idx)-opac(1,idx-1))
       End If

       ! SO4
       wdep(i)%n(1) = opac(2,idx-1)+h*(opac(2,idx)-opac(2,idx-1))
       wdep(i)%k(1) = opac(3,idx-1)+h*(opac(3,idx)-opac(3,idx-1))
       ! BC
       !>>> TvN
       !wdep(i)%n(2) = opac(4,idx-1)+h*(opac(4,idx)-opac(4,idx-1))
       !wdep(i)%k(2) = opac(5,idx-1)+h*(opac(5,idx)-opac(5,idx-1))
       wdep(i)%n(2) = nscale * ( opac(4,idx-1)+h*(opac(4,idx)-opac(4,idx-1)) )
       wdep(i)%k(2) = kscale * ( opac(5,idx-1)+h*(opac(5,idx)-opac(5,idx-1)) )
       !<<< TvN
       ! SS 
       wdep(i)%n(5) = opac(6,idx-1)+h*(opac(6,idx)-opac(6,idx-1))
       wdep(i)%k(5) = opac(7,idx-1)+h*(opac(7,idx)-opac(7,idx-1))

       ! Interpolate ECHAM-HAM data
       findEchamham: Do idx = 1,echamhamdim
          If(echamham(1,idx) .gt. wl) exit findEchamham
       End Do findEchamham
       If(idx .gt. echamhamdim) then
          idx = echamhamdim
          h = 1.0
       Else If(idx .eq. 1) then
          idx = 2
          h = 0.0
       Else
          h = (wl-echamham(1,idx-1))/(echamham(1,idx)-echamham(1,idx-1))
       End If

       ! OC 
       !>>> TvN
       ! The 'ECHAM-HAM' data currently used for POM
       ! are based on the data from Fenn et al. (1985) 
       ! for the 'water soluble' component,
       ! but at a reduced number of wavelengths up to 15 um.
       ! It would be better to use the original table
       ! in the input file.
       !<<< TvN
       wdep(i)%n(3) = echamham(2,idx-1)+h*(echamham(2,idx)-echamham(2,idx-1))
       wdep(i)%k(3) = echamham(3,idx-1)+h*(echamham(3,idx)-echamham(3,idx-1))


       ! SOA
       ! >>TB
       ! For now (March 2017) we use OC indices for SOA
       ! <<TB
       wdep(i)%n(4) = echamham(2,idx-1)+h*(echamham(2,idx)-echamham(2,idx-1))
       wdep(i)%k(4) = echamham(3,idx-1)+h*(echamham(3,idx)-echamham(3,idx-1))

       ! DU
       wdep(i)%n(6) = echamham(4,idx-1)+h*(echamham(4,idx)-echamham(4,idx-1))
       wdep(i)%k(6) = echamham(5,idx-1)+h*(echamham(5,idx)-echamham(5,idx-1))

       ! Interpolate Segelstein data
       findSegelstein: Do idx = 1,segelsteindim
          If(segelstein(1,idx) .gt. wl) exit findSegelstein
       End Do findSegelstein
       If(idx .gt. segelsteindim) then
          idx = segelsteindim
          h = 1.0
       Else If(idx .eq. 1) then
          idx = 2
          h = 0.0
       Else
          h = (wl-segelstein(1,idx-1))/(segelstein(1,idx)-segelstein(1,idx-1))
       End If

       ! H2O
       wdep(i)%n(7) = segelstein(2,idx-1)+h*(segelstein(2,idx)-segelstein(2,idx-1))
       wdep(i)%k(7) = segelstein(3,idx-1)+h*(segelstein(3,idx)-segelstein(3,idx-1))
    enddo


  IF (LHOOK) CALL DR_HOOK('TM5M7_OPTICS_WAVELEN_INIT',1,ZHOOK_HANDLE)
  END SUBROUTINE TM5M7_OPTICS_WAVELEN_INIT
