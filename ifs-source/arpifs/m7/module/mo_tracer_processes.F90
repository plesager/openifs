!! SPDX-License-Identifier: BSD-3-Clause
!! Copyright (c) 2021 hammoz, MPI fuer Meteorologie, FZJ 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE mo_tracer_processes

  USE mo_kind,             ONLY: wp 
  IMPLICIT NONE

  
  PRIVATE

  !! Interface routines         ! purpose                   ! called by
  PUBLIC :: xt_conv_massfix     ! adjust tracer mass in convection
  PUBLIC :: xt_borrow           ! another mass fixer

!#ifdef HAMMOZ
  !--- module variables
  REAL(wp), ALLOCATABLE :: zxtte_old(:,:,:)          ! for xt_convmassfix
  !$OMP THREADPRIVATE(zxtte_old)
!#endif
  CONTAINS

  SUBROUTINE xt_conv_massfix (kproma,        kbdim,             klev,         &
                              klevp1,        ktrac,             krow,         &
                              papp1,         paphp1,            pxtte,        &
                              loini, pxtbound)

    ! *xt_massfix* corrects the tendencies of each column to
    !              conserve mass
    !
    ! Author:
    ! -------
    ! Philip Stier, MPI-Met,         11/2003
    !
    ! Method:
    ! -------
    ! To keep the model as closely as possible to the intended physical
    ! tendency changes, the correction of the tendencies is
    ! imposed proportionally to the absolute values of the
    ! applied tendency in each layer.
    !
    ! Total mass error tendency [kg m-2 s-1]:
    !
    !    dxtdt(jl)=SUM(xtte(jl,jk)*dpg(jl,jk))+xtbound(jl,jk)
    !
    !    where xtbound is a boundary condition, i.e. deposition flux
    !
    ! Correction proportional to tendency [kg m-2 s-1]:
    !
    !                         - ABS(xtte(jl,jk)*dpg(jl,jk))
    !    dxtfix(jl,jk)   = --------------------------------- * dxdt(jl)
    !                       SUM(ABS(xtte(jl,jk)*dpg(jl,jk)))
    !
    ! Resulting in a corrective tendency [kg kg-1 s-1]:
    !
    !    zxttefix(jl,jk) = dxtfix(jl,jk) / dpg(jl,jk)
    !
    !
    ! Arguments:
    ! ----------
    !
    ! pxtte     = cumulative tendency for all tracers [kg kg-1 s-1]
    !
    ! Usage:
    ! ------
    ! Call twice: first with loini=.TRUE. to store the old tendency
    !             then to fix the mass of processes that have modified
    !             pxtte meanwhile with the current pxtte and loini=.FALSE.

  USE mo_tracdef,             ONLY: trlist             ! tracer info variable
  USE mo_physical_constants,  ONLY: grav

    IMPLICIT NONE


    !--- Arguments:
    INTEGER,      INTENT(IN)    :: kproma, kbdim, klev, klevp1, ktrac, krow
    LOGICAL,      INTENT(IN)    :: loini
    REAL(wp),     INTENT(IN)    :: papp1(kbdim,klev),       &
                                   paphp1(kbdim,klevp1)
    REAL(wp),     INTENT(INOUT) :: pxtte(kbdim,klev,ktrac)
!>>SF
    REAL(wp),     INTENT(IN)    :: pxtbound(kbdim,ktrac) ! boundary condition (wet deposition) [kg m-2 s-1]
!<<SF


    !--- Local Variables:
    INTEGER               :: jl, jk, jt
    REAL(wp)              :: zeps, zxttefix
    REAL(wp)              :: zdxtdt(kbdim), zdxtdtsum(kbdim)
    REAL(wp)              :: zxtte(kbdim,klev), zdpg(kbdim,klev)

    !--- 1) Initialization mode:

    IF (.NOT. ALLOCATED(zxtte_old)) ALLOCATE(zxtte_old(kbdim,klev,ktrac))
    IF (loini) THEN
       zxtte_old(1:kproma,:,:) = pxtte(1:kproma,:,:)
       RETURN             ! initialisation done: return from subroutine
    END IF

    !--- 2) Mass fix mode:

    zeps=EPSILON(1.0_wp)

    !--- 2.1) Calculate auxiliary variable dp/g :
    !--- Uppermost level:
    zdpg(1:kproma,1)=2._wp*(paphp1(1:kproma,2)-papp1(1:kproma,1))/grav
    !--- Other levels:
    DO jk=2, klev
       zdpg(1:kproma,jk)=(paphp1(1:kproma,jk+1)-paphp1(1:kproma,jk))/grav
    END DO

    !--- 2.2) Apply mass fixer
    DO jt=1, ktrac
!!     IF(trlist%ti(jt)%nwetdep > 0) THEN
       IF(trlist%ti(jt)%nconvmassfix > 0) THEN

          zdxtdt(1:kproma)   =0.0_wp
          zdxtdtsum(1:kproma)=0.0_wp

          !--- Accumulated tendency since initialization:
          zxtte(1:kproma,:)=pxtte(1:kproma,:,jt)-zxtte_old(1:kproma,:,jt)

          !--- Calculate vertically integrated mass error tendency [kg m-2 s-1]:
          DO jk=1, klev
             DO jl=1, kproma

                zdxtdt(jl)=zdxtdt(jl)      +zxtte(jl,jk)*zdpg(jl,jk)

                zdxtdtsum(jl)=zdxtdtsum(jl)+ABS(zxtte(jl,jk)*zdpg(jl,jk))

             END DO
          END DO

!>>SF
          zdxtdt(1:kproma)=zdxtdt(1:kproma)+pxtbound(1:kproma,jt) !SF restore boundary cond. calc.
!<<SF

          DO jk=1, klev
             DO jl=1, kproma
                IF (ABS(zdxtdtsum(jl)) > zeps) THEN

                   zxttefix=-((ABS(zxtte(jl,jk)*zdpg(jl,jk))/zdxtdtsum(jl))*zdxtdt(jl))/zdpg(jl,jk)

                   pxtte(jl,jk,jt)=pxtte(jl,jk,jt)+zxttefix

                END IF
             END DO ! kproma
          END DO ! klev

       END IF ! nwetdep
    END DO ! ktrac


  END SUBROUTINE xt_conv_massfix
  
  SUBROUTINE xt_borrow(kproma, kbdim,  klev,  klevp1, ktrac,       &
                       papp1,  paphp1,                             &
                       pxtm1,  pxtte                               )

    ! *xt_borrow* borrowing scheme to correct for negative
    !             tracer masses conserveing mass within the
    !             column
    ! Note by SF: actually, the scheme is also good for non-mass tracers, therefore it is
    ! now applied to all transported tracers.
    !
    ! Author:
    ! -------
    ! Philip Stier, MPI-Met,         06/2004
    !
    ! Method:
    ! -------
    ! If negative tracer values occur the scheme is iteratively
    ! borrowing tracer mixing ratios from the grid-box below.
    ! If the lowest layer can not compensate for accumulated remaining
    ! corrections, the procedure is repeated from bottom to top.
    !
    ! Restrictions:
    ! -------------
    ! Columns with a negative total integrated tracer
    ! content are set to zero despite the associated mass
    ! error.

    USE mo_physical_constants,    ONLY: grav
    USE mo_time_control, ONLY: time_step_len
    USE mo_tracdef,      ONLY: trlist, AEROSOLMASS, GAS, AEROSOLNUMBER
    USE mo_advection,    ONLY: no_advection !SF #246
    USE mo_submodel,     ONLY: lmoz  !csld #330

    IMPLICIT NONE

    INTEGER, INTENT(IN)    :: kproma, kbdim, klev, klevp1, ktrac
    REAL(wp),INTENT(IN)    :: papp1(kbdim,klev), paphp1(kbdim,klevp1)
    REAL(wp),INTENT(IN)    :: pxtm1(kbdim,klev,ktrac)
    REAL(wp),INTENT(INOUT) :: pxtte(kbdim,klev,ktrac)

    !--- Local variables:

    INTEGER :: jl, jk, jt

    REAL(wp):: zxtp1, ztmst, zeps

    REAL(wp):: zxtbor(kbdim)

    REAL(wp):: zdpg(kbdim,klev)

    LOGICAL :: lborrtrac !csld #330 : determines if the tracer can be borrowed or not

    !--- 0) Initializations:
    ztmst=time_step_len
    zeps=10._wp*EPSILON(1.0_wp)

    !--- 1) Calculate auxiliary variable dp/g :

    !--- Uppermost level:

    zdpg(1:kproma,1)=2._wp*(paphp1(1:kproma,2)-papp1(1:kproma,1))/grav

    !--- Other levels:

    DO jk=2, klev

       zdpg(1:kproma,jk)=(paphp1(1:kproma,jk+1)-paphp1(1:kproma,jk))/grav

    END DO

    !--- 2) Borrowing scheme:
    
    DO jt=1, ktrac

       !>>csld #330
       ! use xt_borrow on all aerosols species
       lborrtrac = trlist%ti(jt)%nphase == AEROSOLMASS .OR. trlist%ti(jt)%nphase == AEROSOLNUMBER

       IF (.NOT. lmoz) THEN    ! in case of "ham only" runs, use of xt_borrow on gaz species as well
          lborrtrac = lborrtrac .OR. (trlist%ti(jt)%nphase == GAS)
       END IF

       IF ( (trlist%ti(jt)%ntran /= no_advection)       .AND. &  !SF #246 exclude non-advected tracers
            lborrtrac) THEN
       !<<csld #330

          !--- Start borrowing scheme:
    
          !--- Integrate from top to bottom:

          zxtbor(1:kproma)=0.0_wp

          DO jk=1, klev
             DO jl=1, kproma

                !--- Check if corrected mass, including the fix of the layer above,
                !    yields a negative mixing ratio:
                !    (Convert mass correction converted to [kg kg-1] in current layer)

                zxtp1=pxtm1(jl,jk,jt)+pxtte(jl,jk,jt)*ztmst + zxtbor(jl)/zdpg(jl,jk)

                IF ( zxtp1 > 0.0_wp ) THEN

                   !--- Subtract corrected tracer mass from current layer:
                   !    (zxtbor is negative)
       
                   pxtte(jl,jk,jt)=pxtte(jl,jk,jt)+zxtbor(jl)/(zdpg(jl,jk)*ztmst)
    
                   !--- Reset mass correction:
    
                   zxtbor(jl)=0.0_wp

                ELSE

                   !--- Adjust tendency to yield zero:

                   pxtte(jl,jk,jt)=-pxtm1(jl,jk,jt)/ztmst

                   !--- Add correcting mass mixing ratio and convert to [kg m-2]:
                   !    (implicit summation due to the inclusion
                   !     of zxtbor(jk-1) in zxtp1)

                   zxtbor(jl)=zxtp1*zdpg(jl,jk)

                END IF

             END DO
          END DO

          !--- If surface layer cannot compensate accumulated correction:
          !    Iterate from bottom to top

          DO jk=klev, 1, -1
             DO jl=1, kproma
                IF (zxtbor(jl) < -zeps) THEN

                   !--- Check if corrected mass, including the fix of the layer below,
                   !    yields a negative mixing ratio:
                   !    (Convert mass correction converted to [kg kg-1] in current layer)

                   zxtp1=pxtm1(jl,jk,jt)+pxtte(jl,jk,jt)*ztmst + zxtbor(jl)/zdpg(jl,jk)

                   IF ( zxtp1 > 0.0_wp ) THEN

                      !--- Subtract corrected tracer mass from current layer:
                      !    (zxtbor is negative)

                      pxtte(jl,jk,jt)=pxtte(jl,jk,jt)+zxtbor(jl)/(zdpg(jl,jk)*ztmst)

                      !--- Reset mass correction:

                      zxtbor(jl)=0.0_wp

                   ELSE

                      !--- Adjust tendency to yield zero:

                      pxtte(jl,jk,jt)=-pxtm1(jl,jk,jt)/ztmst
                   
                      !--- Add correcting mass mixing ratio and convert to [kg m-2]:
                      !    (implicit summation due to the inclusion
                      !     of zxtbor(jk-1) in zxtp1)
                   
                      zxtbor(jl)=zxtp1*zdpg(jl,jk) 
                   
                   END IF
                   
                END IF
             END DO 
          END DO
       
       END IF 
    END DO
          
  END SUBROUTINE xt_borrow
          




  
!!mgs!!   -------   obsolete code  (perhaps re-use some parts of this to enhance burden diag? ------  
!!mgs!!   SUBROUTINE trastat
!!mgs!! 
!!mgs!!   ! Description:
!!mgs!!   !
!!mgs!!   ! Prints out accumulated mass budgets for tracers at
!!mgs!!   ! the end of a run
!!mgs!! 
!!mgs!!   USE mo_control,        ONLY: ngl
!!mgs!!   USE mo_mpi,            ONLY: p_sum, p_communicator_d, p_pe, p_io
!!mgs!! 
!!mgs!!   !  Local scalars: 
!!mgs!!   REAL(wp) :: zmglob, zmnhk, zmshk, zmstrat, zmtrop, zqcount
!!mgs!!   INTEGER ::  jt
!!mgs!! 
!!mgs!!   !  Local arrays: 
!!mgs!!   REAL(wp) :: zmstratn(ntrac+1), zmstrats(ntrac+1), zmtropn(ntrac+1),    &
!!mgs!!               zmtrops(ntrac+1)
!!mgs!! 
!!mgs!!   !  Intrinsic functions 
!!mgs!!   INTRINSIC SUM
!!mgs!! 
!!mgs!! 
!!mgs!!   !  Executable statements 
!!mgs!! 
!!mgs!!   zqcount = 1.0_wp/(REAL(icount,wp))
!!mgs!!   IF (p_pe == p_io) THEN
!!mgs!!     CALL message('',' Tracer mass budget:')
!!mgs!!     CALL message('',separator)
!!mgs!!     CALL message('',' Averaged mass budgets in [kg] ')
!!mgs!!     CALL message('',' global   n-hem  s-hem  tropo  strat  n-tro s-tro  n-str s-str ')
!!mgs!!   ENDIF
!!mgs!!   DO jt = 1, ntrac + 1
!!mgs!!      zmtropn(jt)  = SUM(tropm(1:ngl/2,jt))  * zqcount
!!mgs!!      zmstratn(jt) = SUM(stratm(1:ngl/2,jt)) * zqcount
!!mgs!!      zmtrops(jt)  = SUM(tropm(ngl/2+1:ngl,jt))  * zqcount
!!mgs!!      zmstrats(jt) = SUM(stratm(ngl/2+1:ngl,jt)) * zqcount
!!mgs!!   END DO
!!mgs!!   zmtropn  = p_sum (zmtropn,  p_communicator_d)
!!mgs!!   zmstratn = p_sum (zmstratn, p_communicator_d)
!!mgs!!   zmtrops  = p_sum (zmtrops,  p_communicator_d)
!!mgs!!   zmstrats = p_sum (zmstrats, p_communicator_d)
!!mgs!! 
!!mgs!!   IF (p_pe == p_io) THEN    
!!mgs!!     DO jt = 1, ntrac + 1
!!mgs!!        zmnhk   = zmtropn(jt)  + zmstratn(jt)
!!mgs!!        zmshk   = zmtrops(jt)  + zmstrats(jt)
!!mgs!!        zmtrop  = zmtropn(jt)  + zmtrops(jt)
!!mgs!!        zmstrat = zmstratn(jt) + zmstrats(jt)
!!mgs!!        zmglob  = zmnhk + zmshk
!!mgs!! 
!!mgs!!        IF (jt <= ntrac) THEN
!!mgs!!          WRITE (message_text,'(a,i2,9e9.2)')                    &
!!mgs!!               ' Tracer: ', jt, zmglob, zmnhk,                   &
!!mgs!!               zmshk, zmtrop, zmstrat, zmtropn(jt), zmtrops(jt), &
!!mgs!!               zmstratn(jt), zmstrats(jt)
!!mgs!!          CALL message('',message_text)
!!mgs!!        ELSE
!!mgs!!          WRITE (message_text,'(a   ,9e9.2)')                    &
!!mgs!!               ' Air mass: ', zmglob, zmnhk,                     &
!!mgs!!               zmshk, zmtrop, zmstrat, zmtropn(jt), zmtrops(jt), &
!!mgs!!               zmstratn(jt), zmstrats(jt)
!!mgs!!          CALL message('',message_text)
!!mgs!!        END IF
!!mgs!!     END DO
!!mgs!!   END IF
!!mgs!!   
!!mgs!!   END SUBROUTINE trastat
!!mgs!! 
  
END MODULE mo_tracer_processes

