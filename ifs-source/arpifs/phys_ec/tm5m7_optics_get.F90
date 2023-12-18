SUBROUTINE TM5M7_OPTICS_GET(m_eff, xg, Cext, a, g, cext_table, a_table, g_table,status)

!*** * TM5M7_OPTICS_GET* 
!
!
!-----------------------------------------------------------------------------
!                    TM5                                                     !
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE:	TM5M7_OPTICS_GET
!
! !DESCRIPTION: Main optical properties routine. Here the interpolated 
!               values for extinction coefficient, single scattering
!               albedo and assymetry parameter are retrieved and returned. 
!
!\\
!\\
  
!     SOURCE.
!     -------
!
!     Taken from TM5-code
!
!     MODIFICATIONS.
!     --------------
!
!
! !REVISION HISTORY: 
!
!   12 Aug 2008 - Michael Kahnert, SMHI
!    6 Feb 2011 - Achim Strunk -
!
!      Sep 2021 - V. Huijnen: first introduction into OpenIFS
!
!
!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMLUN    ,ONLY : NULERR

USE TM5M7_OPTICS_DATA, ONLY : LKVAL, KVAL,N1R,N_RIR,N_RII,XS,N_X

IMPLICIT NONE

!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------



!
! !INPUT PARAMETERS:
!
complex,                           INTENT(IN)    :: m_eff
REAL(KIND=JPRB),                   INTENT(IN)    :: xg
REAL(KIND=JPRB), dimension(:,:,:), INTENT(IN)    :: cext_table, a_table, g_table
!
! !OUTPUT PARAMETERS:
!
REAL(KIND=JPRB),             INTENT(OUT)   :: Cext, a, g
INTEGER(KIND=JPIM),          INTENT(INOUT) :: status
    




!*       0.2   LOCAL VARIABLES
!              ---------------
REAL(KIND=JPRB)    :: n, k, n1, k1, dk1, dk2, lk
REAL(KIND=JPRB)    :: modrad, modrad1, dr, dr1, slope, cext1, a1, g1
INTEGER(KIND=JPIM) :: i
INTEGER(KIND=JPIM) :: i_n, i_k, i_knew


REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TM5M7_OPTICS_GET',0,ZHOOK_HANDLE)



    status = 0
    !FAILS when K=1.00 because kval(15)=1.0 and k.lt.kval(15) false
    !k1=15

    n=real(m_eff)
    k=imag(m_eff)
    if(k.lt.kval(1))then
       k1=kval(1)
       i_knew=1
    elseif(k.gt.kval(15))then
       k1=kval(15)
       i_knew=15
    else
       get_k: do i=2,15
          if(k.lt.kval(i))then
             dk1=k-kval(i-1)
             dk2=kval(i)-k
             if(dk1.le.dk2)then
                k1=kval(i-1)
                i_knew=i-1
             else
                k1=kval(i)
                i_knew=i
             endif
             exit get_k
          elseif (i==15 .and. abs(kval(i)-k) <1e-20)then
             k1=kval(i)
             i_knew=i
          endif
       enddo get_k
    endif
    lk=-log10(k1)
    !#n1=float(int(50*n+0.5))/50.
    !do i_n = 1, n_rir
    !   if (abs(n1 - n1r(i_n)) < 1e-4) exit
    !enddo

    !i_n = 1 + int((n-1.12)*50+0.5)
    !AJS: I guess n is a number on the (increasing) n1r axis; search nearest index:
    i_n = size(n1r)
    do i = 1, size(n1r)
       if ( n <= n1r(i) ) then
          i_n = i
          exit
       end if
    end do
    if ( i_n > 1 ) then
       if ( abs(n-n1r(i_n-1)) < abs(n-n1r(i_n)) ) i_n = i_n - 1
    else if ( i_n < n_rir ) then
       if ( abs(n-n1r(i_n+1)) < abs(n-n1r(i_n)) ) i_n = i_n + 1
    end if

    do i_k = 1, n_rii
       if (abs(lk - lkval(i_k)) < 1e-4) exit
    enddo

    !     ! PLS - test : ik can be determined without the preceding loop "do i_k = 1, n_rii"
    !     if (i_k.ne.i_knew) then
    !        status = 1
    !        write (NULERR,'(" PLSPLS  ik NE iknew = ",2(i2,2x))')i_k,i_knew 
    !     endif

    ! following the "get_k" loop above, the only way to get into this is to
    ! have a NaN for k in the first place ? 
    if (i_k > n_rii) then
       status = 1
       write(NULERR,*)'FATAL ERROR: i_k value outside LUT'
               write (NULERR,'("    lk = ",E16.4)')lk 
               write (NULERR,'("    k1 = ",E16.4)')k1 
               write (NULERR,'("     k = ",E16.4)')k  
               write (NULERR,'("     n = ",E16.4)')n  
               write (NULERR,'("     i_k = ",I6)')i_k  
               write (NULERR,'("     n_rii = ",I6)')n_rii  

               do i_n = 1, 15
                  write (NULERR,'("    lkval(",i2,") = ",E16.4)')i_n,lkval(i_n)
       !           call goPr
                  write (NULERR,'("     kval(",i2,") = ",E16.4)')i_n,kval(i_n)
       !           call goPr
               enddo

       IF(LHOOK) CALL DR_HOOK('TM5M7_OPTICS_GET',1,ZHOOK_HANDLE)
       RETURN
    endif

    ! Added check (15-7-2010 - P. Le Sager)
    if (i_n > n_rir) then
       status = 1
       WRITE(NULERR,*)'FATAL ERROR: i_n value out of range'
       IF(LHOOK) CALL DR_HOOK('TM5M7_OPTICS_GET',1,ZHOOK_HANDLE)
       RETURN
    endif


    !>>> TvN
    ! Since xs(1) equals zero a problem may occur when xg is negative,
    ! because modrad.gt.xg would become true for i=1 in that case,
    ! while modrad1, Cext1, a1 and g1 are not yet set.
    ! It is not clear if negative xg values can ever occur,
    ! but if they do that should be prevented when calculating rg.
    ! 
    ! However, the problem may perhaps already occur
    ! when xg equals zero, because of the finite machine precision.
    !
    ! In any case, it is desired to initialize modrad1, Cext1, a1 and g1.
    ! Cext1, a1 and g1 should be set to their table entries for i=1,
    ! which are all zero:
    Cext1 = cext_table(1,i_n,i_k)
    a1 = a_table(1,i_n,i_k)
    g1 = g_table(1,i_n,i_k)
    ! modrad1 can be initialized to any value different from xs(1),
    ! to prevent division by zero.
    modrad1=xs(1)-9.99e-4
    ! This combination will force slope to zero
    ! and Cext, a and g to the table entries for i=1 (zero),
    ! in case modrad.gt.xg is true for i=1.
    !<<< TvN
    get_values: do i = 1, n_x
       modrad = xs(i)
       a = a_table(i,i_n, i_k)
       g = g_table(i,i_n, i_k)
       cext = cext_table(i,i_n, i_k)
       ! With small concentrations, like in the stratosphere, M7 does not trust its radii.
       ! See m7_averageproperties -> zinsvol. The M7-radius goes to zero
       ! Modrad may never be larger than xg the first time. Modrad1 is not set,
       ! will be zero (or something worse) and dr will be zero (or something worse).
       ! slope is demolished, makeing all values NaN. Therefore, it is modrad .gt. xg
       ! instead of modrad .ge. xg
       !
       ! PLS-TRANSLATION - It boils down to: "on the first iteration of this loop, if modrad=xg and
       !  you go into the if-statement below, you are in trouble, because modrad1 can be
       !  anything. To prevent that, replace GE with GT to avoid going into the if-statement at the
       !  first iteration."
       !
       if(modrad.gt.xg)then
          dr=modrad-modrad1
          dr1=xg-modrad1
          slope=(Cext-Cext1)/dr
          Cext=Cext1+slope*dr1
          slope=(a-a1)/dr
          a=a1+slope*dr1
          slope=(g-g1)/dr
          g=g1+slope*dr1
          exit get_values 
       endif
       modrad1=modrad
       Cext1=Cext
       a1=a
       g1=g
    enddo get_values



IF(LHOOK) CALL DR_HOOK('TM5M7_OPTICS_GET',1,ZHOOK_HANDLE)
END SUBROUTINE TM5M7_OPTICS_GET

