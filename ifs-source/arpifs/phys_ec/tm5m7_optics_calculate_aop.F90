SUBROUTINE TM5M7_OPTICS_CALCULATE_AOP(KIDIA,KFDIA,KLON,KLEV, nwl,NCONTR, wdep, ecearth_units, &
                          &         AOP_IN, &
                          &         PAOP_OUT_EXT, PAOP_OUT_A, PAOP_OUT_G, PAOP_OUT_ADD )

!*** * TM5M7_OPTICS_CALCULATE_AOP* 
!
!
!-----------------------------------------------------------------------------
!                    TM5                                                     !
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE:	OPTICS_CALCULATE_AOP
!
! !DESCRIPTION: This routine writes on PAOP_OUT_* (module wide parameters) 
!		the retrieved aerosol properties. The caller has to ensure
!		that these fields have been allocated properly.
!   IMPORTANT:  OC is actually POM. 
!		Remember converting OC to POM when sending it to this method.
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
!   12 Aug 2008 - Michael Kahnert, SMHI
!    6 Feb 2011 - Achim Strunk -
!
!      Sep 2021 - V. Huijnen: first introduction into OpenIFS
!
!
!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE YOMLUN    ,ONLY : NULOUT,NULERR

USE YOMCST, ONLY : RPI
USE TM5M7_DATA, ONLY : SIGMA, H2SO4_FACTOR,NH4NO3_FACTOR, SO4_DENSITY,MSA_DENSITY, &
  & CMEDR2MMEDR,CARBON_DENSITY,POM_DENSITY,SOA_DENSITY,SS_DENSITY,DUST_DENSITY,&
  & MODAL_DATA,NMOD,NSOL
USE TM5M7_OPTICS_DATA, ONLY: WAVELENDEP, AOPI,NADD, &
  & cext_159, a_159, g_159, cext_200, a_200, g_200

IMPLICIT NONE

!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

INTEGER(KIND=JPIM)              , INTENT(IN) :: KIDIA,KFDIA,KLON,KLEV,NWL,NCONTR
TYPE(WAVELENDEP), DIMENSION(NWL), INTENT(IN) :: WDEP
LOGICAL                         , INTENT(IN) :: ECEARTH_UNITS
TYPE(AOPI), DIMENSION(KLON,KLEV), INTENT(IN) :: AOP_IN 


REAL(KIND=JPRB), INTENT(OUT):: PAOP_OUT_EXT(KLON,KLEV,NWL,NCONTR)
REAL(KIND=JPRB), INTENT(OUT):: PAOP_OUT_A(KLON,KLEV,NWL)
REAL(KIND=JPRB), INTENT(OUT):: PAOP_OUT_G(KLON,KLEV,NWL)
REAL(KIND=JPRB), INTENT(OUT), OPTIONAL :: PAOP_OUT_ADD(KLON,KLEV,NWL,NADD) ! additional parameters



!*       0.2   LOCAL VARIABLES
!              ---------------
REAL(KIND=JPRB), DIMENSION(KLON,KLEV) :: NCsca, incext
REAL(KIND=JPRB)             :: Cexti, ai, gi, NCscai, xg
COMPLEX                     :: m_eff
REAL(KIND=JPRB), DIMENSION(:),allocatable :: lnsigma

INTEGER(KIND=JPIM) :: i, imode, JL,JK
INTEGER(KIND=JPIM) :: statusomp
LOGICAL            :: coarse
REAL(KIND=JPRB)   :: totvoldry, modfrac
REAL(KIND=JPRB), DIMENSION(:,:,:), Pointer :: cext_table, a_table, g_table
REAL(KIND=JPRB)    :: TWOPI

REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-----------------------------------------------------------------------

#include "tm5m7_get_refr_idx.intfb.h"
#include "tm5m7_optics_get.intfb.h"


!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TM5M7_OPTICS_CALCULATE_AOP',0,ZHOOK_HANDLE)
  

    TWOPI=2*RPI

!    allocate( NCsca ( KLON,KLEV ) )
!    allocate( incext( KLON,KLEV ) )
    STATUSOMP=0


    !     Sulphate based on OPAC (Hess et al., 1998):

    !=======================================================================
    !     Get refractive indices of each component at the given wavelength:
    !=======================================================================
    do i = 1, nwl   ! loop over wavelengths

       if( wdep(i)%split .or. wdep(i)%insitu ) then 
          allocate( lnsigma( nmod ) )
          lnsigma = log(sigma)
       end if

       PAOP_OUT_EXT(KIDIA:KFDIA,1:KLEV,i,:)  = 0.0_JPRB
       PAOP_OUT_G  (KIDIA:KFDIA,1:KLEV,i)    = 0.0_JPRB
       PAOP_OUT_A  (KIDIA:KFDIA,1:KLEV,i)    = 0.0_JPRB
       if( wdep(i)%insitu ) PAOP_OUT_ADD(KIDIA:KFDIA,1:KLEV,i,:)  = 0.0_JPRB 
       NCsca      (KIDIA:KFDIA,1:KLEV)      = 0.0_JPRB

       do imode = 1, 7        ! loop over M7 modes
          coarse = (imode .eq. 4 .or. imode .eq. 7)

          If (coarse) then
             cext_table => cext_200
             a_table => a_200
             g_table => g_200
          Else
             cext_table => cext_159
             a_table => a_159
             g_table => g_159
          End if


          !=======================================================================
          !     Compute effective refractive index of internally mixed aerosols
          !     for each grid cell and mode:
          !=======================================================================
          DO JK=1,KLEV
          DO JL=KIDIA,KFDIA
             !write(1001,*) aop_in(JL,JK)%SO4(imode)
             !write(1002,*) aop_in(JL,JK)%SOA(imode)
             !write(1003,*) aop_in(JL,JK)%BC(imode)
             !write(1005,*) aop_in(JL,JK)%OC(imode)
             !write(1004,*) aop_in(JL,JK)%SS(imode)
             !write(1006,*) aop_in(JL,JK)%DU(imode)
             !write(1007,*) aop_in(JL,JK)%NO3(imode)
             call tm5m7_get_refr_idx( wdep(i), &
                  aop_in(JL,JK)%SO4(imode) + aop_in(JL,JK)%NO3(imode), & ! H2-SO4 + NH4NO3
                  aop_in(JL,JK)%BC (imode),  aop_in(JL,JK)%OC (imode), &
                  aop_in(JL,JK)%SOA (imode), &
                  aop_in(JL,JK)%SS (imode),  aop_in(JL,JK)%DU (imode), &
                  aop_in(JL,JK)%h2o(imode),  imode, m_eff)
             !if (statusomp ==1) then
             !   WRITE(NULERR,'(" Problem with GET_REFR_IDX-NAN inputs:  ", 9(E16.4,2x))') &
             !        aop_in(JL,JK)%SO4(imode), aop_in(JL,JK)%NO3(imode), aop_in(JL,JK)%h2o(imode), &
             !        aop_in(JL,JK)%SOA(imode),aop_in(JL,JK)%BC(imode), aop_in(JL,JK)%OC(imode),&
             !        aop_in(JL,JK)%SS(imode), aop_in(JL,JK)%DU(imode)
             !   ! Assume edge case not caught. Set m_eff to default (1.0,1.0e-9), same as in ref_idx
             !   ! (subroutine problems with inputs SO4 and WATER?)
             !   m_eff= Cmplx(1.0,1.0e-9)
             !   ! If all edge cases were handled correctly, that should be:
             !   !CALL ABOR1("Problem with tm5m7_get_refr_idx" )
             !endif

             ! cmk added towpi for new netcdf lookup table
             xg = twopi*aop_in(JL,JK)%rg(imode) / wdep(i)%wl

             !=======================================================================
             !     get aerosol optical properties from data base for each mode
             !=======================================================================
             ! add extra safeguard for negative xg. It is unphysical, but have a safeguard anyway
             ! added for refr_idx_nan problem 

             if (xg .gt. 0.0) then
                call TM5M7_OPTICS_GET(m_eff, xg, Cexti, ai, gi, cext_table, a_table, g_table, statusomp )
             else
                Cexti=0.0
                ai=1.0
                gi=1.0
             end if
             if (statusomp ==1) then
                call ABOR1("ERROR due to TM5M7_OPTICS_GET")
             endif

             ! Multiply Cext with lambda^2 to get the cross section.
             Cexti = Cexti*(wdep(i)%wl**2)
             ! this here is extinction coefficient in this mode 
             incext(JL,JK) = aop_in(JL,JK)%numdens(imode) * Cexti
             ! sum up partial coefficients
             PAOP_OUT_EXT(JL,JK,i,1) = PAOP_OUT_EXT(JL,JK,i,1) + incext(JL,JK)

             ! scattering portion
             NCscai = ai * incext(JL,JK)

             ! sum up weights for average (both albedo and asymmetry)
             NCsca    (JL,JK)   = NCsca    (JL,JK)   + NCscai 
             PAOP_OUT_G(JL,JK,i) = PAOP_OUT_G(JL,JK,i) + NCscai * gi

          END DO
          END DO
    
          ! Split extinction to separate contributions from constituents in modes.
          ! A volume mixing is assumed (in contrast to the explicit mixing in get_refr_ind). 
          if( wdep(i)%split .or. wdep(i)%insitu) then 
             DO JK=1,KLEV
              DO JL=KIDIA,KFDIA

                if (wdep(i)%split) then

                   ! The fine-mode contributions to the extinction
                   ! includes the contributions from particles
                   ! with (wet) diameters smaller than 1 micron.
                   ! For wet particles, only part of the accumulation mode
                   ! should be included, and the coarse mode should be excluded.
                   if (.not.coarse) then
                      ! Currently, the contribution of the accumulation mode
                      ! is approximated using weighting by volume scaling factor modfrac.
                      ! For extinction, area weighting would probably be more appropriate!!
                      if (imode .eq. 3 .or. imode .eq. 6 ) then
                         ! - convert number mean radius to volume mean radius (by cmedr2mmedr(imode))
                         ! - 1 micron diameter --> radius is 0.5 microns (rg is also in microns)
                         modfrac = 0.5 + 0.5 * ERF( log( 0.5 / (aop_in(JL,JK)%rg(imode) * cmedr2mmedr(imode)) ) / &
                              (sqrt(2.0) * lnsigma(imode)) )
                      else
                         ! Include full nucleation and Aitken mode contributions.
                         modfrac = 1.0
                      endif
                      PAOP_OUT_EXT(JL,JK,i,10) = PAOP_OUT_EXT(JL,JK,i,10) + modfrac * incext(JL,JK)
                   endif
                    ! total volume from so4/no3/bc/oc/ss/du in this mode (ATTENTION: DRY!!) 
                ! take no3 as so4
                totvoldry = aop_in(JL,JK)%so4(imode)/so4_density    + aop_in(JL,JK)%no3(imode)/so4_density + &
                            aop_in(JL,JK)%bc (imode)/carbon_density + aop_in(JL,JK)%oc (imode)/pom_density + &
                            aop_in(JL,JK)%soa (imode)/soa_density + &
                            aop_in(JL,JK)%ss (imode)/ss_density     + aop_in(JL,JK)%du (imode)/dust_density
                ! check whether there is some volume available 
                ! otherwise assign zeros to extinction increments
                if( totvoldry < tiny(totvoldry) ) then 
                   write(NULOUT,'("WARNING: no volume in mode (",i3,"). assigning zero extinctions")') imode
                   cycle
                end if
                ! H2-SO4 contribution
                PAOP_OUT_EXT(JL,JK,i,2) = PAOP_OUT_EXT(JL,JK,i,2) + incext(JL,JK) * (aop_in(JL,JK)%so4(imode)/so4_density   ) / totvoldry 
                ! BC contribution
                PAOP_OUT_EXT(JL,JK,i,3) = PAOP_OUT_EXT(JL,JK,i,3) + incext(JL,JK) * (aop_in(JL,JK)%bc (imode)/carbon_density) / totvoldry 
                ! POM contribution
                PAOP_OUT_EXT(JL,JK,i,4) = PAOP_OUT_EXT(JL,JK,i,4) + incext(JL,JK) * (aop_in(JL,JK)%oc (imode)/pom_density   ) / totvoldry 
                ! SOA contribution
                PAOP_OUT_EXT(JL,JK,i,5) = PAOP_OUT_EXT(JL,JK,i,5) + incext(JL,JK) * (aop_in(JL,JK)%soa (imode)/soa_density   ) / totvoldry 
                ! SS contribution
                PAOP_OUT_EXT(JL,JK,i,6) = PAOP_OUT_EXT(JL,JK,i,6) + incext(JL,JK) * (aop_in(JL,JK)%ss (imode)/ss_density    ) / totvoldry 
                ! DU contribution
                PAOP_OUT_EXT(JL,JK,i,7) = PAOP_OUT_EXT(JL,JK,i,7) + incext(JL,JK) * (aop_in(JL,JK)%du (imode)/dust_density  ) / totvoldry 
                ! NH4NO3 contribution
                PAOP_OUT_EXT(JL,JK,i,8) = PAOP_OUT_EXT(JL,JK,i,8) + incext(JL,JK) * (aop_in(JL,JK)%no3(imode)/so4_density   ) / totvoldry
                ! Fine-mode contributions for dust and sea salt 
                if (.not.coarse) then
                  PAOP_OUT_EXT(JL,JK,i,11) = PAOP_OUT_EXT(JL,JK,i,11) + incext(JL,JK) * (aop_in(JL,JK)%du (imode)/dust_density  ) / totvoldry * modfrac

                  PAOP_OUT_EXT(JL,JK,i,12) = PAOP_OUT_EXT(JL,JK,i,12) + incext(JL,JK) * (aop_in(JL,JK)%ss (imode)/ss_density  ) / totvoldry * modfrac
               endif
            endif

                   ! Water contribution: 
                   ! Get optical properties without water, the difference will be extinction due to 
                   ! water existence
                   ! - mis-use gi for this 
                   gi = 0.0 
                   call tm5m7_get_refr_idx( wdep(i), &
                        aop_in(JL,JK)%SO4(imode),&! + aop_in(JL,JK)%NO3(imode), &
                        aop_in(JL,JK)%BC (imode),  aop_in(JL,JK)%OC (imode), &
                        aop_in(JL,JK)%SOA(imode), &
                        aop_in(JL,JK)%SS (imode),  aop_in(JL,JK)%DU (imode), &
                        gi, imode, m_eff)
                 !if (statusomp ==1) then
                 !      write (NULERR,'("Problem with GET_REFR_IDX-NAN inputs:  ", 9(E16.4,2x))') &
                 !           aop_in(JL,JK)%SO4(imode), aop_in(JL,JK)%NO3(imode), aop_in(JL,JK)%h2o(imode), &
                 !           aop_in(JL,JK)%SOA(imode),aop_in(JL,JK)%BC(imode), aop_in(JL,JK)%OC(imode),&
                 !           aop_in(JL,JK)%SS(imode), aop_in(JL,JK)%DU(imode)
                 !      ! Assume edge case not caught. Set m_eff to default (1.0,1.0e-9), same as in ref_idx
                 !      ! (subroutine problems with inputs SO4 and WATER?)
                 !      m_eff= Cmplx(1.0,1.0e-9)
                 !      ! If all edge cases were handled correctly, that should be:
                 !      ! call ABOR1("error in tm5m7_optics_calculate_aop")
                 !  endif
                   
                   ! here we need the dry radius!!! 
                   !>>> TvN
                   ! 2*pi should be included (see comment above)
                   !xg = aop_in(JL,JK)%rgd(imode) / wdep(i)%wl
                   xg = twopi*aop_in(JL,JK)%rgd(imode) / wdep(i)%wl
                   !<<< TvN
                   ! TB
                   ! add extra safeguard for negative xg, unphysical but have a safeguard
                   ! added for refr_idx_nan problem 
                   if (xg .gt. 0.0) then
                      call TM5M7_OPTICS_GET(m_eff, xg, Cexti, ai, gi, cext_table, a_table, g_table, statusomp )  
                   else
                      Cexti=0.0
                      ai=1.0
                      gi=1.0
                   end if
                   
                   if (statusomp ==1) then
                      call ABOR1("error in tm5m7_optics_calculate_aop")
                   endif
                   
                   Cexti = Cexti*(wdep(i)%wl**2)
                   
                   if (wdep(i)%split) then
                      ! add difference to water subarray
                      PAOP_OUT_EXT(JL,JK,i,9) = PAOP_OUT_EXT(JL,JK,i,9) + (incext(JL,JK) - aop_in(JL,JK)%numdens(imode) * Cexti)
                   endif
                   
                   if (wdep(i)%insitu) then
                      ! Surface dry extinction and absorption: 
                      !>>> TvN
                      ! Remove cut off for the total values:
                      !modfrac = 0.5 + 0.5 * ERF( log( 5.0 / (aop_in(JL,JK)%rgd(imode) * cmedr2mmedr(imode)) ) / &
                      !     (sqrt(2.0) * lnsigma(imode)) )
                      ! extinction and absorption (extinction-scattering):
                      !PAOP_OUT_ADD(JL,JK,i,1) = PAOP_OUT_ADD(JL,JK,i,1) + modfrac * aop_in(JL,JK)%numdens(imode) * Cexti
                      !PAOP_OUT_ADD(JL,JK,i,2) = PAOP_OUT_ADD(JL,JK,i,2) + modfrac * aop_in(JL,JK)%numdens(imode) * Cexti * (1. - ai)
                      PAOP_OUT_ADD(JL,JK,i,1) = PAOP_OUT_ADD(JL,JK,i,1) + aop_in(JL,JK)%numdens(imode) * Cexti
                      PAOP_OUT_ADD(JL,JK,i,2) = PAOP_OUT_ADD(JL,JK,i,2) + aop_in(JL,JK)%numdens(imode) * Cexti * (1. - ai)
                      PAOP_OUT_ADD(JL,JK,i,3) = PAOP_OUT_ADD(JL,JK,i,3) + aop_in(JL,JK)%numdens(imode) * Cexti * ai * gi
                      ! Add fine-mode contributions
                      ! For dry aerosol, the fine-mode optical properties include 
                      ! the full accumulation mode, but not the coarse mode:
                      if (.not.coarse) then
                         PAOP_OUT_ADD(JL,JK,i,4) = PAOP_OUT_ADD(JL,JK,i,4) + aop_in(JL,JK)%numdens(imode) * Cexti
                         PAOP_OUT_ADD(JL,JK,i,5) = PAOP_OUT_ADD(JL,JK,i,5) + aop_in(JL,JK)%numdens(imode) * Cexti * (1. - ai)
                      endif
                      !<<< TvN
                   endif

             end do ! JL
             end do ! JK
          end if

       enddo   ! modes

       if (ecearth_units) then
          ! return extinction due to absorption
          ! and asymmetry factor multiplied by extinction due to scattering
          ! and convert um**2/m**3 into 1/m (as for extinction below)
          PAOP_OUT_A(KIDIA:KFDIA,1:KLEV,i) = ( PAOP_OUT_EXT(KIDIA:KFDIA,1:KLEV,i,1) - NCsca(KIDIA:KFDIA,1:KLEV) ) * 1.e-12
          PAOP_OUT_G(KIDIA:KFDIA,1:KLEV,i) = PAOP_OUT_G(KIDIA:KFDIA,1:KLEV,i) * 1.e-12
       else
          ! return single-scattering albedo and asymmetry factor
          ! take into account small extinction values
          where( PAOP_OUT_EXT(KIDIA:KFDIA,1:KLEV,i,1) > tiny(0.0) )
             PAOP_OUT_A(KIDIA:KFDIA,1:KLEV,i) = NCsca(KIDIA:KFDIA,1:KLEV) / PAOP_OUT_EXT(KIDIA:KFDIA,1:KLEV,i,1)
          elsewhere
             ! No extinction -> Fill 1.0 for albedo because it looks clean.
             PAOP_OUT_A(KIDIA:KFDIA,1:KLEV,i) = 1.0
          endwhere

          ! take into account small extinction values
          where( NCsca(KIDIA:KFDIA,1:KLEV) > tiny(0.0) )
             PAOP_OUT_G(KIDIA:KFDIA,1:KLEV,i) = PAOP_OUT_G(KIDIA:KFDIA,1:KLEV,i) / NCsca(KIDIA:KFDIA,1:KLEV)
          elsewhere
             ! No scattering at all -> Fill 1.0 for asymmetry because it looks clean.
             PAOP_OUT_G(KIDIA:KFDIA,1:KLEV,i) = 1.0
          endwhere
       endif

       ! convert um**2/m**3 into 1/m
       PAOP_OUT_EXT(KIDIA:KFDIA,1:KLEV,i,:) = PAOP_OUT_EXT(KIDIA:KFDIA,1:KLEV,i,:) * 1e-12

       if( wdep(i)%insitu) then 
          ! take into account small extinction values
          where( PAOP_OUT_ADD(KIDIA:KFDIA,1:KLEV,i,1) - PAOP_OUT_ADD(KIDIA:KFDIA,1:KLEV,i,2) > tiny(0.0) )
            PAOP_OUT_ADD(KIDIA:KFDIA,1:KLEV,i,3) = PAOP_OUT_ADD(KIDIA:KFDIA,1:KLEV,i,3) / (PAOP_OUT_ADD(KIDIA:KFDIA,1:KLEV,i,1)-PAOP_OUT_ADD(KIDIA:KFDIA,1:KLEV,i,2))
          elsewhere
          ! No scattering at all -> Fill 1.0 for asymmetry because it looks clean.
            PAOP_OUT_ADD(KIDIA:KFDIA,1:KLEV,i,3) = 1.0
          endwhere

          PAOP_OUT_ADD(KIDIA:KFDIA,1:KLEV,i,1) = PAOP_OUT_ADD(KIDIA:KFDIA,1:KLEV,i,1) * 1e-12 
          PAOP_OUT_ADD(KIDIA:KFDIA,1:KLEV,i,2) = PAOP_OUT_ADD(KIDIA:KFDIA,1:KLEV,i,2) * 1e-12 
          PAOP_OUT_ADD(KIDIA:KFDIA,1:KLEV,i,4) = PAOP_OUT_ADD(KIDIA:KFDIA,1:KLEV,i,4) * 1e-12 
          PAOP_OUT_ADD(KIDIA:KFDIA,1:KLEV,i,5) = PAOP_OUT_ADD(KIDIA:KFDIA,1:KLEV,i,5) * 1e-12
       endif

       if( wdep(i)%split .or. wdep(i)%insitu ) then
          deallocate( lnsigma )
       endif

       !=======================================================================
    enddo   ! loop over wavelengths

    

    Nullify(Cext_table)
    Nullify(a_table)
    Nullify(g_table)


    !do imode = 1,7
    !      print *, 'Radius,mode   :', imode, sum(aop_in(KIDIA:KFDIA,1:KLEV)%rg(imode))/((KFDIA-KIDIA+1)*KLEV)
    !      print *, 'numden,mode   :', imode, sum(aop_in(KIDIA:KFDIA,1:KLEV)%numdens(imode))/((KFDIA-KIDIA+1)*KLEV)
    !enddo

    ! deallocate( NCsca, incext )

    !do i = 1,nwl
    !   print *, 'AOD per grid box:', wdep(i)%wl, sum(PAOP_OUT_EXT(KIDIA:KFDIA,1:KLEV,i,1))/((KFDIA-KIDIA+1)*KLEV)
    !enddo



IF(LHOOK) CALL DR_HOOK('TM5M7_OPTICS_CALCULATE_AOP',1,ZHOOK_HANDLE)
END SUBROUTINE TM5M7_OPTICS_CALCULATE_AOP
