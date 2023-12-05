SUBROUTINE m7(KIDIA, KFDIA, KLON,   KLEV,           &  ! TM5  indices
              papp1,  prelhum, ptp1,           &  !   "   thermodynamics
              pso4g,  pelvoc, psvoc,  paerml, paernl,     &       !  M7   tracers
              prhop,  pww,    pm6rp,  pm6dry,  &  !   "   aerosol properties
              ptime )                             ! TM5  time step
  !
  !   ****m7* Aerosol model for the system so4,bc,oc,soa,ss,dust in 7 modes.
  !
  !   Authors:
  !   ---------
  !   E. Vignati, JRC/EI    (original source)					    2000
  !   P. Stier, MPI         (f90-version, changes, comments)		2001 
  !   E. Vignati, JRC/IES   (so2 is not required in this version)   2005
  !   T Bergman FMI         (modification for IFS) 2020
  !
  !   Purpose
  !   ---------
  !   Aerosol model for the system so4,bc,oc,ss,dust in 7 modes.
  !
  !   Externals
  !   ---------
  !
  !   *m7_averageproperties* 
  !       calculates the average mass for all modes and the particle 
  !       dry radius and density for the insoluble modes.
  !    
  !   *m7_equiz*   
  !       calculates the ambient radius of sulphate particles
  !
  !   *m7_equimix* 
  !       calculates the ambient radius of so4,bc,oc (dust) particles
  !
  !   *m7_equil* 
  !       calculates the ambient radius of so4,ss particles 
  !
  !   *m7_dgas*    
  !       calculates the sulfate condensation on existing particles
  !
  !   *m7_dnum*    
  !       calculates new gas phase sulfate and aerosol numbers and masses
  !       after condensation, nucleation and coagulation over one timestep
  !
  !   *m7_dconc*   
  !       repartitions aerosol number and mass between the
  !       the modes to account for condensational growth and the formation
  !       of an accumulation mode from the upper tail of the aitken mode and 
  !       of a coarse mode from the upper tail of the accumulation mode
  !

  !USE mo_aero_m7, ONLY: lsnucl, lscoag, lscond,         &
    USE TM5M7_DATA, ONLY:    nmod,   nss,    nsol,   naermod
USE PARKIND1 , ONLY : JPIM, JPRB
!  use tracer_data,   only : tracer_print
!  use GO,            only : gol, goErr, goPr, goBug
 ! use mo_aero,       only : nsoa      !RM

  IMPLICIT NONE 

  !--- Parameter list:
  !
  !  papp1      = atmospheric pressure at time t+1 [Pa]
  !  prelhum    = atmospheric relative humidity [% (0-1)]
  !  ptp1       = atmospheric temperature at time t+1 [K]
  !  pso4g      = mass of gas phase sulfate [molec. cm-3]
  !  paerml     = total aerosol mass for each compound 
  !               [molec. cm-3 for sulphate and ug m-3 for bc, oc, ss, and dust]
  !  paernl     = aerosol number for each mode [cm-3]
  !  prhop      = mean mode particle density [g cm-3]
  !  pm6rp      = mean mode actual radius (wet radius for soluble modes 
  !               and dry radius for insoluble modes) [cm]
  !  pm6dry     = dry radius for soluble modes [cm]
  !  pww        = aerosol water content for each mode [kg(water) m-3(air)]
  !
  !--- Local variables:
  !
  !  zttn       = average mass for single compound in each mode 
  !               [in molec. for sulphate and in ug for bc, oc, ss, and dust]
  !  zhplus     = number of h+ in mole [???] (kg water)-1
  !  zso4_x     = mass of sulphate condensed on insoluble mode x [molec. cm-3]
  !               (calculated in dgas used in concoag)


  ! Parameters:

  INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA, KFDIA, KLON, KLEV

  REAL(KIND=JPRB),INTENT(IN)    :: ptime

  REAL(KIND=JPRB),INTENT(IN)    :: prelhum(KLON,KLEV),        papp1(KLON,KLEV),           &
             ptp1(KLON,KLEV)  
  REAL(KIND=JPRB),INTENT(INOUT)    ::           pso4g(KLON,KLEV),  pelvoc(KLON,KLEV),  psvoc(KLON,KLEV) !RM
 
  REAL(KIND=JPRB),INTENT(INOUT)    :: paerml(KLON,KLEV,naermod), paernl(KLON,KLEV,nmod),     &
             pm6rp(KLON,KLEV,nmod),     pm6dry(KLON,KLEV,nsol),     &
             prhop(KLON,KLEV,nmod),     pww(KLON,KLEV,nmod)

  ! Local variables:


  REAL(KIND=JPRB)    :: zso4_5(KLON,KLEV),         zso4_6(KLON,KLEV),          &
             zso4_7(KLON,KLEV)

  REAL(KIND=JPRB)    :: zhplus(KLON,KLEV,nss)

  REAL(KIND=JPRB)    :: zttn(KLON,KLEV,naermod)

  !
  !--- 0) Initialisations: -------------------------------------------------
  !
  zhplus(:,:,:) = 0.
  pm6dry(:,:,:) = 0. 
  pm6rp(:,:,:)  = 0.
  zttn(:,:,:)   = 0.
  prhop(:,:,:)  = 0.
  pww(:,:,:)    = 0. 
  zso4_5(:,:)   = 0.
  zso4_6(:,:)   = 0.
  zso4_7(:,:)   = 0.

  

!!$  !
!!$  !--- 1) Calculation of particle properties under ambient conditions: -----
!!$  !
!!$  !--- 1.1) Calculate mean particle mass for all modes 
!!$  !         and dry radius and density for the insoluble modes.
!!$  !
!!$  CALL m7_averageproperties(KIDIA, KFDIA, KLON, KLEV, paernl, paerml, zttn, pm6rp, prhop)
!!$  !
!!$  !--- 1.2) Calculate ambient count median radii and density 
!!$  !         for lognormal distribution of particles.
!!$  !
!!$  !         Sulfate particles:
!!$  !
!!$  CALL m7_equiz(KIDIA, KFDIA,  KLON, KLEV,   &
!!$                papp1,   zttn,  ptp1,   &
!!$                prelhum, pm6rp, pm6dry, &
!!$                prhop,   pww,   paernl  )
!!$  !         
!!$  !         Mixed particles with sulfate, b/o carbon and dust: 
!!$  !
!!$  CALL m7_equimix(KIDIA, KFDIA,  KLON, KLEV,   &
!!$                  papp1,   zttn,  ptp1,   &
!!$                  prelhum, pm6rp, pm6dry, &
!!$                  prhop,   pww,   paernl  )
!!$  !
!!$  !         Accumulation and coarse mode particles in presence of
!!$  !         sea salt particles:
!!$  !
!!$  CALL m7_equil(KIDIA, KFDIA, KLON,  KLEV,   prelhum, paerml, paernl, &
!!$                pm6rp,  pm6dry, zhplus, pww,     prhop,  ptp1    )
!!$  !
!!$  !
!!$  !--- 2) Calculate changes in aerosol mass and gas phase sulfate ----------
!!$  !       due to sulfate condensation:
!!$  !       No change in particle mass/number relationships.
!!$  !
!!$  IF (lscond) CALL m7_dgas(KIDIA, KFDIA, KLON,  KLEV,  pso4g,  paerml, paernl, &
!!$                           ptp1,   papp1,  pm6rp, zso4_5, zso4_6, zso4_7, &
!!$                           ptime) 
!!$
!!$  !--- 2b) Calculate changes in aerosol mass ----------
!!$  !       due to organic condensation:
!!$  !       No change in particle mass/number relationships.
!!$  !
!!$
!!$  IF (nsoa .GT. 0 .AND. lscond) THEN
!!$    CALL m7_dgas_org(KIDIA, KFDIA, KLON,  KLEV,  pelvoc, psvoc, paerml, paernl, &
!!$                             ptp1,   papp1,  pm6rp, &
!!$                             ptime) 
!!$  END IF
!!$  !
!!$  !
!!$  !--- 3) Calculate change in particle number concentrations ---------------
!!$  !       due to nucleation and coagulation:
!!$  !       Change particle mass/number relationships.
!!$  !
!!$  ! JadB: Removed "If (lsnucl .OR. lscoag)".
!!$  ! If only lscond is set, the m7_dnum is required for storing sulfuric acid on insoluble aerosols (making the soluble).
!!$  ! Without m7_dnum, the sulfuric acid condensed on insoluble particles is turned into void.
!!$  ! IF (lsnucl.OR.lscoag) CALL m7_dnum(KIDIA, KFDIA, KLON,   KLEV,          &
!!$  !                                    pso4g,  paerml,  paernl, ptp1,  &
!!$  !                                    papp1,  prelhum, pm6rp,  prhop, &
!!$  !                                    zso4_5, zso4_6,  zso4_7, ptime   )
!!$
!!$  CALL m7_dnum(KIDIA, KFDIA, KLON,   KLEV,          &
!!$               pso4g,  pelvoc,  paerml,  paernl, ptp1,  &
!!$               papp1,  prelhum, pm6rp,  prhop, &
!!$               zso4_5, zso4_6,  zso4_7, ptime   )
!!$  !
!!$  !
!!$  !--- 4) Recalculation of particle properties under ambient conditions: ---
!!$  !
!!$  !--- 4.1) Recalculate mean masses for all modes 
!!$  !         and dry radius and density for the insoluble modes.
!!$  !
!!$  CALL m7_averageproperties(KIDIA, KFDIA, KLON, KLEV, paernl, paerml, zttn, pm6rp, prhop)
!!$  !
!!$  !--- 4.2) Calculate ambient count median radii and density 
!!$  !         for lognormal distribution of particles.
!!$  !
!!$  !         Sulfate particles:
!!$  !
!!$  CALL m7_equiz(KIDIA, KFDIA,  KLON, KLEV,   &
!!$                papp1,   zttn,  ptp1,   &
!!$                prelhum, pm6rp, pm6dry, &
!!$                prhop,   pww,   paernl  )
!!$  !
!!$  !         Mixed particles with sulfate, b/o carbon and dust:
!!$  !
!!$  CALL m7_equimix(KIDIA, KFDIA,  KLON, KLEV,   &
!!$                  papp1,   zttn,  ptp1,   &
!!$                  prelhum, pm6rp, pm6dry, &
!!$                  prhop,   pww,   paernl  )
!!$  !
!!$  !         Accumulation and coarse mode particles in presence of
!!$  !         sea salt particles:  
!!$  ! 
!!$  CALL m7_equil(KIDIA, KFDIA, KLON,  KLEV,   prelhum, paerml, paernl,  &
!!$                pm6rp,  pm6dry, zhplus, pww,     prhop,  ptp1     )
!!$  !
!!$  !--- 5) Repartitition particles among the modes: -------------------------
!!$  !
!!$  IF (lscond.OR.lscoag) THEN
!!$     
!!$     CALL m7_dconc(KIDIA, KFDIA, KLON, KLEV, paerml, paernl, pm6dry)
!!$     
!!$  END IF
!!$  !
!!$  !--- 6) Recalculation of particle properties under ambient conditions: ---
!!$  !
!!$  !--- 6.1) Calculate mean particle mass for all modes 
!!$  !         and dry radius and density for the insoluble modes:
!!$  !
!!$  CALL m7_averageproperties(KIDIA, KFDIA, KLON, KLEV, paernl, paerml, zttn, pm6rp, prhop)
!!$  !
!!$  !--- 6.2) Calculate ambient count median radii and density 
!!$  !         for lognormal distribution of particles.
!!$  !
!!$  !         Sulfate particles:
!!$  !
!!$  CALL m7_equiz(KIDIA, KFDIA,  KLON, KLEV,   &
!!$                papp1,   zttn,  ptp1,   &
!!$                prelhum, pm6rp, pm6dry, &
!!$                prhop,   pww,   paernl  )
!!$  !
!!$  !         Mixed particles with sulfate, b/o carbon and dust: 
!!$  !
!!$  CALL m7_equimix(KIDIA, KFDIA,  KLON, KLEV,   &
!!$                  papp1,   zttn,  ptp1,   &
!!$                  prelhum, pm6rp, pm6dry, &
!!$                  prhop,   pww,   paernl  )
!!$  !
!!$  !         Accumulation and coarse mode particles in presence of
!!$  !         sea salt particles:
!!$  !
!!$  CALL m7_equil(KIDIA, KFDIA, KLON,  KLEV, prelhum, paerml, paernl, &
!!$                pm6rp,  pm6dry, zhplus, pww,   prhop,  ptp1    )
!!$

!  write(*,*) 'ou2', 'h2so4= ', pso4g(2100,1),'num1= ',paernl(2100,1,1)
!  write(*,*) 'ou2', 'BCsol= ', paerml(2100,1,5),'BCins= ', paerml(2100,1,8)
!  write(*,*) 'ou2', 'POsol= ', paerml(2100,1,9),'POins= ', paerml(2100,1,12)
!  write(*,*) 'ou2', 'num2= ', paernl(2100,1,2), 'num5= ', paernl(2100,1,5)
!  write(*,*) 'ou2', 'BCtot= ', paerml(2100,1,5)+ paerml(2100,1,8)+paerml(2100,1,6)+paerml(2100,1,7)
!  write(*,*) 'ou2', 'cond5= ', zso4_5(2100,1), '6= ',zso4_6(2100,1), '7= ',zso4_7(2100,1)
  !
END SUBROUTINE m7
