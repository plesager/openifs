SUBROUTINE hamm7_init(YGFL, YRRIP, CHEM_SCHEME)

! ╭────────────────────────────────────────────────────────────────────────────╮
! │                                                      (updated 14-MAY-2024) │
! │ Purpose :                                                                  │
! │ -------                                                                    │
! │    init routine HAM-M7 aerosol in OpenIFS                                  │
! │                                                                            │
! │                                                                            │
! │ Interface :                                                                │
! │ ---------                                                                  │
! │   *hamm7_init* is called from *CNT4*.                                      │
! │   REMARK: the code assumes that TM5M7_INIT has already been called!        │
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
! │     2020-11-11   - Thomas Kuehn (FMI/UEF)                                  │
! │                                                                            │
! │ Modifications :                                                            │
! │ -------------                                                              │
! │     ?          - Eemeli Holopainen (FMI) : ?                               │
! │     May.  2024 - R. Checa-Garcia (KNMI)  : revision for CY48r1             │
! │                                                                            │
! ╰────────────────────────────────────────────────────────────────────────────╯

!----------------------------------------------------------------------
!*       0.0   USE STATEMENTS
!              ---------------

! IFS/OPENIFS modules ---------------------------------------------------------
USE PARKIND1, ONLY : JPRB
USE YOMHOOK,  ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMRIP,   ONLY : TRIP
USE YOM_YGFL, ONLY : TYPE_GFLD   ! Gives type for YGFL
USE YOMMP0,   ONLY : MYPROC
USE YOMLUN,   ONLY : NULOUT

! --- M7 modules --------------------------------------------------------------
USE MO_TIME_CONTROL, ONLY: init_mo_time_control
USE MO_SPECIES,      ONLY: speclist            ! tracer species in HAM
USE MO_HAM, ONLY:     &
     sizeclass,       & ! aerosol classes in HAM
     nclass,          & ! number of aerosol classes in HAM
     aerocomp,        & ! aerosol compounds by size class in HAM
     naerocomp,       & ! amount of aerosol mass tracers in HAM
     subm_gasspec,    & ! gas phase species in HAM
     subm_ngasspec      ! number of gas phase species in HAM

USE OIFS_to_HAM, ONLY: init_ind_oifs_ham, &  ! init index list for OIFS and HAM
   &                   ind_oifs_ham          ! index list type
     !ind_class_OIFS,        & ! index list for aerosol sizeclasses from OIFS
     !ind_class_HAM,         & ! index list for aerosol sizeclasses from HAM
     !ind_mass_OIFS,        & ! index list for aerosol masses from OIFS
     !ind_mass_HAM,         & ! index list for aerosol masses from HAM
     !ind_gas_OIFS,         & ! index list for gases from OIFS
     !ind_gas_HAM,          &  ! index list for gases from HAM
     !ind_cloud_OIFS,       & ! index list for cloud variables from OIFS
     !ind_cloud_HAM           ! index list for cloud variables from HAM

USE MO_TRACDEF,     ONLY:  GAS, AEROSOL,         & ! species type identifiers
  &                        GAS_OR_AEROSOL, ntrac, trlist

USE MO_HAM_M7_TRAC, ONLY:  idt_cdnc_ham,         & !index for HAM CDNC
  &                        idt_icnc_ham            !index for HAM ICNC

USE MO_HAM_INIT, ONLY:     &
     start_ham,            &
     ham_initialize,       &
     ham_define_tracer       !eehol: added ham_define_tracer

USE MO_ADVECTION, ONLY:    & !eehol: added for advection initialization
     iadvec, tpcore

USE MO_SUBMODEL, &
     ONLY:     & !eehol: added mo_submodel routines 
     lham, id_ham

USE MO_SPECIES, &
     ONLY: speclist, init_splist !eehol: added initialization for species list

USE MO_HAM_SOA, &
     ONLY: soaprop

!eehol: activation initialization
USE MO_ACTIV,   &
     ONLY: activ_initialize, idt_cdnc, idt_icnc

USE MO_PARAM_SWITCHES, &
     ONLY: ncd_activ, nactivpdf, lcdnc_progn


IMPLICIT NONE
TYPE(TYPE_GFLD), INTENT(IN)    :: YGFL
TYPE(TRIP),      INTENT(IN)    :: YRRIP
CHARACTER(len=20), INTENT(IN)  :: CHEM_SCHEME

!----------------------------------------------------------------------
!*       0.5   LOCAL VARIABLES
!              ---------------
INTEGER ::                                   &         ! looping indices
 &  j_yaero, j_ychem,                        &         ! IFS
 &  j_class, j_mass, j_spec, j_gas, j_cloud, &         ! HAM
 &  kt, znclass, znaerocomp, zsubm_ngasspec, zcloudind ! eehol: indices for OIFS to HAM

CHARACTER(len=64) :: int_str, int_str_ham              ! eehol: integer as string
LOGICAL, PARAMETER :: LLDEBUG=.FALSE.                  ! Debug flag

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE                      ! return status

!----------------------------------------------------------------------
!*       1.0   EXECUTABLE CODE
!              ---------------
!ASSOCIATE(NACTAERO=>YGFL%NACTAERO, LAERCHEM=>YGFL%LAERCHEM)

IF (LHOOK) CALL DR_HOOK('HAMM7_INIT',0,ZHOOK_HANDLE)

! getting the correct time step:
CALL init_mo_time_control(YRRIP)

!eehol: set advection scheme
iadvec = tpcore !comes from ECHAM mo_control.f90
!eehol: activation initialization 
nactivpdf = 1 !eehol: using PDF to calculate updraft. Hardcoded for now.. need to check this later (add setphys to oifs?)
ncd_activ = 2 !eehol: Abdul-Razzak and Ghan activation scheme. Hardcoded for now.. need to check this later (add setphys to oifs?)
lcdnc_progn = .TRUE.

!eehol: set submodel parameters and flags
CALL init_splist !eehol: added init for splist

CALL start_ham

id_ham=1
CALL ham_define_tracer

CALL ham_initialize
! ham_init_memory calls:
! ham_nucl_initialize, if nsnucl+nonucl > 0
!CALL ham_init_memory !eehol: this is not needed now

CALL activ_initialize

!eehol: set cdnc and icnc indices for HAM
idt_cdnc = idt_cdnc_ham
idt_icnc = idt_icnc_ham

!<--eehol: initialize index lists for OIFS and HAM
znclass = nclass
znaerocomp = naerocomp
zsubm_ngasspec = subm_ngasspec

! RCHG -> for consistency I would define zcloudind in the MO_HAM like the similar
!         definitions rather than here.
zcloudind = 2                                            !eehol CDNC and ICNC

CALL init_ind_oifs_ham(znclass, znaerocomp, zsubm_ngasspec, zcloudind) ! RCHG -> allocates and initialize to zero.
!-->eehol

! assigning tracer indices as found in YGFL to the
! corresponding metadata fields in HAM
ASSOCIATE(&
     NAERO=>YGFL%NAERO, YAERO=>YGFL%YAERO, & ! aerosol tracer meta-data
     NCHEM=>YGFL%NCHEM, YCHEM=>YGFL%YCHEM, LAERCHEM=>YGFL%LAERCHEM & ! chemistry tracer meta-data
)

! The following code tries to identify the dynamically generated tracers
! in OpenIFS to the independently generated metadata generated in HAM.
! The point of this procedure is to be able to directly access the
! tracer data stored in PCEN and PTENC in the HAM routines without
! having to copy the data into new fields.
! The naming conventions in IFS and HAM are somewhat different, so
! tracer naming is somewhat hard-coded. As a first instance to troubleshooting,
! check the logfile 'fort.2001' in the run directory and see which tracers 
! are not recognised by HAM. Then check the file 
! AC-experiments/ctrl/Table/bins_hamm7ver1.csv and compare the names you get 
! in the files.

! aerosol tracers
LABEL_IFS_AERO: DO j_yaero = 1,NAERO

   ! looping over the number tracers in HAM
   LABEL_CLASS: DO j_class = 1,nclass

      ! sizeclass names in HAM-M7 are two-letter (xy) combinations: 
      ! x: N,K,A,C = nucleation, Aitken, accumulation, and coarse, respectively, while
      ! y: S,I = soluble, insoluble, respectively.
      ! In the IFS-csv table, these combinations are suffixed with '_N'

      IF (TRIM(YAERO(j_yaero)%CNAME) == TRIM(sizeclass(j_class)%shortname)//'_N') THEN

         ! In case of a match, we set the tracer index in the HAM meta data to the
         ! according IFS-index and write a note into the logfile
         !sizeclass(j_class)%idt_no = j_yaero

         kt = sizeclass(j_class)%idt_no
         ind_oifs_ham%ind_class_OIFS(j_class) = j_yaero !eehol: take indices for sizeclass to a vector for OIFS
         ind_oifs_ham%ind_class_HAM(j_class) = kt !eehol: take indices for sizeclass to a vector for HAM
         IF (LLDEBUG .AND. MYPROC == 1) THEN
            WRITE(int_str,*) j_yaero
            WRITE(int_str_ham,*) kt !sizeclass(j_class)%idt_no
            WRITE(2000+MYPROC,'(a)') 'sizeclass: '//TRIM(YAERO(j_yaero)%CNAME)
            WRITE(3000+MYPROC,'(a)') 'OIFS size: '//TRIM(YAERO(j_yaero)%CNAME)
            WRITE(3000+MYPROC,'(a)') 'j_yaero ='//TRIM(int_str)
            WRITE(3000+MYPROC,'(a)') 'HAM size: '//TRIM(trlist%ti(sizeclass(j_class)%idt_no)%fullname)
            WRITE(3000+MYPROC,'(a)') 'sizeclass idt_no ='//TRIM(int_str_ham)
         END IF

         ! once the tracer is identified, we can check the next one
         CYCLE LABEL_IFS_AERO
      END IF
   END DO LABEL_CLASS

   ! looping over the aerosol compound tracers in HAM
   LABEL_MASS: DO j_mass = 1,naerocomp
      ! aerosol compound names in HAM-M7 are combinations of the chemical species
      ! and the sizeclass. The format is zzz_xy, where x and y are as in the sizeclass
      ! and zzz is a 2- or 3-letter abbreviation
      ! standard names are OC, BC, SO4, DU, and SS for organic carbon, black carbon
      ! sulfate, mineral dust and seasalt, respectively. Other species may be
      ! included as well, depending on the HAM setup.

      j_class = aerocomp(j_mass)%iclass  ! index to size class
      j_spec  = aerocomp(j_mass)%spid    ! index to species

      IF (TRIM(YAERO(j_yaero)%CNAME) == (TRIM(speclist(j_spec)%shortname)//'_'//TRIM(sizeclass(j_class)%shortname))) THEN

         ! In case of a match, we set the tracer index in the HAM meta data to the
         ! according IFS-index and write a note into the logfile
         !aerocomp(j_mass)%idt = j_yaero
         !ind_mass_OIFS(j_mass) = j_yaero !eehol: take indices for masses to a vector
         !ind_mass_HAM(j_mass) = aerocomp(j_mass)%idt !eehol: take indices for masses to a vector
         kt = aerocomp(j_mass)%idt
         ind_oifs_ham%ind_mass_OIFS(j_mass) = j_yaero !eehol: take indices for masses to a vector for OIFS
         ind_oifs_ham%ind_mass_HAM(j_mass) = kt !eehol: take indices for masses to a vector for HAM
         IF (LLDEBUG .AND. MYPROC == 1) THEN
            WRITE(int_str,*) j_yaero
            WRITE(int_str_ham,*) kt
            WRITE(2000+MYPROC,'(a)') 'mass: '//TRIM(YAERO(j_yaero)%CNAME)
            WRITE(6000+MYPROC,'(a)') 'OIFS mass: '//TRIM(YAERO(j_yaero)%CNAME)
            WRITE(6000+MYPROC,'(a)') 'j_yaero ='//TRIM(int_str)
            WRITE(6000+MYPROC,'(a)') 'HAM mass: '//TRIM(trlist%ti(aerocomp(j_mass)%idt)%fullname)
            WRITE(6000+MYPROC,'(a)') 'mass idt_no ='//TRIM(int_str_ham)
         END IF
         
         ! once the tracer is identified we can check the next one
         CYCLE LABEL_IFS_AERO

      ELSE IF ( (TRIM(YAERO(j_yaero)%CNAME) == 'POM_'//TRIM(sizeclass(j_class)%shortname)) .AND. (TRIM(speclist(j_spec)%shortname) == 'OC') ) THEN

         ! In case of a match, we set the tracer index in the HAM meta data to the
         ! according IFS-index and write a note into the logfile
         !aerocomp(j_mass)%idt = j_yaero

         kt = aerocomp(j_mass)%idt
         ind_oifs_ham%ind_mass_OIFS(j_mass) = j_yaero !eehol: take indices for masses to a vector for OIFS
         ind_oifs_ham%ind_mass_HAM(j_mass) = kt !eehol: take indices for masses to a vector for HAM

         IF (LLDEBUG .AND. MYPROC == 1) THEN
            WRITE(int_str,*) j_yaero
            WRITE(int_str_ham,*) kt
            WRITE(2000+MYPROC,'(a)') 'mass: '//TRIM(YAERO(j_yaero)%CNAME)
            WRITE(6000+MYPROC,'(a)') 'OIFS mass: '//TRIM(YAERO(j_yaero)%CNAME)
            WRITE(6000+MYPROC,'(a)') 'j_yaero ='//TRIM(int_str)
            WRITE(6000+MYPROC,'(a)') 'HAM mass: '//TRIM(trlist%ti(aerocomp(j_mass)%idt)%fullname)
            WRITE(6000+MYPROC,'(a)') 'mass idt_no ='//TRIM(int_str_ham)
         END IF

         ! once the tracer is identified we can check the next one
         CYCLE LABEL_IFS_AERO
      END IF

   END DO LABEL_MASS

   !eehol: looping over cloud compounds
   LABEL_CLOUD: DO j_cloud = idt_cdnc,idt_icnc !eehol: only CDNC and ICNC
      
      IF (TRIM(YAERO(j_yaero)%CNAME) == TRIM(trlist%ti(j_cloud)%fullname)) THEN

         ! In case of a match, we set the tracer index in the HAM meta data to the
         ! according IFS-index and write a note into the logfile
         !sizeclass(j_class)%idt_no = j_yaero
         kt = j_cloud
         ind_oifs_ham%ind_cloud_OIFS(j_cloud-idt_cdnc+1) = j_yaero !eehol: take indices for sizeclass to a vector for OIFS
         ind_oifs_ham%ind_cloud_HAM(j_cloud-idt_cdnc+1) = kt !eehol: take indices for sizeclass to a vector for HAM
         IF (LLDEBUG .AND. MYPROC == 1) THEN
            WRITE(int_str,*) j_yaero
            WRITE(int_str_ham,*) kt !sizeclass(j_class)%idt_no
            WRITE(2000+MYPROC,'(a)') 'cloud class: '//TRIM(YAERO(j_yaero)%CNAME)
            WRITE(3000+MYPROC,'(a)') 'OIFS size: '//TRIM(YAERO(j_yaero)%CNAME)
            WRITE(3000+MYPROC,'(a)') 'j_yaero ='//TRIM(int_str)
            WRITE(3000+MYPROC,'(a)') 'HAM size: '//TRIM(trlist%ti(kt)%fullname)
            WRITE(3000+MYPROC,'(a)') 'sizeclass idt_no ='//TRIM(int_str_ham)
         END IF

         ! once the tracer is identified, we can check the next one
         CYCLE LABEL_IFS_AERO
      END IF

   END DO LABEL_CLOUD

   ! If we end up here, the tracer is not a sizeclass or mass
   IF (LLDEBUG .AND. MYPROC == 1) THEN
      WRITE(2000+MYPROC,'(a)') 'neither mass nor sizeclass in HAM: '//TRIM(YAERO(j_yaero)%CNAME)
   END IF

END DO LABEL_IFS_AERO

! chemistry tracers
! chemistry tracers in HAM do not have separate metadata, because there is only one 
! tracer per species. Therefore we have to check the HAM species metadata and
! make sure that the metadata actually describes a gas.
IF (LAERCHEM.AND.(TRIM(CHEM_SCHEME)=="tm5".OR.TRIM(CHEM_SCHEME)=="SimChem"))THEN
   LABEL_IFS_CHEM: DO j_ychem = 1,NCHEM

   ! looping over the gas phase tracers in HAM
   LABEL_GAS: DO j_gas = 1,subm_ngasspec
      j_spec = subm_gasspec(j_gas)

      ! If the species does not exist in gas form in HAM, we don't consider it further
      IF ((speclist(j_spec)%nphase == GAS) .OR. (speclist(j_spec)%nphase == GAS_OR_AEROSOL)) THEN
         IF (TRIM(YCHEM(j_ychem)%CNAME) == TRIM(speclist(j_spec)%shortname)) THEN
            ! In case of a match, we set the tracer index in the HAM meta data to the
            ! according IFS-index and write a note into the logfile
            !speclist(j_spec)%idt = j_ychem
            kt = speclist(j_spec)%idt
            ind_oifs_ham%ind_gas_OIFS(j_gas) = j_ychem !eehol: take indices for gases to a vector for OIFS
            ind_oifs_ham%ind_gas_HAM(j_gas) = kt !eehol: take indices for gases to a vector for HAM

            ! turn wetdep off for gases if LAERCHEM=.true.
            trlist%ti(kt)%nwetdep=0
            IF (LLDEBUG .AND. MYPROC == 1) THEN
               WRITE(int_str,*) j_ychem
               WRITE(int_str_ham,*) kt
               WRITE(2000+MYPROC,'(a)') 'gas: '//TRIM(YCHEM(j_ychem)%CNAME)
               WRITE(7000+MYPROC,'(a)') 'OIFS gas: '//TRIM(YCHEM(j_ychem)%CNAME)
               WRITE(7000+MYPROC,'(a)') 'j_ychem ='//TRIM(int_str)
               WRITE(7000+MYPROC,'(a)') 'HAM gas: '//TRIM(trlist%ti(speclist(j_spec)%idt)%fullname)
               WRITE(7000+MYPROC,'(a)') 'gas idt_no ='//TRIM(int_str_ham)
               ! RCHG -> line below printed non string characters. Probably issue with dimensions
               !         strings, so commented. FIXME
               !WRITE(7000+MYPROC,'(a)') 'wetdep',trlist%ti(kt)
            END IF

            ! once the tracer is identified we can check the next one
            CYCLE LABEL_IFS_CHEM

         END IF
         ! sulfuric acid is special, because in HAM it is tagged as H2SO4, while in IFS it
         ! is called SO4 ('tm5' case) 'H2SO4' or ('SimChem' case )
         IF ( (TRIM(speclist(j_spec)%shortname) == 'H2SO4') .AND. &
              & ( ( TRIM(YCHEM(j_ychem)%CNAME) == 'H2SO4' .AND. TRIM(CHEM_SCHEME)=="SimChem" ) &
              & .OR. ( TRIM(YCHEM(j_ychem)%CNAME) == 'SO4' .AND. TRIM(CHEM_SCHEME)=="tm5" )  ) ) THEN
            ! In case of a match, we set the tracer index in the HAM meta data to the
            ! according IFS-index and write a note into the logfile
            kt = speclist(j_spec)%idt
            ind_oifs_ham%ind_gas_OIFS(j_gas) = j_ychem !eehol: take indices for gases to a vector for OIFS
            ind_oifs_ham%ind_gas_HAM(j_gas) = kt !eehol: take indices for gases to a vector for HAM
            ! turn wetdep off for SO4 if LAERCHEM=.true.
            trlist%ti(kt)%nwetdep=0

            IF (LLDEBUG .AND. MYPROC == 1) THEN
               WRITE(int_str,*) j_ychem
               WRITE(int_str_ham,*) kt
               WRITE(2000+MYPROC,'(a)') 'gas: '//TRIM(YCHEM(j_ychem)%CNAME)
               WRITE(7000+MYPROC,'(a)') 'OIFS gas: '//TRIM(YCHEM(j_ychem)%CNAME)
               WRITE(7000+MYPROC,'(a)') 'j_ychem ='//TRIM(int_str)
               WRITE(7000+MYPROC,'(a)') 'HAM gas: '//TRIM(trlist%ti(speclist(j_spec)%idt)%fullname)
               WRITE(7000+MYPROC,'(a)') 'gas idt_no ='//TRIM(int_str_ham)
               WRITE(7000+MYPROC,'(a)') 'wetdep=',trlist%ti(kt)%nwetdep
            END IF

            ! once the tracer is identified we can check the next one
            CYCLE LABEL_IFS_CHEM
         END IF
      END IF
   END DO LABEL_GAS

   ! If we end up here, the tracer is not a sizeclass or mass
   IF (LLDEBUG .AND. MYPROC == 1) THEN
      WRITE(2000+MYPROC,'(a)') 'not a gas in HAM: '//TRIM(YCHEM(j_ychem)%CNAME)
   END IF

   END DO LABEL_IFS_CHEM
 
!ELSE IF (LAERCHEM.AND.TRIM(CHEM_SCHEME)=="SimChem") THEN
!   Write(9191,*)'simple SO4 in development'
!
!   LABEL_IFS_CHEM_SO4: DO j_ychem = 1,NAERO
!
!   ! looping over the gas phase tracers in HAM
!   LABEL_GAS_SO4: DO j_gas = 1,subm_ngasspec
!      j_spec = subm_gasspec(j_gas)
!      ! If the species does not exist in gas form in HAM, we don't consider it further
!      IF ((speclist(j_spec)%nphase == GAS) .OR. (speclist(j_spec)%nphase == GAS_OR_AEROSOL)) THEN
!         IF(MYPROC==1)THEN
!            write(4001,*)TRIM(YAERO(j_ychem)%CNAME),TRIM(speclist(j_spec)%shortname)
!         END IF
!         IF (TRIM(YAERO(j_ychem)%CNAME) == TRIM(speclist(j_spec)%shortname)) THEN
!            ! In case of a match, we set the tracer index in the HAM meta data to the
!            ! according IFS-index and write a note into the logfile
!            !speclist(j_spec)%idt = j_ychem
!            kt = speclist(j_spec)%idt
!            ind_oifs_ham%ind_gas_OIFS(j_gas) = j_ychem !eehol: take indices for gases to a vector for OIFS
!            ind_oifs_ham%ind_gas_HAM(j_gas) = kt !eehol: take indices for gases to a vector for HAM
!            ! turn wetdep on for gases if LAERCHEM=.false.
!            trlist%ti(kt)%nwetdep=1
!            IF (LLDEBUG .AND. MYPROC == 1) THEN
!               WRITE(int_str,*) j_ychem
!               WRITE(int_str_ham,*) kt
!               WRITE(2000+MYPROC,'(a)') 'gas: '//TRIM(YCHEM(j_ychem)%CNAME)
!               WRITE(7000+MYPROC,'(a)') 'OIFS gas: '//TRIM(YCHEM(j_ychem)%CNAME)
!               WRITE(7000+MYPROC,'(a)') 'j_ychem ='//TRIM(int_str)
!               WRITE(7000+MYPROC,'(a)') 'HAM gas: '//TRIM(trlist%ti(speclist(j_spec)%idt)%fullname)
!               WRITE(7000+MYPROC,'(a)') 'gas idt_no ='//TRIM(int_str_ham)
!            END IF
!            ! once the tracer is identified we can check the next one
!            CYCLE LABEL_IFS_CHEM_SO4
!         END IF
!
!         ! sulfuric acid is special, because in HAM it is tagges as H2SO4, while in IFS it
!         ! is called SO4 -- hard-coding for now
!         !IF ( (TRIM(YAERO(j_ychem)%CNAME) == 'SO4_gas') .AND. (TRIM(speclist(j_spec)%shortname) == 'H2SO4') ) THEN
!         !WRITE(*,*)"j_spec:",j_spec
!         !WRITE(*,*)"j_ychem:",j_ychem
!         !WRITE(*,*)"YAERO(j_ychem)%CNAME:",YAERO(j_ychem)%CNAME
!         !WRITE(*,*)"speclist(j_spec)%shortname:",speclist(j_spec)%shortname
!         !IF ( (TRIM(YAERO(j_ychem)%CNAME) == 'SO4_gas') .AND. (TRIM(speclist(j_spec)%shortname) == 'H2SO4') ) THEN
!         IF ( (TRIM(YAERO(j_ychem)%CNAME) == 'SO4') .AND. (TRIM(speclist(j_spec)%shortname) == 'H2SO4') ) THEN
!            ! In case of a match, we set the tracer index in the HAM meta data to the
!            ! according IFS-index and write a note into the logfile
!            kt = speclist(j_spec)%idt
!            ind_oifs_ham%ind_gas_OIFS(j_gas) = j_ychem !eehol: take indices for gases to a vector for OIFS
!            ind_oifs_ham%ind_gas_HAM(j_gas) = kt !eehol: take indices for gases to a vector for HAM
!            ! turn wetdep on for gases if LAERCHEM=.false.
!            trlist%ti(kt)%nwetdep=1
!            IF (LLDEBUG .AND. MYPROC == 1) THEN
!               WRITE(int_str,*) j_ychem
!               WRITE(int_str_ham,*) kt
!               WRITE(2000+MYPROC,'(a)') 'gas: '//TRIM(YCHEM(j_ychem)%CNAME)
!               WRITE(7000+MYPROC,'(a)') 'OIFS gas: '//TRIM(YCHEM(j_ychem)%CNAME)
!               WRITE(7000+MYPROC,'(a)') 'j_ychem ='//TRIM(int_str)
!               WRITE(7000+MYPROC,'(a)') 'HAM gas: '//TRIM(trlist%ti(speclist(j_spec)%idt)%fullname)
!               WRITE(7000+MYPROC,'(a)') 'gas idt_no ='//TRIM(int_str_ham)
!            END IF
!
!            ! once the tracer is identified we can check the next one
!            CYCLE LABEL_IFS_CHEM_SO4
!         END IF
!      END IF
!   END DO LABEL_GAS_SO4
!   END DO LABEL_IFS_CHEM_SO4
   
ELSE
   ! Note that here we are checking that:
   ! - LAERCHEM is true (should always be the case with M7). It was used to distinguish between chem_scheme='tm5' and "no chemistry" in 43r3/AC with M7 activated (and still the case in 48r1 with AER!).
   ! - CHEM_SCHEME is either 'tm5' or 'SimChem'
   CALL ABOR1(" hamm7_init: UNCOUPLED CHEMISTRY SCHEME "//TRIM(CHEM_SCHEME) )
END IF

! -- LOG
WRITE(NULOUT,'("====== HAMM7_INIT ===== ")')

WRITE(NULOUT,'("Number of  size classes:", I3)') znclass
WRITE(NULOUT,'(" class# / IFS id / HM7 id / IFSNAME / M7NAME ")')
DO J_CLASS = 1,NCLASS
  WRITE(NULOUT,'(1x,I6,3x,I6,3x,I6,3x,A,2x,A)') &
       & J_CLASS, &
       & IND_OIFS_HAM%IND_CLASS_OIFS(J_CLASS), &
       & IND_OIFS_HAM%IND_CLASS_HAM(J_CLASS),  &
       & TRIM(YAERO(IND_OIFS_HAM%IND_CLASS_OIFS(J_CLASS))%CNAME), &
       & TRIM(trlist%ti(sizeclass(J_CLASS)%idt_no)%fullname)
ENDDO

WRITE(NULOUT,'("Number of  mass tracers:", I3)') znaerocomp
WRITE(NULOUT,'("  mass# / IFS id / HM7 id / IFSNAME / M7NAME ")')
DO J_MASS = 1,NAEROCOMP
  WRITE(NULOUT,'(1x,I6,3x,I6,3x,I6,3x,A,2x,A)') &
       & J_MASS, &
       & IND_OIFS_HAM%IND_MASS_OIFS(J_MASS), &
       & IND_OIFS_HAM%IND_MASS_HAM(J_MASS),  &
       & TRIM(YAERO(IND_OIFS_HAM%IND_MASS_OIFS(J_MASS))%CNAME), &
       & TRIM(trlist%ti(aerocomp(J_MASS)%idt)%fullname)
ENDDO

WRITE(NULOUT,'("Number of   gas tracers:", I3)') zsubm_ngasspec
WRITE(NULOUT,'("   gas# / IFS id / HM7 id / IFSNAME / M7NAME  ")')
DO J_GAS = 1,SUBM_NGASSPEC
  J_SPEC = SUBM_GASSPEC(J_GAS)
  WRITE(NULOUT,'(1x,I6,3x,I6,3x,I6,3x,A,2x,A)') &
       & J_GAS, &
       & IND_OIFS_HAM%IND_GAS_OIFS(J_GAS), &
       & IND_OIFS_HAM%IND_GAS_HAM(J_GAS),  &
       & TRIM(YCHEM(IND_OIFS_HAM%IND_GAS_OIFS(J_GAS))%CNAME), &
       & TRIM(TRLIST%TI(SPECLIST(J_SPEC)%IDT)%FULLNAME)
ENDDO

WRITE(NULOUT,'("Number of cloud tracers:", I3)') zcloudind
WRITE(NULOUT,'(" cloud# / IFS id / HM7 id / IFSNAME / M7NAME  ")')
DO J_CLOUD = IDT_CDNC,IDT_ICNC
  WRITE(NULOUT,'(1x,I6,3x,I6,3x,I6,3x,A,2x,A)') &
       & J_CLOUD-IDT_CDNC+1, &
       & IND_OIFS_HAM%IND_CLOUD_OIFS(J_CLOUD-IDT_CDNC+1), &
       & IND_OIFS_HAM%IND_CLOUD_HAM(J_CLOUD-IDT_CDNC+1),  &
       & TRIM(YAERO(IND_OIFS_HAM%IND_CLOUD_OIFS(J_CLOUD-IDT_CDNC+1))%CNAME), &
       & TRIM(trlist%ti(IND_OIFS_HAM%IND_CLOUD_HAM(J_CLOUD-IDT_CDNC+1))%fullname)
ENDDO

IF (LLDEBUG .AND. MYPROC == 1) THEN
   WRITE(5001+MYPROC,*) 'HAM class idts =', ind_oifs_ham%ind_class_HAM(:)
   WRITE(5001+MYPROC,*) 'OIFS class idts =', ind_oifs_ham%ind_class_OIFS(:)
   WRITE(5001+MYPROC,*) 'HAM mass idts =', ind_oifs_ham%ind_mass_HAM(:)
   WRITE(5001+MYPROC,*) 'OIFS mass idts =', ind_oifs_ham%ind_mass_OIFS(:)
   WRITE(5001+MYPROC,*) 'HAM gas idts =', ind_oifs_ham%ind_gas_HAM(:)
   WRITE(5001+MYPROC,*) 'OIFS gas idts =', ind_oifs_ham%ind_gas_OIFS(:)
   WRITE(5001+MYPROC,*) 'HAM cloud idts =', ind_oifs_ham%ind_cloud_HAM(:)
   WRITE(5001+MYPROC,*) 'OIFS cloud idts =', ind_oifs_ham%ind_cloud_OIFS(:)
END IF

END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('HAMM7_INIT',1,ZHOOK_HANDLE)
!END ASSOCIATE
END SUBROUTINE hamm7_init
