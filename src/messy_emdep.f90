MODULE messy_emdep

  ! MESSy-SMCL FOR SUBMODEL EMDEP
  !
  ! EMDEP = CONTROL MODULE FOR SUB-SUB-MODELS
  !           EMIS, DRYDEP, XTSURF, LNOX
  !
  ! AUTHOR: Patrick Joeckel, MPICH, Sep 2003

  IMPLICIT NONE
  PRIVATE

  ! GLOBAL PARAMETERS
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'emdep'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '1.0b' ! mz_sw_20040121 (minor udates of 1.0a)

  ! ESS_lg_20120721+
  CHARACTER(LEN=250), PUBLIC :: casename = '' ! ESS_lg_20130402+
  LOGICAL, PUBLIC :: l_emdep = .false.      ! GLOBAL SWITCH
  LOGICAL, PUBLIC :: l_emis = .false.
  LOGICAL, PUBLIC :: l_drydep = .false.
  LOGICAL, PUBLIC :: l_xtsurf = .false.
  LOGICAL, PUBLIC :: l_wetskinRH = .false. ! ESS_lg_20150618+
  LOGICAL, PUBLIC :: l_fvpd = .false.      ! ESS_lg_20130424+
  
  INTEGER, PUBLIC :: ndtgstart = 2010062900
  INTEGER, PUBLIC :: ntimestep = 1440 ! 1440 timesteps for 1 day with dt=60
  INTEGER, PUBLIC :: nwrite = 30

  REAL, PUBLIC    :: dtimestep = 60.
  REAL, PUBLIC    :: zlatitude = 2.5
  REAL, PUBLIC    :: zlongitude = -65. ! ESS_lg_20130817+
  REAL, PUBLIC    :: zpress = 101325.
  REAL, PUBLIC    :: zslf = 1.
  REAL, PUBLIC    :: zvegfrac = 1.
  REAL, PUBLIC    :: zforestfr = 1.
  REAL, PUBLIC    :: zalbedo = 0.12
  REAL, PUBLIC    :: zlai = 6.
  REAL, PUBLIC    :: zcanheight = 20.
  REAL, PUBLIC    :: zMLHmax = 64. 
  REAL, PUBLIC    :: zrheightsl = 42. ! ESS_lg_20130408+ reference height surface layer, half between zcanheight and zMLHmax
  REAL, PUBLIC    :: zroughness = 2.
  REAL, PUBLIC    :: zqm1 = 0.012
  REAL, PUBLIC    :: zws = 0.2
  REAL, PUBLIC    :: zwsmax = 0.3
  REAL, PUBLIC    :: zprc = 0.
  REAL, PUBLIC    :: zprl = 0.
  REAL, PUBLIC    :: zwetskin =0.

  INTEGER, PUBLIC :: ziladprof = 1 
  ! ESS_lg_20120721-

  PUBLIC :: emdep_read_nml_ctrl

CONTAINS

  ! ----------------------------------------------------------------------
  SUBROUTINE emdep_read_nml_ctrl(status, iou)

    ! EMDEP MODULE ROUTINE (CORE)
    !
    ! READ EMDEP NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES
    !
    ! Author: Patrick Joeckel, MPICH, Sep 2003

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    ! ESS_lg_20120721+
    NAMELIST /CTRL/  casename, l_emis, l_drydep, l_xtsurf,                              & ! ESS_lg_20130402+ case
                     ndtgstart, ntimestep, nwrite, dtimestep,                           &
                     zlatitude, zlongitude, zpress, zslf, zvegfrac, zforestfr, zalbedo, & ! ESS_lg_20130817+
                     zlai, zcanheight, zMLHmax, zrheightsl, zroughness, ziladprof,      & ! ESS_lg_20130408+ zrheightsl
                     zqm1, zws, zwsmax, zprc, zprl, l_wetskinRH, zwetskin, l_fvpd         ! ESS_lg_20150618+, lwetskinRH, ESS_lg_20130424+ lfvpd
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='emdep_read_nml_ctrl'
    LOGICAL      :: lex          ! file exists ?
    INTEGER      :: fstat        ! file status

    ! INITIALIZE
    status = 1 ! ERROR

    ! INITIALIZE GLOBAL CONTROL VARIABLES
    ! -> DEFAULT VALUES ARE SET AT DECLARATION ABOVE

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)

    ! ESS_lg_20120721+
    WRITE(*,*)'CTRL:'
    WRITE(*,*)  'Case study    = ',casename ! ESS_lg_20130402+
    WRITE(*,*)  'l_emis        = ',l_emis
    WRITE(*,*)  'l_drydep      = ',l_drydep 
    WRITE(*,*)  'l_xtsurf      = ',l_xtsurf
    WRITE(*,'(1a,i10)')   ' Start reference date/time = ',ndtgstart
    WRITE(*,'(1a,i5)')    ' # of timesteps        = ',ntimestep
    WRITE(*,'(1a,i4)')    ' freq. of writing data = ',nwrite
    WRITE(*,'(1a,f5.0)')  ' length of timestep [s]= ',dtimestep
    WRITE(*,'(1a,f6.1)')  ' latitude [deg]        = ',zlatitude
	WRITE(*,'(1a,f6.1)')  ' longitude [deg]       = ',zlongitude ! ESS_lg_20130817+
    WRITE(*,'(1a,f5.0)')  ' pressure [hPa]        = ',1e-2*zpress
    WRITE(*,'(1a,f4.1)')  ' surf. land fraction   = ',zslf
    WRITE(*,'(1a,f4.1)')  ' vegetation fraction   = ',zvegfrac
    WRITE(*,'(1a,f4.1)')  ' forest fraction       = ',zforestfr
    WRITE(*,'(1a,f4.2)')  ' albedo                = ',zalbedo
    WRITE(*,'(1a,f4.1)')  ' LAI [m2 m-2]          = ',zlai
    WRITE(*,'(1a,i2)')    ' LAD prof. (1=uniform, 2= ~75% in crownlayer) = ',ziladprof
    WRITE(*,'(1a,f4.1)')  ' canopy height [m]     = ',zcanheight
    WRITE(*,'(1a,f6.0)')  ' max. mixed layer depth [m] = ',zMLHmax
    WRITE(*,'(1a,f4.1)')  ' ref. height SL [m]    = ',zrheightsl
    WRITE(*,'(1a,f4.1)')  ' roughness length [m]  = ',zroughness
    WRITE(*,'(1a,f5.3)')  ' SL moist.[g H2O g-1]  = ',zqm1
    WRITE(*,'(1a,f4.1)')  ' soil moisture [m]     = ',zws
    WRITE(*,'(1a,f4.1)')  ' field capacity [m]    = ',zwsmax
    WRITE(*,'(1a,f4.1)')  ' convective prec.[mm]  = ',1e3*zprc
    WRITE(*,'(1a,f4.1)')  ' large scale prec.[mm] = ',1e3*zprl
    WRITE(*,*)'switch to calc. wet skin fraction from RH = ',l_wetskinRH ! ESS_lg_20150618+
    WRITE(*,'(1a,f4.1)')  ' wet skin fraction     = ',zwetskin
    WRITE(*,*)'switch for vap. press. def. = ',l_fvpd ! ESS_lg_20130424+
    ! ESS_lg_20120721-

    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    ! CHECK NAMELIST CONSISTENCY
    l_emdep = .true.

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE emdep_read_nml_ctrl

END MODULE messy_emdep
