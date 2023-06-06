! Author:
! Laurens Ganzeveld,     MPICH, 2005, ESS 2012/2013

MODULE messy_emdep_xtsurf_box

  USE messy_emdep
  USE messy_emdep_mem ! mz_lg_20050721+
  USE messy_emdep_xtsurf
  USE messy_emdep_emis_mem
  USE messy_emdep_emis
  USE messy_main_constants_mem, ONLY: dp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: emdep_xtsurf_initialize
  PUBLIC :: emdep_xtsurf_tracer_init ! ESS_lg_20120722+
  PUBLIC :: emdep_xtsurf_vdiff

CONTAINS

  ! --------------------------------------------------------------------------

  SUBROUTINE emdep_xtsurf_initialize

    IMPLICIT NONE
    INTEGER                     :: status ! error status

    ! read CTRL namelist
    CALL emdep_read_nml_ctrl(status, 99)
    ! mz_lg_20050719+ reading the settings of the xtsurf module
    ! read CTRL_XTSURF namelist
    CALL emdep_xtsurf_read_nml_ctrl(status, 99)
    ! mz_lg_20050719+ emissions are also needed
    ! read CTRL_EMIS namelist
    CALL emdep_emis_read_nml_ctrl(status, 99)
    IF (status /= 0) STOP

  END SUBROUTINE emdep_xtsurf_initialize

  ! --------------------------------------------------------------------------

  SUBROUTINE emdep_xtsurf_vdiff

    USE messy_main_constants_mem,   ONLY: g, dp, i4, pi, &
                                        avo=>N_A, amd=>M_air ! mz_lg_20050719+

    IMPLICIT NONE

    ! general declarations 
    INTEGER             :: i, ii, jt, jk, ip

    ! some constants
    REAL(dp), PARAMETER ::    vkarman=0.41  !von Karman constant [-]

    ! timestepping parameters
    INTEGER,  PARAMETER :: nobs_max  = 100              ! MAQ_lg_20171015+ increased to 100, maximum number of observed parameters

    INTEGER  :: ndtginit,                            &  ! starting date/time [yyyymmddhh]
                icent,                               &  ! century
                Nstep,                               &  ! actual # of timesteps being used
                nstepday,                            &  ! # of timesteps per day
                nstepinit,                           &  ! # of timesteps to correct for NDTGINIT different from 0
                nday,                                &  ! number of simulated days
                nprint                                  ! number of timesteps that data are written to output files
    
    CHARACTER(len=255) dummy, dummy2, dummy3            ! ESS_lg_20140416+ dummy2, dummy3
                                
    REAL(dp) :: delta_time                              ! lenght of timestep in [s]

    REAL(dp) :: DayReal, SoDecli, Latitu, sinpsi, photon, fct, dn, doffset
    REAL(dp), PARAMETER :: dusk = 0.0721347 ! = 5.E-2/LOG(2.) = photon at dusk
    REAL(dp), PARAMETER :: Cancer = 23.45 * (pi/180.)
		
    REAL(dp) :: fAPIN, fBPIN, fSQTERP, fATERP, fLIMO, fMYRC ! MAQ_lg_20170413+ added extra MT's ! ESS_lg_20130119+ added the partitioning of monoterpene and sqterp emissions
	REAL(dp) :: vpair, vpleaf, vpsfc, rbveg             ! ESS_lg_20130424+

    INTEGER  :: ntrac
    REAL(dp) :: recalcmrconc, VdO3_canopy               ! MAQ_lg_20220127+ added diagnostic term to secure also writing < 0 values! term to recalculate mixing ratio to conc.
    REAL(dp) :: recalctend, recalctendppmv              ! terms to recalculate tendencies to ppbv/ppmv hr-1
	INTEGER, PARAMETER :: nunits = 4
        
    ! mz_lg_20040716+ surface cover parameters
    INTEGER             :: iladprof

    INTEGER, DIMENSION(:), ALLOCATABLE ::  &
             ndtgact,                             &  ! actual date/time
             iyear,                               &  ! year
             imon,                                &  ! month
             iday,                                &  ! day
             ihr,                                 &  ! hour
             imin,                                &  ! minutes
             isec,                                &  ! seconds
             Jday                                    ! julian day of simulation
 
    REAL(dp), DIMENSION(:), ALLOCATABLE :: &
        timeday                                  ! time at the day        

    CHARACTER(LEN=20), DIMENSION(:), ALLOCATABLE ::   &
      ldatltime
          
    ! ESS_lg_20140416+ added declarations to write headers and parameters for the changing vertical resolution model 
    CHARACTER(LEN=14), DIMENSION(:,:), ALLOCATABLE ::   & ! maximum length 14 characters that allows 1 space for 15 char. column
	  parname, parunit
    CHARACTER(LEN=14) :: mixratio_unit(nunits)
	DATA mixratio_unit /'1E-3 pptv','pptv','ppbv','ppmv'/
	  
	REAL, DIMENSION(:,:,:),   ALLOCATABLE :: mixratios      ! recalculating/assigning mixing ratios for output
    REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: tendencies     ! recalculating/assigning tendencies for output
    REAL :: recalcconc
	REAL :: recalc_conc_VCD ! ESS_lg_20150714+ recalculation from concentration to Vertical Column Density
    INTEGER,  PARAMETER :: npar_out  = 12, ntend = 5        ! max. number of output parameters that require a flexible parameter name assignment (see veg_mlay.out) 
    INTEGER :: iwind, iKh, ijNO2, iVdO3, iOHr, istomfluxO3, iconc,  &
	           ixteemis, ixtedryd, ixtechem, ixtediff, ixtetend, &
			   npar, nlevels(npar_out),                &   ! counter of number of output parameters and number of levels for each parameter 
               jjt, jjt_emis, jjt_dryd, jjt_chem, jjt_diff, jjt_xte, iunit
    ! ESS_lg_20140416-	
	  
    INTEGER :: iip	
    ! ESS_lg_20140416-
	  
    REAL(dp), DIMENSION(:), ALLOCATABLE :: &
      slf         , & ! land-sea mask [0-1]
      slm         , & ! land-sea mask [0/1]
      vgrat       , & ! vegetation fraction [0-1]
      forestfr    , & ! forest fraction [0-1]
      dm          , & ! foliar density [g m-2]
	  dmbase      , & ! long-term average foliar density [g m-2] ! MAQ_lg_20160817+
      lai         , & ! leaf area index [m m-2]
	  laibase     , & ! long-term average leaf area index [m m-2] ! MAQ_lg_20160621+
      hc          , & ! canopy height [m]         ! mz_lg_20050719+ 
      disp        , & ! displacement height [m]   ! mz_lg_20050719-
      zrefsl          ! ESS_lg_20130408+ 
	  
    LOGICAL, DIMENSION(:), ALLOCATABLE  :: &
      loland          ! land-sea mask, logical

    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: &
      lad         , & ! Leaf Area Density (LAD) profile [0-1]
	  zrefcan     , & ! ESS_lg_20140416+ height of the canopy layers
	  zrefcan_hr  , & ! ESS_lg_20140416+ height of the canopy layers for the high-resolution calculations
      tcan            ! MAQ_20170505+ added Tcanopy to use a temperature gradient inside the canopy
	  
    ! some general input parameters
    REAL(dp), DIMENSION(:), ALLOCATABLE :: &
      scaling_rad , & ! scaling function for radiation
      globrad     , & ! global radiation [W m-2]
      netrad      , & ! net surface radiation [W m-2]
      zen         , & ! zenith angle [0-90]
      cossza      , & ! cosine of the zenith angle [0-1.56]
      tair        , & ! air temperature [K]
      tsurf       , & ! surface/skin/canopy temperature [K]
      tsoil       , & ! soil temperature [K] 
      press       , & ! surface layer pressure [Pa]
      qm1         , & ! surface layer moisture [g H2O g-1 air]
      qs          , & ! saturated surface layer moisture [g H2O g-1 air] ! MAQ_lg_20170720+
      qsam        , & ! saturated moisture at leaf level [g H2O g-1 air] ! MAQ_lg_20170720+
      ws          , & ! soil moisture [m]
      wsmax           ! field capacity [m] 

    ! ESS_lg_20120717+ some extra input to introduce temporal cycle in input parameters
    REAL(dp):: &
      scaling     , & ! scaling function #1
      globradmax  , & ! maximum global radiation
	  netradmax   , & ! maximum net radiation (cloud free)
      albedo      , & ! surface albedo [0-1] 
      tmin        , & ! minimum surface temperature of temporal cycle
      tmax        , & ! maximum surface temperature of temporal cycle
      tamin       , & ! minimum air temperature of temporal cycle
      tamax       , & ! maximum air temperature of temporal cycle
      tslmin      , & ! minimum soil temperature
      tslmax      , & ! maximum soil temperature
      umin        , & ! minimum u-component windspeed
      umax        , & ! maximum u-component windspeed
      vmin        , & ! minimum v-component windspeed
      vmax        , & ! maximum v-component windspeed
	  rhavg       , & ! MAQ_20140414+ average RH
      rhdif           ! maximum difference in RH relative to average 

    ! ESS_lg_20120722+ read-in observations        
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: &
      observ          ! observed parameters

    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: &
      readdata

    ! canopy radiation and turbulence parameters
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: &
      rvd         , & ! diffusive rad. component [W m-2]
      fsl         , & ! fraction of sunlit leaves [0-1]
	  fslbase     , & ! fraction of sunlit leaves [0-1] for long-term average LAI ! MAQ_lg_20160621+
      PAR         , & ! PAR inside the canopy [umol m-2] ! MAQ_lg_20180112+	  
      u_veg       , & ! mz_lg_20050719+ wind speed in canopy [m s-1]
      Kh              ! ESS_lg_20120722+ eddy-diffusivity heat/tracers

    REAL(dp), DIMENSION(:), ALLOCATABLE :: &
      rbvd        , & ! direct beam component [W m-2]
      rvdsl       , & ! surface layer difussive radiation [W m-2] ! mz_lg_20050721+
      Kh_obs_sl   , & ! Surface layer exchange coefficient for heat (and scalars) above canopy [m2 s-1] !MAQ_AV_20200309+
      Kh_obs_cl   , & ! Surface layer exchange coefficient for heat (and scalars) in canopy [m2 s-1]    !MAQ_AV_20200309+,
      zustveg_obs     ! Observed friction velocity over forest [m s-1]                 ! MAQ_AV_20201005+                                      ! a range from 1-5	

    ! mz_LG_20040716+ parameters for soil-biogenic NO emission simulations
    ! some extra timestepping parameters
    INTEGER :: &
      init_step   , & ! initial timestep
      ndaylen     , & ! length of day in seconds
      imonth          ! month

    ! ESS_lg_20130107+ added the assignment of the names of the NO emission classes
    CHARACTER*25, DIMENSION(:), ALLOCATABLE :: &
       NOemclassname  ! names of the NO emission classes
    ! ESS_lg_20130107-
           
    REAL(dp), DIMENSION(:), ALLOCATABLE :: &
      cultiv      , & ! cultivation intensity [old: 0-15, new 2004: 0-1]
      fertil      , & ! fertilizer application
      prc         , & ! convective precipitation [m]
      prl         , & ! large-scale precipitation [m]
      prectot     , & ! monthly accumulated precipitation [m]
      noemis_w    , & ! wet soil emission factor [ng N m-2 s-1]
      noemis_d        ! dry soil emission factor [ng N m-2 s-1]

    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: &
      cpold       , & ! convective precipitation record
      lspold          ! large-scale precipitation record

    REAL(dp), DIMENSION(:), ALLOCATABLE :: &
      cp          , & ! convective precipitation [m]
      lsp         , & ! large-scale precipitation [m]  
      pulsing     , & ! pulsing regime    [index, 1-3]
      plsday      , & ! timing of pulse   [number of timesteps that pulse is active]
      plsdurat    , & ! duration of pulse [days]
      pls         , & ! the actual pulse  [ - ]
      crf         , & ! mz_lg_20050522+ added the Canopy Reduction Fact. [0-1]
      ! ESS_lg_20130117+ added for the HONO and NOx emissions from NO3 photolysis
      fslsum      , & ! canopy intergrated fraction of sunlit leaves [0-1]
      cthru       , & ! throughfall associated with rainfall; for a value of 1 all the canopy is wetted and NO3 being removed
      HONO_jNO3em , & ! "big leaf" HONO emission flux [molecules m-2 s-1]
      NOx_jNO3em  , & ! "big leaf" NOx emission flux [molecules m-2 s-1]
      ! ESS_lg_20130117-
      NO_slflux   , & ! NO soil flux [molecules m-2 s-1] ! mz_lg_20050719+ 
      NO2_slflux  , & ! MAQ_lg_20210727+ added NO2 soil (deposition) flux for extra diagnostics
      NO_emflux   , & ! NO soil-biogenic emission flux [molecules m-2 s-1]
      HONO_slflux , & ! soil HONO emission flux [molecules m-2 s-1]
      NH3_slflux  , & ! MAQ_lg_20230329+ soil NH3 emission flux [molecules m-2 s-1]
	  SQT_slflux      ! MAQ_lg_20230601+ soil SQT emission flux [molecules m-2 s-1]

    ! radon emissions (to test turbulent transport) 
    REAL(dp), DIMENSION(:), ALLOCATABLE :: &
      radon_slflux, &    ! radon emission flux [atoms m-2 s-1]
      co2_slflux         ! ESS_lg_20130503+ and CO2 to test canopy exchange

	! mz_lg_20040921+ added the latitude and longitude
    REAL(dp), DIMENSION(:), ALLOCATABLE :: &
      latitude, &        ! latitude [degrees], SH (-90 - 0), NH (0 -90)
      longitude          ! longitude
	  
    INTEGER, PARAMETER :: ncl_noemis  =    12 ! No. of NO emission classes [YL95: 12  DK97: 17]
    INTEGER  iNOemcl  ! ESS_lg_20130107+ modified to assigned read-in emission class

    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: &
      noemclass1  , & ! fractional coverage for each NO emission class [YL95: 0-12  DK97: 0-17]
      ddNO3       , & ! NO3 dry deposition flux [molecules m-2 s-1]
      NO3s            ! accumulated amount of NO3 on leaf surface

	LOGICAL :: lstart ! switch for starting simulation

    ! mz_LG_20040716+ parameters needed for biogenic VOC emissions
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: &
      VOC_emfact  , & ! VOC emission factors [ug C g-1 hr-1]
      VOC_emflux      ! VOC emission fluxes

    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: &
      ISOP_emflux , & ! isoprene emission fluxes [molecules m-2 s-1]
      MONO_emflux , & ! monoterpene emission fluxes [molecules m-2 s-1]
      OVOC_emflux , & ! other VOC emission fluxes [molecules m-2 s-1]
      HONO_emflux , & ! HONO emissions from HNO3 photolysis [molecules m-2 s-1]
      NOx_emflux  , & ! NOx emissions from HNO3 photolysis [molecules m-2 s-1]
	  NH3_emflux      ! MAQ_lg_20230329+ NH3 emissions by vegetation [molecules m-2 s-1]

    ! parameters used in the calculation aerodynamic resistance
    REAL(dp), DIMENSION(:), ALLOCATABLE    :: &
      cfml        , & ! drag coefficient over land
      cfncl       , & ! scaling term land
      ril         , & ! Richardson number over land [-] 
      cdnl        , & ! neutral drag coefficient over land
      geopot_3d   , & ! geopotential height, BL
      geopot_3d_sl, & ! ESS_20130208+ geopotential height, SL
      tvir        , & ! virtual temperature [K]
      tvl         , & ! virtual temperature over land [K]   
      um1         , & ! windspeed, u component [m s-1]
      vm1         , & ! windspeed, v component [m s-1]
	  u_sl        , & ! ESS_lg_20131125+ surface layer wind speed [m s-1]
      az0         , & ! echam's surface roughness for momentum [m]  
      z0m         , & ! vegetation surface roughness [m]         
      z0mslsn     , & ! snow/ice/bare soil surface roughness [m]   
      zref        , & ! surface layer reference height 
      ! mz_lg_20050719+ added some extra parameters defining the grid      
      grvol       , & ! surface layer grid volume [m3]
      grmass      , & ! surface layer grid volume [kg]
      pdp         , & ! surface layer pressure thickness
      prhoa       , & ! surface layer density [kg m-3]
      MLH             ! ESS_lg_2013014+ added the mixed layer height
      
    REAL(dp), DIMENSION(:), ALLOCATABLE :: &
      ! output (of calculation aerodynamic resistance)
      zrahl       , & ! aerodynamic resistance over land [s m-1]
      zrahveg     , & ! aerodynamic resistance over vegetation [s m-1]
      zrahslsn    , & ! aerodynamic resistance over bare soil/snow-ice [s m-1]
      zustarl     , & ! friction velocity over land [m s-1]
      zustveg     , & ! friction velocity over vegetation [m s-1]
      zustslsn    , & ! friction velocity over bare soil/snow-ice [m s-1]

      ! input (for calculation dry deposition velocity)
      zfrl        , & ! land fraction ! mz_lg_20040619+ added
      zcvbs       , & ! cover with bare soil, bare soil fraction
      zcvwat      , & ! cover with water, water fraction
      zcvs        , & ! cover with snow, snow fraction
      zcvw        , & ! cover with wet vegetation, wet skin fraction
      zvgrat      , & ! vegetation fraction
      rh          , & ! relative humidity [0-1]
	  rco_leaf_AGS, & ! ESS_lg_20130516+ AGS big leaf stomatal resistance
      rmesCO2     , & ! ESS_lg_20130516+ time dependent mesophyllic resistance CO2
	  rcutCO2     , & ! ESS_lg_20130516+ time dependent cuticular resistance CO2
	  ccompCO2    , & ! MAQ_lg_20170720+ added the CO2 compensation point
      fws         , & ! soil moisture stress attenuation function [0-1]
	  fvpd        , & ! ESS_lg_20130424+ added vapor pressure deficit term
	  DS_out      , & ! MAQ_AV_20201023+ added diagnostic VPD output
	  D0_out

    ! ESS_lg_20130516+ for the AGS model
	REAL(dp), DIMENSION(:,:), ALLOCATABLE :: &                          				  
      rco_leaf    , & ! ESS_lg_20130516+ leaf stomatal resistance, changed to make it layer dependent [s m-1]
	  rco_leaf_AGS_ml ! AGS multi layer stomatal resistance considering CO2 concentrations
						  
	REAL(dp), DIMENSION(:,:), ALLOCATABLE :: &
      soilph          ! soil pH, 7 classes

    ! mz_lg_20040720+ added some extra parameters needed to calculate
    !     the leaf stomatal resistance for box model studies that include
    !     the diurnal cycle in the vegetation uptake as a function of
    !     surface net radiation

    REAL(dp)  zwcr, zwpwp  ! soil moisture threshold parameters (e.g., permanent wilting point)

    REAL(dp), PARAMETER :: cva = 5000.  !   constant to define the stomatal resistance
    REAL(dp), PARAMETER :: cvb =   10.  !   constant to define the stomatal resistance
    REAL(dp), PARAMETER :: cvc =  100.  !   minimum stomatal resistance
    REAL(dp), PARAMETER :: cvk =     .9 !
    REAL(dp), PARAMETER :: cvbc= cvb*cvc!   cvb*cvc
    REAL(dp), PARAMETER :: cvkc= cvk*cvc!   cvk*cvc
    REAL(dp), PARAMETER :: cvabc=(cva+cvbc)/cvc   !   (cva+cvbc)/cvc
    REAL(dp), PARAMETER :: cvrad = 0.55 !   fraction of the net s.w radiation contributing to p.a.r
 
    REAL(dp) :: zva, zvabc, zvb, zvbc, zvc
    REAL(dp) :: zvk, zvkc, zvklt, zvrad
    REAL(dp) :: znetrad, zepsr, zabcs

    ! mz_lg_20040721+ and parameters needed to calculate roughness weighted  with all the surface
    !    cover fractions through the drag coeff. for surface cover fraction
    REAL(dp) :: zvisc, zustarn, zustarm, zcdragw, zcdragv, zcdragslsn, zcdrag
        
    ! ESS_lg_20120928+ added to consider a variable mixed layer (boundary layer) depth 
    REAL(dp) :: zrefmin, MLHmax
    
	! MAQ_lg_20230403+ added to calculate qsat and qm
	REAL(dp) :: esat, qsat
	
    CHARACTER(LEN=20), DIMENSION(:), ALLOCATABLE ::   &
      trname 

    ! declaration of tracer properties that are used to calculate their surface uptake resistances
    REAL(dp), DIMENSION(:), ALLOCATABLE :: &
      moleweight  , & ! molecular weight, needed for estimating surface resistances 
      reactivity  , & ! reactivity coefficient [0:non-react., 0.1:semi react., 1:react.]  
      henrycoeff      ! henrycoeff coefficient [mol atm-1]

    LOGICAL      :: & 
      lo_derived  , & ! true whenever SO2 and O3 are defined as tracers
      lexist_O3       ! true whenever O3 is present

    LOGICAL, DIMENSION(:), ALLOCATABLE   :: & 
      lvd_bigl    , & ! true for the tracer when the required parameters are defined
      lexist_GAS      ! true for the tracer when it is a gas

    ! different dry deposition resistance terms
    REAL(dp), DIMENSION(:), ALLOCATABLE  :: &
      diff        , & ! term that corrects for difference in diffusivity, stomatal exchange
      diffrb      , & ! term that corrects for difference in diffusivity, quasi-laminary boundary layer res.
      rsoil       , & ! soil resistance [s m-1]
      rwater      , & ! water resistance
      rws         , & ! wet skin resistance
      rsnow       , & ! snow/ice resistance
      rmes        , & ! mesophyll resistance
      rcut        , & ! cuticular resistance
      ccomp       , & ! mz_lg_20050721+ added the compensation point
      rj_max          ! and maximum surface layer photolysis rates

    ! mz_lg_20040721+ extra parameters of the aerosol dry deposition code.
    INTEGER, PARAMETER :: &
      AIR=             1, & ! index indicating tracer being a gas
      AEROSOL=         2    ! index indicating tracer being an aerosol
      
    INTEGER, DIMENSION(:), ALLOCATABLE :: &
      medium  

    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: &
      ! input (general)
      pxtm1       , & ! surface layer mixing ratio at t=t+1 ! mz_lg_20050719+
      pxtm1_obs       ! observed surface layer mixing ratios

    ! mz_lg_20050719- end declaration dry deposition and turbulence 
    !    parameters and start some additional declarations for the 
    !    atmosphere-biosphere model calculations 

    LOGICAL, DIMENSION(:), ALLOCATABLE   :: & 
      latmbios,       & ! true for atmosphere-biosphere calculations
      latmbios_emveg, & ! true for tracer vegetation emissions
      latmbios_emsoil,& ! true for tracer soil emissions
      latmbios_photo, & ! true for photolysis rates
	  latmbios_output   ! ESS_lg_20140418+ to define output parameters

    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: &
      ! input (general)
      patmbiosflux,   & ! atmosphere-biosphere exchange flux [molec. m-2 s-1]
      pcrf,           & ! Canopy Reduction Factor [ratio]
      psurfflux,      & ! surface exchange flux [molec. m-2 s-1]
	  pstomfluxO3,    & ! ESS_lg_20140509+ O3 stomatal flux
      rj,             & ! surface layer photolysis rates
      OHreact           ! OH reactivity
          
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE:: &
      ! input (general)
      pxtmveg,        & ! canopy mixing ratio at t=t+1 ! mz_lg_20050719+
      rj_veg,         & ! within-canopy photolysis rates
      pvdveg,         & ! within-canopy dry deposition velocities
      xteemis,        & ! emission tendency [molecules g-1 s-1]
      xtedryd,        & ! dry deposition tendency [molecules g-1 s-1]
      xtechem,        & ! chemistry tendency [molecules g-1 s-1]
      xtediff           ! vertical diffusion tendency [molecules g-1 s-1]

    ! ===============================================================================
    ! mz_lg_20050719+ start of initializing the data for actual calculations
    
    ! ESS_lg_20120717+ this can be replaced by 
    !    reading in observed parameters to analyse the model performance under realistic conditions

    ndtginit                   = ndtgstart          ! starting date and time [yyyymmddhh]
    Nstep                      = ntimestep          ! actual # of timesteps being used
    delta_time                 = dtimestep          ! lenght of timestep in [s]
    nstepday                   = 86400/delta_time   ! # of timesteps per day
    nprint                     = nwrite             ! frequency of writing output

    ! ESS_lg_20130104+ assigning the number of tracers dependent on if the chemistry scheme
           !   is used or not, no chemistry: ntrac=16, with chemistry: ntrac=62 
    ntrac=18                             ! MAQ_lg_20190110+ added COS ! ESS_lg_20130503+ added CO2
    IF (l_xtsurf_veg_mlay_chem) ntrac=67 ! MAQ_lg_20190110+ added COS ! MAQ_lg_20170412+ added 3 monoterpene species. ! ESS_lg_20130503+ added CO2

    ! MAQ_lg20170512+ added netradmax initialization
	netradmax=0.

    ! ESS_lg_20130118+ start allocating all the arrays
    ALLOCATE(ndtgact(Nstep))
    ALLOCATE(iyear(Nstep)) 
    ALLOCATE(imon(Nstep)) 
    ALLOCATE(iday(Nstep)) 
    ALLOCATE(ihr(Nstep)) 
    ALLOCATE(imin(Nstep)) 
    ALLOCATE(isec(Nstep)) 
    ALLOCATE(Jday(Nstep))
    ALLOCATE(timeday(Nstep))
    ALLOCATE(ldatltime(Nstep))
        
    ALLOCATE (slf(Nstep))     
    ALLOCATE (slm(Nstep))        
    ALLOCATE (vgrat(Nstep))       
    ALLOCATE (forestfr(Nstep))   
    ALLOCATE (dm(Nstep))  
    ALLOCATE (dmbase(Nstep))  ! MAQ_lg_20160817+   	
    ALLOCATE (lai(Nstep)) 
    ALLOCATE (laibase(Nstep)) ! MAQ_lg_20160621+	
    ALLOCATE (hc(Nstep)) 
    ALLOCATE (disp(Nstep))             
    ALLOCATE (zrefsl(Nstep)) ! ESS_lg_20130408+
    ALLOCATE (loland(Nstep)) 
	ALLOCATE (scaling_rad(Nstep))
    ALLOCATE (globrad(Nstep))
    ALLOCATE (netrad(Nstep))
    ALLOCATE (zen(Nstep)) 
    ALLOCATE (cossza(Nstep))
    ALLOCATE (tair(Nstep))
    ALLOCATE (tsurf(Nstep))
    ALLOCATE (tsoil(Nstep)) 
    ALLOCATE (press(Nstep)) 
    ALLOCATE (qm1(Nstep))
	ALLOCATE (qs(Nstep))   ! MAQ_lg_20170720+
	ALLOCATE (qsam(Nstep)) ! MAQ_lg_20170720+
    ALLOCATE (ws(Nstep))
    ALLOCATE (wsMAX(Nstep))
    ALLOCATE (rbvd(Nstep)) 
    ALLOCATE (rvdsl(Nstep))
    ALLOCATE (NOemclassname(ncl_NOemis))
    ALLOCATE (cultiv(Nstep)) 
    ALLOCATE (fertil(Nstep)) 
    ALLOCATE (prc(Nstep))
    ALLOCATE (prl(Nstep))
    ALLOCATE (prectot(Nstep)) 
    ALLOCATE (noemis_w(Nstep)) 
    ALLOCATE (noemis_d(Nstep))
    ALLOCATE (cp(Nstep))  
    ALLOCATE (lsp(Nstep))  
    ALLOCATE (pulsing(Nstep))
    ALLOCATE (plsday(Nstep)) 
    ALLOCATE (plsdurat(Nstep))
    ALLOCATE (pls(Nstep)) 
    ALLOCATE (crf(Nstep)) 
    ALLOCATE (NO_slflux(Nstep))
	ALLOCATE (NO2_slflux(Nstep)) ! MAQ_lg_20210727+ added NO2 soil (deposition) flux for extra diagnostics
    ALLOCATE (NO_emflux(Nstep))
    ALLOCATE (RADON_slflux(Nstep))  
    ALLOCATE (CO2_slflux(Nstep)) ! ESS_lg_20130503+	
    ALLOCATE (latitude(Nstep))
	ALLOCATE (longitude(Nstep))  ! ESS_lg_20130817+
    ALLOCATE (cfml(Nstep))  
    ALLOCATE (cfncl(Nstep)) 
    ALLOCATE (ril(Nstep)) 
    ALLOCATE (cdnl(Nstep)) 
    ALLOCATE (geopot_3d(Nstep))
    ALLOCATE (geopot_3d_sl(Nstep)) ! ESS_20130208+ 	
    ALLOCATE (tvir(Nstep)) 
    ALLOCATE (tvl(Nstep)) 
    ALLOCATE (um1(Nstep)) 
    ALLOCATE (vm1(Nstep))
    ALLOCATE (u_sl(Nstep)) ! ESS_lg_20131125+
    ALLOCATE (az0(Nstep)) 
    ALLOCATE (z0m(Nstep))               
    ALLOCATE (z0mslsn(Nstep))     
    ALLOCATE (zref(Nstep))        
    ALLOCATE (grvol(Nstep))   
    ALLOCATE (grmass(Nstep)) 
    ALLOCATE (pdp(Nstep)) 
    ALLOCATE (prhoa(Nstep))
    ALLOCATE (MLH(Nstep)) ! ESS_20130104+ added the mixed layer height        
    ALLOCATE (zrahl(Nstep)) 
    ALLOCATE (zrahveg(Nstep)) 
    ALLOCATE (zrahslsn(Nstep))
    ALLOCATE (zustarl(Nstep))
    ALLOCATE (zustveg(Nstep)) 
    ALLOCATE (zustslsn(Nstep))
    ALLOCATE (zfrl(Nstep)) 
    ALLOCATE (zcvbs(Nstep)) 
    ALLOCATE (zcvwat(Nstep))
    ALLOCATE (zcvs(Nstep))
    ALLOCATE (zcvw(Nstep)) 
    ALLOCATE (zvgrat(Nstep))
    ALLOCATE (rh(Nstep)) 
    ALLOCATE (rco_leaf_AGS(Nstep)) ! ESS_lg_20130516+	
    ALLOCATE (rmesCO2(Nstep))
    ALLOCATE (rcutCO2(Nstep))  ! ESS_lg_20130516-	
    ALLOCATE (ccompCO2(Nstep)) ! MAQ_lg_20170720+ added the CO2 compensation point

    ALLOCATE (DS_out(Nstep))   ! MAQ_AV_20201023+ added diagnostic VPD output
    ALLOCATE (D0_out(Nstep))   ! MAQ_AV_20201023+ added diagnostic VPD output

    ALLOCATE (fws(Nstep)) 
	ALLOCATE (fvpd(Nstep)) ! ESS_lg_20130424+ added vapor pressure deficit term
    ! ESS_lg_20130117+ added for HONO/NOx source from nitrate photolysis
    ALLOCATE (fslsum(Nstep)) 
    ALLOCATE (ddNO3(Nstep,nveglay))   
    ALLOCATE (NO3s(Nstep,nveglay))
    ALLOCATE (cthru(Nstep))
    ALLOCATE (HONO_jNO3em(Nstep))
    ALLOCATE (NOx_jNO3em(Nstep))
    ALLOCATE (HONO_slflux(Nstep)) ! MAQ_lg_20170517+ soil HONO flux
    ALLOCATE (NH3_slflux(Nstep))  ! MAQ_lg_20230329+ soil NH3 flux
    ALLOCATE (SQT_slflux(Nstep))  ! MAQ_lg_20230601+ soil SQT flux
    ! ESS_lg_20130117- 
 
    ! allocation of 2D arrays
    ALLOCATE (lad(Nstep,nveglay_hr)) 
    ALLOCATE (zrefcan(Nstep,nveglay))       ! ESS_lg_20140416+   	
    ALLOCATE (zrefcan_hr(Nstep,nveglay_hr)) ! ESS_lg_20140416+ 
    ALLOCATE (tcan(Nstep,nveglay_hr))       ! MAQ_lg_20170505+ added Tcanopy
    ALLOCATE (observ(Nstep,nobs_max))
    ALLOCATE (readdata(Nstep,nobs_max))
    ALLOCATE (rvd(Nstep,nveglay_hr))  
    ALLOCATE (fsl(Nstep,nveglay_hr))  
    ALLOCATE (fslbase(Nstep,nveglay_hr))    ! MAQ_lg_20160621+ 
	ALLOCATE (PAR(Nstep,nveglay_hr))        ! MAQ_lg_20180112+ added for extra output
    ALLOCATE (u_veg(Nstep,nveglay_hr)) 
    ALLOCATE (Kh(Nstep,nveglay)) 

    ALLOCATE (Kh_obs_sl(Nstep))        ! MAQ_AV_20200309+  
    ALLOCATE (Kh_obs_cl(Nstep))        ! MAQ_AV_20200309+
    ALLOCATE (zustveg_obs(Nstep))      ! MAQ_AV_20201005+

	ALLOCATE (cpold(Nstep,ndrydays))  
    ALLOCATE (lspold(Nstep,ndrydays))
    ALLOCATE (noemclass1(Nstep,ncl_noemis))  
    ALLOCATE (VOC_emfact(Nstep,nspec_vocemis))  
    ALLOCATE (VOC_emflux(Nstep,nspec_vocemis))     
    ALLOCATE (ISOP_emflux(Nstep,nveglay))
    ALLOCATE (MONO_emflux(Nstep,nveglay)) 
    ALLOCATE (OVOC_emflux(Nstep,nveglay)) 
    ALLOCATE (HONO_emflux(Nstep,nveglay))
    ALLOCATE (NOx_emflux(Nstep,nveglay))
	ALLOCATE (NH3_emflux(Nstep,nveglay))      ! MAQ_lg_20230329+ added NH3
    ALLOCATE (soilph(Nstep,ncl_soilph)) 
    ALLOCATE (rco_leaf(Nstep,nveglay))        ! ESS_lg_20130516+ moved to here; made it layer dependent
	ALLOCATE (rco_leaf_AGS_ml(Nstep,nveglay)) ! ESS_lg_20130516+ 
    
    ! allocation of tracer paramters, 1D
    ALLOCATE (trname(ntrac))
    ALLOCATE (moleweight(ntrac))  
    ALLOCATE (reactivity(ntrac)) 
    ALLOCATE (henrycoeff(ntrac)) 
    ALLOCATE (lvd_bigl(ntrac))
    ALLOCATE (lexist_GAS(ntrac))  
    ALLOCATE (diff(ntrac))  
    ALLOCATE (diffrb(ntrac))
    ALLOCATE (rsoil(ntrac))
    ALLOCATE (rwater(ntrac))
    ALLOCATE (rws(ntrac)) 
    ALLOCATE (rsnow(ntrac))
    ALLOCATE (rmes(ntrac))
    ALLOCATE (rcut(ntrac))
    ALLOCATE (ccomp(ntrac))
    ALLOCATE (rj_MAX(ntrac)) 
    ALLOCATE (medium(ntrac))  
    ALLOCATE (latmbios(ntrac))
    ALLOCATE (latmbios_emveg(ntrac))
    ALLOCATE (latmbios_emsoil(ntrac))
    ALLOCATE (latmbios_photo(ntrac))
	ALLOCATE (latmbios_output(ntrac)) ! ESS_lg_20140418+

    ! 2D arrays
    ALLOCATE (pxtm1(Nstep,ntrac))    
    ALLOCATE (pxtm1_obs(Nstep,ntrac))  
    ALLOCATE (patmbiosflux(Nstep,ntrac))   
    ALLOCATE (pcrf(Nstep,ntrac))
    ALLOCATE (psurfflux(Nstep,ntrac))
	ALLOCATE (pstomfluxO3(Nstep,nveglay))
    ALLOCATE (rj(Nstep,ntrac))
    ALLOCATE (OHreact(Nstep,nveglay+1)) ! ESS_lg_20130113+

    ! 3D arrays
    ALLOCATE (pxtmveg(Nstep,nveglay,ntrac))  
    ALLOCATE (rj_veg(Nstep,nveglay,ntrac))  
    ALLOCATE (pvdveg(Nstep,nveglay,ntrac))
    ! ESS_lg_2013016+ process tendencies
    ALLOCATE (xteemis(Nstep,nveglay+1,ntrac))
    ALLOCATE (xtedryd(Nstep,nveglay+1,ntrac))
    ALLOCATE (xtechem(Nstep,nveglay+1,ntrac))
    ALLOCATE (xtediff(Nstep,nveglay+1,ntrac))
    ! ESS_lg_2013016- 

    ! ====================================================================================================
    ! start assigning parameter values
        
    ! ESS_lg_20130118+ calculation of century and initial year, month, day and hour from ndtginit
    icent      = (ndtginit/1e8)
    iyear(1)   = (ndtginit/1e6)-icent*1e2
    imon(1)    = (ndtginit/1e4)-icent*1e4-iyear(1)*1e2
    iday(1)    = (ndtginit/1e2)-icent*1e6-iyear(1)*1e4-imon(1)*1e2
    ihr(1)     = (ndtginit/1)  -icent*1e8-iyear(1)*1e6-imon(1)*1e4-iday(1)*1e2 ! ESS_lg_20130415+ modified calculation of ihr

    ! ESS_lg_20150310+ to be used in the calculation of DayReal (photolysis etc)
    ! MAQ_lg_20230405+ found out, based on the simulations for the ATTO site, starting at 22.00 that this calculation of nstepinit 
	!  resulted in wrong timing of calculated parameter DayReal used for the first-order estimates of the j-val values. This
	!  specific case for ihr(1) > 12  and < 24 has been outcomemmented!

    !IF (ihr(1).gt.12.AND.ihr(1).LT.24) THEN
	!  nstepinit=(24-ihr(1))*(3600./dtimestep)
	!ELSE 

	nstepinit=ihr(1)*(3600./dtimestep)

	!ENDIF
    ! MAQ_lg_20230405- 
    ! ESS_lg_20150310-
	
    ! ESS_lg_20130118+ calculation of initial Julian day
    Jday(1)=0
    DO i=1,imon(1)-1
      Jday(1)=Jday(1)+ndaymonth(i)
    ENDDO
    Jday(1)=Jday(1)+INT(iday(1))

	! ESS_lg_20130228+ added the calculation of the actual NDTG, also used to process the observations
    ndtgact(1)=ndtginit
	
    ! calculating the year, month, day and hour for each timestep
    nday=0
  	  
 
	! ESS_lg_20130118+ the parameter LDATLTIME to use for the time column in the output files
    ! MAQ_lg_20210330+ filling ldatltime for i=1 
    i=1
    isec(i)=INT(i*delta_time)
    ihr(i)=ihr(1)+INT(isec(i)/3600.)-INT(nday*24)
    imin(i)=INT(isec(i)/60.)-INT((ihr(i)-ihr(1))*60)-INT(nday*24*60) ! ESS_lg_20130415+ modified 

    ldatltime(i)=''
    ip=index(ldatltime(i),' ')
    write(dummy,'(i2.2,a1)')  &
     +    iday(i),'-'
    ldatltime(i)=ldatltime(i)(1:ip-1)//dummy
    ip=index(ldatltime(i),' ')
    write(dummy,'(i2.2,a1)')  &
     +    imon(i),'-'
    ldatltime(i)=ldatltime(i)(1:ip-1)//dummy
    ip=index(ldatltime(i),' ')
    write(dummy,'(i2.2,a1)')  &
     +    iyear(i),'-'
    ldatltime(i)=ldatltime(i)(1:ip-1)//dummy
    ip=index(ldatltime(i),' ')

    write(dummy,'(i2.2,a1)') &
     +    ihr(i),':'
    ldatltime(i)=ldatltime(i)(1:ip-1)//dummy
    ip=index(ldatltime(i),' ')
    write(dummy,'(i2.2)')  &
     +    imin(i)
    ldatltime(i)=ldatltime(i)(1:ip-1)//dummy
    ldatltime(i)(9:9)=' '        ! introduce the space between date and time
    ! MAQ_lg_20210330-

    DO i=2,Nstep
      isec(i)=INT(i*delta_time)
      ihr(i)=ihr(1)+INT(isec(i)/3600.)-INT(nday*24)
      imin(i)=INT(isec(i)/60.)-INT((ihr(i)-ihr(1))*60)-INT(nday*24*60) ! ESS_lg_20130415+ modified 
      iday(i)=iday(i-1)
      Jday(i)=Jday(i-1)
      imon(i)=imon(i-1)
      iyear(i)=iyear(i-1)

      ! MAQ_lg_20230412+ modification of ndaymonth for leapyear
	  IF (MOD(iyear(i),4) == 0) THEN
	     ndaymonth(2)=29
	  ELSE 
	     ndaymonth(2)=28
	  ENDIF
      ! MAQ_lg_20230412-

      timeday(i)=REAL(ihr(i))+REAL(imin(i))/100.

      ! Changing the day
      IF (ihr(i)>=24.AND.ihr(i-1)<24) THEN
        nday=nday+1
        iday(i)=iday(i)+1
        Jday(i)=Jday(i-1)+1
        ihr(i)=0
      ENDIF
      ! Changing the month
      IF (iday(i)>ndaymonth(imon(i))) THEN
         imon(i)=imon(i)+1
         iday(i)=1
      ENDIF
      ! Changing the year
      IF (Jday(i)>365) THEN
        Jday(i)=1
        imon(i)=1
        iyear(i)=iyear(i)+1
      ENDIF

      ! ESS_lg_20130228+ added the calculation of the actual NDTG, also used to process the observations
	  ! ESS_lg_20131124+ modified to secure proper calculation of NDTGACT
      ndtgact(i)=INT(INT(icent*1e8)+INT(iyear(i)*1e6)+INT(imon(i)*1e4)+ &
	             INT(iday(i)*1e2)+INT(ihr(i)))
				
  	  ! ESS_lg_20130118+ the parameter LDATLTIME to use for the time column in the output files
      ldatltime(i)=''
      ip=index(ldatltime(i),' ')
      write(dummy,'(i2.2,a1)')  &
     +    iday(i),'-'
      ldatltime(i)=ldatltime(i)(1:ip-1)//dummy
      ip=index(ldatltime(i),' ')
      write(dummy,'(i2.2,a1)')  &
     +    imon(i),'-'
      ldatltime(i)=ldatltime(i)(1:ip-1)//dummy
      ip=index(ldatltime(i),' ')
      write(dummy,'(i2.2,a1)')  &
     +    iyear(i),'-'
      ldatltime(i)=ldatltime(i)(1:ip-1)//dummy
      ip=index(ldatltime(i),' ')

      write(dummy,'(i2.2,a1)') &
     +    ihr(i),':'
      ldatltime(i)=ldatltime(i)(1:ip-1)//dummy
      ip=index(ldatltime(i),' ')
      write(dummy,'(i2.2)')  &
     +    imin(i)
      ldatltime(i)=ldatltime(i)(1:ip-1)//dummy
      ldatltime(i)(9:9)=' '        ! introduce the space between date and time
    ENDDO

    ! ESS_lg_20130118-
        
    ! ESS_lg_20120718+ initialization of location, surface pressure and surface cover parameters
    ! surface cover fractions
    latitude(:)     = zlatitude              ! latitude in degrees
	longitude(:)    = zlongitude             ! ESS_lg_20130817+ longitude in degrees
    press(:)        = zpress                 ! surface pressure [Pa]
    slf(:)          = zslf                   ! surface land fraction [0-1]
    slm(:)          = 1.                     ! land-sea mask, integer
    loland(:)       = .true.                 ! land-sea mask, logical 
    vgrat(:)        = zvegfrac               ! vegetation fraction
    forestfr(:)     = zforestfr              ! forest fraction

    ! Canopy structure parameters
    lai(:)            = zlai          ! LAI
    ! MAQ_lg_20170929+ modified definition/explanation of laibase
	laibase(:)        = lai(:)        ! set this value to e.g., a reference LAI of 3 when conducting a sensitivity analysis for LAI for 
                                      ! a range from 1-5	
	
    ! introduced definition of various Leaf Area Density (LAD) profiles for four equidistant canopy layers
	! ESS_lg_20140416+ modified, to deal with the possibility to conduct simulations with > 2 (nveglay)/4 (nveglay_hr) canopy layers
    iladprof          = ziladprof     ! 1 for uniform profile, iladprof=2: prescribed profile 
    IF (iladprof == 1) THEN
      lad(:,:)        = 1._dp/nveglay_hr ! uniform profile
      write(*,*)'Initialization LAD profile: uniform profile'
	  write(*,*)'ENTER to continue'
	  read (*,*)
	ELSEIF (iladprof == 2) THEN 
      write(*,*)'Initialization LAD profile: top profile, check emdep_xtsurf_vdiff for details'
	  write(*,*)'ENTER to continue'
	  read (*,*)
	  IF (nveglay_hr == 4) THEN         ! with ~75% in crownlayer 
        lad(:,1)        = 0.334         ! LAD first crown-layer
        lad(:,2)        = 0.447         ! LAD second crown-layer
        lad(:,3)        = 0.203         ! LAD first understory-layer
        lad(:,4)        = 0.016         ! LAD second understory-layer
	  ELSE IF (nveglay_hr == 6) THEN    ! with ~75% in crownlayer 
        lad(:,1)        = 0.25          ! LAD first crown-layer
        lad(:,2)        = 0.40          ! ...
        lad(:,3)        = 0.131         ! MAQ_lg_20170509+ with these numbers for the top layers the total LAD is the same as that for 4 layers
        lad(:,4)        = 0.073         ! ...
        lad(:,5)        = 0.073         ! ...
        lad(:,6)        = 0.073         ! LAD lowest understory-layer
      ElSE
        write(*,'(1a,i3,1a)') &
		  ' Top LAD profile is defined for nveglay_hr=4 and 6 but not for: ',nveglay_hr,' layers'
        write(*,'(1a,i3,1a)') &
		  ' Change nveglay_hr to 4/6 OR change ILADPROF to 1 OR add LAD profile for: ',nveglay_hr,' layers'
		write(*,*)'STOP called in messy_emdep_xtsurf'
	    STOP
      ENDIF		
    ENDIF
    ! ESS_lg_20140416-
	
    hc(:)           = zcanheight      ! canopy height
	zrefsl(:)       = zrheightsl      ! ESS_lg_20130408+ reference height surface layer
    z0m(:)          = zroughness      ! roughness length for momentum
    disp(:)         = (2./3.)*hc(:)   ! displacement height
    
	! ESS_lg_20140416+ definition of height of canopy layers, assuming equidistant layers inside the canopy
    DO i=1,nveglay
	  zrefcan(:,i)   = hc(:)/nveglay*(nveglay+1-i)-0.5*hc(:)/nveglay
      WRITE(*,'(1a,i2,1x,f4.1,1a)')' Reference height layer # ',i,zrefcan(1,i),' m'
	ENDDO

	DO i=1,nveglay_hr
	  zrefcan_hr(:,i)   = hc(:)/nveglay_hr*(nveglay_hr+1-i)-0.5*hc(:)/nveglay_hr
      WRITE(*,'(1a,i2,1x,f4.1,1a)')' Reference height layers for high resolution calc. # ',i,zrefcan_hr(1,i),' m'
	ENDDO
	! ESS_lg_20140416-

    ! surface layer properties   
    zref(:)         = hc(:)+zrefsl(:) ! ESS_lg_20130208+ surface layer reference height, ref. height above the canopy
    geopot_3d_sl(:) = zref(:)*g       ! ESS_lg_20130208+ surface layer geopotential height 
    geopot_3d(:)    = zref(:)*g       ! ESS_lg_20130208+ geopotential height 
    prhoa(:)        = 1.225_dp        ! density in kg m-3
    grvol(:)        = 2._dp*zref*1.e10_dp ! grid size of 100 by 100 km 
    grmass(:)       = grvol(:)*prhoa(:)   ! grid mass in kg
    pdp(:)          = 2._dp*zref(:)*prhoa(:)*g*1.e3_dp  ! ESS_lg_20130220+ zref(:)! density difference surface layer
    qm1(:)          = zqm1            ! surface layer moisture in g H2O g-1 air
    qs(:)           = zqm1            ! saturated surface layer moisture in g H2O g-1 air ! MAQ_lg_20170720+
    qsam(:)         = zqm1            ! saturated moisture at leaf-level in g H2O g-1 air ! MAQ_lg_20170720+
    ! mz_lg_20050719-
 
    ! ESS_lg_20120717+ some parameters used to introduce a diurnal cycle in input data
    globradmax=1000.   ! maximum global radiation
    albedo=zalbedo     ! typical albedo of surface cover

    ! a diurnal cycle in some basic micrometeorological drivers of atmosphere-biosphere exchanges
    DO i=1,Nstep

       ! Calculate radiation as a function of latitude, season and time of the day; 
       ! see also emdep_update_jval for these calculations as well as:
       ! http://education.gsfc.nasa.gov/experimental/all98invproject.site/pages/science-briefs/ed-stickler/ed-irradiance.html
       ! Day as a real value (at start of spring DayReal = 0):
       !DayReal = real(i-nstepinit)/real(nstepday)+Jday(1)  ! ESS_lg_20150310+ added nstepinit
	   DayReal = real(nstepinit+i)/real(nstepday)+Jday(1)  ! MAQ_lg_20211103+ bugfix, nstepinit+1 ! ESS_lg_20150310+ added nstepinit

       ! seasonal cycle, Solar declination angle:
       SoDecli = -Cancer * COS (2.*PI*(DayReal+10.)/365.25)
       ! to correct for latitude
       Latitu = latitude(i) * (pi/180.)                           
       zen(i)    =   ACOS(SIN (Latitu) * SIN (SoDecli) &        ! Zenith angle
   &                   +COS (Latitu) * COS (SoDecli) * &
   &                    COS ((15./360.)*(timeday(i)-12)*2.*PI))  
       cossza(i) =   COS(zen(i))                                ! cosine of zenith angle
       globrad(i)=   MAX(0.,globradmax*MAX(0.-dp,cossza(i)))    ! global radiation, and scaled with photon to consider lat, time and season
       netrad(i) =   globrad(i)*(1-albedo)                      ! net radiation calculated from the global radiation and albedo

       ! ESS_lg_20130118+ some scaling functions to introduce diurnal cycles in some input parameters,
       doffset=1. ! 1-hr delay in temperature, windspeed and Mixed layer increases
       scaling =   COS ((15./360)*(timeday(i)-(12+doffset))*2.*PI)  ! scaling function for (surface) temperature and windspeed
                                                                    ! the term doffset results in an delay in the response by doffset hr's

	   ! ESS_lg_20130120+ minimum temperature of tmin K and maximum of tmax, 
       ! with the selected numbers being used to consider the larger tmin and tmax 
       ! as a function of latitude (25.*COS(Latitu)) as well as the larger seasonal 
       ! cycles in higher latitudes (15.*SIN(Latitu)*SIN(SoDecli))
       tmin=275.+25.*COS(Latitu)+15.*SIN(Latitu)*SIN(SoDecli)
       tmax=275.+30.*COS(Latitu)+15.*SIN(Latitu)*SIN(SoDecli)
       ! temporal cycle in Tair with a dampened temporal cycle compared to Tsurf
       tamin=tmin+0.1*(tmax-tmin) ! MAQ_lg_20170413+ the amplitude was determined by a constant of 0.25; for ATTO this has been reduced to 0.1
       tamax=tmax-0.1*(tmax-tmin)
       ! temporal cycle in Tsoil with a dampened temporal cycle compared to Tsurf
       tslmin=tmin+0.1*(tmax-tmin)
       tslmax=tmax-0.1*(tmax-tmin)
       ! temporal cycle in u-v windspeed component
       umin=2.0
       umax=5.0
       vmin=0.0
       vmax=0.0

       Kh_obs_sl(i) = 0.   ! MAQ_AV_20200309+ added initialization of Kh_obs_sl
       Kh_obs_cl(i) = 0.   ! MAQ_AV_20200309+ added initialization of Kh_obs_cl
       zustveg_obs(i) = 0. ! MAQ_AV_20201005+ added initialization of zustveg_obs

	   ! MAQ_20130414+ temporal cycle in RH
	   rhavg=0.8
	   rhdif=0.15
       ! MAQ_20130414- 
	   
       tsurf(i)  =   tmin   +(tmax-tmin)*scaling      ! surface/skin/canopy temperature
       tair(i)   =   tamin  +(tamax-tamin)*scaling    ! air temperature
       tsoil(i)  =   tslmin +(tslmax-tslmin)*scaling  ! soil temperature

	   tcan(i,:) = tsurf(i) ! MAQ_lg_20170505+ added Tcanopy
	   
       ! mz_lg_20040715+ virtual temperature; temporarily set to tsurf-1
       tvir(i)   =   tsurf(i)-1. 

       ! wind speed reflected by u and v component
       ! MAQ_20130414+ modified to arrive at more realistic calculations of u and v 
       um1(i)    =   (umax+umin)/2+((umax-umin)/2)*scaling
       vm1(i)    =   (vmax+vmin)/2+((vmax-vmin)/2)*scaling

       ! MAQ_20170414+ included an estimate of the diurnal cycle in RH which is also important to
	   ! determine the diurnal cycle in wet skin fraction and, thus dry deposition
       rh(i)    =   rhavg - rhdif*scaling	   
       ! MAQ_20170414-
	   
       ! ESS_lg_20120928+ introduced a daily varying depth of the surface mixed layer resembling the BL        
       zrefmin         = zcanheight+zrheightsl ! ESS_lg_20130208+ minimum reference height surface/mixed layer
       MLHmax          = zMLHmax         ! maximum depth of surface/mixed layer
       MLH(i)          = MAX(2*(zrefmin-zcanheight)+zcanheight,MLHmax * scaling) ! diurnal cycle in mixed layer depth
       zref(i)         = (MLH(i)-zcanheight)/2 ! ESS_lg_20130104+ the reference height is in the middle of the mixed layer
       geopot_3d(i)    = zref(i)*g       ! geopotential height 
       prhoa(i)        = 1.225_dp        ! density in kg m-3
       grvol(i)        = 2._dp*zref(i)*1.e10_dp ! grid size of 100 by 100 km 
       grmass(i)       = grvol(i)*prhoa(i)   ! grid mass in kg
       pdp(i)          = 2._dp*zref(i)*prhoa(i)*g*1.e3_dp  ! ESS_lg_20130220+ removed zcanheight! density difference surface layer

       ! calculated atmospheric stability from temperature gradients and wind speed
       ril(i)=(g*(tair(i)-tsurf(i))*zrefsl(i))/(((tair(i)+tsurf(i))/2.)*  &  ! ESS_lg_20131125; reference height surface layer!
                MAX(1e-10,um1(i)**2+vm1(i)**2))

    END DO

    ! some other general micromet input data
    ws(:)       =   zws ! soil moisture [m]
    wsMAX(:)    =zwsmax ! field capacity [m]
	
    ! prepare specific input data for biogenic VOC emissions
    VOC_emfact(:,:)= 0._dp
    fAPIN          = 0._dp      
    fBPIN          = 0._dp      
    fSQTERP        = 0._dp
    ! MAQ_lg_20170413+ added
	fATERP         = 0._dp      
    fLIMO          = 0._dp      
    fMYRC          = 0._dp
    ! MAQ_lg_20170413-	
	
    IF (l_emis_bio_VOC) THEN
        dm(:)               = lai(:)*125      ! foliar density [g m-2], =LAI*specific leaf weight
        dmbase(:)           = laibase(:)*125. ! MAQ_lg_20160817+ use the long-term average LAI for recalculation of VOC_EMFACT
        VOC_emfact(:,iisop) = zisopemfact ! isoprene emission factor [ug C g-1 hr-1], 
        VOC_emfact(:,imono) = zmonoemfact ! monoterpene emission factor [ug C g-1 hr-1]
        VOC_emfact(:,iovoc) = zovocemfact ! OVOC's emission factor [ug C g-1 hr-1]
        fAPIN               = zfAPIN      ! fraction of monoterpene emission emitted as a-pinene
        fBPIN               = zfBPIN      ! fraction of monoterpene emission emitted as b-pinene
        fSQTERP             = zfSQTERP    ! fraction of monoterpene emission emitted as sesquiterpenes

        ! MAQ_lg_20170413+
        fATERP              = zfAPIN*0.6  ! fraction of monoterpene emission emitted as a-terpinene
        fLIMO               = zfAPIN*5.   ! fraction of monoterpene emission emitted as limonene
        fMYRC               = zfAPIN*1.   ! fraction of monoterpene emission emitted as myrcene
        ! MAQ_lg_20170413-
		
	ENDIF

    ! prepare specific input data for soil biogenic NOx emissions
    IF (l_emis_bio_NO) THEN
        init_step   =     0 ! initial timestep
        ndaylen     = 86400 ! length of day in seconds
        imonth      = imon(1) ! month

        itrop=1             ! ESS_lg_20130118+ bug fix; defining itrop

        cultiv(:)   = zcult ! cultivation intensity [old: 0-15, new 2004: 0-1]
        fertil(:)   = zfert ! fertilizer application (synthetic/manure)
        prc(:)      =  zprc ! convective rainfall [m]
        prl(:)      =  zprl ! large-scale rainfall [m]
        prectot(:)  =     0.! monthly accumulated precipitation [m]
        noemis_w(:) =     0.! mz_lg_20050719+
        noemis_d(:) =     0.! mz_lg_20050719+
        noemclass1(:,:)=  0.! mz_lg_20050719+

        IF (ncl_noemis.eq.12) THEN
           
           ! ESS_lg_2013017+ modified to assign the read-in NO emission class
           NOemclassname(1)='Water'
           NOemclassname(2)='-'
           NOemclassname(3)='-'
           NOemclassname(4)='-'
           NOemclassname(5)='Tundra'
           NOemclassname(6)='Grass/shrub'
           NOemclassname(7)='Woodland'
           NOemclassname(8)='Decid.-forest'
           NOemclassname(9)='Conif.-forest'
           NOemclassname(10)='Dry-decid.-forest'
           NOemclassname(11)='Rainforest'
           NOemclassname(12)='Agric.-land'
           iNOemcl=iNOemclass         ! NO emission class

           itrop=11 ! ESS_lg_20130118+ bug fix; defining itrop

           ! mz_lg_20040716+ assigning the fractional coverage to each NO emission class.
           !     Normally this info is being read in from a preprocessed input file but 
           !     here in the boxmodel it a coverage of 1 is assigned to the selected 
           !     emission class (info in papers by Yienger and Levy, 1995, 
           !     1=Water, 2=- , 3=- , 4=-, 5=Tundra, 6=Grass/shrub,
           !     7=Woodland, 8=Decid. forest, 9=Conif. forest, 10=Dry-decid. forest
           !     11=Rain forest, 12=agric. land)

           NOemclass1(:,iNOemcl)= 1. ! fractional coverage of 1 for the selected NO emission class
           ! ESS_lg_20130107-

           DO i=1,ncl_noemis
              ! mz_lg_20040407+ assigning the emission factors
              noemis_w(:)=noemis_w(:)+noemfact_wet(i)*noemclass1(:,i)
              noemis_d(:)=noemis_d(:)+noemfact_dry(i)*noemclass1(:,i)
           ENDDO

        ENDIF
        lstart      =  .TRUE. ! switch for starting simulation    
    ENDIF

    ! assigning the radon emission flux; a typical flux for tropical soils (more recent compared to
	! the initially used emanation rate of 3000 atoms m-2 s-1 by Trumbore et al.) is 30 m Bq m-2 s-1 
	! this resembles 1.43 atoms cm-2 s-1, this 14300 atoms m-2 s-1.  
 
    RADON_slflux(:)=zradonemis ! [atoms m-2 s-1]

    ! ESS_lg_20130503+ assigning the CO2 respiration/emission flux
	CO2_slflux(:)=zco2emis*avo
	
	! MAQ_lg_20230329+ assigning a soil NH3 emission flux
	NH3_slflux(:)=znh3emis*avo
	
    ! mz_lg_20050719+ definition of a large selection of parameters normally
    !    being used in the dry deposition calculations. However, since these
    !    are also mostly used for the calculations of atmosphere-biosphere
    !    exchanges, e.g., tracer properties and turbulent exchanges the whole
    !    set of parameters is initialized, including some extra settings for
    !    the subroutine emdep_xtsurf_veg_mlay

    cfml(:)         = 1000.           ! momemtum drag coefficient over tropical forest
    cfncl(:)        = 2500.
    cdnl(:)         = (1./LOG(100./z0m(:)))**2   ! neutral drag coefficient
    z0mslsn(:)      = 0.005           ! snow/ice/bare surface roughness of 0.5 cm 
                                      ! (might be smaller, first order estimate)
    ! others
    soilph(:,:)     = 1./7.           ! similar contribution by 7 pH classes

	! ESS_lg_20150618+ calculation or assignment of wet skin fraction
    zcvw(:)=zwetskin ! assigning constant wet skin fraction
	! ESS_lg_20150618-

    ! ESS_lg_20130208+ tracer initialization and reading of input is moved to here to secure that the read data 
	!    are used to calculate the stomatal resistance, soil moisture stress function and J-values

    ! ESS_lg_20120722+ tracer initialization dependent on the call of the chemistry scheme
    PRINT *,'Tracer initialization in emdep_xtsurf_tracer_init'

    call emdep_xtsurf_tracer_init( l_xtsurf_veg_mlay_chem,                   &
          Nstep, ntrac, trname, moleweight, reactivity, henrycoeff, medium,  &
          lo_derived, lexist_O3, lvd_bigl, lexist_GAS,                       &
          latmbios, latmbios_emveg, latmbios_emsoil, latmbios_photo,         &
          latmbios_output, rj_max, pxtm1, pxtmveg)                             ! ESS_lg_20140418+

    ! ESS_lg_20120722-

        ! ESS_lg_20120722+ having the option to read in observations to use these to constrain
    !    the model simulations with these data and to use the model to support analysis of
    !    these observations

    pxtm1_obs(:,:)=-9999.9999  ! give default a value to NA tracer mixing ratios

    IF (l_readdata) THEN

       PRINT *,'Reading a file with observations to constrain the simulations'
       PRINT *,'Check carefully for the assignment of these observations to model parameters'

       call emdep_xtsurf_readdata(infilename,Nstep,nprint,delta_time,ndtgact,ldatltime,observ)  ! ESS_lg_20130228+ ndtgact!

       ! assigning the read-in observations to the specific model input data

       DO i=1,Nstep

	     !==================================================================================================
         ! Note that the code below is specific to assigning the read-in the observed surface layer 
         ! NO, NO2, O3 and CO mixing ratios and precipitation from the example input file. The observed 
         ! surface layer mixing ratios are assigned to the parameter pxtm1_obs which is then used in the model
         ! to "nudge" the calculate the surface layer mixing ratio (pxtm1) (see Ganzeveld et al., 2006) on the
         ! nudging in such stand-alone/1-D model approaches. This nudging of the tracer concentrations/
         ! mixing ratios is done to simulate chemical conditions close to the observed conditions, 
         ! as such implicitly adding the role of advection of air masses with different composition
         ! =================================================================================================
           
         ! only assigned when not NA but also >= 0 (not using small negative concentrations)
         IF (casename.EQ.'example_run') THEN
           IF (observ(i,1) > 0.) pxtm1_obs(i,idt_NO)=1e-9*observ(i,1)  ! only assigning NO when NO > 0
           IF (observ(i,2) > 0.) pxtm1_obs(i,idt_NO2)=1e-9*observ(i,2) ! only assigning NO2 when NO2 > 0

           ! only assigning a value when not NA
           IF (observ(i,3) > -9999.) pxtm1_obs(i,idt_O3)=1e-9*observ(i,3)  ! observed O3 is given in ppbv, only assigned when not NA
           IF (observ(i,4) > -9999.) pxtm1_obs(i,idt_CO)=1e-9*observ(i,4)  ! observed CO is given in ppbv, only assigned when not NA

           ! observed precipitation
           IF (observ(i,6) > -9999.) prc(i)=1e-3*observ(i,6) ! assigning the read-in precipitation from the example file
                                                           ! observ(:,6) is in mm whereas prc is in 
         ENDIF

         ! MAQ_lg_20230506+ added another read-in dataset showing a more complete input dataset, including more details
		 !    micromet measurements besides measurements on tracer concentrations. This dataset is that of Bosco Fontana, 
		 !    which has been used as input for a study by Auke Visser on the impact of constraining such a canopy exchange 
         !    model with the observation inferred eddy diffusivities and also also checking on the role of high soil NO 
         !    emissions in O3 deposition. In order to use this example dataset, change the emdep.nml file to that reflecting
		 !    the use of MLC-CHEM for this Bosco Fontana case study
		 
	     ! MAQ_AV_20190913+ added Bosco Fontana (data courtesy of Giacomo Gerosa, UniCatt, ITA, see Finco et al., 2018, ACP)
	     IF (casename.EQ.'BoscoFontana') THEN
	 
	       IF (observ(i,2) > -999.99) netrad(i) = MAX(0.,observ(i,2))
	       IF (observ(i,6) > -999.99) tair(i) = observ(i,6) + 273.16 !32 m
	       IF (observ(i,8) > -999.99) prc(i) = 1.e-3 * observ(i,8)
	       IF (observ(i,12) > -999.99) rh(i) = observ(i,12)/100. !32 m
	    
	       !Derive Tsurf from Tair at lowest measured level:
	       IF (observ(i,5) > -999.99 .AND. netrad(i).gt.0.) THEN
	         tsurf(i) = observ(i,5) + 273.16 + netrad(i) * 0.0022 !Tsurf > Tair during daytime
	       ELSEIF (observ(i,3) > 4000. .AND. netrad(i).lt.1e-5) THEN !MAQ_AV_20200407: note the bug in this line: 4000 should be -999.99! Test effect later!
	         tsurf(i) = observ(i,3) + 273.16 -1. !stable conditions at night
	       ENDIF
           
		   ! MAQ_lg_20210401+ added estimate of the Tsoil being set at a constant value of 300K
		   tsoil(i)=300.
	    
	       IF (observ(i,21) > -999.99) um1(i) = observ(i,21)
	       IF (observ(i,29) > -999.99) pxtm1_obs(i,idt_O3) = observ(i,29) * 1.e-9
	    
	       !MAQ_AV_20200407+ added observed NO and NO2 concentrations
	       IF (observ(i,35) > -999.99) pxtm1_obs(i,idt_NO2) = observ(i,35) * 1.e-9
	       IF (observ(i,36) > -999.99) pxtm1_obs(i,idt_NO) = observ(i,36) * 1.e-9

	       !MAQ_AV_20200306+ added observation-derived K_h values
	       IF (observ(i,38) > -999.99) Kh_obs_sl(i) = observ(i,38)
	       IF (observ(i,39) > -999.99) Kh_obs_cl(i) = observ(i,39)
	       !MAQ_AV_20200306-
	    
           ! ESS_GK_20150306+ calculating the soil CO2 respiration flux from the soil temperature
	       CO2_slflux(i)=1e-6*0.4281*exp(0.1369*(tsoil(i)-273.15))*avo ! molecules m-2 s-1
	       ! ESS_GK_20150306-

   	     ENDIF
         ! MAQ_lg_20230506-

		 ! ======================================================================================================
         ! the next lines of code can be activated (remvoving the "!" in front of the lines: IF (observ(i,*)...) 
         ! to use also the observed net radiation, temperatures, wind speed and moisture whenever available in the 
         ! input file; to activate assign the seocond index of the array 'observ' that resembles the column 
         ! from the input file with the respective param.) and be sure that the that units are correctly defined

         ! observed net radiation (global radiation corrected for albedo)
         ! ESS_lg_20130228+ introduced a system that secures that the J-value estimates are corrected
         ! for the difference between the actual net radiation (clouds) and the estimated max. net 
         ! radiation (no-clouds)
         scaling_rad(i)=1.

         !IF (observ(i,8) > -9999.) THEN
		    !netrad(i)=observ(i,8)*(1.-albedo)  ! observ(i,8) = Downward irradiance
            ! -------------------------------------------------------------------------------
            ! ESS_lg_20130228+ whenever the observed net radiation is being used then the
			! radiation dependent parameters should be corrected for the difference between the
			! maximum (cloud-free) net radiation and the actual net radiation)
            ! -------------------------------------------------------------------------------
			!netradmax=globrad(i)*(1.-albedo)
            !IF (netradmax>0.) scaling_rad(i)=MIN(1.,MAX(0.,netrad(i)/netradmax))
		 !ENDIF
         ! ESS_lg_20130228- 

         ! observed air temperature (temperature measured in atmospheric surface layer) 
         !IF (observ(i,*) > -9999.) tair(i)=observ(i,*) ! assigning the read-in air temperature [K]
                 
         ! observed surface/skin/canopy temperature (temperature measured at the top of the canopy) 
         !IF (observ(i,*) > -9999.) tsurf(i)=observ(i,*) ! assigning the read-in surface temperature [K]
                 
         ! observed soil temperature (temperature measured in soil, e.g., in the top few cm of the soil) 
         !IF (observ(i,*) > -9999.) tsoil(i)=observ(i,*) ! assigning the read-in soil temperature [K]

         ! observed windspeed (U-component, measured in atmospheric surface layer); when only the windspeed is known
         ! not having the individual U (x- or West-East direction) and V (y- or North-South direction component) then
         ! assign the observed windspeed to U and set vm1 to 0!           
         !IF (observ(i,*) > -9999.) um1(i)=observ(i,*) ! assigning the read-in windspeed, U-component [m s-1]

         ! observed windspeed (V-component, measured in atmospheric surface layer) 
         !IF (observ(i,*) > -9999.) vm1(i)=observ(i,*) ! assigning the read-in windspeed, V-component [m s-1]

         ! observed moisture, measured in atmospheric surface layer) 
         !IF (observ(i,*) > -9999.) qm1(i)=observ(i,*) ! assigning the read-in surface moisture g H2O g-1 air
         
		 ! ESS_lg_20130503+ temporarily fixing the surface layer CO2 concentrations
		 pxtm1_obs(i,idt_CO2)=pxtm1(1,idt_CO2)

         ! ESS_lg_20131125+ calculate the surface layer wind speed from um and vm
		 u_sl(i)=(um1(i)**2+vm1(i)**2)**0.5

         ! update the calculated atmospheric stability from the observed temperature gradients and windspeed
         ril(i)=(g*(tair(i)-tsurf(i))*zrefsl(i))/(((tair(i)+tsurf(i))/2.)*  &
                MAX(1e-10,um1(i)**2+vm1(i)**2))
 
		 ! ESS_lg_20130208+ and update the virtual temperature used in the calculation of Monin-obukhov lenght!
         ! NOTE that these calculations should be still updated using the actual formula to calculate Tvirtual
		 tvir(i)=tair(i)
		 tvl(i)=tsurf(i)
         ! ESS_lg_20130208-

	  ENDDO
       
    ENDIF
    ! ESS_lg_20120722-
    ! ESS_lg_20130208-

    ! ESS_lg_20150618+ moved calculation of the surface cover fractions after the definition of of zcvw
    !--- Calculate bare soil fraction:
    !    It is calculated as residual term over land
    DO i=1,Nstep
      zfrl(i)=slm(i)

	  ! ESS_lg_20150618+ calculation or assignment of wet skin fraction
      !   This is now simply based on a direct scaling using the RH relative to a minimum RH;
	  !   This should be replaced by the parameterization proposed by G. Lammel (Linda Voss, personal communications)
      IF (l_wetskinRH) THEN
        ! ESS_lg_20150618+ see table A2 in Lammel et al., 1999, report # 286 of MPI-M, on the formation of nitrous
		! acid: Parameterizations and comparison with observations
		! calculation of the wet skin fraction from the formula on SAIdry-vegetation
		! zcvw=1.-(SAIdry vegetation/SAIvegetation)
	    IF (rh(i) >= 0.55 .AND. rh(i) < 0.9) THEN
		  zcvw(i)=1.-(1.-(rh(i)-0.55)/0.35)
	    ELSEIF (rh(i) >=0.9) THEN
		  zcvw(i)=1.
		ELSE
		  zcvw(i)=0.
		ENDIF
	  ENDIF

	  ! ESS_lg_20150618-

      IF (zfrl(i) .GT. 0.) THEN
        zvgrat(i)=(1.-zcvw(i))*vgrat(i)     ! MAQ_20170208+ bug fix to avoid having a total surface cover fraction > 1
        zcvbs(i)=(1.-zcvw(i))*(1.-vgrat(i)) ! bare soil fraction
      ELSE
        zvgrat(i)=0.
        zcvbs(i)=0.
      ENDIF

 	ENDDO
    ! ESS_lg_20150618-
	
    DO i=1,Nstep

	   ! -------------------------------------------------------------------------------
       ! ESS_lg_20150311+ determining the scaling function to estimate the non-observed
	   ! radiation dependent terms (photolysis); netradmax, scaled with the cos zenith angle 
	   ! and the actual netrad are used to calculate this radiation scaling function
       ! -------------------------------------------------------------------------------

	   IF (netradmax*MAX(0.,COS ((15./360.)*(timeday(i)-12)*2.*PI))>0.) &
	     scaling_rad(i)= &
		   MIN(1.,MAX(0.,netrad(i)/MAX(0.,netradmax*MAX(0.,COS ((15./360.)*(timeday(i)-12)*2.*PI)))))

	   ! calculation of soil moisture stress function from the actual and maximum soil moisture
       zwcr=0.75*wsmax(i)   
       zwpwp=0.35*wsmax(i)
       fws(i)=MAX(0.,MIN(1.,(ws(i)-zwpwp)/(zwcr-zwpwp)))  ! ESS_lg_20130127+ to avoid Fws > 1, soil moisture stress function

       ! mz_lg_20040720+ including the assignment/calculation of parameters 
       !    needed to calculate the leaf stomatal resistance according to 
       !    Sellers, 1986, as included in ECHAM. This parameter can be replaced
       !    by any alternative representation of stomatal exchanges. 
       zva=cva
       zvb=cvb
       zvc=cvc
       zvbc=cvbc
       zvk=cvk
       zvkc=cvkc
       zvabc=cvabc
       zvrad=cvrad
       zepsr=1.e-10
       znetrad=MAX(zepsr,netrad(i)*zvrad)
       zabcs=(zva+zvbc)/(zvc*znetrad)
       ! calculation of leaf stomatal resistance using the echam
       ! equation applying an LAI of 1, included by Laurens Ganzeveld, 18/10/01 

       IF (znetrad > 1e-10) THEN
         rco_leaf(i,:)=1./((zvb*(LOG((zabcs*EXP(zvk)+1.)/(zabcs+1.))) & ! ESS_lg_20130516+ layer dependent, not yet correct for rad. extinction
               /zvabc-(LOG((zabcs+EXP(-zvk))/(zabcs+1.))))/zvkc)
       ELSE
         rco_leaf(i,:)=1e10
       ENDIF

       ! -----------------------------------------------------------------

       ! mz_lg_20040721+ and calculation of the surface roughness from the 
       !    specific drag coefficients over each surface cover and the fractions,
       !    assuming a blending height of 100 m

       zcdragv=(1./LOG(100./z0m(i)))**2
       zcdragslsn=(1./LOG(100./z0mslsn(i)))**2
       zcdrag=(zvgrat(i)+zcvw(i))*zcdragv+            &  ! drag over (wet/dry) vegetation 
              (zcvbs(i))+zcdragslsn                      ! drag over bare soil
       az0(i)=100./EXP(1._dp/SQRT(MAX(0.001_dp,zcdrag)))

    END DO

    ! calculate diurnal cycle in photolysis rates
    rj(:,:)=0._dp
    IF (l_xtsurf_veg_mlay.AND.l_xtsurf_veg_mlay_chem.AND.l_xtsurf_veg_mlay_chem_photolysis) THEN
      DO i=1,Nstep

        ! estimating J values based on assigned maximum values (see emdep_xtsurf_tracer_init)
        ! and scaling functions for time and latitude
        DO jt=1,ntrac
          IF (rj_MAX(jt) > 0._dp) rj(i,jt)=rj_MAX(jt)*  &
             MAX(0.,COS(REAL(0.5*pi*latitude(1)/90.)))* &     ! scaling with latitude
             MAX(0.,SIN(-0.5*Pi+REAL(2*Pi*i)/nstepday))       ! scaling with time
        ENDDO

        ! call of subroutine also used in MECCA's (Rolf Sander) box model
        CALL emdep_update_jval(i, nstepday, nstepinit, Jday(i), latitude(i), rj) ! ESS_lg_20150310+ nstepinit
		
        rj(i,:)=rj(i,:)*scaling_rad(i)  ! ESS_lg_20130228+ correcting the rj value estimates for "cloud cover"

      ENDDO 
    ENDIF

    ! ==========================================================================================
    ! ESS_lg_20120721+ actual start of simulations by calling the different subroutines

    ! mz_lg_20040716+ biogenic VOC emissions
    IF (l_emis_bio_VOC.OR.l_xtsurf_veg_mlay_chem.OR.l_xtsurf_AGS) THEN

        PRINT *,'Calculation of canopy radiation profiles (emdep_xtsurf_calcprof)'

        ! mz_LG_20020115 calling of subroutine in which the vertical profiles
        !     of radiation and windspeed in the canopy are calculated.
        !     The radiation profiles are used to calculate the VOC emission
        !     fluxes and are also used to calculate the photodissociation
        !     rates in the multi-layer vegetation model.

        ! MAQ_lg_20160621+
        CALL emdep_xtsurf_calcprof( Nstep, & ! ESS_lg_20120722+
           nveglay_hr, globrad, cossza, laibase, lad,   &! MAQ_lg_20201016+ globrad instead of netrad
           rbvd, rvdsl, rvd, fslbase)  ! mz_lg_20050721+ added SL diffusive rad.
        ! MAQ_lg_20160621-
		   
        CALL emdep_xtsurf_calcprof( Nstep, & ! ESS_lg_20120722+
           nveglay_hr, globrad, cossza, lai, lad,   & ! MAQ_lg_20201016+ globrad instead of netrad
           rbvd, rvdsl, rvd, fsl)  ! mz_lg_20050721+ added SL diffusive rad.

        ! show some results; canopy radiation parameters
        OPEN(unit=2,file='output/Radiation_canopy.out',status='UNKNOWN')

        ! ESS_20140416+ modified writing of radiation/canopy structure data for higher vertical resolution model 
        WRITE(2,'(1a)') &
           'Canopy radiation/structure data for high resolution of number of vegetation layers (21 max)'
        WRITE(2,'(1a)') &
           'LAD profile, glob. rad., diffusive beam comp., frac. sunlit leaves (FSL) per layer and sum, and scaled with LAI'  ! MAQ_lg_20201016+ globrad instead of netrad, ! MAQ_20180112+
		   
		! to assign parname and parunit for the header of the file
		! MAQ_20180112+ added number of output parameters
		ALLOCATE (parname(4,nveglay_hr))
		ALLOCATE (parunit(4,nveglay_hr))
 
        ! assigning some of this header info and some further manipulations for output
        fslsum(:)=0.
        DO i=1,nveglay_hr

          ! Leaf Area Density profile
		  dummy='LAD('
          iip=INDEX(dummy,' ')
          IF (zrefcan_hr(1,i) >= 10.) THEN
            write(dummy2,'(f4.1)')zrefcan_hr(1,i)
	      ELSE
            write(dummy2,'(f3.1)')zrefcan_hr(1,i) 
	      ENDIF 
          dummy=dummy(1:iip-1)//dummy2
          iip=INDEX(dummy,' ')
 		  dummy3=')'
          parname(1,i)=dummy(1:iip-1)//dummy3
          parunit(1,i)='fraction'

          ! Diffusive radiation profile
		  dummy='RVD('
          iip=INDEX(dummy,' ')
          IF (zrefcan_hr(1,i) >= 10.) THEN
            write(dummy2,'(f4.1)')zrefcan_hr(1,i)
	      ELSE
            write(dummy2,'(f3.1)')zrefcan_hr(1,i) 
	      ENDIF 
          dummy=dummy(1:iip-1)//dummy2
          iip=INDEX(dummy,' ')
 		  dummy3=')'
          parname(2,i)=dummy(1:iip-1)//dummy3
          parunit(2,i)='W m-2'
          
          ! MAQ_20180112+ added output
          ! PAR for each layer
		  dummy='PAR('
          iip=INDEX(dummy,' ')
          IF (zrefcan_hr(1,i) >= 10.) THEN
            write(dummy2,'(f4.1)')zrefcan_hr(1,i)
	      ELSE
            write(dummy2,'(f3.1)')zrefcan_hr(1,i) 
	      ENDIF 
          dummy=dummy(1:iip-1)//dummy2
          iip=INDEX(dummy,' ')
 		  dummy3=')'
          parname(3,i)=dummy(1:iip-1)//dummy3
          parunit(3,i)='umol m-2'

          ! fraction of sunlit leaves per layer
		  dummy='FSL('
          iip=INDEX(dummy,' ')
          IF (zrefcan_hr(1,i) >= 10.) THEN
            write(dummy2,'(f4.1)')zrefcan_hr(1,i)
	      ELSE
            write(dummy2,'(f3.1)')zrefcan_hr(1,i) 
	      ENDIF 
          dummy=dummy(1:iip-1)//dummy2
          iip=INDEX(dummy,' ')
 		  dummy3=')'
          parname(4,i)=dummy(1:iip-1)//dummy3
          parunit(4,i)='[-]'
		  ! calculation of PAR expressed in umol m-2 per layer
		  PAR(:,i)=4.405*(rbvd(:)*fsl(:,i)+rvd(:,i)*(1-fsl(:,i)))

          ! MAQ_20180112-

		  ! calculation of integrated fraction of sunlit leaves
          fslsum(:)=fslsum(:)+fsl(:,i)*lad(:,i)		  
		ENDDO

        ! MAQ_20180112+ added output parameters
		WRITE(2,'(A10,A16,50A15)') &
           'istep','ldatltime',(parname(1,i),i=1,nveglay_hr),'Glob. rad.',& ! MAQ_lg_20201016+ globrad instead of netrad
		                       (parname(2,i),i=1,nveglay_hr),             &
							   (parname(3,i),i=1,nveglay_hr),             &
							   (parname(4,i),i=1,nveglay_hr),             &
							       'FSL-sum','FSL*LAI'
        WRITE(2,'(A10,A16,50A15)') &
           '#','dd-mm-yy hh:min',(parunit(1,i),i=1,nveglay_hr),'W m-2', &
		                         (parunit(2,i),i=1,nveglay_hr),         &
								 (parunit(3,i),i=1,nveglay_hr),         &
								 (parunit(4,i),i=1,nveglay_hr),         &
								    '[-]','[m2 m-2]'
        DO i=1,Nstep
          IF (mod(i,nprint) == 0) THEN
            WRITE(2,'(I10,2x,A14,50(1x,E14.6))') &
                i,ldatltime(i),(LAD(i,ii),ii=1,nveglay_hr),globrad(i),& ! MAQ_lg_20201016+ globrad instead of netrad
				               (rvd(i,ii),ii=1,nveglay_hr),           &
							   (PAR(i,ii),ii=1,nveglay_hr),           &
							   (fsl(i,ii),ii=1,nveglay_hr),           &
							       fslsum(i),fslsum(i)*lai(i)
          ENDIF
        END DO
        ! MAQ_20180112-
		
        DEALLOCATE (parname)
		DEALLOCATE (parunit)
        ! ESS_20140416-

        CLOSE(unit=2)

    ENDIF

    ! mz_lg_20040716+ biogenic VOC emissions
    IF (l_emis_bio_VOC) THEN

        PRINT *,'Calculation of biogenic VOC emissions (emdep_emis_bio_VOC***) based on: '

        ! mz_LG_20020115 calling of subroutine in which the biogenic VOC emissions
        !     are calculated using a modified version of the Guenther et al. 1995
        !     algorithm or MEGAN

        IF (.NOT.l_emis_bio_VOC_MEGAN) THEN               ! ESS_lg_20130817+ G95 or MEGAN

          PRINT *,'Guenther et al., 1995'
          ! MAQ_lg_20160817+ added to compare the emission factor with that of MEGAN
          WRITE(*,'(1a,f8.1)') &
            ' The G95 canopy-scale C5H8 emission factor [ug C m-2 hr-1] is: ', &
              VOC_emfact(1,iisop)*DMBASE(1) ! MAQ_lg_20160817+ use the long-term average LAI!
	      print *,'ENTER TO CONTINUE'
	      read (*,*)
          ! MAQ_lg_20160817-		  

		  CALL emdep_emis_bio_VOC ( Nstep,              & ! ESS_lg_20120722+
            nveglay_hr, nveglay, iISOP, iMONO, iOVOC,   &
            dm, lad, VOC_emfact, tsurf, tcan,           & ! MAQ_lg_20170505+ added Tcanopy
            rbvd, rvd, fsl, VOC_emflux,                 & ! mz_lg_20040423+
            ISOP_emflux, MONO_emflux, OVOC_emflux)

		ELSE  ! ESS_lg_20130817+ calling of subroutine in which the biogenic VOC emissions
              !     are calculated according to MEGAN algorithm

          PRINT *,'MEGAN'
			  
		  CALL emdep_emis_bio_VOC_MEGAN( Nstep,           &
		    delta_time, latitude, longitude,              & ! ESS_lg_20120722+
            nveglay_hr, nveglay, iISOP, iMONO, iOVOC,     &
            dm, lai, lad, VOC_emfact, tsurf, tcan,        & ! MAQ_lg_20170505+ added Tcanopy
            rbvd, rvd, fsl, fslbase, dmbase,              & ! MAQ_lg_20160817+ dmbase MAQ_20160621+ added fslbase! mz_lg_20040423+
            fws, VOC_emflux, ISOP_emflux, MONO_emflux, OVOC_emflux) ! MAQ_20160621+ added fws  

		ENDIF                                             ! ESS_lg_20130817-
		  
        ! ESS_20140416+ modified writing of radiation data for higher vertical resolution model 
	    ! show some results; VOC emissions
        OPEN(unit=2,file='output/VOCemis.out',status='UNKNOWN')
        WRITE(2,'(1a)') &
           'VOC emissions calculated according to Guenther et al., 1995; isoprene emission fluxes for nveglay'

 		! to assign parname and parunit for the header of the file
		ALLOCATE (parname(1,nveglay))
		ALLOCATE (parunit(1,nveglay))
 
        ! assigning some of this header info and some further manipulations for output
        DO i=1,nveglay
		  dummy='ISOPfl('
          iip=INDEX(dummy,' ')
          IF (zrefcan(1,i) >= 10.) THEN
            write(dummy2,'(f4.1)')zrefcan(1,i)
	      ELSE
            write(dummy2,'(f3.1)')zrefcan(1,i) 
	      ENDIF
          dummy=dummy(1:iip-1)//dummy2
          iip=INDEX(dummy,' ')
 		  dummy3=')'
          parname(1,i)=dummy(1:iip-1)//dummy3
          parunit(1,i)='molec m-2 s-1'
        ENDDO

		WRITE(2,'(1a)') &
           'Net. rad., Tsurf, C5H8: all canopy layers, C10H16 and OVOCs only the top-layer emissions'
        WRITE(2,'(A10,A16,50A15)') &
           'istep','ldatltime','Net. rad.','Tsurf',(parname(1,i),i=1,nveglay),'MONOflux','OVOCflux'
        WRITE(2,'(A10,A16,50A15)') &
           '#','dd-mm-yy hh:min','W m-2','K',(parunit(1,i),i=1,nveglay),'molec m-2 s-1','molec m-2 s-1'
        DO i=1,Nstep
          IF (mod(i,nprint) == 0) THEN
            WRITE(2,'(I10,2x,A14,50(1x,E14.6))')      &
                i,ldatltime(i),netrad(i),tsurf(i),    &
               (ISOP_emflux(i,ii),ii=1,nveglay),      &
                MONO_emflux(i,1),OVOC_emflux(i,1)
          ENDIF
        END DO

	    DEALLOCATE (parname)
		DEALLOCATE (parunit)
        ! ESS_20140416-

        CLOSE(unit=2)				


	ELSE
        PRINT *,'Calculation of biogenic VOC emissions switched off'
    END IF

    ! mz_lg_20040716+ soil-biogenic NO emissions
    IF (l_emis_bio_NO) THEN

        PRINT *,'Calculation of soil-biogenic NO emissions (emdep_emis_bio_NO)'

        CALL emdep_emis_bio_NO(                                              &
           latitude, Nstep, init_step, ndaylen, delta_time, imonth,          & ! mz_lg_20040921+
           ncl_noemis, itrop, lstart, lcrfyl95, l_xtsurf_veg_mlay,           &
           l_emis_bio_NO_pls,                                                & ! ESS_lg_20120717+ removed lemis_bio_NO_mm
           cultiv,  fertil,                    &
           tsoil,   ws,                        &
           prc,     prl,                       &
           prectot,                            &
           noemis_w,noemis_d,                  &
           noemclass1, iNOemcl,                &
           lai,     slf,                       & ! mz_lg_20050522+ added slf
           cpold,                              &
           lspold,                             &
           pulsing, plsday,                    &        
           plsdurat,cp,                        &
           lsp,     pls,                       &
           crf,     NO_slflux,                 & ! mz_lg_20050719+, NO_slflux, mz_lg_20050522+ added CRF
           NO_emflux) ! mz_lg_20040426+

           ! show some results; soil NOx emissions
           OPEN(unit=2,file='output/soilNOxemis.out',status='UNKNOWN')
           WRITE(2,'(1a)') &
             'Soil NOx emissions calculated according to an implementation of Yienger and Levy 1995'
           WRITE(2,'(1a)') &
             'Tsoil, soil wetness, soil NO emission flux and canopy top flux (when) using CRF'
           ! ESS_lg_20130107+
           WRITE(2,'(2a)') &
             'Soil NO emission fluxes are calculated for Yienger and Levy emission class: ',NOemclassname(iNOemclass)
           WRITE(2,'(1a,f6.2,1a,f6.2)') &
             'Implying a wet soil NO emission factor: ',noemis_w(1),' [ngN m-2 s-1] and a dry soil NO emission factor: ',noemis_d(1)
           ! ESS_lg_20130107-
           WRITE(2,'(A10,A16,4A15)') &
             'istep','ldatltime','Tsoil','ws','NO_slflux','NO_emflux'
           WRITE(2,'(A10,A16,4A15)') &
             '#','dd-mm-yy hh:min','K','m','molec m-2 s-1','molec. m-2 s-1'
           DO i=1,Nstep
             IF (mod(i,nprint) == 0) THEN
               WRITE(2,'(I10,2x,A14,4(1x,E14.6))') &
                  i,ldatltime(i),tsoil(i),ws(i),   &
                  NO_slflux(i),NO_emflux(i)
             ENDIF
           END DO
           CLOSE(unit=2)

    ELSE
        PRINT *,'Calculation of soil-biogenic NO emissions switched off'
    END IF

    PRINT *,'Calculation of surface layer turbulent exchanges (emdep_xtsurf_calcra)'

    !--- Calculate the aerodynamic resistance:

    CALL emdep_xtsurf_calcra( Nstep,& ! ESS_lg_20120722+
         vkarman, loland,           &
	     cfml,cfncl,                &
         ril,cdnl,                  &
         geopot_3d_sl,tvir,tvl,     & ! ESS_lg_20130208+
         um1, vm1,                  &
         az0, z0m,                  &
         zustveg_obs,               & ! MAQ_20230605+ read-in ustar 
         zrahl,zrahveg,zrahslsn,    &
         zustarl,zustveg,zustslsn)

    ! mz_lg_20060719+ calling of subroutine to calculate the surface
    !   uptake resistances from the reactivity coefficient and the 
    !   henry coeff. by scaling with the SO2 and O3 uptake resistances,
    !   which are defined specifically in this subroutine in addition
    !   to a selection of other gases of the Ganzeveld et al., JGR, 1995 
    !   and 1998 papers on O3, NOx and HNO3 and SOx dry deposition     

    PRINT *,'Calculation of surface uptake resistances (emdep_xtsurf_calcrs)'

    CALL emdep_xtsurf_calc_rs(                               &
       lo_derived, lvd_bigl, trname, moleweight, reactivity, &
       henrycoeff, idt_SO2, diff, diffrb, rsoil, rwater, rws,& ! ESS_lg_20150619+
       rsnow, rmes, rcut)

	! mz_lg_20040716+ biogenic HONO and NOx (NO2) emissions from nitrate photolysis:
           ! NOTE that this subroutine is called after the calculation of the surface layer 
           ! turbulence parameters since this subroutine needs as input an estimate of the 
           ! HNO3 dry deposition flux on the canopy which uses the HNO3 dry deposition velocity 
           ! calculated as 1/(ra+rb)

	! ESS_lg_20130424+ added the calculation of the vapor pressure deficit effect
    !    calculation of the water vapour pressure of the air and within the leaf
    fvpd(:)=1. ! MAQ_lg_20170713+ here initialized at 1 to secure that a value of 1 is default used/written to output when not activated
    IF (l_fvpd) THEN
      DO i=1,nstep
 	    fvpd(i)=1.
        ! calculation of quasi-laminar boundary layer resistance
        rbveg=(2./(zustveg(i)*0.40))
        vpair=6.108*10**(7.5*(tsurf(i)-273.15)/(237.3+(tsurf(i)-273.15)))*(rh(i)*100.)/100.
        vpleaf=6.108*10**(7.5*(tsurf(i)-273.15)/(237.3+(tsurf(i)-273.15))) 
        vpsfc =(vpleaf-vpair)*rbveg/(rco_leaf(i,1)+rbveg)+vpair ! ESS_lg_20130516+ using rco_leaf of top-layer

        IF (vpsfc.gt.vpleaf) THEN
           fvpd(i)=1.
        ELSE
          fvpd(i) = 1 - 0.02*(vpleaf-vpsfc)
          if(fvpd(i).lt.0.5) fvpd(i) = 0.5
        ENDIF

		!print *,'vpsfc: ',i,fvpd(i),vpair,vpleaf,vpsfc,tsurf(i),rh(i)

        rco_leaf(i,:)=rco_leaf(i,:)/fvpd(i)
	  ENDDO
    ENDIF
	! ESS_lg_20130424-
		   		   		   
    IF (l_emis_bio_jNO3) THEN

      PRINT *,'Calculation of HONO/NOx emissions from NO3 photolysis (emdep_emis_bio_jNO3)'
 
      ! calculating the (H)NO3 dry deposition flux from the HNO3 mixing ratio, recalculated to 
      ! to concentration in molecules m-3 and multiplied with the VdHNO3 in m s-1, calculated here
      ! as 1/Ra with Ra in s m-1                
      ddNO3(:,1)=pxtm1(:,idt_HNO3)*(1e6*1.e-3_dp*prhoa(:)/(amd/avo))*1./(zrahveg(:))

      CALL emdep_xtsurf_jNO3( Nstep,                                         & ! ESS_lg_20120722+
        delta_time, l_xtsurf_veg_mlay, nveglay_hr, nveglay, zcvw, zvgrat,    &
        prc, prl, fsl, netrad, netradmax, lad, lai, ddNO3, fslsum, cthru,    & ! MAQ_lg_20170524+ added netrad + netradmax
        NO3s, HONO_jNO3em, NOx_jNO3em, HONO_emflux, NOx_emflux) 
 
      ! MAQ_lg_20170526+ all the output is now written after the call of veg_mlay since the emissions are
	  !   calculated in MLC-CHEM using the updated concentrations within that subroutine

    ELSE
      PRINT *,'Calculation of HONO/NOx emissions from NO3 photolysis switched off'
    END IF

    ! MAQ_lg_20170517+ included a subroutine to calculate the soil HONO emission flux
    IF (l_emis_soil_HONO) THEN

      !PRINT *,'Calculation of soil HONO emission flux'
      !CALL emdep_emis_soilHONO( Nstep, prhoa, ws, tsoil, pxtmveg(:,nveglay,idt_HONO), HONO_slflux) 
	ELSE
      PRINT *,'Calculation of soil HONO emissions switched off'
    END IF	
    ! MAQ_lg_20170517+
	
    ! ESS_lg_20130516+ calculation of stomatal resistance as a function of radiation, temperature, 
	!   moisture and CO2 concentrations; the AGS model
	IF (l_xtsurf_AGS) THEN
      PRINT *,'Calc. of stomatal resistance also considering role of CO2 (emdep_xtsurf_AGS)'
	  PRINT *,'The selected vegetation type (see namelist) is: ',Agstype
	  PRINT *,'See the code for the details on which vegetation type this number represents'
	  PRINT *,'ENTER to continue'
	  read (*,*)
	  
      ! MAQ_lg_20170607+ to test the implementation 
      !	  tsurf(:)=280.
      !	  netrad(:)=600.

	  CALL emdep_xtsurf_AGS(nstep, nveglay, Agstype,          & ! Ags vegetation type, C3; Agstype=1, C4; Agstype=2, conifereous forest; Agstype=3, tropical forest; Agstype=4
                          tsurf, lai, press,                  & ! lai, pressure
                          ws, wsmax, 1e6*pxtm1(:,idt_CO2),    & ! CO2 in ppmv!
                          qm1, qs, qsam,                      & ! MAQ_lg_20170720+ modified; was qm1, qm1, qm1 ! qsam, saturation point at the leaf level	 
                          netrad, rbvd, rvd, fsl,             & ! radiation parameters						  
	                      rco_leaf_AGS,  rco_leaf_AGS_ml,     & ! big leaf and canopy layer stomatal resistance
						  rmesCO2, rcutCO2, ccompCO2)	        ! MAQ_lg_20170720+ added the assignment of the CO2 compensation point in ppmv  ! CO2 exchange resistances

      ! ESS_20140416+ modified writing of radiation data for higher vertical resolution model 
	  ! show some results; canopy stomatal resistance profiles and CO2 resistances
      OPEN(unit=2,file='output/AGS.out',status='UNKNOWN')
      WRITE(2,'(1a)') &
         'AGS: big-leaf and canopy layer stomatal conductance considering the role of CO2, mesophyllic CO2 resist.'

	  ! to assign parname and parunit for the header of the file
      ALLOCATE (parname(1,nveglay))
	  ALLOCATE (parunit(1,nveglay))
 
      ! assigning some of this header info and some further manipulations for output
      DO i=1,nveglay
	    dummy='/rco_ml('
        iip=INDEX(dummy,' ')
        IF (zrefcan(1,i) >= 10.) THEN
          write(dummy2,'(f4.1)')zrefcan(1,i)
	    ELSE
          write(dummy2,'(f3.1)')zrefcan(1,i) 
	    ENDIF 
        dummy=dummy(1:iip-1)//dummy2
        iip=INDEX(dummy,' ')
 		dummy3=')'
        parname(1,i)=dummy(1:iip-1)//dummy3
        parunit(1,i)='cm s-1'
      ENDDO	 

      WRITE(2,'(A10,A16,50A15)') &
         'istep','ldatltime','/rco_AGS_bigl',(parname(1,i),i=1,nveglay),'rmesCO2','ccompCO2' ! MAQ_lg_20170720+ added the CO2 compensation point 
      WRITE(2,'(A10,A16,50A15)') &
         '#','dd-mm-yy hh:min','[cm s-1]',(parunit(1,i),i=1,nveglay),'[s m-1]','[ppmv]'      ! MAQ_lg_20170720+ added the CO2 compensation point 
      DO i=1,Nstep
        IF (mod(i,nprint) == 0) THEN
          WRITE(2,'(I10,2x,A14,50(1x,E14.6))') &
              i,ldatltime(i),100./rco_leaf_AGS(i),(100./rco_leaf_AGS_ml(i,ii),ii=1,nveglay),rmesCO2(i),ccompCO2(i)  ! MAQ_lg_20170720+ added the CO2 compensation point 
        ENDIF
      END DO

	  DEALLOCATE (parname)
      DEALLOCATE (parunit)
      ! ESS_20140416-

      CLOSE(unit=2)

	  ! assigning the AGS multi-layer stomatal resistance to rco_leaf used in xtsurf_veg_mlay
  	  DO i=1,Nstep
        rco_leaf(i,:)=rco_leaf_AGS_ml(i,:)

        ! MAQ_lg_20190905+ deactivated since the VPD is already considered in Ags
		! ESS_lg_20160303+ added the calculation of the vapor pressure deficit effect
        !    calculation of the water vapour pressure of the air and within the leaf
        !IF (l_fvpd) THEN
        !  rco_leaf(i,:)=rco_leaf_AGS_ml(i,:)/fvpd(i)
		!ENDIF
        ! ESS_lg_20160303-
        ! MAQ_lg_20190905-
		
	  ENDDO
    ENDIF
	! ESS_lg_20130424-

	! ESS_lg_20130516-

    ! mz_lg_20050721+ and calculation of atmosphere-biosphere exchanges 

    IF (l_xtsurf_veg_mlay) THEN

       PRINT *,'Preparing calculation of atmosphere-biosphere exchanges'

       ! mz_LG_20020115 calling of subroutine in which the vertical 
       !     windspeed profile in the canopy is calculated. These      
       !     profiles are being used in the multi-layer model to
       !     calculate the turbulent exchange between the layers.

       PRINT *,'Calculation of canopy wind speed profile (emdep_xtsurf_wndprof)'

       CALL emdep_xtsurf_wndprof( Nstep, & ! ESS_lg_20120722+
         nveglay, zustveg, lai, hc, forestfr, z0m, disp, u_veg)

       ! mz_lg_20050721+ calculation of photolysis profiles in canopy
       !    needed for the calculation of within-canopy photochemistry 

       IF (l_xtsurf_veg_mlay_chem) THEN

         rj_veg(:,:,:)=0._dp
         IF (l_xtsurf_veg_mlay_chem_photolysis) THEN
           PRINT *,'Calculation of canopy photolysis profiles (emdep_xtsurf_rjveg)'

           CALL emdep_xtsurf_rjveg( Nstep,           & ! ESS_lg_20120722+
                nveglay, nveglay_hr, latmbios_photo, &
                lai, hc, rvdsl, rvd, fsl, rj, rj_veg)  
         ENDIF

       ENDIF

       ! ESS_lg_20120718+ call of subroutine to assign/calculate tracer compensation point
       ! ESS_lg_20150807+ important change of the location of the initialization of ccomp,
	   !   Before it was done only when l_xtsurf_veg_mlay_ccomp was true but this resulted
	   !   in not defined values for ccomp when the switch was false
	   
       ccomp(:)        = 0._dp ! resetting the tracer compensation points
       IF (l_xtsurf_veg_mlay_ccomp) THEN

         PRINT *,'Calculation/assignment of tracer compensation points (emdep_xtsurf_ccomp)'         
         CALL emdep_xtsurf_ccomp(Nstep,prhoa,vgrat,tsurf,pxtmveg(:,1,idt_CO2),idt_CO2,&
		                         idt_NO2,idt_NH3,ccomp)

       ENDIF

       ! mz_lg_20050719+ the actual call to the subroutine that calculates
       ! the atmosphere-biosphere exchanges

       PRINT *,'Calculation of atmosphere-biosphere exchanges (emdep_xtsurf_veg_mlay)'
       IF (l_xtsurf_veg_mlay_chem) PRINT *,'Including the simulations of canopy chemistry'

       CALL emdep_xtsurf_veg_mlay( casename, nstep,                    & ! MAQ_lg_20210330+ added casename ! ESS_lg_20120722+
            latmbios, latmbios_emveg, latmbios_emsoil,                 & !
            lvd_bigl, lexist_GAS, trname, ccomp, delta_time,           & !
            l_emis, l_emis_bio_NO, l_emis_bio_VOC, l_emis_bio_jNO3,    & !
            l_emis_soil_HONO, l_drydep, l_xtsurf_Ags, loland,          & ! MAQ_20201105+ l_xtsurf_Ags, MAQ_lg_20170517+ added l_emis_soil_HONO
			zfrl, zcvbs, zcvw, zvgrat, zrefsl,                         & ! ESS_lg_20131125+ zref, reference height surface layer
			zrahveg, zustveg, u_sl, u_veg,                             & ! ESS_lg_20131125+ u_sl 
			netrad, rj, rj_veg, tsurf, tsoil, ws, press, qm1, rh,      & ! MAQ_lg_20170517+ tsoil and ws, ESS_lg_20120725+
            lai, lad, hc, rco_leaf, fws, diff, diffrb, rcut, rmes,     & !
            rws, rsoil, rmesCO2, rcutCO2, ccompCO2, grvol, grmass,     & ! MAQ_lg_20170720+ added the CO2 compensation point ! ESS_lg_20130516+ added rmesCO2 and rcutCO2
            pdp, prhoa, prc, prl, fsl, fslsum, netradmax, cthru,       & ! MAQ_lg_20170526+ extra input to constrain HONO/NOx emission calculations
            NO3s, HONO_jNO3em, NOx_jNO3em,                             & !
			NO_slflux, NO2_slflux, ISOP_emflux, MONO_emflux, OVOC_emflux, & ! MAQ_lg_20210727+ added NO2 soil (deposition) flux for extra diagnostics
            HONO_emflux, HONO_slflux, ddNO3, NOx_emflux,               & ! MAQ_lg_20170517+ added HONO soil flux and ddNO3 
            NH3_emflux, NH3_slflux,                                    & ! MAQ_lg_20230329+ added NH3 vegetation/soil emissions			
			RADON_slflux, CO2_slflux, SQT_slflux,                      & ! MAQ_lg_20230601+ added SQT soil flux, ESS_lg_20130503+ added CO2
			fAPIN, fBPIN, fSQTERP, fATERP, fLIMO, fMYRC,               & ! MAQ_20170413+ added extra MT's ! ESS_lg_20130119+ added the partitioning of monoterpene and sqterp emissions
            pxtm1, pxtmveg, pxtm1_obs,                                 & ! ESS_lg_20120722+ added pxtm1_obs
            weight_pxtm1_obs, Kh, Kh_obs_sl, Kh_obs_cl, zustveg_obs,   & ! ESS_lg_20120903+ weight_pxtm1_obs ! &  ESS_lg_20120722+ added Kh for diagnostics ! MAQ_AV_20200309+ added Kh_obs_sl & Kh_obs_cl, MAQ_AV_20201005+ added zustveg_obs
            patmbiosflux, pcrf, pvdveg, psurfflux,                     & !
            pstomfluxO3, xteemis, xtedryd, xtechem, xtediff, OHreact)    ! ESS_lg_20140509+ pstomflux

       PRINT *,'End calculation of atmosphere-biosphere exchanges (emdep_xtsurf_veg_mlay)'

       ! show some results; canopy concentrations, fluxes, etc.
       OPEN(unit=2,file='output/veg_mlay.out',status='UNKNOWN')
       WRITE(2,'(1a)') &
         'Canopy exchange model output including tracer concentrations, fluxes and CRFs'
       WRITE(2,'(1a)') &
         'CRF reflects here the explicitly calculated ratio of canopy-top to emission flux'

	   ! ESS_20140416+ modified writing of radiation data for higher vertical resolution model 

	   ! to assign parname and parunit for the header of the file
	   ALLOCATE (parname(npar_out,nveglay*ntrac))
	   ALLOCATE (parunit(npar_out,nveglay*ntrac))
	   ALLOCATE (mixratios(Nstep,nveglay+1,ntrac))
 
       ! assigning some of this header info and some further manipulations for output

	   ! ---wind speed
	   npar=1
	   iwind=npar
	   nlevels(npar)=0
 	   nlevels(npar)=nlevels(npar)+1
	   dummy='u('
       iip=INDEX(dummy,' ')
       IF (zref(1) < 10.) THEN
		 write(dummy2,'(f3.1)')zref(1)
	   ELSE
		 write(dummy2,'(f4.1)')zref(1)
	   ENDIF
       dummy=dummy(1:iip-1)//dummy2
       iip=INDEX(dummy,' ')
 	   dummy3='m)'
       parname(npar,nlevels(npar))=dummy(1:iip-1)//dummy3
	   parunit(npar,nlevels(npar))='m s-1'
	   ! Canopy layers
	   DO jk=1,nveglay
         nlevels(npar)=nlevels(npar)+1
		 dummy='u('
         iip=INDEX(dummy,' ')
         IF (zrefcan(1,jk) < 10.) THEN
		   write(dummy2,'(f3.1)')zrefcan(1,jk)
		 ELSE
		   write(dummy2,'(f4.1)')zrefcan(1,jk)
		 ENDIF
         dummy=dummy(1:iip-1)//dummy2
         iip=INDEX(dummy,' ')
 		 dummy3='m)'
         parname(npar,nlevels(npar))=dummy(1:iip-1)//dummy3
		 parunit(npar,nlevels(npar))='m s-1'
       ENDDO

	   ! ---Kh, between surface layer and top-canopy layer and between canopy layers
	   npar=npar+1
	   iKh=npar
	   nlevels(npar)=0
       nlevels(npar)=nlevels(npar)+1
	   dummy='Kh(~'
       iip=INDEX(dummy,' ')
       ! Kh(1) is representative from zref to reference height top-canopy layer
       IF (zref(1)+zrefcan(1,1) < 10.) THEN
		 write(dummy2,'(f3.1)') (zref(1)+zrefcan(1,1))/2.
	   ELSE
		 write(dummy2,'(f4.1)')(zref(1)+zrefcan(1,1))/2.
	   ENDIF
       dummy=dummy(1:iip-1)//dummy2
       iip=INDEX(dummy,' ')
 	   dummy3='m)'
       parname(npar,nlevels(npar))=dummy(1:iip-1)//dummy3
	   parunit(npar,nlevels(npar))='m s-1'
	   ! Canopy layers
	   DO jk=1,nveglay-1
         nlevels(npar)=nlevels(npar)+1
		 dummy='Kh(~'
         iip=INDEX(dummy,' ')
         IF ((zrefcan(1,jk)+zrefcan(1,jk+1))/2 < 10.) THEN
		   write(dummy2,'(f3.1)') (zrefcan(1,jk)+zrefcan(1,jk+1))/2. 
		 ELSE
		   write(dummy2,'(f4.1)')(zrefcan(1,jk)+zrefcan(1,jk+1))/2. 
		 ENDIF
         dummy=dummy(1:iip-1)//dummy2
         iip=INDEX(dummy,' ')
 		 dummy3='m)'
         parname(npar,nlevels(npar))=dummy(1:iip-1)//dummy3
		 parunit(npar,nlevels(npar))='m s-1'
       ENDDO	

       ! ---jNO2, surface and canopy layers
	   npar=npar+1
	   ijNO2=npar
	   nlevels(npar)=0
 	   nlevels(npar)=nlevels(npar)+1
	   dummy='jNO2('
       iip=INDEX(dummy,' ')
       IF (zref(1) < 10.) THEN
		 write(dummy2,'(f3.1)')zref(1)
	   ELSE
		 write(dummy2,'(f4.1)')zref(1)
	   ENDIF
       dummy=dummy(1:iip-1)//dummy2
       iip=INDEX(dummy,' ')
 	   dummy3='m)'
       parname(npar,nlevels(npar))=dummy(1:iip-1)//dummy3
	   parunit(npar,nlevels(npar))='s-1'
	   ! Canopy layers
	   DO jk=1,nveglay
         nlevels(npar)=nlevels(npar)+1
		 dummy='jNO2('
         iip=INDEX(dummy,' ')
         IF (zrefcan(1,jk) < 10.) THEN
		   write(dummy2,'(f3.1)')zrefcan(1,jk)
		 ELSE
		   write(dummy2,'(f4.1)')zrefcan(1,jk)
		 ENDIF
         dummy=dummy(1:iip-1)//dummy2
         iip=INDEX(dummy,' ')
 		 dummy3='m)'
         parname(npar,nlevels(npar))=dummy(1:iip-1)//dummy3
		 parunit(npar,nlevels(npar))='s-1'
       ENDDO		
		 
       ! ---VdO3, only canopy layers
	   npar=npar+1
	   iVdO3=npar
	   nlevels(npar)=0
	   ! Canopy layers
	   DO jk=1,nveglay
         nlevels(npar)=nlevels(npar)+1
		 dummy='VdO3('
         iip=INDEX(dummy,' ')
         IF (zrefcan(1,jk) < 10.) THEN
		   write(dummy2,'(f3.1)')zrefcan(1,jk)
		 ELSE
		   write(dummy2,'(f4.1)')zrefcan(1,jk)
		 ENDIF
         dummy=dummy(1:iip-1)//dummy2
         iip=INDEX(dummy,' ')
 		 dummy3='m)'
         parname(npar,nlevels(npar))=dummy(1:iip-1)//dummy3
		 parunit(npar,nlevels(npar))='cm s-1'
       ENDDO		

       ! ---OH reactivity, surface and canopy layers
	   npar=npar+1
	   iOHr=npar
	   nlevels(npar)=0
 	   nlevels(npar)=nlevels(npar)+1
	   dummy='OHr('
       iip=INDEX(dummy,' ')
       IF (zref(1) < 10.) THEN
		 write(dummy2,'(f3.1)')zref(1)
	   ELSE
		 write(dummy2,'(f4.1)')zref(1)
	   ENDIF
       dummy=dummy(1:iip-1)//dummy2
       iip=INDEX(dummy,' ')
 	   dummy3='m)'
       parname(npar,nlevels(npar))=dummy(1:iip-1)//dummy3
	   parunit(npar,nlevels(npar))='s-1'
	   ! Canopy layers
	   DO jk=1,nveglay
         nlevels(npar)=nlevels(npar)+1
		 dummy='OHr('
         iip=INDEX(dummy,' ')
         IF (zrefcan(1,jk) < 10.) THEN
		   write(dummy2,'(f3.1)')zrefcan(1,jk)
		 ELSE
		   write(dummy2,'(f4.1)')zrefcan(1,jk)
		 ENDIF
         dummy=dummy(1:iip-1)//dummy2
         iip=INDEX(dummy,' ')
 		 dummy3='m)'
         parname(npar,nlevels(npar))=dummy(1:iip-1)//dummy3
		 parunit(npar,nlevels(npar))='s-1'
       ENDDO			   
	   	   
       ! ---O3 stomatal flux, only canopy layers
	   npar=npar+1
	   istomfluxO3=npar
	   nlevels(npar)=0
	   ! Canopy layers
	   DO jk=1,nveglay
         nlevels(npar)=nlevels(npar)+1
		 dummy='O3stfl('
         iip=INDEX(dummy,' ')
         IF (zrefcan(1,jk) < 10.) THEN
		   write(dummy2,'(f3.1)')zrefcan(1,jk)
		 ELSE
		   write(dummy2,'(f4.1)')zrefcan(1,jk)
		 ENDIF
         dummy=dummy(1:iip-1)//dummy2
         iip=INDEX(dummy,' ')
 		 dummy3='m)'
         parname(npar,nlevels(npar))=dummy(1:iip-1)//dummy3
		 parunit(npar,nlevels(npar))='molec. m-2 s-1'
       ENDDO
	   
	   ! ESS_lg_20140418+ tracer concentrations (only a selection, see tracer_init)
       npar=npar+1
	   iconc=npar
	   nlevels(npar)=0
	   jjt=0
	   recalcmrconc=1.e-3_dp*prhoa(1)/(amd/avo)
	   DO jt=1,ntrac
 	  	 IF (latmbios_output(jt)) THEN
           ! Surface layer concentrations
 		   nlevels(npar)=nlevels(npar)+1
		   jjt=jjt+1
	       dummy=trname(jt)
           iip=INDEX(dummy,' ')
 		    dummy2='('
           dummy=dummy(1:iip-1)//dummy2
           iip=INDEX(dummy,' ')
           IF (zref(1) < 10.) THEN
			 write(dummy2,'(f3.1)')zref(1)
		   ELSE
			 write(dummy2,'(f4.1)')zref(1)
		   ENDIF
           dummy=dummy(1:iip-1)//dummy2
           iip=INDEX(dummy,' ')
 		   dummy3='m)'
           parname(npar,nlevels(npar))=dummy(1:iip-1)//dummy3
 
           ! recalculation of mixing ratio to pptv/ppbv/ppmv
           recalcconc=1e15 ! starting at 1e-3 pptv level
           iunit=1
		   DO i=1,Nstep
             IF (recalcconc*pxtm1(i,jt) > 1.e3 .AND. iunit < nunits) THEN ! ESS_lg_20150611: iunit should not exceed nunits
			   recalcconc=recalcconc*1e-3
			   iunit=iunit+1
			 ENDIF
		   ENDDO
		   parunit(npar,nlevels(npar))=mixratio_unit(iunit)

		   ! Canopy layers
		   DO jk=1,nveglay
             nlevels(npar)=nlevels(npar)+1
			 dummy=trname(jt)
             iip=INDEX(dummy,' ')
 		     dummy2='('
             dummy=dummy(1:iip-1)//dummy2
             iip=INDEX(dummy,' ')
             IF (zrefcan(1,jk) < 10.) THEN
			   write(dummy2,'(f3.1)')zrefcan(1,jk)
			 ELSE
			   write(dummy2,'(f4.1)')zrefcan(1,jk)
			 ENDIF
             dummy=dummy(1:iip-1)//dummy2
             iip=INDEX(dummy,' ')
 		     dummy3='m)'
             parname(npar,nlevels(npar))=dummy(1:iip-1)//dummy3

             ! recalculation of mixing ratio to pptv/ppbv/ppmv
	  		 DO i=1,Nstep
               IF (recalcconc*pxtmveg(i,jk,jt) > 1.e3 .AND. iunit < nunits) THEN ! ESS_lg_20150611: iunit should not exceed nunits
			     recalcconc=recalcconc*1e-3
				 iunit=iunit+1
			   ENDIF
			 ENDDO
 			 parunit(npar,nlevels(npar))=mixratio_unit(iunit)

		   ENDDO
           ! MAQ_lg_20160608+ whenever the concentrations inside the canopy are higher then those in the 
		   ! surface layer (e.g. ppbv vs pptv level) then, because we use the final recalculation term for the
		   ! the lowest vegetation, then all the concentration units are corrected
		   parunit(npar,nlevels(npar)-nveglay:nlevels(npar))=mixratio_unit(iunit)

		   ! special cases (OH, radon, etc)
		   IF (jt == idt_OH) THEN
			 recalcconc=1e-6*recalcmrconc
			 parunit(npar,nlevels(npar)-nveglay:nlevels(npar))='1e6 molec cm-3'
		   ENDIF
		   IF (jt == idt_rad) THEN
			 recalcconc=recalcmrconc
			 parunit(npar,nlevels(npar)-nveglay:nlevels(npar))='atoms cm-3'
		   ENDIF

		   mixratios(:,1,jjt)=recalcconc*pxtm1(:,jt)
		   mixratios(:,2:nveglay+1,jjt)=recalcconc*pxtmveg(:,:,jt)

		 ENDIF
	   ENDDO		 
		 
       ! ESS_lg_20140418+ and writing the file 
	   IF (l_xtsurf_veg_mlay_chem) THEN
         WRITE(2,'(A10,A16,350A15)') &
               'istep','ldatltime'                                                        &
			  ,parname(iwind,1:nlevels(iwind))                                            & 
			  ,'u*-veg','netrad','Tair','Tsurf','RH','RiB','ML-depth'                     & ! ESS_lg_20120928+ added mixed layer depth
              ,parname(iKh,1:nlevels(iKh))                                                & 
              ,'LAI','stom. cond-cl.','fvpd','ws','f-ws','scale_rad','f-wetskin'          & ! MAQ_20170606+ fvpd ! MAQ_2060621+ added LAI ! ESS_lg_20150406+ f-wetskin ! ESS_lg_20150324+ add scaling_rad ! ESS_lg_20130127+ added soil moisture
			  ,parname(ijNO2,1:nlevels(ijNO2))                                            & 
			  ,parname(iconc,1:nlevels(iconc))                                            &
              ,parname(iVdO3,1:nlevels(iVdO3))                                            &
              ,'VdNO2_crownl','VdISOP_soill','VdCH2O_crownl'                              & ! ESS_lg_20150406+ VdCH2O
              ,'VdHONO_cl','VdHONO_sll'                                                   & ! MAQ_lg_20170524+ VdHONO
              ,'VdNH3_cl','VdNH3_sll'                                                     & ! MAQ_lg_20230329+ VdNH3
              ,'VdNO_cl','VdNO_sll'                                                       & ! MAQ_lg_20180730+ VdNO
              ,'VdCO2_cl','VdCO2_sll'                                                     & ! ESS_lg_20130811+ VdCO2
              ,'VdCOS_cl','VdCOS_sll'                                                     & ! MAQ_lg_20190111+ VdCOS
              ,'VdAPIN_cl','VdAPIN_sll'                                                   &

!              ,'VdAPINP1a_cl','VdAPINP1a_sll'                                             &

              ,'NO_slflux','NO2_slflux','emflx-NOx_cl','emflx-NOx_sll'                    & ! MAQ_lg_20210727+ added NO2 soil (deposition) flux for extra diagnostics
              ,'emflx-HONO_cl','emflx-HONO_sll','emflx-NH3_cl','emflx-NH3_sll'            & ! MAQ_lg_20230329+ added NH3 emission flux
              ,'emflx-ISOP_cl','emflx-ISOP_sll','emflx-APIN_cl','emflx-APIN_sll'          &
              ,'emflx-BPIN_cl','emflx-BPIN_sll','emflx-SQTP_cl','emflx-SQTP_sll'          &
			  ,'Rad_slflux','CO2_slflux'                                                  & ! MAQ_lg_20170830+ CO2 soilflux ! ESS_lg_20150601+ Radon emission flux
              ,'atmbioflx-O3','atmbioflx-NO','atmbioflx-NO2','atmbioflx-HONO'             & ! ESS_lg_20150823+ added HONO
              ,'atmbioflx-NH3','atmbioflx-ISOP','atmbioflx-APIN','atmbioflx-SQT'          & ! MAQ_lg_20230329+ added NH3
              ,'atmbioflx-CO2','atmbioflx-COS','atmbioflx-Rad'                            & ! MAQ_lg_20190917+ COS ! ESS_lg_20150611+ Radon ! ESS_lg_20130619+
              ,'CRF_NOx','CRF_HONO','CRF_ISOP','CRF_APIN','CRF_SQTERP'                    & ! ESS_lg_20150823+ added CRF HONO
              ,'CRF_aTERP','CRF_LIMO','CRF_MYRC'                                          & ! MAQ_lg_20170413+
              ,parname(iOHr,1:nlevels(iOHr))                                              & ! ESS_lg_20130113+
              ,parname(istomfluxO3,1:nlevels(istomfluxO3))                                & ! ESS_lg_20140508+
			  ,'VdO3-canopy'                                                                ! ESS_lg_20150624+ Canopy-scale VdO3
		 WRITE(2,'(A10,A16,350A15)') &
               '#','dd:mm:yy hh:min'                                                      &
			  ,parunit(iwind,1:nlevels(iwind))                                            &
			  ,'m s-1','W m-2','K','K','0-1','-','m'                                      & ! ESS_lg_20120928+ 
			  ,parunit(iKh,1:nlevels(iKh))                                                &
              ,'m2 m-2','cm s-1','0-1','m','0-1','0-1','0-1'                              & ! MAQ_lg_20170606+ fvpd ! MAQ_20160621+ added LAI ! ESS_lg_20150406+ f-wetskin ! ESS_lg_20150324+ added scaling_rad ! ESS_lg_20130127+ added soil moisture
			  ,parunit(ijNO2,1:nlevels(ijNO2))                                            &
			  ,parunit(iconc,1:nlevels(iconc))                                            & 
			  ,parunit(iVdO3,1:nlevels(iVdO3))                                            &
              ,'cm s-1','cm s-1','cm s-1'                                                 & ! ESS_lg_20150406+ VdCH2O
              ,'cm s-1','cm s-1'                                                          & ! MAQ_lg_20170524+ VdHONO
			  ,'cm s-1','cm s-1'                                                          & ! MAQ_lg_20230329+ VdNH3 
			  ,'cm s-1','cm s-1'                                                          & ! MAQ_lg_20180730+ VdNO
              ,'cm s-1','cm s-1'                                                          & ! ESS_lg_20130811+ VdCO2
              ,'cm s-1','cm s-1'                                                          & ! MAQ_lg_20190111+ VdCOS
              ,'cm s-1','cm s-1'                                                          &
              ,'molec m-2 s-1','molec m-2 s-1','molec m-2 s-1','molec m-2 s-1'            & ! MAQ_lg_20210727+ added NO2 soil (deposition) flux for extra diagnostics
              ,'molec m-2 s-1','molec m-2 s-1','molec m-2 s-1','molec m-2 s-1'            & ! MAQ_lg_20230329+ added NH3 emission fluxes
              ,'molec m-2 s-1','molec m-2 s-1','molec m-2 s-1','molec m-2 s-1'            &
              ,'molec m-2 s-1','molec m-2 s-1','molec m-2 s-1','molec m-2 s-1'            &
              ,'atoms m-2 s-1','molec m-2 s-1'                                            & ! MAQ_lg_20170830+ CO2 soilflux! ESS_lg_20150611+ Radon
              ,'molec m-2 s-1','molec m-2 s-1','molec m-2 s-1','molec m-2 s-1'            & ! ESS_lg_20150823+ added HONO
              ,'molec m-2 s-1','molec m-2 s-1','molec m-2 s-1','molec m-2 s-1'            & ! MQ_lg_20230329+ NH3
              ,'molec m-2 s-1','molec m-2 s-1','atoms m-2 s-1'                            & ! MAQ_lg_20190917+ COS ! ESS_lg_20150611+ Radon ! ESS_lg_20130619+
              ,'ratio','ratio','ratio','ratio','ratio'                                    & ! ESS_lg_20150823+ added CRF HONO
              ,'ratio','ratio','ratio'                                                    & ! MAQ_lg_20170413+
              ,parunit(iOHr,1:nlevels(iOHr))                                              & ! ESS_lg_20130113+
              ,parunit(istomfluxO3,1:nlevels(istomfluxO3))                                & ! ESS_lg_20140508+
			  ,'cm s-1'                                                                     ! ESS_lg_20150624+ Canopy-scale VdO3
		   DO i=1,Nstep
           IF(mod(i,nprint) == 0) THEN
             recalcmrconc=1.e-3_dp*prhoa(i)/(amd/avo)
             VdO3_canopy=-1e-4*patmbiosflux(i,idt_O3)/(recalcmrconc*pxtm1(i,idt_O3))        ! MAQ_lg_20220127+ calculating the VdO3_canopy for writing properly to output file
             WRITE(2,'(I10,2x,A14,350(1x,E14.6))')                                         &
                 i,ldatltime(i)                                                            &
				,sqrt(um1(i)**2+vm1(i)**2),u_veg(i,1:nveglay)                              &
				,zustveg(i),netrad(i),Tair(i),Tsurf(i),rh(i),ril(i),MLH(i)                 &
				,Kh(i,1:nveglay)                                                           & 
                ,lai(i),100./rco_leaf(i,1),fvpd(i),ws(i),fws(i),scaling_rad(i),zcvw(i)     & ! MAQ_lg_20170606+ fvpd ! MAQ_lg_20160621+ added LAI ! ESS_lg_20150324+ added scaling_rad ! ESS_lg_20130127+ added soil moisture
				,rj(i,idt_NO2),rj_veg(i,1:nveglay,idt_NO2)                                 & 
                ,mixratios(i,1:nveglay+1,1:jjt)                                            &
                ,pvdveg(i,1:nveglay,idt_O3)                                                &
				,pvdveg(i,1,idt_NO2),pvdveg(i,nveglay,idt_ISOP),pvdveg(i,1,idt_CH2O)       & ! ESS_lg_20150406+ VdCH2O (top layer NO2/CH2O and soil layer isoprene)
				,pvdveg(i,1,idt_HONO),pvdveg(i,nveglay,idt_HONO)                           & ! MAQ_lg_20170524+ VdHONO (top layer and soil layer)
				,pvdveg(i,1,idt_NH3),pvdveg(i,nveglay,idt_NH3)                             & ! MAQ_lg_20230329+ VdNH3 (top layer and soil layer)
				,pvdveg(i,1,idt_NO),pvdveg(i,nveglay,idt_NO)                               & ! MAQ_lg_20180730+ VdNO (top layer and soil layer)
				,pvdveg(i,1,idt_CO2),pvdveg(i,nveglay,idt_CO2)                             & ! ESS_lg_20130811+ VdCO2 (top layer and soil layer)
				,pvdveg(i,1,idt_COS),pvdveg(i,nveglay,idt_COS)                             & ! MAQ_lg_20190111+ VdCOS (top layer and soil layer)
                ,pvdveg(i,1,idt_apin),pvdveg(i,nveglay,idt_apin)                           &    

!                ,pvdveg(i,1,idt_APINP1a),pvdveg(i,nveglay,idt_APINP1a)                     &

                ,NO_slflux(i),NO2_slflux(i),NOx_emflux(i,1),NOx_emflux(i,nveglay)          & ! MAQ_lg_20210727+ added NO2 soil (deposition) flux for extra diagnostics
                ,HONO_emflux(i,1),HONO_emflux(i,nveglay)                                   &
                ,NH3_emflux(i,1),NH3_emflux(i,nveglay)                                     & ! MAQ_lg_20230329+ NH3 emission flux 
                ,ISOP_emflux(i,1),ISOP_emflux(i,nveglay)                                   &
                ,fAPIN*MONO_emflux(i,1),fAPIN*MONO_emflux(i,nveglay)                       & 
                ,fBPIN*MONO_emflux(i,1),fBPIN*MONO_emflux(i,nveglay)                       & 
                ,fSQTERP*MONO_emflux(i,1),fSQTERP*MONO_emflux(i,nveglay)                   &
                ,radon_slflux(i), CO2_slflux(i)                                            & ! MAQ_lg_20170830+ CO2 soilflux ! ESS_lg_20150611+ radon soil flux				
                ,patmbiosflux(i,idt_O3),patmbiosflux(i,idt_NO),patmbiosflux(i,idt_NO2)     &
                ,patmbiosflux(i,idt_HONO),patmbiosflux(i,idt_NH3)                          & ! MAQ_lg_20230329+ NH3 ! ESS_lg_20150823+ added HONO
                ,patmbiosflux(i,idt_ISOP),patmbiosflux(i,idt_APIN),patmbiosflux(i,idt_SQTERP) &                      
                ,patmbiosflux(i,idt_CO2),patmbiosflux(i,idt_COS),patmbiosflux(i,idt_RAD)   & ! MAQ_lg_20190917+ COS ! ESS_lg_20150611+ Radon flux ! ESS_lg_20130619+ added CO2
                ,pcrf(i,idt_NO2),pcrf(i,idt_HONO)                                          & ! ESS_lg_20150823+ added CRF HONO
				,pcrf(i,idt_ISOP),pcrf(i,idt_APIN),pcrf(i,idt_SQTERP)                      &
                ,pcrf(i,idt_ATERP),pcrf(i,idt_LIMO),pcrf(i,idt_MYRC)                       & ! MAQ_lg_20170413+ added MT's 
                ,OHreact(i,1:nveglay+1),pstomfluxO3(i,1:nveglay)                           & ! ESS_lg_20140508+ pstomfluxO3 ESS_lg_20130113+
                ,VdO3_canopy                                                                 ! MAQ_lg_20220127+ now writing VdO3_canopy to secure also writing values < 0 !ESS_lg_20150624+ Canopy-scale VdO3
			 ENDIF
		 END DO

       ELSE ! non-chemistry output

         WRITE(2,'(A10,A16,100A15)') &
              'istep','ldatltime','u-sl','u_crownl','u_soill','u*-veg'               &
             ,'netrad','Tair','Tsurf','RH','RiB','ML-depth','Kh_sl','Kh_soill'       & ! ESS_lg_20120928+ added mixed layer depth
             ,'stom. cond-cl.','ws','f-ws','jNO2_sl','jNO2_crownl','jNO2_soill'      & ! ESS_lg_20130127+ added soil moisture
             ,'O3_sl','O3_crownl','O3_soill','NO_sl','NO_crownl','NO_soill'          &
             ,'NO2_sl','NO2_crownl','NO2_soill','ISOP_sl','ISOP_crownl','ISOP_soill' &
             ,'APIN_sl','APIN_crownl','APIN_soill','Rn_sl','Rn_crownl','Rn_soill'    &
			 ,'CO2_sl','CO2_crownl','CO2_soill'                                      & ! ESS_lg_20130503+ added CO2  
             ,'VdO3_crownl','VdO3_soill','VdNO2_crownl','VdISOP_soill'               &
             ,'VdCO2_cl','VdCO2_sll'                                                 & ! ESS_lg_201308511+ VdCO2
             ,'VdAPINP1a_cl','VdAPINP1a_sll'                                         &
             ,'atmbioflx-O3','atmbioflx-NO','atmbioflx-NO2'                          &
             ,'atmbioflx-ISOP','atmbioflx-APIN','atmbioflx-SQT'                      &
             ,'CRF_NOx','CRF_HONO','CRF_ISOP','CRF_APIN','CRF_SQTERP'                  ! ESS_lg_20150823+ added CRF HONO
         WRITE(2,'(A10,A16,100A15)') &
              '#','dd:mm:yy hr:min','m s-1','m s-1','m s-1','m s-1'                  &
             ,'W m-2','K','K','0-1','-','m','m2 s-1','m2 s-1'                        &
             ,'cm s-1','m','0-1','s-1','s-1','s-1'                                   & ! ESS_lg_20130127+ added soil moisture
             ,'ppbv','ppbv','ppbv','pptv','pptv','pptv'                              &
             ,'pptv','pptv','pptv','ppbv','ppbv','ppbv'                              &
             ,'pptv','pptv','pptv','atoms cm-3','atoms cm-3','atoms cm-3'            &
	         ,'ppmv','ppmv','ppmv'                                                   & ! ESS_lg_20130503+ CO2
             ,'cm s-1','cm s-1','cm s-1','cm s-1'                                    &
			 ,'cm s-1','cm s-1'                                                      & ! ESS_lg_201308511+ VdCO2
			 ,'cm s-1','cm s-1'                                                      &
             ,'molec m-2 s-1','molec m-2 s-1','molec m-2 s-1'                        &
             ,'molec m-2 s-1','molec m-2 s-1','molec m-2 s-1'                        &
             ,'ratio','ratio','ratio','ratio','ratio'                                  ! ESS_lg_20150823+ added CRF HONO
         DO i=1,Nstep
           IF(mod(i,nprint) == 0) THEN
             WRITE(2,'(I10,2x,A14,100(1x,E14.6))') &
               i,ldatltime(i),sqrt(um1(i)**2+vm1(i)**2),u_veg(i,1),u_veg(i,2),zustveg(i) &
              ,netrad(i),Tair(i),Tsurf(i),rh(i),ril(i),MLH(i),Kh(i,1),Kh(i,2)            & ! ESS_lg_20130127+ added soil moisture
              ,100./rco_leaf(i,1),ws(i),fws(i),rj(i,idt_NO2),rj_veg(i,1,idt_NO2),rj_veg(i,2,idt_NO2) &
              ,1e9*pxtm1(i,idt_O3),1e9*pxtmveg(i,1,idt_O3),1e9*pxtmveg(i,2,idt_O3)       &
              ,1e12*pxtm1(i,idt_NO),1e12*pxtmveg(i,1,idt_NO),1e12*pxtmveg(i,2,idt_NO)    &
              ,1e12*pxtm1(i,idt_NO2),1e12*pxtmveg(i,1,idt_NO2),1e12*pxtmveg(i,2,idt_NO2) &
              ,1e9*pxtm1(i,idt_ISOP),1e9*pxtmveg(i,1,idt_ISOP),1e9*pxtmveg(i,2,idt_ISOP) &
              ,1e12*pxtm1(i,idt_APIN),1e12*pxtmveg(i,1,idt_APIN),1e12*pxtmveg(i,2,idt_APIN) &
              ,recalcmrconc*pxtm1(i,idt_RAD),recalcmrconc*pxtmveg(i,1,idt_RAD),recalcmrconc*pxtmveg(i,2,idt_RAD) &
              ,1e6*pxtm1(i,idt_CO2),1e6*pxtmveg(i,1,idt_CO2),1e6*pxtmveg(i,2,idt_CO2)    & ! ESS_lg_20130503+ CO2 
              ,pvdveg(i,1,idt_O3),pvdveg(i,2,idt_O3),pvdveg(i,1,idt_NO2)                 &
              ,pvdveg(i,2,idt_ISOP),pvdveg(i,1,idt_CO2),pvdveg(i,2,idt_CO2)              & ! ESS_lg_201308511+ VdCO2
              ,pvdveg(i,1,idt_APINP1a),pvdveg(i,2,idt_APINP1A)                           & 
              ,patmbiosflux(i,idt_O3),patmbiosflux(i,idt_NO),patmbiosflux(i,idt_NO2)     &
              ,patmbiosflux(i,idt_ISOP),patmbiosflux(i,idt_APIN),patmbiosflux(i,idt_SQTERP)  &
              ,pcrf(i,idt_NO2),pcrf(i,idt_HONO),pcrf(i,idt_ISOP),pcrf(i,idt_APIN),pcrf(i,idt_SQTERP) ! ESS_lg_20150823+ added CRF HONO
           ENDIF 
         ENDDO
       ENDIF

	   DEALLOCATE (parname)
  	   DEALLOCATE (parunit)
	   DEALLOCATE (mixratios)
		 
	   ! ESS_20140416- 

	   CLOSE(unit=2)

       ! show some results; process tendencies of a selection of tracers
       ! ESS_lg_20130106+ all process tendencies are calculated in molecules g-1 s-1.
       ! They are recalculated to ppbv hr-1 using the following recalculation term:
       recalctend    =3600*amd*1e9/avo
       recalctendppmv=3600*amd*1e6/avo
	   
       OPEN(unit=2,file='output/xttend.out',status='UNKNOWN')
       WRITE(2,'(1a)') &
         'Canopy exchange model output: process tendencies'
       WRITE(2,'(1a)') &
         'Process tendencies in ppbv/ppmv hr-1, em=emissions, dd=dry deposition, ch=chemistry and df=diffusion'

	   ! ESS_20140416+ modified writing of radiation data for higher vertical resolution model 

	   ! to assign parname and parunit for the header of the file
	   ALLOCATE (parname(npar_out,nveglay*ntrac))
	   ALLOCATE (parunit(npar_out,nveglay*ntrac))
	   ALLOCATE (tendencies(Nstep,ntend,nveglay+1,ntrac))
	   
       ! assigning some of this header info and some further manipulations for output
 
	   ! ESS_lg_20140418+ tracer tendencies (only a selection, see tracer_init)
	   ! emission tendencies, only canopy layers
       npar=npar+1
       ixteemis=npar
	   nlevels(npar)=0
	   jjt_emis=0
	   DO jt=1,ntrac
     	 IF (latmbios_output(jt)) THEN
           IF (jt == idt_NO2 .or. jt == idt_ISOP .or. jt == idt_CO2 .or.  &
		       jt == idt_APIN .or. jt == idt_MYRC .or. jt == idt_HONO) THEN  ! MAQ_20170413+ added some of the terpenes) THEN
             jjt_emis=jjt_emis+1
		     ! Only canopy layers
		     DO jk=1,nveglay
			   nlevels(npar)=nlevels(npar)+1
			   dummy=trname(jt)
			   
			   ! special case; NOx is considered
			   IF (jt == idt_NO2) dummy='NOx'
			   
               iip=INDEX(dummy,' ')
 		       dummy2='em('
               dummy=dummy(1:iip-1)//dummy2
               iip=INDEX(dummy,' ')
               IF (zrefcan(1,jk) < 10.) THEN
			     write(dummy2,'(f3.1)')zrefcan(1,jk)
			   ELSE
			     write(dummy2,'(f4.1)')zrefcan(1,jk)
			   ENDIF
               dummy=dummy(1:iip-1)//dummy2
               iip=INDEX(dummy,' ')
 		       dummy3='m)'
               parname(npar,nlevels(npar))=dummy(1:iip-1)//dummy3
			   parunit(npar,nlevels(npar))='ppbv hr-1'
               IF (jt == idt_CO2) parunit(npar,nlevels(npar))='ppmv hr-1'
			 ENDDO

			 ! Recalculation of tendencies; special cases, NOx, and CO2 in ppmv
			 IF (jt == idt_NO2) THEN
  	   	       tendencies(:,1,2:nveglay+1,jjt_emis)= &
			     recalctend*(xteemis(:,2:nveglay+1,jt)+xteemis(:,2:nveglay+1,idt_NO))
			 ELSE IF (jt == idt_CO2) THEN
               tendencies(:,1,2:nveglay+1,jjt_emis)= &
			     recalctendppmv*xteemis(:,2:nveglay+1,jt)
			 ELSE
			   tendencies(:,1,2:nveglay+1,jjt_emis)=recalctend*xteemis(:,2:nveglay+1,jt)
			 ENDIF

		   ENDIF
		 ENDIF
	   ENDDO	

	   ! dry deposition tendencies, only canopy layers
       npar=npar+1
       ixtedryd=npar
	   nlevels(npar)=0
	   jjt_dryd=0
	   DO jt=1,ntrac
 	  	 IF (latmbios_output(jt)) THEN
           IF (jt == idt_O3 .or. jt == idt_HNO3 .or. jt == idt_NO2 .or. jt == idt_CO2 .or. jt == idt_CH2O .or.  &
		       jt == idt_APIN .or. jt == idt_MYRC .or. jt == idt_HONO) THEN  ! MAQ_20170413+ added some of the terpenes) THEN
             jjt_dryd=jjt_dryd+1
		     ! Only canopy layers
		     DO jk=1,nveglay
			   nlevels(npar)=nlevels(npar)+1
			   dummy=trname(jt)

			   ! special case; NOx is considered
			   IF (jt == idt_NO2) dummy='NOx'
			   
               iip=INDEX(dummy,' ')
 		       dummy2='dd('
               dummy=dummy(1:iip-1)//dummy2
               iip=INDEX(dummy,' ')
               IF (zrefcan(1,jk) < 10.) THEN
			     write(dummy2,'(f3.1)')zrefcan(1,jk)
			   ELSE
			     write(dummy2,'(f4.1)')zrefcan(1,jk)
			   ENDIF
               dummy=dummy(1:iip-1)//dummy2
               iip=INDEX(dummy,' ')
 		       dummy3='m)'
               parname(npar,nlevels(npar))=dummy(1:iip-1)//dummy3
               parunit(npar,nlevels(npar))='ppbv hr-1'
               IF (jt == idt_CO2) parunit(npar,nlevels(npar))='ppmv hr-1'
			 ENDDO

			 ! Recalculation of tendencies; special cases, NOx, and CO2 in ppmv
			 IF (jt == idt_NO2) THEN
  	   	       tendencies(:,2,2:nveglay+1,jjt_dryd)= &
			     recalctend*(xtedryd(:,2:nveglay+1,jt)+xtedryd(:,2:nveglay+1,idt_NO))
			 ELSE IF (jt == idt_CO2) THEN
 	   	       tendencies(:,2,2:nveglay+1,jjt_dryd)= &
			     recalctendppmv*xtedryd(:,2:nveglay+1,jt)
			 ELSE
  	  	       tendencies(:,2,2:nveglay+1,jjt_dryd)=recalctend*xtedryd(:,2:nveglay+1,jt)
             ENDIF
		   ENDIF
		 ENDIF
	   ENDDO

	   ! chemistry, all layers
       npar=npar+1
       ixtechem=npar
	   nlevels(npar)=0
	   jjt_chem=0
	   DO jt=1,ntrac
 	  	 IF (latmbios_output(jt)) THEN
           IF (jt == idt_O3 .or. jt == idt_ISOP .or. jt == idt_NO2 .or. jt == idt_CH2O .or.  &
		       jt == idt_APIN .or. jt == idt_MYRC .or. jt == idt_HONO) THEN  ! MAQ_20170413+ added some of the terpenes) THEN
 		     nlevels(npar)=nlevels(npar)+1
		     jjt_chem=jjt_chem+1
	         dummy=trname(jt)

	         ! special case; NOx is considered
			 IF (jt == idt_NO2) dummy='NOx'
			 
             iip=INDEX(dummy,' ')
 		     dummy2='ch('
             dummy=dummy(1:iip-1)//dummy2
             iip=INDEX(dummy,' ')
             IF (zref(1) < 10.) THEN
			   write(dummy2,'(f3.1)')zref(1)
		     ELSE
			   write(dummy2,'(f4.1)')zref(1)
		     ENDIF
             dummy=dummy(1:iip-1)//dummy2
             iip=INDEX(dummy,' ')
 		     dummy3='m)'
             parname(npar,nlevels(npar))=dummy(1:iip-1)//dummy3
             parunit(npar,nlevels(npar))='ppbv hr-1' ! MAQ_lg_20210824+ was ppbv s-1, should be in ppbv hr-1
			 
		     ! Canopy layers
		     DO jk=1,nveglay
               nlevels(npar)=nlevels(npar)+1
			   dummy=trname(jt)

			   ! special case; NOx is considered
			   IF (jt == idt_NO2) dummy='NOx'
			   
               iip=INDEX(dummy,' ')
 		       dummy2='ch('
               dummy=dummy(1:iip-1)//dummy2
               iip=INDEX(dummy,' ')
               IF (zrefcan(1,jk) < 10.) THEN
			     write(dummy2,'(f3.1)')zrefcan(1,jk)
			   ELSE
			     write(dummy2,'(f4.1)')zrefcan(1,jk)
			   ENDIF
               dummy=dummy(1:iip-1)//dummy2
               iip=INDEX(dummy,' ')
 		       dummy3='m)'
               parname(npar,nlevels(npar))=dummy(1:iip-1)//dummy3
               parunit(npar,nlevels(npar))='ppbv hr-1'
			 ENDDO

			 ! Recalculation of tendencies; special case, NOx, 
			 IF (jt == idt_NO2) THEN
  	   	       tendencies(:,3,1:nveglay+1,jjt_chem)= &
			     recalctend*(xtechem(:,1:nveglay+1,jt)+xtechem(:,1:nveglay+1,idt_NO))
			 ELSE
   		       tendencies(:,3,1:nveglay+1,jjt_chem)=recalctend*xtechem(:,1:nveglay+1,jt)
			 ENDIF
		   ENDIF
		 ENDIF
	   ENDDO

	   ! diffusion tendencies, all layers
       npar=npar+1
       ixtediff=npar
	   nlevels(npar)=0
	   jjt_diff=0
	   DO jt=1,ntrac
 	  	 IF (latmbios_output(jt)) THEN
           IF (jt == idt_O3 .or. jt == idt_ISOP .or. jt == idt_NO2 .or. jt == idt_CO2 .or. jt == idt_CH2O .or.  &
		       jt == idt_APIN .or. jt == idt_MYRC .or. jt == idt_HONO) THEN  ! MAQ_20170413+ added some of the terpenes) THEN
 		     nlevels(npar)=nlevels(npar)+1
		     jjt_diff=jjt_diff+1
	         dummy=trname(jt)

			 ! special case; NOx is considered
			 IF (jt == idt_NO2) dummy='NOx'

			 iip=INDEX(dummy,' ')
 		     dummy2='df('
             dummy=dummy(1:iip-1)//dummy2
             iip=INDEX(dummy,' ')
             IF (zref(1) < 10.) THEN
			   write(dummy2,'(f3.1)')zref(1)
		     ELSE
			   write(dummy2,'(f4.1)')zref(1)
		     ENDIF
             dummy=dummy(1:iip-1)//dummy2
             iip=INDEX(dummy,' ')
 		     dummy3='m)'
             parname(npar,nlevels(npar))=dummy(1:iip-1)//dummy3
             parunit(npar,nlevels(npar))='ppbv hr-1'                    ! MAQ_lg_20210824+ was ppbv s-1, should be in ppbv hr-1 
             IF (jt == idt_CO2) parunit(npar,nlevels(npar))='ppmv hr-1' ! MAQ_lg_20210824+ was ppmv s-1, should be in ppmv hr-1

		     ! Canopy layers
		     DO jk=1,nveglay
               nlevels(npar)=nlevels(npar)+1
			   dummy=trname(jt)

			   ! special case; NOx is considered
			   IF (jt == idt_NO2) dummy='NOx'

               iip=INDEX(dummy,' ')
 		       dummy2='df('
               dummy=dummy(1:iip-1)//dummy2
               iip=INDEX(dummy,' ')
               IF (zrefcan(1,jk) < 10.) THEN
			     write(dummy2,'(f3.1)')zrefcan(1,jk)
			   ELSE
			     write(dummy2,'(f4.1)')zrefcan(1,jk)
			   ENDIF
               dummy=dummy(1:iip-1)//dummy2
               iip=INDEX(dummy,' ')
 		       dummy3='m)'
               parname(npar,nlevels(npar))=dummy(1:iip-1)//dummy3
               parunit(npar,nlevels(npar))='ppbv hr-1'
               IF (jt == idt_CO2) parunit(npar,nlevels(npar))='ppmv hr-1'
			 ENDDO

			 ! Recalculation of tendencies; special cases, NOx, and CO2 in ppmv
			 IF (jt == idt_NO2) THEN
  	   	       tendencies(:,4,1:nveglay+1,jjt_diff)= &
			     recalctend*(xtediff(:,1:nveglay+1,jt)+xtediff(:,1:nveglay+1,idt_NO))
		     ELSE IF (jt == idt_CO2) THEN
 	   	       tendencies(:,4,1:nveglay+1,jjt_diff)= &
			     recalctendppmv*xtediff(:,1:nveglay+1,jt)
			 ELSE
			   tendencies(:,4,1:nveglay+1,jjt_diff)=recalctend*xtediff(:,1:nveglay+1,jt)
			 ENDIF
		   ENDIF
		 ENDIF
	   ENDDO		   

	   ! total tracer tendencies, all layers
       npar=npar+1
       ixtetend=npar
	   nlevels(npar)=0
	   jjt_xte=0
	   DO jt=1,ntrac
 	  	 IF (latmbios_output(jt)) THEN
           IF (jt == idt_O3 .or. jt == idt_ISOP .or. jt == idt_NO2 .or. jt == idt_CO2 .or. jt == idt_CH2O .or.  &
		       jt == idt_APIN .or. jt == idt_MYRC .or. jt == idt_HONO) THEN  ! MAQ_20170413+ added some of the terpenes
 		     nlevels(npar)=nlevels(npar)+1
		     jjt_xte=jjt_xte+1
	         dummy=trname(jt)

			 ! special case; NOx is considered
			 IF (jt == idt_NO2) dummy='NOx'

			 iip=INDEX(dummy,' ')
 		     dummy2='tot('
             dummy=dummy(1:iip-1)//dummy2
             iip=INDEX(dummy,' ')
             IF (zref(1) < 10.) THEN
			   write(dummy2,'(f3.1)')zref(1)
		     ELSE
			   write(dummy2,'(f4.1)')zref(1)
		     ENDIF
             dummy=dummy(1:iip-1)//dummy2
             iip=INDEX(dummy,' ')
 		     dummy3='m)'
             parname(npar,nlevels(npar))=dummy(1:iip-1)//dummy3
             parunit(npar,nlevels(npar))='ppbv hr-1'                    ! MAQ_lg_20210824+ was ppbv s-1, should be in ppbv hr-1
             IF (jt == idt_CO2) parunit(npar,nlevels(npar))='ppmv hr-1' ! MAQ_lg_20210824+ was ppmv s-1, should be in ppmv hr-1

		     ! Canopy layers
		     DO jk=1,nveglay
               nlevels(npar)=nlevels(npar)+1
			   dummy=trname(jt)

			   ! special case; NOx is considered
			   IF (jt == idt_NO2) dummy='NOx'

               iip=INDEX(dummy,' ')
 		       dummy2='tot('
               dummy=dummy(1:iip-1)//dummy2
               iip=INDEX(dummy,' ')
               IF (zrefcan(1,jk) < 10.) THEN
			     write(dummy2,'(f3.1)')zrefcan(1,jk)
			   ELSE
			     write(dummy2,'(f4.1)')zrefcan(1,jk)
			   ENDIF
               dummy=dummy(1:iip-1)//dummy2
               iip=INDEX(dummy,' ')
 		       dummy3='m)'
               parname(npar,nlevels(npar))=dummy(1:iip-1)//dummy3
               parunit(npar,nlevels(npar))='ppbv hr-1'
               IF (jt == idt_CO2) parunit(npar,nlevels(npar))='ppmv hr-1'
			 ENDDO

			 ! Recalculation of tendencies; special cases, NOx, and CO2 in ppmv
			 IF (jt == idt_NO2) THEN
  	   	       tendencies(:,5,1:nveglay+1,jjt_xte)=                                  &
			     recalctend*(xtediff(:,1:nveglay+1,jt)+xtediff(:,1:nveglay+1,idt_NO) &
	                        +xtechem(:,1:nveglay+1,jt)+xtechem(:,1:nveglay+1,idt_NO) &
	                        +xtedryd(:,1:nveglay+1,jt)+xtedryd(:,1:nveglay+1,idt_NO) &
	                        +xteemis(:,1:nveglay+1,jt)+xteemis(:,1:nveglay+1,idt_NO))
		     ELSE IF (jt == idt_CO2) THEN
 	   	       tendencies(:,5,1:nveglay+1,jjt_xte)=        &
			     recalctendppmv*(xtediff(:,1:nveglay+1,jt) &
	                            +xtedryd(:,1:nveglay+1,jt) &
	                            +xteemis(:,1:nveglay+1,jt))
			 ELSE
			   tendencies(:,5,1:nveglay+1,jjt_xte)=       &
			     recalctend*(xtediff(:,1:nveglay+1,jt)    &
				              +xtechem(:,1:nveglay+1,jt)  &
							  +xtedryd(:,1:nveglay+1,jt)  &
                              +xteemis(:,1:nveglay+1,jt)) 							  
			 ENDIF
		   ENDIF
		 ENDIF
	   ENDDO	
	   
       ! ESS_lg_20140423+ and writing the file 
	   WRITE(2,'(A10,A16,250A15)') &
             'istep','ldatltime'                        &
            ,parname(ixteemis,1:nlevels(ixteemis))      &
            ,parname(ixtedryd,1:nlevels(ixtedryd))      & 
            ,parname(ixtechem,1:nlevels(ixtechem))      &
            ,parname(ixtediff,1:nlevels(ixtediff))      & 
            ,parname(ixtetend,1:nlevels(ixtetend)) 			
       WRITE(2,'(A10,A16,250A15)') &
             '#','dd:mm:yy hh:min'                      &
            ,parunit(ixteemis,1:nlevels(ixteemis))      &
            ,parunit(ixtedryd,1:nlevels(ixtedryd))      & 
            ,parunit(ixtechem,1:nlevels(ixtechem))      &
            ,parunit(ixtediff,1:nlevels(ixtediff))      & 
            ,parunit(ixtetend,1:nlevels(ixtetend)) 			
       DO i=1,Nstep
         IF(mod(i,nprint) == 0) THEN
            WRITE(2,'(I10,2x,A14,250(1x,E14.6))')       &
              i,ldatltime(i)                            &

              ! emissions, only canopy layers
			 ,tendencies(i,1,2:nveglay+1,1:jjt_emis)    &

			 ! dry deposition, only canopy layers
             ,tendencies(i,2,2:nveglay+1,1:jjt_dryd)    &

             ! chemistry, all layers!
             ,tendencies(i,3,1:nveglay+1,1:jjt_chem)    &   

             ! vertical transport, also all layers!
             ,tendencies(i,4,1:nveglay+1,1:jjt_diff)    & 

             ! total tracer tendencies, also all layers!
             ,tendencies(i,5,1:nveglay+1,1:jjt_xte)  			 
         ENDIF
      END DO
      CLOSE(unit=2)
	   
	  DEALLOCATE (parname)
	  DEALLOCATE (parunit)
	  DEALLOCATE (tendencies)
	  
    ENDIF !     IF (l_xtsurf_veg_mlay) THEN

    ! MAQ_lg_20170526+ moved to here the writing of the HONO/NOx emission diagnostics since these terms are
	!   recalculated in the subroutine veg_mlay to also use the actual updated information about the 
	!   tracer concentrations
    IF (l_emis_bio_jNO3) THEN

      PRINT *,'Writing of output on HONO/NOx emissions from NO3 photolysis (emdep_emis_bio_jNO3)'
 
      ! ESS_20140416+ modified writing of radiation data for higher vertical resolution model 
      ! show some results; HONO and NOx emissions from nitrate photolysis
      OPEN(unit=2,file='output/HONO-NOxem_jNO3.out',status='UNKNOWN')
      WRITE(2,'(1a)') &
         'HONO and NOx emissions due to nitrate photolysis on the leaf surface'
      WRITE(2,'(1a)') &
         'Dry deposition of NO3, surface nitrate, throughfall, fraction of sunlit leaves and HONO and NOx emission fluxes'
      WRITE(2,'(1a)') &
         'For a throughfall of 1 all the accumulated nitrate is being removed'

 	  ! to assign parname and parunit for the header of the file
      ALLOCATE (parname(2,nveglay))
	  ALLOCATE (parunit(2,nveglay))
 
      ! assigning some of this header info and some further manipulations for output
      DO i=1,nveglay
	    dummy='HONOem('
        iip=INDEX(dummy,' ')
        IF (zrefcan(1,i) >= 10.) THEN
          write(dummy2,'(f4.1)')zrefcan(1,i)
	    ELSE
          write(dummy2,'(f3.1)')zrefcan(1,i) 
	    ENDIF 
        dummy=dummy(1:iip-1)//dummy2
        iip=INDEX(dummy,' ')
 		dummy3=')'
        parname(1,i)=dummy(1:iip-1)//dummy3
        parunit(1,i)='molec m-2 s-1'

	    dummy='NO2em('
        iip=INDEX(dummy,' ')
        IF (zrefcan(1,i) >= 10.) THEN
          write(dummy2,'(f4.1)')zrefcan(1,i)
	    ELSE
          write(dummy2,'(f3.1)')zrefcan(1,i) 
	    ENDIF 
        dummy=dummy(1:iip-1)//dummy2
        iip=INDEX(dummy,' ')
 		dummy3=')'
        parname(2,i)=dummy(1:iip-1)//dummy3
		parunit(2,i)='molec m-2 s-1'
      ENDDO

	  IF (l_xtsurf_veg_mlay) THEN
        WRITE(2,'(A10,A16,50A15)') &
          'istep','ldatltime','ddNO3','NO3s','prec.','throughfall','Fslsum',           &
          (parname(1,i),i=1,nveglay),(parname(2,i),i=1,nveglay)
        WRITE(2,'(A10,A16,50A15)') &
          '#','dd-mm-yy hh:min','molec m-2 s-1','molec m-2','mm','[0-1]','[0-1]'   &
          ,(parunit(1,i),i=1,nveglay),(parunit(2,i),i=1,nveglay)
        DO i=1,Nstep
          IF (mod(i,nprint) == 0) THEN
            WRITE(2,'(I10,2x,A14,50(1x,E14.6))') &
              i,ldatltime(i),ddNO3(i,1),NO3s(i,1),(prc(i)+prl(i))*1e3,cthru(i),fslsum(i),  &
             (HONO_emflux(i,ii),ii=1,nveglay),(NOx_emflux(i,ii),ii=1,nveglay)
          ENDIF
        END DO
      ELSE  ! "big-leaf" HONO and NOx emissions from NO3 photolysis
        WRITE(2,'(A10,4A15)') &
          'istep','ldatltime','ddNO3','NO3s','HONO_emflux','NOx_emflux'
        WRITE(2,'(A10,4A15)') &
          '#','dd-mm-yy hh:min','molec m-2 s-1','molec m-2','molec m-2 s-1','molec. m-2 s-1'
        DO i=1,Nstep
           IF (mod(i,nprint) == 0) THEN
             WRITE(2,'(I10,2x,A14,8(1x,E14.6))') &
               i,ldatltime(i),ddNO3(i,1),NO3s(i,1),cthru(i),fslsum(i),HONO_jNO3em(i),NOx_jNO3em(i)
           ENDIF
        END DO
      ENDIF

	  DEALLOCATE (parname)
	  DEALLOCATE (parunit)
      ! ESS_20140416- 

      CLOSE(unit=2)

    ELSE
      PRINT *,'Calculation of HONO/NOx emissions from NO3 photolysis switched off'
    END IF

    ! MAQ_lg_20170517+ included a subroutine to calculate the soil HONO emission flux
    IF (l_emis_soil_HONO) THEN

      ! ESS_20140416+ modified writing of HONO soil emission flux 
      OPEN(unit=2,file='output/soilHONOem.out',status='UNKNOWN')
      WRITE(2,'(1a)') &
         'soil HONO emission flux'
 	  WRITE(2,'(A10,A16, A15)') &
          'istep','ldatltime','HONO_slemflux'
      WRITE(2,'(A10,2x,A14,4A15)') &
          '#','dd-mm-yy hh:min','molec m-2 s-1'
      DO i=1,Nstep
        IF (mod(i,nprint) == 0) THEN
           WRITE(2,'(I10,2x,A14,1(1x,E14.6))') &
             i,ldatltime(i),HONO_slflux(i)
        ENDIF
      END DO
      CLOSE(unit=2)

    ELSE
      PRINT *,'Calculation of soil HONO emissions switched off'
    END IF	
    ! MAQ_lg_20170517+
    ! MAQ_lg_20170526+ end writing of all output for diagnostics
    ! ---------------------------------------------------------------------------------
	
    ! deallocation of 1D arrays
    DEALLOCATE(ndtgact)
    DEALLOCATE(iyear) 
    DEALLOCATE(imon) 
    DEALLOCATE(iday)
    DEALLOCATE(ihr)        
    DEALLOCATE(imin)        
    DEALLOCATE(isec)
    DEALLOCATE(Jday)
    DEALLOCATE(timeday)
    DEALLOCATE(ldatltime)

    DEALLOCATE (slf)     
    DEALLOCATE (slm)        
    DEALLOCATE (vgrat)       
    DEALLOCATE (forestfr)   
    DEALLOCATE (dm)             
    DEALLOCATE (dmbase)  ! MAQ_lg_20160817+    
    DEALLOCATE (lai) 
    DEALLOCATE (laibase) ! MAQ_lg_20160621+	
    DEALLOCATE (hc)             
    DEALLOCATE (disp)
    DEALLOCATE (zrefsl) ! ESS_20130408+	
    DEALLOCATE (loland)
    DEALLOCATE (scaling_rad)
    DEALLOCATE (globrad)
    DEALLOCATE (netrad)
    DEALLOCATE (zen) 
    DEALLOCATE (cossza)
    DEALLOCATE (tair)
    DEALLOCATE (tsurf)
    DEALLOCATE (tsoil) 
    DEALLOCATE (tcan)  ! MAQ_lg_20170505+ added Tcanopy
    DEALLOCATE (press) 
    DEALLOCATE (qm1)
    DEALLOCATE (qs)    ! MAQ_lg_20170720+
    DEALLOCATE (qsam)  ! MAQ_lg_20170720+
    DEALLOCATE (ws)
    DEALLOCATE (wsmax) 
    DEALLOCATE (rbvd) 
    DEALLOCATE (rvdsl)
    DEALLOCATE (NOemclassname)
    DEALLOCATE (cultiv) 
    DEALLOCATE (fertil) 
    DEALLOCATE (prc)
    DEALLOCATE (prl)
    DEALLOCATE (prectot) 
    DEALLOCATE (noemis_w) 
    DEALLOCATE (noemis_d)
    DEALLOCATE (cp)  
    DEALLOCATE (lsp)  
    DEALLOCATE (pulsing)
    DEALLOCATE (plsday) 
    DEALLOCATE (plsdurat)
    DEALLOCATE (pls) 
    DEALLOCATE (crf) 
    DEALLOCATE (NO_slflux)
	DEALLOCATE (NO2_slflux) ! MAQ_lg_20210727+ added NO2 soil (deposition) flux for extra diagnostics
    DEALLOCATE (NO_emflux)
    DEALLOCATE (RADON_slflux)
    DEALLOCATE (CO2_slflux) ! ESS_lg_20130503+ 	
    DEALLOCATE (latitude)   ! ESS_lg_20130817+
	DEALLOCATE (longitude)  ! ESS_lg_20130817+
    DEALLOCATE (cfml)  
    DEALLOCATE (cfncl) 
    DEALLOCATE (ril) 
    DEALLOCATE (cdnl) 
    DEALLOCATE (geopot_3d)  
    DEALLOCATE (geopot_3d_sl) ! ESS_lg_20130208+  
    DEALLOCATE (tvir) 
    DEALLOCATE (tvl) 
    DEALLOCATE (um1) 
    DEALLOCATE (vm1)
    DEALLOCATE (u_sl) ! ESS_lg_20131125+
    DEALLOCATE (az0) 
    DEALLOCATE (z0m)               
    DEALLOCATE (z0mslsn)     
    DEALLOCATE (zref)        
    DEALLOCATE (grvol)   
    DEALLOCATE (grmass) 
    DEALLOCATE (pdp) 
    DEALLOCATE (prhoa) 
    DEALLOCATE (MLH) ! ESS_lg_20130104+ 
    DEALLOCATE (zrahl) 
    DEALLOCATE (zrahveg) 
    DEALLOCATE (zrahslsn)
    DEALLOCATE (zustarl)
    DEALLOCATE (zustveg) 
    DEALLOCATE (zustslsn)
    DEALLOCATE (zfrl) 
    DEALLOCATE (zcvbs) 
    DEALLOCATE (zcvwat)
    DEALLOCATE (zcvs)
    DEALLOCATE (zcvw) 
    DEALLOCATE (zvgrat)
    DEALLOCATE (rh) 
    DEALLOCATE (rco_leaf)
    DEALLOCATE (rco_leaf_AGS) ! ESS_lg_20130516+	
	DEALLOCATE (rmesCO2)
	DEALLOCATE (rcutCO2) ! ESS_lg_20130516-
    DEALLOCATE (ccompCO2)! MAQ_lg_20170720+ added the CO2 compensation point
    DEALLOCATE (fws)  
    DEALLOCATE (fvpd) ! ESS_lg_20130424+	

    ! ESS_lg_20130117+ added for HONO/NOx source from nitrate photolysis
    DEALLOCATE (fslsum) 
    DEALLOCATE (ddNO3)   
    DEALLOCATE (NO3s)
    DEALLOCATE (cthru)
    DEALLOCATE (HONO_jNO3em)
    DEALLOCATE (NOx_jNO3em)
    ! ESS_lg_20130117- 

    ! 2D arrays
    DEALLOCATE (lad)
	DEALLOCATE (zrefcan)    ! ESS_lg_20140416+
    DEALLOCATE (zrefcan_hr) ! ESS_lg_20140416+
    DEALLOCATE (observ)
    DEALLOCATE (readdata)     
    DEALLOCATE (rvd)  
    DEALLOCATE (fsl)
    DEALLOCATE (fslbase) ! MAQ_lg_20160621+	
    DEALLOCATE (PAR)     ! MAQ_lg_20180112+			
    DEALLOCATE (u_veg) 
    DEALLOCATE (Kh)
    DEALLOCATE (cpold)  
    DEALLOCATE (lspold)
    DEALLOCATE (noemclass1)  
    DEALLOCATE (VOC_emfact)  
    DEALLOCATE (VOC_emflux)     
    DEALLOCATE (ISOP_emflux)
    DEALLOCATE (MONO_emflux) 
    DEALLOCATE (OVOC_emflux) 
    DEALLOCATE (HONO_emflux)
    DEALLOCATE (HONO_slflux)  ! MAQ_lg_20170517+
    DEALLOCATE (NOx_emflux)
    DEALLOCATE (NH3_emflux)
    DEALLOCATE (NH3_slflux)  ! MAQ_lg_20230329+ NH3
    DEALLOCATE (SQT_slflux)  ! MAQ_lg_20230601+ SQT
    DEALLOCATE (soilph)
    DEALLOCATE (rco_leaf_AGS_ml) ! ESS_lg_20130516+	

    ! DEALLocation of tracer paramters, 1D
    DEALLOCATE (trname)
    DEALLOCATE (moleweight)  
    DEALLOCATE (reactivity) 
    DEALLOCATE (henrycoeff) 
    DEALLOCATE (lvd_bigl)
    DEALLOCATE (lexist_GAS)  
    DEALLOCATE (diff)  
    DEALLOCATE (diffrb)
    DEALLOCATE (rsoil)
    DEALLOCATE (rwater)
    DEALLOCATE (rws) 
    DEALLOCATE (rsnow)
    DEALLOCATE (rmes)
    DEALLOCATE (rcut)
    DEALLOCATE (ccomp)
    DEALLOCATE (rj_max) 
    DEALLOCATE (medium)  
    DEALLOCATE (latmbios)
    DEALLOCATE (latmbios_emveg)
    DEALLOCATE (latmbios_emsoil)
    DEALLOCATE (latmbios_photo)
    DEALLOCATE (latmbios_output) ! ESS_lg_20140418+

    ! 2D arrays
    DEALLOCATE (pxtm1)    
    DEALLOCATE (pxtm1_obs)  
    DEALLOCATE (patmbiosflux)   
    DEALLOCATE (pcrf)
    DEALLOCATE (psurfflux)
	DEALLOCATE (pstomfluxO3) ! ESS_lg_20140509+
    DEALLOCATE (rj)
    DEALLOCATE (OHreact) ! ESS_lg_20130113+

    ! 3D arrays
    DEALLOCATE (pxtmveg)  
    DEALLOCATE (rj_veg)  
    DEALLOCATE (pvdveg)
    DEALLOCATE (xteemis)
    DEALLOCATE (xtedryd)
    DEALLOCATE (xtechem)
    DEALLOCATE (xtediff)        

  END SUBROUTINE emdep_xtsurf_vdiff

  ! --------------------------------------------------------------------------------- 
  
  SUBROUTINE emdep_xtsurf_tracer_init( l_xtsurf_veg_mlay_chem,           &
      Nstep, ntrac, trname, moleweight, reactivity, henrycoeff, medium,  &
      lo_derived, lexist_O3, lvd_bigl, lexist_GAS,                       &
      latmbios, latmbios_emveg, latmbios_emsoil, latmbios_photo,         &
      latmbios_output, rj_max, pxtm1, pxtmveg)                            ! ESS_lg_20140418+

    ! ESS_lg_20120718+ declaration of tracers considered in the emissions, 
    !   dry deposition and canopy exchanges simulations. 
    !   To extend the calculations to other trace gases/aerosols, these tracers and 
    !   aerosols and their properties should be included consistently here. 
    !   modifying: ntrac, trname, idt_trc, etc. 

    IMPLICIT NONE
    INTEGER :: i, ip
    INTEGER, INTENT(in) :: Nstep, ntrac
    LOGICAL, INTENT(in) :: l_xtsurf_veg_mlay_chem

    CHARACTER(LEN=20), INTENT(out), DIMENSION(ntrac) ::   &
      trname 

    ! declaration of tracer properties that are used to calculate their surface uptake resistances
    REAL(dp), INTENT(out), DIMENSION(ntrac) :: &
      moleweight  , & ! molecular weight, needed for estimating surface resistances 
      reactivity  , & ! reactivity coefficient [0:non-react., 0.1:semi react., 1:react.]  
      henrycoeff      ! henrycoeff coefficient [mol atm-1]

    LOGICAL, INTENT(out) :: & 
      lo_derived  , & ! true whenever SO2 and O3 are defined as tracers
      lexist_O3       ! true whenever O3 is present

    LOGICAL, INTENT(out), DIMENSION(ntrac) :: & 
      lvd_bigl    , & ! true for the tracer when the required parameters are defined
      lexist_GAS      ! true for the tracer when it is a gas

    ! mz_lg_20040721+ extra parameters of the aerosol dry deposition code.
    INTEGER, PARAMETER :: &
      AIR=             1, & ! index indicating tracer being a gas
      AEROSOL=         2    ! index indicating tracer being an aerosol
      
    INTEGER, INTENT(out),DIMENSION(ntrac) :: &
      medium  

    LOGICAL, INTENT(out),DIMENSION(ntrac) :: & 
      latmbios,       & ! true for atmosphere-biosphere calculations
      latmbios_emveg, & ! true for tracer vegetation emissions
      latmbios_emsoil,& ! true for tracer soil emissions
      latmbios_photo, & ! true for photolysis rates
	  latmbios_output   ! ESS_lg_20140418+

    REAL(dp), INTENT(out), DIMENSION(ntrac) :: &
      rj_max            ! maximum surface layer photolysis rates [s-1]

    REAL(dp), INTENT(out), DIMENSION(Nstep,ntrac) :: &
      pxtm1             ! surface layer mixing ratio ! mz_lg_20050719+

    REAL(dp), INTENT(out), DIMENSION(Nstep,nveglay,ntrac) :: &
      pxtmveg           ! canopy mixing ratio ! mz_lg_20050719+

    IF (.NOT.l_xtsurf_veg_mlay_chem) THEN
      idt_O3=1         ! tracer indices for ozone, etc.
      idt_HNO3=2         
      idt_NO=3           
      idt_NO2=4          
      idt_CO=5          
      idt_SO2=6          
      idt_H2O2=7         
      idt_ISOP=8         
      idt_APIN=9        
      idt_BPIN=10        
      idt_SQTERP=11     
      idt_PAN=12          
      idt_HCOOH=13      
      idt_NH3=14        
      idt_APINP1A=15   
      idt_RAD=16       ! adding more tracer indices when required, consistent 
                       ! with ntrac and tracer names above!      
      idt_CO2=17       ! ESS_lg_20130503+

      idt_COS=18       ! MAQ_lg_20190110+ added COS for Hyytiala tower analysis
	  
    ELSE  ! chemistry scheme based on the implementation of the CBM4 isoprene chemistry in the SCM

      idt_NOX=1
      idt_O3=2
      idt_CH4=3
      idt_CO=4
      idt_HNO3=5
      idt_H2O2=6
      idt_CH3O2H=7
      idt_CH2O=8
      idt_ALD2=9
      idt_PAR=10
      idt_OLE=11
      idt_ETH=12
      idt_PAN=13
      idt_ACET=14
      idt_ISOP=15
      idt_MGLY=16
      idt_ISOPRD=17
      idt_METHAC=18
      idt_MVK=19
      idt_MEK=20
      idt_MPAN=21
      idt_NTR=22
      idt_DMS=23
      idt_SO2=24
      idt_SO4=25
      idt_RAD=26
      idt_ISONTR=27
      idt_HCOOH=28
      idt_CH3CO2H=29
      idt_NH2=30
      idt_NH3=31
      idt_NH4=32
      idt_APIN=33
      idt_BPIN=34
      idt_SQTERP=35
      idt_HONO=36
      idt_CH2OHO2H=37
      idt_RCHOHO2H=38
      idt_CH3OH=39
      idt_CH3CN=40
      idt_SQTERP2B=41
      idt_SQTERP1B=42
      idt_MTTERP=43
      idt_HEXANE=44
      idt_BUTADIENE=45
      idt_TMBENZENE=46
      idt_NO=47
      idt_NO2=48
      idt_NO3=49
      idt_N2O5=50
      idt_HNO4=51
      idt_OH=52
      idt_HO2=53
      idt_CH3O2=54
      idt_C2O3=55
      idt_XO2=56
      idt_ROR=57
      idt_XO2N=58
      idt_RXPAR=59
      idt_BXO2N=60
      idt_MC3O3=61

      idt_APINP1A=62 ! added the a-pinene aerosol product to test in-canopy cycling of this aerosol  
      idt_CO2=63     ! ESS_lg_20130503+

      idt_ATERP=64   ! MAQ_lg_20170412+ added a-terpenine for the ATTO tower analysis
	  idt_LIMO =65   ! MAQ_lg_20170412+ added limonene for the ATTO tower analysis
      idt_MYRC =66   ! MAQ_lg_20170412+ added myrcene for the ATTO tower analysis

      idt_COS  =67   ! MAQ_lg_20190110+ added COS for Hyytiala tower analysis
	  
	ENDIF

    trname(:)='' ! ESS_lg_20130104+ to initialize trname

    ! the next tracers are included in the non-chemistry as well as the CBM4 tracer list
    ! when adding a tracer name be sure that the tracer index is defined
    trname(idt_O3)      ='O3    '
    trname(idt_HNO3)    ='HNO3  '
    trname(idt_NO)      ='NO    '
    trname(idt_NO2)     ='NO2   '
    trname(idt_CH4)     ='CH4   ' ! ESS_lg_20130417+ added to secure also initialization of methane
    trname(idt_CO)      ='CO    '
    trname(idt_SO2)     ='SO2   '
    trname(idt_H2O2)    ='H2O2  '
    trname(idt_ISOP)    ='ISOP  ' 
    trname(idt_APIN)    ='APIN  '
    trname(idt_BPIN)    ='BPIN  '
    trname(idt_SQTERP)  ='SQTERP'
    trname(idt_PAN)     ='PAN   '
    trname(idt_HCOOH)   ='HCOOH '
    trname(idt_NH3)     ='NH3   '
    trname(idt_APINP1A) ='APINP1A'
    trname(idt_RAD)     ='RADON '
    trname(idt_CO2)     ='CO2   ' ! ESS_lg_20130503+

    trname(idt_ATERP)     ='aTERP   ' ! MAQ_lg_20170412+
    trname(idt_LIMO)      ='LIMO    ' ! MAQ_lg_20170412+
    trname(idt_MYRC)      ='MYRC    ' ! MAQ_lg_20170412+
    trname(idt_COS)       ='COS     ' ! MAQ_lg_20190110+
	
    IF (l_xtsurf_veg_mlay_chem) THEN
      trname(idt_MVK)      = 'MVK'
      trname(idt_METHAC)   = 'METHAC'
      trname(idt_HONO)     = 'HONO'
	  trname(idt_CH2O)     = 'HCHO' ! ESS_lg_20150318+ added HCHO
	  trname(idt_NO3)      = 'NO3'  ! MAQ_lg_20170607+ added NO3
    ENDIF

    ! defining the tracer properties used in the calculation of the surface uptake resistances
    moleweight(:)=0._dp ! initialize

    moleweight(idt_O3)=48    ! molecular weight
    moleweight(idt_HNO3)=63  
    moleweight(idt_NO)=30    
    moleweight(idt_NO2)=46

    moleweight(idt_HONO)=47  ! MAQ_lg_20170524+  
	
    moleweight(idt_SO2)=64   
    moleweight(idt_H2O2)=34  
    moleweight(idt_ISOP)=68
    moleweight(idt_PAN)=121  
    moleweight(idt_HCOOH)=46
    moleweight(idt_NH3)=17

	moleweight(idt_CO2)=44  ! ESS_lg_20130503+
    moleweight(idt_CH2O)=32 ! ESS_lg_20150324+

    moleweight(idt_APIN)=136 ! MAQ_lg_20160610+
	moleweight(idt_BPIN)=136 ! MAQ_lg_20160610+

	moleweight(idt_ATERP)=136 ! MAQ_lg_20170412+
	moleweight(idt_LIMO)=136  ! MAQ_lg_20170412+
    moleweight(idt_MYRC)=136  ! MAQ_lg_20170412+
    moleweight(idt_COS)=60    ! MAQ_lg_20180717+ added COS
	
    reactivity(idt_O3)=1.    ! reactivity coefficient [0:non-react., 0.1:semi react., 1:react.]  
    reactivity(idt_HNO3)=1.  
    reactivity(idt_NO)=1.    
    reactivity(idt_NO2)=1.   

    reactivity(idt_HONO)=0.1  ! MAQ_lg_20170524+ (personal communications Putian Zhou 2017, see Wesely) 

    reactivity(idt_SO2)=0.   
    reactivity(idt_H2O2)=1.  
    reactivity(idt_ISOP)=0.  
    reactivity(idt_PAN)=0.1  
    reactivity(idt_HCOOH)=0.
    reactivity(idt_NH3)=0.  

    reactivity(idt_CO2)=0.  ! ESS_lg_20130503+
    reactivity(idt_CH2O)=1. ! ESS_lg_20150324+
	
    reactivity(idt_APIN)=0. ! MAQ_lg_20160610+
	reactivity(idt_BPIN)=0. ! MAQ_lg_20160610+

    reactivity(idt_ATERP)=0. ! MAQ_lg_20170412+
	reactivity(idt_LIMO)=0.  ! MAQ_lg_20170412+
    reactivity(idt_MYRC)=0.  ! MAQ_lg_20170412+
    reactivity(idt_COS)=0.   ! MAQ_lg_20180725+ 
	
	! ESS_lg_20150619+ modifications based on further feedback by Putian Zhou, having updated 
	! Wesely's values for Henry constant for a selection of tracers with the constants provided by
	! Rolf Sander, 2015, see also document Wesely_MESSy_Sander.pdf provided by Putian Zhou

    henrycoeff(idt_O3)=1.0e-2   ! ESS_lg_20150619+ modifications; henrycoeff coefficient [M atm-1]      
    henrycoeff(idt_HNO3)=8.9e4  ! ESS_lg_20150619+ was 1.E14  
    henrycoeff(idt_NO)=1.9e-3   ! ESS_lg_20150619+ was 2.E-3    
    henrycoeff(idt_NO2)=1.0e-2  ! ESS_lg_20150619+ was 2.E-3   

    henrycoeff(idt_HONO)=4.9e1  ! MAQ_lg_20170524+ personal communications Putian Zhou (H for 298.15) 

    henrycoeff(idt_SO2)=1.4     ! ESS_lg_20150619+ was 1.E5    
    henrycoeff(idt_H2O2)=9.2e4  ! ESS_lg_20150619+ was 1.E5   
    henrycoeff(idt_ISOP)=1.3e-2 ! MAQ_lg_20170413+ see draft paper Putian Zhou (Hyytiala BVOC exchange analysis)    
    henrycoeff(idt_PAN)=3.0     ! ESS_lg_20150619+ was 3.6     
    henrycoeff(idt_HCOOH)=8.9e3 ! ESS_lg_20150619+ was 4.E6
    henrycoeff(idt_NH3)=6.0e1   ! ESS_lg_20150619+ was 2.E4     

    henrycoeff(idt_CO2)=0.01    ! ESS_lg_20130503+
    henrycoeff(idt_CH2O)=3.2e3  ! ESS_lg_20150619+ was 1.e5 ! ESS_lg_20150324+ just an assigned high 

    henrycoeff(idt_APIN)=3.0e-2 ! MAQ_lg_20160610+
	henrycoeff(idt_BPIN)=1.6e-2 ! MAQ_lg_20160610+

    ! MAQ_lg_20170413+ see spreadsheet provided by Ana Serrano
    henrycoeff(idt_ATERP)=2.3e-2! MAQ_lg_20170413+ see draft paper Putian Zhou (Hyytiala BVOC exchange analysis);other minor MT
	henrycoeff(idt_LIMO)=4.9e-2 ! MAQ_lg_20170413+ see draft paper Putian Zhou (Hyytiala BVOC exchange analysis)
    henrycoeff(idt_MYRC)=8.9e-2 ! MAQ_lg_20170413+ see draft paper Putian Zhou (Hyytiala BVOC exchange analysis)
    henrycoeff(idt_COS)=2.1E-2  ! MAQ_lg_20180725+ 2.1E-02*EXP(3300.*(1./298.15-1./T)) ! [M/atm] type: L, ref: 2828
	
    ! definition of some tracer properties to distinguish between aerosols and tracers and 
    ! needed for the calculation of the surface uptake resistances for gaseous dry deposition

    medium(:)=AIR ! default assigning AIR (gas-phase) to all tracers
    medium(idt_APINP1A)=AEROSOL ! ESS_lg_20120718+ except of a-pinene aerosol product

    ! determining if O3 exist, if the tracer is a gas and if the Vd should be calculated
    lo_derived=.true.
    lexist_O3=.false.
    lexist_GAS(:)=.false.
    lvd_bigl(:)=.false.

    ! initialization of some tracer concentrations in mixing ratios
    pxtm1(:,:)=0.1e-15_dp
    DO i=1,ntrac
       ! mz_lg_20050719+ some tracer initialization
       IF (trname(i) == 'O3      ') pxtm1(:,i)=15.e-9_dp
       IF (trname(i) == 'CO      ') pxtm1(:,i)=80.e-9_dp
       IF (trname(i) == 'CH4     ') pxtm1(:,i)=1700.e-9_dp
       IF (trname(i) == 'NO      ') pxtm1(:,i)=2.e-12_dp
       IF (trname(i) == 'NO2     ') pxtm1(:,i)=0.1e-9_dp
       IF (trname(i) == 'ISOP    ') pxtm1(:,i)=5.0e-9_dp
       IF (trname(i) == 'NH3     ') pxtm1(:,i)=0.5e-9_dp
       IF (trname(i) == 'RADON   ') pxtm1(:,i)=0.4e-19_dp

       IF (trname(i) == 'CO2     ') pxtm1(:,i)=380.e-6_dp  ! ESS_lg_20130503+
       IF (trname(i) == 'COS     ') pxtm1(:,i)=400.e-12_dp ! MAQ_lg_20190111+
	   
       ! mz_lg_20050719-
    ENDDO

    ! mz_lg_20050719+ added setting of tracer settings for 
    ! atmosphere-biosphere exchanges
    latmbios(:)=.false. ! mz_lg_20050719+
    latmbios_emveg(:)=.false. 
    latmbios_emsoil(:)=.false. 
    latmbios_photo(:)=.false. 
	latmbios_output(:)=.false.

    ! defining first-order estimates of maximum (tropical) photolysis rates
    IF (l_xtsurf_veg_mlay_chem) THEN 
      rj_MAX(idt_H2O2)=1e-15 ! still to be defined
      rj_MAX(idt_HNO3)=1e-15 ! still to be defined
      rj_MAX(idt_NO2)=1.e-6_dp ! some first estimate
      rj_MAX(idt_N2O5)=1e-15 ! still to be defined
      rj_MAX(idt_CH2O)=1e-15 ! still to be defined
      rj_MAX(idt_O3)=1e-6    ! some first estimate
      rj_MAX(idt_CH3O2H)=1e-15 ! still to be defined
      rj_MAX(idt_HNO4)=1e-15 ! still to be defined
      rj_MAX(idt_NO3)=1e-15 ! still to be defined
      rj_MAX(idt_PAN)=1e-15 ! still to be defined
      rj_MAX(idt_ALD2)=1e-15 ! still to be defined
      rj_MAX(idt_ACET)=1e-15 ! still to be defined
      rj_MAX(idt_MGLY)=1e-15 ! still to be defined
      rj_MAX(idt_HONO)=1e-15 ! still to be defined
    ENDIF

    ! determining which tracers and aerosols should be considered in atmosphere-biosphere calculations
    print *,'Atmosphere-biosphere exchanges of the following species is simulated'
    DO i=1,ntrac

       IF (trname(i) == 'O3      ') lexist_O3 = .true.
       IF (medium(i) == AIR)        lexist_GAS(i) = .true.
       IF ((lexist_GAS(i)) .AND. (moleweight(i) > 0. ) .AND. &
            henrycoeff(i) > 0.) lvd_bigl(i)=.true.

       ip=INDEX(trname(i),' ')
       IF (ip > 1) THEN ! ESS_lg_20130104+ to secure only writing those tracers that have been assigned a name
         print *,'tracer # and name: ',i,trname(i)
             latmbios(i)=.true.
			 latmbios_output(i)=.true.  ! ESS_lg_20140416+
       ENDIF
       ! ESS_lg_20140418+ excluding some tracers in the output file(s)
	   IF (trname(i).eq.'CH4'.OR.trname(i).eq.'CO'.OR.    &
	       trname(i).eq.'SO2') latmbios_output(i)=.false.
		   
	   ! ESS_lg_20140418+ adding some extra ones
	   IF (i == idt_OH) THEN
	     trname(i)='OH'
	     latmbios_output(i)=.TRUE.
	   ENDIF
	   ! ESS_lg_20140418-
	   
       IF (latmbios(i)) pxtmveg(:,1:nveglay,i)=pxtm1(1,i)       

       IF (rj_MAX(i) > 0.) latmbios_photo(i)=.true.
	   
    ENDDO

  END SUBROUTINE emdep_xtsurf_tracer_init

  ! -------------------------------------------------------------------------------------

  SUBROUTINE emdep_update_jval(i, nstepday, nstepinit, jday0, latitude, rj) ! ESS_lg_20150310+ nstepinit

    ! calculation of J values based on implementation of update_jval subroutine
    ! of MECCA'a box model (Rolf Sander)

    IMPLICIT NONE

    INTEGER,  INTENT(in)    :: i, nstepday, nstepinit, jday0 ! timestep and # of timesteps per day
    REAL(dp), INTENT(in)    :: latitude
    REAL(dp), INTENT(inout) :: rj(:,:)

    INTRINSIC SIN, COS, EXP, LOG

    REAL(dp) :: DayReal, SoDecli, Latitu, SINPSI, PHOTON, FCT, DN
    REAL(dp), PARAMETER :: DUSK = 0.0721347 ! = 5.E-2/LOG(2.) = PHOTON at dusk
    REAL(dp), PARAMETER :: Cancer = 23.45 * (pi/180.)

    ! ESS_lg_20130119+ modified calculation of radiation properties
        ! Calculate radiation as a function of latitude, season and time of the day; 
    ! see also emdep_update_jval for these calculations as well as:
    ! http://education.gsfc.nasa.gov/experimental/all98invproject.site/pages/science-briefs/ed-stickler/ed-irradiance.html
    ! Day as a real value (at start of spring DayReal = 0):
    !DayReal = real(i-nstepinit)/real(nstepday)+Jday0  ! ESS_lg_20150310+ nstepinit
	DayReal = real(nstepinit+i)/real(nstepday)+Jday0  ! MAQ_lg_20211103+ bugfix, nstepinit+1 ! ESS_lg_20150310+ added nstepinit

    ! seasonal cycle, Solar declination angle:
    SoDecli = -Cancer * COS (2.*PI*(DayReal+10.)/365.25)
    ! ESS_lg_20130119-

    ! diurnal cycle of psi, the solar elevation angle
    Latitu = latitude * (pi/180.)
    SINPSI = SIN (Latitu) * SIN (SoDecli) &
   &            - COS (Latitu) * COS (SoDecli) * COS (2 * PI * DayReal)
    ! PHOTON is approximately the positive part of SINPSI but
    !   - avoid sharp switching on and off
    !   - add some light at dawn before sunrise and at dusk after sunset
    PHOTON = DUSK*LOG(1.+EXP(SINPSI/DUSK))
    ! DN is a day/night switch to set photolyses to about 0 at night
    FCT=50.
    DN = EXP(FCT*SINPSI)/(EXP(-FCT*SINPSI)+EXP(FCT*SINPSI))

    ! J values from PAPER model by Landgraf et al.
    ! J value as a function of PHOTON: J = A*exp(-B/(C+PHOTON))
    ! PHOTON=sin(PSI)=COS(THETA)
    ! THETA=zenith angle, PSI=90-THETA
    ! PHOTON/COS(LATITU) reaches about 1 on noon of first day of spring
    ! temp profile of atmos. is: data/prof.AFGL.midl.sum
    ! surface albedo :  0.00
    ! ozone column (Dobson) :300.00
    ! cloud cover: OFF
    ! paper model was started on: 00-01-18
    rj(i,idt_O3)     = DN*4.4916E-04*EXP(-3.5807E+00 /(3.2382E-01+PHOTON)) !J01
    rj(i,idt_H2O2)   = DN*1.8182E-05*EXP(-1.2565E+00 /(1.7904E-01+PHOTON)) !J03
    rj(i,idt_NO2)    = DN*1.2516E-02*EXP(-4.5619E-01 /(6.9151E-02+PHOTON)) !J04
    rj(i,idt_NO3)    = DN*2.2959E-02*EXP(-1.0873E-01 /(2.1442E-02+PHOTON)) !J05
    rj(i,idt_N2O5)   = DN*1.1375E-04*EXP(-1.1092E+00 /(1.7004E-01+PHOTON)) !J07
    rj(i,idt_HNO3)   = DN*3.4503E-06*EXP(-2.1412E+00 /(2.6143E-01+PHOTON)) !J08
    rj(i,idt_CH3O2H) = DN*1.2858E-05*EXP(-1.1739E+00 /(1.7044E-01+PHOTON)) !J09
    rj(i,idt_CH2O)   = DN*8.6933E-05*EXP(-1.3293E+00 /(1.6405E-01+PHOTON)) !J10
    rj(i,idt_HNO4)   = DN*2.1532E-05*EXP(-1.9648E+00 /(2.1976E-01+PHOTON)) !J91
    rj(i,idt_HONO)   = DN*2.9165E-03*EXP(-5.1317E-01 /(7.4940E-02+PHOTON)) !J92

    !print *,i,nstepinit,nstepday,Jday0,rj(i,idt_NO2)
    !read (*,*)

  END SUBROUTINE emdep_update_jval

END MODULE messy_emdep_xtsurf_box

!*****************************************************************************

PROGRAM emdep_xtsurf

  USE messy_emdep_xtsurf_box, ONLY: emdep_xtsurf_initialize, emdep_xtsurf_vdiff

  IMPLICIT NONE

  CALL emdep_xtsurf_initialize ! read CTRL namelist
  CALL emdep_xtsurf_vdiff      ! calculate online emissions and dry deposition

END PROGRAM emdep_xtsurf

!*****************************************************************************
