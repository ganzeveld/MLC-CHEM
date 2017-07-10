MODULE messy_emdep_xtsurf

  USE messy_emdep,        ONLY: modstr, modver

  ! mz_lg_20040430+ added the messy file with some general constants
  USE messy_main_constants_mem,   ONLY: g, dp, i4, pi, &
                                        avo=>N_A, amd=>M_air ! mz_lg_20050719+
  ! mz_lg_20040430- 

  IMPLICIT NONE

  ! mz_pj_20030903+
  CHARACTER(len=*), PARAMETER :: SUBMODSTR='MLC_CHEM'
  CHARACTER(len=*), PARAMETER :: SUBmodver = 'Summer 2015'
  ! mz_pj_20030903-

  LOGICAL :: lxtsurf                 ! global switch (internal)
  LOGICAL :: l_xtsurf_veg_mlay       ! global switch (namelist)
  LOGICAL :: l_xtsurf_veg_mlay_chem  ! global switch (namelist) ! mz_lg_20050725+
  LOGICAL :: l_xtsurf_veg_mlay_chem_reactions ! global switch (namelist) ESS_lg_20120920+
  LOGICAL :: l_xtsurf_veg_mlay_chem_photolysis ! global switch (namelist) ESS_lg_20120920+
  LOGICAL :: l_xtsurf_veg_mlay_ccomp ! global switch (namelist) ! ESS_lg_20120721+
  LOGICAL :: l_xtsurf_AGS            ! global switch (namelist) ! ESS_lg_20130516+
  
  LOGICAL, PUBLIC :: l_readdata = .false.
  
  CHARACTER(LEN=75), PUBLIC :: infilename

  REAL(dp), PUBLIC :: weight_pxtm1_obs

  REAL(dp), PUBLIC, PARAMETER :: ptmst_veg=7200._dp    ! mz_lg_20050719+ timestep canopy exch.
  REAL(dp), PUBLIC, PARAMETER :: dtime_mixing=43000._dp ! mz_lg_20050719+ turnover time
  REAL(dp), PUBLIC, PARAMETER :: hcmin=1._dp           ! mz_lg_20050719+ min. canopy height  

  INTEGER, PUBLIC, PARAMETER :: ncl_soilph=7    ! number of soil pH classes
  ! mz_LG_20021118+ declaration of the the maximum number of vegetation
  !     layers. 
  INTEGER, PUBLIC :: nveglay_hr=4  ! ESS_lg_20140401+ maximum number of layers for high res. calc. (radiation, VOC emis., etc.)
  INTEGER, PUBLIC :: nveglay=2     ! mz_lg_20050718+ number of canopy layers for exchange calc. (<=nveglay_hr)  
  ! Ags vegetation type: for the simulations of stomatal exchange as function of also the CO2 concentration
  INTEGER, PUBLIC :: Agstype=1     ! Ags vegetation type, C3; Agstype=1, C4; Agstype=2, coniferous forest; Agstype=3, tropical forest; Agstype=4
 
CONTAINS

  SUBROUTINE emdep_xtsurf_read_nml_ctrl(status, iou)

    ! XTSURF MODULE ROUTINE (CORE)
    !
    ! READ xtsurf NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES
    !
    ! Author: Patrick Joeckel, MPICH, Feb 2002
    ! Modified: Laurens Ganzeveld, MPICH, 31-10-2002

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER  :: substr='emdep_xtsurf_read_nml_ctrl'
    LOGICAL                      :: lex          ! file exists ?
    INTEGER                      :: fstat        ! file status

    NAMELIST /CTRL_XTSURF/     l_xtsurf_veg_mlay,      &
                               l_xtsurf_veg_mlay_chem, &  ! mz_lg_20050725+
                               l_xtsurf_veg_mlay_chem_reactions, & !  ESS_lg_20120920+
                               l_xtsurf_veg_mlay_chem_photolysis, & !  ESS_lg_20130120+
                               l_xtsurf_veg_mlay_ccomp,&  ! ESS_lg_20120721+
							   l_xtsurf_AGS,           &  ! ESS_lg_20130516+
                               Agstype,                &  ! ESS_lg_20150623+
                               l_readdata,             &
                               infilename,             &
                               weight_pxtm1_obs,       &
							   nveglay,                &  ! ESS_lg_20140808+
							   nveglay_hr
                                                           
    status = 1 ! ERROR ON RETURN

    ! INITIALIZE GLOBAL CONTROL VARIABLES
    lxtsurf=.false.
    l_xtsurf_veg_mlay=.false.
    l_xtsurf_veg_mlay_chem=.false.           ! mz_lg_20050725+ 
    l_xtsurf_veg_mlay_chem_reactions=.true.  ! ESS_lg_20120920+ default the model will calculate the gas-phase chemistry
                                             ! whenever l_xtsurf_veg_mlay_chem = .TRUE.
    l_xtsurf_veg_mlay_chem_photolysis=.true. ! ESS_lg_20120920+ default the model will calculate photolysis
                                             ! whenever l_xtsurf_veg_mlay_chem = .TRUE.
    l_xtsurf_veg_mlay_ccomp=.false.          ! ESS_LG_20120721+
	l_xtsurf_AGS=.false.                     ! ESS_LG_20130516+

	Agstype=1                                ! ESS_lg_20150623+
	
    l_readdata=.false.
    infilename=''
    weight_pxtm1_obs=1.
	nveglay=2                                ! ESS_lg_20140808+ default version with 2 canopy layers for exchange calc.
	nveglay_hr=4                             ! and 4 high resolution canopy layers

    CALL read_nml_open(lex, substr, iou, 'CTRL_XTSURF', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL_XTSURF, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL_XTSURF', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    lxtsurf = .true.

    CALL read_nml_close(substr, iou, modstr)
    
    WRITE(*,*)  'l_xtsurf_veg_mlay          = ',l_xtsurf_veg_mlay
    WRITE(*,*)  'l_xtsurf_veg_mlay_ccomp    = ',l_xtsurf_veg_mlay_ccomp        
    WRITE(*,*)  'l_xtsurf_veg_mlay_chem     = ',l_xtsurf_veg_mlay_chem
    WRITE(*,*)  'l_xtsurf_veg_mlay_chem_reactions = ',l_xtsurf_veg_mlay_chem_reactions
    WRITE(*,*)  'l_xtsurf_veg_mlay_chem_photolysis = ',l_xtsurf_veg_mlay_chem_photolysis
    WRITE(*,*)  'l_xtsurf_AGS               = ',l_xtsurf_AGS  
	WRITE(*,*)  'Agstype                    = ',Agstype  
    WRITE(*,*)  'l_readdata                 = ',l_readdata
    WRITE(*,*)  'Input file name            = ',infilename
    WRITE(*,*)  'weight_pxtm1_obs           = ',weight_pxtm1_obs
    WRITE(*,*)  'nveglay                    = ',nveglay	
    WRITE(*,*)  'nveglay_hr                 = ',nveglay_hr	
	
    status = 0 ! no ERROR

  END SUBROUTINE emdep_xtsurf_read_nml_ctrl

  !===========================================================================

  SUBROUTINE emdep_xtsurf_calcra ( nstep,  & ! ESS_lg_20120722+
               ckap, loland,               &
               cfml, cfncl,                &
               ril,  cdnl,                 &
               geopot_3d,  tvir, tvl,      &
               um1,  vm1,  az0, z0m,       &      
               prahl,   prahveg, prahslsn, &
               pustarl, pustveg, pustslsn)

    ! LG- Subroutine in which the aerodynamic resistance (Ra) is being
    !     calculated This resistance is used in the calculation of the
    !     dry deposition velocity.Ra (and also the friction velocity u*
    !     are being calculated for each surface cover fraction to consider
    !     pronounced differences in the surface roughness
    !
    !     Laurens Ganzeveld, 1998 (see for references the papers mentioned
    !     in the dry deposition model Vd_bl), modified for implementation
    !     in echam5, October, 2001
    ! --------------------------------------------------------------------

    ! mz_lg_20040430+ 
    ! Interface:
    ! ----------
    ! input 
    ! nstep     : # of timesteps ! ESS_lg_20120722+
    ! ckap      : von Karman constant [-]
    ! loland    : land-sea mask, logical
    ! cfml      : drag coefficient over land
    ! cfncl     : scaling term land
    ! ril       : Richardson number over land [-] 
    ! cdnl      : neutral drag coefficient over land
    ! geopot_3d : geopotential height
    ! tvir      : virtual temperature [K]
    ! tvl       : virtual temperature over land [K]   
    ! um1       : windspeed, u component [m s-1]
    ! vm1       : windspeed, v component [m s-1]
    ! az0       : echam's surface roughness for momentum [m]  
    ! z0m       : vegetation surface roughness [m]         
    ! 
    ! output 
    ! prahl     : aerodynamic resistance over land [s m-1]
    ! prahveg   : aerodynamic resistance over vegetation [s m-1]
    ! prahslsn  : aerodynamic resistance over bare soil/snow-ice [s m-1]
    ! pustarl   : friction velocity over land [m s-1]
    ! pustveg   : friction velocity over vegetation [m s-1]
    ! pustslsn  : friction velocity over bare soil/snow-ice [m s-1]

    IMPLICIT NONE 

    ! I/O
    INTEGER,  INTENT(in)  :: nstep ! ESS_lg_20120722+
    REAL(dp), INTENT(in)  :: ckap
    LOGICAL,  INTENT(in)  :: loland(:)
    REAL(dp), INTENT(in)  :: cfml(:),  cfncl(:),  &
      ril(:),  cdnl(:), geopot_3d(:), tvir(:),    &
      tvl(:),  um1(:), vm1(:), az0(:), z0m(:)    
    REAL(dp), INTENT(out) ::   prahl(:),  &
      prahveg(:), prahslsn(:), pustarl(:), pustveg(:), pustslsn(:)

    ! mz_lg_20040423-

    ! LG- internal constants

    REAL(dp), PARAMETER :: zkmkh=0.74, &
                           pepdu2=1.e-3 ! check with the value zepdu2 in vdiff!

    ! LG- local variables

    INTEGER  :: jl, klon  ! mz_lg_20040430+ modified
    REAL(dp) :: zsurf,zmonin,zoverl,zxzsurf,zxzref

    REAL(dp),  DIMENSION(:)  , ALLOCATABLE  :: zcdnveg, zcdnslsn, &
          cmveg, cmslsn, cml, psih

    ! INITIALIZATION
    ! mz_lg_20040423+ modified
    klon=nstep ! SIZE(loland) ! ESS_lg_20120722+
    ALLOCATE(zcdnveg(klon))
    ALLOCATE(zcdnslsn(klon))
    ALLOCATE(cmveg(klon))
    ALLOCATE(cmslsn(klon))
    ALLOCATE(cml(klon))
    ALLOCATE(psih(klon))

    DO jl=1,klon

       ! LG- calculation of drag coefficient and u*. A change with
       !     the previous calculations in echam4 is that the exchange
       !     coefficients, calculated in echam5 over land, ice and water,
       !     is being used to explicitly calculate the aerodynamic
       !     resistances over these surface
       ! ==================================================================
       ! mz_lg_20030102+ over land

       cml(jl)=cdnl(jl)*cfml(jl)/cfncl(jl)
       pustarl(jl)=SQRT(cml(jl))*SQRT(MAX(pepdu2,um1(jl)**2+ &
            vm1(jl)**2))

       ! LG- Computation of stability correction term, The stability
       !     correction functions are taken from Stull (page 383-385)
       !     (08-11-98) and are slightly different from those by Williams
       !     and Hicks et al., which were originaly being used in the dry
       !     deposition scheme.

       zsurf=(geopot_3d(jl)/g)
       if(zsurf.le.0.1e-10_dp) zsurf=1.e-10_dp ! avoid zero due to initialization ! mz_lg_20040923
       if (ril(jl).gt.0.) then

          ! LG- calculating the Monin-Obukhov lenght directly applying the
          !     formula given by Stull, 9.7.5k, page 386
          zmonin=(pustarl(jl)*((tvir(jl)+tvl(jl))/2.)*                      &
               SQRT(MAX(pepdu2,um1(jl)**2+vm1(jl)**2)))/(ckap*g*(tvir(jl)-  &
               tvl(jl)))
		  zoverl=zsurf/zmonin
          psih(jl)=-4.7*zoverl

		  !print *,'pri > 0',jl,pustarl(jl),tvir(jl),tvl(jl),&
          !    SQRT(MAX(pepdu2,um1(jl)**2+vm1(jl)**2)),zmonin,psih(jl)

       elseif (ril(jl).lt.0_dp) then
          zmonin=zsurf/ril(jl)
          zoverl=zsurf/zmonin
          zxzsurf=zkmkh*(1.-9.*(zoverl))**(0.5)
          zxzref=zkmkh
          psih(jl)= &
               (2.*LOG((1.+zxzsurf)/2.)+LOG((1.+zxzsurf**2.)/2.)-  &
               2.*ATAN(zxzsurf))-  & ! primitive function value for z
               (2.*LOG((1.+zxzref)/2.)+LOG((1.+zxzref**2.)/2.)- &
               2.*ATAN(zxzref))      ! primitive function value for zz
          ! mz_lg_20040923+
       else
          ! avoid division by zero in case ril(jl) = 0
          psih(jl)=1.0_dp
          ! mz_lg_20040923-
       endif

       prahl(jl)=MAX(1._dp,(1._dp/(pustarl(jl)*ckap))* &
            (LOG((geopot_3d(jl)/g)/REAL(az0(jl)))-psih(jl)))

       if (z0m(jl).gt.0..or.loland(jl)) then
          zcdnveg(jl)=(ckap/LOG(1.+geopot_3d(jl)/ &
               (g*MAX(0.02_dp,z0m(jl)))))**2
          zcdnslsn(jl)=(ckap/LOG(1.+geopot_3d(jl)/ &
               (g*0.005)))**2
          cmveg(jl)=zcdnveg(jl)*cfml(jl)/cfncl(jl)
          cmslsn(jl)=zcdnslsn(jl)*cfml(jl)/cfncl(jl)
          pustveg(jl)=SQRT(cmveg(jl))*sqrt(MAX(pepdu2, &
               um1(jl)**2+vm1(jl)**2))
  	      pustslsn(jl)=SQRT(cmslsn(jl))*sqrt(MAX(pepdu2, &
               um1(jl)**2+vm1(jl)**2))
          prahveg(jl)=MAX(1._dp,(1./(pustveg(jl)*ckap))* &
               (LOG((geopot_3d(jl)/g)/ &
               MAX(0.02_dp,z0m(jl)))-psih(jl)))
		  prahslsn(jl)=MAX(1._dp,(1./(pustslsn(jl)*ckap))* &
               (LOG((geopot_3d(jl)/g)/ &
               0.005)-psih(jl)))
       else
          cmveg(jl)=cml(jl)
          cmslsn(jl)=cml(jl)
          pustveg(jl)=pustarl(jl)
          pustslsn(jl)=pustarl(jl)
          prahveg(jl)=prahl(jl)
          prahslsn(jl)=prahl(jl)
       endif

    ENDDO

    DEALLOCATE(zcdnveg)
    DEALLOCATE(zcdnslsn)
    DEALLOCATE(cmveg)
    DEALLOCATE(cmslsn)
    DEALLOCATE(cml)
    DEALLOCATE(psih)

  END SUBROUTINE emdep_xtsurf_calcra

  ! mz_lg_20050718+ included a subroutine for the atmosphere-biosphere
  !      exchange calculations  
  !==============================================================================

  SUBROUTINE emdep_xtsurf_veg_mlay( nstep,                     & ! ESS_lg_20120722+
    latmbios, latmbios_emveg, latmbios_emsoil,                 &
    lvd_bigl, lexist_GAS, trname, ccomp, ptmst,                &
    l_emis, l_emis_bio_NO, l_emis_bio_VOC, l_emis_bio_jHNO3,   &
    l_drydep, loland, pfrl, pcvbs, pcvw, pvgrat, zrefsl,       & ! ESS_lg_20131125+ zrefls
	rahveg, ustveg, u_sl, u_veg,                               & ! ESS_lg_20131125+ u_sl
	netrad, rj, rj_veg, tsurf, press, qm1, rh_2m,              & ! ESS_lg_20120725+ added rj, rj_veg, press and qm1 for chemistry
    lai, lad, hc, rco_leaf, fws, diff, diffrb, rcut, rmes,     &
    rws, rsoil, rmesCO2, rcutCO2, grvol, grmass, pdp, prhoa,   & ! ESS_lg_20130516+ added rmesCO2 and rcutCO2
    no_slflux, isop_emflux, mono_emflux, ovoc_emflux, hono_emflux,        &
    nox_emflux, radon_slflux, co2_slflux,                      & ! ESS_lg_20130503+ added CO2
	fAPIN, fBPIN, fSQTERP,                                     & ! ESS_lg_20130119+ added the partitioning of monoterpene and sqterp emissions
    pxtm1, pxtmveg, pxtm1_obs,                                 & ! ESS_lg_20120722+ added pxtm1_obs
    weight_pxtm1_obs, Kh, patmbiosflux, pcrf, pvdveg, psurfflux, & ! ESS_lg_20120903+ weight_pxtm1_obs, ESS_lg_20120722+ Kh for diagnostics
    pstomfluxO3, xteemis, xtedryd, xtechem, xtediff, OHreact)      ! ESS_lg_20140509+ added pstomfluxO3 ESS_lg_20130113+
        
    !=====================================================================
    !    Program to calculate the atmosphere-biosphere trace gas and
    !    aerosol exchanges with a two-layer representation of the canopy
    !    In this subroutine the emission within the canopy and the dry 
    !    deposition process as a function of plant activity are considered. 
    !    Then the canopy top flux is calculated from the gradient between 
    !    the vegetation layers and the surface layer concentrations and 
    !    the eddy diffusivity. The diagnostic vegetation layer concentrations 
    !    are then updated within the chemistry routine to correct for 
    !    chemical production or destruction. 
    !
    !    See for more information about the code the papers by
    !    Ganzeveld, L. et al., Atmosphere-biosphere trace gas exchanges 
    !    simulated with a single-column model, J. Geophys. Res., 107, 2002.
    !    Ganzeveld, L., et al., The influence of soil-biogenic NOx 
    !    emissions on the global distribution of reactive trace gases: 
    !    the role of canopy processes, J. Geophys. Res., 107, 2002.
    !=====================================================================

    ! mz_lg_20050721+ 
    ! Interface:
    ! ----------
    ! input 
    ! nstep           : # of timesteps ! ESS_lg_20120722+
    ! latmbios        : logical indicating atmosphere-biosphere exchange calc.
    ! latmbios_emveg  : logical indicating vegetation emis. in this subroutine
    ! latmbios_emsoil : logical indicating soil emissions in this subroutine
    ! lvd_bigl        : logical indicating "big-leaf" dry deposition calculations
    ! lexist_GAS      : logical indicating if tracer is gas
    ! trname          : tracer name
    ! ccomp           : tracer compensation point
    ! ptmst           : timestep lenght
    ! l_emis          : general swith for emission calculations
    ! l_emis_bio_NO   : switch for soil-biogenic NO emissions
    ! l_emis_bio_VOC  : switch for biogenic VOC emissions
    ! l_emis_bio_jHNO3: switch for biogenic NOx emissions due to HNO3 photolysis
    ! l_drydep        : general swith for dry deposition calculations
    ! loland   : land-sea mask, logical
    ! pfrl     : land fraction
    ! pcvbs    : bare soil fraction
    ! pcvw     : wet skin fraction
    ! pvgrat   : vegetation fraction
	! zrefsl   : reference height surface layer [m]
    ! rahveg   : aerodynamic resistance over vegetation [s m-1]
    ! ustveg   : friction velocity over vegetation [m s-1]
	! u_sl     : wind speed in surface layer [m s-1]
    ! u_veg    : windspeed in canopy [m s-1]
    ! netrad   : net surface radiation [W m-2]
    ! rj       : surface layer photolysis [s-1]
    ! rj_veg   : canopy photolysis [s-1]
    ! tsurf    : surface temperature [K]
    ! press    : surface layer pressure ! ESS_lg_20120725+
    ! qm1      : surface layer moisture ! ESS_lg_20120725+
    ! rh_2m    : relative humidity [0-1]
    ! lai      : Leaf Area Index [m2 m-2]
    ! lad      : Leaf Area Density profile [m2 m-3]
    ! hc       : canopy height [m]
    ! rco_leaf : H2O leaf stomatal resistance [s m-1]
    ! diff     : diffusivity term to get tracer stomatal resistance
    ! diffrb   : term that corrects for difference in diff., Rb
    ! rcut     : cuticular resistance [s m-1]
    ! rmes     : mesophyll resistance [s m-1]
    ! rws      : wet skin resistance [s m-1]
    ! rsoil    : soil resistance [s m-1]
	! rmesCO2  : mesophyllic resistance CO2 ! ESS_lg_20130516+
	! rcutCO2  : cuticular resistance CO2   ! ESS_lg_20130516+
    ! fws      : soil moisture stress attenuation function [0-1]
    ! grvol    : volume surface layer [m^3]
    ! grmass   : mass in surface layer [g]
    ! pdp      : pressure thickness surface layer
    ! prhoa    : density surface layer [kg m-3]
    ! no_slflux  : soil NO emission flux [molec. m-2 s-1]
    ! isop_emflux : Isoprene emission fluxes [molec. m-2 s-1]
    ! mono_emflux : monoterpene emission fluxes [molec. m-2 s-1]
    ! ovoc_emflux : other VOC emission fluxes [molec. m-2 s-1]
    ! hono_emflux : HONO emissions due to HNO3 photodiss [molec. m-2 s-1]
    ! hono_emflux : NOx emissions due to HNO3 photodiss [molec. m-2 s-1]
    ! radon_slflux: Radon emission flux [atoms m-2 s-1]
	! co2_slflux  : CO2 soil respiration flux [molec. m-2 s-1]
    ! fAPIN       : fraction of monoterpene emissions being emitted in the form of alpha-pinene
    ! fBPIN       : fraction of monoterpene emissions being emitted in the form of beta-pinene
    ! fSQTERP     : the the amount of the sesquiterpene emissions relatve to the monoterpene emissions
    ! pxtm1       : gas and aerosol mixing ratios
    ! pxtm1_obs   : observed (surface layer) mixing ratios
    ! weight_pxtm1_obs: weight of considering observed tracer concentrations
    !
    ! input/output
    ! pxtmveg     : canopy- gas and aerosol mixing ratios
    ! psurfflux   : surface exchange flux [molec m-2 s-1]
	! pstomfluxO3 : O3 stomatal uptake flux  ! ESS_lg_20140509+
    !
    ! output
    ! Kh          : Eddy diffusivity for heat (and tracers) [m2 s-1]
    ! patmbiosflux: atmosphere-biosphere exchange flux [molec m-2 s-1]
    ! pcrf        : Canopy Reduction Factor [ratio]
    ! pvdveg      : canopy dry deposition velocities [cm s-1]
    ! xteemis     : emission tendency [molecules g-1 s-1]
    ! xtedryd     : dry deposition tendency [molecules g-1 s-1]
    ! xtechem     : chemical tendency [molecules g-1 s-1]
    ! xtediff     : vertical diffusion tendency [molecules g-1 s-1]
    ! OHreact     : OH reactivity [s-1]

    USE messy_emdep_mem

    IMPLICIT NONE 

    ! I/O
    INTEGER,  INTENT(in)  :: nstep ! ESS_lg_20120722+
    REAL(dp), INTENT(in)  :: fAPIN, fBPIN, fSQTERP
    CHARACTER(LEN=20), INTENT(in) :: trname(:)
    LOGICAL,  INTENT(in)  :: l_emis, l_emis_bio_NO, l_emis_bio_VOC, &
                             l_emis_bio_jHNO3, l_drydep
    LOGICAL,  INTENT(in)  :: latmbios(:), lvd_bigl(:), lexist_GAS(:),  &
                             loland(:)
    LOGICAL,  INTENT(inout) :: latmbios_emveg(:), latmbios_emsoil(:)

    REAL(dp), INTENT(inout)  ::                                   & ! ESS_lg_20130516+ changed to inout to allow overwriting rmes
      ptmst, pfrl(:), pcvbs(:), pcvw(:), pvgrat(:), zrefsl(:),    & ! ESS_lg_20131125+ zref
	  rahveg(:), ustveg(:), u_sl(:), u_veg(:,:),                  & ! ESS_lg_20131125+ u_sl
	  netrad(:), tsurf(:), rh_2m(:), lai(:),                      & 
      lad(:,:), hc(:), rco_leaf(:,:), diff(:), diffrb(:), rcut(:),& ! ESS_lg_20130516+ rco_leaf layer dependent
      rmes(:), rws(:), rsoil(:), rmesCO2(:), rcutCO2(:), fws(:),  & ! ESS_lg_20130516+
	  grvol(:), grmass(:), pdp(:), prhoa(:), no_slflux(:), isop_emflux(:,:), &
      mono_emflux(:,:), ovoc_emflux(:,:), hono_emflux(:,:),       &
      radon_slflux(:), co2_slflux(:), ccomp(:), pxtm1_obs(:,:),   & ! ESS_lg_20130503+ CO2, ESS_lg_20120722+
      rj(:,:), rj_veg(:,:,:), press(:), qm1(:), weight_pxtm1_obs    ! ESS_lg_20120903+ weight_pxtm1_obs ! ESS_lg_20120725+

    REAL(dp), INTENT(inout)  :: pxtm1(:,:), pxtmveg(:,:,:), psurfflux(:,:), & ! ESS_20120718+ pxtm1 (is updated)
                                pstomfluxO3(:,:), nox_emflux(:,:)             ! ESS_lg_20140509+

    REAL(dp), INTENT(out)    :: Kh(:,:), patmbiosflux(:,:), pcrf(:,:), & ! ESS_lg_20120722+ Kh
                                pvdveg(:,:,:), xteemis(:,:,:),         & ! ESS_lg_20130106+ added the proces tendecies
                                xtedryd(:,:,:), xtechem(:,:,:), xtediff(:,:,:), &
                                OHreact(:,:) ! ESS_lg_20130113+
    ! more declarations

    REAL(dp),  DIMENSION(:,:,:), ALLOCATABLE  :: &
         emisflux
    REAL(dp),  DIMENSION(:,:)  , ALLOCATABLE  :: &
         rstomx,     &
         rleaf,      &
         rleafw,     &
         rbveg,      &
         lad_veglay, &
         emisflux_tot, &
         fluxctop,   &
         zxtm1,      &
         zxtmveg,    &
         dc,         &
         Ccanopy        ! ESS_lg_20120721+
    REAL(dp),  DIMENSION(:)    , ALLOCATABLE  :: &
         zxtmveg_old, &
         emifacveg,  &
         grvol_veg,  &
         grmass_veg, &
         pdp_veg,    &
         prhoa_veg,  &
         fluxveg,    &
         terma,      &
         termb,      &
         termc,      &
         termd,      &
         terme,      &
         em,         &
         pz,         &
         rahcan,     &
         rmesophyll

	! ESS_lg_20140423+ declarations for new numerical system for > 2 layers
    REAL(dp),  DIMENSION(:)    , ALLOCATABLE  :: &
	     termf,        &
         termw,        &
         termv,        &
         termwa,       &
         termvb,       &
	     pz_edge,      &
         zxtm,         & 
         ts,           &
	     trid_a,       &
	     trid_b,       &
         trid_b2,      &
	     trid_c,       &
	     trid_d,       &
	     trid_x,       &
	     NN,		   &
	     NM,		   &
	     MN,		   &
	     MM,		   &
	     pre_a,		   &
	     pre_b,		   &
	     pre_c	
		 		 
    ! mz_lg_20020115 local declarations

    INTEGER :: jl, jt, jk, jjk, klon, ntrac, it, nstep_sub, wp_sw,wp_sw_max, wp_int ! ESS_lg_20140423+ added terms
                   
    REAL(dp):: dz, dz1, dz2, dz3, ptmst_sub, ts_dd, ts_em,   &
               ts_turb, mass_old, mass_new, dmass, dmass_em, &
               dmass_dd, zxtm1_old, rtot, Khmin, ddfrac,     &
               rsveg_lay, rswet, dzsl

    REAL(dp), parameter :: stomblock=0._dp ! no stomatal blocking due to waterfilm/droplets

    ! mz_lg_20061021+ parameters to calculate aerosol dry deposition
    LOGICAL :: lvdaer(SIZE(latmbios))
    INTEGER,  parameter :: nmod=1, jmod=1
    REAL(dp):: rx1,rx2,um10,s,sc,cunning,sedspeed,vb_veg,vim_veg,st_veg,re,    &
         eff,rdrop,zdrop,qdrop,vkd_veg,vkc_veg,diffc,relax,eps,phi,alpha1,     &
         alpharat,vk1,vk2,zm,zmvd,zn,zdm,zdn,zdlnr
    REAL(dp):: zalphae(SIZE(loland)),zbeta(SIZE(loland)),zalpha(SIZE(loland))
    ! mz_lg_20031110+ added for some modifications of the vegetation surface
    !     resistance based on the paper by Gallagher et al., JGR, 2002. The
    !     paper shows some more details about the parameterization by 
    !     Slinn, Atmos. Environment, Vol 16, 1785-1794, 1982
    REAL(dp)    :: vin_veg,zrebound
    REAL(dp)    :: zrint
    REAL(dp)    :: zr(SIZE(loland))               ! Aerosol Radius of the respective tracer at klev
    REAL(dp)    :: densaer(SIZE(loland),nmod)
    !--- Assign values to used constants:
    REAL(dp), PARAMETER :: fln10=2.302585
    REAL(dp), PARAMETER :: w2pi=2.506638
    REAL(dp), PARAMETER :: zg=9.8e2         ! ESS_lg_20120721+ g cm s-2; changed to zg; important to avoid 
                                            ! use of wrong g in other calculations; cm s-2 at sea level
    REAL(dp), PARAMETER :: dynvisc=1.789e-4 ! g cm-1 s-1
    REAL(dp), PARAMETER :: cl=0.066*1e-4    ! mean free path [cm] (particle size also in cm)
    REAL(dp), PARAMETER :: bc= 1.38e-16     ! boltzman constant [g cm-2 s-1 K-1] (1.38e-23 J deg-1)
    REAL(dp), PARAMETER :: kappa=1.         ! shapefactor
    REAL(dp), PARAMETER :: visc=0.15        ! molecular viscocity [cm2 s-1]
    REAL(dp), PARAMETER :: vkar=0.40        ! von karman constant
    ! mz_lg_20031014+ added for the Vd aerosol over vegetation
    ! mz_lg_20040602+ modified
    REAL(dp), PARAMETER :: ZAS    = 10.E-6*1.E2  ! um -> CM, see paper Gallagher and Slinn, 1982
                ! here the smallest collector size is set at 10 um
                ! Those are the values for Slinn's 82 model (see Table 1)
    REAL(dp), PARAMETER :: zAL    = 1.E-3*1.E2   ! mm -> CM, see paper Gallagher and Slinn, 1982
                ! here the largest collector size is set at 1 mm  
    ! mz_lg_20061021

    ! INITIALIZATION
    klon=nstep ! SIZE(loland) ! ESS_lg_20120722+
    ntrac=SIZE(latmbios)

    ALLOCATE(emisflux(klon,nveglay,ntrac))

    ALLOCATE(rstomx(klon,ntrac))
    ALLOCATE(rleaf(klon,ntrac))
    ALLOCATE(rleafw(klon,ntrac))
    ALLOCATE(lad_veglay(klon,nveglay))
    ALLOCATE(fluxctop(klon,ntrac))
    ALLOCATE(emisflux_tot(klon,ntrac))
    ALLOCATE(zxtm1(klon,ntrac))

    ALLOCATE(zxtmveg(nveglay,ntrac))

    ALLOCATE(dc(klon,ntrac))

    ALLOCATE(rbveg(nveglay,ntrac))
    ALLOCATE(Ccanopy(nveglay,ntrac)) ! ESS_lg_20120721+

    ALLOCATE(rmesophyll(ntrac))

    ALLOCATE(zxtmveg_old(nveglay))
    ALLOCATE(emifacveg(nveglay))
    ALLOCATE(grvol_veg(nveglay))
    ALLOCATE(grmass_veg(nveglay))
    ALLOCATE(pdp_veg(nveglay))
    ALLOCATE(prhoa_veg(nveglay))

    ALLOCATE(fluxveg(nveglay+1))
    ALLOCATE(terma(nveglay+1))
    ALLOCATE(termb(nveglay+1))
    ALLOCATE(termc(nveglay+1))
	ALLOCATE(termd(nveglay+1))
	ALLOCATE(terme(nveglay+1))
    ALLOCATE(termf(nveglay+1))
    ALLOCATE(em(nveglay+1))
    ALLOCATE(pz(nveglay+1))

    ALLOCATE(rahcan(klon))

	! ESS_lg_20140423+ 
    ALLOCATE(pz_edge(nveglay+1))
	ALLOCATE(trid_a(nveglay+1))
	ALLOCATE(trid_b(nveglay+1))
	ALLOCATE(trid_b2(nveglay+1))
	ALLOCATE(trid_c(nveglay+1))
	ALLOCATE(trid_d(nveglay+1))
	ALLOCATE(trid_x(nveglay+1))
	ALLOCATE(pre_a(nveglay+1))
	ALLOCATE(pre_b(nveglay+1))	
	ALLOCATE(pre_c(nveglay+1))
	! ESS_lg_20140423-
	
    ! initializations
    patmbiosflux=0._dp
    pcrf=0._dp
    prhoa_veg(:)=prhoa(1)

    DO jt = 1, ntrac
      IF ( latmbios(jt) ) THEN
         ! mz_lg_20050725+ assigning the mesophyllic resistance
         !    (rmes = intent(in)!)
         rmesophyll(jt)=rmes(jt)

         ! mz_lg_20061021+ if the tracer is not a gas, it must be aerosol
         ! defining the aerosol dry deposition switch
         lvdaer(jt)=.FALSE. ! ESS_lg_20120718+
         IF (l_drydep .AND. .NOT. lexist_GAS(jt)) lvdaer(jt)=.TRUE.
      ENDIF
    ENDDO

    ! mz_lg_20050725+ re-assigning the tracer mesophyllic resistance
    !     of NO2 and NO that were calculated in the big-leaf deposition
    !     scheme from the ozone stomatal resistance
    IF (idt_NO2 > 0) &             
       rmesophyll(idt_NO2)=1._dp    ! see emdep_xtsurf_calc_rs
    IF (idt_NO > 0) &             
       rmesophyll(idt_NO)=500._dp   ! see emdep_xtsurf_calc_rs

    ! mz_lg_20050725+ determining the LAD for the nveglay layers
    DO jk=1,nveglay  
       lad_veglay(:,jk)=0._dp
       ! summing the LAD
       DO jjk=(jk-1)*nveglay_hr/nveglay+1,jk*nveglay_hr/nveglay
         lad_veglay(:,jk)=lad_veglay(:,jk)+lad(:,jjk)
       END DO
    END DO
    ! mz_lg_20050725-

    ! mz_lg_20050718+ assigning canopy emission fluxes

    IF (l_emis) THEN

      emisflux=0._dp
      emisflux_tot=0._dp
      xteemis=0._dp  ! ESS_lg_20130106+ emission tendency
          
      DO jl=1,klon

         IF (idt_NO > 0 .AND. l_emis_bio_NO) THEN
           latmbios_emsoil(idt_NO)=.TRUE.                
                    emisflux(jl,nveglay,idt_NO)=no_slflux(jl)  ! NO in molec. m-2 s-1

! MAQ_lg_20160922+ soil emission fluxes in surface layer to show difference in big-leaf versus ml approach
!                    emisflux(jl,1,idt_NO)=no_slflux(jl)  ! NO in molec. m-2 s-1

		 END IF

         ! ESS_lg_20120726+ added radon to test turbulent transport
         IF (idt_RAD > 0) THEN
           latmbios_emsoil(idt_RAD)=.TRUE.                
                    emisflux(jl,nveglay,idt_RAD)=radon_slflux(jl)  ! Radon in atoms m-2 s-1
         END IF
         ! ESS_lg_20120726-

		 ! ESS_lg_20130503+ added CO2 to test canopy exchange
         IF (idt_CO2 > 0) THEN
           latmbios_emsoil(idt_CO2)=.TRUE.                
              emisflux(jl,nveglay,idt_CO2)=co2_slflux(jl)  ! CO2 respiration flux in molec. m-2 s-1
		 END IF
         ! ESS_lg_20130503-
		 
         ! lg-   since the isoprene emission flux is already a grid average
         !       emission flux since all the controlling parameters are grid 
         !       averages, e.g. dry matter, for the bulkveg or veg_mlay routine, 
         !       the grid average isoprene emission flux is scaled with the
         !       relative fraction of vegetation and wet skin fraction, in      
         !       order to get an overall emission flux which resembles that
         !       calculated in the subroutine emdep_vocemis.f 
    
         ! lg-   isop in molec. m-2 s-1, see for the recalculation from the
         !       basic flux to proper units the subroutine vocemis

         IF ((pvgrat(jl)+pcvw(jl)) > 0._dp) THEN
            IF (idt_ISOP > 0 .AND. l_emis_bio_VOC) THEN
              latmbios_emveg(idt_ISOP)=.TRUE.  
              emisflux(jl,1:nveglay,idt_ISOP)=isop_emflux(jl,1:nveglay)/ &
                    (pvgrat(jl)+pcvw(jl))
            END IF

            ! mz_lg_20050810+ added alpha and beta pinene taking an even
            ! partioning
            IF (idt_APIN > 0 .AND. l_emis_bio_VOC) THEN
              latmbios_emveg(idt_APIN)=.TRUE.  
              emisflux(jl,1:nveglay,idt_APIN)=fAPIN*mono_emflux(jl,1:nveglay)/     &
                        (pvgrat(jl)+pcvw(jl))
            END IF

            IF (idt_BPIN > 0 .AND. l_emis_bio_VOC) THEN
              latmbios_emveg(idt_BPIN)=.TRUE.  
              emisflux(jl,1:nveglay,idt_BPIN)=fBPIN*mono_emflux(jl,1:nveglay)/     &
                        (pvgrat(jl)+pcvw(jl))
            END IF

            ! mz_lg_20050810+ sesquiterpene fluxes 
            IF (idt_SQTERP > 0 .AND. l_emis_bio_VOC) THEN
              latmbios_emveg(idt_SQTERP)=.TRUE.  
              emisflux(jl,1:nveglay,idt_SQTERP)=fSQTERP*mono_emflux(jl,1:nveglay)/ &
                       (pvgrat(jl)+pcvw(jl))
            END IF

            ! mz_lg_20050807+ added the NOx and HONO source through photolysis
            ! deposited HNO3. It is assumed that the NOx is emitted in the
            ! form of NO2

            ! NOx/NO2
            IF (idt_NO2 > 0 .AND. l_emis_bio_jHNO3) THEN
              latmbios_emsoil(idt_NO2)=.TRUE.                
              emisflux(jl,1:nveglay,idt_NO2)=nox_emflux(jl,1:nveglay)/ &
                       (pvgrat(jl)+pcvw(jl))
            END IF

            ! HONO
            IF (idt_HONO > 0 .AND. l_emis_bio_jHNO3) THEN
              latmbios_emsoil(idt_HONO)=.TRUE.                
              emisflux(jl,1:nveglay,idt_HONO)=hono_emflux(jl,1:nveglay)/ &
                       (pvgrat(jl)+pcvw(jl))
            END IF

         END IF ! IF ((pvgrat(jl)+p...
      END DO ! DO jl=1,klon
    ENDIF ! IF (l_emis)
 
    ! mz_lg_20050718+ calculating dry deposition from the dry deposition 
    ! velocity within the canopy as a function of uptake resistances of 
    ! the vegetation and wet skin fraction 

    IF (l_drydep) THEN
 
       xtedryd=0._dp  ! ESS_lg_20130106+ dry deposition tendency

       pvdveg=1.e-20_dp ! mz_lg_20050721+ initializing pvdveg also since 
                        ! it is needed to calculate the dry deposition timescale
       DO jk=1,nveglay
         DO jt=1,ntrac
           ! mz_lg_20061021+ gas-phase dry deposition velocities
           IF ( (latmbios(jt)) .AND. (lvd_bigl(jt)) ) THEN
             DO jl=1,klon
               rtot=0._dp

               ! ESS_lg_20130516+
               IF (l_xtsurf_AGS.AND.jt.EQ.idt_CO2) THEN
			      rmesophyll(jt)=rmesCO2(jl)
				  rcut(jt)=1e5 ! rcutCO2(jl) ! ESS_lg_20130619+ large rcut to avoid too large uptake (also at night)
			   ENDIF
               ! ESS_lg_20130516-			   

               IF (hc(jl) > hcmin .AND. &
                  (pvgrat(jl)+pcvw(jl)) >= 0._dp) THEN

                  rstomx(jl,jt)=diff(jt)*rco_leaf(jl,jk)/fws(jl) ! ESS_lg_20130516+ rco_leaf layer dependent

                  ! lg-    calculation of quasi-laminar boundary layer resistance
                  !        for canopy layers 

                  rbveg(jk,jt)=diffrb(jt)*180._dp*   &
                       (0.07_dp/MAX(1.e-10_dp,u_veg(jl,jk)))**0.5 
                               ! mz_lg_20050725+ layer index changed: u_veg!

                  rleaf(jl,jt)=(1._dp/((2._dp/rcut(jt)) +              &  ! dry vegetation, both sides
                        (1._dp/(rstomx(jl,jt)+rmesophyll(jt)))))
 						
 				  !     LG- 112003, modified by introducing the calculation of the leaf resistance
                  !     for the wet vegetation including the possibility to study the role of the
                  !     stomatal blocking for the exchanges                  

                  rleafw(jl,jt)=                                         &
                         (1._dp/((2._dp/rws(jt))+                        &
                                 (1._dp/(rws(jt)))*stomblock+            &
                                 (1._dp/(rstomx(jl,jt)+rmesophyll(jt)))* &
                                 (1.-stomblock)))

                  ! ESS_lg_20131125+ added to consider role of compensation point
				  IF (ccomp(jt)>0.) THEN
				    rleaf(jl,jt)=-99.99
					rleafw(jl,jt)=-99.99
				  ENDIF
				  ! ESS_lg_20131125-

				  ! calculation of in-canopy turbulent resistance needed to 
                  ! consider soil deposition
                  rahcan(jl)=MAX(1._dp,14.*lai(jl)*hc(jl)/ &
                           ustveg(jl))

                  ! understorey
                  IF (jk == nveglay) THEN

                     ! LG-   determining the reference height 
                     pz(jk)=0.25_dp*hc(jl)

                     IF (rleaf(jl,jt) == 0._dp) cycle

                     ! lg-    calculation of the total uptake resistance of the vegetated and
                     !        wet skin fraction from the specific land cover type resistances
                     !        and the land cover fractions, this gives the total uptake
                     !        velocities for these fractions, reflecting solely the within
                     !        canopy uptake processes. the turbulent transport is considered
                     !        in the calculations of the canopy top fluxes. the bulk dry 
                     !        deposition velocities are scaled with the LAD and radiation 
                     !        profile for the two layers in order to scale the "big leaf" 
                     !        dry deposition velocity accounting for the biomass distribution 
                     !        and the larger removal rates for the sunlit leaves


                     ! LG-    for the lowest canopy layer, the limiting turbulent transport from
                     !        the reference height to the soil surface is also considered.

                     IF (rleaf(jl,jt) > 0._dp) THEN  ! dry vegetation
                        rsveg_lay=1._dp/((1._dp/                          & 
                            (rsoil(jt)+rahcan(jl)*pz(jk)/hc(jl)))+        &
                            (1._dp/(rbveg(jk,jt)+rleaf(jl,jt))/           & ! ESS_lg_20140703+ bug fix, added ) after rleaf(jl,jt)
                             MAX(1.e-5_dp,lad_veglay(jl,jk)*lai(jl))))      ! ESS_lg_20140703+ bug fix, remove here ) after lai(jl)
					 ELSE
                        rsveg_lay=1.e10
                     END IF
                     IF (rleafw(jl,jt) > 0._dp) THEN   ! ESS_LG_20120721+
                       rswet=1._dp/((1._dp/                               & 
                            (rsoil(jt)+rahcan(jl)*pz(jk)/hc(jl)))+        &
                            (1._dp/(rbveg(jk,jt)+rleafw(jl,jt))/          & ! ESS_lg_20140703+ bug fix, added ) after rleaf(jl,jt)
                             MAX(1.e-5_dp,lad_veglay(jl,jk)*lai(jl))))      ! ESS_lg_20140703+ bug fix, remove here ) after lai(jl)
                     ELSE
                       rswet=1.e10
                     ENDIF

                     ! LG-     end

                  ELSE ! crown- and other canopy layers

                      IF (rleaf(jl,jt) > 0._dp) THEN
                        rsveg_lay=(rbveg(jk,jt)+rleaf(jl,jt))/      &
                                 MAX(1.e-5_dp,lad_veglay(jl,jk)*lai(jl))
                      ELSE
                        rsveg_lay=1.e10
                      ENDIF
                      IF (rleafw(jl,jt) > 0._dp) THEN  ! ESS_lg_20120721+
                        rswet=(rbveg(jk,jt)+rleafw(jl,jt))/           &
                                 MAX(1.e-5_dp,lad_veglay(jl,jk)*lai(jl))
                      ELSE
                        rswet=1.e10
                      ENDIF

				  END IF ! IF (jk == nveglay) THEN

                  ! mz _lg_20050721+ calculation of surface uptake resistance 
                  rtot=1._dp/         &
                     ((pvgrat(jl)/(pvgrat(jl)+pcvw(jl)))*(1./(rsveg_lay))+  &
                      (pcvw(jl)/(pvgrat(jl)+pcvw(jl)))*(1./(rswet)))

                  pvdveg(jl,jk,jt)=100./rtot ! cm s-1

				  !if (jt.eq.idt_CO2) pvdveg(jl,jk,jt)=0.
				  
				END IF ! IF (hc(jl) > hcmin

             END DO ! jl=1,klon

           END IF ! IF ( (latmbios(jt))

           ! mz_lg_20061021+ aerosol dry deposition inside the canopy

           IF ( (latmbios(jt)) .AND. lvdaer(jt) ) THEN
             DO jl=1,klon

               !--- possible correction for large humidity close to surface/in canopy
               ! mz_lg_20061021+ currently no effect has been considered
               zalphae(jl)=1._dp  
               zbeta(jl)=1._dp

               zrint=1.e-8*1.E2  ! um -> CM, mz_lg_20061021+ assumed particle radius 
                                 ! of 1e-2 um
                      densaer(jl,jmod)=1.7e3 ! mz_lg_20061021+ assumed density of 1.7 kg cm-3
              
               IF (hc(jl) > hcmin .AND. &
                  (pvgrat(jl)+pcvw(jl)) >= 0._dp) THEN

                  ! calculation of in-canopy turbulent resistance needed to 
                  ! consider soil deposition
                  rahcan(jl)=MAX(10._dp,14.*lai(jl)*hc(jl)/ &  ! mz_lg_20061021+ minimum
                           ustveg(jl))                         ! Rahcan of 10 s m-1

                  ! Cunningham factor
                  cunning=1.+(cl/(zalphae(jl)*zrint*1.e2_dp)**zbeta(jl))*   &
                         (2.514+0.800*EXP(-0.55*(zrint*1.e2_dp)/cl))
                  ! Diffusivity:
                  diffc =(bc*tsurf(jl)*cunning)/(3.*pi*dynvisc*  &
                         (zalphae(jl)*zrint*1.e2_dp)**zbeta(jl))
                  ! Relaxation coefficient
                  relax=(densaer(jl,jmod)*1.e-3*                         &
                      (((zalphae(jl)*zrint)**zbeta(jl))**2.)*            &
                         cunning)/(18.*dynvisc*kappa)
                  ! Sedimentation:
                  sedspeed=(((((zalphae(jl)*zrint)**zbeta(jl))**2.)* &
                       densaer(jl,jmod)*1.e-3*zg*cunning)/ &  ! ESS_lg_20120721+ g to zg; important modification (see above)
                       (18.*dynvisc))
                  ! Calculation of schmidt and stokes number
                  sc      =visc/diffc
                  ! mz_lg_20031014+, modified calculation of the Stokes number
                  !     over vegetated surfaces, see paper by Gallagher et al., 
                  !     JGR 2002. zAL is a characteristic radius for the
                  !     largest collectors comprising the surface
                  st_veg=MAX((relax*(100.*ustveg(jl))**2.)/(zg*zAL),1.E-1_dp) ! ESS_lg_20120721+ g to zg; important modification (see above)
                  ! Brownian diffusion
                  vb_veg   =(1./vkar)*((ustveg(jl)/u_veg(jl,jk))**2)*100.*  &
                             u_veg(jl,jk)*(sc**(-2./3.))
                  ! mz_lg_20031014+, modified calculation of the impaction over
                  !     vegetated surfaces, see paper by Gallagher et al., JGR 2002. 
                  !     We have applied here the parameterization by Slinn [1982] 
                  !     over vegetated surfaces
                  vim_veg   =(1./vkar)*((ustveg(jl)/u_veg(jl,jk))**2)*100.* &
                             u_veg(jl,jk)*(st_veg**2/(1.+st_veg**2))
                  ! mz_lg_20031014+, modified calculation of the surface resistance over 
                  !     vegetated surfaces, see paper by Gallagher et al., JGR 2002. 
                  !     The calculation includes the interception collection efficiency 
                  !     vim and a rebound correction factor R
                  vin_veg  =(1./vkar)*((ustveg(jl)/u_veg(jl,jk))**2)*100.*  &
                            u_veg(jl,jk)*(1./2.)*(zrint/zAS)**2       ! equation 18
                  zrebound =exp(-st_veg**0.5)                         ! equation 19
                  vkd_veg  =zrebound*(vb_veg+vim_veg+vin_veg)
                  vkc_veg     =(100./rahcan(jl))
                  pvdveg(jl,jk,jt)  =1./((1./vkc_veg)+(1./vkd_veg))   ! cm s-1

                  ! mz_lg_20061021-

               END IF ! IF (hc(jl) > hcmin

             END DO ! jl=1,klon

           END IF ! IF ( (latmbios(jt) .AND. lvdaer(jt))

         END DO ! jt=1,ntrac

       END DO ! jk=1,nveglay

    ENDIF

    ! mz_lg_20050719+ inclusion of the actual exchanges calculations

    ! LG- ===second num. solver coupling vert. diff., dry depos. and emiss==
    !     calculation of concentrations as a function of vertical transport
    !     dry deposition and emissions within one equation, so without the 
    !     operator splitting as it has been used in the default version. 

    IF (nveglay >= 2) THEN

       ! LG-  calculation of canopy top flux from gradient between surface 
       !      layer concentration and the average canopy concentration.

       DO jl=1,klon

	      print *,'veg_mlay: timestep',jl

          ! LG-   only calculating a canopy top flux for the areas with a 
          !       canopy height > HCMIN and with a vegetation and wet skin 
          !       fraction> 0, in order to prevent the scheme to calculate fluxes
          !       for the snow covered surfaces with a canopy height > HCMIN 

          IF (hc(jl) > hcmin .AND.   &
             (pvgrat(jl)+pcvw(jl)) > 0._dp) THEN

             ! LG-   determining the reference height of the canopy layers

             ! ESS_lg_20140423+ modified calculation to consider > 2 layers
             DO jk=1,nveglay
               pz(jk+1)=hc(jl)-((jk-1)*(hc(jl)/nveglay)+0.5*(hc(jl)/nveglay))
			   pz_edge(jk)=hc(jl)-(jk-1)*(hc(jl)/nveglay)
			 ENDDO
             pz_edge(nveglay+1)=0.
			 ! LG-    and volume, the volume for the two layers is similar

             grvol_veg(1:nveglay)=((hc(jl)/nveglay)/          &
               (pdp(jl)/(prhoa(jl)*g*1.e3_dp)))*      &
               (pvgrat(jl)+pcvw(jl))*grvol(jl)
             grmass_veg(1:nveglay)=(grvol_veg(1:nveglay)/grvol(jl))*  &
               grmass(jl)

             pdp_veg(:)=(hc(jl)/nveglay)*(prhoa(jl)*g*1.e3_dp)

             ! LG-    and that of the surface layer
             pz(1)=(pdp(jl)/(prhoa(jl)*g*1.e3_dp))/2._dp+    &
                hc(jl)

			 ! LG-    calculating the eddy diffusivity from the Dz between the 
             !        reference height of the first and the second canopy layer. 
             !        There are two approaches which both calculate this resistance 
             !        for the canopy height. Therefore the term HC/DZ has been 
             !        introduced to correct for the difference between the vertical 
             !        extent that has been represented

             dz=pz(1)-pz(2)
             Khmin=(dz**2)/dtime_mixing
             Kh(jl,1)=MAX(Khmin,(dz/rahveg(jl)))

             rahcan(jl)=MAX(1._dp,14.*lai(jl)*hc(jl)/ &
                 ustveg(jl))
             ! ESS_lg_20140423+ modified calculation of Kh for > 2 layers
		     DO jk=1,nveglay-1   		! ESS_km_20130221 added loop 
               dz=pz(jk+1)-pz(jk+2)
               Khmin=(dz**2)/dtime_mixing
               ! ESS_lg_20131125+ calculation from Rahcan (which is a parameterization based on measurements within Maize!)
			   Kh(jl,jk+1)=MAX(Khmin,(dz/(rahcan(jl)/(hc(jl)/dz))))  
               ! ESS_lg_20131125+ estimate of in-canopy Kh based on scaling SL Kh (calc. from Rahveg) with wind speed profile
               Kh(jl,jk+1)=MAX(Khmin,(zrefsl(jl)-pz(jk+1))/rahveg(jl)* &
			       ((u_veg(jl,jk)+u_veg(jl,jk+1))/2)/              & ! average wind speed in canopy at canopy layer edge
			       ((u_sl(jl)+u_veg(jl,1))/2))                       ! estimate of wind speed around the canopy top
			 ENDDO
             ! ESS_lg_20140423-			
			
             DO jt=1,ntrac

                ! mz_lg_20050719+ all the concentrations are recalculated
                ! from mixing ratio to molecules cm-3; kg m-3 -> g cm-3 1e-3
                IF (pxtm1_obs(jl,jt) > 0._dp) THEN ! -9999.999) THEN
                   zxtm1(jl,jt)=(pxtm1(jl,jt)+weight_pxtm1_obs*(pxtm1_obs(jl,jt)-pxtm1(jl,jt))) &
                     *1.e-3_dp*prhoa(jl)/(amd/avo) ! using the observed (surface layer) mixing ratios
                ELSE
                   zxtm1(jl,jt)=pxtm1(jl,jt)*1.e-3_dp*prhoa(jl)/(amd/avo)
                ENDIF
                zxtmveg(:,jt)=pxtmveg(jl,:,jt)*         &
                      1.e-3_dp*prhoa_veg(:)/(amd/avo)

                ! ESS_lg_20130106+ -------initializing the vertical diffusion tendency--------------------
                            
                DO jk=1,nveglay
                  jjk=jk+1                           
                  xtediff(jl,jjk,jt)=zxtmveg(jk,jt)/(1e-3_dp*prhoa_veg(jk))   
                ENDDO
                xtediff(jl,1,jt)=zxtm1(jl,jt)/(1e-3_dp*prhoa(jl)) 
                ! ESS_lg_20130106-

                IF (.NOT.latmbios(jt)) cycle

                ! LG-     setting the soil surface flux to zero

                fluxveg(nveglay+1)=0._dp

                ! ESS_lg_20120721+ calculation of compensation point and resulting emission flux; 
                ! normally done with the dry deposition calculations (see above). However, since canopy 
                ! concentrations are needed for these calculations, which are updated within this 
                ! jl=1,klon loop, this calculation has been moved to this location

                IF (ccomp(jt) > 0._dp) THEN

				  DO jk=1,nveglay
                     zxtmveg(jk,jt)=pxtmveg(jl,jk,jt)*1.e-3_dp*prhoa(jl)/(amd/avo)
                     ! LG-     considering a compensation point 
                     ! ESS_lg_20070803+   calculation of canopy concentration (at leaf level)
                     Ccanopy(jk,jt)=((zxtmveg(jk,jt)/(0.01*rbveg(jk,jt)))+                   & ! Ccanopy in [molec. cm-3]
                       (ccomp(jt)/(0.01*rstomx(jl,jt))))/                                    & ! [molec. cm-2 s-1]
                       (1./(0.01*rbveg(jk,jt))+1./(0.01*rstomx(jl,jt))+1./(0.01*rcut(jt)))     ! [molec. cm-2 s-1 (s-1 cm)-1

                     pvdveg(jl,jk,jt)=0._dp                                         ! setting Vd to 0.
                     ! ESS_lg_20131125+ if Ccanopy < zxtmveg we calculate a negative emission flux; a deposition term
					 ! This term should still be use to calculate properly the dry deposition/emission tendency
					 ! For a negative emission flux, this flux/divided by the zxtmveg gives then inferred pvdveg used 
					 ! to calculate later on the dry deposition tendency
                     emisflux(jl,jk,jt)=emisflux(jl,jk,jt)+                      &  ! adding to already def. flux.
                         ((pvgrat(jl)/(pvgrat(jl)+pcvw(jl)))*                    &
                         (Ccanopy(jk,jt)-zxtmveg(jk,jt))/(0.01*rbveg(jk,jt)))    &  ! molecules cm-2 s-1
                         *1.e4*lad_veglay(jl,jk)*lai(jl)                            ! molecules m-2 s-1 for whole canopy layer

                  ENDDO

                END IF ! IF (ccomp(jt) > 0._dp) THEN
				
				! ESS_lg_20140508+ calculation of O3 stomatal flux to diagnose the role of the stomatal versus 
                !      non-stomatal flux in the total ozone removal

                IF (jt == idt_O3) THEN
                  DO jk=1,nveglay
				    pstomfluxO3(jl,jk)=                                           &
                           -zxtmveg(jk,jt)*1e6*                                   & ! [molec. m-3], negative for deposition flux  
                           (1./(rbveg(jk,jt)+rstomx(jl,jt)+rmes(jt)))                               
                  ENDDO
				ENDIF
 				! ESS_lg_20140508-

                ! ESS_lg_20120721- 

                ! LG-     definition of length of subtimestep used for coupled 
                !         turbulence dry deposition and emission calculations. 
                !         The number of required subtimesteps is a function of the 
                !         timescale of the dry deposition and the emission process 
                !         calculated from the thickness of the canopy layers and the dry
                !         deposition velocity/emission flux. 

                ! LG-     dry deposition timescale

                ts_dd=(hc(jl)/2._dp)/   &
                   MAX(1.e-2_dp,1.e-2_dp*pvdveg(jl,1,jt))

                ! LG-     turbulence timescale

                ts_turb=((hc(jl)/2._dp)**2._dp)/MAX(1.e-10_dp,Kh(jl,1))

                ! LG-     the applied subtimestep is determined from the dry deposition 
                !         timescale, the minimum value of the timestep is determined by 
                !         the value of the parameter PTMST_VEG

                ptmst_sub=MAX(ptmst_veg,MIN(ptmst,0.1_dp*MIN(ts_dd,ts_turb))) 
                nstep_sub=MAX(1,INT(ptmst/ptmst_sub))

                ! LG-     to get exactly the same total emission fluxes, the applied
                !         length of the subtimestep needs to be a connected to the
                !         determined number of subtimesteps

                ptmst_sub=ptmst/nstep_sub

                ! LG-     end

                DO jk=nveglay+1,1,-1
                  jjk=jk-1

                  ! EMIFAC in m2 s cm-3

                  IF (jjk > 0) THEN
                    emifacveg(jjk)=0._dp
                    emifacveg(jjk)=1.e-6_dp*ptmst_sub/(hc(jl)/nveglay) ! ESS_lg_20140730+ nveglay instead of 2_dp                  
                  ENDIF

                  ! mz_lg_20050719+ all the concentrations are recalculated
                  ! from mixing ratio to molecules cm-3

                  IF (JK == nveglay+1) THEN
                    dz1=(hc(jl)/nveglay)
                    dz2=(hc(jl)/nveglay)
                    dz3=0._dp
                    terma(jk)=Kh(jl,jjk)*ptmst_sub/(pz(jk-1)*dz2- &
                      pz(jk)*dz2)
                    termb(jk)=0._dp
                    termc(jk)=1._dp+1.e-2_dp*(pvdveg(jl,jjk,jt)/dz2)*ptmst_sub+ &
                      terma(jk)+termb(jk)
                    em(jk)=emifacveg(jjk)*emisflux(jl,jjk,jt) 
                    zxtmveg_old(jjk)=zxtmveg(jjk,jt)
                  ELSEIF(jk > 1) THEN
                    dz1=(pdp(jl)/(prhoa(jl)*g*1.e3_dp))
                    dz2=(hc(jl)/nveglay)
                    dz3=(hc(jl)/nveglay)
                    terma(jk)=Kh(jl,jjk)*ptmst_sub/(pz(jk-1)*dz2-  &
                      pz(jk)*dz2)
                    termb(jk)=Kh(jl,jjk+1)*ptmst_sub/(pz(jk)*dz2-  &
                      pz(jk+1)*dz2)
                    termc(jk)=1._dp+1.e-2_dp*(pvdveg(jl,jjk,jt)/dz2)*ptmst_sub+ &
                      terma(jk)+termb(jk)
                    em(jk)=emifacveg(jjk)*emisflux(jl,jjk,jt) 
                    zxtmveg_old(jjk)=zxtmveg(jjk,jt)
                  ELSEIF(jk == 1) THEN
                    dz1=0._dp
                    dz2=(pdp(jl)/(prhoa(jl)*g*1.e3_dp))
                    dz3=(hc(jl)/nveglay)
                    terma(jk)=0._dp
                    termb(jk)=Kh(jl,jjk+1)*ptmst_sub/(pz(jk)*dz2- &
                      pz(jk+1)*dz2)
                    termc(jk)=1._dp+terma(jk)+termb(jk)
                    em(jk)=0._dp 
                    zxtm1_old=zxtm1(jl,jt)
				  ENDIF
                  termd(jk)=terma(jk)/termc(jk)
                  terme(jk)=termb(jk)/termc(jk)

                ENDDO

                ! LG-     Start loop for performing emission, dry deposition, and 
                !         turbulence calculations multiple times to deal with 
                !         subtimestep scale process timescales such as for HNO3 and 
                !         O3 dry deposition. 

                dmass_em=0.
                dmass_dd=0.

				! ESS_lg_20140423+ calculation for the original 2-layer model version
                !     ------------------------------------------------------------------
                !     ---> This code can only be used for a 2-layer canopy since the 
                !     implicit solver is restricted to a maximum of three layers, the 
                !     surface layer and the canopy layers <-----
                !     ------------------------------------------------------------------

				IF (nveglay == 2) THEN

	   			  DO it=1,nstep_sub

                    ! LG-      top canopy layer

                    jk=2
                    jjk=1

                    ! LG-      budget calculations for dry deposition

                    ddfrac=(pvdveg(jl,jjk,jt)*ptmst_sub/100._dp)/(hc(jl)/2._dp)
                    dmass_dd=dmass_dd-pxtmveg(jl,jjk,jt)*ddfrac* &
                        grvol_veg(jjk)

                    ! LG-      for emissions within the canopy and the surface layer, 

                    dmass_em=dmass_em+em(jk)*grvol_veg(jjk)
                                   
                    ! LG-      end

                    zxtmveg(jjk,jt)=                                      &
                       ((zxtmveg(jjk,jt)+em(jk))/termc(jk)+               &
                        (termd(jk)*zxtm1(jl,jt))/termc(jk-1)+             &
                        (terme(jk)*zxtmveg(jjk+1,jt))/termc(jk+1))/       &
                        (1._dp-termd(jk)*terme(jk-1)-terme(jk)*termd(jk+1))
				
                    ! ESS_lg_20130106+ calculating the dry deposition tendency which is used 
                    !              with the total tendency to arrive at the vertical diffusion tendency

                    ! ESS_lg_20131125+ calculation of pvdveg from a negative emission (= deposition) flux
				    IF (emisflux(jl,jjk,jt)<0._dp) pvdveg(jl,jjk,jt)=-emisflux(jl,jjk,jt)/zxtmveg(jjk,jt)

                    IF (pvdveg(jl,jjk,jt).GT.0._dp)  &
                      xtedryd(jl,jk,jt)=-1.e-2*pvdveg(jl,jjk,jt)*(zxtmveg(jjk,jt)/  &  
                        (1.e-3_dp*prhoa_veg(jjk)))/(hc(jl)/2._dp)  

				    ! ESS_lg_20131125+ and resetting pvdveg to zero again (for correct calculation of atmosphere-biosphere flux)
				    IF (emisflux(jl,jjk,jt)<0._dp) pvdveg(jl,jjk,jt)=0._dp
					   
				    ! ESS_lg_20130106+ calculating the emission tendency
                    IF (em(jk).gt.0._dp) &
                       xteemis(jl,jk,jt)=em(jk)/(1.e-3_dp*prhoa_veg(jjk))/ptmst_sub 
                    ! ESS_lg_20130106+   
                                           
                    ! LG-      lower canopy layer
 
                    jk=3
                    jjk=2

                    ! LG-      budget calculations for dry deposition

                    ddfrac=(pvdveg(jl,jjk,jt)*ptmst_sub/100._dp)/(hc(jl)/2._dp)
                    dmass_dd=dmass_dd-zxtmveg(jjk,jt)*ddfrac* &
                       grvol_veg(jjk)

                    ! LG-      for emissions within the canopy and the surface layer, 

                    dmass_em=dmass_em+em(jk)*grvol_veg(jjk)

                    ! LG-      end

                    zxtmveg(jjk,jt)=                                 &
                       (zxtmveg(jjk,jt)+em(jk))/termc(jk)+           &
                        termd(jk)*zxtmveg(jjk-1,jt)

				    ! ESS_lg_20130106+ calculating the dry deposition tendency which is used 
                    !              with the total tendency to arrive at the vertical diffusion tendency

                    ! ESS_lg_20131125+ calculation of pvdveg from a negative emission (= deposition) flux
				    IF (emisflux(jl,jjk,jt)<0._dp) pvdveg(jl,jjk,jt)=-emisflux(jl,jjk,jt)/zxtmveg(jjk,jt)

                    IF (pvdveg(jl,jjk,jt).GT.0._dp)  &
                      xtedryd(jl,jk,jt)=-1.e-2*pvdveg(jl,jjk,jt)*(zxtmveg(jjk,jt)/  &  
                         (1.e-3_dp*prhoa_veg(jjk)))/(hc(jl)/2._dp)  

				    ! ESS_lg_20131125+ and resetting pvdveg to zero again (for correct calculation of atmosphere-biosphere flux)
				    IF (emisflux(jl,jjk,jt)<0._dp) pvdveg(jl,jjk,jt)=0._dp
					   
                    ! ESS_lg_20130106+ calculating the emission tendency
                    IF (em(jk).gt.0._dp) &
                      xteemis(jl,jk,jt)=em(jk)/(1.e-3_dp*prhoa_veg(jjk))/ptmst_sub 
                    ! ESS_lg_20130106+   

                    ! LG-      surface layer

                    jk=1
                    jjk=0

                    zxtm1(jl,jt)=                                             &
                      (1.-(pvgrat(jl)+pcvw(jl)))*zxtm1(jl,jt)+                &
                      (pvgrat(jl)+pcvw(jl))*((zxtm1(jl,jt)+em(jk))/termc(jk)+ &
                       terme(jk)*zxtmveg(jjk+1,jt))

					! LG-      end subtimestep loop DO IT=1,NSTEP_SUB
           
                  END DO

				ELSE

  				  ! ESS_lg_20140423+ otherwise the calculations for > 2 layers using Crank-Nicolson
                  ! This implementation is based on the code that was implemented in the multi-layer version
                  ! of the snowpack chemical exchange modelling system (MLSN-CHEM) by Keenan Murray, 2013.
                  ! This code also uses the contribution by advection (wind pumping) which is relevant for
                  ! snow-ice exchange but this has been de-activated in this particular implementation of 
				  ! the code for the multi-layer vegetation model system. 				  

                  DO it=1,nstep_sub
			
			        DO wp_int=1,1,1 ! -1,1,2; ESS_km_20130221+ this loop switches the velocity direction

    	  		      pre_a(:)=0.0
			          trid_a(:)=0.0
			          pre_b(:)=0.0
			          trid_b(:)=0.0
			          trid_c(:)=0.0
			          pre_c(:)=0.0
			          trid_d(:)=0.0
			
			          ! ESS_km_20130225+ This loop steps through depths calculating variables used to define contributions of
			          ! advection, diffusion, emissions, and deposition in the different layers.
			          ! These variables will be used to create a Crank-Nicholson discretization to be solved by a Thomas solver
			          ! Brief variable definitions:
			          ! terma - diffusion from above layer
			          ! termb - diffusion from below layer
			          ! termf - deposition from current layer

  		      	      DO jk=1,nveglay+1

   	      		        jjk=jk-1

			            IF (JK == nveglay+1) THEN
        	        
        	              dz3=0._dp
			              dz1=pz(jk-1)-pz(jk)
			              dz2=pz_edge(jjk)-pz_edge(jjk+1)

			              terma(jk)=Kh(jl,jjk)*ptmst_sub/(dz2*dz1)
			              termb(jk)=0._dp						
                	      termc(jk)=1._dp+1.e-2_dp*(pvdveg(jl,jjk,jt)/dz2)*ptmst_sub+ &
                             terma(jk)+termb(jk)

  	  		              termf(jk)=1.e-2_dp*(pvdveg(jl,jjk,jt)/dz2)*ptmst_sub

			              pre_a(jk)=-terma(jk)
			              pre_b(jk)=(2._dp+termf(jk)+terma(jk)+termb(jk))
			              pre_c(jk)=-termb(jk)

			              trid_c(jk)=pre_c(jk)
 			              trid_b(jk)=pre_b(jk)
   			              trid_a(jk)=pre_a(jk)		

                          em(jk)=emifacveg(jjk)*emisflux(jl,jjk,jt) 
 
        	              zxtmveg_old(jjk)=zxtmveg(jjk,jt)

  			            ELSEIF(jk > 1) THEN

        	              dz1=pz(jk-1)-pz(jk)
        	              dz2=pz_edge(jjk)-pz_edge(jjk+1)
        	              dz3=pz(jk)-pz(jk+1)

			 	          terma(jk)=Kh(jl,jjk)*ptmst_sub/(dz2*dz1)			
 			              termb(jk)=Kh(jl,jjk+1)*ptmst_sub/(dz2*dz3)			
        	              termc(jk)=1._dp+1.e-2_dp*(pvdveg(jl,jjk,jt)/dz2)*ptmst_sub+ &
                            terma(jk)+termb(jk)
			              termf(jk)=1.e-2_dp*(pvdveg(jl,jjk,jt)/dz2)*ptmst_sub

			              pre_a(jk)=-terma(jk)
			              pre_b(jk)=(2._dp+termf(jk)+terma(jk)+termb(jk))
			              pre_c(jk)=-termb(jk)

			              trid_c(jk)=pre_c(jk)
			              trid_b(jk)=pre_b(jk)
			              trid_a(jk)=pre_a(jk)

        	              em(jk)=emifacveg(jjk)*emisflux(jl,jjk,jt) 
        	              zxtmveg_old(jjk)=zxtmveg(jjk,jt)
						
	  		            ELSEIF(jk == 1) THEN
					  
        	              dz1=0._dp
        	              dz2=2.0*abs(pz(jk)-pz_edge(jk))
			              dz3=pz(jk)-pz(jk+1)

        	              terma(jk)=0._dp
			              termb(jk)=Kh(jl,jjk+1)*ptmst_sub/(dz2*dz3)
			              termc(jk)=1._dp+terma(jk)+termb(jk)
		   	              termf(jk)=0.0_dp !1.e-2_dp*(pvdveg(jl,jjk,jt)/dz2)*ptmst_sub

			              pre_a(jk)=-terma(jk)
			              pre_b(jk)=(2._dp+termf(jk)+terma(jk)+termb(jk))
			              pre_c(jk)=-termb(jk)

			              trid_c(jk)=pre_c(jk)
			              trid_b(jk)=pre_b(jk)
			              trid_a(jk)=pre_a(jk)

        	              em(jk)=0._dp 
        	              zxtm1_old=zxtm1(jl,jt)

	  	    	        ENDIF
		                termd(jk)=terma(jk)/termc(jk)
		                terme(jk)=termb(jk)/termc(jk)
			          ENDDO   ! ESS_km_20130220  end of do loop through depths

                      ! ESS_lg_20140423+ modified budget calculations for dry deposition 
					  DO jk=2,nveglay+1
					    jjk=jk-1
			            ddfrac=(pvdveg(jl,jjk,jt)*ptmst_sub/100._dp)/(hc(jl)/nveglay)
                        dmass_dd=dmass_dd-pxtmveg(jl,jjk,jt)*ddfrac* &
                          grvol_veg(jjk)

                        ! LG-      for emissions within the canopy and the surface layer, 

                        dmass_em=dmass_em+em(jk)*grvol_veg(jjk)
                                                       
                        ! ESS_lg_20130106+ calculating the dry deposition tendency which is used 
                        !              with the total tendency to arrive at the vertical diffusion tendency

                        ! ESS_lg_20131125+ calculation of pvdveg from a negative emission (= deposition) flux
				        IF (emisflux(jl,jjk,jt)<0._dp) pvdveg(jl,jjk,jt)=-emisflux(jl,jjk,jt)/zxtmveg(jjk,jt)

                        IF (pvdveg(jl,jjk,jt).GT.0._dp)  &
                          xtedryd(jl,jk,jt)=-1.e-2*pvdveg(jl,jjk,jt)*(zxtmveg(jjk,jt)/  &  
                            (1.e-3_dp*prhoa_veg(jjk)))/(hc(jl)/nveglay)  

				        ! ESS_lg_20131125+ and resetting pvdveg to zero again (for correct calculation of atmosphere-biosphere flux)
				        IF (emisflux(jl,jjk,jt)<0._dp) pvdveg(jl,jjk,jt)=0._dp
					   
				        ! ESS_lg_20130106+ calculating the emission tendency
                        IF (em(jk).gt.0._dp) &
                          xteemis(jl,jk,jt)=em(jk)/(1.e-3_dp*prhoa_veg(jjk))/ptmst_sub 
                        ! ESS_lg_20130106+   
 
					  ENDDO

					  ! ESS_lg_20140423-
                    
			          DO jk=1,nveglay+1
			            jjk=jk-1

			            IF (jk == nveglay+1) THEN
			              trid_d(jk)=-pre_a(jk)*zxtmveg(jjk-1,jt)+    &
			                (4._dp-pre_b(jk))*zxtmveg(jjk,jt)+2._dp*em(jk)
			            ELSEIF (jk > 2) THEN
			              trid_d(jk)=-pre_a(jk)*zxtmveg(jjk-1,jt)+    &
				            (4._dp-pre_b(jk))*zxtmveg(jjk,jt)         &
				            -pre_c(jk)*zxtmveg(jjk+1,jt)+2._dp*em(jk)
			            ELSEIF (jk == 2) THEN
			              trid_d(jk)=-pre_a(jk)*zxtm1(jl,jt)+         &
			                (4._dp-pre_b(jk))*zxtmveg(jjk,jt)         &
				            -pre_c(jk)*zxtmveg(jjk+1,jt)+2._dp*em(jk)
			            ELSEIF (jk == 1) THEN
			  	          trid_d(jk)= (4._dp-pre_b(jk))*zxtm1(jl,jt)  &
				            -pre_c(jk)*zxtmveg(jjk+1,jt)+2._dp*em(jk)
			            ENDIF

			          ENDDO

		  	          CALL solve_tridiag(trid_a,trid_b,trid_c,trid_d,trid_x,nveglay+1)

  			          DO jk=2,nveglay+1
  			            zxtmveg(jk-1,jt)=trid_x(jk)
 			          ENDDO
			          zxtm1(jl,jt)=trid_x(1)

				    ENDDO  ! wp_int loop

		            ! LG-      end subtimestep loop DO IT=1,NSTEP_SUB
				              
                  END DO

				ENDIF
				
                ! LG-     calculation of the fluxes at the top of the canopy
                !         layer from the dC/dt, the deposition/emission flux and the 
                !         flux between the layer and the underlying layer
  
                ! mz_lg_20050721+ fluxveg is in molecules cm-2 s-1 !
                DO jk=nveglay,1,-1
                  fluxveg(jk)=                                      & 
                     -((1.e2_dp*(hc(jl)/nveglay)*                   &  ! ESS_lg_20140730+ nveglay instead of 2._dp
                       (zxtmveg(jk,jt)-zxtmveg_old(jk)))/ptmst+     &
                        pvdveg(jl,jk,jt)*zxtmveg(jk,jt)-            &
                        1.e-4_dp*emisflux(jl,jk,jt)-                &
                        fluxveg(jk+1))

                  ! ESS_lg_20150617+ activate to check the calculated canopy flux profiles
                  ! For Radon there should be a very small/no flux divergence with its limited sinks 				  
				  !IF (jt.eq.idt_rad) THEN
				  !  print *,'jk,iradon,flux_veg and emission flux: ',jk,fluxveg(jk),emisflux(jl,jk,jt)
				  !ENDIF
				  ! ESS_lg_20150617-

				END DO

                ! diagnostics
                fluxctop(jl,jt)=1.e-11*fluxveg(1)

                ! mz_lg_20050721+ patmbiosflux is in molecules m-2 s-1 !
                patmbiosflux(jl,jt)=1.e4_dp*fluxveg(1) ! mz_lg_20050719+ assigning flux
     
                ! mz_lg_20050721+ Canopy Reduction Factor, which is the ratio
                !     of the canopy top flux to the emission flux

                ! ESS_lg_20120718+
                ! mz_lg_20050725+ determining the canopy intergrated 
                ! emission flux
                DO jk=1,nveglay
                  emisflux_tot(jl,jt)=emisflux_tot(jl,jt)+emisflux(jl,jk,jt)
                ENDDO
                ! ESS_lg_20120718-

                IF (emisflux_tot(jl,jt) > 0._dp) THEN
                   pcrf(jl,jt)=patmbiosflux(jl,jt)/emisflux_tot(jl,jt)
                ELSE
                   pcrf(jl,jt)=1._dp ! setting it default to 1. (also for zero
                                     ! nocturnal emissions of VOC's
                ENDIF
       
                ! LG-     calculation of total surface flux considering all the 
                !         land cover fractions, this flux must be defined in 
                !         molecules m-2 s-1

                psurfflux(jl,jt)=psurfflux(jl,jt)+           &
                  (pvgrat(jl)+pcvw(jl))*patmbiosflux(jl,jt)

                ! diagnostics
                dc(jl,jt)=(zxtmveg(1,jt)-zxtmveg_old(1))/(1.e-3_dp*prhoa(jl))

             END DO ! DO jt=1,ntrac

             ! ESS_lg_20120722+ call of chemistry routine(s)
             IF (l_xtsurf_veg_mlay_chem) THEN

               ! ESS_lg_20130106+ -------prepare chemistry tendency--------------------
               DO jt=1,ntrac                                     
                 DO jk=1,nveglay
                   jjk=jk+1
                   xtechem(jl,jjk,jt)=zxtmveg(jk,jt)/(1.e-3_dp*prhoa_veg(jk))
                 ENDDO
                 xtechem(jl,1,jt)=zxtm1(jl,jt)/(1.e-3_dp*prhoa(jl))
               ENDDO
               ! ESS_lg_20130106-          

               CALL emdep_xtsurf_chem(jl, Nstep, ntrac,             &
                       ptmst, prhoa, tsurf, press, qm1, rj, rj_veg, &
                       zxtm1, zxtmveg, OHreact) ! ESS_lg_20130113+                

               ! ESS_lg_20130106+ -------calculating the chemistry tendency--------------------
               DO jt=1,ntrac                                     
                 DO jk=1,nveglay
                   jjk=jk+1
                   xtechem(jl,jjk,jt)=(zxtmveg(jk,jt)/(1.e-3_dp*prhoa_veg(jk))- &
                      xtechem(jl,jjk,jt))/ptmst
                 ENDDO
                 xtechem(jl,1,jt)=(zxtm1(jl,jt)/(1.e-3_dp*prhoa(jl))- &
                    xtechem(jl,1,jt))/ptmst
               ENDDO
               ! ESS_lg_20130106-          

             ENDIF
             ! ESS_lg_20120722-

             ! ESS_lg_20130106+ -------calculating the vertical diffusion tendency--------------------
             DO jt=1,ntrac                                     
               DO jk=1,nveglay
                 jjk=jk+1                           
                 xtediff(jl,jjk,jt)=(zxtmveg(jk,jt)/(1e-3_dp*prhoa_veg(jk))-  &           
                   xtediff(jl,jjk,jt))/ptmst-                                 &
                   xtedryd(jl,jjk,jt)-xteemis(jl,jjk,jt)-xtechem(jl,jjk,jt)   
			   ENDDO
               xtediff(jl,1,jt)=(zxtm1(jl,jt)/(1e-3_dp*prhoa(jl))-            &           
                      xtediff(jl,1,jt))/ptmst-                                &
                      xtedryd(jl,1,jt)-xteemis(jl,1,jt)-xtechem(jl,1,jt)  
 			 ENDDO
             ! ESS_lg_20130106-
                         
             ! ESS_lg_20120718+ updating pxtmveg and pxtm1
             ! The concentrations are recalculated from molecules cm-3 to mixing ratio;
             ! g cm-3 -> kg m-3 1e3 
             DO jk=1,nveglay
               IF (jl.LT.klon) & ! ESS_lg_20120924+ modified to overcome boundary problems
                 pxtmveg(jl+1,jk,:)=zxtmveg(jk,:)/  &    ! ESS_lg_20120721+ to secure updating the concentrations; jl+1
                     (1.e-3_dp*prhoa_veg(jk)/(amd/avo))
			 ENDDO
                
		     IF (jl.LT.klon) & ! ESS_lg_20120924+ modified to overcome boundary problems
               pxtm1(jl+1,:)=zxtm1(jl,:)/  &             ! ESS_lg_20120721+ to secure updating the concentrations; jl+1
                  (1.e-3_dp*prhoa(jl)/(amd/avo))

             ! ESS_lg_20120718-

             ! mz_lg_20050721+ NOx Canopy Reduction Factor: special case
             ! This will be modified whenever also here the tracer family 
             ! concept will be introduced
             IF ( (emisflux_tot(jl,idt_NO) > 0._dp) .AND. &
                  (idt_NO > 0 .AND. idt_NO2 > 0) ) &
               pcrf(jl,idt_NO2)=   &
                  (patmbiosflux(jl,idt_NO)+patmbiosflux(jl,idt_NO2))/ &
                   emisflux_tot(jl,idt_NO)

               ! ESS_lg_20130124+ added the updated NOx_emflux to also include the leaf emissions
               ! associated with the role of the compensation point in addition to the potential source
               ! of NOx in the form of NO2 due to nitrate photolysis
               nox_emflux(jl,1)=emisflux(jl,1,idt_NO2)
               nox_emflux(jl,2)=emisflux(jl,2,idt_NO2)
               ! ESS_lg_20130124-

             END IF ! LG-    endif  IF (HC(JL).GT.HCMIN.AND.....

       END DO ! DO jl=1,klon

    END IF ! LG-  endif IF (nveglay == 2)

    ! mz_lg_20050719+ deallocate 3D arrays
    DEALLOCATE(emisflux)

    ! mz_lg_20050719+ deallocate 2D arrays
    DEALLOCATE(rstomx)
    DEALLOCATE(rleaf)
    DEALLOCATE(rleafw)
    DEALLOCATE(rbveg)
    DEALLOCATE(lad_veglay)
    DEALLOCATE(emisflux_tot)
    DEALLOCATE(fluxctop)
    DEALLOCATE(zxtm1)
    DEALLOCATE(zxtmveg)
    DEALLOCATE(dc)

    ! mz_lg_20050719+ deallocate 1D arrays
    DEALLOCATE(zxtmveg_old)
    DEALLOCATE(emifacveg)
    DEALLOCATE(grvol_veg)
    DEALLOCATE(grmass_veg)
    DEALLOCATE(pdp_veg)
    DEALLOCATE(prhoa_veg)

    DEALLOCATE(fluxveg)
    DEALLOCATE(terma)
    DEALLOCATE(termb)
    DEALLOCATE(termc)
    DEALLOCATE(termd)
    DEALLOCATE(terme)
    DEALLOCATE(termf)
    DEALLOCATE(em)
    DEALLOCATE(pz)
    DEALLOCATE(rahcan)
    DEALLOCATE(rmesophyll)
    ! mz_lg_20050719-

	! ESS_lg_20140423+ 
	DEALLOCATE(pz_edge)
	DEALLOCATE(trid_a)
	DEALLOCATE(trid_b)
	DEALLOCATE(trid_b2)
	DEALLOCATE(trid_c)
	DEALLOCATE(trid_d)
	DEALLOCATE(trid_x)
	DEALLOCATE(pre_a)
	DEALLOCATE(pre_b)	
	DEALLOCATE(pre_c)
	! ESS_lg_20140423-
	
    RETURN

  END SUBROUTINE emdep_xtsurf_veg_mlay

  ! mz_lg_20050718-

  !===========================================================================

  SUBROUTINE emdep_xtsurf_calc_rs( &
    lo_derived, lvd_bigl, trname, moleweight, dryreac, henry, &
    idt_SO2, diff, diffrb, rsoil, rwater, rws, rsnow, rmes, rcut) ! ESS_lg_20150619+
    !-----------------------------------------------------------------------------
    ! subroutine emdep_xtsurf_calc_rs, to calculate the values of the uptake resistances
    ! required to calculate the trace gas dry deposition velocity. this routine
    ! is based on an approach by wesely, 1989, in which the uptake resistances of
    ! trace gases, for which the dry deposition velocities have not been observed,
    ! are estimated based on the henry coefficient and a reactivity coefficient and
    ! the uptake resistances of so2 and o3, of the "big leaf" dry deposition scheme
    ! by ganzeveld and j. lelieveld j. geophys. res., 100, 20,999-21,012,1995,
    ! ganzeveld et al.,j. geophys. res., 103, 5679-5694, 1998 and ganzeveld et al,
    ! submitted to j. geophys. res., 2001. for more information of the wesely
    ! approach see atmospheric environment vol 23, no 6, 1293-1304.
    !
    ! the program needs as input data the molecular mass of the defined trace
    ! gases, the henry coefficient [M atm-1] and an estimated reactivity
    ! coefficient which has 3 distinct values: 0 for non-reactive species,
    ! (e.g, so2, acetaldehyde), 0.1 for moderately reactive species (e.g., pan),
    ! and 1 for reactive species (e.g., o3, hno3). these values are defined in
    ! the module mo_*_request_tracer
    !-----------------------------------------------------------------------------

    ! mz_lg_20040503+ 
    ! Interface:
    ! ----------
    ! input 
    ! lo_derived : true whenever SO2 and O3 are defined as tracers
    ! lvd_bigl   : true for the tracer when the required parameters are defined
    ! trname     : tracer name
    ! moleweight : molecular weight 
    ! dryreac    : reactivity coefficient [0:non-react., 0.1:semi react., 1:react.]        
    ! henry      : henry coefficient [mol atm-1]
    ! 
    ! output
    ! 

    IMPLICIT NONE 

    ! I/O
    LOGICAL,  INTENT(in)   :: lo_derived, lvd_bigl(:)
    CHARACTER(LEN=20), INTENT(in) :: trname(:)
    REAL(dp), INTENT(in)   :: moleweight(:), dryreac(:), henry(:)
    REAL(dp), INTENT(inout):: diff(:), diffrb(:), &
            rsoil(:), rwater(:), rws(:), rsnow(:), rmes(:), rcut(:)

    ! declarations

    INTEGER :: jt, ntrac, idt_SO2 ! ESS_lg_2015061

    REAL(dp):: diffrb_so2, rsoil_so2, rwater_so2, rws_so2,  &
               rsnow_so2,  rmes_so2,  rcut_so2,   diff_so2, &
               diffrb_o3,  rsoil_o3,  rwater_o3,  rws_o3,   &
               rsnow_o3,   rmes_o3,   rcut_o3,    diff_o3

    ! INITIALIZATION
    ! mz_lg_20040423+ added
    ntrac=SIZE(lvd_bigl)

    ! mz_lg_20050721+ initializing
    diff(:)=0._dp 
    diffrb(:)=0._dp
    rsoil(:)=0._dp
    rwater(:)=0._dp
    rws(:)=0._dp
    rsnow(:)=0._dp
    rmes(:)=0._dp
    rcut(:)=0._dp
    ! mz_lg_20050721-

    ! mz_lg_20020115 definition of the different resistances for the 6
    !     species of the original dry deposition scheme:
    !     ------------o3, hno3, no, no2, so2 and so4--------------

    ! mz_lg_20020115 definition of terms which correct for the diffusivity
    !     for the computation of the boundary layer resistance, (v/dx)**2/3
    !     v = 0.189 cm2 s-1 (heat) and dx has been calculated according to:
    !     dx=0.212*sqrt(mh2o/mx) with 0.212 being the diffusivity of water
    !     vapour and m is the molar mass of the component.

    !--- attribute specific parameters in the following order:
    !
    !    - soil resistance
    !    - sea water resistance, which is generally similar to the wet skin
    !    - wet skin reservoir resistance
    !    - snow resistance
    !    - mesophyll resistance
    !    - cuticle resistance
    !    - diffusivity coefficient, to correct stomatal resistance for
    !      differences in diffusivity between water vapour and the
    !      specific trace gas (sqrt(molmass trace gas)/sqrt(molmass h2o))

    !--- Values for SO2 and O3 are required in any case to estimate the
    !    resistances of the other species considered in the deposition scheme.
    !    => Define hard-linked resistances  for SO2 and O3 even if they
    !       are not present as tracer.

    ! mz_lg_20030811+, SO2
    diffrb_so2=1.6_dp
    rsoil_so2=250._dp
    rwater_so2=1._dp
    rws_so2=100._dp
    rsnow_so2=1._dp
    rmes_so2=1._dp
    rcut_so2=1.e5_dp
    diff_so2=1.9_dp

    ! mz_lg_20030811+, O3
    diffrb_o3=1.2_dp
    rsoil_o3=400._dp
    rwater_o3=2000._dp
    rws_o3=2000._dp
    ! ESS_lg_20150811+ modified for sensitivity analyses
    ! rws_o3=750._dp

    rsnow_o3=2000._dp

    ! mz_lg_20030111+ modified for sensitivity analyses
    ! rsnow_o3=750._dp

    rmes_o3=1._dp
    rcut_o3=1.e5_dp
    diff_o3=1.6_dp

    ! mz_lg_20020115 the resistances are only calculated for species that
    !     have not been included in the original dry deposition scheme.

    DO jt=1, ntrac

     ! mz_lg_20030811+, modified to deal with the fact that there might
     !     multiple tracers starting with the same basename but with
     !     with different subnames, e.g., SO2 (e.g., SO2, SO2_GM7)

     SELECT CASE(TRIM(trname(jt))) ! mz_pj_20040330

     CASE('SO2')

       diffrb(jt)=diffrb_so2
       rsoil(jt)=rsoil_so2
       rwater(jt)=rwater_so2
       rws(jt)=rws_so2
       rsnow(jt)=rsnow_so2
       rmes(jt)=rmes_so2
       rcut(jt)=rcut_so2
       diff(jt)=diff_so2

     CASE('O3')

       diffrb(jt)=diffrb_o3
       rsoil(jt)=rsoil_o3
       rwater(jt)=rwater_o3
       rws(jt)=rws_o3
       rsnow(jt)=rsnow_o3
       rmes(jt)=rmes_o3
       rcut(jt)=rcut_o3
       diff(jt)=diff_o3

     CASE('SO4')

       diffrb(jt)=1.8_dp
       rsoil(jt)=1.e5_dp
       rwater(jt)=1.e5_dp
       rws(jt)=1.e5_dp
       rsnow(jt)=1.e5_dp
       rmes(jt)=1.e5_dp
       rcut(jt)=1.e5_dp
       diff(jt)=2.7_dp

     CASE('HNO3')

       diffrb(jt)=1.4_dp
       rsoil(jt)=1._dp
       rwater(jt)=1._dp
       rws(jt)=1._dp
       rsnow(jt)=1._dp
       rmes(jt)=1._dp
       rcut(jt)=1._dp
       diff(jt)=1.9_dp

     CASE('NO')

       diffrb(jt)=1.1_dp
       rsoil(jt)=1.e5_dp
       rwater(jt)=1.e5_dp
       rws(jt)=1.e5_dp
       rsnow(jt)=1.e5_dp
       rmes(jt)=500._dp
       rcut(jt)=1.e5_dp
       diff(jt)=1.3_dp

     CASE('NO2')

       diffrb(jt)=1.2_dp
       rsoil(jt)=600._dp
       rwater(jt)=1.e5_dp
       rws(jt)=1.e5_dp

	   ! ESS_lg_20150823+
       !rws(jt)=750._dp
	   
       rsnow(jt)=1.e5_dp
       rmes(jt)=1._dp
       rcut(jt)=1.e5_dp
       diff(jt)=1.6_dp

	 ! ESS_lg_20130503+  
	 CASE('CO2')

       diffrb(jt)=1.2_dp
       rsoil(jt)=1.e5_dp
       rwater(jt)=1.e5_dp
       rws(jt)=1.e5_dp
       rsnow(jt)=1.e5_dp
       rmes(jt)=1._dp
       rcut(jt)=1.e5_dp
       diff(jt)=1.6_dp

	 ! ESS_lg_20130503-  
	   
     CASE DEFAULT

       ! mz_lg_20030411+ modified: only O3 and SO2 need to be defined for estimating
       !     the various resistances, not the other species of the default scheme

       ! mz_lg_20040504+ lo_derived is true whenever O3 is present. SO2 is now
       !     longer required because of the above defined specific resistances for
       !     SO2 and O3

       IF (lo_derived .AND. (lvd_bigl(jt))) THEN

          ! calculation of term which is used to correct the stomatal resistance
          ! for differences in the diffusitivy (see also equation 4).

          diff(jt)=SQRT(moleweight(jt)/18.)     ! mz_pj_20040330

          ! calculation of the term to correct for differences in diffusivity
          ! between water vapor and the trace gas. it is calculated from:
          ! diff bl=(v/dx)**2/3, with v=0.189 sm-1 and dx= dh2o/sqrt(mh2o/mx),
          ! with dh2o=0.212

          diffrb(jt)=(0.189/(0.212/diff(jt)))**(2./3.)

          ! calculation of rmx, the mesophyll resistance

          rmes(jt)=1./(henry(jt)/3000.+100.*dryreac(jt))

          ! calculation of rlux, the cuticular resistance, equation 7 of wesely's
          ! paper

          ! ESS_lg_20150619+ modified calculation of Rcut scaling with 1/henry(idt_SO2) instead of
		  ! scaling with 1e-5 (which resembled the ...)

          !rcut(jt)=1./(henry(jt)/1e5+dryreac(jt))*rcut_o3

          rcut(jt)=1./(henry(jt)/henry(idt_SO2)+dryreac(jt))*rcut_o3

          ! calculation of rgsx, the soil resistance, equation 9 of wesely's
          ! paper
          
		  !rsoil(jt)=1./(henry(jt)/(1e5*rsoil_so2)+ &
          !     dryreac(jt)/rsoil_o3)    
          rsoil(jt)=1./(henry(jt)/(henry(idt_SO2)*rsoil_so2)+ &
               dryreac(jt)/rsoil_o3)          

          ! the snow resistance is similar as the soil resistance
          ! mz_lg_20050628+ modified to assure calculations of rsnow for
          ! for species with a Henry and reactivity of that of ozone similar
          ! to the rsnowO3 of 2000 s m-1 (based on feedback from A. Pozzer)

          !rsnow(jt)=1./(henry(jt)/(1e5*rsnow_so2)+ &
          !     dryreac(jt)/rsnow_o3)
          rsnow(jt)=1./(henry(jt)/(henry(idt_SO2)*rsnow_so2)+ &
               dryreac(jt)/rsnow_o3)

          ! calculation of rlux-wet, the wet skin resistance, equation 14 of
          ! wesely's paper

          rws(jt)=1./(1./(3.*rws_so2)+1.e-7*henry(jt)+ &
               dryreac(jt)/rws_o3)          

          ! calculation of sea uptake resistance, using equation 9 of wesely's
          ! paper

          !rwater(jt)=1./(henry(jt)/(1.e5*rwater_so2)+ &
          !     dryreac(jt)/rwater_o3)
          rwater(jt)=1./(henry(jt)/(henry(idt_SO2)*rwater_so2)+ &
               dryreac(jt)/rwater_o3)			   
			   
       ENDIF

     END SELECT

    ENDDO

  END SUBROUTINE emdep_xtsurf_calc_rs

  !=============================================================================

  SUBROUTINE emdep_xtsurf_calcprof( nstep, & ! ESS_lg_20120722+
    nveglay_hr, netrad, cossza_2d, lai, lad, rbvd, rvdsl, rvd, fsl) ! mz_lg_20050721+
    ! ---------------------------------------------------------------
    !     Calculation of distribution of PAR (diffuse and direct)
    !     within the canopy as a function of the solar zenith angle
    !     the LAI and the radiation above the canopy. The code
    !     is taken from the DDIM model (Dry deposition Inferential
    !     Model) and developed by Norman and Weiss, 1985. The
    !     windprofile is also calculated acc. to Cionco.
    !
    !     Laurens Ganzeveld, 1998, modified for implementation
    !     in echam5, October, 2001
    ! --------------------------------------------------------------------
    ! mz_sw_20040206 slightly rewritten (also removed all subroutine input parameters)

    ! mz_lg_20040503+ 
    ! Interface:
    ! ----------
    ! input 
    ! nstep     : # of timesteps ! ESS_lg_20120722+
    ! nveglay_hr: no. of canopy layers; high resolution ! mz_lg_20050721+
    ! netrad      : net surface radiation [W M-2]
    ! cossza_2d : cosine of the zenith angle [0-1]
    ! lai       : Leaf Area Index [m2 m-2]
    ! lad       : leaf area density profile
    !
    ! output
    ! rbvd      : direct beam irradiance [W m-2]
    ! rvdsl     : diffusive irradiance in surface layer [W m-2] ! mz_lg_20050721+
    ! rvd       : diffusive irradiance in canopy [W m-2]
    ! fsl       : fraction of sunlit leaves [0-1]

    IMPLICIT NONE 

    ! I/O
    INTEGER,  INTENT(in)  :: nstep ! ESS_lg_20120722+
    INTEGER,  INTENT(in)  :: nveglay_hr ! mz_lg_20050721+

    REAL(dp), INTENT(in)  :: netrad(:), cossza_2d(:), lai(:), lad(:,:)
    REAL(dp), INTENT(out) :: rbvd(:), rvdsl(:), rvd(:,:), fsl(:,:) 
    ! mz_lg_20040503-

    INTEGER :: klon, i, ii, iday, jl

    REAL(dp) :: zrvd(nveglay_hr),  zfsl(nveglay_hr), &
                mlai(nveglay_hr,2),mlaitot(nveglay_hr)

    REAL(dp) :: rg,zen,parbeam,zrbvd,zrvdsl,twopi

    !
    !  *********************************************************************
    !  *                                                                   *
    !  *        K E Y         V A R I A B L E S                            *
    !  *                                                                   *
    !  *                                                                   *
    !  *   RG = GLOBAL RADIATION (WATTS/M^2)                               *
    !  *                                                                   *
    !  *   LAI = LEAF AREA INDEX                                           *
    !  *                                                                   *
    !  *   JDAY = JULIAN DAY                    LAIW = WINTER LAI          *
    !  *                                                                   *
    !  *   H! = CANOPY HEIGHT (M)                                          *
    !  *                                                                   *
    !  *   FSL = FRACTION OF SUNLIT LEAVES    CSNL = COSINE OF LEAF NORMAL *
    !  *                                                                   *
    !  *   RVD = DIFFUSE VISIBLE RADIATION    RBVD = VISIBLE BEAM RADIATION*
    !  *                                                                   *
    !  *   PAR = PHOTOSYNTHETICALLY ACTIVE RADIATION (.4 - .7 MICRONS)     *
    !  *                                                                   *
    !  *   FVIS = FRACTION OF THE GLOBAL RADIATION THAT IS PAR             *
    !  *                                                                   *
    !  *   PCNTLF = PERCENTAGE OF MAX LEAF AREA FOR NON-CONIFERS           *
    !  *                                                                   *
    !  *                                                                   *
    !  *   ** UNLESS SPECIFIED, UNITS ARE SI                               *
    !  *                                                                   *
    !  *********************************************************************
    !

    ! mz_sw_20040206+
    ! INITIALIZATION
    klon=nstep ! SIZE(netrad) ! ESS_lg_20120722+
    ! mz_sw_20040206-

    twopi = 2.*pi

    ! added
    rbvd=0._dp
    rvdsl=0._dp
    rvd=0._dp
    fsl=0._dp

    ! LG- loop longitude

    DO jl=1,klon
       rg=netrad(jl)
       zen=ACOS(cossza_2d(jl))
           
       ! LG-  Be carefull with LAD vertical profile, see also the routine
       !      vegetation.f, LAD(1) is the LAD value for the top vegetation
       !      layer (in constrast to the single column model)

       ! LG-  normally the LAD profile for one ecosystem type (mlai(:,1) is
       !      considered but this can be extended to more vegetation types

       DO ii=1,nveglay_hr
          mlai(ii,1) = lad(jl,ii)*lai(jl)
          mlaitot(ii) = mlai(ii,1)
       ENDDO

       iday = 1

       IF((rg.LT.10).OR.(zen.GT.(twopi/4))) then
          iday = 0
          zrbvd = 0._dp
          zrvdsl= 0._dp
          DO ii=1,nveglay_hr
             zrvd(ii) = 0._dp
             zfsl(ii) = 0._dp
          ENDDO
       ENDIF

       ! LG- call of routine in which the radiation profiles in the
       !     canopy are being calculated

       IF(iday.GT.0)  &
            CALL emdep_xtsurf_canrad2(nveglay_hr,mlaitot,zen,rg, &
                                zfsl,zrvd,zrbvd,zrvdsl,parbeam)
                                                                
       DO ii=1,nveglay_hr
          ! mz_lg_20050725+ note the change in the index with the value
          ! for i=1 reflecting the value in the highest canopy layer! 
          i=nveglay_hr+1-ii
          fsl(jl,i)= zfsl(ii)
          rvd(jl,i)= zrvd(ii)
       ENDDO

       rbvd(jl)=zrbvd
       rvdsl(jl)=zrvdsl

       ! LG- end longitudinal loop

    ENDDO

  END SUBROUTINE emdep_xtsurf_calcprof

  !=============================================================

  SUBROUTINE emdep_xtsurf_canrad2(n,pai,zen,rg,fsl,rvd,rbvd,rvdsl,parbeam)
    !
    !    ***************************************************
    !    *  S U B R O U T I N E    C A N R A D 2           *
    !    *                                                 *
    !    *  THIS SUBROUTINE COMPUTES THE VERTICAL PROFILE  *
    !    *  OF VISIBLE RADIATION (PAR), BOTH BEAM AND      *
    !    *  DIFFUSE COMPONENTS FOR A MULTILAYER CANOPY     *
    !    *  NLEV IS THE TOPLAYER ABOVE THE CANOPY!!!       *
    !    ***************************************************

    !  _________________________________________________________
    !
    !    It turned out that there is an error in this code
    !    RBVD is calculated from the RDV by dividing through cos(ZEN),
    !    which is incorrect!!!!. It likely has to do with the definition
    !    of the zenith angle in DDIM
    !  _________________________________________________________
    !

    IMPLICIT NONE

    INTEGER :: N

    REAL(dp) :: PAI(N),RG,ZEN,FVD,RBV(N),PI, &
         ALPHA,BETA,KXM,KX(9),IB(N),MIB(N),AV(N), &
         TV,PV,X1(11),KXD(11),ID(N),MID(N),Z1,Z2,Z3,FSL(N), &
         RVD(N),RVU(N),ORVU,CKVU,CPAI

    INTEGER :: I,J,K,M,ITER

    REAL(dp) :: OT,RDVIS,RFVIS,WA,RDIR,RFIR,RVT,RIRT,FVIS, &
         FIR,RATIO,FVB,PARBEAM,PARDIFF,THETA,RBVD,RVDSL

    !
    !    ******************************************************************
    !    *                                                              *
    !    *   THE NEXT FEW STATEMENTS DETERMINE THE BEAM AND DIFFUSE       *
    !    *   COMPONENTS OF VISIBLE RADIATION (FOR DETAILS SEE WEISS AND   *
    !    *   NORMAN, 1985, AGRICULTURAL METEOROLOGY.......                *
    !    *                                                              *
    !    ******************************************************************
    !

    ZEN=MIN(ZEN,1.56_dp)
    OT=35./(1224.*COS(ZEN)**2.+1)**.5
    RDVIS=600.*2.7182**(-.185*OT)*COS(ZEN)
    RFVIS=0.4*(600-RDVIS)*COS(ZEN)
    WA= 1320*.077*(2.*OT)**0.3
    RDIR=(720.*2.7182**(-0.06*OT)-WA)*COS(ZEN)
    RFIR=0.60*(720.-WA-RDIR)*COS(ZEN)
    RVT=RDVIS+RFVIS
    RIRT=RDIR+RFIR
    FVIS=RVT/(RIRT+RVT)
    FIR=RIRT/(RIRT+RVT)
    RATIO=RG/(RVT+RIRT)
    IF(RATIO.GE..9) THEN
       RATIO=0.899_dp
       RG=RVT+RIRT
    ENDIF

    FVB=RDVIS/RVT*(1.-((.9-RATIO)/7.)**0.67)
    FVD=1.-FVB

    ! LG-  calculation of RVDSL, the diffusive radiation in the surface
    !      layer (see subroutine RJVEG and CALCRAD)

    RVDSL = MAX(1.E-5_dp,RG*FVIS*FVD)

    ! LG-  end

    PARBEAM = MAX(1.E-5_dp,RG*FVIS*FVB*COS(ZEN))
    PARDIFF = MAX(1.E-5_dp,RG*FVIS*FVD)

    !
    !     **************************
    !     *  INITIALIZE CONSTANTS  *
    !     **************************
    !

    PI=3.1415926_dp

    !
    !     *******************************************************
    !     *  SET FRACTION OF TOTAL THAT IS BEAM RADIATION AND   *
    !     *  THEN SEPARATE INTO INTO IR AND VISIBLE COMPONENTS  *
    !     *******************************************************
    !
    !

    RBV(N)=RG*FVIS*(1.-FVD)

    !
    !     *****************************************
    !     * COMPUTATION OF EXTINCTION COEFFICIENT *
    !     *****************************************
    !

    ALPHA=PI/2. - ZEN
    BETA=PI/36.
    KXM = 0.5/SIN(ALPHA)
    THETA = PI/36
    DO I=1,9
       KX(I)=0._dp
       THETA = THETA+PI/18
    ENDDO

    !
    !     ****************************************************
    !     *  COMPUTE PROBABILTIY FUNCTION FOR PENETRATION    *
    !     *  OF THE BEAM COMPONENT , FROM NORMAN, 1979       *
    !     *  MODIFICATION OF THE AERIAL ENVIRONMENT OF CROPS *
    !     ****************************************************
    !

    DO I=1,N-1
       K=N - I

       ! LG-       originally it is PAI(K+1) but this has been replaced by
       !           PAI(I+1)

       IB(K+1)=2.7182818**(-KXM*PAI(I+1))
       MIB(K+1)=1.0 - IB(K+1)
       RBV(K)=RBV(K+1)*IB(K+1)
    ENDDO

    !
    !     ************************************
    !     *  SET SOIL VISIBLE AND IR ALBEDO  *
    !     ************************************
    !

    AV(1)=0.10_dp
    TV=0.01_dp
    PV=0.08_dp
    ALPHA=PI/20.

    !
    !     **************************************************
    !     *   LOOP FOR COMPUTING OFTENLY USED FACTORS      *
    !     **************************************************
    !

    DO I=2,11
       X1(I)=SIN(ALPHA)*COS(ALPHA)
       KXD(I)=.5/SIN(ALPHA)
       ALPHA=ALPHA+PI/20.
    ENDDO

    !
    !     *******************************************************
    !     * LOOP FOR COMPUTING ID (DIFFUSE RADIATION PENETRATION*
    !     * FUNCTION) AND A (R UP/R DOWN)                       *
    !     *******************************************************
    !

    DO I=2,N
       ID(I)=0.0_dp
       Z1=0.0_dp
       K=1
       DO J=1,5

          ! LG-          originally it is PAI(I) but this has been replaced by
          !              PAI(N+1-I)

          Z2=2.718282**(-PAI(N+1-I)*KXD(K+1))*X1(K+1)

          ! LG-          originally it is PAI(I) but this has been replaced by
          !              PAI(N+1-I)

          Z3=2.718282**(-PAI(N+1-I)*KXD(K+2))*X1(K+2)
          ID(I)=ID(I)+PI*(Z1+4.*Z2+Z3)/(20.*3.)
          Z1=Z3
          K=K+2
       ENDDO
       ID(I)=ID(I)*2.0
       IF(ID(I).GT.1) ID(I) = 1
       MID(I)=1.0 - ID(I)
       AV(I)=AV(I-1)*(TV*MID(I)+ID(I))*(TV*MID(I)+ID(I))/  &
            (1.0 - AV(I-1)*PV*MID(I)) + PV*MID(I)
    ENDDO

    !
    !     *************************************************
    !     *  INITIALIZE DOWNWARD DIFFUSE COMPONENTS OF    *
    !     *  VISIBLE RADIATION THE CANOPY                 *
    !     *************************************************
    !

    RVD(N)=RG*FVD*(1.0 - FIR)

    !
    !    ****************************************************
    !    * COMPUTE DOWNWARD DIFFUSE RADIATION AT EACH LEVEL *
    !    ****************************************************
    !
    DO I=1,N-1
       K=N-I
       RVD(K)=RVD(K+1)*(TV*MID(K+1)+ID(K+1))/(1.0-AV(K)*PV &
            *MID(K+1))
    ENDDO

    !
    !     *************************************************
    !     *  COMPUTE UPWARD DIFFUSE FLUXES USING THE SOIL *
    !     *  ALBEDO                                     *
    !     *************************************************
    !

    RVU(1)=AV(1)*(RVD(1)+RBV(1))

    DO I=2,N
       RVU(I)=RVU(I-1)*(TV*MID(I)+ID(I))*AV(I)/(AV(I)-PV*MID(I))
    ENDDO

    !
    !     *************************************
    !     *  START MAIN LOOP FOR ITERATIONS   *
    !     *************************************
    !
    !

    ORVU=0.0_dp
    ITER=0
    DO M=1,100
       ITER=ITER + 1
       DO I=1,N-1
          K = N - I
          RVD(K)=RVD(K+1)*(TV*MID(K+1)+ID(K+1))+RVU(K)*PV*MID(K+1) &
               +RBV(K+1)*MIB(K+1)*TV
       ENDDO

       RVU(1)=AV(1)*(RVD(1)+RBV(1))

       !
       !     ************************************************
       !     *   COMPUTE THE UPWARD DIFFUSE FLUXES FOR      *
       !     *   THE VISIBLE WAVELENGTHS                    *
       !     ************************************************
       !

       DO I=2,N
          RVU(I)=RVU(I-1)*(TV*MID(I)+ID(I))+RVD(I)*PV*MID(I)+ &
               RBV(I)*MIB(I)*PV
       ENDDO

       CKVU=ABS(RVU(1)-ORVU)

       !
       !     *****************************************************
       !     *  CHECK FOR CONVERGENCE OF VALUES                  *
       !     *  CRITERIA FOR CONVERGENCE: DIFFERENCE BETWEEN     *
       !     *  THE OLD AND NEW VALUES AT THE TENTH/FOURTH LEVEL *
       !     *  CHANGES BY NO MORE THAN 2 WATTS/METER SQ.        *
       !     *****************************************************
       !
       !

       IF((CKVU.LE.0.01)) GO TO 100
       ORVU=RVU(1)
    ENDDO

100 CPAI = 0.0_dp

    DO I=1,N
       K = N+1-I

       ! LG-    in original code PAI(K) with top layer=N, replaced by PAI(I)

       CPAI=CPAI+PAI(I)
       FSL(K)=2.7183**(-KXM*CPAI)
    ENDDO
    RBVD=RBV(N)*COS(ZEN)

    RETURN

  END SUBROUTINE emdep_xtsurf_canrad2

  ! ESS_lg_20120721++ implementation of subroutine to calculate/assign compensation points
  !=======================================================================================

  SUBROUTINE emdep_xtsurf_ccomp(nstep,prhoa,vgrat,tsurf,pxtmveg_CO2,idt_CO2, &
                                idt_NO2,idt_NH3,ccomp)

    ! ----------------------------------------------------------
    ! Calculation/definition of the stomatal compenstion
    ! point for use in the dry deposition calculations  
    ! ----------------------------------------------------------
    ! Interface:
    ! ----------
    ! input 
    ! nstep     : # of timesteps
    ! prhoa     : air density in kg m-3
    ! vgrat     : vegetation fraction
    ! tsurf     : surface temperature
	! pxtmveg_CO2 : crown layer CO2 concentration
    ! idt_CO2   : tracer index CO2
    ! idt_NO2   : tracer index NO2
    ! idt_NH3   : tracer index NH3
    !
    ! output
    ! ccomp     : tracer compensation points

    IMPLICIT NONE 

    ! I/O
    INTEGER,  INTENT(in)  :: nstep, idt_CO2, idt_NO2, idt_NH3
    REAL(dp), INTENT(in)  :: prhoa(:),vgrat(:),tsurf(:),pxtmveg_CO2(:)
    REAL(dp), INTENT(out) :: ccomp(:)

    ! local parameters
    INTEGER :: jl,klon

    klon=nstep ! SIZE(prhoa)

    DO jl=1,klon

      !ccomp(idt_CO2)=365*1e-3*prhoa(jl)/(amd*1e6/avo) ! 350 ppmv to molucules cm-3 ! alternative 0.99*pxtmveg_CO2(jl)   

	  ccomp(idt_NO2)=0.5*1e-3*prhoa(jl)/(amd*1e9/avo) ! 0.5 ppbv to molucules cm-3

      ! LG- compensation point of NH3, 10-2001 it is set to a constant value.
      !     However, there is more information in the paper by Asman et al., 1998
      !     (submitted that time, LGP475), showing that the stomatal compensation
      !     point Xs can be calculated as Xs=f(T)[NH4+]/[H+], with f(T) expressing 
      !     the temperature dependence (leaf temperature) and [NH4+] and [H+] being
      !     the apoplastic concentrations, [NH4+] is very sensitive to the leaf   
      !     N status and the external N supply. Typical concentrations, given in
      !     mM range between 0.04 and 2.3 mM, dependent on the external NH4+/NH3 
      !     concentrations. Apoplastic pH (and thus -log[H+]) ranges between 
      !     5.5 and 6.5. For f(T), maybe that papers by Farquhar show some specific
      !     functions, but also the results in Figure 6 of the Asman paper can be
      !     used to define f(T)

      ccomp(idt_NH3)=0.5*1e-3*prhoa(jl)/(amd*1.e9/avo) ! 0.5 ppbv
      ! ESS_lg_20070806+ see excercise NH3 compensation model, Cc model approach 
      ccomp(idt_NH3)=0.
      IF (vgrat(jl) > 0.) THEN ! ESS_lg_20070801+ included the vegfrac
         ccomp(idt_NH3)=0.05*                                       & ! scalings factor, typical concentration should be 1 ug m-3
            ((2749644./tsurf(jl))*EXP(-10378./tsurf(jl))*1e9*2300.) & ! ESS_lg_20070706+ ug m-3
              *(1e-6*1e-6*avo/17.)                                    ! molecules cm-3     
      ENDIF

    ENDDO

  END SUBROUTINE emdep_xtsurf_ccomp
  ! ESS_lg_20120721-

  ! ESS_lg_20120721++ implementation of subroutine to calculate/assign compensation points
  !=======================================================================================

  SUBROUTINE emdep_xtsurf_chem( istep, nstep, ntrac,          &  
                ptmst, prhoa, temp, press, qm1, rj, rj_veg,   &
                pxtm1, pxtmveg, OHreact) ! ESS_lg_20130113+ 

    ! -------------------------------------------------------------------------------
    ! Calculation of chemistry based on the CBM4 isoprene chemistry scheme of the SCM 
    ! -------------------------------------------------------------------------------
    ! Interface:
    ! ----------
    ! input 
    ! istep     : index of actual timestep 
    ! nstep     : # of timesteps
    ! ptmst     : lenght of the timestep
    ! ntrac     : number of tracers
    ! prhoa     : air density in kg m-3
    ! temp      : temperature
    ! press     : pressure
    ! qm1       : moisture in g H2O g-1 air
    ! rj        : surface layer photolysis [s-1]
    ! rj_veg    : canopy photolysis [s-1]
    ! pxtm1     : surface layer concentrations
    ! pxtmveg   : canopy layer concentrations
    ! OHreact   : OH reactivity [s-1]

    use messy_emdep_mem

    IMPLICIT NONE 

    ! I/O
    INTEGER,  INTENT(in)  :: istep, nstep, ntrac ! ESS_lg_20120722+
    REAL(dp), INTENT(in)  :: ptmst, prhoa(:),temp(:),press(:), qm1(:), rj(:,:), rj_veg(:,:,:)
    REAL(dp), INTENT(inout) :: pxtm1(:,:), pxtmveg(:,:)
    REAL(dp), INTENT(out) :: OHreact(:,:) ! ESS_lg_20130113+

    ! local parameters
    INTEGER :: jl,klon
    INTEGER, PARAMETER :: NLON=1, MAE=15, NBINREAC=82
    REAL,    PARAMETER :: ZMAIR=28.970
    LOGICAL, PARAMETER :: LSULFCHEM=.TRUE.,LNH3CHEM=.TRUE., &
	                      l_chem_canopy=.true.  ! MAQ_lg_20160923+ adding a switch off chemistry inside the canopy

    INTEGER JK,JT,ITER,NSTOP
    REAL PM(NLON,nveglay+1,NTRAC), &
   &     PZ(NLON,nveglay+1),PT(NLON,nveglay+1), &
   &     PQ(NLON,nveglay+1),PDP(NLON,nveglay+1),PP(NLON,nveglay+1), &
   &     PPHOTCHEM(NLON,nveglay+1,MAE)

    REAL ZPM(NLON,nveglay+1,NTRAC), &
   &     ZPRHOA(NLON,nveglay+1),ZPT(NLON,nveglay+1),ZPQ(NLON,nveglay+1), &
   &     ZPP(NLON,nveglay+1),ZPPHOTCHEM(NLON,nveglay+1,NTRAC)

    REAL RJNO2(NLON),RJHNO3(NLON),RJO3D(NLON),RJH2O2(NLON),            &
   &     RJMEPE(NLON),RJBCH2O(NLON),RJACH2O(NLON),RJN2O5(NLON),        &
   &     RJANO3(NLON),RJBNO3(NLON),RJHNO4(NLON),RJHONO(NLON),          &
   &     RNOO3(NLON),RHO2NO(NLON),RMO2NO(NLON),RNO2OH(NLON),           &
   &     ROHHNO3(NLON),RNO2O3(NLON),RNONO3(NLON),RNO2NO3(NLON),        &
   &     RN2O5(NLON),RHNO4OH(NLON),RHNO4M(NLON),RNO2HO2(NLON),         &
   &     RODM(NLON),RH2OOD(NLON),RMOO2(NLON),RO3HO2(NLON),             &
   &     RCOOH(NLON),RO3OH(NLON),RHPOH(NLON),RFRMOH(NLON),             &
   &     RCH4OH(NLON),ROHPCAT(NLON),ROHPFRM(NLON),RMO2HO2(NLON),       &
   &     RHO2OH(NLON),RHO2HO2(NLON),RN2O5AQ(NLON),ROHSO2(NLON),        &
   &     ROHDMS(NLON),RNO3DMS(NLON)

    REAL AIR(NLON),O2(NLON),H2O(NLON),                                 &
   &     CH40(NLON),CO0(NLON),HNO30(NLON),H2O20(NLON),CH3O2H0(NLON),   &
   &     ZNO0(NLON),ZNO20(NLON),ZNO30(NLON),ZN2O50(NLON),HNO40(NLON),  &
   &     OH0(NLON),HO20(NLON),O30(NLON),OD0(NLON),CH3O0(NLON),         &
   &     CH3O20(NLON),CH2O0(NLON),ODDN0(NLON),                         &
   &     DMS0(NLON),SO20(NLON),SAER0(NLON),RADON0(NLON),               &
   &     CH4(NLON),CO(NLON),HNO3(NLON),H2O2(NLON),CH3O2H(NLON),        &
   &     ZNO(NLON),ZNO2(NLON),ZNO3(NLON),ZN2O5(NLON),HNO4(NLON),       &
   &     OH(NLON),HO2(NLON),O3(NLON),OD(NLON),CH3O(NLON),              &
   &     CH3O2(NLON),CH2O(NLON),                                       &
   &     DMS(NLON),SO2(NLON),SAER(NLON),RADON(NLON)                    

    REAL RJPAN(NLON),RJALD2(NLON),RJACET(NLON),RJMGLY(NLON),           &
   &     RJMEK(NLON),RJNITR(NLON),RH2OH(NLON),ROHOH(NLON),             &
   &     RFRMNO3(NLON),RMO2MO2(NLON),RALD2OH(NLON),RALD2NO3(NLON),     &
   &     RC23NO(NLON),RC23NO2(NLON),EQPAN,RPAN(NLON),RC23C23(NLON),    &
   &     RC23HO2(NLON),RC23MO2(NLON),RETHOH(NLON),RETHO3(NLON),        &
   &     RPAROH(NLON),RRORA(NLON),RRORB(NLON),RRORNO2(NLON),           &
   &     RRXPAR(NLON),ROLEOH(NLON),ROLEO3(NLON),RMGLYOH(NLON),         &
   &     RISOPOH(NLON),RISOPO3(NLON),RISOPNO2(NLON),RISPDOH(NLON),     &
   &     RISPDO3(NLON),RISPNO3(NLON),RMTHCOH(NLON),RMTHCO3(NLON),      &
   &     RMTHCNO3(NLON),RMVKOH(NLON),RMVKO3(NLON),RMC23NO(NLON),       &
   &     RMC23HO2(NLON),RMC23NO2(NLON),RMPAN(NLON),RMPANO3(NLON),      &
   &     RMPANOH(NLON),RMEKOH(NLON),RACETOH(NLON),                     &
   &     RXO2NO(NLON),ROLENO3(NLON),RISOPNO3(NLON),RISPDNO3(NLON),     &
   &     RXO2XO2(NLON),RXO2HO2(NLON),RXO2NNO(NLON),RXO2N(NLON),        &
   &     RXO2NHO2(NLON),RXO2NXO2(NLON),RBXO2NNO(NLON),RBXO2N(NLON),    &
   &     RBXO2NHO2(NLON),RBXO2NXO2(NLON),RBXO2NXO2N(NLON),             &
   &     RISONTROH(NLON),RHCOOHOH(NLON),RCH3CO2HOH(NLON),              &
   &     RC23MO2A(NLON),RC23MO2B(NLON),ROHNH3(NLON),RNONH2(NLON),      &
   &     RNO2NH2(NLON),RHO2NH2(NLON),RO2NH2(NLON),RO3NH2(NLON),        &
   &     RSO4NH3(NLON),RMATERPOH(NLON),RMBTERPOH(NLON),RSQTERPOH(NLON),&
   &     RMATERPO3(NLON),RMBTERPO3(NLON),RSQTERPO3(NLON),              &
   &     RMATERPNO3(NLON),RMBTERPNO3(NLON),RSQTERPNO3(NLON),           &
   &     RNOOH(NLON),RNO2H2O(NLON),                                    &
   &     RCH2OHO2H(NLON),RRCHOHO2H(NLON),                              &
   &     RSQTERP2BOH(NLON),RSQTERP2BO3(NLON),                          & ! mz_lg_20060330+
   &     RSQTERP2BNO3(NLON),RSQTERP1BOH(NLON),RSQTERP1BO3(NLON),       &
   &     RSQTERP1BNO3(NLON),                                           &
   &     RMTTERPOH(NLON),RMTTERPO3(NLON),RMTTERPNO3(NLON),             & ! mz_lg_20060403+
   &     RISPDHO2(NLON),                                               & ! ESS_lg_20080711+
   &     RHEXANEOH(NLON),RBUTADIENEOH(NLON),RTMBENZENEOH(NLON),        & ! ESS_lg_20100628+ selection of anthropogenic alkenes
   &     RHEXANEO3(NLON),RBUTADIENEO3(NLON),RTMBENZENEO3(NLON),        &
   &     RHEXANENO3(NLON),RBUTADIENENO3(NLON),RTMBENZENENO3(NLON)

    REAL ALD20(NLON),PAR0(NLON),OLE0(NLON),ETH0(NLON),                 &
   &     PAN0(NLON),ACET0(NLON),ISOP0(NLON),MGLY0(NLON),ISOPRD0(NLON), &
   &     METHAC0(NLON),MVK0(NLON),MEK0(NLON),MPAN0(NLON),NITR0(NLON),  &
   &     C2O30(NLON),XO20(NLON),ROR0(NLON),XO2N0(NLON),RXPAR0(NLON),   &
   &     BXO2N0(NLON),MC3O30(NLON),ALD2(NLON),                         &
   &     PAR(NLON),OLE(NLON),ETH(NLON),PAN(NLON),ACET(NLON),           &
   &     ISOP(NLON),MGLY(NLON),ISOPRD(NLON),METHAC(NLON),MVK(NLON),    &
   &     MEK(NLON),MPAN(NLON),NITR(NLON),C2O3(NLON),XO2(NLON),         &
   &     ROR(NLON),XO2N(NLON),RXPAR(NLON),BXO2N(NLON),MC3O3(NLON),     &
   &     ISONTR0(NLON),ISONTR(NLON),HCOOH0(NLON),HCOOH(NLON),          &
   &     CH3CO2H0(NLON),CH3CO2H(NLON),NH20(NLON),NH2(NLON),            &
   &     NH30(NLON),NH3(NLON),NH40(NLON),NH4(NLON),ACID0(NLON),        &
   &     ACID(NLON),MATERP0(NLON),MBTERP0(NLON),SQTERP0(NLON),         &
   &     HONO0(NLON),CH2OHO2H0(NLON),RCHOHO2H0(NLON),                  &
   &     MATERP(NLON),MBTERP(NLON),SQTERP(NLON),HONO(NLON),            &
   &     CH2OHO2H(NLON),RCHOHO2H(NLON),                                &
   &     SQTERP2B0(NLON),SQTERP1B0(NLON),                              & ! mz_lg_20060330+
   &     SQTERP2B(NLON),SQTERP1B(NLON),                                &
   &     MTTERP0(NLON),MTTERP(NLON),                                   &
   &     HEXANE0(NLON),BUTADIENE0(NLON),TMBENZENE0(NLON),              & ! ESS_lg_20100628+
   &     HEXANE(NLON),BUTADIENE(NLON),TMBENZENE(NLON)

    REAL PCH3O2,XLCH3O2,PC2O3,XLC2O3,PXO2,XLXO2,PROR,XLROR,PRXPAR,     &
   &     XLRXPAR,PXO2N,XLXO2N,PBXO2N,XLBXO2N,POLE,XLOLE,PMGLY,XLMGLY,  &
   &     PISOP,XLISOP,PISOPRD,XLISOPRD,PMETHAC,XLMETHAC,PMVK,XLMVK,    &
   &     PMEK,XLMEK,PMC3O3,XLMC3O3,PCH3O2H,XLCH3O2H,PPAN,XLPAN,        &
   &     PMPAN,XLMPAN,PETH,XLETH,PACET,XLACET,PNTR,XLNTR,              &
   &     PALD2,XLALD2,PPAR,XLPAR,                                      &
   &     PISONTR,XLISONTR,PHCOOH,XLHCOOH,PCH3CO2H,XLCH3CO2H,           &
   &     PNH2,XLNH2,XLNH3,PMSA,PMATERP,XLMATERP,PMBTERP,XLMBTERP,      &
   &     PSQTERP,XLSQTERP,PHONO,XLHONO,PCH2OHO2H,XLCH2OHO2H,           &
   &     PRCHOHO2H,XLRCHOHO2H,PSQTERP2B,XLSQTERP2B,PSQTERP1B,XLSQTERP1B, & ! mz_lg_20060330+
   &     PMTTERP,XLMTTERP,                                             &   ! mz_lg_20060403+
   &     PHEXANE,XLHEXANE,PBUTADIENE,XLBUTADIENE,PTMBENZENE,XLTMBENZENE    ! ESS_lg_20100628+ 

    REAL ZFARR,ZF3BOD2,RX1,RX2,ER,ZTREC,ZF3BOD,CT,C0,XP,XL,DLT,        &
   &     ZMH2O,PKL,DT,DT2,DGHNO3,DGAIR,GAMICE,GAMWAT,ZRHOA,ZN2,        &
   &     ZT3REC,RX3,RX4,EQN2O5,EQHNO4,P1,R12,R21,XL1,P2,               &
   &     XL2,P3,XL3,X1,X2,X3,C1,C2,C3,Y2,XJT,R21T,R12T,R12TC,          &
   &     R21TC,XJTC,ACUB,BCUB,CCUB,CUBDET,DNO2,R57,R56,                &
   &     R65,R75,P5,XL5,R66,X5,P6,XL6,X6,C6,XL7,C7,Y1,R89,             &
   &     P8,XL8,X4,C5,R98,XL9,R1011,R1012,C10,R1211,R1112,             &
   &     P11,X11,C11,C12,C1112,XLOD,PHNO3,XLHNO3,PH2O2,XLH2O2,         &
   &     PCH2O,XLCH2O,PCO,XLSO2,PSAER,XLDMS,ZFAC,                      &
   &     RTERPO3                                             ! ESS_lg_20080714+

    REAL RN2O5L(NLON)

    ! LG- extra declarations 

    REAL POHHO2(NLON,nveglay+1),PHO2OH(NLON,nveglay+1), &
         H2OLTR,IS(NBINREAC,nveglay+1),OLDNOY,NEWNOY,DODDN, &
         RSP_LTIME,FXO2_ISOP

    ! LG- declaration of functions which are used in the EBI scheme

    ZFARR(RX1,ER,ZTREC)=RX1*EXP(ER*ZTREC)
    ZF3BOD(RX1,RX2)=RX1/(1+RX1/RX2)*0.6**  &
        (1./(1+ALOG10(RX1/RX2)**2))
    ZF3BOD2(RX1,RX2)=RX1/(1+RX1/RX2)*0.3** &
        (1./(1+ALOG10(RX1/RX2)**2))

    ! further initialization
    DT=PTMST
    DT2=DT**2
    ZMH2O=18.

    ! start the chemistry calculations

    DO 99 JK=1,nveglay+1
    DO 99 JL=1,NLON

      ! LG- assigning the tracer concentrations 

      DO JT=1,NTRAC
        IF (JK.EQ.1) THEN 
          ZPM(JL,JK,JT)=MAX(1.E-30,pxtm1(istep,JT)) ! avoiding numerical problems
        ELSE
          ZPM(JL,JK,JT)=MAX(1.E-30,pxtmveg(JK-1,JT))
        ENDIF
      ENDDO

      ! ESS_lg_20140814+ initialization of intensity of segregation term, this didn't provide a problem 
	  ! with the default 2-layer version but resulted in NaN when the number of layers was made flexible
	  DO JT=1,NBINREAC
        IS(JT,JK)=0.
	  ENDDO
	  ! ESS_lg_20140814-

      ! LG-  assigning some other parameters, e.g. temperature etc.
      ZPRHOA(JL,JK)=1e-3*prhoa(JL)
      ZPQ(JL,JK)=qm1(JL)
      ZPT(JL,JK)=temp(JL)
      ZPP(JL,JK)=press(JL)

      ! LG-  assigning the photodissociation rates 
      DO JT=1,NTRAC
        IF (JK.EQ.1) THEN
           ZPPHOTCHEM(JL,JK,JT)=rj(istep,JT)
        ELSE
           ZPPHOTCHEM(JL,JK,JT)=rj_veg(istep,JK-1,JT)
        ENDIF
      ENDDO

   99 CONTINUE

      DO 101 JK=1,nveglay+1

      DO 1102 JL=1,NLON
        ZRHOA=ZPRHOA(JL,JK)
        AIR(JL)=ZRHOA*AVO/ZMAIR
        O2(JL)=0.209476*AIR(JL)
        ZN2=0.78084*AIR(JL)
        H2O(JL)=ZPQ(JL,JK)*ZRHOA*AVO/ZMH2O  ! molecules cm-3
        H2O(JL)=MAX(H2O(JL),0.1)

        ! LG-   added the calculation of the available water in liters!

        H2OLTR=ZPQ(JL,JK)*ZRHOA*1.E-3 ! zpq in g H2O g-1 air and zrhoa 
             ! in g air cm-3 gives g H2O cm-3. 1 g = 1e-3 kg = 1 liter
 1102 CONTINUE

      ! initialization of photolysis rates
      ! rjach2o: ch2o -> co+h2
      ! rjbch2o: ch2o -> co+2ho2
      ! rjano3:  no3 -> no2+o3 
      ! rjbno3:  no3 -> no

      DO 1103 JL=1,NLON
        IF (ZPPHOTCHEM(JL,JK,idt_NO2).GT.1e-20) THEN
          RJH2O2(JL)=ZPPHOTCHEM(JL,JK,idt_H2O2)
          RJHNO3(JL)=ZPPHOTCHEM(JL,JK,idt_HNO3)
          RJNO2(JL)=ZPPHOTCHEM(JL,JK,idt_NO2)
          RJN2O5(JL)=ZPPHOTCHEM(JL,JK,idt_N2O5)
          RJACH2O(JL)=ZPPHOTCHEM(JL,JK,idt_CH2O)
          RJBCH2O(JL)=ZPPHOTCHEM(JL,JK,idt_CH2O)
          RJO3D(JL)=ZPPHOTCHEM(JL,JK,idt_O3)
          RJMEPE(JL)=ZPPHOTCHEM(JL,JK,idt_CH3O2H)
          RJHNO4(JL)=ZPPHOTCHEM(JL,JK,idt_HNO4)
          RJANO3(JL)=ZPPHOTCHEM(JL,JK,idt_NO3)
          RJBNO3(JL)=ZPPHOTCHEM(JL,JK,idt_NO3)
          RJPAN(JL)=ZPPHOTCHEM(JL,JK,idt_PAN)
          RJALD2(JL)=ZPPHOTCHEM(JL,JK,idt_ALD2)
          RJACET(JL)=ZPPHOTCHEM(JL,JK,idt_ACET)
          RJMGLY(JL)=ZPPHOTCHEM(JL,JK,idt_MGLY)
        ELSE
          RJNO2(JL)=0.
          RJHNO3(JL)=0.
          RJO3D(JL)=0.
          RJH2O2(JL)=0.
          RJMEPE(JL)=0.
          RJBCH2O(JL)=0.
          RJACH2O(JL)=0.
          RJN2O5(JL)=0.
          RJANO3(JL)=0.
          RJBNO3(JL)=0.
          RJHNO4(JL)=0.
          RJPAN(JL)=0.
          RJALD2(JL)=0.
          RJACET(JL)=0.
          RJMGLY(JL)=0.
        ENDIF
        RJMEK(JL)=3.6E-4*RJNO2(JL)
        RJNITR(JL)=4.8*RJHNO3(JL)

        ! LG-   12-2003, added the photodissociation of HONO (see personal
        !       communication by Yvonne Trebs:
        !       J(HONO)=0.189*J(NO2)+8.433*10^(-2)*J(NO2)^2  
        !       Kraus and Hofzumahaus (1998)

        RJHONO(JL)=0.189*RJNO2(JL)+8.433E-2*RJNO2(JL)**2
 1103 CONTINUE

      ! initializing the reaction rates
      DO 1104 JL=1,NLON
        ZTREC=1./ZPT(JL,JK)
        ZT3REC=300./ZPT(JL,JK)

        ! ESS_lg_20080722+ modified to include the intensity of segregation

        RNOO3(JL)=(1.+IS(1,JK))*ZFARR(2.E-12,-1400.,ZTREC)
        RHO2NO(JL)=(1.+IS(2,JK))*ZFARR(3.5E-12,250.,ZTREC)
        RMO2NO(JL)=(1.+IS(3,JK))*ZFARR(3.0E-12,280.,ZTREC)
           RX1=2.5E-30*ZT3REC**4.4*AIR(JL)
           RX2=1.6E-11*ZT3REC**1.7
        RNO2OH(JL)=(1.+IS(4,JK))*ZF3BOD(RX1,RX2)
           RX1=ZFARR(7.2E-15,785.,ZTREC)
           RX3=ZFARR(1.9E-33,725.,ZTREC)
           RX4=ZFARR(4.1E-16,1440.,ZTREC)
        ROHHNO3(JL)=(1.+IS(5,JK))*RX1+RX3*AIR(JL)/(1.+RX3*AIR(JL)/RX4)
        RNO2O3(JL)=(1.+IS(6,JK))*ZFARR(1.2E-13,-2450.,ZTREC)
        RNONO3(JL)=(1.+IS(7,JK))*ZFARR(1.5E-11,170.,ZTREC)
           RX1=2.2E-30*ZT3REC**3.9*AIR(JL)
           RX2=1.5E-12*ZT3REC**0.7
        RNO2NO3(JL)=(1.+IS(8,JK))*ZF3BOD(RX1,RX2)
          EQN2O5=4.E-27*EXP(10930.*ZTREC)
        RN2O5(JL)=RNO2NO3(JL)/EQN2O5
        RHNO4OH(JL)=(1.+IS(9,JK))*ZFARR(1.3E-12,380.,ZTREC)
           RX1=1.8E-31*ZT3REC**3.2*AIR(JL)
           RX2=4.7E-12*ZT3REC**1.4
        RNO2HO2(JL)=(1.+IS(10,JK))*ZF3BOD(RX1,RX2)
          EQHNO4=2.1E-27*EXP(10900.*ZTREC)
        RHNO4M(JL)=RNO2HO2(JL)/EQHNO4

        RODM(JL)=0.2094*ZFARR(3.2E-11,70.,ZTREC) &
             +0.7808*ZFARR(1.8E-11,110.,ZTREC)
        RH2OOD(JL)=(1.+IS(11,JK))*2.2E-10
        RH2OH(JL)=0.5E-6*ZFARR(5.5E-12,-2000.,ZTREC)
        ROHOH(JL)=ZFARR(4.3E-12,-240.,ZTREC)
        RO3HO2(JL)=(1.+IS(12,JK))*ZFARR(1.1E-14,-500.,ZTREC)
        RCOOH(JL)=(1.+IS(13,JK))*1.5E-13*(1.+0.6*ZPP(JL,JK)/101325.)
        RO3OH(JL)=(1.+IS(14,JK))*ZFARR(1.6E-12,-940.,ZTREC)
        RHPOH(JL)=(1.+IS(15,JK))*ZFARR(2.9E-12,-160.,ZTREC)
        RFRMOH(JL)=(1.+IS(16,JK))*1.0E-11
        RFRMNO3(JL)=(1.+IS(24,JK))*ZFARR(3.4E-13,-1900.,ZTREC)  ! ESS_lg_20080722+ note the number
        RCH4OH(JL)=(1.+IS(17,JK))*ZFARR(2.45E-12,-1775.,ZTREC)
        ROHPCAT(JL)=(1.+IS(18,JK))*0.7*ZFARR(3.8E-12,200.,ZTREC)
        ROHPFRM(JL)=(1.+IS(19,JK))*0.3*ZFARR(3.8E-12,200.,ZTREC)
        RMO2HO2(JL)=(1.+IS(20,JK))*ZFARR(3.8E-13,800.,ZTREC)
        RMO2MO2(JL)=(1.+IS(21,JK))*ZFARR(9.1E-14,416.,ZTREC)
        RHO2OH(JL)=(1.+IS(22,JK))*ZFARR(4.8E-11,250.,ZTREC)
        RHO2HO2(JL)=(1.+IS(23,JK))*(ZFARR(2.3E-13,600.,ZTREC)+ &
     &               ZFARR(1.7E-33*AIR(JL),1000.,ZTREC))* &
     &              (1.+H2O(JL)*ZFARR(1.4E-21,2200.,ZTREC))

        ! aerosol N2O5 -> 2 HNO3; set to zero
        RN2O5AQ(JL)=0. 
        RN2O5L(JL)=0. 

        ! O1D steady state
        RJO3D(JL)=RJO3D(JL)*H2O(JL)*RH2OOD(JL)/ &
     &     (H2O(JL)*RH2OOD(JL)+AIR(JL)*RODM(JL))

        ! organic reaction rates;
        ! taken from Duncan and Chameides (1998); Stockwell ea (1997);
        ! Houweling ea (1998)
        RALD2OH(JL)=(1.+IS(25,JK))*ZFARR(7.E-12,250.,ZTREC)
        RALD2NO3(JL)=(1.+IS(26,JK))*ZFARR(1.4E-12,-1900.,ZTREC)
        RC23NO(JL)=(1.+IS(27,JK))*ZFARR(3.5E-11,-180.,ZTREC)
          RX1=9.7E-29*ZT3REC**5.6*AIR(JL)
          RX2=9.3E-12*ZT3REC**1.5
        RC23NO2(JL)=(1.+IS(28,JK))*ZF3BOD2(RX1,RX2)
          EQPAN=8.62E-29*EXP(13954.*ZTREC)
        RPAN(JL)=RC23NO2(JL)/EQPAN
        RC23C23(JL)=(1.+IS(29,JK))*ZFARR(2.5E-12,500.,ZTREC)
        RC23HO2(JL)=(1.+IS(30,JK))*6.5E-12
        RC23MO2(JL)=(1.+IS(52,JK))*6.5e-12  ! ESS_lg-20080722+ added reaction for Is
          RX1=1.E-28*ZT3REC**0.8*AIR(JL)
          RX2=8.8E-12
        RETHOH(JL)=(1.+IS(35,JK))*ZF3BOD(RX1,RX2)
        RETHO3(JL)=(1.+IS(36,JK))*ZFARR(9.14e-15,-2580.,ZTREC)
        RPAROH(JL)=(1.+IS(31,JK))*8.1E-13
        RRORA(JL)=ZFARR(1.E15,-8000.,ZTREC)
        RRORB(JL)=1.6E3

        ! ESS_20080722+ some added Is indices for the next selection of reactions

        RRORNO2(JL)=(1.+IS(52,JK))*1.5E-11   ! ESS_lg-20080722+ added reaction for Is
        RRXPAR(JL)=8E-11
        ROLEOH(JL)=(1.+IS(32,JK))*ZFARR(5.2E-12,504.,ZTREC)
        ROLEO3(JL)=(1.+IS(33,JK))*ZFARR(1.4E-14,-2100.,ZTREC)
        ROLENO3(JL)=(1.+IS(34,JK))*7.7E-15
        RMGLYOH(JL)=(1.+IS(37,JK))*1.7E-11
        RISOPOH(JL)=(1.+IS(38,JK))*ZFARR(2.54E-11,410.,ZTREC)
        RISOPO3(JL)=(1.+IS(39,JK))*ZFARR(12.3E-15,-2013.,ZTREC)
        RISOPNO3(JL)=(1.+IS(40,JK))*ZFARR(4.E-12,-446.,ZTREC)
        RISOPNO2(JL)=(1.+IS(53,JK))*1.5E-19
        RISPDOH(JL)=(1.+IS(54,JK))*6.1E-11
        RISPDO3(JL)=(1.+IS(55,JK))*4.2E-18
        RISPDNO3(JL)=(1.+IS(56,JK))*1.E-13
        RMTHCOH(JL)=(1.+IS(57,JK))*ZFARR(1.9E-11,176.,ZTREC)
        RMTHCO3(JL)=(1.+IS(58,JK))*ZFARR(1.4E-15,-2114.,ZTREC)
        RMTHCNO3(JL)=(1.+IS(59,JK))*ZFARR(1.5E-12,-1726.,ZTREC)
        RMVKOH(JL)=(1.+IS(60,JK))*ZFARR(4.1E-12,453.,ZTREC)
        RMVKO3(JL)=(1.+IS(61,JK))*ZFARR(7.5E-16,-1520.,ZTREC)
        RMC23NO(JL)=(1.+IS(62,JK))*RC23NO(JL)
        RMC23HO2(JL)=(1.+IS(63,JK))*RC23HO2(JL)
        RMC23NO2(JL)=(1.+IS(64,JK))*RC23NO2(JL)
        RMPAN(JL)=RPAN(JL)
        RMPANO3(JL)=(1.+IS(65,JK))*8.2E-18
        RMPANOH(JL)=(1.+IS(66,JK))*3.6E-12
        RMEKOH(JL)=(1.+IS(67,JK))*ZFARR(2.9E-13,413.,ZTREC) 

        ! mz_lg_20060215+ check the mecca code for the MEK+OH reaction rate which 
        !      appears to be a factor 10 larger compared to values of the CBM4 scheme
        !      RMEKOH(JL)=ZFARR(1.3e-12,-25.,ZTREC) {&1207}

        RACETOH(JL)=(1.+IS(68,JK))*ZFARR(1.7E-12,-600.,ZTREC)

        ! ESS_lg_20080722+ and back to some Is indices of TM3 scheme

        RXO2NO(JL)=(1.+IS(41,JK))*8.1E-12
        RXO2XO2(JL)=(1.+IS(42,JK))*ZFARR(1.7E-14,1300.,ZTREC)
        RXO2HO2(JL)=(1.+IS(43,JK))*ZFARR(7.7E-14,1300.,ZTREC)
        RXO2NNO(JL)=(1.+IS(44,JK))*8.1E-12
        RXO2N(JL)=ZFARR(1.7E-14,1300.,ZTREC)
        RXO2NHO2(JL)=(1.+IS(45,JK))*ZFARR(7.7E-14,1300.,ZTREC)
        RXO2NXO2(JL)=(1.+IS(46,JK))*ZFARR(3.5E-14,1300.,ZTREC) 

        ! ESS_lg_20080722+ skipped for the time being the BXO2N Is terms
        RBXO2NNO(JL)=8.1E-12 ! Duncan and Chameides, Table 2, 104

        ! WP-  changed RBXO2NNO

        RBXO2NNO(JL)=ZFARR(0.5*4.2E-12,180.,ZTREC)

        ! WP-  end

        RBXO2N(JL)=ZFARR(1.7E-14,1300.,ZTREC)
        RBXO2NHO2(JL)=ZFARR(7.7E-14,1300.,ZTREC)
        RBXO2NXO2(JL)=8.1E-12
        RBXO2NXO2N(JL)=ZFARR(3.4E-14,1300.,ZTREC)
        ! ESS_lg_20080722-

        ! LG-  sulfur chemistry reaction rates

           RX1=3.0E-31*ZT3REC**3.3*AIR(JL)
           RX2=1.5E-12
        ROHSO2(JL)=(1.+IS(48,JK))*ZF3BOD(RX1,RX2)
        ROHDMS(JL)=(1.+IS(49,JK))*ZFARR(1.2E-11,-260.,ZTREC)
        RNO3DMS(JL)=(1.+IS(50,JK))*ZFARR(1.9E-13,500.,ZTREC)

        ! WP-  added isontr+oh

        RISONTROH(JL)=(1.+IS(47,JK))*3.E-11

        ! WP-  end

        ! LG-  added hcooh+oh, and ch3coh2+oh

        RHCOOHOH(JL)=(1.+IS(74,JK))*4.5E-13
        RCH3CO2HOH(JL)=(1.+IS(75,JK))*ZFARR(4.E-13,200.,ZTREC)

        ! LG- updated reaction rates for c2o3 reactions

        RC23HO2(JL)=(1.+IS(30,JK))*ZFARR(1.3E-13,1040.,ZTREC)
           RX1=9.8E-12
           RX2=ZFARR(2.2E6,-3870.,ZTREC)
        RC23MO2A(JL)=(1.+IS(71,JK))*RX1/(1.+RX2)
        RC23MO2B(JL)=(1.+IS(71,JK))*RX1*(1.+1./RX2)

        ! LG-  ammonia reaction rates

        ROHNH3(JL)=ZFARR(1.7E-12,-710.,ZTREC) !1.56e-13 at 298K
        RNONH2(JL)=ZFARR(3.8E-12,+450.,ZTREC) !1.72e-11
        RNO2NH2(JL)=ZFARR(2.1E-12,650.,ZTREC) !1.86e-11
        RHO2NH2(JL)=3.4E-11
        RO2NH2(JL)=6.0E-21
        RO3NH2(JL)=ZFARR(4.3E-12,-930.,ZTREC) !1.89e-13 at 298K

        ! LG-  see TM3 code for definition of the parameters

        ! knh3so4 is uptake coefficient on H2SO4. 1 uptake of NH3 consumes 1 acid molecule.  
        ! 
        !  rr(jl,knh3so4)=het_nh3(jl)/1e-9/y(jl,iair) 

        RSO4NH3(JL)=0.

        ! LG-  added the monoterpenes destruction rates (Kostas Tsigaridis, personal 
        !      communication, March, 2002)

        RMATERPOH(JL)=(1.+IS(72,JK))*ZFARR(12.1E-12,444.,ZTREC)
        RMBTERPOH(JL)=(1.+IS(73,JK))*ZFARR(23.8E-12,357.,ZTREC)
        RMATERPO3(JL)=(1.+IS(74,JK))*ZFARR(1.01E-15,-732.,ZTREC)
        RMBTERPO3(JL)=(1.+IS(75,JK))*1.5E-17
        RMATERPNO3(JL)=(1.+IS(76,JK))*ZFARR(1.19E-12,490.,ZTREC)
        RMBTERPNO3(JL)=(1.+IS(77,JK))*2.51E-12

        RMTTERPOH(JL)=0.        ! (1.+IS(78,JK))*2.1E-10  ! mz_lg_20060403+ based on ratio terpinolene/
                               ! humulene of 38/28 minutes lifetime, check table 
        RMTTERPO3(JL)=(1.+IS(79,JK))*39.E-17   !180.E-17, 39e-17 for ~1 hour lifetime ! mz_lg_20060403+ based on ratio terpinolene/
                               ! humulene of 13/2 minutes lifetime
        RMTTERPNO3(JL)=0.       ! still needs to be defined

        ! LG-   added the sesquiterpenes (personal communications with Boris Bonn, April, 2003)

        RSQTERPOH(JL)=(1.+IS(80,JK))*2.9E-10
        RSQTERPO3(JL)=(1.+IS(81,JK))*1170.E-17
        RSQTERPNO3(JL)=(1.+IS(82,JK))*3.5E-11

        ! mz_lg_20060330+ added some first-order estimates of the sequiterpene oxidation
        !       products having 2- and 1- double bonds

        RSQTERP2BOH(JL)=2.9E-11
        RSQTERP2BO3(JL)=1.5E-16   ! mz_lg_20060530+ new reaction rates based on 
                                  ! numbers provided by R. Winterhalter, May 2006
        RSQTERP2BNO3(JL)=3.5E-12

        RSQTERP1BOH(JL)=2.9E-12
        RSQTERP1BO3(JL)=1.E-17    ! mz_lg_20060530+ new reaction rates based on 
                                  ! numbers provided by R. Winterhalter, May 2006
        RSQTERP1BNO3(JL)=3.5E-13
        ! LG-  end

        ! ESS_lg_20100628+ added a selection of anthropogenic alkenens, Vinayak Sinha, June 2010

        RHEXANEOH(JL)=5.8E-12
        RHEXANEO3(JL)=0.
        RHEXANENO3(JL)=0.

        RBUTADIENEOH(JL)=6.7E-11
        RBUTADIENEO3(JL)=0.
        RBUTADIENENO3(JL)=0.

        RTMBENZENEOH(JL)=7.7E-11
        RTMBENZENEO3(JL)=0.
        RTMBENZENENO3(JL)=0.        

        ! ESS_lg_20100628+ 

        ! LG-   12-2003, added some of the HONO chemistry
        !       Rate constant: NO+OH+M --> HONO+M
        !       9*10^(-12) cm^3 molecule^(-1) s^(-1) Seinfeld&Pandis (1998)

        RNOOH(JL)=9.E-12

        !       Rate constant: 2NO2+H2O--> HONO+HNO3 (1st order reaction)
        !       5.6*10^(-6) *(100/mixing height [m]) *s^(-1) Harrison et al. (1996)

        RNO2H2O(JL)=ZFARR(8.4E7,-2900.,ZTREC)/(AVO*H2OLTR)

        ! LG-   added the decomposition rate of HMHP into HCOOH and H2O (see Peter Neebs
        !       thesis, Tabel 3.15, pp 98) and that of HAHP into H2O2 and RCHO, see page 108
        !       (Chapter 6, summary and comparison..., Valverde's thesis), k1, k2

        RCH2OHO2H(JL)=6.E-4
        RRCHOHO2H(JL)=3.5E-3 ! 0. for sensitivity analysis, nocturnal H2O2 exchanges 

        ! ESS_lg_20080711+ added a reaction rate for isoprene products and HO2 to mimic the ISO2+HO2 reaction
        !       forming OH and ISOOH of the MIM scheme 
        RISPDHO2(JL)=ZFARR(2.22E-13,1300.,ZTREC)  ! see gas.eqn in subdirectory messy/mecca of ECHAM5/MESSy
        ! ESS_lg_20080711-       

 1104 CONTINUE

      ! - - -
      !      Set concentrations for time 0
      !

      DO 1105 JL=1,NLON
        O30(JL)=ZPM(JL,jk,idt_o3)
        CH40(JL)=ZPM(JL,jk,idt_ch4)
        CO0(JL)=ZPM(JL,jk,idt_co)
        HNO30(JL)=ZPM(JL,jk,idt_hno3)
        H2O20(JL)=ZPM(JL,jk,idt_h2o2)
        CH3O2H0(JL)=ZPM(JL,jk,idt_ch3o2h)
        CH2O0(JL)=ZPM(JL,jk,idt_ch2o)
        ALD20(JL)=ZPM(JL,jk,idt_ald2)
        PAR0(JL)=ZPM(JL,jk,idt_par)
        OLE0(JL)=ZPM(JL,jk,idt_ole)
        ETH0(JL)=ZPM(JL,jk,idt_eth)
        PAN0(JL)=ZPM(JL,jk,idt_pan)
        ACET0(JL)=ZPM(JL,jk,idt_acet)
        ISOP0(JL)=ZPM(JL,jk,idt_isop)
        MGLY0(JL)=ZPM(JL,jk,idt_mgly)
        ISOPRD0(JL)=ZPM(JL,jk,idt_isoprd)
        METHAC0(JL)=ZPM(JL,jk,idt_methac)
        MVK0(JL)=ZPM(JL,jk,idt_mvk)
        MEK0(JL)=ZPM(JL,jk,idt_mek)
        MPAN0(JL)=ZPM(JL,jk,idt_mpan)
        NITR0(JL)=ZPM(JL,jk,idt_ntr)

        ! LG-   sulfur chemistry

        DMS0(JL)=ZPM(JL,jk,idt_dms)
        SO20(JL)=ZPM(JL,jk,idt_so2)
        SAER0(JL)=ZPM(JL,jk,idt_so4)

        ! LG-   radon

            RADON0(JL)=ZPM(JL,jk,idt_rad)

        ! WP-   added isontr

        ISONTR0(JL)=ZPM(JL,jk,idt_isontr)

        ! LG-   added formic acid and acetic acid

        HCOOH0(JL)=ZPM(JL,jk,idt_hcooh)
        CH3CO2H0(JL)=ZPM(JL,jk,idt_ch3co2h)

        ! LG-   added ammonia chemistry

        NH20(JL)=ZPM(JL,jk,idt_nh2)
        NH30(JL)=ZPM(JL,jk,idt_nh3)
        NH40(JL)=ZPM(JL,jk,idt_nh4)

        ! LG-   added the monoterpenes as a group

        MATERP0(JL)=ZPM(JL,jk,idt_apin)
        MBTERP0(JL)=ZPM(JL,jk,idt_bpin)
        MTTERP0(JL)=ZPM(JL,jk,idt_mtterp) ! mz_lg_20060403+ extra monoterpene

        ! LG-   added the sesquiterpenes as a group

        SQTERP0(JL)=ZPM(JL,jk,idt_sqterp)

        ! mz_lg-20060330+ added the 2- and 1- double bond oxidation products

        SQTERP2B0(JL)=ZPM(JL,jk,idt_sqterp2b)
        SQTERP1B0(JL)=ZPM(JL,jk,idt_sqterp1b)

        ! LG-   added HONO

        HONO0(JL)=ZPM(JL,jk,idt_hono)

        ! LG-   added hydroxy-methyl/alkylhydroxyperoxide

        CH2OHO2H0(JL)=ZPM(JL,jk,idt_ch2oho2h)
        RCHOHO2H0(JL)=ZPM(JL,jk,idt_rchoho2h)

        ! ESS_LG_20100628+ added a selection of anthropogenic alkenes

        HEXANE0(JL)=ZPM(JL,jk,idt_hexane)
        BUTADIENE0(JL)=ZPM(JL,jk,idt_butadiene)
        TMBENZENE0(JL)=ZPM(JL,jk,idt_tmbenzene)

        ! ESS_LG_20100628-

        ZNO0(JL)=ZPM(JL,jk,idt_no)
        ZNO20(JL)=ZPM(JL,jk,idt_no2)
        ZNO30(JL)=MAX(ZPM(JL,jk,idt_no3),1.E-30)
        ZN2O50(JL)=MAX(ZPM(JL,jk,idt_n2o5),1.E-30)
        HNO40(JL)=MAX(ZPM(JL,jk,idt_hno4),1.E-30)
        OH0(JL)=ZPM(JL,jk,idt_oh)
        HO20(JL)=ZPM(JL,jk,idt_ho2)
        CH3O20(JL)=ZPM(JL,jk,idt_ch3o2)
        C2O30(JL)=ZPM(JL,jk,idt_c2o3)
        XO20(JL)=ZPM(JL,jk,idt_xo2)
        ROR0(JL)=ZPM(JL,jk,idt_ror)
        XO2N0(JL)=ZPM(JL,jk,idt_xo2n)
        RXPAR0(JL)=ZPM(JL,jk,idt_rxpar)
        RXPAR0(JL)=MIN(RXPAR0(JL),PAR0(JL))
        BXO2N0(JL)=ZPM(JL,jk,idt_bxo2n)
        MC3O30(JL)=ZPM(JL,jk,idt_mc3o3)
   
        ! LG-   end

        ! - - -
        !     Set concentrations for time t
        !
        O3(JL)=O30(JL)
        CH4(JL)=CH40(JL)
        CO(JL)=CO0(JL)
        HNO3(JL)=HNO30(JL)
        H2O2(JL)=H2O20(JL)
        CH3O2H(JL)=CH3O2H0(JL)
        CH2O(JL)=CH2O0(JL)
        ALD2(JL)=ALD20(JL)
        PAR(JL)=PAR0(JL)
        OLE(JL)=OLE0(JL)
        ETH(JL)=ETH0(JL)
        PAN(JL)=PAN0(JL)
        ACET(JL)=ACET0(JL)
        ISOP(JL)=ISOP0(JL)
        MGLY(JL)=MGLY0(JL)
        ISOPRD(JL)=ISOPRD0(JL)
        METHAC(JL)=METHAC0(JL)
        MVK(JL)=MVK0(JL)
        MEK(JL)=MEK0(JL)
        MPAN(JL)=MPAN0(JL)
        NITR(JL)=NITR0(JL)
        ZNO(JL)=ZNO0(JL)
        ZNO2(JL)=ZNO20(JL)
        ZNO3(JL)=ZNO30(JL)
        ZN2O5(JL)=ZN2O50(JL)
        HNO4(JL)=HNO40(JL)
        OH(JL)=OH0(JL)
        HO2(JL)=HO20(JL)
        CH3O2(JL)=CH3O20(JL)
        C2O3(JL)=C2O30(JL)
        XO2(JL)=XO20(JL)
        ROR(JL)=ROR0(JL)
        XO2N(JL)=XO2N0(JL)
        RXPAR(JL)=RXPAR0(JL)
        BXO2N(JL)=BXO2N0(JL)
        MC3O3(JL)=MC3O30(JL)

        ODDN0(JL)=ZNO(JL)+ZNO2(JL)+ZNO3(JL)+2*ZN2O5(JL)+HNO4(JL)+ &
     &            PAN(JL)+MPAN(JL)+NITR(JL)

        ! LG-   sulfur chemistry

        DMS(JL)=DMS0(JL)
        SO2(JL)=SO20(JL)
        SAER(JL)=SAER0(JL)

        ! LG-   radon

            RADON(JL)=RADON0(JL)

        ! WP-   added isontr

        ISONTR(JL)=ISONTR0(JL)

        ! WP-   added formic acid and acetic acid

        HCOOH(JL)=HCOOH0(JL)
            CH3CO2H(JL)=CH3CO2H0(JL)

        ! LG-   ammonia chemistry
    
        NH2(JL)=NH20(JL)
        NH3(JL)=NH30(JL)
        NH4(JL)=NH40(JL)

        ! LG-   added the monoterpenes as a group

        MATERP(JL)=MATERP0(JL)
        MBTERP(JL)=MBTERP0(JL)
        MTTERP(JL)=MTTERP0(JL) ! mz_lg_20060403+ extra monoterpene

        ! LG-   added the sesquiterpenes as a group

        SQTERP(JL)=SQTERP0(JL)

        ! mz_lg_20060330+

        SQTERP2B(JL)=SQTERP2B0(JL)
        SQTERP1B(JL)=SQTERP1B0(JL)

        ! LG-   added HONO

        HONO(JL)=HONO0(JL)

        ! LG-   added hydroxy-methyl/alkylhydroxyperoxide

        CH2OHO2H(JL)=CH2OHO2H0(JL)
        RCHOHO2H(JL)=RCHOHO2H0(JL)

        ! ESS_LG_20100628+ added a selection of anthropogenic alkenes

        HEXANE(JL)=HEXANE0(JL)
        BUTADIENE(JL)=BUTADIENE0(JL)
        TMBENZENE(JL)=TMBENZENE0(JL)

 1105 CONTINUE

      !
      !  **** CHEMISTRY STARTS HERE *****
      !

      DO 2104 JL=1,NLON

          ! LG- different iteration mechanism compared to the original ECHAM/CBM4
          !     code, the number of iterations is determined by the relative
          !     difference between the new and old concentration of NOy 

          ! ESS_lg_20120920+ added the option to skip gas-phase chemistry
          ! this allows to compare the simulations with the extensive tracer list
          ! with and without the role of gas-phase chemistry

          IF (.NOT.l_xtsurf_veg_mlay_chem_reactions) GOTO 2104

		  ! MAQ_2016923+ skip the chemistry inside the canopy
		  IF (.NOT.l_chem_canopy.AND.jk>1) GOTO 2104
		  
          ITER=0
          OLDNOY=ODDN0(JL)

 9999     CONTINUE

          ! --- First group: NO NO2 O3
          P1=RJBNO3(JL)*ZNO3(JL) &
            ! LG- added the production of NO through the photolysis of HONO
     &      +RJHONO(JL)*HONO(JL)

          R12=0.2*RISOPNO2(JL)*ISOP(JL)
          R21=RHO2NO(JL)*HO2(JL)+RMO2NO(JL)*CH3O2(JL)+ &
     &      RC23NO(JL)*C2O3(JL)+RXO2NO(JL)*XO2(JL)+ &
     &      RMC23NO(JL)*MC3O3(JL)
          XL1=RNONO3(JL)*ZNO3(JL)+RXO2NNO(JL)*XO2N(JL)+ &
     &       RBXO2NNO(JL)*BXO2N(JL) &

             ! LG- added the production of HONO from NO+OH
     &       +RNOOH(JL)*OH(JL)

          P2=RJHNO3(JL)*HNO3(JL)+RJN2O5(JL)*ZN2O5(JL)+ &
     &      RN2O5(JL)*ZN2O5(JL)+HNO4(JL)* &
     &      (RJHNO4(JL)+RHNO4M(JL)+RHNO4OH(JL)*OH(JL))+ &
     &      2*RNONO3(JL)*ZNO3(JL)*ZNO(JL)+ &
     &      PAN(JL)*(RPAN(JL)+RJPAN(JL))+ &
     &      MPAN(JL)*(RMPAN(JL)+RJPAN(JL)+0.70*RMPANO3(JL)*O3(JL))+ &
     &      ZNO3(JL)*(RJANO3(JL)+ROLENO3(JL)*OLE(JL)+ &
     &       0.2*RISOPNO3(JL)*ISOP(JL))+ &
     &      RJNITR(JL)*NITR(JL) &

            ! WP- added isontr:isontr+oh=no2
     &     +RISONTROH(JL)*ISONTR(JL)*OH(JL)

          XL2=RNO2OH(JL)*OH(JL)+RNO2NO3(JL)*ZNO3(JL)+ &
     &     RNO2HO2(JL)*HO2(JL)+RNO2O3(JL)*O3(JL)+ &
     &     RC23NO2(JL)*C2O3(JL)+RRORNO2(JL)*ROR(JL)+ &
     &     0.8*RISOPNO2(JL)*ISOP(JL)+RMC23NO2(JL)*MC3O3(JL)

          P3=RJANO3(JL)*ZNO3(JL)+ROHOH(JL)*OH(JL)*OH(JL) &
         
            ! LG- added acetic acid:o3 formation in reaction 
            !     c2o3+ch3o2->0.3*(ch3co2h+o3) (the ratio is already accounted
            !     for in the reaction rate)

     &      +RC23HO2(JL)*C2O3(JL)*HO2(JL)
         
          XL3=RO3HO2(JL)*HO2(JL)+RO3OH(JL)*OH(JL)+ &
     &      RNO2O3(JL)*ZNO2(JL)+RJO3D(JL)+RETHO3(JL)*ETH(JL)+ &
     &      RISOPO3(JL)*ISOP(JL)+ROLEO3(JL)*OLE(JL)+ &
     &      RISPDO3(JL)*ISOPRD(JL)+RMTHCO3(JL)*METHAC(JL)+ &
     &      RMVKO3(JL)*MVK(JL)+RMPANO3(JL)*MPAN(JL)+ &
     
            ! mz_lg_20060814+ added the ozone sinks through the reactions with the 
            !    terpenes

     &      RMATERPO3(JL)*MATERP(JL)+RMBTERPO3(JL)*MBTERP(JL)+ &
     &      RSQTERPO3(JL)*SQTERP(JL)+ &
     &      RSQTERP2BO3(JL)*SQTERP2B(JL)+RSQTERP1BO3(JL)*SQTERP1B(JL)

          X1=ZNO0(JL)+P1*DT
          X2=ZNO20(JL)+P2*DT
          X3=O30(JL)+P3*DT
          C1=1.+XL1*DT
          C2=1.+XL2*DT
          C3=1.+XL3*DT
          Y2=RNOO3(JL)*DT/(C1*C3)
          XJT=RJNO2(JL)*DT
          R21T=R21*DT
          R12T=R12*DT
          R12TC=R12T/C2
          R21TC=R21T/C1
          XJTC=XJT/C2
          ! --- resolving for unknown x
          ACUB=-1.*Y2*(1.+R12TC+R21TC)
          BCUB=1.+R12TC+XJTC+R21TC+ &
     &           Y2*(R12TC*(X1-X2)+2.*R21TC*X1+X1+X3)
          CCUB=X2*(R12TC+XJTC)-X1*R21TC+Y2*X1* &
     &           (X2*R12TC-X3-R21TC*X1)
          CUBDET=BCUB*BCUB-4.*ACUB*CCUB
          CUBDET=MAX(CUBDET,1.E-20)
          DNO2=(-1.*BCUB+SQRT(CUBDET))/(2.*ACUB)
          ZNO2(JL)=(X2+DNO2)/C2
          ZNO(JL)=(X1-DNO2)/C1
          O3(JL)=(X3+XJTC*(X2+DNO2))/(C3+Y2*C3*(X1-DNO2))

          !  --- Second group: HO2 OH HNO4
          R57=RJHNO4(JL)+RHNO4M(JL)
          R56=RCOOH(JL)*CO(JL)+RO3OH(JL)*O3(JL)+ &
     &       RHPOH(JL)*H2O2(JL)+RFRMOH(JL)*CH2O(JL)+ &
     &       RH2OH(JL)*AIR(JL)+RETHOH(JL)*ETH(JL)+ &
     &       ROLEOH(JL)*OLE(JL)+ &
     &       0.91*RISOPOH(JL)*ISOP(JL)+ &
     &       0.68*RISPDOH(JL)*ISOPRD(JL)+0.5*RMTHCOH(JL)*METHAC(JL)+ &
     &       0.3*RMVKOH(JL)*MVK(JL)+0.11*RPAROH(JL)*PAR(JL)+ &
     &       0.4*RMPANOH(JL)*MPAN(JL)

          ! LG- 12092004+ the HO2 production term from the OH reactions

          PHO2OH(JL,JK)=R56

          ! LG- March 2002, adding the production of OH in the ozonolysis reactions of 
          !     the terpenes, the yields are 0.76 and 0.33 for the ozonolysis of alpha 
          !     and beta pinenes respectively (Kostas Tsigaridis, personal communication)
          !     These reactions might be an important nocturnal source of OH

          R65=RHO2NO(JL)*ZNO(JL)+RO3HO2(JL)*O3(JL)

          ! LG- modified based on personal communications with Kostas, 23 July
          !     2002, see also the change in the OH production term P6!

          ! LG- 12092004+ modified by adding an extra term that resembles the OH
          !     production related to the terpene ozonolysis. This term was originally
          !     included in the R65 term but the latter is also used to calculate the
          !     HO2 production/destruction and the terpene ozonolysis reactions do not
          !     result in an HO2 destruction, as it was initially suggested by including
          !     the reactions in the R65 term!

          RTERPO3= &
     &         0.76*RMATERPO3(JL)*O3(JL) &
     &        +0.33*RMBTERPO3(JL)*O3(JL) &
     &        +1.0 *RMTTERPO3(JL)*O3(JL) & ! mz_lg_20060403+ 

          ! LG- further modified based on personal communications with Boris Bonn, 
          !     April, 2003. The yield of 0.22 is representative for the species alpha-
          !     humulene, which reactions rates with OH, O3 and NO3 are also applied

     &        +0.22*RSQTERPO3(JL)*O3(JL) & ! check out this yield, a 50%
                       ! uncertainty would require a substantially lower
                       ! emission for a significant OH production 

     &        +0.1*RSQTERP2BO3(JL)*O3(JL) & ! what is the yield of 2-double bond sesquiterpene-O3 product ?
     &        +0.1*RSQTERP1BO3(JL)*O3(JL)   ! what is the yield of 1-double bond sesquiterpene-O3 product ?


          ! LG- assigning the OH production rate for interpretation of OH 
          !     production rates  

          ! LG- 12092004+ the OH production term from the HO2 reactions
          POHHO2(JL,JK)=R65

          R75=RNO2HO2(JL)*ZNO2(JL)
          P5=2*RJBCH2O(JL)*CH2O(JL)+RMO2NO(JL)*CH3O2(JL)*ZNO(JL)+ &
     &      0.66*RMO2MO2(JL)*CH3O2(JL)*CH3O2(JL)+ &
     &      RJMEPE(JL)*CH3O2H(JL)+2*RJALD2(JL)*ALD2(JL)+ &
     &     C2O3(JL)* &
     &       (RC23NO(JL)*ZNO(JL)+2*RC23C23(JL)*C2O3(JL)+ &
     &       RC23MO2(JL)*CH3O2(JL))+ &
     &     O3(JL)* &
     &       (0.44*ROLEO3(JL)*OLE(JL)+0.12*RETHO3(JL)*ETH(JL)+ &
     &       0.23*RISPDO3(JL)*ISOPRD(JL)+0.1*RMTHCO3(JL)*METHAC(JL)+ &
     &       0.08*RMPANO3(JL)*MPAN(JL))+ &
     &     ZNO3(JL)* &
     &       (0.8*RISOPNO3(JL)*ISOP(JL)+RISPDNO3(JL)*ISOPRD(JL)+ &
     &       0.5*RMTHCNO3(JL)*METHAC(JL)+RFRMNO3(JL)*CH2O(JL))+ &
     &     0.8*RISOPNO2(JL)*ISOP(JL)*ZNO2(JL)+ &
     &     RRORB(JL)*ROR(JL)+RJMEK(JL)*MEK(JL)+ &
     &     1.17*RRORA(JL)*ROR(JL)+ &
     &     RJNITR(JL)*NITR(JL) &
     
           ! LG- added formic acid:hcooh+oh -> ho2
     &     +RHCOOHOH(JL)*HCOOH(JL)*OH(JL) &
     
           ! LG- added acetic acid:production ho2 in reaction c2o3+ch3o2->ch2o+ho2

     &     +RC23MO2A(JL)*C2O3(JL)*CH3O2(JL) 
     
          XL5=R65+R75+RMO2HO2(JL)*CH3O2(JL)+RHO2OH(JL)*OH(JL)+ &
     
          ! LG- added acetic acid: change of constant 0.21 in 1 (see paper Duncan),
          !     There is production of 0.79 ho2 from 1 ho2, thus a net loss of
          !     0.21 ho2, due to the removal of reaction 53 this yield is changed 
          !     to 1
     &     1.*RC23HO2(JL)*C2O3(JL)+ &

     &      RXO2HO2(JL)*XO2(JL)+RXO2NHO2(JL)*XO2N(JL)+ &
     &      RBXO2NHO2(JL)*BXO2N(JL)+RMC23HO2(JL)*MC3O3(JL)

          R66=2.*RHO2HO2(JL)
          X5=HO20(JL)+P5*DT

          ! ESS_lg_20080714+ added the calculated of contribution of isoprene oxidation to the production of
          !         XO2, which is needed to estimate the production of OH from the isoprene related XO2+HO2 
          !         reaction.

          FXO2_ISOP=( &
     &      OH(JL)*(0.99*RISOPOH(JL)*ISOP(JL)+0.68*RISPDOH(JL)*ISOPRD(JL))+ &
     &      O3(JL)*(0.2*RISPDO3(JL)*ISOPRD(JL)+0.2*RISOPO3(JL)*ISOP(JL))+ &
     &      ZNO3(JL)*(RISPDNO3(JL)*ISOPRD(JL)+RISOPNO3(JL)*ISOP(JL))+ &
     &            ISOP(JL)*RISOPNO2(JL)*ZNO2(JL) )/ &
     &           ( RJALD2(JL)*ALD2(JL)+C2O3(JL)*(RC23NO(JL)*ZNO(JL)+2*RC23C23(JL)*C2O3(JL))+ &
     &      OH(JL)*(RACETOH(JL)*ACET(JL)+ROLEOH(JL)*OLE(JL)+ &
     &      RETHOH(JL)*ETH(JL)+0.87*RPAROH(JL)*PAR(JL)+ &
     &      RMGLYOH(JL)*MGLY(JL)+0.99*RISOPOH(JL)*ISOP(JL)+ &
     &      0.68*RISPDOH(JL)*ISOPRD(JL)+0.5*RMTHCOH(JL)*METHAC(JL)+ &
     &      1.5*RMEKOH(JL)*MEK(JL)+RMVKOH(JL)*MVK(JL)+ &
     &      RMPANOH(JL)*MPAN(JL))+ &
     &      O3(JL)* &
     &       (0.22*ROLEO3(JL)*OLE(JL)+0.2*RISPDO3(JL)*ISOPRD(JL)+ &
     &        0.2*RISOPO3(JL)*ISOP(JL)+0.1*RMTHCO3(JL)*METHAC(JL)+ &
     &        0.05*RMVKO3(JL)*MVK(JL))+ &
     &      ZNO3(JL)* &
     &       (0.91*ROLENO3(JL)*OLE(JL)+RISPDNO3(JL)*ISOPRD(JL)+ &
     &        RISOPNO3(JL)*ISOP(JL)+0.5*RMTHCNO3(JL)*METHAC(JL))+ &
     &      ISOP(JL)*RISOPNO2(JL)*ZNO2(JL)+0.76*RRORA(JL)*ROR(JL)+ &
     &      RJMEK(JL)*MEK(JL) )

          ! ESS_lg_20080714-

          ! LG- modified based on personal communications with Kostas, July 2002

          P6=RJHNO3(JL)*HNO3(JL)+2.*RJO3D(JL)*O3(JL)+ &
     &     2.*RJH2O2(JL)*H2O2(JL)+RJMEPE(JL)*CH3O2H(JL)+ &
     &     O3(JL)* &
     &       (0.4*RISPDO3(JL)*ISOPRD(JL)+0.1*ROLEO3(JL)*OLE(JL)+ &
     &       0.2*RISOPO3(JL)*ISOP(JL)+0.1*RMTHCO3(JL)*METHAC(JL)+ &
     &       0.05*RMVKO3(JL)*MVK(JL)+0.04*RMPANO3(JL)*MPAN(JL)) &
     &      +0.76*RMATERPO3(JL)*MATERP(JL)*O3(JL) &
     &      +0.33*RMBTERPO3(JL)*MBTERP(JL)*O3(JL) &
     &      +1.0 *RMTTERPO3(JL)*MTTERP(JL)*O3(JL) & ! mz_lg_20060403+ 

          ! LG- futher modified based on personal communications with Boris Bonn, 
          !     April 2003, to include the ozonolysis of the sesquiterpenes

     &      +0.22*RSQTERPO3(JL)*SQTERP(JL)*O3(JL) &
     &      +0.35*RSQTERP2BO3(JL)*SQTERP2B(JL)*O3(JL) & ! yield of 2-double bond sesquiterpene-O3 product
                                                        ! taken as middle of range given by R. Winterhalter, May 2006, 
     &      +0.35*RSQTERP1BO3(JL)*SQTERP1B(JL)*O3(JL) & ! yield of 1-double bond sesquiterpene-O3 product
                                                        ! taken as middle of range given by R. Winterhalter, May 2006, 

          ! LG- and added the production of OH due to HONO photodissociation
     &      +RJHONO(JL)*HONO(JL) &

          ! mz_lg_20060614+ added the production of OH from the photolysis of Glycolaldehyde,
          !     paper by Magneron et al., 2005. The implementation as it is suggests that we assume
          !     that all the higher aldehydes consist of Glycolaldehyde. The competing reaction of
          !     Glycolaldehyde with OH has been considered in through the reaction ALD2*OH. What has not
          !     been considered yet are the oxidation products, such as methanol

          !     &      +0.25*RJALD2(JL)*ALD2(JL)  ! assuming a yield of 25%

          ! mz_lg_200601128+ added the production of OH from the XO2+HO2 reaction

     &      +0.4*RXO2HO2(JL)*XO2(JL)*HO2(JL) &
     &      +FXO2_ISOP*4.5*RXO2HO2(JL)*XO2(JL)*HO2(JL) 

          ! ESS_lg_20080711+ recycling of OH in isoprene oxidation an extra yield of 
          !   3 OH for for an FXO2_ISOP of 1. The typical FXO2_ISOP ranges between 0.5-1

          XL6=RHNO4OH(JL)*HNO4(JL)+RHO2OH(JL)*HO2(JL)+ &
     &       RNO2OH(JL)*ZNO2(JL)+ROHHNO3(JL)*HNO3(JL)+ &
     &       R56+RCH4OH(JL)*CH4(JL)+ROHPCAT(JL)*CH3O2H(JL)+ &
     &       2*ROHOH(JL)*OH(JL)+0.89*RPAROH(JL)*PAR(JL)+ &
     &       RMGLYOH(JL)*MGLY(JL)+0.09*RISOPOH(JL)*ISOP(JL)+ &
     &       0.32*RISPDOH(JL)*ISOPRD(JL)+ &
     &       0.5*RMTHCOH(JL)*METHAC(JL)+RMEKOH(JL)*MEK(JL)+ &
     &       0.7*RMVKOH(JL)*MVK(JL)+RALD2OH(JL)*ALD2(JL)+ &
     &       0.6*RMPANOH(JL)*MPAN(JL)+RACETOH(JL)*ACET(JL) &

          ! LG- including oxidation of monoterpenes

     &      +RMATERPOH(JL)*MATERP(JL) &
     &      +RMBTERPOH(JL)*MBTERP(JL) &
     &      +RMTTERPOH(JL)*MTTERP(JL) & ! mz_lg_20060403+ 

          ! LG- futher modified based on personal communications with Boris Bonn, 
          !     April 2003, to include the ozonolysis of the sesquiterpenes

     &      +RSQTERPOH(JL)*SQTERP(JL) &
     &      +RSQTERP2BOH(JL)*SQTERP2B(JL) & ! mz_lg_20060330+
     &      +RSQTERP1BOH(JL)*SQTERP1B(JL) &

          ! LG- including HONO production through NO+OH

     &      +RNOOH(JL)*ZNO(JL)

          X6=OH0(JL)+P6*DT
          C6=1.+XL6*DT
          XL7=R57+RHNO4OH(JL)*OH(JL)
          C7=1.+XL7*DT
          Y1=R57/C7     ! R57=RJHNO4(JL)+RHNO4M(JL)
          Y2=R56/C6     ! R56=RCOOH(JL)*CO(JL)+RO3OH(JL)*O3(JL)+ etc.
          ACUB=R66*DT   ! R66=2.*RHO2HO2(JL)
          BCUB=1.+XL5*DT-DT2*(Y1*R75+Y2*R65) ! R75=RNO2HO2(JL)*ZNO2(JL)
                                             ! R65=RHO2NO(JL)*ZNO(JL)+RO3HO2(JL)*O3(JL)
          CCUB=-1.*X5-DT*(Y1*HNO40(JL)+Y2*X6)
          CUBDET=BCUB*BCUB-4.*ACUB*CCUB
          CUBDET=MAX(CUBDET,1.E-20)
          HO2(JL)=(-1.*BCUB+SQRT(CUBDET))/(2.*ACUB)
          OH(JL)=(X6+R65*HO2(JL)*DT)/C6
          HNO4(JL)=(HNO40(JL)+R75*DT*HO2(JL))/C7 ! R75=RNO2HO2(JL)*ZNO2(JL)

          ! LG- added the production of HONO including the reaction of NO and OH
          !     and the reaction of NO2 and H2O. H2O is not included in the 
          !     production term since this is already implictly considered (the
          !     units of RNO2H2O is s-1 instead of cm3 molecules-1 s-1)

          ! LG- it still must be find out if the twice the NO2 concentration or
          !     only once the NO2 concentration should be used (check the
          !     parameterization). The HNO3 being formed is not gas-phase HNO3
          !     but HNO3 present in the water/aerosol: 2NO2+H2O -> HONO+HNO3.
          !     See also paper Harrisson et al, JGR, 1996. 

          PHONO=RNOOH(JL)*ZNO(JL)*OH(JL)+RNO2H2O(JL)*2.*ZNO2(JL) &
     &          +0.8*RBXO2NNO(JL)*BXO2N(JL)*ZNO(JL) ! mz_lg_20060405+ 
          XLHONO=RJHONO(JL)
          HONO(JL)=(HONO0(JL)+PHONO*DT)/(1.+XLHONO*DT)

          ! --- Third group: NO3 N2O5

          R89=RJN2O5(JL)+RN2O5(JL)
          P8=ROHHNO3(JL)*HNO3(JL)*OH(JL)+ &
     &        RNO2O3(JL)*ZNO2(JL)*O3(JL)+ &
     &        RMPANOH(JL)*MPAN(JL)*OH(JL)
          XL8=RJBNO3(JL)+RJANO3(JL)+RNONO3(JL)*ZNO(JL)+ &
     &        RNO2NO3(JL)*ZNO2(JL)+RALD2NO3(JL)*ALD2(JL)+ &
     &        RFRMNO3(JL)*CH2O(JL)+ROLENO3(JL)*OLE(JL)+ &
     &        RISOPNO3(JL)*ISOP(JL)+RISPDNO3(JL)*ISOPRD(JL)+ &
     &        RMTHCNO3(JL)*METHAC(JL)
          X4=ZNO30(JL)+P8*DT
          C5=1.+XL8*DT
          R98=RNO2NO3(JL)*ZNO2(JL)
          XL9=RJN2O5(JL)+RN2O5(JL)+RN2O5AQ(JL)+RN2O5L(JL)
          C6=1.+XL9*DT
          C7=(C5*C6-R89*R98*DT2)
          ZN2O5(JL)=(C5*ZN2O50(JL)+R98*DT*X4)/C7
          ZNO3(JL)=(C6*X4+R89*DT*ZN2O50(JL))/C7
          PCH3O2=RCH4OH(JL)*CH4(JL)*OH(JL)+ &
     &      OH(JL)*(ROHPCAT(JL)*CH3O2H(JL))+RJACET(JL)*ACET(JL)+ &
     &      0.5*RJMGLY(JL)*MGLY(JL)

          ! rc23mo2 not considered since P=L
          XLCH3O2=RMO2NO(JL)*ZNO(JL)+RMO2HO2(JL)*HO2(JL)+ &
     &      2*RMO2MO2(JL)*CH3O2(JL) &
     
           ! LG- added acetic acid: destruction of ch3o2 in reaction with c2o3

     &     +RC23MO2A(JL)*C2O3(JL)+RC23MO2B(JL)*C2O3(JL) 
     
          CH3O2(JL)=(CH3O20(JL)+PCH3O2*DT)/(1.+XLCH3O2*DT)

          ! --- Some others C2O3 XO2 etc
          ! C2O3

          PC2O3=ALD2(JL)*(RALD2OH(JL)*OH(JL)+RALD2NO3(JL)*ZNO3(JL))+ &
     &     OH(JL)* &
     &     (RMGLYOH(JL)*MGLY(JL)+ &
     &     0.7*RMVKOH(JL)*MVK(JL)+0.5*RMEKOH(JL)*MEK(JL))+ &
     &     RJACET(JL)*ACET(JL)+PAN(JL)*(RPAN(JL)+RJPAN(JL))+ &
     &     RJMEK(JL)*MEK(JL)+0.5*RJMGLY(JL)*MGLY(JL)+ &
     &     O3(JL)* &
     &      (0.17*RISPDO3(JL)*ISOPRD(JL)+0.1*RMTHCO3(JL)*METHAC(JL)+ &
     &       0.70*RMPANO3(JL)*MPAN(JL))
          XLC2O3=RC23NO(JL)*ZNO(JL)+RC23NO2(JL)*ZNO2(JL)+ &
     &     2*RC23C23(JL)*C2O3(JL)+RC23HO2(JL)*HO2(JL)+ &
     
           ! LG- added acetic acid: using both rc23mo2a and b     
    
     &     RC23MO2A(JL)*CH3O2(JL)+RC23MO2B(JL)*CH3O2(JL)

          C2O3(JL)=(C2O30(JL)+PC2O3*DT)/(1.+XLC2O3*DT)
          ! XO2
          PXO2=RJALD2(JL)*ALD2(JL)+ &
     &     C2O3(JL)* &
     &     (RC23NO(JL)*ZNO(JL)+2*RC23C23(JL)*C2O3(JL))+ &

           ! LG- removed acetic acid: removal of reaction 53 (see paper Duncan)

           !     &    0.79*RC23HO2(JL)*HO2(JL))+

     &     OH(JL)* &
     &      (RACETOH(JL)*ACET(JL)+ROLEOH(JL)*OLE(JL)+ &
     &       RETHOH(JL)*ETH(JL)+0.87*RPAROH(JL)*PAR(JL)+ &
     &       RMGLYOH(JL)*MGLY(JL)+0.99*RISOPOH(JL)*ISOP(JL)+ &
     &       0.68*RISPDOH(JL)*ISOPRD(JL)+0.5*RMTHCOH(JL)*METHAC(JL)+ &
     &       1.5*RMEKOH(JL)*MEK(JL)+RMVKOH(JL)*MVK(JL)+ &
     &       RMPANOH(JL)*MPAN(JL))+ &
     &     O3(JL)* &
     &      (0.22*ROLEO3(JL)*OLE(JL)+0.2*RISPDO3(JL)*ISOPRD(JL)+ &
     &       0.2*RISOPO3(JL)*ISOP(JL)+0.1*RMTHCO3(JL)*METHAC(JL)+ &
     &       0.05*RMVKO3(JL)*MVK(JL))+ &
     &     ZNO3(JL)* &
     &      (0.91*ROLENO3(JL)*OLE(JL)+RISPDNO3(JL)*ISOPRD(JL)+ &
     &       RISOPNO3(JL)*ISOP(JL)+0.5*RMTHCNO3(JL)*METHAC(JL))+ &
     &     ISOP(JL)*RISOPNO2(JL)*ZNO2(JL)+0.76*RRORA(JL)*ROR(JL)+ &
     &     RJMEK(JL)*MEK(JL)

          XLXO2=RXO2NO(JL)*ZNO(JL)+2*RXO2XO2(JL)*XO2(JL)+ &
     &      RXO2HO2(JL)*HO2(JL)+RXO2NXO2(JL)*XO2N(JL)+ &
     &     RBXO2NXO2(JL)*BXO2N(JL)
          XO2(JL)=(XO20(JL)+PXO2*DT)/(1.+XLXO2*DT)

          ! shorties
          PROR=0.76*RPAROH(JL)*PAR(JL)*OH(JL)+0.03*RRORA(JL)*ROR(JL)
          XLROR=RRORA(JL)+RRORB(JL)+RRORNO2(JL)*ZNO2(JL)
          ROR(JL)=(ROR0(JL)+PROR*DT)/(1.+XLROR*DT)
          PRXPAR=3.39*RRORA(JL)*ROR(JL)+ROLEO3(JL)*OLE(JL)*O3(JL)+ &
     &      ROLENO3(JL)*OLE(JL)*ZNO3(JL)+ &
     &      OH(JL)* &
     &        (0.1*RPAROH(JL)*PAR(JL)+ROLEOH(JL)*OLE(JL))
          XLRXPAR=RRXPAR(JL)*PAR(JL)
          RXPAR(JL)=(RXPAR0(JL)+PRXPAR*DT)/(1.+XLRXPAR*DT)
          PXO2N=0.13*RPAROH(JL)*PAR(JL)*OH(JL)+0.03*RRORA(JL)*ROR(JL)+ &
     &      0.09*ROLENO3(JL)*OLE(JL)*ZNO3(JL)
          XLXO2N=RXO2NNO(JL)*ZNO(JL)+2*RXO2N(JL)*XO2N(JL)+ &
     &      RXO2NHO2(JL)*HO2(JL)+RXO2NXO2(JL)*XO2(JL)+ &
     &      RBXO2NXO2N(JL)*BXO2N(JL)
          XO2N(JL)=(XO2N0(JL)+PXO2N*DT)/(1.+XLXO2N*DT)
          PBXO2N=1.*RISOPOH(JL)*ISOP(JL)*OH(JL) ! mz_lg_20060406+ originally 0.09, set to 1 for tests
          XLBXO2N=RBXO2NNO(JL)*ZNO(JL)+2*RBXO2N(JL)*BXO2N(JL)+ &
     &      RBXO2NHO2(JL)*HO2(JL)+RBXO2NXO2(JL)*XO2(JL)+ &
     &      RBXO2NXO2N(JL)*XO2N(JL)
          BXO2N(JL)=(BXO2N0(JL)+PBXO2N*DT)/(1.+XLBXO2N*DT)
          POLE=0.
          XLOLE=ROLEOH(JL)*OH(JL)+ROLEO3(JL)*O3(JL)+ROLENO3(JL)*ZNO3(JL)
          OLE(JL)=(OLE0(JL)+POLE*DT)/(1.+XLOLE*DT)
          PMGLY= &
     &      OH(JL)* &
     &       (0.15*RISPDOH(JL)*ISOPRD(JL)+0.08*RMTHCOH(JL)*METHAC(JL)+ &
     &        0.3*RMVKOH(JL)*MVK(JL)+RACETOH(JL)*ACET(JL))+ &
     &      O3(JL)* &
     &       (0.65*RISPDO3(JL)*ISOPRD(JL)+0.9*RMTHCO3(JL)*METHAC(JL)+ &
     &        0.95*RMVKO3(JL)*MVK(JL))
          XLMGLY=RMGLYOH(JL)*OH(JL)+RJMGLY(JL)
          MGLY(JL)=(MGLY0(JL)+PMGLY*DT)/(1.+XLMGLY*DT)
          PISOP=0.
          XLISOP=RISOPOH(JL)*OH(JL)+RISOPO3(JL)*O3(JL)+ &
     &      RISOPNO3(JL)*ZNO3(JL)+RISOPNO2(JL)*ZNO2(JL)
          ISOP(JL)=(ISOP0(JL)+PISOP*DT)/(1.+XLISOP*DT)
          PISOPRD=ISOP(JL)* &
     &     (0.36*RISOPOH(JL)*OH(JL)+0.1*RISOPO3(JL)*O3(JL)+ &
     &      0.2*RISOPNO3(JL)*ZNO3(JL)+0.2*RISOPNO2(JL)*ZNO2(JL))
          XLISOPRD=RISPDOH(JL)*OH(JL)+RISPDO3(JL)*O3(JL)+ &
     &      RISPDNO3(JL)*ZNO3(JL)
          ISOPRD(JL)=(ISOPRD0(JL)+PISOPRD*DT)/(1.+XLISOPRD*DT)
          PMETHAC=ISOP(JL)* &
     &     (0.23*RISOPOH(JL)*OH(JL)+0.39*RISOPO3(JL)*O3(JL))
          XLMETHAC=RMTHCOH(JL)*OH(JL)+RMTHCO3(JL)*O3(JL)+ &
     &      RMTHCNO3(JL)*ZNO3(JL)
          METHAC(JL)=(METHAC0(JL)+PMETHAC*DT)/(1.+XLMETHAC*DT)
          PMVK=ISOP(JL)*(0.32*RISOPOH(JL)*OH(JL)+0.16*RISOPO3(JL)*O3(JL))
          XLMVK=RMVKOH(JL)*OH(JL)+RMVKO3(JL)*O3(JL)
          MVK(JL)=(MVK0(JL)+PMVK*DT)/(1.+XLMVK*DT)
          PMEK=0.42*RMTHCOH(JL)*OH(JL)*METHAC(JL)+ &
     &      ISOPRD(JL)*(0.48*RISPDOH(JL)*OH(JL)+0.28*RISPDO3(JL)*O3(JL))
          XLMEK=RMEKOH(JL)*OH(JL)+RJMEK(JL)
          MEK(JL)=(MEK0(JL)+PMEK*DT)/(1.+XLMEK*DT)
          PMC3O3=ISOP(JL)*0.2*RISOPO3(JL)*O3(JL)+ &
     &      OH(JL)* &
     &       (0.31*RISPDOH(JL)*ISOPRD(JL)+0.5*RMTHCOH(JL)*METHAC(JL))+ &
     &      0.5*RMTHCNO3(JL)*METHAC(JL)*ZNO3(JL)+ &
     &      MPAN(JL)*(RMPAN(JL)+RJPAN(JL))
          XLMC3O3=RMC23NO(JL)*ZNO(JL)+RMC23NO2(JL)*ZNO2(JL)+ &
     &      RMC23HO2(JL)*HO2(JL)
          MC3O3(JL)=(MC3O30(JL)+PMC3O3*DT)/(1.+XLMC3O3*DT)

          PCH3O2H=RMO2HO2(JL)*CH3O2(JL)*HO2(JL)
          XLCH3O2H=(ROHPCAT(JL)+ROHPFRM(JL))*OH(JL)+RJMEPE(JL)
          CH3O2H(JL)=(CH3O2H0(JL)+PCH3O2H*DT)/(1.+XLCH3O2H*DT)

          ! LG-  082003 added the production and destruction terms of hydroxy-
          !      methyl/alkylhydroxyperoxide

          PCH2OHO2H=0.105*RMBTERPO3(JL)*MBTERP(JL)*O3(JL)
          XLCH2OHO2H=RCH2OHO2H(JL)  ! LG- HMHP goes into HCOOH and H2O: 
                      ! The decomposition rate is taken from Peter Neebs thesis
          CH2OHO2H(JL)=(CH2OHO2H0(JL)+PCH2OHO2H*DT)/(1.+XLCH2OHO2H*DT)

          ! LG-  Table 2.7, Valverde's Thesis, alpha- and beta pinene destruction to mimic 
          !      differences between endo and exocyclic alkene ozonolyis 

          PRCHOHO2H=0.05*RMATERPO3(JL)*MATERP(JL)*O3(JL)+ &
     &           0.23*RMBTERPO3(JL)*MBTERP(JL)*O3(JL)
          XLRCHOHO2H=RRCHOHO2H(JL) ! LG- all the HAHP goes to H2O2 (see below)
                                   ! LG- correct???? adding production of RCHO
          RCHOHO2H(JL)=(RCHOHO2H0(JL)+PRCHOHO2H*DT)/(1.+XLRCHOHO2H*DT)

          PHNO3=RNO2OH(JL)*ZNO2(JL)*OH(JL)+ &
     &      2.*(RN2O5AQ(JL)+RN2O5L(JL))*ZN2O5(JL)+ &
     &      ZNO3(JL)* &
     &        (RALD2NO3(JL)*ALD2(JL)+RFRMNO3(JL)*CH2O(JL)+ &
     &        RMTHCNO3(JL)*METHAC(JL)+0.8*RISOPNO3(JL)*ISOP(JL)+ &
     &        RISPDNO3(JL)*ISOPRD(JL))+ &
     &      0.8*RISOPNO2(JL)*ZNO2(JL)*ISOP(JL)

            ! WP- removed isontr:

            !     &   +RBXO2NNO(JL)*BXO2N(JL)*ZNO(JL)
        
          XLHNO3=RJHNO3(JL)+OH(JL)*ROHHNO3(JL)
          HNO3(JL)=(HNO30(JL)+PHNO3*DT)/(1.+XLHNO3*DT)
          PH2O2=RHO2HO2(JL)*HO2(JL)*HO2(JL) &

            ! LG-  added the production of H2O2 through the ozonolysis of the terpenes/
            !      alkenes. The key production term is the further destruction of
            !      RCHOHO2H

     &      +RRCHOHO2H(JL)*RCHOHO2H(JL) ! LG- assuming a 100% conversion of HAHP to H2O2

          XLH2O2=RJH2O2(JL)+RHPOH(JL)*OH(JL)
          H2O2(JL)=(H2O20(JL)+PH2O2*DT)/(1.+XLH2O2*DT)
          PCH2O=RMO2NO(JL)*CH3O2(JL)*ZNO(JL)+RJMEPE(JL)*CH3O2H(JL)+ &
     &      ROHPFRM(JL)*CH3O2H(JL)*OH(JL)+RJALD2(JL)*ALD2(JL)+ &
     &      1.33*RMO2MO2(JL)*CH3O2(JL)*CH3O2(JL)+ &
     &      C2O3(JL)*(RC23NO(JL)*ZNO(JL)+2*RC23C23(JL)*C2O3(JL)+ &

            ! LG- removed acetic acid: removal of reaction 53 (see paper Duncan)

            !     &     0.79*RC23HO2(JL)*HO2(JL)+ &

            ! LG- added acetic acid: updating rc23mo2a and b

     &      RC23MO2A(JL)*CH3O2(JL)+RC23MO2B(JL)*CH3O2(JL))+ &

     &      OH(JL)* &
     &       (1.56*RETHOH(JL)*ETH(JL)+ROLEOH(JL)*OLE(JL)+ &
     &        0.63*RISOPOH(JL)*ISOP(JL)+0.14*RISPDOH(JL)*ISOPRD(JL)+ &
     &        0.08*RMTHCOH(JL)*METHAC(JL)+0.3*RMVKOH(JL)*MVK(JL)+ &
     &        0.4*RMPANOH(JL)*MPAN(JL)+0.5*RMEKOH(JL)*MEK(JL))+ &
     &      O3(JL)* &
     &       (RETHO3(JL)*ETH(JL)+0.74*ROLEO3(JL)*OLE(JL)+ &
     &        0.6*RISOPO3(JL)*ISOP(JL)+0.39*RISPDO3(JL)*ISOPRD(JL)+ &
     &        0.7*RMPANO3(JL)*MPAN(JL)+0.2*RMTHCO3(JL)*METHAC(JL)+ &
     &        0.1*RMVKO3(JL)*MVK(JL))+ &
     &      ZNO3(JL)* &
     &        (ROLENO3(JL)*OLE(JL)+0.33*RISPDNO3(JL)*ISOPRD(JL))+ &
     &      MC3O3(JL)*(RMC23NO(JL)*ZNO(JL)+2*RMC23HO2(JL)*HO2(JL))
            XLCH2O=RJACH2O(JL)+RJBCH2O(JL)+OH(JL)*RFRMOH(JL)+ &
     &        ZNO3(JL)*RFRMNO3(JL)
          CH2O(JL)=(CH2O0(JL)+PCH2O*DT)/(1.+XLCH2O*DT)
          PPAN=ZNO2(JL)*RC23NO2(JL)*C2O3(JL)+ &
     &     MPAN(JL)*(0.3*RMPANO3(JL)*O3(JL)+0.4*RMPANOH(JL)*OH(JL))
          XLPAN=RPAN(JL)+RJPAN(JL)
          PAN(JL)=(PAN0(JL)+PPAN*DT)/(1.+XLPAN*DT)
          PMPAN=RMC23NO2(JL)*MC3O3(JL)*ZNO2(JL)
          XLMPAN=RMPAN(JL)+RJPAN(JL)+RMPANO3(JL)*O3(JL)+ &
     &      RMPANOH(JL)*OH(JL)
          MPAN(JL)=(MPAN0(JL)+PMPAN*DT)/(1.+XLMPAN*DT)
          PETH=0.
          XLETH=RETHOH(JL)*OH(JL)+RETHO3(JL)*O3(JL)
          ETH(JL)=(ETH0(JL)+PETH*DT)/(1.+XLETH*DT)
          PACET=0.64*RRORA(JL)*ROR(JL)+0.6*RMPANOH(JL)*MPAN(JL)*OH(JL)+ &
     &      0.8*RJNITR(JL)*NITR(JL)
          XLACET=RACETOH(JL)*OH(JL)+RJACET(JL)
          ACET(JL)=(ACET0(JL)+PACET*DT)/(1.+XLACET*DT)
          PNTR=RRORNO2(JL)*ROR(JL)*ZNO2(JL)+ &
     &      RXO2NNO(JL)*XO2N(JL)*ZNO(JL)
          XLNTR=RJNITR(JL)
          NITR(JL)=(NITR0(JL)+PNTR*DT)/(1.+XLNTR*DT)

          ! WP- added isontr:

          PISONTR=0.2*RBXO2NNO(JL)*BXO2N(JL)*ZNO(JL)
          XLISONTR=RJNITR(JL)+RISONTROH(JL)*OH(JL)
          ISONTR(JL)=(ISONTR0(JL)+PISONTR*DT)/(1.+XLISONTR*DT)

          ! WP- end

          ! LG- March 2002, added the destruction of the monoterpenes

          PMATERP=0.
          XLMATERP=RMATERPOH(JL)*OH(JL)+RMATERPO3(JL)*O3(JL)+ &
     &      RMATERPNO3(JL)*ZNO3(JL)
          MATERP(JL)=(MATERP0(JL)+PMATERP*DT)/(1.+XLMATERP*DT)
          PMBTERP=0.
          XLMBTERP=RMBTERPOH(JL)*OH(JL)+RMBTERPO3(JL)*O3(JL)+ &
     &      RMBTERPNO3(JL)*ZNO3(JL)
          MBTERP(JL)=(MBTERP0(JL)+PMBTERP*DT)/(1.+XLMBTERP*DT)

          ! mz_lg_20060403+ added extra monoterpene

          PMTTERP=0.
          XLMTTERP=RMTTERPOH(JL)*OH(JL)+RMTTERPO3(JL)*O3(JL)+ &
     &      RMTTERPNO3(JL)*ZNO3(JL)
          MTTERP(JL)=(MTTERP0(JL)+PMTTERP*DT)/(1.+XLMTTERP*DT)

          ! LG- April 2002, added the destruction of the sesquiterpenes

          PSQTERP=0.
          XLSQTERP=RSQTERPOH(JL)*OH(JL)+RSQTERPO3(JL)*O3(JL)+ &
     &      RSQTERPNO3(JL)*ZNO3(JL)
          SQTERP(JL)=(SQTERP0(JL)+PSQTERP*DT)/(1.+XLSQTERP*DT)

          ! mz_lg_20060330+ added the 2- and 1-double bond sesquiterpenes

          PSQTERP2B=RSQTERPO3(JL)*SQTERP(JL)*O3(JL)
          XLSQTERP2B=RSQTERP2BOH(JL)*OH(JL)+RSQTERP2BO3(JL)*O3(JL)+ &
     &      RSQTERP2BNO3(JL)*ZNO3(JL)
          SQTERP2B(JL)=(SQTERP2B0(JL)+PSQTERP2B*DT)/(1.+XLSQTERP2B*DT)

          PSQTERP1B=RSQTERP2BO3(JL)*SQTERP2B(JL)*O3(JL)
          XLSQTERP1B=RSQTERP1BOH(JL)*OH(JL)+RSQTERP1BO3(JL)*O3(JL)+ &
     &      RSQTERP1BNO3(JL)*ZNO3(JL)
          SQTERP1B(JL)=(SQTERP1B0(JL)+PSQTERP1B*DT)/(1.+XLSQTERP1B*DT)

          ! LG- added formic/acetic acid:the formation and destruction of 
          !    formic acid and acetic acid

          PHCOOH=0.56*0.5*RISOPO3(JL)*ISOP(JL)*O3(JL) &
     &     +RCH2OHO2H(JL)*CH2OHO2H(JL)  ! LG- added the production through the
                                        ! decomposition of HMHP
          XLHCOOH=RHCOOHOH(JL)*OH(JL)
          HCOOH(JL)=(HCOOH0(JL)+PHCOOH*DT)/(1.+XLHCOOH*DT)
          PCH3CO2H=RC23HO2(JL)*C2O3(JL)*HO2(JL)+ &
     &           RC23MO2B(JL)*C2O3(JL)*CH3O2(JL)
          XLCH3CO2H=RCH3CO2HOH(JL)*OH(JL)
          CH3CO2H(JL)=(CH3CO2H0(JL)+PCH3CO2H*DT)/(1.+XLCH3CO2H*DT)

          ! ESS_lg_20100628+, added a selection of anthropogenic alkenes

          PHEXANE=0.
          XLHEXANE=RHEXANEOH(JL)*OH(JL)+RHEXANEO3(JL)*O3(JL)+ &
     &      RHEXANENO3(JL)*ZNO3(JL)
          HEXANE(JL)=(HEXANE0(JL)+PHEXANE*DT)/(1.+XLHEXANE*DT)

          PBUTADIENE=0.
          XLBUTADIENE=RBUTADIENEOH(JL)*OH(JL)+RBUTADIENEO3(JL)*O3(JL)+ &
     &      RBUTADIENENO3(JL)*ZNO3(JL)
          BUTADIENE(JL)=(BUTADIENE0(JL)+PBUTADIENE*DT)/(1.+XLBUTADIENE*DT)

          PTMBENZENE=0.
          XLTMBENZENE=RTMBENZENEOH(JL)*OH(JL)+RTMBENZENEO3(JL)*O3(JL)+ &
     &      RTMBENZENENO3(JL)*ZNO3(JL)
          TMBENZENE(JL)=(TMBENZENE0(JL)+PTMBENZENE*DT)/(1.+XLTMBENZENE*DT)

          ! ESS_lg_20100628- 

          CH4(JL)=CH40(JL)/(1.+RCH4OH(JL)*OH(JL)*DT)
          PCO=CH2O(JL)*XLCH2O+RJALD2(JL)*ALD2(JL)+ &
     &      O3(JL)* &
     &       (0.44*RETHO3(JL)*ETH(JL)+0.33*ROLEO3(JL)*OLE(JL)+ &
     &        0.13*RMPANO3(JL)*MPAN(JL))+ &
     &      MGLY(JL)*(OH(JL)*RMGLYOH(JL)+1.5*RJMGLY(JL))
          CO(JL)=(CO0(JL)+PCO*DT)/(1.+RCOOH(JL)*OH(JL)*DT)
          XLALD2=RALD2OH(JL)*OH(JL)+RALD2NO3(JL)*ZNO3(JL)+RJALD2(JL)
          PALD2= &
     &      OH(JL)*(0.22*RETHOH(JL)*ETH(JL)+0.11*RPAROH(JL)*PAR(JL)+ &
     &       ROLEOH(JL)*OLE(JL)+0.7*RMVKOH(JL)*MVK(JL)+ &
     &       RMEKOH(JL)*MEK(JL)+0.19*RISPDOH(JL)*ISOPRD(JL))+ &
     &      O3(JL)*(0.1*RMTHCO3(JL)*METHAC(JL)+0.5*ROLEO3(JL)*OLE(JL)+ &
     &       0.15*RISOPO3(JL)*ISOP(JL)+0.06*RISPDO3(JL)*ISOPRD(JL))+ &
     &      ZNO3(JL)*(ROLENO3(JL)*OLE(JL)+0.8*RISOPNO3(JL)*ISOP(JL)+ &
     &       0.33*RISPDNO3(JL)*ISOPRD(JL))+ &
     &      RJMEK(JL)*MEK(JL)+ &
     &      1.1*RRORA(JL)*ROR(JL)+0.8*RISOPNO2(JL)*ZNO2(JL)*ISOP(JL)+ &
     &      0.2*RJNITR(JL)*NITR(JL)
          ALD2(JL)=(ALD20(JL)+PALD2*DT)/(1.+XLALD2*DT)
          PPAR=0.
          XLPAR=RPAROH(JL)*OH(JL)+RRXPAR(JL)*RXPAR(JL)
          PAR(JL)=(PAR0(JL)+PPAR*DT)/(1.+XLPAR*DT)

          ! LG-  added sulfur chemistry:

          IF (LSULFCHEM) THEN
            XLSO2=ROHSO2(JL)*OH(JL)
            SO2(JL)=(SO20(JL)+(ROHDMS(JL)*OH(JL)+ &
     &          RNO3DMS(JL)*ZNO3(JL))*DMS(JL)*DT)/(1.+XLSO2*DT)
            PSAER=ROHSO2(JL)*OH(JL)*SO2(JL)
            SAER(JL)=SAER0(JL)+PSAER*DT
            XLDMS=ROHDMS(JL)*OH(JL)+RNO3DMS(JL)*ZNO3(JL)
            DMS(JL)=DMS0(JL)/(1.+XLDMS*DT)
          ENDIF

          ! LG-  and ammonia chemistry

          IF (LNH3CHEM) THEN
            ACID(JL)=(ACID0(JL)+(2.*XLSO2*SO2(JL))*DT)/ &
     &          (1.+RSO4NH3(JL)*NH3(JL)*DT)
            NH4(JL)=NH40(JL)+ACID(JL)*RSO4NH3(JL)*NH3(JL)*DT
            PNH2=OH(JL)*ROHNH3(JL)

            ! LG-    note that the destruction of NH3 by the reaction with OH is not
            !        included here. In the paper by Dentener and Crutzen [1994], this
            !        reaction is mentioned being an important source of N2O but also a
            !        very uncertain source due to the uncertainty in the rate constants

            XLNH3=ACID(JL)*RSO4NH3(JL)+PNH2
            NH3(JL)=NH30(JL)/(1.+XLNH3*DT)
            XLNH2=RNONH2(JL)*ZNO(JL)+RNO2NH2(JL)*ZNO2(JL)+ &
     &      RHO2NH2(JL)*HO2(JL)+RO2NH2(JL)+RO3NH2(JL)*O3(JL)
            NH2(JL)=(NH20(JL)+NH3(JL)*PNH2*DT)/(1.+XLNH2*DT)
          ENDIF

          ! LG-  end

          NEWNOY=ZNO(JL)+ZNO2(JL)+ZNO3(JL)+2*ZN2O5(JL)+HNO4(JL)+ &
     &        PAN(JL)+MPAN(JL)+NITR(JL)

          ! LG-  radon decay, Radioactive decay of 222^Ra; life time of 5.5 (or 3.8??) day

          RSP_LTIME=1./(5.5*86400.)
          RADON(JL)=RADON0(JL)*EXP(-PTMST*RSP_LTIME)

          ! LG-  iteration counter

          ITER=ITER+1

          IF (ABS((NEWNOY-OLDNOY)/NEWNOY).LT.1.E-15.OR.ITER.GT.250)  &
     &      GOTO 999
          OLDNOY=NEWNOY

          GOTO 9999

  999   CONTINUE  

 2104 CONTINUE 

      ! **** --- Finish all
      DO 106 JL=1,NLON

        ZPM(JL,JK,idt_NOX)=ZNO(JL)+ZNO2(JL)+ZNO3(JL)+2*ZN2O5(JL)+HNO4(JL)
        DODDN=ZPM(JL,JK,idt_NOX)+PAN(JL)+MPAN(JL)+NITR(JL)-ODDN0(JL)

        ZPM(JL,JK,idt_O3)=O3(JL)
        ZPM(JL,JK,idt_CO)=CO(JL)
        ZPM(JL,JK,idt_HNO3)=ZPM(JL,JK,idt_HNO3)-DODDN
        ZPM(JL,JK,idt_H2O2)=H2O2(JL)
        ZPM(JL,JK,idt_CH3O2H)=CH3O2H(JL)
        ZPM(JL,JK,idt_CH2O)=CH2O(JL)
        ZPM(JL,JK,idt_ALD2)=ALD2(JL)
        ZPM(JL,JK,idt_PAR)=PAR(JL)
        ZPM(JL,JK,idt_OLE)=OLE(JL)
        ZPM(JL,JK,idt_ETH)=ETH(JL)
        ZPM(JL,JK,idt_PAN)=PAN(JL)
        ZPM(JL,JK,idt_ACET)=ACET(JL)
        ZPM(JL,JK,idt_ISOP)=ISOP(JL)
        ZPM(JL,JK,idt_MGLY)=MGLY(JL)
        ZPM(JL,JK,idt_ISOPRD)=ISOPRD(JL)
        ZPM(JL,JK,idt_METHAC)=METHAC(JL)
        ZPM(JL,JK,idt_MVK)=MVK(JL)
        ZPM(JL,JK,idt_MEK)=MEK(JL)
        ZPM(JL,JK,idt_MPAN)=MPAN(JL)
        ZPM(JL,JK,idt_NTR)=NITR(JL)
        ZPM(JL,JK,idt_DMS)=DMS(JL)
        ZPM(JL,JK,idt_SO2)=SO2(JL)
        ZPM(JL,JK,idt_SO4)=SAER(JL)
        ZPM(JL,JK,idt_RAD)=RADON(JL)
        ZPM(JL,JK,idt_ISONTR)=ISONTR(JL)
        ZPM(JL,JK,idt_HCOOH)=HCOOH(JL)
        ZPM(JL,JK,idt_CH3CO2H)=CH3CO2H(JL)
        ZPM(JL,JK,idt_NH2)=NH2(JL)
        ZPM(JL,JK,idt_NH3)=NH3(JL)
        ZPM(JL,JK,idt_NH4)=NH4(JL)
        ZPM(JL,JK,idt_APIN)=MATERP(JL)
        ZPM(JL,JK,idt_BPIN)=MBTERP(JL)
        ZPM(JL,JK,idt_MTTERP)=MTTERP(JL) 
        ZPM(JL,JK,idt_SQTERP)=SQTERP(JL)
        ZPM(JL,JK,idt_SQTERP2B)=SQTERP2B(JL) 
        ZPM(JL,JK,idt_SQTERP1B)=SQTERP1B(JL)
        ZPM(JL,JK,idt_HONO)=HONO(JL)
        ZPM(JL,JK,idt_CH2OHO2H)=CH2OHO2H(JL)
        ZPM(JL,JK,idt_RCHOHO2H)=RCHOHO2H(JL)
        ZPM(JL,JK,idt_HEXANE)=HEXANE(JL)
        ZPM(JL,JK,idt_BUTADIENE)=BUTADIENE(JL) 
        ZPM(JL,JK,idt_TMBENZENE)=TMBENZENE(JL)
        ZPM(JL,JK,idt_NO)=ZNO(JL)
        ZPM(JL,JK,idt_NO2)=ZNO2(JL)
        ZPM(JL,JK,idt_NO3)=ZNO3(JL)
        ZPM(JL,JK,idt_N2O5)=ZN2O5(JL)
        ZPM(JL,JK,idt_HNO4)=HNO4(JL)
        ZPM(JL,JK,idt_OH)=OH(JL)
        ZPM(JL,JK,idt_HO2)=HO2(JL)
        ZPM(JL,JK,idt_CH3O2)=CH3O2(JL)
        ZPM(JL,JK,idt_C2O3)=C2O3(JL)
        ZPM(JL,JK,idt_XO2)=XO2(JL)
        ZPM(JL,JK,idt_ROR)=ROR(JL)
        ZPM(JL,JK,idt_XO2N)=XO2N(JL)
        ZPM(JL,JK,idt_RXPAR)=RXPAR(JL)
        ZPM(JL,JK,idt_BXO2N)=BXO2N(JL)          
        ZPM(JL,JK,idt_MC3O3)=MC3O3(JL)

        ! ESS_lg_2013013+ added OH reactivity

        OHreact(istep,JK)=(                            &
               RHNO4OH(JL)*ZPM(JL,JK,idt_hno4)+        &
               RHO2OH(JL)*ZPM(JL,JK,idt_ho2)+          &
               RNO2OH(JL)*ZPM(JL,JK,idt_no2)+          &
               ROHHNO3(JL)*ZPM(JL,JK,idt_hno3)+        &
               RCOOH(JL)*ZPM(JL,JK,idt_co)+            &
               RO3OH(JL)*ZPM(JL,JK,idt_o3)+            &
               RHPOH(JL)*ZPM(JL,JK,idt_h2o2)+          &
               RFRMOH(JL)*ZPM(JL,JK,idt_ch2o)+         &
               RH2OH(JL)*ZPRHOA(JL,JK)*6.022045E23/ZMAIR+ &
               RETHOH(JL)*ZPM(JL,JK,idt_eth)+          &
               ROLEOH(JL)*ZPM(JL,JK,idt_ole)+          &
               0.91*RISOPOH(JL)*ZPM(JL,JK,idt_isop)+   &
               0.68*RISPDOH(JL)*ZPM(JL,JK,idt_isoprd)+ &
               0.5*RMTHCOH(JL)*ZPM(JL,JK,idt_methac)+  &
               0.3*RMVKOH(JL)*ZPM(JL,JK,idt_mvk)+      &
               0.11*RPAROH(JL)*ZPM(JL,JK,idt_par)+     &
               0.4*RMPANOH(JL)*ZPM(JL,JK,idt_mpan)+    &
               RCH4OH(JL)*ZPM(JL,JK,idt_ch4)+          &
               ROHPCAT(JL)*ZPM(JL,JK,idt_ch3o2h)+      &
               ROHOH(JL)*ZPM(JL,JK,idt_oh)+            &
               0.89*RPAROH(JL)*ZPM(JL,JK,idt_par)+     &
               RMGLYOH(JL)*ZPM(JL,JK,idt_mgly)+        &
               0.09*RISOPOH(JL)*ZPM(JL,JK,idt_isop)+   &
               0.32*RISPDOH(JL)*ZPM(JL,JK,idt_isoprd)+ &
               0.5*RMTHCOH(JL)*ZPM(JL,JK,idt_methac)+  &
               RMEKOH(JL)*ZPM(JL,JK,idt_mek)+          &
               0.7*RMVKOH(JL)*ZPM(JL,JK,idt_mvk)+      &
               RALD2OH(JL)*ZPM(JL,JK,idt_ald2)+        & 
               0.6*RMPANOH(JL)*ZPM(JL,JK,idt_mpan)+    &
               RACETOH(JL)*ZPM(JL,JK,idt_acet)+        &
               RMATERPOH(JL)*ZPM(JL,JK,idt_apin)+      &
               RMBTERPOH(JL)*ZPM(JL,JK,idt_bpin)+      &
               RMTTERPOH(JL)*ZPM(JL,JK,idt_mtterp)+    &
               RSQTERPOH(JL)*ZPM(JL,JK,idt_sqterp)+    &
               RSQTERP2BOH(JL)*ZPM(JL,JK,idt_sqterp2b)+&
               RSQTERP1BOH(JL)*ZPM(JL,JK,idt_sqterp1b)+&
               RHEXANEOH(JL)*ZPM(JL,JK,idt_hexane)+    &  ! ESS_lg_20100628+ &added some anthropogenic alkenes
               RBUTADIENEOH(JL)*ZPM(JL,JK,idt_butadiene)+ &
               RTMBENZENEOH(JL)*ZPM(JL,JK,idt_tmbenzene)+ &! ESS_lg_20100628- 
               RNOOH(JL)*ZPM(JL,JK,idt_no))

        ! ESS_lg_2013013-

  106 CONTINUE

  101 CONTINUE

      DO 199 JK=1,nveglay+1
      DO 199 JL=1,NLON
      ! LG- assigning the tracer concentrations 

      DO JT=1,NTRAC
        IF (JK.EQ.1) THEN 
          pxtm1(istep,JT)=ZPM(JL,JK,JT)
        ELSE
          pxtmveg(JK-1,JT)=ZPM(JL,JK,JT)
        ENDIF
      ENDDO

  199 CONTINUE

  END SUBROUTINE emdep_xtsurf_chem
  ! ESS_lg_20120721-

  ! mz_lg_20050719+ implementation of subroutine to calculate wind profile
  !     in the canopy
  !=============================================================================

  SUBROUTINE emdep_xtsurf_wndprof( nstep, & ! ESS_lg_20120722+
    nveglay, ustar, lai, hc, forestfr, z0m, disp, u_veg)

    ! ---------------------------------------------------------------
    !     Subroutine wndprof determines the mean wind speed profile 
    !     within the canopy using a form of the exponential wind profile 
    !     (cionco). The code is taken from the DDIM model (Dry deposition 
    !     Inferential Model) 
    !
    !     Laurens Ganzeveld, 1998, modified for implementation
    !     in echam5, July 2005
    ! --------------------------------------------------------------------
    ! Interface:
    ! ----------
    ! input 
    ! nstep     : # of timesteps ! ESS_lg_20120722+
    ! nveglay   : no. of canopy layers
    ! ustar     : friction velocity [m s-1]
    ! lai       : leaf area index [m2 m-2]
    ! hc        : canopy height [m]
    ! forestfr  : forest fraction [0-1]
    ! z0m       : surface roughness [m]
    ! disp      : displacement height [m]
    !
    ! output

    IMPLICIT NONE 

    ! I/O
    INTEGER,  INTENT(in)  :: nstep ! ESS_lg_20120722+
    INTEGER,  INTENT(in)  :: nveglay
    REAL(dp), INTENT(in)  :: ustar(:), lai(:), hc(:),   &
                             forestfr(:), z0m(:), disp(:)
    REAL(dp), INTENT(out) :: u_veg(:,:)

    ! local parameters
    INTEGER :: klon, i, ii, k, jl, ptype
    REAL(dp):: beta, alpha, d0, z0

    ! INITIALIZATION
    klon=nstep ! SIZE(ustar) ! ESS_lg_20120722+

    DO jl=1,klon
 
      u_veg(jl,:)=1.e-2_dp*ustar(jl) ! minimum windspeed to avoid
                                     ! dividing by zero (resistances)

      IF (hc(jl) > hcmin .AND. lai(jl) > 1.e-10_dp) THEN

        ! LG-   there are two different LAD profiles distinguished, one for
        !       forest canopies and one for non-forested canopies, these 
        !       are assigned based on the fraction of forested area in each
        !       grid square. For the forests, PTYPE=2 whereas for other canopies
        !       PTYPE=4 (August 2000) (this should be consistent with the
        !       preprocessed LAD profile fields!)

        ptype=2
      
        IF (forestfr(jl) < 0.5_dp) ptype=4     
      
        !
        !       ************************
        !       * SET ALPHA AND BETA   *
        !       ************************
        !

        beta = 1._dp-(ptype-1._dp)*.25_dp
        IF (hc(jl) > 10._dp) THEN
          alpha = lai(jl)
          IF (alpha > 4.0_dp) alpha = 4.0_dp
        ELSE
          alpha = .65_dp*lai(jl)
          IF (alpha > 3.0_dp) alpha = 3.0_dp
        ENDIF

        !
        !      ******************************************************
        !      *  COMPUTATION OF D0 AND Z0 AS A FUNCTION            *
        !      *  OF LAI BASED ON MODEL COMPUTATIONS OF             *
        !      *  SHAW AND PERIERA (1982) AGRICULTURAL METEOROLOGY  *
        !      ******************************************************
        !

        d0 = hc(jl)*(.05_dp+0.5_dp*lai(jl)**0.20_dp +  &
           (ptype-1._dp)*.05_dp)

        IF (lai(jl) < 1._dp) THEN
           z0 = hc(jl)*0.1_dp
        ELSE
           z0 = hc(jl)*(0.23_dp - 0.1_dp*lai(jl)**0.25_dp - &
              (ptype-1._dp)*.015_dp)
        ENDIF

        IF (disp(jl) > 0._dp) d0=disp(jl)
        IF (z0m(jl) > 0._dp) z0=z0m(jl)

        ! mz_lg_20050725+ note the change in layer index: the top
        ! layer is not indicated by index number 1: the soil-canopy
        ! layer is indicated by nveglay
        u_veg(jl,1) = 2.5_dp*ustar(jl)*LOG((hc(jl)-d0+z0)/z0)

        !
        !     *****************************************
        !     *  ESTIMATE WITHIN CANOPY WIND PROFILE  *
        !     *  USING A FORM OF CIONCO MODEL         *
        !     *****************************************
        !
        DO i=1,nveglay
           ! mz_lg_20050725+ note the change in layer index:
           k=nveglay+1-i
           u_veg(jl,i)=u_veg(jl,1)*EXP(-alpha*((1._dp-FLOAT(k)/ &
              FLOAT(nveglay))**beta))
        END DO

      ENDIF ! IF (hc(jl) > hcmin

    END DO ! end DO jk=1,klon

  END SUBROUTINE emdep_xtsurf_wndprof
  ! mz_lg_20050719-

  ! mz_lg_20120722+ implementation of subroutine to read in a file that contains observations
  ! =========================================================================================

  SUBROUTINE emdep_xtsurf_readdata(infilename,nstop,nprint,dtime,ndtgact,ldatltime,readdata)  ! ESS_lg_20130228+ ndtgact
     ! -------------------------------------------------------------------------------------
     !  This subroutine reads in a file containing observed parameter
     !  values, the format of reading needs to be adapted to the specific
     !  file structure. An example of the format of the file for which this      
     !  code has been written can be found in the directory: 
     !  input/data/   
     !
     !  Written by Laurens Ganzeveld, 26-01-99, updated for use in box model, 22-07-2012
     ! -------------------------------------------------------------------------------------

     IMPLICIT NONE

     INTEGER,  INTENT(in)   :: nstop,nprint,ndtgact(:) ! ESS_lg_20130228+ 

     REAL(dp), INTENT(in)   :: dtime

     CHARACTER(LEN=75), INTENT(in) :: infilename
     CHARACTER(LEN=20), INTENT(in) :: ldatltime(:) ! ESS_lg_20130415+
	 
     INTEGER I,II,IY,N,J,NDATA,NDATA_MAX,NPARAM,NPARAM_MAX,NSTEP,     &
             NDATE,NODAY,NYMDIN,NTBASE,NDIN,NYEAR,NMONTH,NDAY,DAY1,   &
             NYLEN,NCBASE,IDAT2C,JDAY_OBS,NOBS,IP,IIP,NDTG,NDTG_INIT, & ! ESS_lg_20101216+
             NUNDATA,NUNOBS

     PARAMETER (NDATA_MAX=100000,NPARAM_MAX=30,NUNDATA=2,NUNOBS=3) ! ESS_lg_20100625+ increased to deal with larger files

     REAL STARTTIME, ENDTIME, TIME(NDATA_MAX), HOUR, MIN, SEC, DT, DTIME_INT_SEC, DTIME_OBS    

     REAL(dp) :: OBSERV(NDATA_MAX,NPARAM_MAX)
     REAL(dp) :: OBSERVM(NDATA_MAX,NPARAM_MAX)
     REAL(dp) :: READDATA(:,:)

     CHARACTER*250 DUMMY
     CHARACTER*12 PARNAME(0:NPARAM_MAX)

     INTEGER IDMAX(12)
     DATA IDMAX/31,28,31,30,31,30,31,31,30,31,30,31/

     ! ------------------PROGRAM STARTS HERE----------------------
     DO NSTEP=1,NSTOP

       IF (NSTEP.EQ.1) THEN
         WRITE(*,'(2a)') &
           ' Start reading file with observations: ',infilename

         OPEN(NUNDATA,FILE=infilename,STATUS='unknown')
         READ(NUNDATA,'(a100)') DUMMY

         ! LG-   reading the begin and end time of the observations

         READ(NUNDATA,*) NDATE,STARTTIME,ENDTIME,NODAY,DTIME_OBS

         ! LG-   calculation of some time parameters of the observations

         NYMDIN = NDATE/100
         NTBASE = NDATE - 100*NYMDIN
         NYEAR  = NYMDIN/10000
         NMONTH =(NYMDIN-NYEAR*10000)/100
         NDAY   = NYMDIN-NYEAR*10000-NMONTH*100
         NDTG   = NDTGACT(1) ! ESS_lg_20130228+ 
         NDTG_INIT=NDATE ! ESS_lg_20120102+ removed STARTTIME ! ESS_lg_20101216+ 

         DAY1=0
         DO 101 I=1,NMONTH
           IF (I.GT.1) DAY1=DAY1+IDMAX(I-1)
   101   CONTINUE
         JDAY_OBS=DAY1+NDAY

         ! ESS_lg_20120112+ introduced to avoid using the first data points for NSTEP=0 for 
         !       a STARTTIME < the NDTG_INIT
         IF (NDTG.GT.NDTG_INIT) THEN  ! ESS_lg_20130228+ 
           WRITE(*,'(1a)') &
             ' The start time of the read-in observations is before the start time of the simulation'
           WRITE(*,'(1a)') &
             ' To avoid an unwanted shift in the timing of read-in observations; modify the init. NDTG'
           WRITE(*,'(1a)') &
             ' STOP called in emdep_xtsurf_readdata'
           STOP
         ENDIF
         ! ESS_lg_20120112-

         ! LG-   reading the header with the names of the parameters

         READ(NUNDATA,'(A250)') DUMMY
         READ(NUNDATA,'(A250)') DUMMY

         N=0
         IIP=0
   77    CONTINUE        
         IP=INDEX(DUMMY(IIP+1:200),' ')
         PARNAME(N)=DUMMY(IIP+1:IIP+IP-1)
         IF (DUMMY(IIP+1:IIP+IP+1).EQ.' ') GOTO 88
         IIP=IIP+IP
         N=N+1
         GOTO 77
   88    CONTINUE

         NPARAM=N-1

         WRITE(*,'(1a,i3,1a)') &
           ' The read datafile contains data of: ',NPARAM,' parameters'
         DO N=1,NPARAM
           WRITE(*,'(1a,i2,1x,1a,1x,1a)') ' No: ',N,'Name:',PARNAME(N)
         ENDDO
         WRITE(*,'(1A,A10)'),' as a function of ',PARNAME(0)

         WRITE(*,'(1a,4i4)') &
           ' The observations represent the year/month/day/Julian day ', &
             NYEAR,NMONTH,NDAY,JDAY_OBS

         WRITE(*,'(1a,f6.2,1a,f6.2,1a,i3,1a,f6.0,1a)') &
           ' The observations start at ',STARTTIME,' GMT and end at ', &
             ENDTIME,' GMT after ',NODAY,' day(s), timestep: ',DTIME_OBS,' [s]'

         WRITE(*,'(1A)')' Press enter to continue'
         READ (*,*)

         I=1

         ! LG-   opening of file for checking the average observations 

         OPEN(UNIT=NUNOBS,FILE='output/avg_obs.out', &
           FORM='FORMATTED',STATUS='UNKNOWN')
         WRITE(NUNOBS,'(a10,a10,a16,100(1x,a12))') &
           'nstep','time [s]','ldatltime',(PARNAME(N),N=1,NPARAM)

     100 CONTINUE

         READ(NUNDATA,*,ERR=200,END=300) &
            TIME(I),(OBSERV(I,N),N=1,NPARAM)
        
! ESS_lg_20120921+ to check reading-in data uncomment the following print and read statement
!         print *,'reading in observations: ',I,TIME(i),(OBSERV(I,N),N=1,NPARAM) 
!         read (*,*)

         I=I+1  ! ESS_lg_20100625+ 

         GOTO 100      

   200   PRINT *,'error reading data file, check file'  ! in case of error
         STOP

   300   NDATA=I-1                            ! IN CASE OF EOF
         II=1

         CLOSE(NUNDATA,STATUS='keep')

       ENDIF

       ! LG- calculation of the time in second since the start of the 
       !     simulation in order to apply the proper observations, which
       !     should be defined in default format in seconds 

       DTIME_INT_SEC=NSTEP*DTIME

       ! LG- 200301, in case of a database with the starting time being later 
       !     then the initial reference time of the simulation then the time array of the 
       !     observations is contineously corrected for the evolving time of the simulation by
       !     adding the term DTIME such that when the first data point of the observation dataset
       !     is reached that then the parameters DTIME_INT_SEC and TIME have similar value and where
       !     the time difference between observation(I) and observation(I+1) are used for the interpolation 

       IF (NDTGACT(NSTEP).LT.NDTG_INIT) THEN ! ESS_lg_20130228+ ! ESS_lg_20101216+, further modified, ESS_lg_20100625+ GT instead of GE, check if this
                                             ! also works properly for a dT of 3600s or longer
         DO I=1,NDATA
           TIME(I)=TIME(I)+DTIME
         ENDDO
       ENDIF

       ! LG- resetting the model applied value
       DO N=1,NPARAM
         OBSERVM(NSTEP,N)=-9999.999 ! ESS_lg_20120726+
       ENDDO

       ! LG- calculation of the timestep average parameter value

       ! LG- 200301, in case of an initial time of available input data that
       !     is larger then the actual time in the integration, the 
       !     interpolation is not being done.

       IF (NDTGACT(NSTEP).LT.NDTG_INIT) GOTO 402 ! ESS_lg_20130228+ ! ESS_lg_20101216+, further modified, ESS_lg_20100625+ GT instead of GE, check if this
                                                 ! also works properly for a dT of 3600s or longer

	   DO 401 N=1,NPARAM

         IF (II.GT.NDATA) GOTO 400

         IF (DTIME_INT_SEC.GE.TIME(II)) THEN
    	   DTIME_OBS=TIME(II+1)-TIME(II)
           IF (OBSERV(II,N).LT.-9999..OR.OBSERV(II+1,N).LT.-9999.) & ! ESS_lg_20110309+ LT -9999. NAs!
             GOTO 399                                       
		   IF (II+1.LE.NDATA.AND.DTIME_OBS.GT.0.) THEN
 	  	     OBSERVM(NSTEP,N)=OBSERV(II,N)* & ! ESS_lg_20120726+
               (1.-((DTIME_INT_SEC-TIME(II))/DTIME_OBS))+ &
               OBSERV(II+1,N)*(1.-((TIME(II+1)-DTIME_INT_SEC)/DTIME_OBS))
		   ELSE
			   OBSERVM(NSTEP,N)=OBSERV(II,N) ! ESS_lg_20120726+
           ENDIF
   399     CONTINUE
           IF (DTIME_INT_SEC+DTIME.GE.TIME(II+1).AND.N.EQ.NPARAM) II=II+1
         ENDIF

   400   CONTINUE
   401 CONTINUE
   402 CONTINUE

       IF (NSTEP.EQ.0.OR.MOD(NSTEP,NPRINT).EQ.0)  &
         WRITE(NUNOBS,'(i10,1x,f9.0,2x,A14,100(1x,f12.4))')   &
     &      NSTEP,DTIME_INT_SEC,LDATLTIME(NSTEP),(OBSERVM(NSTEP,N),N=1,NPARAM) ! ESS_lg_20120726+

       DO N=1,NPARAM
         READDATA(NSTEP,N)=OBSERVM(NSTEP,N)
       ENDDO

    ENDDO
    CLOSE(NUNOBS)

  END SUBROUTINE emdep_xtsurf_readdata
  ! ESS_lg_20120722- 

  ! mz_lg_20050721+ implementation of subroutine to calculate photolysis rates
  !     in canopy
  !=============================================================================

  SUBROUTINE emdep_xtsurf_rjveg( nstep, & ! ESS_lg_20120722+
    nveglay, nveglay_hr, latmbios_photo, lai, hc, rvdsl, rvd, fsl, &
    rj, rj_veg)

    ! ----------------------------------------------------------
    !     Calculation of vertical profiles within the canopy
    !     of the photodissociation rates using the vertical 
    !     profiles of direct and diffusive PAR calculated in the
    !     subroutines xtsurf_calcprof   
    !
    !     Laurens Ganzeveld, 1998, modified for implementation
    !     in echam5, July, 2005
    ! --------------------------------------------------------------------
    ! Interface:
    ! ----------
    ! input 
    ! nstep     : # of timesteps ! ESS_lg_20120722+
    ! nveglay   : no. of canopy layers
    ! nveglay_hr: no. of canopy layers; high resolution ! mz_lg_20050721+
    ! latmbios_photo : switch indicating tracer photolysis rates
    ! lai       : leaf area index [m2 m-2]
    ! hc        : canopy height [m]
    ! rvdsl     : diffusive irradiance in surface layer [W m-2] ! mz_lg_20050721+
    ! rvd       : diffusive irradiance in canopy [W m-2]
    ! fsl       : fraction of sunlit leaves [0-1]
    ! rj        : photolysis rates [s-1]
    !
    ! output
    ! rj_veg    : photolysis rates in canopy [s-1]

    IMPLICIT NONE 

    ! I/O
    INTEGER,  INTENT(in)  :: nstep ! ESS_lg_20120722+
    INTEGER,  INTENT(in)  :: nveglay, nveglay_hr
    LOGICAL,  INTENT(in)  :: latmbios_photo(:)
    REAL(dp), INTENT(in)  :: lai(:), hc(:),     &
                             rvdsl(:), rvd(:,:), fsl(:,:), rj(:,:)  
    REAL(dp), INTENT(out) :: rj_veg(:,:,:)

    ! local parameters
    INTEGER :: nphoto, klon, jl, jt, jk, jjk 
    REAL(dp):: fsl_avg(nveglay),rvd_avg(nveglay)

    ! INITIALIZATION
    klon=nstep ! SIZE(rvdsl) ! ESS_lg_20120722+
    rj_veg=0._dp
    nphoto=SIZE(latmbios_photo)

    DO jt=1,nphoto
      rj_veg(:,:,jt)=0._dp    ! mz_lg_20050721+ initialization
      IF (latmbios_photo(jt)) THEN
        DO jl=1,klon
          IF (hc(jl) > hcmin .AND. lai(jl) > 1.e-10_dp) THEN
            DO jk=1,nveglay

              ! mz_lg_20050725+ summing the parameters and determining the
              ! the average of the low resolution fields

              fsl_avg(jk)=0._dp
              rvd_avg(jk)=0._dp
              DO jjk=(jk-1)*nveglay_hr/nveglay+1,jk*nveglay_hr/nveglay
                fsl_avg(jk)=fsl_avg(jk)+fsl(jl,jjk)
                rvd_avg(jk)=rvd_avg(jk)+rvd(jl,jjk)
              END DO
              fsl_avg(jk)=fsl_avg(jk)/(nveglay_hr/nveglay)
              rvd_avg(jk)=rvd_avg(jk)/(nveglay_hr/nveglay)

              rj_veg(jl,jk,jt)= &
                 MAX(0._dp,((1.-fsl_avg(jk))*rvd_avg(jk)/ &
                 MAX(1.e-5_dp,rvdsl(jl))+fsl_avg(jk))*rj(jl,jt))

            END DO ! DO jk=1,nveglay
          END IF
        END DO ! end DO jk=1,klon
      END IF ! IF (latmbios_photo(jt)) THEN
    END DO ! DO jt=1,nphoto

  END SUBROUTINE emdep_xtsurf_rjveg
  ! mz_lg_20050719-

  ! ESS_lg_20130503+ subroutine AGS to calculate the stomatal resistance as a function of CO2 and 
  !   water parameters, etc.
  !=============================================================================
  
  SUBROUTINE emdep_xtsurf_ags(nstep, nveglay, Agstype,     & ! Ags vegetation type, C3; Agstype=1, C4; Agstype=2, conifereous forest; Agstype=3, tropical forest; Agstype=4
                              PTSM1M, PVLTM, PAPHM1,       & ! Radiation, LAI, pressure
                              PWSM1M, PWSMAX,PXTMCO2,      &
                              PQM1, PQS,   PQSAM,          & ! PQSAM, saturation point at the leaf level	 
	                          PSRFL, RBVD, RVD, FSL,       &
							  PRS0,  PRS0_AGSML,           &
							  RMESCO2, RCUTCO2)	  
      IMPLICIT NONE 																	
 
      INTEGER,  INTENT(in)  :: nstep 
      INTEGER,  INTENT(in)  :: nveglay, Agstype

      REAL(dp), INTENT(in)  :: PTSM1M(:),PVLTM(:),PAPHM1(:),       &
                               PWSM1M(:), PWSMAX(:), PXTMCO2(:),   &
							   PQM1(:), PQS(:),PQSAM(:),           &
							   PSRFL(:),RBVD(:),RVD(:,:),FSL(:,:)

	  REAL(dp), INTENT(out) :: PRS0_AGSML(:,:),PRS0(:),RMESCO2(:),RCUTCO2(:)

	  !*    LOCAL STORAGE
      !     ----- -------
      !

      INTEGER klon,klevs
      INTEGER jl, jk, j, k, MCATY, MSOTY
      PARAMETER (MCATY=4,MSOTY=1,KLEVS=1) ! ESS_lg_20150623+ MSOTY was 12, reduced now to 1	  

      REAL(dp) ZDSMAX,ZA1

      REAL(dp) ZGAMMA, ZGM, ZWROOT,ZSOIL,ZFMIN,ZDS,ZF, &
           ZCI,ZEPSILON,ZAMAX,ZAM,ZRD,ZPAR,ZAN,    &
           ZCMIN,ZAMIN,ZGSC,ZGS,ZGSCS,             &
           ZRO,ZCS,ZB,ZDSTAR

	  REAL(dp) VTMPC1,RD,RV,TMELT,CVRAD  ! ESS_lg_20130516+ extra declarations

      REAL(dp) ZEPSR,ZEPSW, E1

      REAL(dp)  CGAMMA25(MCATY),CGAMMAQ10(MCATY),CGM25(MCATY), &
            CGMQ10(MCATY),CGMT1(MCATY),CGMT2(MCATY),       &
            CGC(MCATY),CF0(MCATY),CEPSILON0(MCATY),        &
            CAMAX25(MCATY),CAMAXQ10(MCATY),CAMAXT1(MCATY), &
            CAMAXT2(MCATY),CRICO(MCATY),                   &
			CWPWPT(MSOTY),CQWEVAPT(MSOTY)
      REAL(dp)  PWSAM1M(KLEVS), CVROOTSA(KLEVS) ! ESS_lg_20130516+ temporarily local parameters for multi-layer soil properties   
						
      REAL(dp) C25,C10,CDOT3,C9,C1DOT6,CMA,CMV,CMC02,C1000
      REAL(dp) CCS,CEXPAR

	  ! INITIALIZATION
      ! mz_lg_20040423+ modified
      klon=nstep ! SIZE(loland) ! ESS_lg_20120722+

      ZEPSR=1.E-10 
      ZEPSW=1.E-3 

      ! ESS_lg_20130516+ normally done in the 1-D model system in iniags.f

      ! LG- 05-2001, added the vegetation type number 3 to represent a special
      !     class for coniferous trees in terms of some applied constants for the
      !     physiological model

      CGAMMA25(1)=45.
      CGAMMA25(2)=2.8
      CGAMMA25(3)=45.
      CGAMMA25(4)=45. ! ESS_lg_20150624+ tropical rainforest, taking similar values as C3 vegetation
	                  ! except of specific modifications given below
	  
      CGAMMAQ10(1)=1.5
      CGAMMAQ10(2)=1.5
      CGAMMAQ10(3)=1.5
      CGAMMAQ10(4)=1.5! ESS_lg_20150624+ tropical rainforest, taking similar values as C3 vegetation
	                  ! except of specific modifications given below
	  
      C25=25.
      C10=10.

      CGM25(1)=7.0      ! normally 7.0, reset to smaller value to check the Rstom
                        ! Reinder Ronda has proposed value of 3.8 for tropical 
			            ! forest !
      CGM25(2)=17.5
      CGM25(3)=7.0      ! ESS_lg_20110330+ see Steeneveld et al, 2002, table 3.5, gm=3.5, was 7.0,
      CGM25(4)=7.0      ! ESS_lg_20150624+ tropical rainforest, taking similar values as C3 vegetation
	                    ! except of specific modifications given below
	  
      CGMQ10(1)=2.0
      CGMQ10(2)=2.0
      CGMQ10(3)=2.0
      CGMQ10(4)=1.25!   ! ESS_GK_20150407 Lowered from 2.0 to 1.25 to simulate a more tropical forest plant behavior
	  
      CDOT3=.3

      CGMT1(1)=5.
      CGMT1(2)=13.
      CGMT1(3)=5.
      CGMT1(4)=10.    ! ESS_lg_20150624+ tropical rainforest, see BSc thesis study Gijs Koetsenruijter, July 2015
	  
      CGMT2(1)=28.
      CGMT2(2)=36.
      CGMT2(3)=28.
      CGMT2(4)=38.    ! ESS_lg_20150624+ tropical rainforest, see BSc thesis study Gijs Koetsenruijter, July 2015
	  
      CGC(1)=.25	  ! gmin,c, see report by van de Kassteele
      CGC(2)=.25
      CGC(3)=.25	  ! 1., ESS_lg_20110325+ the value for coniferous forest was 1 but this resulted in way too
                      ! large resistance and small CO2 flux for Hyytiala ! 
					  ! a value of 0 gives too small uptake rates for tropical forest canopies !!!
      CGC(4)=.25      ! ESS_lg_20150624+ tropical rainforest, taking similar values as C3 vegetation
	                  ! except of specific modifications given below

	  CF0(1)=.89	  
      CF0(2)=.85
      CF0(3)=.90      ! ESS_lg_20110530+ 0.6 based on personal communications with Cor Jacobs, May 2010, 
                      ! ESS_lg_20110330+ see Steeneveld et al, 2002, table 3.5, CF0(3)=0.9, F0, was 0.4 ! 0.4 selected for coniferous trees (Kassteele)
      CF0(4)=.89      ! ESS_lg_20150624+ tropical rainforest, taking similar values as C3 vegetation
	                  ! except of specific modifications given below


	  CRICO(1)=-0.07
      CRICO(2)=-0.15
      CRICO(3)=-0.15  ! ESS_lg_20110330+ see Steeneveld et al, 2002, table 3.5, Ad, was -0.15 different CRICO for vegetation type C3
      CRICO(4)=-0.07  ! ESS_lg_20150624+ tropical rainforest, taking similar values as C3 vegetation
	                  ! except of specific modifications given below
	  
      CCS=340.

      CEPSILON0(1)=.017
      CEPSILON0(2)=.014
      CEPSILON0(3)=.017 ! ESS_lg_20110330+ see Steeneveld et al, 2002, table 3.5, alpha0, not changed,
      CEPSILON0(4)=.017 ! ESS_lg_20150624+ tropical rainforest, taking similar values as C3 vegetation
	                    ! except of specific modifications given below
	  
      CAMAX25(1)=2.2  ! normally 2.2, reset to smaller value to check the Rstom
                      ! Reinder Ronda has proposed value of 1 for tropical forest !
      CAMAX25(2)=1.7
      CAMAX25(3)=0.45 ! ESS_lg_20110330+ see Steeneveld et al, 2002, table 3.5, Am,max=1.1 ! was 0.45, different CAMAX25 for vegetation type 3
      CAMAX25(4)=2.2  ! ESS_lg_20150624+ tropical rainforest, taking similar values as C3 vegetation
	                  ! except of specific modifications given below
	  
      CAMAXQ10(1)=2.
      CAMAXQ10(2)=2.
      CAMAXQ10(3)=2.
	  CAMAXQ10(4)=2.  ! ESS_lg_20150624+ tropical rainforest, taking similar values as C3 vegetation
	                  ! except of specific modifications given below

      CAMAXT1(1)=8. 
      CAMAXT1(2)=13.
      CAMAXT1(3)=8.
	  CAMAXT1(4)=8.   ! ESS_lg_20150624+ tropical rainforest, taking similar values as C3 vegetation
	                  ! except of specific modifications given below

      CAMAXT2(1)=38.
      CAMAXT2(2)=38.
      CAMAXT2(3)=38.
	  CAMAXT2(4)=38.  ! ESS_lg_20150624+ tropical rainforest, taking similar values as C3 vegetation
	                  ! except of specific modifications given below

      C9=9.
      C1DOT6=1.6
      CMA=28.9
      CMC02=44.
      C1000=1000.
      CMV=18.
      CEXPAR=.7

      ! ESS_lg_20130516+ extra constant assignments
      RD=287.05
	  RV=461.51
	  TMELT=273.16
	  CVRAD=0.55  ! ECHAM4, ECMWF =0.5
	  ! ESS_lg_20130516-
	  
      DO JL=1,klon
	    VTMPC1=RV/RD-1.
	    ZRO = PAPHM1(JL)/  &
               (RD*PTSM1M(JL)*(1.+VTMPC1*PQSAM(JL)))
	    K=Agstype
		
        ! LG- the parameter relevant to the CO2 exchange flux is CGAMMA25(K), which
        !     represents the CO2 compensation point for the C3/C4 species. Its value 
        !     is defined in the initialisation file INIAGS.f in [mg m-3] and must be 
        !     corrected for the air density to arrive at mixing ratios.

        ZGAMMA=CGAMMA25(K)*CGAMMAQ10(K)**(((PTSM1M(JL)-TMELT)-C25)/C10)
        ZGAMMA=ZGAMMA*(CMC02/CMA)*ZRO
       
        ZGM=(CGM25(K)*CGMQ10(K)**(((PTSM1M(JL)-TMELT)-C25)/C10)) / &
            ((1+EXP(CDOT3*(CGMT1(K)-(PTSM1M(JL)-TMELT))))          &
            *(1+EXP(CDOT3*((PTSM1M(JL)-TMELT)-CGMT2(K)))))

		! LG- assigning the calculated CO2 mesophyllic resistance to the parameter
        !     RMES and the cuticular resistance for the different plant types

        RMESCO2(JL)=(1./ZGM)*1000     ! ESS_lg_20140311+ bug fix; check! s m-1
        RCUTCO2(JL)=(1./CGC(K))*1000. ! s m-1

        ! LG- end

        ZAMAX=  &
         (CAMAX25(K)*CAMAXQ10(K)**(((PTSM1M(JL)-TMELT)-C25)/C10)) / &
           ((1+EXP(CDOT3*(CAMAXT1(K)-(PTSM1M(JL)-TMELT))))          &
           *(1+EXP(CDOT3*((PTSM1M(JL)-TMELT)-CAMAXT2(K)))))

		J=1  ! soiltype

        ! LG- assigning the soil moisture parameters of EC4_VDIFF to the parameters
        !     used within this routine, which are originally defined for the 
        !     ECMWF surface scheme 

        CWPWPT(J)=0.35*PWSMAX(JL) 
        CQWEVAPT(J)=1./(0.75*PWSMAX(JL)-0.35*PWSMAX(JL))
        ! ESS_lg_20130516+ modified statement to use code for multilayer soil moisture content to arrive at 
        ! ZWROOT=WSM1M		

		ZWROOT=0.
        DO JK=1,KLEVS
		  PWSAM1M(JK)=PWSM1M(JL)           ! temporarily assignment of WSM1M to PWSAM1M
		  CVROOTSA(JK)=1./FLOAT(KLEVS)     ! temporarily set to one to arrive at ZWROOT=WSM1M
          ZWROOT=ZWROOT+PWSAM1M(JK)*CVROOTSA(JK)
        ENDDO
        ! ESS_lg_20130516+

        ZSOIL=MAX(ZEPSW,AMIN1(1.,(ZWROOT-CWPWPT(J))*CQWEVAPT(J)))
        ZSOIL=2.*ZSOIL-ZSOIL**2.
        ZGM=ZGM/C1000 

        ZB=CGC(K)-(1./9.)*ZGM*C1000                                ! ESS_lg_20110623+ gmin,w/1.6-0.11*gm term from Ronda et al., 2001
        ZFMIN=(-ZB+SQRT(ZB**2.+4.*CGC(K)*ZGM*C1000))/(2*ZGM*C1000) ! ESS_lg_20110623+ equation A9 from Ronda et al., 2001
      
        ZDS=(PQS(JL)-PQSAM(JL))*C1000
        ZDS=PQS(JL)/(0.622+0.378*PQS(JL))  &
         -PQSAM(JL)/(0.622+0.378*PQSAM(JL))
		ZDS=ZDS*PAPHM1(JL)/C1000
        ZDSMAX=(ZFMIN-CF0(K))/CRICO(K)
        ZDS=MIN(ZDSMAX,MAX(ZDS,0.))
        ZF=CF0(K)*(1.-ZDS/ZDSMAX)+ZFMIN*(ZDS/ZDSMAX) ! ESS_lg_20110623+ including the VPD effect in fmin, from Ronda et al., 2001

        ! LG- originally, the formula ZCS=CCS*(CMC02/CMA)*ZRO is being used to 
        !     recalculate the predefined CO2 concentration in ppmv,to a concentration
        !     in mg m-3, by the introduction of the explicitly resolved CO2 
        !     concentrations it is possible to consider the feedback between the 
        !     assimilation rate and the CO2 concentration in the atmosphere and 
        !     the canopy.

        ZCS=PXTMCO2(JL)*(CMC02/CMA)*ZRO ! ESS_lg_20101221+ PXTMCO2 already in ppmv

        ! LG- end

        ZCI=ZF*ZCS+(1-ZF)*ZGAMMA
        ZEPSILON=CEPSILON0(K)*(ZCS-ZGAMMA)/(ZCS+2.*ZGAMMA) ! ESS_lg_20110623+ equation A3 from Ronda et al., 2001
        ZAM=ZAMAX*(1-EXP((-ZGM*(ZCI-ZGAMMA))/ZAMAX))       ! ESS_lg_20110623+ equation A3, Am = Am,max() from Ronda et al., 2001
        ZRD=ZAM/C9                                         ! ESS_lg_20110623+ dark respiration, equation A6 from Ronda et al., 2001
        ZPAR=MAX(ZEPSR,PSRFL(JL)*CVRAD) 
        ZB=(ZEPSILON*CEXPAR*ZPAR)/(ZAM+ZRD)                ! ESS_lg_20110623+ the term used in e-(x) in equation A2 from Ronda et al., 2001

        !        ZAN=ZAM+ZRD &                             ! ESS_lg_20110623+ the term Am+Rd in equation A2 from Ronda et al., 2001
        !          -((ZAM+ZRD)/(CEXPAR*PVLTM(JL))) &
        !          *(E1(ZB*EXP(-CEXPAR*PVLTM(JL)))-E1(ZB))

        ! LG- 07-2002, modified calculation of ZAN according to a modification by 
        !     Reiner Ronda, personal communication, July 2002. The correction is meant
        !     to deal with a large overestimation of the CO2 and water vapor fluxes
        !     from tropical rainforest. The correction includes the dependence of the
        !     plant-physiological processes on the nitrogen content, which should 
        !     increase with height within the canopy

        ZAN=((ZAM+ZRD)/(CEXPAR*PVLTM(JL)))*(1.-EXP(-ZB))* &
          (1.-EXP(-CEXPAR*PVLTM(JL)))
	    ZA1=1./(1.-CF0(K))
        ZDSTAR=(1./(ZA1-1.))*ZDSMAX
        ZGSC=PVLTM(JL)*(CGC(K)/C1000  &
          +ZA1*ZSOIL*ZAN/((ZCS-ZGAMMA)*(1.+ZDS/ZDSTAR)))
        ZGS=C1DOT6*ZGSC

        PRS0(JL)=1./(ZGS)
        IF (PQM1(JL).GT.PQS(JL)) PRS0(JL)=1.E5 	! original code, PWET(JL)=0. (ECMWF)

        ! LG- setting the nocturnal stomatal resistance to a arbitrarly selected
        !     large value

        IF (PSRFL(JL).LT.0.1) PRS0(JL)=1.E5

        ! LG- determining the vertical profiles of the stomatal resistance using
        !     vertical profile of radiation within the canopy

        IF (l_xtsurf_veg_mlay) THEN

          DO JK=1,nveglay

            ! LG-     introduction of the explicitly resolved CO2 concentrations 
 
            ! LG-     the formula ZCS=CCS*(CMC02/CMA)*ZRO is being used to 
            !         recalculate the predefined CO2 concentration in ppmv,to a concentration
            !         in mg m-3.

            ZCS=PXTMCO2(JL)*(CMC02/CMA)*ZRO ! ESS_lg_20101221+ PXTMCO2 in ppmv
            ZCI=ZF*ZCS+(1-ZF)*ZGAMMA
            ZEPSILON=CEPSILON0(K)*(ZCS-ZGAMMA)/(ZCS+2.*ZGAMMA)
            ZAM=ZAMAX*(1-EXP((-ZGM*(ZCI-ZGAMMA))/ZAMAX))
            ZRD=ZAM/C9
            ZPAR=MAX(ZEPSR,(RBVD(JL)*FSL(JL,JK)+  &
			                RVD(JL,JK)*(1-FSL(JL,JK)))*CVRAD) ! ESS_lg_20140313+ added (1-FSL(JL,JK))
            ZB=(ZEPSILON*CEXPAR*ZPAR)/(ZAM+ZRD)

            !            ZAN=ZAM+ZRD &
            !              -((ZAM+ZRD)/(CEXPAR*1.)) &	        ! multi-layer approach, an LAI of 1 is being
            !              *(E1(ZB*EXP(-CEXPAR*1.))-E1(ZB))	    ! used to calculate the leaf stomatal resistance

            ! LG- 07-2002, modified calculation of ZAN according to a modification by 
            !     Reinder Ronda, personal communication, July 2002. The correction is meant
            !     to deal with a large overestimation of the CO2 and water vapor fluxes
            !     from tropical rainforest. The correction includes the dependence of the
            !     plant-physiological processes on the nitrogen content, which should 
            !     increase with height within the canopy

            ZAN=((ZAM+ZRD)/(CEXPAR*1.))*(1.-EXP(-ZB))* &
              (1.-EXP(-CEXPAR*1.))

			ZA1=1./(1.-CF0(K))
            ZDSTAR=(1./(ZA1-1.))*ZDSMAX
            ZGSC=1.*(CGC(K)/C1000  &	! LAI of 1 (PVLTM)
              +ZA1*ZSOIL*ZAN/((ZCS-ZGAMMA)*(1.+ZDS/ZDSTAR)))
            ZGS=C1DOT6*ZGSC

            ! LG-     in the original code the PWET is set to zero for dew formation 
            !         conditions. This doesn't result in a large evapotranspiration since the
            !         gradient is downward. However, since the stomatal resistance is being
            !         used also for CO2 exchange (and other species), it should be set to a
            !         large value, assuming that during dew formation, uptake (photosynthesis)
            !         is limited 

            IF (PQM1(JL).GT.PQS(JL)) THEN
              PRS0_AGSML(JL,JK)=1.E5 ! original code, PWET(JL)=0. 
            ELSE
              PRS0_AGSML(JL,JK)=1./(ZGS)
            ENDIF

            ! LG-     setting the nocturnal stomatal resistance to a arbitrarly selected
            !         large value

            IF (PSRFL(JL).LT.0.1.OR.PVLTM(JL).LT.0.01) PRS0_AGSML(JL,JK)=1.E5
			
          ENDDO

		ENDIF

        ! LG- end

      ENDDO 

  END SUBROUTINE emdep_xtsurf_ags	
  
  ! ============================================================================

  SUBROUTINE solve_tridiag(a,b,c,d,x,n)
    implicit none

    !        a - sub-diagonal (means it is the diagonal below the main diagonal)
    !        b - the main diagonal
    !        c - sup-diagonal (means it is the diagonal above the main diagonal)
    !        d - right part
    !        x - the answer
    !        n - number of equations
 
    integer,intent(in) :: n
    real(8),dimension(n),intent(in) :: a,b,c,d
    real(8),dimension(n),intent(out) :: x
    real(8),dimension(n) :: cp,dp
    real(8) :: m
    integer i
 
! initialize c-prime and d-prime
    cp(1) = c(1)/b(1)
    dp(1) = d(1)/b(1)
! solve for vectors c-prime and d-prime
    do i = 2,n
       m = b(i)-cp(i-1)*a(i)
       cp(i) = c(i)/m
       dp(i) = (d(i)-dp(i-1)*a(i))/m
    enddo
! initialize x
    x(n) = dp(n)
! solve for x from the vectors c-prime and d-prime
    do i = n-1, 1, -1
      x(i) = dp(i)-cp(i)*x(i+1)
    end do
 	
  END SUBROUTINE solve_tridiag
    
END MODULE messy_emdep_xtsurf
