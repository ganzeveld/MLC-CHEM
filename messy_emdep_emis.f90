MODULE messy_emdep_emis

  ! mz_sw_20040121+
  ! Purpose:
  ! ---------
  ! This modules contains all subroutines to calculate the emission flux of aerosols and 
  ! gases that can be used in echam5 as well as in boxmodels. ! mz_lg_20040501+ added
  ! The off-line emission fluxes are being read in the subroutine emdep_emis_global_start.
  ! The  on-line emission fluxes are calculated in additional routines that are 
  ! called from emdep_emis_vdiff, depending on the existence of the associated tracers.
  !
  ! Authors:
  ! ----------
  ! Laurens Ganzeveld, MPI-CHEM
  ! Swen Metzger    (metzger@mpch-mainz.mpg.de), MPI-CHEM, Jan 2004
  ! (Included SO2 OC/BC, Sea salt and Dust according to the MESSy structure.
  ! mz_sw_20040121-

  ! mz_sw_20040121+
  ! mz_lg_20040430+ modified messy_emdep_emis_mem
  USE messy_emdep_emis_mem
  ! mz_sw_20040121-

  ! mz_lg_20040429+ added the messy file with some general constants
  USE messy_main_constants_mem,   ONLY: pi, N_A, g, dp, i4
  ! mz_lg_20040429- 

  IMPLICIT NONE
  PRIVATE

  ! mz_pj_20030903+
  INTEGER :: status
  CHARACTER(len=*) ,PUBLIC, PARAMETER :: SUBMODSTR='emdep_emis'
  CHARACTER(len=*) ,PUBLIC, PARAMETER :: SUBmodver = '0.9' ! mz_lg_20040429+ (major udate of 0.3)
  ! mz_pj_20030903-

  ! mz_sw_20040121+
  ! Make subroutines accessable:
  ! mz_lg_20040428+ modified for removing some of the subroutines
  PUBLIC :: emdep_emis_read_nml_ctrl
  PUBLIC :: emdep_emis_bio_VOC
  PUBLIC :: emdep_emis_bio_VOC_MEGAN ! ESS_lg_20130817+
  PUBLIC :: emdep_emis_bio_NO
  PUBLIC :: emdep_emis_jNO3 ! mz_lg_20050408+ added
  ! mz_sw_20040121-

  ! mz_sw_20040121+
  ! define some internal dummy switches
  LOGICAL, PUBLIC :: l_emis_bio_NO_dummy    = .false.  ! global switch (namelist)
  LOGICAL, PUBLIC :: l_emis_bio_VOC_dummy   = .false.  ! global switch (namelist)
  LOGICAL, PUBLIC :: l_emis_bio_jNO3_dummy  = .false.  ! global switch (namelist) ! mz_lg_20050408+
  ! mz_lg_20040921-

  ! mz_sw_20040121+
  ! removed from messy_emdep_emis.f90 and modified
  LOGICAL, PUBLIC :: l_emis_bio_NO      = .false.! global switch (namelist)
  LOGICAL, PUBLIC :: l_emis_bio_VOC     = .false.! global switch (namelist)
  LOGICAL, PUBLIC :: l_emis_bio_VOC_MEGAN = .false. ! ESS_lg_20130817+ global switch (namelist)
  LOGICAL, PUBLIC :: l_emis_bio_jNO3    = .false.! global switch (namelist) ! mz_lg_20050408+

  LOGICAL, PUBLIC :: l_emis_bio_NO_pls  = .false.! global switch (namelist)
  LOGICAL, PUBLIC :: lcrfyl95           = .false.! global switch (namelist)
  INTEGER, PUBLIC :: iNOemclass         =    8   ! soil NO emission class, 8 resembles Decid. forest
  REAL, PUBLIC    :: zcult              =    0.0 ! cultivation index [0-1] 
  REAL, PUBLIC    :: zfert              =    0.0 ! application of fertilizer
  REAL, PUBLIC    :: zisopemfact        =    8.0 ! isoprene emission factor [ugC g-1 hr-1]
  REAL, PUBLIC    :: zmonoemfact        =    0.4 ! monoterpene emission factor [ugC g-1 hr-1]
  REAL, PUBLIC    :: zovocemfact        =    0.2 ! Other VOC emission factor [ugC g-1 hr-1]
  REAL, PUBLIC    :: zfAPIN             =    0.45! fraction of monoterpene emission flux emitted as alpha-pinene
  REAL, PUBLIC    :: zfBPIN             =    0.45! fraction of monoterpene emission flux emitted as beta-pinene
  REAL, PUBLIC    :: zfSQTERP           =    0.1 ! sesquiterpene emissions relative to total monoterpene emission flux
  REAL, PUBLIC    :: zradonemis         =  0.3e4 ! Radon emission flux [atoms m-2 s-1]
  REAL, PUBLIC    :: zco2emis           =  7.e-6 ! CO2 emission flux [umol CO2 m-2 s-1]
       
  ! mz_sw_20040121-

CONTAINS

  SUBROUTINE emdep_emis_read_nml_ctrl(status, iou)

    ! EMIS MODULE ROUTINE (CORE)
    !
    ! READ EMIS NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES
    !
    ! Author: Patrick Joeckel, MPICH, Feb 2002
    ! Modified: Laurens Ganzeveld, MPICH, 14-11-2002
    !           Swen Metzger    (metzger@mpch-mainz.mpg.de), MPI-CHEM, Jan 2004

    ! mz_lg_20040430+ added
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    ! I/O
    INTEGER, INTENT(OUT)        :: status
    INTEGER, INTENT(IN)         :: iou   ! logical I/O unit

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='emdep_emis_read_nml_ctrl'
    LOGICAL                     :: lex          ! file exists ?
    INTEGER                     :: fstat        ! file status

    NAMELIST /CTRL_EMIS/  l_emis_bio_NO,    &
                          l_emis_bio_VOC,   &
						  l_emis_bio_VOC_MEGAN, & ! ESS_20130817+
                          l_emis_bio_jNO3,  & ! mz_lg_20050408+ added
                          l_emis_bio_NO_pls,&
                          lcrfyl95,         &
                          iNOemclass,       &
                          zcult,            &
                          zfert,            &
                          zisopemfact,      &
                          zmonoemfact,      &
                          zovocemfact,      &
                          zfAPIN,           &
                          zfBPIN,           &
                          zfSQTERP,         &
                          zradonemis,       &
                          zco2emis            ! ESS_lg_20130503+						  

    status = 1 ! ERROR ON RETURN

    CALL read_nml_open(lex, substr, iou, 'CTRL_EMIS', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL_EMIS, IOSTAT=fstat)

    WRITE(*,*)'CTRL_EMIS:'
    WRITE(*,*)  'l_emis_bio_NO     = ',l_emis_bio_NO
    WRITE(*,*)  'l_emis_bio_VOC    = ',l_emis_bio_VOC
	WRITE(*,*)  'l_emis_bio_VOC_MEGAN = ',l_emis_bio_VOC_MEGAN ! ESS_lg_20130817+
    WRITE(*,*)  'l_emis_bio_jNO3   = ',l_emis_bio_jNO3 ! mz_lg_20050408+ added
    WRITE(*,*)  'l_emis_bio_NO_pls = ',l_emis_bio_NO_pls
    WRITE(*,*)  'lcrfyl95          = ',lcrfyl95
    WRITE(*,'(1a,i2)')    ' soil NO emission class [0-12] = ',iNOemclass
    WRITE(*,'(1a,f6.1)')  ' Cultivation index [0-1]    = ',zcult
    WRITE(*,'(1a,f6.1)')  ' C5H8 emission fact. [ugC g-1 hr-1] = ',zisopemfact
    WRITE(*,'(1a,f6.1)')  ' Monoterpene emission fact. [ugC g-1 hr-1] = ',zmonoemfact
    WRITE(*,'(1a,f6.1)')  ' Other VOC emission fact. [ugC g-1 hr-1]   = ',zovocemfact
    WRITE(*,'(1a,f6.2)')  ' fraction of mono. terpene emiss: A-pinene = ',zfAPIN
    WRITE(*,'(1a,f6.2)')  ' fraction of mono. terpene emiss: B-pinene = ',zfBPIN
    WRITE(*,'(1a,f6.2)')  ' Sesquiterpene emis. relative to mono. emis= ',zfSQTERP
    WRITE(*,'(1a,f6.1)')  ' radon emission flux [atoms m-2 s-1]       = ',zradonemis
    WRITE(*,'(1a,e8.2)')  ' CO2 emis/respiration flux [mol m-2 s-1]   = ',zco2emis
	
    CALL read_nml_check(fstat, substr, iou, 'CTRL_EMIS', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR
    lemis=.true.

  END SUBROUTINE emdep_emis_read_nml_ctrl

  !==============================================================================

  PURE SUBROUTINE emdep_emis_bio_VOC( nstep,                & ! ESS_lg_20120722+
    nveglay_hr, nveglay, iisop, imono, iovoc,               & ! mz_lg_20050725+
    dm, lad, voc_emfact, tslm1, rbvd, rvd, fsl, voc_emflux, & ! mz_lg_20040423+
    isop_emflux, mono_emflux, ovoc_emflux)                    ! mz_lg_20050725+

    ! ----------------------------------------------------------------------
    !     This program calculates the emission of volatile organic compounds
    !     from vegetation as a function of: biome (Leaf Area Index) and
    !     temperature and Photosynthetically Active Radiation (PAR). The model
    !     considers the extinction of PAR within the canopy as a function
    !     of the Leaf Area Index which is derived from the Olson ecosystems
    !     database (1992), which discerns 72 ecosystems and
    !     their characteristics (see CDROM and paper by Guenther et al., 1995).
    !     This Olson database is also applied to distinguish between different
    !     biomes which show distinct different standard emission factors
    !     (the emission rate taken at a standard temperature and for a standard
    !     amount of PAR, e.g. 30 degrees C and 1000 umol m-2 s-1). The
    !     units of the GEIA database which contains monthly average
    !     emission fluxes is mg C m-2 month-1 and in order to compare
    !     this model with these data the same units are applied where
    !     possible. This model version applies the model of Weiss and Norman
    !     (1985) to calculate the extinction of PAR as a function of the the
    !     Leaf Area Index, the distribution of the LAI (Leaf Area Density),
    !     the fraction of leaves and the orientation of these leaves. This in
    !     contrast with the original model used for the GEIA emission inventory
    !     which applies the formulas by Norman, 1982.
    !
    !     Laurens Ganzeveld 1997, modified October 2001 for implementation in
    !     ECHAM4/5 f90 versions !
    ! ----------------------------------------------------------------------
    ! mz_sw_20040206 slightly rewritten (also removed all subroutine input parameters)

    ! mz_lg_20040423+ 
    ! Interface:
    ! ----------
    ! input 
    ! nstep     : # of timesteps ! ESS_lg_20120722+
    ! nveglay_hr: number of canopy layers; high resolution
    ! nveglay:  : number of canopy layers ! mz_lg_20050725+     
    ! iisop     : isoprene index, and similar for monoterpenes and other VOC's
    ! dm        : foliar density [g m-2]
    ! lad       : leaf area density profiles [fraction]
    ! voc_emfact: VOC emission factors [ug C g-1 hr-1]
    ! tslm1     : surface temperature [K]
    ! rbvd      : direct incoming radiation [ W m-2]
    ! rvd       : diffusive radiation [W m-2]
    ! fls       : fraction of sunlit leaves [-]
    ! 
    ! output 
    ! voc_emflux: VOC emission flux [molecules m-2 s-1] 
    ! isop_emflux: Isoprene emission flux per canopy layer [molecules m-2 s-1] ! mz_lg_20050725+ 
    ! mono_emflux: Monoterpenes
    ! ovoc_emflux: other VOC's
    !

    IMPLICIT NONE 

    ! I/O
    INTEGER,  INTENT(in)  :: nstep ! ESS_lg_20120722+
    INTEGER,  INTENT(in)  :: nveglay_hr, nveglay, iisop, imono, iovoc ! mz_lg_20050725+
    REAL(dp), INTENT(in)  :: dm(:), lad(:,:), voc_emfact(:,:), &
                             tslm1(:), rbvd(:), rvd(:,:), fsl(:,:)   
    REAL(dp), INTENT(out) :: voc_emflux(:,:), &
                             isop_emflux(:,:), & ! mz_lg_20050725+
                             mono_emflux(:,:), &
                             ovoc_emflux(:,:)    ! mz_lg_20050725-

    ! mz_lg_20040423-

    ! LG- local parameters

    INTEGER :: jl, ii

    REAL :: fluxshade(SIZE(rbvd),nveglay_hr), fluxsun(SIZE(rbvd),nveglay_hr),    &
            clshade(SIZE(rbvd),nveglay_hr), foldenslay(SIZE(rbvd),nveglay_hr),   &
            pardif(SIZE(rbvd),nveglay_hr)
    REAL :: clsun(SIZE(rbvd)), pardir(SIZE(rbvd)), ct(SIZE(tslm1))

    REAL, PARAMETER :: xmc=12.        ! molar weight of C

    !  -- assigning of values of the used constants, (see Guenther et al.,
    !     1993, JGR). TSC is the leaf temperature at standard conditions,
    !     The term RECALC is the recalculation factor for getting the
    !     net short wave radiation/PAR im umol m-2 s-1 instead of W m-2.
    !     This term is taken as the average of the recalculation factor
    !     for clear sky (4.24) and diffuse conditions (4.57).
    !     See Ecological Physics by J. Hage, D96-6, IMAU and the official
    !     reference is: Grace, J., Plant-Atmosphere relationships, Chapman &
    !     Hall
    !

    REAL, PARAMETER :: alpha=0.0027
    REAL, PARAMETER :: cl1=1.066
    REAL, PARAMETER :: ct1=95000.
    REAL, PARAMETER :: ct2=230000.
    REAL, PARAMETER :: tsc=303.
    REAL, PARAMETER :: tm=314.
    REAL, PARAMETER :: r=8.314
    REAL, PARAMETER :: el=0.42
    REAL, PARAMETER :: beta=0.09
    REAL, PARAMETER :: recalc=4.405

    ! mz_lg_20030702+ added
    INTEGER :: jt, i ! mz_lg_200507025+ i

    ! LG- start calculation of emissions

    !  --  initialisation of emission flux of sunlit and shaded leaves

    fluxsun=0.
    fluxshade=0.
    voc_emflux(:,iisop)=0.
    voc_emflux(:,imono)=0.
    voc_emflux(:,iovoc)=0.

    !  --  make flux voc_emfact dependent on the radiation (PAR)
    !      and temperature (see Guenther et al., 1993, JGR)
    !      The PAR is calculated from the net short wave radiation
    !      at the surface

    ! LG-  In contrast to Guenther et al. 1995 who used the monthly
    !      mean air temperature (Leemans and Cramer), we use the surface
    !      temperature. This surface temperature represents the temperature
    !      of all the four surface cover fractions. Especially for the
    !      semi-arid regions with some vegetation, this can introduce some
    !      bias since the surface temperature will be mainly controlled
    !      by the bare soil temperature. One solution is to introduce some
    !      diagnostically derived leaf temperature from the energy balance
    !      parameters and resistances (21-09-1999)

    ! LG- calculation of temperature attunation function

    ct(:)=exp((ct1*(tslm1(:)-tsc))/ &
           (r*tsc*tslm1(:)))/ &
           (1.+exp((ct2*(tslm1(:)-tm))/ &
           (r*tsc*tslm1(:))))

    !  --  Four layers within the canopy are distinguished
    !      and for each layer the extinction of PAR is determined.
    !      The amount of total biomass is distributed over these canopy
    !      layers, expressed by the Leaf Area Density (LAD) and
    !      combined with the LAI to calculate the emission from each
    !      layer and the total emission from the biome.

    !      mz_lg_20050725+ NOTE that there has been a change in the canopy
    !      layer index with in all subroutines now a consistent layer index
    !      use with no. 1 always being the crown-layer whereas the soil-canopy
    !      layer index = nveglay or nveglay_hr (dependent on resolution of
    !      calculations

    !      The direct PAR is calculated from the direct visible
    !      radiation and the zenith angle and combined with the fraction
    !      of sunlit leaves and the total biomass yielding the emission flux
    !      of the fraction directly effected by the sun. The diffuse PAR is
    !      a function of the location within the canopy and the fraction of
    !      shaded leaves (1.-FSL)

    pardir(:)=rbvd(:)*recalc      ! direct incoming PAR

    DO ii=1,nveglay_hr               ! loop vertical layers (radiation profile)
       pardif(:,ii)=rvd(:,ii)*recalc     ! diffuse PAR
       foldenslay(:,ii)=dm(:)*lad(:,ii)    ! foliar density profile
       clsun(:)=(alpha*cl1*pardir(:))/  & ! light attenuation function sunlit leaves
            (sqrt(1.+alpha**2*pardir(:)**2))
       clshade(:,ii)=(alpha*cl1*pardif(:,ii))/ & ! light att. funct. shaded leaves
            (sqrt(1.+alpha**2*pardif(:,ii)**2))
       fluxsun(:,ii)=voc_emfact(:,iisop)*clsun(:)*ct(:)* &
            foldenslay(:,ii)*fsl(:,ii) ! flux from sunlit leaves
	   fluxshade(:,ii)=voc_emfact(:,iisop)*clshade(:,ii)*ct(:)* &
            foldenslay(:,ii)*(1.-fsl(:,ii)) ! flux from shaded leaves

	   ! LG- calculation of integrated isoprene emission rate for
       !     bulk approach, the emission is in ug C m-2 hr-1 and is
       !     recalculated to the emission in [kg C m-2 s-1] (1.E9/3600)
       !     and from that to molecules m-2 s-1. The term 1.E3 it
       !     to recalculate from kg to g, 1/XMC to recalculate to mol C
       !     and the term 1/5 is to correct for the 5 C molecules.
       !     In order to recalc from mol isoprene m-2 s-1 to molecules its
       !     multiplied with the avogadro number

       voc_emflux(:,iisop)=voc_emflux(:,iisop)+         &
          (fluxsun(:,ii)+fluxshade(:,ii))*1.e-9/(3600.)*      &
           1.e3*(1./xmc)*(1./5.)*N_A         ! total flux for use in "big leaf" model

    ENDDO ! end loop vertical layers

    ! mz_lg_20050725+ determining the emission flux for the default 2-layer
    !     canopy model needed for the subroutine xtsurf_veg_mlay
    DO i=1,nveglay  
       isop_emflux(:,i)=0._dp
       ! summing the flux!
       DO ii=(i-1)*nveglay_hr/nveglay+1,i*nveglay_hr/nveglay
         isop_emflux(:,i)=isop_emflux(:,i)+ &
               (fluxsun(:,ii)+fluxshade(:,ii))*1.e-9/(3600.)*      &
                1.e3*(1./xmc)*(1./5.)*N_A  
       END DO
    END DO
    ! mz_lg_20050725-

    ! LG- emission of monoterpenes and OVOC's

    ! mz_lg_20030811+, modified by including the recalculation to the units
    !     molecules m-2 s-1, with for monoterpenes assuming 10 C and for
    !     OVOC's 15 C per molecule

    voc_emflux(:,imono)=voc_emfact(:,imono)*exp(beta*(tslm1(:)-tsc))* &
         dm(:)*1.e-9/(3600.)*1.e3*(1./xmc)*(1./10.)*N_A
    voc_emflux(:,iovoc)=voc_emfact(:,iovoc)*exp(beta*(tslm1(:)-tsc))* &
         dm(:)*1.e-9/(3600.)*1.e3*(1./xmc)*(1./15.)*N_A

    ! mz_lg_20050725+ determining the emission flux for the default 2-layer
    !     canopy model needed for the subroutine xtsurf_veg_mlay
    DO i=1,nveglay  
       mono_emflux(:,i)=0._dp
       ovoc_emflux(:,i)=0._dp
       ! summing the flux!
       DO ii=(i-1)*nveglay_hr/nveglay+1,i*nveglay_hr/nveglay
         mono_emflux(:,i)=mono_emflux(:,i)+ &
            voc_emflux(:,imono)*lad(:,ii)  ! mz_lg_20050725+ simply scaling with LAD 
         ovoc_emflux(:,i)=ovoc_emflux(:,i)+ &
            voc_emflux(:,iovoc)*lad(:,ii) 
       END DO
    END DO
    ! mz_lg_20050725-

  END SUBROUTINE emdep_emis_bio_VOC

  !==============================================================================
  ! ESS_lg_20130817+ added MEGAN model
  SUBROUTINE emdep_emis_bio_VOC_MEGAN( nstep,                  &
    delta_time, latitude, longitude,                           & ! ESS_lg_20120722+
    nveglay_hr, nveglay, iisop, imono, iovoc,                  & ! mz_lg_20050725+
    dm, lai, lad, voc_emfact, tslm1, rbvd, rvd, fsl, fslbase, dmbase,  & ! MAQ_lg_20160817+ dmbase  MAQ_20160621+ added fslbase ! mz_lg_20040423+
    fws, voc_emflux, isop_emflux, mono_emflux, ovoc_emflux)      ! MAQ_20160621+ added fws ! mz_lg_20050725+

    ! -----------------------------------------------------------------------
    !     LG 122005: see emdep_emis_bio_VOC for more details. This code calculates the 
    !     VOC emissions but then using the updated MEGAN emission algorithm
    ! -----------------------------------------------------------------------

    ! Interface:
    ! ----------
    ! input 
    ! nstep     : # of timesteps ! ESS_lg_20120722+
    ! delta_time: timestep
	! latitude  : actual latitude (SH, -90 - 0, NH, 0-90) ! mz_lg_20040921+
	! longitude : actual longitude
    ! nveglay_hr: number of canopy layers; high resolution
    ! nveglay:  : number of canopy layers ! mz_lg_20050725+     
    ! iisop     : isoprene index, and similar for monoterpenes and other VOC's
    ! dm        : foliar density [g m-2]
	! dmbase    : long-term average foliar density [g m-2] ! MAQ_lg_20160817+
	! lai       : leaf area index [m2 m-2]
    ! lad       : leaf area density profiles [fraction]
    ! voc_emfact: VOC emission factors [ug C g-1 hr-1]
    ! tslm1     : surface temperature [K]
    ! rbvd      : direct incoming radiation [ W m-2]
    ! rvd       : diffusive radiation [W m-2]
    ! fls       : fraction of sunlit leaves [-]
    ! flsbase   : fraction of sunlit leaves for long-term average LAI ! MAQ_lg_20160621+
    ! fws       : soil moisture stress term ! MAQ_lg_20160621+
    ! 
    ! output 
    ! voc_emflux: VOC emission flux [molecules m-2 s-1] 
    ! isop_emflux: Isoprene emission flux per canopy layer [molecules m-2 s-1] ! mz_lg_20050725+ 
    ! mono_emflux: Monoterpenes
    ! ovoc_emflux: other VOC's
    !

    IMPLICIT NONE 

    ! I/O
    INTEGER,  INTENT(in)  :: nstep ! ESS_lg_20120722+
    INTEGER,  INTENT(in)  :: nveglay_hr, nveglay, iisop, imono, iovoc ! mz_lg_20050725+
    REAL(dp), INTENT(in)  :: delta_time, latitude(:), longitude(:)    ! mz_lg_20040921+
    REAL(dp), INTENT(in)  :: dm(:), lai(:), lad(:,:), &
                             tslm1(:), rbvd(:), rvd(:,:), fsl(:,:), fslbase(:,:), fws(:), & ! MAQ_lg_20160621+ fslbase & fws  
							 dmbase(:) ! MAQ_lg_20160817+ dmbase  
    REAL(dp), INTENT(inout)  :: voc_emfact(:,:)
    REAL(dp), INTENT(out) :: voc_emflux(:,:), &
                             isop_emflux(:,:), & ! mz_lg_20050725+
                             mono_emflux(:,:), &
                             ovoc_emflux(:,:)    ! mz_lg_20050725-

    ! mz_lg_20040423-

    ! LG- local parameters
    INTEGER :: klon, jl, ii
	INTEGER :: istep(SIZE(rbvd))

    REAL :: fluxshade(SIZE(rbvd),nveglay_hr,3), fluxsun(SIZE(rbvd),nveglay_hr,3),    & ! MAQ_lg_20160615+ 3 BVOC classes 
            foldenslay(SIZE(rbvd),nveglay_hr),  pardif(SIZE(rbvd),nveglay_hr)
    REAL :: pardir(SIZE(rbvd)), P240(SIZE(rbvd)), P24(SIZE(rbvd)),            & ! ESS_lg_20150409+ time dependent P240/24
            T240(SIZE(rbvd)), T24(SIZE(rbvd)),                                & ! ESS_lg_20150409+ time dependent T240/24
	        CPsun(SIZE(rbvd)), CPshade(SIZE(rbvd)), alpha(SIZE(rbvd))           ! ESS_lg_20150409+ time dependent CPsun, CPshade, alpha
			
    REAL, PARAMETER :: xmc=12.        ! molar weight of C

    LOGICAL, PARAMETER :: lemterp_light=.false. ! MAQ_lg_20160615+
	
    !  -- assigning of values of the used constants, (see Guenther et al.,
    !     1993, JGR). TSC is the leaf temperature at standard conditions,
    !     The term RECALC is the recalculation factor for getting the
    !     net short wave radiation/PAR im umol m-2 s-1 instead of W m-2.
    !     This term is taken as the average of the recalculation factor
    !     for clear sky (4.24) and diffuse conditions (4.57).
    !     See Ecological Physics by J. Hage, D96-6, IMAU and the official
    !     reference is: Grace, J., Plant-Atmosphere relationships, Chapman &
    !     Hall
    !
    REAL, PARAMETER :: cl1=1.066
    REAL, PARAMETER :: tsc=303.
    REAL, PARAMETER :: tm=314.
    REAL, PARAMETER :: r=8.314
    REAL, PARAMETER :: el=0.42
    REAL, PARAMETER :: beta=0.09
    REAL, PARAMETER :: recalc=4.405

    ! mz_lg_20030702+ added
    INTEGER :: jt, i, istart ! mz_lg_200507025+ i

    ! ESS_lg_20130817+ to read the MEGAN input data and MEGAN parameter initialization
    CHARACTER*70 chdum
    INTEGER, PARAMETER :: ngrid_megan_veg=62640, nlonm=720,nlatm=360,   &
	  iisopem=1,imonoem=2,iovocem=3,ich3ohem=4,iacetem=5,iacetaldem=6,  &
      iformaldem=7,iacetacidem=8,iformacidem=9

	! ESS_lg_20070803-

	INTEGER :: grid_megan_veg(ngrid_megan_veg), indx, indx_mlc_chem
  
    REAL ::  &
        brt(ngrid_megan_veg),fet(ngrid_megan_veg),fdt(ngrid_megan_veg), &
        shr(ngrid_megan_veg),grs(ngrid_megan_veg),crp(ngrid_megan_veg), & 
        emis_btr1(1,ngrid_megan_veg),emis_fet1(1,ngrid_megan_veg),      &
        emis_fdt1(1,ngrid_megan_veg),emis_shr1(1,ngrid_megan_veg),      & 
        emis_grs1(1,ngrid_megan_veg),emis_crp1(1,ngrid_megan_veg),      & 
        emis_ave1(1,ngrid_megan_veg)

    REAL ::  &
	      res, MEA, DEA, HEAsun, HEAshade, Clai, GAMMAa, Fnew, Anew,    &
          Fgro, Agro, Fmat, Amat, Fsen, Asen, LAIc, LAIp,               &
          T, Tg, Ti, zTm, FSL_OPTIMUM, P0sun, P0shade,                  & ! ESS_lg_20150409+ removed P240/24, CPsun, Cpshade
          x, Ct1, Ct2, Topt, T0, Eopt, fsl_avg, fsl_avgb,               & ! MAQ_lg_20160621+ ! ESS_lg_20150409+ removed T240/24
          La, CLsun, CLshade, CT                                          ! ESS_lg_20150409+ removed alpha, ESS_lg_20070803+ 

    PARAMETER (res=0.5, ANEW=0.01, AGRO=0.5, AMAT=1., ASEN=0.33, Ti=12, zTm=28, &
               Ct1=95., Ct2=230., T0=297., P0sun=200., P0shade=50.) ! ESS_lg_20120307+			   
			   
    ! INITIALIZATION
    ! 
    klon=nstep

    ! ESS_lg_20150409+ calculation of the 240 and 24 hour average PPFD and T
	! 24 hours
	DO jl=1,klon
      ! ESS_lg_20150409+ calculation of array index that determines the starting point over which the 
	  ! the time-average PPFD and T should be determined
   	  istart=MAX(1,INT(jl-INT(24*3600/delta_time)))
      P24(jl)=0.
	  T24(jl)=0.
      DO i=istart,jl
	    istep(jl)=MAX(1,INT(jl+1-istart))
        P24(jl)=P24(jl)+rbvd(i)*recalc 
        T24(jl)=T24(jl)+tslm1(i)
	  ENDDO
	  P24(jl)=MAX(1.,P24(jl)/istep(jl))
	  T24(jl)=T24(jl)/istep(jl)
	ENDDO

	! 240 hours, 10 days
	DO jl=1,klon
 	  istart=MAX(1,INT(jl-INT(240*3600/delta_time)))
      P240(jl)=0.
	  T240(jl)=0.
      DO i=istart,jl
	    istep(jl)=MAX(1,INT(jl+1-istart))
        P240(jl)=P240(jl)+rbvd(i)*recalc 
        T240(jl)=T240(jl)+tslm1(i)
	  ENDDO
	  P240(jl)=MAX(1.,P240(jl)/istep(jl))
	  T240(jl)=T240(jl)/istep(jl)

      CPsun(jl)=0.0468*EXP(0.0005*(P24(jl)-P0sun))*P240(jl)**0.6       ! equation 7, ACPD paper
      CPshade(jl)=0.0468*EXP(0.0005*(P24(jl)-P0shade))*P240(jl)**0.6   ! ESS_lg_20120307+, equation 7, ACPD paper
      alpha(jl)=0.004-0.0005*LOG(P240(jl))                             ! equation 6,  " 
 	ENDDO

    ! ESS_lg_20150409-  

    ! ESS_lg_20130817+ added reading in the MEGAN input files on plant functional types and emission factors

    ! LG- 122005, added the reading in of the file containing MEGAN input files
    !     including grid info, plant functional types and emission factors

    ! LG- 122005, plant functional types
    OPEN (10,FILE=    &                                          
        'input/MEGAN/PFT2000m302a.txt', STATUS='OLD')

    READ (10,'(70A)') chdum

    ! BTR1: broadleaf trees 
    ! FET1: Fineleaf evergreen trees 
    ! FDT1: Fineleaf deciduous trees 
    ! SHR1: shrub 
    ! GRS1: grass, non-vascular plants and other ground cover 
    ! CRP1: crop 

    DO i=1,ngrid_megan_veg
      READ (10,*) grid_megan_veg(i),brt(i),fet(i),fdt(i), &
                  shr(i),grs(i),crp(i)
    ENDDO

    CLOSE(10)

    ! LG- 122005, reading in the file with the MEGAN isoprene emission factors
	! Emission factor files: EFtccccyyyyrrrvv.txt
    ! where t is emission type which indicates which emission algorithms are used 
	! cccc is chemical species (e.g., isop= isoprene)
	! yyyy = year (Emission factors can vary for different years because plant species composition can change)
    ! rrr= grid resolution: m30 is 30 minute spatial resolution
    ! vv = file version number. Version 2a was released in November 2005 and is described by Guenther et al. (2005).
    ! The file includes a column for each emission factor that varies spatially. This differes for different compounds so the header row is different for different compounds. For isoprene, which only has type 1 emission, the header row is
	! m30,BTR1,FET1,FDT1,SHR1,GRS1,CRP1,AVE1
    ! for alpha-pinene, which has type 1 and type 2 emission, the header row is
	! m30,BTR1,FET1,FDT1,SHR1,GRS1,CRP1,AVE1, BTR2,FET2,FDT2,SHR2,GRS2,CRP2,AVE2
    ! m30 is the grid cell number (see above to get lat/lon)
	! BTR1 is the broadleaf tree emission factor (micrograms compound m-2 h-1) for type 1 emission (i.e., use type n emission algorithms to estimate variations from these emission factors) 
	! FET1 is the Fineleaf evergreen tree emission factor for type 1 emission
	! FDT1 is the Fineleaf deciduous tree emission factor for type 1 emission
    ! SHR1 is the shrub emission factor for type 1 emission
	! GRS1 is the grass, non-vascular plants and other ground cover emission factor for type 1 emission
    ! CRP1 is the crop emission factor for type 1 emission
	! AVE1 is the landscape weighted average emission factor for type 1 emission that accounts for bare ground and water and so is the average for the entire grid cell (not just for the vegetated surface)
	! BTR2 is the broadleaf tree emission factor for type 2 emission
	! and so on ...

    OPEN (10,FILE=  &                                            
        'input/MEGAN/EFisop2000m302a.txt',STATUS='OLD')

	READ (10,'(70A)') chdum
    DO i=1,ngrid_megan_veg
      READ (10,*) grid_megan_veg(i),       &
            emis_btr1(1,i),emis_fet1(1,i), & 
            emis_fdt1(1,i),emis_shr1(1,i), & 
            emis_grs1(1,i),emis_crp1(1,i), & 
            emis_ave1(1,i)
	ENDDO

    CLOSE(10)

    ! mz_lg_20051230+ calculation of MEGAN index that resembles the
    !    location of the model indicated by latitude and longitude in degrees

    indx=-999
    indx_mlc_chem=nlonm*((90.-latitude(1))/res)
    indx_mlc_chem=INT(nlonm*INT(indx_mlc_chem/nlonm)+ABS(-180-longitude(1))/res)

    DO i=1,ngrid_megan_veg
      IF (grid_megan_veg(i).EQ.indx_mlc_chem) indx=i
    ENDDO
    IF (indx.EQ.-999) THEN 
      WRITE(*,'(1a)') &
            ' emis_bio_voc_megan: data not available for this grid'
      WRITE(*,'(1a)') &
            ' emission factors are set to zero'
	  voc_emfact(:,iisopem)=0.
    ENDIF

    ! The emission factor in MEGAN is in ug compound m-2 hr-1 in
    ! contrast to the emission factor given in the G95 algorithm 
    ! which is ug C g-1 hr-1. Consequently, we have to scale the
    ! the MEGAN emission factors with foliar density [g m-2] to
    ! arrive at comparable units. 

    ! MAQ_lg_20160817+ calculation of this term below, normally used in the calculation
	! of the MEA (see below), could be applied to correc the canopy-scale emission factor to
	! to the leaf-scale emission factor but before doing so; THIS SHOULD BE CAREFULLY CHECKED
    ! voc_emfact * Clai /dm?
	
    ! calculation of the canopy environment (radiation & T) impact on effective emissions
    ! Clai(:)=0.49*LAI(:)/((1.+0.2*LAI(:)*LAI(:))**0.5)           ! equation 15 ACP-6-3181 paper, in equation 3, 3.2.1
    ! MAQ_lg_20160817-
	
    IF (dmbase(1).GT.0..AND.indx.GT.0)          & ! MAQ_lg_20160817+ dmbase rather than dm                    
       voc_emfact(:,iisopem)=                   & 
     	   ((brt(indx)*emis_btr1(1,indx) +      &
             fet(indx)*emis_fet1(1,indx) +      &
             fdt(indx)*emis_fdt1(1,indx) +      &
             shr(indx)*emis_shr1(1,indx) +      &
             grs(indx)*emis_grs1(1,indx) +      &
             crp(indx)*emis_crp1(1,indx))/100.) & ! to correct the percentage cover
     			      /dmbase(1)   ! MAQ_lg_20160817+ use the long-term average LAI/DM to scale the emission factor ! 

    WRITE(*,'(1a,f8.1)') &
      ' The MEGAN canopy-scale C5H8 emission factor [ug C m-2 hr-1] is: ', &
        voc_emfact(1,iisopem)*DMBASE(1) ! MAQ_lg_20160817+ use the long-term average LAI!
	print *,'ENTER TO CONTINUE'
	read (*,*)

    ! LG- start calculation of emissions

    !  --  initialisation of emission flux of sunlit and shaded leaves

    fluxsun=0.
    fluxshade=0.
    voc_emflux(:,iisop)=0.
    voc_emflux(:,imono)=0.
    voc_emflux(:,iovoc)=0.

	! NOTE that Guenther et al. 1995 defined an leaf level emission factor 
	! (i.e., as the emission you would get if the entire leaf were exposed to 
	! PPFD of 1000 micromol m-2 s-1) while MEGAN defines a canopy level 
	! emission factor (i.e. this is the actual canopy scale emission that you 
	! get for the standard conditions).
	! So, for comparison with the Guenther et al. 1995 EF, this must be 
	! multiplied by a factor that accounts for the fact that not all leaves 
	! are in sunlight (this is usually ~0.4 for an LAI of 5-6), or the MEGAN
	! emission factor used in the calculation of the actual emissions found	
	! below must be divided by the total fraction of sunlit leaves

	fsl_optimum=0.15  ! the fraction of sunlit leaves at 1500 umol m-2 s-1
                      ! this depends on the amount and distribution of biomass 

    ! ==========================================================================
    ! ESS_lg_20130820+ start of longitude loop

    DO jl=1,klon
	
	  ! LG- initialisation of integrated emission rates 
      ! which are being used in the hydrocarbon chemistry mode without
      ! considering the interactions within the biosphere

      Topt=313.+(0.6*(T240(jl)-T0))                             ! equation 8, ACP manuscript on MEGAN
      Eopt=2.038*EXP(0.05*(T24(jl)-T0))*EXP(0.05*(T240(jl)-T0)) ! equation 9, ACP manuscript on MEGAN
      x=((1./Topt)-(1./tslm1(jl)))/0.00831
      CT=Eopt*(Ct2*EXP(Ct1*x)/(Ct2-Ct1*(1.-EXP(Ct2*x))))        ! equation 15 ACP-6-3181 paper, equation 5, ACPD manuscript on MEGAN
  
      !  -- A number of layers within the canopy are distinguished
      !     and for each layer the extinction of PAR is determined.
      !     The amount of total biomass is distributed over these
      !     layers, expressed by the Leaf Area Density (LAD) and 
      !     combined with the LAI/DM to calculate the emission from each 
      !     layer and the total emission from the biome.  
	  !     The direct PAR is calculated from the direct visible
      !     radiation and the zenith angle and combined with the fraction
      !     of sunlit leaves and the total biomass yielding the emission flux 
      !     of the fraction directly effected by the sun. The diffuse PAR is
      !     a function of the location within the canopy and the fraction of
      !     shaded leaves (1.-FSL)  

      ! MAQ_lg_20160817+ calculation of the canopy average fraction of sunlit leaves;
	  !     first the actual one (for the LAI) and then below for the long-term average LAI
	  !     These terms were initially used to correct the emission flux for the fact
	  !     that MEGAN gives the canopy scale emission factor whereas here we deal with
	  !     the leaf-scale emission factor then being upscaled. Rather then using fsl_avg
	  !     it seems more obvious to potentially correct for this canopy environment factors
	  !     in the recalculation of the canopy scale emission factor to the leaf-scale (see above)
      !     The code below here is now obsolete but leave it in for the time being
	  
      fsl_avg=0.          
      DO ii=1,nveglay_hr
        ! mz_lg_2006011+ determining the average fraction of sunlit leaves
        fsl_avg=fsl_avg+fsl(jl,ii)
      ENDDO
      fsl_avg=fsl_avg/nveglay_hr

      ! MAQ_lg_20160621+ determining the average fraction of sunlit leaves for the long-term LAI
	  fsl_avgb=0.
      DO ii=1,nveglay_hr
        ! mz_lg_2006011+ determining the average fraction of sunlit leaves
        fsl_avgb=fsl_avgb+fslbase(jl,ii)
      ENDDO
      fsl_avgb=fsl_avgb/nveglay_hr
      ! MAQ_lg_20160621-	  
      ! MAQ_lg_20160817+

	  pardir(jl)=rbvd(jl)*recalc            ! direct incoming PAR

      DO ii=1,nveglay_hr                    ! loop vertical layers (radiation profile)
        pardif(jl,ii)=rvd(jl,ii)*recalc     ! diffuse PAR
        foldenslay(jl,ii)=dm(jl)*lad(jl,ii) ! foliar density profile

        ! mz_lg_20060111+ notice the changes in the calculations of the light and temperature
        !       attenuation factors

        CLsun=CPsun(jl)*((alpha(jl)*pardir(jl))/  &         ! equation 4, ACP manuscript on MEGAN
           (SQRT(1.+alpha(jl)**2*pardir(jl)**2)))

		! ESS_lg_20120307+
        CLshade=CPshade(jl)*((alpha(jl)*pardif(jl,ii))/  &  ! equation 4, ACP manuscript on MEGAN
           (SQRT(1.+(alpha(jl)**2)*pardif(jl,ii)**2)))

        ! mz_lg_20051230+ consistent with MEGAN documentation, calculation of 
        !       monthly, daily and hourly emission activity factors, where we 
        !       distinguish the sunlit and shade biomass

		! MEA:
        jt=30.*86400./delta_time ! number of timesteps between current LAI (LAIc) and past LAI (LAIp)
        T=jt*delta_time/86400.   ! length of the time step expressed in days between current LAI (LAIc) and past LAI (LAIp) 
        ! ESS_lg_20130820+ modify!!!
        LAIc=LAI(jl)
		LAIp=LAI(jl)
        ! MAQ_lg_20160816+ 
		IF (jl > jt) THEN
		  LAIp=LAI(jl-jt)
		ENDIF
        ! MAQ_lg_20160816- 

        IF (T.GT.TM) THEN
          Tg=zTm
        ELSE
	      Tg=T
        ENDIF

        IF (LAIc.EQ.LAIp) THEN
          Fnew=0.
          Fgro=0.
          Fsen=0.
          Fmat=1.
        ELSEIF (LAIc.LT.LAIp) THEN
          Fnew=0.
          Fgro=0.
          Fsen=(LAIp-LAIc)/LAIp
          Fmat=1.-Fnew 
	    ELSEIF (LAIc.GT.LAIp) THEN  
          Fsen=0.
          IF (T.LE.Ti) THEN
	        Fnew=1.-(LAIp/LAIc)                           ! equation 5a
            Fgro=0.                                       ! equation 5b
	      ELSEIF (T.GT.Ti) THEN
            Fnew=(T/Ti)*(1.-(LAIp/LAIc))                  ! equation 5a
            Fgro=((Tg-Ti)/T)*(1.-(LAIp/LAIc))             ! equation 5b
          ELSEIF (T.LE.zTm) THEN
            Fmat=(LAIp/LAIc)                              ! equation 5c
          ELSEIF (T.GT.zTm) THEN
            Fmat=(LAIp/LAIc)+((T-zTm)/T)*(1.-(LAIp/LAIc)) ! equation 5c
          ENDIF
        ENDIF

        GAMMAa=Fnew*Anew+Fgro*Agro+Fmat*Amat+Fsen*Asen    ! equation 4, GAMMAa should have a value of 1 for evergreen forest

        ! ESS_lg_20120308+ based on the interpretation of the MEGAN model together with Tom Pugh
        !       check if this calculation of the MEA is correctly done. According to Tom this should only be
        !       included in the big leaf model set-up and also the correct implementation of this equation
        !       should be checked in more detail. 

        ! calculation of the canopy environment (radiation & T) impact on effective emissions
        ! Clai=0.49*LAI(jl)/((1.+0.2*LAI(jl)*LAI(jl))**0.5)           ! equation 15 ACP-6-3181 paper, in equation 3, 3.2.1
        Clai=1. ! MAQ_lg_20160817+ default the canopy model is used; so no use of canopy environment parameterization in MEA

        MEA=Clai*GAMMAa                                   ! Monthly emission activity factor; equation 2

        ! Daily emission activity; set to 1
        DEA=1.                   	                      ! Daily emission activity factor
 
        ! Hourly emission activity
        HEAsun=CLsun*CT            ! ESS_lg_20120307+, Hourly emission activity factor sunlit leaves
        HEAshade=CLshade*CT        ! ESS_lg_20120307+, Hourly emission activity factor shaded leaves
		
        !  --   make flux EMISFACT dependent on the radiation (PAR)
        !       and temperature (see Guenther et al., 1993, JGR)
        !       The PAR is calculated from the net short wave radiation
        !       at the surface (See subroutine CALCPAR) 
        fluxsun(jl,ii,iisopem)=voc_emfact(jl,iisopem)*MEA*DEA*HEAsun* &
               foldenslay(jl,ii)*fsl(jl,ii)            ! assigning the layer flux considering
                                                       ! the significantly larger fraction of
                                                       ! of sunlit leaves in top of canopy
        ! ESS_lg_20120307+ modified was set at zero
        fluxshade(jl,ii,iisopem)=voc_emfact(jl,iisopem)*MEA*DEA*HEAshade* &
               foldenslay(jl,ii)*(1.-fsl(jl,ii)) 

	    ! MAQ_lg_20160615+ light dependent monoterpene emissions
		IF (lemterp_light) THEN ! ESS_LG_20100614+ included an extra switch
          IF (jl.EQ.1 .AND. ii.EQ.1) THEN
	        WRITE(*,'(1a)')' emdep_emis_bio_VOC_MEGAN: modified repres. of terpene emiss.!'
            WRITE(*,'(1a)')' Also including the light dependence'
            READ (*,*)
          ENDIF

          fluxsun(jl,ii,imonoem)=voc_emfact(jl,imonoem)*MEA*DEA*HEAsun* &
               foldenslay(jl,ii)*fsl(jl,ii)
          fluxshade(jl,ii,imonoem)=voc_emfact(jl,imonoem)*MEA*DEA*HEAshade* &
               foldenslay(jl,ii)*(1.-fsl(jl,ii))
        ENDIF ! MAQ_lg_20160615-

	  ENDDO

	ENDDO  ! end DO jl=1,klon

	! mz_lg_20050725+ determining the emission flux for the default 2-layer
    !     canopy model needed for the subroutine xtsurf_veg_mlay
    DO i=1,nveglay  
       isop_emflux(:,i)=0._dp
       mono_emflux(:,i)=0._dp
       ! summing the flux!
       DO ii=(i-1)*nveglay_hr/nveglay+1,i*nveglay_hr/nveglay
         isop_emflux(:,i)=isop_emflux(:,i)+ &
               (fluxsun(:,ii,iisopem)+fluxshade(:,ii,iisopem))*1.e-9/(3600.)*      &
                1.e3*(1./xmc)*(1./5.)*N_A  
         mono_emflux(:,i)=mono_emflux(:,i)+ &
               (fluxsun(:,ii,imonoem)+fluxshade(:,ii,imonoem))*1.e-9/(3600.)*      &
                1.e3*(1./xmc)*(1./10.)*N_A  
	   END DO
    END DO
    ! mz_lg_20050725-
	
	! LG- emission of monoterpenes and OVOC's

    ! mz_lg_20030811+, modified by including the recalculation to the units
    !     molecules m-2 s-1, with for monoterpenes assuming 10 C and for
    !     OVOC's 15 C per molecule

	IF (.NOT.lemterp_light) &  ! MAQ_lg_20160615+
      voc_emflux(:,imono)=voc_emfact(:,imono)*exp(beta*(tslm1(:)-tsc))* &
         dm(:)*1.e-9/(3600.)*1.e3*(1./xmc)*(1./10.)*N_A
		 
    voc_emflux(:,iovoc)=voc_emfact(:,iovoc)*exp(beta*(tslm1(:)-tsc))* &
         dm(:)*1.e-9/(3600.)*1.e3*(1./xmc)*(1./15.)*N_A

    ! mz_lg_20050725+ determining the emission flux for the default 2-layer
    !     canopy model needed for the subroutine xtsurf_veg_mlay
    DO i=1,nveglay  
       IF (.NOT.lemterp_light) &  ! MAQ_lg_20160615+
	     mono_emflux(:,i)=0._dp
		 
       ovoc_emflux(:,i)=0._dp
       ! summing the flux!
       DO ii=(i-1)*nveglay_hr/nveglay+1,i*nveglay_hr/nveglay
	     IF (.NOT.lemterp_light) &  ! MAQ_lg_20160615+
           mono_emflux(:,i)=mono_emflux(:,i)+ &
              voc_emflux(:,imono)*lad(:,ii)  ! mz_lg_20050725+ simply scaling with LAD 
			  
         ovoc_emflux(:,i)=ovoc_emflux(:,i)+ &
            voc_emflux(:,iovoc)*lad(:,ii) 
       END DO
    END DO
    ! mz_lg_20050725-

  END SUBROUTINE emdep_emis_bio_VOC_MEGAN  

  ! ESS_lg_20130817-
  !==============================================================================

  SUBROUTINE emdep_emis_bio_NO(                                         &
    latitude,  nstep, init_step, ndaylen, delta_time, month,            & ! mz_lg_20040921+ 
    ncl_noemis, itrop, lstart, lcrfyl95, l_veg_mlay,                    &
    lemis_bio_NO_pls,                                                   & ! ESS_lg_20120717+ removed lemis_bio_NO_mm
    cultiv, fertil, tsoil, ws, prc, prl, prectot,                       &
    noemis_w, noemis_d, noemclass, iNOemcl, lai, slf, cpold, lspold,    & ! ESS_lg_20150204+ added iNOemcl to properly assign frc_trop ! mz_lg_20050522+ added slf
    pulsing, plsday, plsdurat, cp, lsp, pls,                            & ! mz_lg_20040426+ 
    crfyl95, no_slflux, no_emflux)                                        ! mz_lg_20050614+ further modified
    ! ---------------------------------------------------------------------
    !     This program calculates the soil-biogenic NO-emission
    !     as a function of: biome, soil wetness, soil temperature,
    !     the distribution of cultivation/agriculture, N-fertilizer loss,
    !     canopy reduction, the pulsing (which is the enhanced emission
    !     due to rainfall), and rice-emission-reduction. The program
    !     is originally developed by Peter van den Broek, 1995, and
    !     edited by Laurens Ganzeveld 1996/1998.
    !
    !     10-2001, Modified for including the code in echam4/echam5 f90. More
    !     information about the model can be found in the paper by Yienger
    !     and Levy, "Empirical model of global soil-biogenic NOx emissions"
    !     JGR 100, 1995 and Ganzeveld et al., "The influence of soil-biogenic
    !     NOx emissions on the global distribution of reactive trace gases:
    !     the role of canopy processes", submitted to JGR, 2001. There is
    !     also the option to use an alternative emission inventory by
    !     Davidson, E., and W. Kingerlee, "A global inventory of nitric oxide
    !     emissions from soils", Nutrient Cycling in Agroecosystems, 48, 37-50,
    !     1997. The two different inventories are referred to in this routine
    !     by YL95 and DK97
    !
    !     Laurens Ganzeveld, October, 2001
    ! ---------------------------------------------------------------------
    ! mz_sw_20040206 slightly rewritten (also removed all subroutine input parameters)

    ! mz_lg_20040426+ 
    ! Interface:
    ! ----------
    ! input 
    ! latitude  : actual latitude (SH, -90 - 0, NH, 0-90) ! mz_lg_20040921+
    ! nstep     : timestep
    ! init_step : initial timestep
    ! ndaylen   : length of day in seconds
    ! delta_time: timestep
    ! month     : month
    ! ncl_noemis: number of NO emission classes [12, YL95, 17, DK97]
    ! itrop     : index number that resembles the tropical forest emis. class 
    ! cultiv    : cultivation intensity [old: 0-15, new 2004: 0-1]
    ! fertil    : fertilizer application (synthetic/manure)
    ! tsoil     : soil temperature  [K]
    ! ws        : soil moisture     [m]
    ! prc       : convective rainfall [m]
    ! prl       : large-scale rainfall [m]
    ! prectot   : monthly accumulated precipitation [m]
    ! noemis_w  : wet soil emission factor [ng N m-2 s-1]
    ! noemis_d  : dry soil emission factor [ng N m-2 s-1]
    ! noemclass : NO emission class, fractional coverage [0-12, YL95, 0-17 DK97]
    ! iNOemcl,  : NO emission class index, also needed to properly assign frc_trop
    ! lai       : leaf area index [m2 m-2]
    ! slf       : land sea mask [0-1] ! mz_lg_20050522+ added
    ! lstart    : switch to indicate the start of simulation
    ! lcrfyl95  : switch to use YL95's Canopy Reduction Factor (CRF)
    ! l_veg_mlay: switch to indicate use of explicit canopy model
    ! lemis_bio_NO_pls: switch to use the pulsing subroutine
    ! 
    ! input/output 
    ! cpold     : convective rainfall record for ndrydays
    ! lspold    : large-scale rainfall record for ndrydays
    ! pulsing   : pulsing regime    [index, 1-3]
    ! plsday    : timing of pulse   [number of timesteps that pulse is active]
    ! plsdurat  : duration of pulse [days]
    !
    ! output  
    ! crfyl95   : Canopy Reduction Factor [0-1] ! mz_lg_20050614+
    ! pls       : the actual pulse  [ - ]
    ! no_slflux : NO soil emission flux [molecules m-2 s-1] ! mz_lg_20050614
    ! no_emflux : NO soil-biogenic emission flux [molecules m-2 s-1]
    !

    IMPLICIT NONE 

    ! I/O
    INTEGER,  INTENT(in)  :: nstep, init_step,                    & ! mz_lg_20040921+ removed klat, jrow
                             ndaylen, month, ncl_noemis, itrop, iNOemcl  ! ESS_lg_20150204+ iNOemcl
    LOGICAL,  INTENT(in)  :: lstart, lcrfyl95,                    &
                             l_veg_mlay, lemis_bio_NO_pls           ! ESS_lg_20120717+ removed lemis_bio_NO_mm
                             
    REAL(dp), INTENT(in)  :: delta_time
    REAL(dp), INTENT(in)  :: cultiv(:),  fertil(:), tsoil(:),     &
                             ws(:),      prc(:),    prl(:),       &
                             prectot(:),                          &
                             noemis_w(:),noemis_d(:),             &
                             noemclass(:,:), lai(:),              &  
                             slf(:) ! mz_lg_20050522+ added  

    REAL(dp), INTENT(in)  :: latitude(:)     ! mz_lg_20040921+

    REAL(dp), INTENT(inout) :: cpold(:,:), lspold(:,:),           &
                             pulsing(:), plsday(:), plsdurat(:)
    REAL(dp), INTENT(out) :: cp(:),      lsp(:),     pls(:),      &
                             crfyl95(:), no_slflux(:),            & ! mz_lg_20050614+ added
                             no_emflux(:)

    ! mz_lg_20040426-

    REAL :: fert,    fertseas,soiltemp, soilws,                   &
            prec,    prect,   fwd,      fwdagri,   cult,          &
            frc_trp, frc_oth, crf,      ks,        kc,      sai   ! mz_lg_20040921+ removed echres

    REAL, PARAMETER :: wtd=0.10         ! threshold value to distinguish wet and dry soils
    ! mz_lg_20040617+ there has been a change in the assumed fraction of NO loss through the
    !     fertilizer application. The value is reduced from 2.5 to 0.7%, which is partly 
    !     compensated for by not only considering the synthetic fertilizer application but also
    !     the animal manure which increase the total N contribution for 1995 from about 70 Tg N
    !     to about 170 Tg N!
    REAL, PARAMETER :: fertloss=0.007   ! mz_lg_20040617+ modified fraction of NO loss
                                        ! based on paper by Lex Bouwman, GBC: 2002
                                        ! "Modeling of N2O and NO emissions...": was 0.025)
    REAL, PARAMETER :: tmp1=10.         ! threshold temperature
    REAL, PARAMETER :: tmp2=30.         ! threshold temperature
    REAL, PARAMETER :: fctr=0.103

    INTEGER :: klon, i

    !--- Local Variables: ! mz_sw_20040121

    ! mz_lg_20040921+ changed declarations related to assigment of fertilizer
    !     seasonality
    INTEGER :: ic, jl, vegtype

    REAL(dp) :: mtsh, mtnh
    ! mz_lg_20040921-

    ! mz_lg_20040428+ added
    REAL(dp),PARAMETER :: amn   = 14.00_dp       ! molecular weight of N

    ! mz_LG_20020115 external arrays

    ! mz_sw_20040121+

    ! INITIALIZATION
    ! mz_lg_20040426+ modified
    klon=nstep ! SIZE(tsoil) ! ESS_lg_20120722+

    ! mz_lg_20050725+ added
    crfyl95=0._dp
    no_slflux=0._dp
    no_emflux=0._dp
    ! mz_lg_20050725-

    ! --  Definition of threshold latitude to distinguish midlatitudes
    !     from the tropics, mtsh and mtnh (MidlatitudesTropicsSH and NH)
    !     the border is in this study at 30 N and 30 S

    ! mz_lg_20040921+ modified definition of latitudes used to determine
    !     seasonality in fertilizer application, echres and trp have been 
    !     removed: not used
    mtsh=-30.0_dp   ! INT(120./echres)+1     !  = 30 s
    mtnh=30.0_dp    ! INT(60./echres)+1      !  = 30 n
    ! mz_lg_20040921- 

    ! ======================================================================
    ! mz_LG_20020115 calling of subroutine pulse in which the pulsing effect
    !     on the the NO emission is determined from the precipitation data,
    !     moreover, the total precipiation of the month is determined from
    !     the accumulated precipitation calculated in a previous simulation
    !     and this is being used to determine if there is a dry or wet
    !     season month in the tropics
    ! ==================================================================

    IF (lemis_bio_NO_pls) THEN

       ! mz_lg_20040605+ the ELSE has been removed since pls is only
       !    added as a stream element whenever lemis_bio_NO_pls

       CALL emdep_emis_bio_NOpulse(                                      &
            klon, nstep, init_step, ndaylen, delta_time,                 & ! mz_lg_20040921+ removed klat, jrow
            ndrydays, lstart, prc, prl, cpold, lspold, pulsing, plsday,  &
            plsdurat, cp, lsp, pls)    
    ENDIF

    ! ==========================================================================

    ! mz_LG_20020115 start of longitude loop

    DO jl=1,klon

       ! explicit calculation of NO emission flux each timestep

       ! mz_LG_20020115 assigning the field to local parameters

       cult=cultiv(jl)
       fert=fertil(jl)
       soiltemp=tsoil(jl)-273.15 ! in echam4, td3 was being used
       soilws=ws(jl)

       ! mz_LG_20020115 prc and prl are the convective and large scale
       !     precipation in m calculated at the end of the subroutine
       !     physc after the call to the subroutine vdiff and this
       !     routine, which implies that its not the actual rainfall
       !     but the rainfall of the previous timestep being used here
       !     in this routine.

       prec=prc(jl)+prl(jl) ! in m
       prect=prectot(jl)*1000. ! ESS_lg_20150204+, precipitation in mm

       !  --   initialisation of flux

       fwd=0.0
       fwdagri=0.0

       ic=cult

       IF (ic.EQ.1.OR.ic.EQ.-9999) cult=0.

       ! --    Bouwman data

       IF (ic.EQ.2.OR.ic.EQ.12) cult=0.20
       IF (ic.EQ.3.OR.ic.EQ.13) cult=0.50
       IF (ic.EQ.4.OR.ic.EQ.14) cult=0.75
       IF (ic.EQ.5.OR.ic.EQ.15) cult=1.

       !  --   fill flux (fwd) with biome-related emission-factor
       !       also dependent on soil wetness (A(w/d)), non-agriculture
       !

       IF (soilws.GT.wtd) THEN              !wet
          fwd=noemis_w(jl)
       ELSE                                 !dry
          fwd=noemis_d(jl)
       ENDIF

       ! mz_LG_20020115   Not performing futher corrections of emission
       !       factors ofthe DK97 inventory since initial tests show that
       !       this yields a much too large emission flux (15-10-2000)

       IF (ncl_noemis.EQ.ncl_yl95) THEN

          !  --   make flux (fwd) dependent of temperature, soil wetness
          !       This is not done for rain forest (f(w/d))

          ! mz_LG_20020115   determining the fraction of coverage with rain forest
          !       and from that the surface cover of the other ecosystems

          frc_trp=0.                                         ! ESS_lg_20150204+ 
          IF (iNOemcl.eq.itrop) frc_trp=noemclass(jl,itrop)  ! ESS_lg_20150204+ 
          frc_oth=(1.-frc_trp)

          ! mz_LG_20020115   only correction for the YL95 inventory, however still assigning
          !       the temperature correction function for interpretation

          IF (soilws.GT.wtd) then   ! wet
             IF (soiltemp.LT.0.) THEN
                fwd=frc_trp*fwd+frc_oth*0.
             ELSE IF (soiltemp.ge.0.and.soiltemp.le.tmp1) THEN
                fwd=frc_trp*fwd+frc_oth*(0.28*fwd*soiltemp)
             ELSE IF (soiltemp.GT.tmp1.and.soiltemp.le.tmp2) THEN
                fwd=frc_trp*fwd+frc_oth*(fwd*exp(fctr*soiltemp))
             ELSE IF (soiltemp.GT.tmp2) THEN
                fwd=frc_trp*fwd+frc_oth*(21.97*fwd)
             ENDIF
          ELSE                         ! dry
             IF (soiltemp.LT.0.) THEN
                fwd=frc_trp*fwd+frc_oth*0.
             ELSE IF (soiltemp.ge.0..and.soiltemp.le.tmp2) THEN
                fwd=frc_trp*fwd+frc_oth*(fwd*soiltemp/30.)
             ELSE IF (soiltemp.GT.tmp2) THEN
                fwd=frc_trp*fwd+frc_oth*fwd
             ENDIF
          ENDIF

          !  --   fill flux (fwd) for rain forest, depending on month
          !       and biome (f(w/d) rain forest), non-agriculture

          !  --   distinguishing between wet and dry season based on the total
          !       amount of precipitation of 150 mm/month for any grid square in
          !       the tropics

          ! MAQ_20160923+
          IF (jl.EQ.1.and.frc_trp.GT.0) THEN
            print *,'emdep_emis_bio_NO: Note that the soil NO emission factor for tropical forest is corrected for precipitation'
            print *,'The read-in/initialized precipitation is: ',prect,' mm (when > 150mm, than a wet soil emission factor is used)'
	        print *,'ENTER to continue'
            read (*,*)
          ENDIF
          ! MAQ_20160923-
		  
          IF (prect.LT.150.) THEN ! 5 driest months
             fwd=frc_trp*8.6+frc_oth*fwd
          ELSE
             fwd=frc_trp*2.6+frc_oth*fwd
          ENDIF

          ! mz_LG_20020115   explicit calculation of canopy reduction factor acc. to Yienger
          !       and Levy

          crf=1.
          ks=8.75
          kc=0.24
          sai=0.

          ! mz_LG_20020115   determining the grid average Stomatal Area Index (SAI) from the
          !       individual SAI values of the twelve ecosystems

          DO vegtype=1,ncl_yl95-1
             sai=sai+sai_veg(vegtype)*noemclass(jl,vegtype)
          ENDDO

          crf=(exp(-ks*sai)+exp(-kc*lai(jl)))/2.

          ! mz_lg_20050719+ assigning the soil emission flux

          !  --   fwd is in [ng N m-2 s-1] whereas the required input is in
          !       molecules NO m-2 s-1, so multiplying with 1e-9 to get g,
          !       dividing by the molecular mass of N to get moles N and
          !       then times avogadro to get molecules
          no_slflux(jl)=MAX(0.0_dp,((fwd*1.e-9)/amn)*N_A)

          ! mz_LG_20020115   Yienger and Levy canopy reduction factor, only applying this
          !       term when the bigleaf approach is being used and the switch
          !       LCRFYL95 is set to TRUE

          IF (lcrfyl95) THEN
             IF (.not.l_veg_mlay) fwd=fwd*crf
          ENDIF

          crfyl95(jl)=crf ! mz_lg_20050520+ assigning natural CRF

          !  ---------------------------------------------------------------
          !  --    start of calculations for agricultural areas
          !

          !  --   fill flux (fwdagri with biome-factor and fert,
          !       dependent of soilwetness and month (A(w/d)),agriculture

          IF (cult.ge.0.20) THEN
             fwdagri=0.36
          ENDIF

          !  --   makes flux (fwdagri dependent on temperature, soil wetness,
          !       agriculture (f(w/d))

          IF (soiltemp.LT.0.) THEN
             fwdagri=0.
          ELSE IF (soiltemp.ge.0..and.soiltemp.le.tmp1) THEN
             fwdagri=0.28*fwdagri*soiltemp
          ELSE IF (soiltemp.GT.tmp1.and.soiltemp.le.tmp2) THEN
             fwdagri=fwdagri*exp(fctr*soiltemp)
          ELSE IF (soiltemp.GT.tmp2) THEN
             fwdagri=21.97*fwdagri
          ENDIF

          !  --   add effect of fertilizer in the periods of application

          fertseas=0.

          ! mz_lg_20040921+ modified definition of dependence of 
          !     seasonality in fertilizer application which now uses
          !     in the box model simulations a predefined latitude whereas
          !     in echam5 the latitude is calculated from philat instead of
          !     using jrow, which resulted in an erroneous calculation of
          !     seasonality with the kproma system

          IF (latitude(jl).GT.mtnh) THEN
             IF (cult.GE.0.20) THEN
                IF (month.GT.4.and.month.LT.9) THEN
                   fertseas=fert
                ENDIF
             ENDIF
          ENDIF
          IF (latitude(jl).LE.mtnh.and.latitude(jl).GE.mtsh) THEN
             IF (cult.GE.0.20) THEN
                fertseas=fert
             ENDIF
          ENDIF
          IF (latitude(jl).LT.mtsh) THEN
             IF (cult.GE.0.20) THEN
                IF (month.LT.3.OR.month.GT.10) THEN
                   fertseas=fert
                ENDIF
             ENDIF
          ENDIF
          ! mz_lg_20040921-

          fwdagri=fwdagri+fertseas

          ! mz_lg_20050719+ assigning the total soil emission flux,
          ! consisting of sum of the non- and agricultural flux and 
          ! weighted with the cultivation intensity 

          !  --   see above for unit recalculations
          no_slflux(jl)=no_slflux(jl)*(1.-cult)+ &
               MAX(0.0_dp,((fwdagri*1.e-9)/amn)*N_A)*cult

          !  --   correct flux (fwdagri for canopy reduction, month
          !       agriculture

          sai=0.

          DO vegtype=12,12
             sai=sai+sai_veg(vegtype)*noemclass(jl,vegtype)
          ENDDO

          crf=(exp(-ks*sai)+exp(-kc*lai(jl)))/2.

          ! mz_LG_20020115   Yienger and Levy canopy reduction factor, only applying this
          !       term when the bigleaf approach is being used

          IF (lcrfyl95) THEN
             IF (.not.l_veg_mlay) fwdagri=fwdagri*crf
          ENDIF

          !  --   combine flux from biomes and flux from agriculture
          !       to overall flux, depending of cultivation-fraction

          fwd=fwd*(1.-cult)+fwdagri*cult

          ! mz_lg_20050520+ determining the overall CRF
          crfyl95(jl)=crfyl95(jl)*(1.-cult)+crf*cult ! mz_lg_20050520+ added crf(jl)

          !  --   add effect of pulsing to overall flux at each timestep

          IF (lemis_bio_NO_pls) THEN
             fwd=fwd*pls(jl) ! mz_lg_20040605+ modified
             no_slflux(jl)=no_slflux(jl)*pls(jl) ! mz_lg_20050719+
          ENDIF

          ! mz_LG_20020115   end IF (ncl_noemis.EQ.ncl_yl95)

       ENDIF

       !=========================================================================
       ! mz_LG_20020115 the option of the using the canopy reduction factor also
       !     for DK97 still needs to be included. This has not be done so far
       !     since the DK97 inventory has only been used in the multi-layer model,
       !     where the role of the canopy processes is explicitly calculated
       !======================================================================

       !  --   fwd is in [ng N m-2 s-1] whereas the required input is in
       !       molecules NO m-2 s-1, so multiplying with 1e-9 to get g,
       !       dividing by the molecular mass of N to get moles N and
       !       then times avogadro to get molecules

       no_emflux(jl)=MAX(0.0_dp,((fwd*1.e-9)/amn)*N_A)

       ! ESS_lg_20150811+ sensitivity analysis
       !no_slflux(jl)=0.5*no_slflux(jl)
	   
       ! mz_LG_20020115  end longitudinal loop

    ENDDO

  END SUBROUTINE emdep_emis_bio_NO

  !=============================================================================

  SUBROUTINE emdep_emis_bio_NOpulse(                             &
    klon, nstep, init_step, ndaylen, delta_time,                 & ! mz_lg_20040921+ removed klat, jrow
    ndrydays, lstart, prc, prl, cpold, lspold, pulsing, plsday,  &
    plsdurat, cp, lsp, pls)      ! mz_sw_20040121 (deleted kbdim)

    ! ------------------------------------------------------------------
    !     This program adds the convective and large scale precipitation
    !     and calculates the pulsing of the NO emission occuring after
    !     a rainfall event after a period of drought. For the pulsing,
    !     a history must be recorded and therefore the precipitation data
    !     of the previous month are also incorporated in the calculation
    !     of the NO emission pulse.
    !
    !     edited by Laurens Ganzeveld, 1998, modified for implementation
    !     in echam4/echam5 f90, Laurens Ganzeveld, October, 2001
    ! ------------------------------------------------------------------

    ! mz_lg_20040426+ 
    ! Interface:
    ! ----------
    ! input 
    ! klon      : number of longitudes
    ! mz_lg_20040921+ removed klat, jrow
    ! nstep     : timestep
    ! init_step : initial timestep
    ! ndaylen   : length of day in seconds 
    ! delta_time: timestep
    ! ndrydays  : number of dry days needed to get a pulse (default 14)
    ! prc       : convective rainfall [m]
    ! prl       : large-scale rainfall [m]
    ! 
    ! input/output 
    ! cpold     : convective rainfall record for ndrydays
    ! lspold    : large-scale rainfall record for ndrydays
    ! pulsing   : pulsing regime    [index, 1-3]
    ! plsday    : timing of pulse   [number of timesteps that pulse is active]
    ! plsdurat  : duration of pulse [days]
    ! 
    ! output  
    ! cp        : daily accumulated convective rainfall
    ! lsp       : daily accumulated large-scale rainfall
    ! cpold     : convective rainfall record previous 14 days
    ! lspold    : large-scale rainfall record previous 14 days
    ! pls       : the actual pulse  [ - ]

    IMPLICIT NONE 

    ! I/O
    INTEGER,  INTENT(in)  :: klon, nstep, init_step, & ! mz_lg_20040921+ removed klat, jrow
                             ndaylen, ndrydays
    LOGICAL,  INTENT(in)  :: lstart
    REAL(dp), INTENT(in)  :: delta_time
    REAL(dp), INTENT(in)  :: prc(:), prl(:)
    REAL(dp), INTENT(inout) :: cpold(:,:), lspold(:,:), &
                               pulsing(:), plsday(:), plsdurat(:)
    REAL(dp), INTENT(out) :: cp(:),  lsp(:), pls(:)

    ! mz_lg_20040426-

    INTEGER :: nstepday, istep

    !--- Local Variables: ! mz_sw_20040121

    INTEGER :: k, jl

    !  The parameter wetdry is the threshold for the amount
    !  of precipitation during 2 weeks below which pulsing will occur,
    !  daily average parameter values are applied for the calculation of
    !  the pulsing

    INTEGER, PARAMETER :: wetdry=10

    REAL :: sum, prect, zdayl, zdtime

    !--- 0) Initialisations: -------------------------------------------------------------

    ! mz_LG_20020115 determining the actual step of the simulation and the
    !     daylength in seconds

    zdayl=REAL(ndaylen)
    zdtime=delta_time

    ! mz_LG_20020115 determining the number of timesteps of one day

    nstepday=int(zdayl/zdtime)

    ! mz_LG_20020115   reseting the precipitation record to zero
    ! mz_LG_20040427+  Note that throughout the code modifications have
    !     have been introduced with respect of keeping track of the 
    !     of the precipitation record. The second array now resembles
    !     the maximum number of dry days needed to get a pulse and the
    !     accumulated amount of precipitation for the actual days is
    !     is now written to the second array element ndryday after the
    !     old values have all been copied in the array element - 1

    IF (lstart) THEN  ! mz_lg_20040428+ modified

       cp(:)=0._dp ! mz_lg_20040621+ moved to here
       lsp(:)=0._dp

       plsday(:)=0._dp
       pulsing(:)=0._dp
       plsdurat(:)=0._dp
       cpold(:,:)=4e-4_dp 
       lspold(:,:)=4e-4_dp  ! for initialization a small value has been 
          ! selected to avoid the calculation of pulse for first timesteps.
          ! The value 4e-4 results to an accumulated amount of precipitation
          ! of ndrydays*(4e-3+4e-3)*1000.=11.2 mm of rainfall
                         
    ENDIF

    ! mz_LG_20020115 start loop longitude

    DO jl=1,klon

       istep=nstep-init_step+1

       !  --  start calculation of pulsing as a function of total
       !      precipitation, the summed precipitation during a period of 14
       !      days

       ! mz_LG_20020115  CP AND LSP are the daily cumulative precipitation,
       !      so the actual precipitation rates are added

       cp(jl)=cp(jl)+prc(jl)
       lsp(jl)=lsp(jl)+prl(jl)

       ! mz_lg_20030118+ modified to reduce the size of the arrays, the
       !      maximum second array dimension is now resembling the duration
       !      of the drought in days. Whenever the day changes, the values
       !      of the aray are copied such that the accumulated precipitation
       !      of day no. ndrydays (default 14) are copied to array number 13
       !      ands forth. Then at the end the newly accumulated precipitation
       !      record is copied to the aray number no. ndrydays. Then the
       !      daily accumulated precipitation rates are set to zero

       IF (MOD(istep,nstepday).eq.0) THEN
          DO k=1,ndrydays-1  ! mz_lg-20040427+ modified
             cpold(jl,k)=cpold(jl,k+1)
             lspold(jl,k)=lspold(jl,k+1)
          ENDDO
          cpold(jl,ndrydays)=cp(jl)    ! mz_lg-20040427+ modified
          lspold(jl,ndrydays)=lsp(jl)  ! mz_lg-20040427+ modified
          cp(jl)=0._dp
          lsp(jl)=0._dp
       ENDIF

       !   -- determining the total precipitation of the previous
       !      ndrydays

       sum=0.
       DO k=1,ndrydays   ! mz_lg-20040427+ modified
          sum=sum+(cpold(jl,k)+lspold(jl,k))*1000.
       ENDDO

       ! mz_LG_20020115  calculating the total amount of precipitation in
       !      mm day-1 (see Yienger and Levy to determine the pulsing regimes).
       !      The daily cumulative rainfall is considered in mm

       prect=(cp(jl)+lsp(jl))*1000.

       !   -- the parameter n is the no. of pulsing in each specific class
       !      1 < prec < 5, 5 < prec < 15 and prec > 15 (n5, n10 and
       !      n15 respect.) Three pulsing regimes are discerned as a function
       !      of the intensity of the precipitation (pulsing=1,2,or 3). The
       !      pulse last for 3,7 or 14 days for the three regimes (plsdurat)
       !      starting at day=1

       ! mz_LG_20020115  all units in mm !!!

       IF (prect .GT. 1. .AND. prect .LE. 5.) THEN
          IF (sum.LT.wetdry.AND.plsday(jl).EQ.0.) THEN
             pulsing(jl)=1._dp
             plsday(jl)=1./nstepday
             plsdurat(jl)=3._dp
          ENDIF
       ELSE IF (prect .GT. 5. .AND. prect .LE. 15.) THEN
          IF (sum.LT.wetdry.AND.plsday(jl).EQ.0.) THEN
             pulsing(jl)=2._dp
             plsday(jl)=1./nstepday
             plsdurat(jl)=7._dp
          ENDIF
       ELSE IF (prect .GT. 15.) THEN
          IF (sum.LT.wetdry.AND.plsday(jl).EQ.0.) THEN
             pulsing(jl)=3._dp
             plsday(jl)=1./nstepday
             plsdurat(jl)=14._dp
          ENDIF
       ENDIF

       IF (pulsing(jl) .EQ. 1. .AND. &
            plsday(jl) .LE. plsdurat(jl)) THEN
          pls(jl)=MAX(1._dp,11.19*EXP(-0.805*plsday(jl)))
       ELSE IF (pulsing(jl) .EQ.2. .AND. &
            plsday(jl) .LE. plsdurat(jl)) THEN
          pls(jl)=MAX(1._dp,14.68*EXP(-0.384*plsday(jl)))
       ELSE IF (pulsing(jl) .EQ. 3. .AND. &
            plsday(jl) .LE. plsdurat(jl)) THEN
          pls(jl)=MAX(1._dp,18.46*EXP(-0.208*plsday(jl)))
       ELSE
          pls(jl)=1._dp
       ENDIF

       ! mz_LG_20020115  increasing of the day, so this needs to be corrected
       !      for the number of timesteps of each day

       IF (plsday(jl) .GT. 0.) plsday(jl)=plsday(jl)+1./nstepday
       IF (plsday(jl) .GT. plsdurat(jl)) THEN
          pulsing(jl)=0._dp
          plsday(jl)=0._dp
          plsdurat(jl)=0._dp
       ENDIF

       ! mz_LG_20020115 end loop longitude

    ENDDO

    ! mz_lg_20040427+ removed the whole code on writing the rainfall record.
    !     with storing the information into the streams this information is
    !     stored, even in the case of a rerun

  END SUBROUTINE emdep_emis_bio_NOpulse

  ! mz_lg_20050408+ added a new subroutine to calculate the emissions of HONO and
  !     NO2 by the photolysis of deposited HNO3 on the leaves

  !==============================================================================

  SUBROUTINE emdep_emis_jNO3( Nstep,                                    & ! ESS_lg_20120722+
    delta_time, l_xtsurf_veg_mlay, nveglay_hr, nveglay, pcvw, pvgrat,   &
    prc, prl, fsl, lad, lai, ddNO3, fslsum, cthru, NO3s, HONO_jNO3em,   &
    NOx_jNO3em, HONO_emflux, NOx_emflux) 

    ! ----------------------------------------------------------------------------
    !   added an subroutine to estimate the HONO and NO2 emission flux related
    !   to the photo-dissociation of HNO3 accumulated at the surface through
    !   dry deposition (Paper Zhou et al., GRL, 2003/2004). 
    !
    !   Laurens Ganzeveld, April, 2005
    ! ---------------------------------------------------------------------------
    ! Interface:
    ! ----------
    ! input 
    ! nstep     : # of timesteps ! ESS_lg_20120722+
    ! delta_time: timestep
    ! l_xtsurf_veg_mlay: switch to indicate use of explicit canopy model
    ! nveglay_hr: number of canopy layers; high resolution
    ! nveglay   : number of canopy layers: multi-layer model
    ! pcvw      : wet skin fraction
    ! pvgrat    : vegetation fraction
    ! prc       : convective rainfall [m]
    ! prl       : large-scale rainfall [m]
        ! fsl       : fraction of sunlit leaves per layer [0-1]
    ! lad       : leaf area density [fraction]
    ! lai       : leaf area index [m2 m-2]
    ! ddhno3    : deposition flux of HNO3 [molecules m-2 s-1]
    !
    ! input/output 
    ! NO3s     : amount of NO3 at surface due to dry deposition
    ! 
    ! output  
    ! fslsum         : total fraction of sunlit leaves [0-1]
    ! cthru          : throughfall; for a value of 1 all the accumulated nitrate is removed
    ! HONO_jNO3em    : bulk HONO biogenic emission flux [molecules m-2 s-1]
    ! NOx_jNO3em     : bulk NO2 biogenic emission flux [molecules m-2 s-1]
    ! HONO_emflux    : vegetation HONO biogenic emission flux [molecules m-2 s-1]
    ! NOx_emflux     : vegetation NO2 biogenic emission flux [molecules m-2 s-1]

    IMPLICIT NONE 

    ! I/O
    REAL(dp), INTENT(in)    :: delta_time
    LOGICAL,  INTENT(in)    :: l_xtsurf_veg_mlay
    INTEGER,  INTENT(in)    :: nstep,nveglay_hr, nveglay  ! ESS_lg_20120722+
    REAL(dp), INTENT(in)    :: pvgrat(:), pcvw(:),                  &
                               prc(:), prl(:), fsl(:,:),            &
                               lad(:,:), lai(:), ddNO3(:)  
    REAL(dp), INTENT(inout) :: NO3s(:)
    REAL(dp), INTENT(out)   :: HONO_jNO3em(:), NOx_jNO3em(:),     &
                               HONO_emflux(:,:), NOx_emflux(:,:), &
                               cthru(:)
    REAL, PARAMETER ::                                            &
        jNO3s_HONO     =  2.5e-5, & ! production rate of HONO from NO3 photolsysis, [s-1], Table 2: 50% humidity
        jNO3s_NOx      =  2.2e-5    ! Table 2: 50% humidity     

    !--- Local Variables:

    REAL(dp)  :: fslsum(SIZE(prc))
    REAL(dp)  :: ztprcp
    INTEGER   :: jl, jk, jjk, klon

    ! INITIALIZATION
    klon=nstep ! SIZE(prc) ! ESS_lg_20120722+

    ! mz_lg_20050410+ determining the total fraction of sunlit leaves
    fslsum(:)=0._dp
    DO jk=1,nveglay_hr
      fslsum(:)=fslsum(:)+fsl(:,jk)
    ENDDO

    ! mz_LG_20020115 start of loop

    DO jl=1,klon

      HONO_emflux(jl,:)=0._dp
      NOx_emflux(jl,:)=0._dp

      IF (jl>1) THEN ! ESS_lg_20130516+ modified based on bounds check
        IF (l_xtsurf_veg_mlay) THEN ! multilayer model approach
  
          ! mz_lg_20050725+ determining the emission flux for the default 2-layer
          !     canopy model needed for the subroutine xtsurf_veg_mlay
          DO jk=1,nveglay  
             ! summing the flux!
             DO jjk=(jk-1)*nveglay_hr/nveglay+1,jk*nveglay_hr/nveglay
               HONO_emflux(jl,jk)=    &
                   lad(jl,jjk)*lai(jl)*fsl(jl,jjk)*jNO3s_HONO*NO3s(jl-1)
               NOx_emflux(jl,jk)=     &
                   lad(jl,jjk)*lai(jl)*fsl(jl,jjk)*jNO3s_NOx*NO3s(jl-1)
             END DO
          END DO
          ! mz_lg_20050725-

        ELSE ! bulk model approach
  
          DO jk=1,nveglay_hr
            fslsum(:)=fslsum(:)+fsl(:,jk)
          ENDDO
          HONO_jNO3em(jl)=(1.-pcvw(jl))*pvgrat(jl)*   & ! only emission from dry vegetation
       &       lai(jl)*fslsum(jl)*jNO3s_HONO*NO3s(jl-1)
          NOx_jNO3em(jl)=(1.-pcvw(jl))*pvgrat(jl)*    &
               lai(jl)*fslsum(jl)*jNO3s_NOx*NO3s(jl-1)

          ! mz_LG_20020115  end longitudinal loop

        ENDIF ! end IF (l_veg_lay) THEN
      ENDIF ! (jl>1) THEN ! ESS_lg_20130516-
		
      ! determining the accumulation of nitrate on the leaf surface using the dry deposition flux as a source
      ! and the losses by the photolysis as a sink term
      IF (jl > 1) THEN
        IF (l_xtsurf_veg_mlay) THEN
            NO3s(jl)=NO3s(jl-1)                            &
               +ddNO3(jl)*delta_time                       &  ! ESS_lg_20130117+ to accumulate, added by dry deposition
                   -(HONO_emflux(jl,1)+HONO_emflux(jl,2)+  &  ! and lost by the photolysis process
                       NOx_emflux(jl,1)+NOx_emflux(jl,2))*delta_time
        ELSE
          NO3s(jl)=NO3s(jl)                            &
           +ddNO3(jl)*delta_time                       &  ! ESS_lg_20130117+ to accumulate, added by dry deposition
           -(HONO_jNO3em(jl)+NOx_jNO3em(jl))*delta_time   ! and lost by the photolysis process
        ENDIF
      ENDIF
    
      ! mz_lg_20050408+ determining the fraction of gridsquare wetted by
      ! rainfall. This is a time-independent constant u where a new parameterization
      ! which is also used by Meryem, gives different results. This approach should be
      ! implemented too to make it consistent with her work
      ztprcp=prc(jl)+prl(jl)
      cthru(jl)=0._dp                                  ! default, no rain, no throughfall
      IF (ztprcp.gt.0.)                             &
       &   cthru(jl)=0.2_dp*(prc(jl)/ztprcp)+       &  ! convective rainfall
       &             1.0_dp*(1.-(prc(jl)/ztprcp))      ! large scale rainfall

      ! mz_lg_20050408+ determine the accumulated amount of HNO3 at surface

      IF (ztprcp.GT.1.e-10_dp) THEN

        ! It is assumed that the all the accumulated HNO3 at that part of
        ! surface that interceps the rain is removed. The parameter CTHRU,
        ! determined in surf.f90 as a function of the relative contribution
        ! of convective rainfall to the total rainfall, is applied to
        ! estimate this HNO3 removal. The effect of rainfall interception is
        ! already partly accounted for by only calculating the emission flux
        ! for the dry vegetation fraction but it can be assumed that the
        ! rainfall is falling randomly within the grid square at an area that
        ! was dry previously. To avoid the contineous accumulation of HNO3
        ! at the dry vegetation fraction, whenever it rains it is assumed 
        ! that for a fraction resembling CTHRU the HNO3 is washed off  

        NO3s(jl)=(1.-cthru(jl))*NO3s(jl)  
      ENDIF
          
    ENDDO

  END SUBROUTINE emdep_emis_jNO3

  ! mz_lg_20050408-

  !=============================================================================

END MODULE messy_emdep_emis
