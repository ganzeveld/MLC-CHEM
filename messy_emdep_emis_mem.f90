! mz_sw_20040121+
MODULE messy_emdep_emis_mem

! AUTHOR:
! Swen Metzger    (metzger@mpch-mainz.mpg.de), MPI-CHEM, Jan 2004
! (messy_emdep_emis_mem.f90 introduced for consistency with MESSy structure)

! mz_lg_2004043+ modified to make it MESSy conform: This file contains
!     those parameters being used in the messy_emdep_emis.f90 file which
!     are completely independent of echam5

! This module provides all variables for the EMDEP EMIS submodel

  USE messy_emdep,            ONLY: modstr, modver

  ! mz_lg_20040430+ added the messy file with some general constants
  USE messy_main_constants_mem,   ONLY: dp
  ! mz_lg_20040430-
 
  IMPLICIT NONE
  SAVE
  !PUBLIC is default

  ! mz_LG_20021118+ declaration of species considered in the VOC emission algorithm
  INTEGER, PARAMETER :: nspec_vocemis=3 ! number of emitted VOC species
  ! mz_sw_20040121-

  INTEGER, PARAMETER :: iisop=1,imono=2,iovoc=3

  ! mz_LG_20020131+ soil-biogenic NO emissions
  INTEGER            :: ncl_noemis         ! number of NO emission classes
  INTEGER            :: itrop              ! tropical forest class index
  INTEGER            :: nday_month         ! the total number of days in the month
  INTEGER, PARAMETER :: ncl_yl95=12        ! 12 emission classes in YL95 inventory
  INTEGER, PARAMETER :: ncl_dk97=18        ! 18 emission classes in DK97 inventory
  INTEGER, PARAMETER :: ndrydays=14        ! number of dry days needed to get pulsing
  INTEGER, PARAMETER :: nday_month_max=31  ! maximum no. of days in month

  ! mz_LG_20020115 emission factor for the twelve ecosystems, wet conditions

  REAL, DIMENSION(ncl_yl95) :: &
       noemfact_wet = (/0.,0.,0.,0.,0.05,0.36,0.17,0.03,0.03,0.06,2.6,0./)

  ! mz_LG_20020115 emission factor for the twelve ecosystems, dry conditions

  REAL, DIMENSION(ncl_yl95) :: &
       noemfact_dry = (/0.,0.,0.,0.,0.37,2.65,1.44,0.22,0.22,0.40,8.6,0./)

  ! mz_LG_20020115 Stomatal Area Index (SAI) for twelve ecosystems, we have
  !     not distinguished the different SAI values for temperate and tropical
  !     regions, since this only makes a difference for grassland
  !     (0.018 vs 0.020, very small difference), and woodland
  !     (0.020 vs 0.040, of which we have taken the average)

  REAL, DIMENSION(ncl_yl95) :: &
       sai_veg = (/0.,0.,0.,0.010,0.019,0.,0.030,0.025,0.036,0.075,0.120,0.032/)

  ! mz_LG_20020115 emission fluxes for the 18 ecosystems distinguished by
  !     Davidson and Kingerlee, 1997, recalculated from kg N ha-1 yr-1 to
  !     ng N m-2 s-1 (which means multiplying it with 3.17)

  REAL, DIMENSION(ncl_dk97) :: &
       noemfact_dk97 = (/0.,11.42,12.68,0.13,0.,0.32,8.56,3.49,1.59, &
       0.63,3.80,18.70,9.83,5.71,0.13,0.95,0.,0./)

  ! mz_lg_20021118- soil-biogenic NO emission parameters

  ! mz_LG_20021118+ declaration and setting of some switches

  ! mz_sw_20040121+
  ! removed from messy_emdep_emis.f90 and modified/extended
  LOGICAL, PUBLIC :: lemis           = .false.   ! global switch (internal)
  LOGICAL, PUBLIC :: l_emis_antrop   = .false.   ! global switch (internal)
  LOGICAL, PUBLIC :: l_emis_bio      = .false.   ! global switch (internal)
  ! mz_sw_20040121-
  ! mz_sw_20040206+
  LOGICAL, PUBLIC :: lemis_annual    = .false.   ! global switch (internal)
  LOGICAL, PUBLIC :: lemis_monthly   = .false.   ! global switch (internal)
  LOGICAL, PUBLIC :: l_emis_SO2      = .false.   ! global switch (internal)
  ! mz_sw_20040206-

END MODULE messy_emdep_emis_mem
! mz_sw_20040121+
