!****************************************************************************
!                Time-stamp: <2005-03-11 14:31:11 sander>
!****************************************************************************

! Definitions of machine precision constants as Fortran PARAMETERs for MESSy
! Definitions of physical constants as Fortran PARAMETERs for MESSy

! Authors:
! Rolf Sander,    MPICH, 2004: original code
! Patrick Jöckel, MPICH, 2004: preprocessor-directives removed; the 
!                              BASEMODEL now may use the constants of the
!                              Modular Earth Submodel System ...

MODULE messy_main_constants_mem

  IMPLICIT NONE
  ! PUBLIC is already default

  CHARACTER(LEN=*), PARAMETER :: modstr = 'MESSy'
  CHARACTER(LEN=*), PARAMETER :: modver = '1.0'

  INTEGER, PARAMETER :: nout = 6     ! standard output stream 
  INTEGER, PARAMETER :: nerr = 6     ! send also to standard output stream

  ! MACHINE PRECISION CONSTANTS
  INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(6,37)
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12,307)
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)
  INTEGER, PARAMETER :: wp = dp

  ! PHYSICAL CONSTANTS
  REAL(dp), PARAMETER :: pi      = 3.14159265358979323846_dp
  REAL(dp), PARAMETER :: R_gas   = 8.314409_dp  ! R [J/K/mol]
  REAL(dp), PARAMETER :: N_A     = 6.022045E23_dp ! Avogadro constant [1/mol]
  REAL(dp), PARAMETER :: g       = 9.80665_dp   ! gravity acceleration [m/s2]
  REAL(dp), PARAMETER :: T0      = 298.15_dp    ! standard temperature [K]
  REAL(dp), PARAMETER :: atm2Pa  = 101325._dp   ! conversion from [atm] to [Pa]
  REAL(dp), PARAMETER :: cal2J   = 4.1868_dp    ! conversion from [cal] to [J]
  REAL(dp), PARAMETER :: k_B     = 1.380662E-23_dp ! Boltzmann constant [J/K]

  REAL(dp), PARAMETER :: c_vKar  = 0.4          !  Karman constant [?]

  ! DRY AIR AND WATER VAPOUR THERMODYNAMIC CONSTANTS
  REAL(dp), PARAMETER :: M_H2O   = 18.0154_dp   ! molar mass of H2O [g/mol]
  REAL(dp), PARAMETER :: rho_H2O = 999.97_dp    ! density of H2O [kg/m3]
  REAL(dp), PARAMETER :: M_air   = 28.970_dp    ! molar mass of dry air [g/mol]
  REAL(dp), PARAMETER :: cp_air  = 1005.46_dp   ! specific heat of dry air at
                                                ! constant pressure in J/K/kg

  ! MXXX = molar mass of element XXX [g/mol]
  REAL(dp), PARAMETER :: MH  =   1.01
  REAL(dp), PARAMETER :: MC  =  12.01
  REAL(dp), PARAMETER :: MN  =  14.01
  REAL(dp), PARAMETER :: MF  =  19.00
  REAL(dp), PARAMETER :: MNa =  22.99
  REAL(dp), PARAMETER :: MO  =  16.00
  REAL(dp), PARAMETER :: MS  =  32.07
  REAL(dp), PARAMETER :: MCl =  35.45
  REAL(dp), PARAMETER :: MBr =  79.90
  REAL(dp), PARAMETER :: MI  = 126.90

  ! PLANETARY PARAMETERS
  REAL(dp), PARAMETER :: radius_earth = 6371000.0_dp ! radius of the Earth [m]

END MODULE messy_main_constants_mem

!*****************************************************************************
