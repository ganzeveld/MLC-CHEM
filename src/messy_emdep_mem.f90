! mz_lg_20050721+
MODULE messy_emdep_mem

! (messy_emdep_mem.f90 introduced for consistency with MESSy structure)

! mz_lg_2004043+ modified to make it MESSy conform: This file contains
!     those parameters being used in the messy_emdep*.f90 files which
!     are completely independent of echam5

! This module provides all variables for the EMDEP submodels

  IMPLICIT NONE
  SAVE
  !PUBLIC is default

  INTEGER, PUBLIC :: & 
      idt_O3       , & ! tracer indices for ozone, etc.
      idt_HNO3     , & 
      idt_NO       , & 
      idt_NO2      , & 
      idt_CO       , &
      idt_SO2      , & 
      idt_H2O2     , & 
      idt_ISOP     , & 
      idt_APIN     , &
      idt_BPIN     , & 
      idt_SQTERP   , &
      idt_PAN      , &  
      idt_HCOOH    , &
      idt_NH3      , &
      idt_APINP1A  , & 
      idt_RAD      , & ! tracers for the non-chemistry set-up
 	  
	  idt_CO2      , & ! ESS_lg_20130503+
      idt_COS      , & ! MAQ_lg_20190110+	 
	  
      idt_NOX      , & ! and added extra tracers for the CBM4 chemistry scheme of the SCM      idt_O3=2
      idt_CH4      , &
      idt_CH3O2H   , &
      idt_CH2O     , &
      idt_ALD2     , &
      idt_PAR      , &
      idt_OLE      , &
      idt_ETH      , &
      idt_ACET     , &
      idt_MGLY     , &
      idt_ISOPRD   , &
      idt_METHAC   , &
      idt_MVK      , &
      idt_MEK      , &
      idt_MPAN     , &
      idt_NTR      , &
      idt_DMS      , &
      idt_SO4      , &
      idt_ISONTR   , &
      idt_CH3CO2H  , &
      idt_NH2      , &
      idt_NH4      , &
      idt_HONO     , &
      idt_CH2OHO2H  ,&
      idt_RCHOHO2H  ,&
      idt_CH3OH    , &
      idt_CH3CN    , &
      idt_SQTERP2B , & 
      idt_SQTERP1B , &
      idt_MTTERP   , &
      idt_HEXANE   , &
      idt_BUTADIENE, &
      idt_TMBENZENE, &
      idt_NO3      , &
      idt_N2O5     , &
      idt_HNO4     , &
      idt_OH       , & 
      idt_HO2      , &
      idt_CH3O2    , &
      idt_C2O3     , &
      idt_XO2      , &
      idt_ROR      , &
      idt_XO2N     , &
      idt_RXPAR    , &
      idt_BXO2N    , &
      idt_MC3O3    , & 

      idt_ATERP    , & ! MAQ_lg_20170412+
      idt_LIMO     , &
      idt_MYRC         ! MAQ_lg_20170412-	 
 
  ! ESS_lg_20130118+ 
  INTEGER, DIMENSION(12) :: &
       ndaymonth = (/31,28,31,30,31,30,31,31,30,31,30,31/)

END MODULE messy_emdep_mem
! mz_lg_20050721-
