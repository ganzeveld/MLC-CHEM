! ************************************************************************
MODULE messy_main_tools
! ************************************************************************

  ! MESSy-SMCL tools
  !
  ! Author: Patrick Joeckel, MPICH, Mainz, Sep 2003

  USE messy_main_constants_mem, ONLY: DP, I4

  IMPLICIT NONE
  PRIVATE

  TYPE PTR_1D_ARRAY
     REAL(DP), DIMENSION(:), POINTER :: PTR
  END TYPE PTR_1D_ARRAY

  TYPE PTR_2D_ARRAY
     REAL(DP), DIMENSION(:,:), POINTER :: PTR
  END TYPE PTR_2D_ARRAY

  TYPE PTR_3D_ARRAY
     REAL(DP), DIMENSION(:,:,:), POINTER :: PTR
  END TYPE PTR_3D_ARRAY

  TYPE PTR_4D_ARRAY
     REAL(DP), DIMENSION(:,:,:,:), POINTER :: PTR
  END TYPE PTR_4D_ARRAY

  TYPE PTR_5D_ARRAY
     REAL(DP), DIMENSION(:,:,:,:,:), POINTER :: PTR
  END TYPE PTR_5D_ARRAY

  PUBLIC :: PTR_1D_ARRAY, PTR_2D_ARRAY, PTR_3D_ARRAY &
          , PTR_4D_ARRAY, PTR_5D_ARRAY

  ! mz_ht_20040414+
  ! FOR LOOKUP TABLE INITIALIZATION (CONVECTION at al.)
  INTEGER,  PARAMETER, PUBLIC :: jptlucu1 =  50000  ! lookup table lower bound
  INTEGER,  PARAMETER, PUBLIC :: jptlucu2 = 400000  ! lookup table upper bound
  ! table - e_s*Rd/Rv
  REAL(dp), TARGET, SAVE,  PUBLIC :: tlucua(jptlucu1:jptlucu2)
  ! table - for derivative calculation
  REAL(dp), TARGET, SAVE,  PUBLIC :: tlucub(jptlucu1:jptlucu2)
  ! table - l/cp
  REAL(dp), TARGET, SAVE,  PUBLIC :: tlucuc(jptlucu1:jptlucu2)
  ! table
  REAL(dp), TARGET, SAVE,  PUBLIC :: tlucuaw(jptlucu1:jptlucu2)
  ! mz_ht_20040414-


  ! SUBROUTINES
  PUBLIC :: read_nml_open       ! Utilities ...
  PUBLIC :: read_nml_check      ! ... to simplify ...
  PUBLIC :: read_nml_close      ! ... namelist input
  PUBLIC :: start_message       ! standard messages for start and end ...
  PUBLIC :: end_message         ! .. of submodel-specific MESSy-routines
  PUBLIC :: iso2ind             ! find index
  PUBLIC :: ind2val             ! find value at index level
  PUBLIC :: int2str             ! convert integer to string
  PUBLIC :: strcrack            ! cracking strings into parts
  ! mz_ht_20040414+
  PUBLIC :: init_convect_tables ! lookup table for convection et al. 
  ! mz_ht_20040414-
CONTAINS

! -----------------------------------------------------------------------
  SUBROUTINE read_nml_open(lex, substr, iou, nmlstr, modstr)

    IMPLICIT NONE

    ! I/O
    LOGICAL, INTENT(OUT)                   :: lex       ! file exists ?
    CHARACTER(LEN=*), INTENT(IN)           :: substr    ! calling routine
    INTEGER,          INTENT(IN)           :: iou       ! unit
    CHARACTER(LEN=*), INTENT(IN)           :: nmlstr    ! namelist
    CHARACTER(LEN=*), INTENT(IN)           :: modstr    ! module name


    CALL start_message(TRIM(modstr), 'INITIALISATION', substr)

    ! CHECK IF FILE EXISTS
    INQUIRE(file=TRIM(modstr)//'.nml', exist=lex)
    IF (.NOT.lex) THEN
       WRITE(*,*) '*** WARNING: FILE '''//TRIM(modstr)//'.nml'&
            &//'''  NOT FOUND !'
       WRITE(*,*) ' '//TRIM(modstr)//' SWITCHED OFF !'
       RETURN
    END IF

    ! OPEN FILE
    OPEN(iou,file=TRIM(modstr)//'.nml')
    WRITE(*,*) 'Reading namelist '''//TRIM(nmlstr)//''''//&
         &' from '''//TRIM(modstr)//'.nml',''' (unit ',iou,') ...'

  END SUBROUTINE read_nml_open
! -----------------------------------------------------------------------

! -----------------------------------------------------------------------
  SUBROUTINE read_nml_check(fstat, substr, iou, nmlstr, modstr)

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(IN)                    :: fstat     ! file status
    CHARACTER(LEN=*), INTENT(IN)           :: substr    ! calling routine
    INTEGER,          INTENT(IN)           :: iou       ! unit
    CHARACTER(LEN=*), INTENT(IN)           :: nmlstr    ! namelist
    CHARACTER(LEN=*), INTENT(IN)           :: modstr    ! module name

    IF (fstat /= 0) THEN
       WRITE(*,*) '*** ERROR: READ ERROR in NAMELIST '''//TRIM(nmlstr)//''''&
            &//'in FILE '''//TRIM(modstr)//'.nml'//''' !'
       CLOSE(iou)
    ELSE
       WRITE(*,*) ' ... OK !'
    END IF

  END SUBROUTINE read_nml_check
! -----------------------------------------------------------------------

! -----------------------------------------------------------------------
  SUBROUTINE read_nml_close(substr, iou, modstr)

    IMPLICIT NONE

    ! I/O
    CHARACTER(LEN=*), INTENT(IN)           :: substr    ! calling routine
    INTEGER,          INTENT(IN)           :: iou       ! unit
    CHARACTER(LEN=*), INTENT(IN)           :: modstr    ! module name

    CLOSE(iou)

    CALL end_message(TRIM(modstr), 'INITIALISATION', substr)

  END SUBROUTINE read_nml_close
! -----------------------------------------------------------------------

! --------------------------------------------------------------------
  SUBROUTINE START_MESSAGE(modstr, str, substr)
    
    IMPLICIT NONE
    
    ! I/O
    CHARACTER(LEN=*), INTENT(IN) :: modstr, str, substr

    WRITE (*,'(75("*"))')
    WRITE(*,*) '*** START ',TRIM(modstr),': ',TRIM(str), &
         ' (',TRIM(substr),')'
    
  END SUBROUTINE START_MESSAGE
! --------------------------------------------------------------------

! --------------------------------------------------------------------
  SUBROUTINE END_MESSAGE(modstr, str, substr)
    
    IMPLICIT NONE
    
    ! I/O
    CHARACTER(LEN=*), INTENT(IN) :: modstr, str, substr
    
    WRITE(*,*) '*** END ',TRIM(modstr),': ',TRIM(str), &
         ' (',TRIM(substr),')'
    WRITE (*,'(75("*"))')

  END SUBROUTINE END_MESSAGE
  ! --------------------------------------------------------------------
  
  ! -----------------------------------------------------------------------
  SUBROUTINE iso2ind(field, iso, k, f, lrev)
    
    ! ISOSURFACE TO INDEX
    ! OUTPUT IS THE LEVEL INDEX FOR A GIVEN ISOSURFACE
    ! NOTE:
    !   THIS ROUTINE IS WORKING PROPERLY ONLY IF THE 'FIELD'
    !   IS MONOTONIC
    !   (e.g., POTENTIAL TEMPERATURE, PRESSURE, POTENTIAL VORTICITY, ...)
    ! METHOD:
    !  lrev = .false. (default) -> SEARCH FROM LOW INDEX TO HIGH INDEX
    !  lrev = .true.            -> SEARCH FROM HIGH INDEX TO LOW INDEX
    !  ADJUST INDEX
    !  OPTIONAL: FRACTION OF BOX 'BELOW' ISO-SURFACE
    
    IMPLICIT NONE

    INTRINSIC :: NINT

    ! I/O
    REAL, DIMENSION(:), INTENT(IN)            :: field  ! input field
    REAL,               INTENT(IN)            :: iso    ! isosurface value
    INTEGER ,           INTENT(OUT)           :: k      ! isosurface level
    REAL,               INTENT(OUT), OPTIONAL :: f      ! fraction in layer
    LOGICAL,            INTENT(IN),  OPTIONAL :: lrev   ! reverse order ?

    ! LOCAL
    INTEGER :: nk, jk
    REAL    :: zf, dk
    LOGICAL :: llrev
    INTEGER :: jstart, jstop, jstep

    k = 0
    nk = SIZE(field)    

    IF (PRESENT(lrev)) THEN
       llrev = lrev
    ELSE
       llrev = .false. ! default
    END IF

    IF (llrev) THEN
       jstart = nk
       jstop  = 2
       jstep = -1
    ELSE
       jstart = 2
       jstop  = nk
       jstep = 1
    END IF

    DO jk = jstart, jstop, jstep
       IF ( ( (iso >= field(jk-1)) .AND. (iso <= field(jk)) ) .OR. &
            ( (iso <= field(jk-1)) .AND. (iso >= field(jk)) ) ) THEN
          k=jk
          EXIT
       END IF
    END DO

    IF ( k == 0 ) THEN   ! NOT FOUND
       IF (llrev) THEN
          k = 2
          IF (PRESENT(f)) f = 1.0
       ELSE
          k = nk
          IF (PRESENT(f)) f = 0.0
       END IF
       RETURN
    END IF

    ! ADJUST INDEX
    ! CALCULATE FRACTION OF BOX 'BELOW' ISO-SURFACE
    !
    ! METHOD: LINEAR INTERPOLATION
    !
    ! THE FOLLOWING CONDITION MUST ALWAYS BE .TRUE.,
    ! SINCE THE FIRST LEVEL WITH 
    !    FIELD(k-1) <= ISO <= FLIELD(k)
    ! OR
    !    FIELD(k-1) >= ISO >= FLIELD(k)
    ! IS SEARCHED
    !
    IF ( ABS( field(k) - field(k-1) ) > TINY(0.0) ) THEN
       zf = ABS( (iso-field(k-1)) / (field(k)-field(k-1)) )    ! e [0,1)
    ELSE
       zf = 0.5  ! SHOULD BE NEVER REACHED !!!
    END IF
    
    zf = MIN(1.,zf)
    zf = MAX(0.,zf)
    
    ! dk = INT(zf+0.5)
    dk = NINT(zf)
    ! zf e [0,0.5] -> dk = 0 -> ONE LEVEL ABOVE (-1)
    ! zf e (0.5,1) -> dk = 1 -> KEEP LEVEL
    
    k = k - 1 + dk
    
    ! CALCULATE FRACTION OF BOX 'BELOW' ISO-SURFACE
    ! ONE LEVEL ABOVE (dk = 0) -> zf e [0, 0.5]
    !                          -> f  e [0.5, 0]
    !       EXAMPLE: zf  = 0   -> ISO AT BOX MID         -> FRACT. = 0.5
    !                zf  = 0.5 -> ISO AT LOWER INTERFACE -> FRACT. = 0.0
    ! KEEP LEVEL      (dk = 1) -> zf e (0.5, 1)
    !                          -> f  e (1, 0.5)
    !       EXAMPLE: zf  = 0.5 -> ISO AT UPPER INTERFACE -> FRACT. = 1.0
    !                zf  = 1   -> ISO AT BOX MID         -> FRACT. = 0.5
    
    !                            FOR dk=1            FOR dk=0
    IF (PRESENT(f)) f = (1.5-zf)*REAL(dk) + (0.5-zf)*REAL(1-dk)
     
  END SUBROUTINE iso2ind
! -----------------------------------------------------------------------

! -----------------------------------------------------------------------
  SUBROUTINE ind2val(val, field, k, f)

    ! CONVERT INDEX (AND FRACTION 'BELOW') IN A MONOTONOUS FIELD
    ! TO THE VALUE
    ! METHOD:
    !   - PICK OUT INDEX
    !   - LINEAR INTERPOLATION BETWEEN NEIGHBOURS, IF f IS PRESENT
    
    IMPLICIT NONE
    
    ! I/O
    REAL,                   INTENT(OUT)          :: val   ! value
    REAL,     DIMENSION(:), INTENT(IN)           :: field ! field
    INTEGER,                INTENT(IN)           :: k     ! level
    REAL,                   INTENT(IN), OPTIONAL :: f     ! fraction
    
    ! LOCAL
    INTEGER :: nk
    REAL    :: ri, gf
    
    nk = SIZE(field)
    
    IF (PRESENT(f)) THEN
       ri  = 0.5 - f           ! e (-0.5,0.5) -> (top, bottom) of box
       IF (ri >= 0.0) THEN
          IF (k == nk) THEN
             gf  = (field(k)-field(k-1))
          ELSE
             gf  = (field(k+1)-field(k))
          END IF
       ELSE
          IF (k == 1) THEN
             gf  = (field(k+1)-field(k))
          ELSE
             gf  = (field(k)-field(k-1))
          END IF
       END IF
       !
       val = field(k) + ri * gf
    ELSE
       val = field(k)
    END IF
    
  END SUBROUTINE ind2val
! -----------------------------------------------------------------------

! ---------------------------------------------------------------------
  SUBROUTINE int2str(str, ii, cpad, cerr)
    
    IMPLICIT NONE
    
    INTRINSIC :: MOD
    
    ! I/O
    CHARACTER(LEN=*), INTENT(OUT)          :: str   ! STRING
    INTEGER,          INTENT(IN)           :: ii    ! INTEGER
    CHARACTER,        INTENT(IN), OPTIONAL :: cpad  ! CHAR FOR PADDING
    CHARACTER,        INTENT(IN), OPTIONAL :: cerr  ! CHAR FOR ERROR
    
    ! LOCAL
    INTEGER              :: n, zi, zn, k
    INTEGER              :: rest
    CHARACTER            :: zcpad
    
    IF (PRESENT(cpad)) THEN
       zcpad = cpad
    ELSE
       zcpad = '0'      ! DEFAULT PADDING
    END IF
    
    n  = LEN(str)
    zi = ii
    zn = n
    
    DO
       rest = MOD(zi, 10)
       zi   = zi/10
       WRITE(str(zn:zn),'(i1)') rest
       zn = zn - 1
       IF (zi == 0) EXIT
       IF (zn == 0) EXIT
    END DO
    
    IF (PRESENT(cerr)) THEN
       IF ( (zn == 0) .AND. (zi /= 0) ) THEN
          DO k = 1, n
             WRITE(str(k:k),'(a1)') cerr
          END DO
       END IF
    END IF
    
    DO k = zn, 1, -1
       WRITE(str(k:k),'(a1)') zcpad
    END DO
    
  END SUBROUTINE int2str
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
  SUBROUTINE strcrack(str, ch, el, n)

    IMPLICIT NONE
    
    ! I/O
    CHARACTER(LEN=*),               INTENT(IN)  :: str
    CHARACTER,                      INTENT(IN)  :: ch
    CHARACTER(LEN=*), DIMENSION(:), POINTER     :: el       
    INTEGER(I4),                    INTENT(OUT) :: n
    
    ! LOCAL
    INTEGER :: idx1, idx2
    
    ! INIT
    IF (ASSOCIATED(el)) DEALLOCATE(el)
    NULLIFY(el)
    n = 0
    
    ! EMPTY STRING
    IF ( (TRIM(str) == '') .OR. (TRIM(str) == ch) ) RETURN
    
    idx1 = 0
    idx2 = 0
    DO 
       idx1 = idx2 + 1
       IF (idx1 > LEN_TRIM(str(:))) EXIT
       IF (INDEX(TRIM(str(idx1:)), ch) == 0) THEN
          idx2 = LEN_TRIM(str(:)) + 1
       ELSE
          idx2 = idx2 + INDEX(TRIM(str(idx1:)), ch)
       END IF
       IF (idx1 == idx2) CYCLE
       
       n = n + 1
       
    END DO
    
    ! ALLOCATE SPACE
    ALLOCATE(el(n))
    
    n = 0
    idx1 = 0
    idx2 = 0
    DO     
       idx1 = idx2 + 1
       IF (idx1 > LEN_TRIM(str(:))) EXIT
       IF (INDEX(TRIM(str(idx1:)), ch) == 0) THEN
          idx2 = LEN_TRIM(str(:)) + 1
       ELSE
          idx2 = idx2 + INDEX(TRIM(str(idx1:)), ch)
       END IF
       IF (idx1 == idx2) CYCLE
       
       n = n + 1
       
       el(n) = str(idx1:idx2-1)
       
    END DO

  END SUBROUTINE strcrack
! ---------------------------------------------------------------------  

! ---------------------------------------------------------------------  
  ! mz_ht_20042510+
  SUBROUTINE init_convect_tables

    ! Lookup tables for convective adjustment code
    !
    ! D. Salmond, CRAY (UK), August 1991, original code

    USE messy_main_constants_mem, ONLY: R_gas, M_H2O, M_air

    IMPLICIT NONE
                       
    REAL(dp), PARAMETER :: zavl1 = -6096.9385_dp
    REAL(dp), PARAMETER :: zavl2 =    21.2409642_dp
    REAL(dp), PARAMETER :: zavl3 =    -2.711193_dp
    REAL(dp), PARAMETER :: zavl4 =     1.673952_dp
    REAL(dp), PARAMETER :: zavl5 =     2.433502_dp 

    REAL(dp), PARAMETER :: zavi1 = -6024.5282_dp
    REAL(dp), PARAMETER :: zavi2 =    29.32707_dp
    REAL(dp), PARAMETER :: zavi3 =     1.0613868_dp
    REAL(dp), PARAMETER :: zavi4 =    -1.3198825_dp
    REAL(dp), PARAMETER :: zavi5 =    -0.49382577_dp   
    
    REAL(dp), PARAMETER :: cpd   = 1005.46_dp  ! specific heat of dry air at
                                               ! constant pressure in J/K/kg
    REAL(dp), PARAMETER :: alv   = 2.5008e6_dp ! latent heat for
    !                                          ! vaporisation in J/kg
    REAL(dp), PARAMETER :: als   = 2.8345e6_dp ! latent heat for
    !                                          ! sublimation in J/kg
    REAL(dp), PARAMETER :: tmelt = 273.15_dp   ! melting temperature of
    !                                          ! ice/snow
!qqq REAL(dp), PARAMETER :: rd    = 287.05_dp   ! gas constant for
!qqq !                                          ! dry air in J/K/kg
    REAL(dp),PARAMETER  :: rd    = 1000. * R_gas / M_air ! gas constant for
    !                                          ! water vapour in J/K/kg
!qqq REAL(dp), PARAMETER :: rv    = 461.51_dp   ! gas constant for
!qqq                                            ! water vapour in J/K/kg
    REAL(dp),PARAMETER :: rv     = 1000. * R_gas / M_H2O
    ! Constants used for computation of saturation mixing ratio
    ! over liquid water (*c_les*) or ice(*c_ies*)
    REAL(dp),PARAMETER :: c3les = 17.269_dp           ! 
    REAL(dp),PARAMETER :: c3ies = 21.875_dp           ! 
    REAL(dp),PARAMETER :: c4les = 35.86_dp            ! 
    REAL(dp),PARAMETER :: c4ies = 7.66_dp             ! 
    REAL(dp),PARAMETER :: c5les = c3les*(tmelt-c4les) ! 
    REAL(dp),PARAMETER :: c5ies = c3ies*(tmelt-c4ies) ! 
   
    REAL(dp) :: z5alvcp, z5alscp, zalvdcp, zalsdcp
    REAL(dp) :: ztt, zldcp
    REAL(dp) :: zcvm3, zcvm4, zcvm5
    REAL(dp) :: zavm1, zavm2, zavm3, zavm4, zavm5

    INTEGER :: it

    z5alvcp = c5les*alv/cpd
    z5alscp = c5ies*als/cpd

    zalvdcp = alv/cpd
    zalsdcp = als/cpd

    DO it = jptlucu1, jptlucu2
      ztt = 0.001_dp*it
      IF ((ztt-tmelt) > 0.0_dp) THEN
        zcvm3 = c3les
        zcvm4 = c4les
        zcvm5 = z5alvcp
        zldcp = zalvdcp
        zavm1 = zavl1
        zavm2 = zavl2
        zavm3 = zavl3
        zavm4 = zavl4
        zavm5 = zavl5
      ELSE
        zcvm3 = c3ies
        zcvm4 = c4ies
        zcvm5 = z5alscp
        zldcp = zalsdcp
        zavm1 = zavi1
        zavm2 = zavi2
        zavm3 = zavi3
        zavm4 = zavi4
        zavm5 = zavi5
      END IF
      tlucuc(it)  = zldcp
      tlucua(it)  = EXP((zavm1/ztt+zavm2+zavm3*0.01_dp* &
           ztt+zavm4*ztt*ztt*1.e-5_dp+zavm5*LOG(ztt)))*rd/rv
      tlucub(it)  = zcvm5*(1.0_dp/(ztt-zcvm4))**2
      tlucuaw(it) = EXP((zavl1/ztt+zavl2+zavl3*0.01_dp* &
           ztt+zavl4*ztt*ztt*1.e-5_dp+zavl5*LOG(ztt)))*rd/rv
    END DO
    
  END SUBROUTINE init_convect_tables
  ! mz_ht_20042510-
! ---------------------------------------------------------------------  

! ************************************************************************
END MODULE messy_main_tools
! ************************************************************************
