!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \note
!>      If parallel mode is distributed certain combination of
!>      "in_use" and "in_space" can not be used.
!>      For performance reasons it would be better to have the loops
!>      over g-vectros in the gather/scatter routines in new subprograms
!>      with the actual arrays (also the adressing) in the parameter list
!> \par History
!>      JGH (29-Dec-2000) : Changes for parallel use
!>      JGH (13-Mar-2001) : added timing calls
!>      JGH (26-Feb-2003) : OpenMP enabled
!>      JGH (17-Nov-2007) : Removed mass arrays
!>      JGH (01-Dec-2007) : Removed and renamed routines
!>      03.2008 [tlaino] : Splitting pw_types into pw_types and pw_methods
!> \author apsi
! **************************************************************************************************
MODULE pw_types

   USE kinds,                           ONLY: dp
   USE pw_grid_types,                   ONLY: pw_grid_type

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: pw_type, pw_p_type
   PUBLIC :: pw_create

! **************************************************************************************************
   TYPE pw_type
      REAL(KIND=dp), DIMENSION(:), POINTER :: cr
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER :: cr3d
      COMPLEX(KIND=dp), DIMENSION(:), POINTER :: cc
      COMPLEX(KIND=dp), DIMENSION(:, :, :), POINTER :: cc3d

      INTEGER :: in_use ! Which data is used [r1d/c1d/r3d/c3d]
      INTEGER :: in_space ! Real/Reciprocal space
      INTEGER :: id_nr ! unique identifier
      INTEGER :: ref_count ! reference count

      TYPE(pw_grid_type), POINTER :: pw_grid
   END TYPE pw_type

! **************************************************************************************************
   TYPE pw_p_type
      TYPE(pw_type), POINTER :: pw
   END TYPE pw_p_type

   ! Flags for the structure member 'in_use'
   INTEGER, PARAMETER, PUBLIC :: REALDATA1D = 301, COMPLEXDATA1D = 302
   INTEGER, PARAMETER, PUBLIC :: REALDATA3D = 303, COMPLEXDATA3D = 304, NODATA = 305

   ! Flags for the structure member 'in_space'
   INTEGER, PARAMETER, PUBLIC :: NOSPACE = 371, REALSPACE = 372, RECIPROCALSPACE = 373
   INTEGER, PUBLIC, PARAMETER :: SQUARE = 391, SQUAREROOT = 392

   ! to generate unique id_nr
   INTEGER, SAVE, PRIVATE :: last_pw_id_nr = 0
   INTEGER, SAVE, PRIVATE :: allocated_pw_count = 0

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pw_types'
   LOGICAL, PARAMETER, PRIVATE :: debug_this_module = .FALSE.

CONTAINS

! **************************************************************************************************
!> \brief allocates and initializes pw_type
!> \param pw the type that will bw allocated and initialized
!> \param pw_grid ...
!> \param use_data which kind of data will be used
!> \param in_space in which space the pw is (real or reciprocal)
!> \param cr3d_ptr pointer with the cr3d data (make sense only if
!>        use_data==REALDATA3D)
!> \par History
!>      11.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE pw_create(pw, pw_grid, use_data, in_space, cr3d_ptr)
      TYPE(pw_type), POINTER                             :: pw
      TYPE(pw_grid_type), POINTER                        :: pw_grid
      INTEGER, INTENT(in)                                :: use_data
      INTEGER, INTENT(in), OPTIONAL                      :: in_space
      REAL(KIND=dp), DIMENSION(:, :, :), OPTIONAL, &
         POINTER                                         :: cr3d_ptr

      CHARACTER(len=*), PARAMETER :: routineN = 'pw_create', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      INTEGER, DIMENSION(:, :), POINTER                  :: bounds

      CALL timeset(routineN, handle)
      ALLOCATE (pw)
      
      last_pw_id_nr = last_pw_id_nr+1
      pw%id_nr = last_pw_id_nr
      pw%ref_count = 1
      NULLIFY (pw%pw_grid)
      pw%in_use = use_data
      pw%pw_grid => pw_grid
      pw%in_space = NOSPACE
      bounds => pw%pw_grid%bounds_local
      
      allocated_pw_count = allocated_pw_count+1
      
      NULLIFY (pw%cr, pw%cc, pw%cr3d, pw%cc3d)
      
      SELECT CASE (use_data)
      CASE (REALDATA1D)
         ALLOCATE (pw%cr(pw%pw_grid%ngpts_cut_local))
         
      CASE (COMPLEXDATA1D)
         ALLOCATE (pw%cc(pw%pw_grid%ngpts_cut_local))
         
      CASE (REALDATA3D)
         IF (PRESENT(cr3d_ptr)) THEN
            IF (ASSOCIATED(cr3d_ptr)) THEN
               pw%cr3d => cr3d_ptr
            END IF
         END IF
         IF (.NOT. ASSOCIATED(pw%cr3d)) THEN
            ALLOCATE (pw%cr3d( &
                      bounds(1, 1):bounds(2, 1), &
                      bounds(1, 2):bounds(2, 2), &
                      bounds(1, 3):bounds(2, 3)))
         END IF

      CASE (COMPLEXDATA3D)
         ALLOCATE (pw%cc3d( &
                   bounds(1, 1):bounds(2, 1), &
                   bounds(1, 2):bounds(2, 2), &
                   bounds(1, 3):bounds(2, 3)))
      CASE (NODATA)
      END SELECT
      IF (PRESENT(in_space)) pw%in_space = in_space
      CALL timestop(handle)
   END SUBROUTINE pw_create

END MODULE pw_types
