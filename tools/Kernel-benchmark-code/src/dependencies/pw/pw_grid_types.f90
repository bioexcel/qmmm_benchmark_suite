!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \par History
!>      JGH (20-12-2000) : Parallel data layout
!> \author APSI
! **************************************************************************************************
MODULE pw_grid_types

   USE kinds,                           ONLY: dp,&
                                              int_8
   IMPLICIT NONE

   PRIVATE
   PUBLIC :: pw_grid_type

   ! (only for reciprocal grid:) fill in half or full space
   INTEGER, PARAMETER, PUBLIC :: HALFSPACE = 211, FULLSPACE = 212
   INTEGER, PARAMETER, PUBLIC :: PW_MODE_LOCAL = 0, PW_MODE_DISTRIBUTED = 1


! **************************************************************************************************
   TYPE pw_grid_type
      INTEGER(int_8) :: ngpts ! # grid points
      INTEGER(int_8) :: ngpts_cut ! # grid points within cutoff
      INTEGER, DIMENSION(2, 3) :: bounds ! lower and upper bounds
      INTEGER, DIMENSION(3) :: npts ! # point in all directions
      INTEGER :: ngpts_local ! # grid points
      INTEGER :: ngpts_cut_local ! # grid points within cutoff
      INTEGER, DIMENSION(2, 3) :: bounds_local ! bounds on local process
      INTEGER, DIMENSION(3) :: npts_local ! local version of npts
      REAL(KIND=dp), DIMENSION(3) :: dr ! grid spacing
      REAL(KIND=dp), DIMENSION(3, 3) :: dh ! incremental cell matrix
      REAL(KIND=dp), DIMENSION(3, 3) :: dh_inv ! inverse incremental cell matrix
      LOGICAL :: orthorhombic ! cell symmetry
      REAL(KIND=dp) :: dvol, vol ! volume element, volume
      REAL(KIND=dp) :: cutoff ! cutoff in a.u.
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: g ! grid point vectors
      REAL(KIND=dp), DIMENSION(:), POINTER :: gsq ! squared vector lengths
      INTEGER, DIMENSION(:, :), POINTER :: g_hat ! grid point indices (Miller)
      INTEGER, DIMENSION(:, :), POINTER :: g_hatmap ! mapped grid point indices (Miller) [CUDA]
      INTEGER :: grid_span ! type HALFSPACE/FULLSPACE
      LOGICAL :: have_g0 ! whether I have G = [0,0,0]
      INTEGER :: first_gne0 ! first g index /= 0 [1/2]
      INTEGER :: id_nr ! tag of this grid
      INTEGER :: reference ! reference grid identifier
      INTEGER, DIMENSION(:), POINTER :: gidx ! ref grid index
      INTEGER :: ref_count ! reference count
      LOGICAL :: spherical ! spherical cutoff?
      COMPLEX(KIND=dp), DIMENSION(:, :), POINTER :: grays ! used by parallel 3D FFT routine
   END TYPE pw_grid_type

END MODULE pw_grid_types

