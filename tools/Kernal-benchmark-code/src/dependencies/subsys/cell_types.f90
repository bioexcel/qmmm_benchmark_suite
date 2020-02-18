!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Handles all functions related to the CELL
!> \par History
!>      11.2008 Teodoro Laino [tlaino] - deeply cleaning cell_type from units
!>      10.2014 Moved many routines from cell_types.F here.
!> \author Matthias KracK (16.01.2002, based on a earlier version of CJM, JGH)
! **************************************************************************************************
MODULE cell_types
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: degree,&
                                              sqrt3
   USE mathlib,                         ONLY: angle,&
                                              det_3x3,&
                                              inv_3x3
   IMPLICIT NONE

   ! Impose cell symmetry
   INTEGER, PARAMETER, PUBLIC               :: cell_sym_none = 0, &
                                               cell_sym_triclinic = 1, &
                                               cell_sym_monoclinic = 2, &
                                               cell_sym_monoclinic_gamma_ab = 3, &
                                               cell_sym_orthorhombic = 4, &
                                               cell_sym_tetragonal_ab = 5, &
                                               cell_sym_tetragonal_ac = 6, &
                                               cell_sym_tetragonal_bc = 7, &
                                               cell_sym_rhombohedral = 8, &
                                               cell_sym_hexagonal = 9, &
                                               cell_sym_cubic = 10


! **************************************************************************************************
!> \brief   Type defining parameters related to the simulation cell
!> \version 1.0
! **************************************************************************************************
   TYPE cell_type
      INTEGER                           :: id_nr, ref_count, symmetry_id
      LOGICAL                           :: orthorhombic ! actually means a diagonal hmat
      REAL(KIND=dp)                     :: deth
      INTEGER, DIMENSION(3)             :: perd
      REAL(KIND=dp), DIMENSION(3, 3)    :: hmat, h_inv
   END TYPE cell_type

   TYPE cell_p_type
      TYPE(cell_type), POINTER :: cell
   END TYPE cell_p_type

   ! Public data types
   PUBLIC :: cell_type, cell_p_type

   ! Public subroutines
   PUBLIC :: get_cell, cell_create
             
   ! Public functions
   PUBLIC :: pbc, scaled_to_real

   INTERFACE pbc
      MODULE PROCEDURE pbc1, pbc2, pbc3, pbc4
   END INTERFACE

CONTAINS



! **************************************************************************************************
!> \brief   Get informations about a simulation cell.
!> \param cell ...
!> \param alpha ...
!> \param beta ...
!> \param gamma ...
!> \param deth ...
!> \param orthorhombic ...
!> \param abc ...
!> \param periodic ...
!> \param h ...
!> \param h_inv ...
!> \param id_nr ...
!> \param symmetry_id ...
!> \date    16.01.2002
!> \author  Matthias Krack
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE get_cell(cell, alpha, beta, gamma, deth, orthorhombic, abc, periodic, &
                       h, h_inv, id_nr, symmetry_id)

      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: alpha, beta, gamma, deth
      LOGICAL, INTENT(OUT), OPTIONAL                     :: orthorhombic
      REAL(KIND=dp), DIMENSION(3), INTENT(OUT), OPTIONAL :: abc
      INTEGER, DIMENSION(3), INTENT(OUT), OPTIONAL       :: periodic
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(OUT), &
         OPTIONAL                                        :: h, h_inv
      INTEGER, INTENT(out), OPTIONAL                     :: id_nr, symmetry_id

      IF (PRESENT(deth)) deth = cell%deth ! the volume
      IF (PRESENT(orthorhombic)) orthorhombic = cell%orthorhombic
      IF (PRESENT(periodic)) periodic(:) = cell%perd(:)
      IF (PRESENT(h)) h(:, :) = cell%hmat(:, :)
      IF (PRESENT(h_inv)) h_inv(:, :) = cell%h_inv(:, :)

      ! Calculate the lengths of the cell vectors a, b, and c
      IF (PRESENT(abc)) THEN
         abc(1) = SQRT(cell%hmat(1, 1)*cell%hmat(1, 1)+ &
                       cell%hmat(2, 1)*cell%hmat(2, 1)+ &
                       cell%hmat(3, 1)*cell%hmat(3, 1))
         abc(2) = SQRT(cell%hmat(1, 2)*cell%hmat(1, 2)+ &
                       cell%hmat(2, 2)*cell%hmat(2, 2)+ &
                       cell%hmat(3, 2)*cell%hmat(3, 2))
         abc(3) = SQRT(cell%hmat(1, 3)*cell%hmat(1, 3)+ &
                       cell%hmat(2, 3)*cell%hmat(2, 3)+ &
                       cell%hmat(3, 3)*cell%hmat(3, 3))
      END IF

      ! Angles between the cell vectors a, b, and c
      ! alpha = <(b,c)
      IF (PRESENT(alpha)) alpha = angle(cell%hmat(:, 2), cell%hmat(:, 3))*degree
      ! beta = <(a,c)
      IF (PRESENT(beta)) beta = angle(cell%hmat(:, 1), cell%hmat(:, 3))*degree
      ! gamma = <(a,b)
      IF (PRESENT(gamma)) gamma = angle(cell%hmat(:, 1), cell%hmat(:, 2))*degree
      IF (PRESENT(id_nr)) id_nr = cell%id_nr
      IF (PRESENT(symmetry_id)) symmetry_id = cell%symmetry_id

   END SUBROUTINE get_cell




! **************************************************************************************************
!> \brief   Apply the periodic boundary conditions defined by a simulation
!>          cell to a position vector r.
!> \param r ...
!> \param cell ...
!> \return ...
!> \date    16.01.2002
!> \author  Matthias Krack
!> \version 1.0
! **************************************************************************************************
   FUNCTION pbc1(r, cell) RESULT(r_pbc)

      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: r
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), DIMENSION(3)                        :: r_pbc

      REAL(KIND=dp), DIMENSION(3)                        :: s

      IF (cell%orthorhombic) THEN
         r_pbc(1) = r(1)-cell%hmat(1, 1)*cell%perd(1)*ANINT(cell%h_inv(1, 1)*r(1))
         r_pbc(2) = r(2)-cell%hmat(2, 2)*cell%perd(2)*ANINT(cell%h_inv(2, 2)*r(2))
         r_pbc(3) = r(3)-cell%hmat(3, 3)*cell%perd(3)*ANINT(cell%h_inv(3, 3)*r(3))
      ELSE
         s(1) = cell%h_inv(1, 1)*r(1)+cell%h_inv(1, 2)*r(2)+cell%h_inv(1, 3)*r(3)
         s(2) = cell%h_inv(2, 1)*r(1)+cell%h_inv(2, 2)*r(2)+cell%h_inv(2, 3)*r(3)
         s(3) = cell%h_inv(3, 1)*r(1)+cell%h_inv(3, 2)*r(2)+cell%h_inv(3, 3)*r(3)
         s(1) = s(1)-cell%perd(1)*ANINT(s(1))
         s(2) = s(2)-cell%perd(2)*ANINT(s(2))
         s(3) = s(3)-cell%perd(3)*ANINT(s(3))
         r_pbc(1) = cell%hmat(1, 1)*s(1)+cell%hmat(1, 2)*s(2)+cell%hmat(1, 3)*s(3)
         r_pbc(2) = cell%hmat(2, 1)*s(1)+cell%hmat(2, 2)*s(2)+cell%hmat(2, 3)*s(3)
         r_pbc(3) = cell%hmat(3, 1)*s(1)+cell%hmat(3, 2)*s(2)+cell%hmat(3, 3)*s(3)
      END IF

   END FUNCTION pbc1

! **************************************************************************************************
!> \brief   Apply the periodic boundary conditions defined by a simulation
!>          cell to a position vector r subtracting nl from the periodic images
!> \param r ...
!> \param cell ...
!> \param nl ...
!> \return ...
!> \date    16.01.2002
!> \author  Matthias Krack
!> \version 1.0
! **************************************************************************************************
   FUNCTION pbc2(r, cell, nl) RESULT(r_pbc)

      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: r
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, DIMENSION(3), INTENT(IN)                  :: nl
      REAL(KIND=dp), DIMENSION(3)                        :: r_pbc

      REAL(KIND=dp), DIMENSION(3)                        :: s

      IF (cell%orthorhombic) THEN
         r_pbc(1) = r(1)-cell%hmat(1, 1)*cell%perd(1)* &
                    REAL(NINT(cell%h_inv(1, 1)*r(1))-nl(1), dp)
         r_pbc(2) = r(2)-cell%hmat(2, 2)*cell%perd(2)* &
                    REAL(NINT(cell%h_inv(2, 2)*r(2))-nl(2), dp)
         r_pbc(3) = r(3)-cell%hmat(3, 3)*cell%perd(3)* &
                    REAL(NINT(cell%h_inv(3, 3)*r(3))-nl(3), dp)
      ELSE
         s(1) = cell%h_inv(1, 1)*r(1)+cell%h_inv(1, 2)*r(2)+cell%h_inv(1, 3)*r(3)
         s(2) = cell%h_inv(2, 1)*r(1)+cell%h_inv(2, 2)*r(2)+cell%h_inv(2, 3)*r(3)
         s(3) = cell%h_inv(3, 1)*r(1)+cell%h_inv(3, 2)*r(2)+cell%h_inv(3, 3)*r(3)
         s(1) = s(1)-cell%perd(1)*REAL(NINT(s(1))-nl(1), dp)
         s(2) = s(2)-cell%perd(2)*REAL(NINT(s(2))-nl(2), dp)
         s(3) = s(3)-cell%perd(3)*REAL(NINT(s(3))-nl(3), dp)
         r_pbc(1) = cell%hmat(1, 1)*s(1)+cell%hmat(1, 2)*s(2)+cell%hmat(1, 3)*s(3)
         r_pbc(2) = cell%hmat(2, 1)*s(1)+cell%hmat(2, 2)*s(2)+cell%hmat(2, 3)*s(3)
         r_pbc(3) = cell%hmat(3, 1)*s(1)+cell%hmat(3, 2)*s(2)+cell%hmat(3, 3)*s(3)
      END IF

   END FUNCTION pbc2

! **************************************************************************************************
!> \brief   Apply the periodic boundary conditions defined by the simulation
!>          cell cell to the vector pointing from atom a to atom b.
!> \param ra ...
!> \param rb ...
!> \param cell ...
!> \return ...
!> \date    11.03.2004
!> \author  Matthias Krack
!> \version 1.0
! **************************************************************************************************
   FUNCTION pbc3(ra, rb, cell) RESULT(rab_pbc)

      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: ra, rb
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), DIMENSION(3)                        :: rab_pbc

      INTEGER                                            :: icell, jcell, kcell
      INTEGER, DIMENSION(3)                              :: periodic
      REAL(KIND=dp)                                      :: rab2, rab2_pbc
      REAL(KIND=dp), DIMENSION(3)                        :: r, ra_pbc, rab, rb_image, rb_pbc, s2r

      CALL get_cell(cell=cell, periodic=periodic)

      ra_pbc(:) = pbc(ra(:), cell)
      rb_pbc(:) = pbc(rb(:), cell)

      rab2_pbc = HUGE(1.0_dp)

      DO icell = -periodic(1), periodic(1)
         DO jcell = -periodic(2), periodic(2)
            DO kcell = -periodic(3), periodic(3)
               r = REAL((/icell, jcell, kcell/), dp)
               CALL scaled_to_real(s2r, r, cell)
               rb_image(:) = rb_pbc(:)+s2r
               rab(:) = rb_image(:)-ra_pbc(:)
               rab2 = rab(1)*rab(1)+rab(2)*rab(2)+rab(3)*rab(3)
               IF (rab2 < rab2_pbc) THEN
                  rab2_pbc = rab2
                  rab_pbc(:) = rab(:)
               END IF
            END DO
         END DO
      END DO

   END FUNCTION pbc3

   !if positive_range == true, r(i) (or s(i)) in range [0, hmat(i,i)],
   !else, r(i) (s(i)) in range [-hmat(i,i)/2, hmat(i,i)/2]
! **************************************************************************************************
!> \brief ...
!> \param r ...
!> \param cell ...
!> \param positive_range ...
!> \return ...
! **************************************************************************************************
   FUNCTION pbc4(r, cell, positive_range) RESULT(r_pbc)

      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: r
      TYPE(cell_type), POINTER                           :: cell
      LOGICAL                                            :: positive_range
      REAL(KIND=dp), DIMENSION(3)                        :: r_pbc

      REAL(KIND=dp), DIMENSION(3)                        :: s

      IF (positive_range) THEN
         IF (cell%orthorhombic) THEN
            r_pbc(1) = r(1)-cell%hmat(1, 1)*cell%perd(1)*FLOOR(cell%h_inv(1, 1)*r(1))
            r_pbc(2) = r(2)-cell%hmat(2, 2)*cell%perd(2)*FLOOR(cell%h_inv(2, 2)*r(2))
            r_pbc(3) = r(3)-cell%hmat(3, 3)*cell%perd(3)*FLOOR(cell%h_inv(3, 3)*r(3))
         ELSE
            s(1) = cell%h_inv(1, 1)*r(1)+cell%h_inv(1, 2)*r(2)+cell%h_inv(1, 3)*r(3)
            s(2) = cell%h_inv(2, 1)*r(1)+cell%h_inv(2, 2)*r(2)+cell%h_inv(2, 3)*r(3)
            s(3) = cell%h_inv(3, 1)*r(1)+cell%h_inv(3, 2)*r(2)+cell%h_inv(3, 3)*r(3)
            s(1) = s(1)-cell%perd(1)*FLOOR(s(1))
            s(2) = s(2)-cell%perd(2)*FLOOR(s(2))
            s(3) = s(3)-cell%perd(3)*FLOOR(s(3))
            r_pbc(1) = cell%hmat(1, 1)*s(1)+cell%hmat(1, 2)*s(2)+cell%hmat(1, 3)*s(3)
            r_pbc(2) = cell%hmat(2, 1)*s(1)+cell%hmat(2, 2)*s(2)+cell%hmat(2, 3)*s(3)
            r_pbc(3) = cell%hmat(3, 1)*s(1)+cell%hmat(3, 2)*s(2)+cell%hmat(3, 3)*s(3)
         END IF
      ELSE
         r_pbc = pbc1(r, cell)
      END IF

   END FUNCTION pbc4

! **************************************************************************************************
!> \brief   Transform scaled cell coordinates real coordinates.
!>          r=h*s
!> \param r ...
!> \param s ...
!> \param cell ...
!> \date    16.01.2002
!> \author  Matthias Krack
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE scaled_to_real(r, s, cell)

      REAL(KIND=dp), DIMENSION(3), INTENT(OUT)           :: r
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: s
      TYPE(cell_type), POINTER                           :: cell

      IF (cell%orthorhombic) THEN
         r(1) = cell%hmat(1, 1)*s(1)
         r(2) = cell%hmat(2, 2)*s(2)
         r(3) = cell%hmat(3, 3)*s(3)
      ELSE
         r(1) = cell%hmat(1, 1)*s(1)+cell%hmat(1, 2)*s(2)+cell%hmat(1, 3)*s(3)
         r(2) = cell%hmat(2, 1)*s(1)+cell%hmat(2, 2)*s(2)+cell%hmat(2, 3)*s(3)
         r(3) = cell%hmat(3, 1)*s(1)+cell%hmat(3, 2)*s(2)+cell%hmat(3, 3)*s(3)
      END IF

   END SUBROUTINE scaled_to_real

! **************************************************************************************************
!> \brief allocates and initializes a cell
!> \param cell the cell to initialize
!> \param hmat the h matrix that defines the cell
!> \param periodic periodicity of the cell
!> \par History
!>      09.2003 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE cell_create(cell,&
        id_nr, ref_count, symmetry_id,&
        orthorhombic,&
        deth, hmat, h_inv, periodic)

     TYPE(cell_type), POINTER                            :: cell
     INTEGER                                             :: id_nr, ref_count, symmetry_id
     LOGICAL                                             :: orthorhombic
     REAL(KIND=dp)                                       :: deth
     INTEGER, DIMENSION(3), INTENT(IN)                   :: periodic
     REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)          :: hmat, h_inv

     ALLOCATE (cell)
     cell%id_nr = id_nr
     cell%ref_count = ref_count
     cell%symmetry_id = symmetry_id
     cell%orthorhombic = orthorhombic
     cell%deth = deth
     cell%perd = periodic
     cell%hmat = hmat
     cell%h_inv = h_inv
     
   END SUBROUTINE cell_create


END MODULE cell_types
