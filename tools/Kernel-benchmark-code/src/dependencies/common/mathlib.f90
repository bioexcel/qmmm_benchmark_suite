!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Collection of simple mathematical functions and subroutines
!> \par History
!>      FUNCTION angle updated and FUNCTION dihedral angle added; cleaned
!>      (13.03.2004,MK)
!> \author MK (15.11.1998)
! **************************************************************************************************
MODULE mathlib

   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE mathconstants,                   ONLY: euler,&
                                              fac
   IMPLICIT NONE

   PRIVATE

   REAL(KIND=dp), PARAMETER             :: eps_geo = 1.0E-6_dp

   ! Public functions

   PUBLIC :: angle, &
             det_3x3, &
             inv_3x3

   INTERFACE det_3x3
      MODULE PROCEDURE det_3x3_1, det_3x3_2
   END INTERFACE


CONTAINS


! **************************************************************************************************
!> \brief  Calculation of the angle between the vectors a and b.
!>         The angle is returned in radians.
!> \param a ...
!> \param b ...
!> \return ...
!> \date    14.10.1998
!> \author  MK
!> \version 1.0
! **************************************************************************************************
   FUNCTION angle(a, b) RESULT(angle_ab)

      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: a, b
      REAL(KIND=dp)                                      :: angle_ab

      REAL(KIND=dp)                                      :: length_of_a, length_of_b
      REAL(KIND=dp), DIMENSION(SIZE(a, 1))               :: a_norm, b_norm

      length_of_a = SQRT(DOT_PRODUCT(a, a))
      length_of_b = SQRT(DOT_PRODUCT(b, b))

      IF ((length_of_a > eps_geo) .AND. (length_of_b > eps_geo)) THEN
         a_norm(:) = a(:)/length_of_a
         b_norm(:) = b(:)/length_of_b
         angle_ab = ACOS(MIN(MAX(DOT_PRODUCT(a_norm, b_norm), -1.0_dp), 1.0_dp))
      ELSE
         angle_ab = 0.0_dp
      END IF

   END FUNCTION angle

! **************************************************************************************************
!> \brief   Returns the determinante of the 3x3 matrix a.
!> \param a ...
!> \return ...
!> \date    13.03.2004
!> \author  MK
!> \version 1.0
! **************************************************************************************************
   FUNCTION det_3x3_1(a) RESULT(det_a)
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: a
      REAL(KIND=dp)                                      :: det_a

      det_a = a(1, 1)*(a(2, 2)*a(3, 3)-a(2, 3)*a(3, 2))+ &
              a(1, 2)*(a(2, 3)*a(3, 1)-a(2, 1)*a(3, 3))+ &
              a(1, 3)*(a(2, 1)*a(3, 2)-a(2, 2)*a(3, 1))

   END FUNCTION det_3x3_1

! **************************************************************************************************
!> \brief   Returns the determinante of the 3x3 matrix a given by its columns.
!> \param a1 ...
!> \param a2 ...
!> \param a3 ...
!> \return ...
!> \date    13.03.2004
!> \author  MK
!> \version 1.0
! **************************************************************************************************
   FUNCTION det_3x3_2(a1, a2, a3) RESULT(det_a)
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: a1, a2, a3
      REAL(KIND=dp)                                      :: det_a

      det_a = a1(1)*(a2(2)*a3(3)-a3(2)*a2(3))+ &
              a2(1)*(a3(2)*a1(3)-a1(2)*a3(3))+ &
              a3(1)*(a1(2)*a2(3)-a2(2)*a1(3))

   END FUNCTION det_3x3_2


! **************************************************************************************************
!> \brief   Returns the inverse of the 3 x 3 matrix a.
!> \param a ...
!> \return ...
!> \date    13.03.2004
!> \author  MK
!> \version 1.0
! **************************************************************************************************
   FUNCTION inv_3x3(a) RESULT(a_inv)

      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: a
      REAL(KIND=dp), DIMENSION(3, 3)                     :: a_inv

      REAL(KIND=dp)                                      :: det_a

      det_a = 1.0_dp/det_3x3(a)

      a_inv(1, 1) = (a(2, 2)*a(3, 3)-a(3, 2)*a(2, 3))*det_a
      a_inv(2, 1) = (a(2, 3)*a(3, 1)-a(3, 3)*a(2, 1))*det_a
      a_inv(3, 1) = (a(2, 1)*a(3, 2)-a(3, 1)*a(2, 2))*det_a

      a_inv(1, 2) = (a(1, 3)*a(3, 2)-a(3, 3)*a(1, 2))*det_a
      a_inv(2, 2) = (a(1, 1)*a(3, 3)-a(3, 1)*a(1, 3))*det_a
      a_inv(3, 2) = (a(1, 2)*a(3, 1)-a(3, 2)*a(1, 1))*det_a

      a_inv(1, 3) = (a(1, 2)*a(2, 3)-a(2, 2)*a(1, 3))*det_a
      a_inv(2, 3) = (a(1, 3)*a(2, 1)-a(2, 3)*a(1, 1))*det_a
      a_inv(3, 3) = (a(1, 1)*a(2, 2)-a(2, 1)*a(1, 2))*det_a

   END FUNCTION inv_3x3

END MODULE mathlib
