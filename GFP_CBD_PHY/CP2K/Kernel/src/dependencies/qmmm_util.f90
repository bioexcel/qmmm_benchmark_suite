!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \par History
!>      09.2004 created [tlaino]
!> \author Teodoro Laino
! **************************************************************************************************
MODULE qmmm_util
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .FALSE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qmmm_util'
   PUBLIC :: spherical_cutoff_factor

CONTAINS

! **************************************************************************************************
!> \brief Computes a spherical cutoff factor for the QMMM interactions
!> \param spherical_cutoff ...
!> \param rij ...
!> \param factor ...
!> \par History
!>      08.2008 created
!> \author Teodoro Laino
! **************************************************************************************************
   SUBROUTINE spherical_cutoff_factor(spherical_cutoff, rij, factor)
      REAL(KIND=dp), DIMENSION(2), INTENT(IN)            :: spherical_cutoff
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rij
      REAL(KIND=dp), INTENT(OUT)                         :: factor

      CHARACTER(len=*), PARAMETER :: routineN = 'spherical_cutoff_factor', &
         routineP = moduleN//':'//routineN

      REAL(KIND=dp)                                      :: r, r0

      r = SQRT(DOT_PRODUCT(rij, rij))
      r0 = spherical_cutoff(1)-20.0_dp*spherical_cutoff(2)
      factor = 0.5_dp*(1.0_dp-TANH((r-r0)/spherical_cutoff(2)))

   END SUBROUTINE spherical_cutoff_factor

END MODULE qmmm_util
