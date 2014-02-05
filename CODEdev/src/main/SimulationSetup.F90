!> \file SimulationSetup.F90
!> \author Sayop Kim

MODULE SimulationSetup_m
   USE Parameters_m, ONLY: wp
   IMPLICIT NONE

   PUBLIC :: InitializeCommunication, UniformSpacing

CONTAINS

!-----------------------------------------------------------------------------!
   SUBROUTINE InitializeCommunication()
!-----------------------------------------------------------------------------!
      USE Parameters_m, ONLY: CODE_VER_STRING
      IMPLICIT NONE
     
      WRITE(*,'(a)') "" 
      WRITE(*,'(a)') "CFD code Version: ", CODE_VER_STRING
   END SUBROUTINE InitializeCommunication


!-----------------------------------------------------------------------------!
   FUNCTION UniformSpacing(xmin,xmax,indx,indxMax) RESULT(outcome)
!-----------------------------------------------------------------------------!
      !Distribute interior grid points based on edge points' coordinates.
      !Linear Interpolateion is made by referring to (i,j,k) indices
      IMPLICIT NONE
      REAL(KIND=wp), INTENT(IN) :: xmin, xmax
      INTEGER, INTENT(IN) :: indx, indxMax
      REAL(KIND=wp) :: outcome, coef
      coef = REAL(indx - 1) / REAL(indxMax - 1)
      outcome = xmin + coef * (xmax - xmin)
   END FUNCTION UniformSpacing

END MODULE SimulationSetup_m
