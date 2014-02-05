!> \file: SimulationVars.F90
!> \author: Sayop Kim

MODULE SimulationVars_m
   USE parameters_m, ONLY : wp
   IMPLICIT NONE

   INTEGER :: imax, jmax, kmax, nmax
   REAL(KIND=wp), ALLOCATABLE, DIMENSION(:,:,:,:) :: xp
   REAL(KIND=wp), ALLOCATABLE, DIMENSION(:,:) :: BOTedge
   REAL(KIND=wp) :: cy
   REAL(KIND=wp), ALLOCATABLE, DIMENSION(:,:,:) :: inverseJacobian
   REAL(KIND=wp), DIMENSION(3,8) :: xblkV    ! x,y,z points at 8 vertices of block
   REAL(KIND=wp) :: error, ErrMax
END MODULE SimulationVars_m
