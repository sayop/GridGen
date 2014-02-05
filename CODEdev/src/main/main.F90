!> \file: main.F90
!> \author: Sayop Kim

PROGRAM main
   USE SimulationSetup_m, ONLY: InitializeCommunication
   USE GridSetup_m, ONLY: InitializeGrid
   USE GridTransform_m, ONLY: GridTransform
   USE io_m, ONLY: WriteTecPlot, filenameLength
   USE Parameters_m, ONLY: wp

   IMPLICIT NONE

   CHARACTER(LEN=filenameLength) :: outputfile = 'output.tec'

   CALL InitializeCommunication
   ! Make initial condition for grid point alignment
   ! Using Algebraic method
   CALL InitializeGrid
   ! Use Elliptic grid points
   CALL GridTransform
   CALL WriteTecPlot(outputfile,'"I","J","K","Jacobian"')
END PROGRAM main
