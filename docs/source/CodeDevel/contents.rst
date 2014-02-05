Code development
================

1. GridGen Code summary
---------------------------

The present project is to make a grid-generator for 3-D computational domain around a modified NACA 00xx series airfoil in a channel. The assigned project is inherently aimed at 2-D grid. However, the currently built GridGen code has a capability of 3-D grid generation.

The source code contains two directories, 'io', and 'main', for input/output related sources and grid-setup related sources, respectively. 'CMakeLists.txt' file is also included for cmake compiling.

::

   $ cd GridGen/CODEdev/src/
   $ ls
   $ CMakeLists.txt  io  main

The **io** folder has **io.F90** file which contains **ReadGridInput()** and **WriteTecPlot()** subroutines. It also includes **input** directory which contains default **input.dat** file.

The **main** folder is only used for containing grid-setup related source files. The main routine is run by **main.F90** which calls important subroutines from **main** folder itself and **io** folder when needed. All the fortran source files **main** folder contains are listed below::

   > GridSetup.F90
   > GridTransform.F90
   > GridTransformSetup.F90
   > main.F90
   > Parameters.F90
   > SimulationSetup.F90
   > SimulationVars.F90


2. Code details
---------------

The source code shown below is **main.F90** and it calls skeletal subroutines for generating grid structure. The main features of the main code is to make initialized variable arrays and read input and write output files::

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

The **main.F90** file first refers to **InitializeGrid** subroutine defined in **GridSetup.F90** file. The main function of this routine is to call again multiple subroutines defined in same file. The subroutine definition shown below summarizes the how the code runs for the grid initialization::

  !-----------------------------------------------------------------------------!
     SUBROUTINE InitializeGrid()
  !-----------------------------------------------------------------------------!
     USE io_m, ONLY: ReadGridInput
     USE SimulationVars_m, ONLY: imax, jmax, kmax,&
                               xblkV, cy
     IMPLICIT NONE
 
     ! Create Bottom Edge coordinate values
     CALL ReadGridInput
     CALL InitializeGridArrays
     CALL CreateBottomEdge
     CALL SetEdgePnts
     CALL GridPntsAlgbra
     CALL GenerateInteriorPoints

     END SUBROUTINE

* **CALL ReadGridInput**: Reads important user defined variables and parameters for grid configuration.

* **CALL InitializeGridArrays**: Initialize the single- and multi-dimensional arrays and set their size with input parameters(for example, imax, jmax, kmax).

* **CALL CreateBottomEdge**: Generate point values for airfoil geometry.

* **CALL SetEdgePnts**: Generate grid points along 8 edges of the computational domain.

* **CALL GridPntsAlgbra**: Based on the edge points, this routine will distribute grid points located on each 6 surfaces of the computational domain.

* **CALL GenerateInteriorPoints**: Based on grid points along the edges and surfaces, this routine will create interior grid points that are aligned with user-defined grid point interpolations.



