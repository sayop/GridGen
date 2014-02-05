How to run the code
===================

1. Code setup
-------------

The GridGen source code has been developed with version management tool, GIT. The git repository was built on 'github.com'. Thus, the source code as well as related document files can be cloned into user's local machine by following command::

   $ git clone http://github.com/sayop/GridGen.git

If you open the git-cloned folder **GridGen**, you will see two different folder and README file. The **CODEdev** folder contains again **bin** folder and **src** folder. In order to run the code, user should run **setup.sh** script in the **bin** folder. It may contain **build** folder, which might have been created in the different platform. Thus it is recommended that user should remove **build** folder before setting up the code. Note that the **setup.sh** script will run **cmake** command. Thus, make sure to have cmake installed on your system::

  $ rm -rf build
  $ ./setup.sh
  -- The C compiler identification is GNU 4.6.3
  -- The CXX compiler identification is GNU 4.6.3
  -- Check for working C compiler: /usr/bin/gcc
  -- Check for working C compiler: /usr/bin/gcc -- works
  -- Detecting C compiler ABI info
  -- Detecting C compiler ABI info - done
  -- Check for working CXX compiler: /usr/bin/c++
  -- Check for working CXX compiler: /usr/bin/c++ -- works
  -- Detecting CXX compiler ABI info
  -- Detecting CXX compiler ABI info - done
  -- The Fortran compiler identification is Intel
  -- Check for working Fortran compiler: /opt/intel/composer_xe_2011_sp1.11.339/bin/intel64/ifort
  -- Check for working Fortran compiler: /opt/intel/composer_xe_2011_sp1.11.339/bin/intel64/ifort  -- works
  -- Detecting Fortran compiler ABI info
  -- Detecting Fortran compiler ABI info - done
  -- Checking whether /opt/intel/composer_xe_2011_sp1.11.339/bin/intel64/ifort supports Fortran 90
  -- Checking whether /opt/intel/composer_xe_2011_sp1.11.339/bin/intel64/ifort supports Fortran 90 -- yes
  -- Configuring done
  -- Generating done
  -- Build files have been written to: /data/ksayop/GitHub.Clone/GridGen/CODEdev/bin/build
  Scanning dependencies of target cfd.x
  [ 12%] Building Fortran object CMakeFiles/cfd.x.dir/main/Parameters.F90.o
  [ 25%] Building Fortran object CMakeFiles/cfd.x.dir/main/SimulationVars.F90.o
  [ 37%] Building Fortran object CMakeFiles/cfd.x.dir/io/io.F90.o
  [ 50%] Building Fortran object CMakeFiles/cfd.x.dir/main/SimulationSetup.F90.o
  [ 62%] Building Fortran object CMakeFiles/cfd.x.dir/main/GridSetup.F90.o
  [ 75%] Building Fortran object CMakeFiles/cfd.x.dir/main/GridTransformSetup.F90.o
  [ 87%] Building Fortran object CMakeFiles/cfd.x.dir/main/GridTransform.F90.o
  [100%] Building Fortran object CMakeFiles/cfd.x.dir/main/main.F90.o
  Linking Fortran executable cfd.x
  [100%] Built target cfd.x
  $ ls
  $ build  cfd.x  input.dat  setup.sh

If you run this, you will get executable named **cfd.x** and **input.dat** files. The input file is made by default. You can quickly change the required options.


2. Input file setup
-------------------

Under construction
