### Windows10
## compile as a single command
#FC = mpif90 -fopenmp -cpp -DMPI
#FFLAGS = -O2 -Wall
#
## compile as a part of GUI
#FC = i686-w64-mingw32-gfortran -fopenmp -cpp
#FFLAGS = -O2 -Wall

### Linux
## compile as a single command(using Intel fortran compiler)
#FC = mpif90 -qopenmp -cpp -DMPI
#FFLAGS = -O2 -warn
#
## compile as a part of GUI
#FC = gfortran -fopenmp -cpp -g -ffpe-trap=zero,invalid,overflow,underflow
#FFLAGS = -O2 -Wall

### Mac
#FC = /usr/local/Cellar/gcc/11.1.0_1/bin/gfortran -fopenmp -cpp -fPIC
FFLAGS = -O2 -Wall

# as single command
#FC = mpif90 -fopenmp -cpp -fPIC -DMPI
#FFLAGS = -O2 -Wall

# # part of GUI
FC = gfortran -fopenmp -cpp -fPIC -g -ffpe-trap=zero,invalid,overflow,underflow 
FFLAGS = -Og -Wall  -fbacktrace

