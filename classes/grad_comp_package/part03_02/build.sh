# build.sh
# Builds and executes the code

FFLAGS="-fdefault-real-8 -fallow-argument-mismatch"

# Compile
gfortran -c $FFLAGS mymodule2.f90
gfortran -c $FFLAGS main2.f90

# Link
gfortran mymodule2.o main2.o -o main2.exe

# Execute
./main2.exe

# Clean up
rm *.o
rm *.mod
