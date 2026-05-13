# build2.sh
# Builds and executes the code

FFLAGS="-fdefault-real-8 -fallow-argument-mismatch"

# Compile most used files from the ADFirstAidKit
gcc -c ADFirstAidKit/adStack.c
gfortran -c $FFLAGS ADFirstAidKit/adBuffer.f

# Compile
gfortran -c $FFLAGS mymodule2.f90
gfortran -c $FFLAGS TapenadeResults_d/mymodule2_d.f90
gfortran -c $FFLAGS TapenadeResults_b/mymodule2_b.f90
gfortran -c $FFLAGS main2.f90

# Link
gfortran adStack.o adBuffer.o mymodule2_d.o mymodule2_b.o mymodule2.o main2.o -o main2.exe

# Execute
./main2.exe

# Clean up
rm *.o
rm *.mod
