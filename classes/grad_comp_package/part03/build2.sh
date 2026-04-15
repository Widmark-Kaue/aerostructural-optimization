# build2.sh
# Builds and executes the code

# Compile most used files from the ADFirstAidKit
gcc -c ADFirstAidKit/adStack.c
gfortran -c -fdefault-real-8 ADFirstAidKit/adBuffer.f

# Compile
gfortran -c -fdefault-real-8 mymodule2.f90
gfortran -c -fdefault-real-8 TapenadeResults_d/mymodule2_d.f90
gfortran -c -fdefault-real-8 TapenadeResults_b/mymodule2_b.f90
gfortran -c -fdefault-real-8 mymodule2.f90
gfortran -c -fdefault-real-8 main2.f90

# Link
gfortran adStack.o adBuffer.o mymodule2_d.o mymodule2_b.o mymodule2.o main2.o -o main2.exe

# Execute
./main2.exe

# Clean up
rm *.o
rm *.mod
