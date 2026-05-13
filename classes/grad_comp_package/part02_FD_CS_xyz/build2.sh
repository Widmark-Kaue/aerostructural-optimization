# build2.sh
# Builds and executes the code

# Compile
gfortran -c -fdefault-real-8 mymodule2.f90 main2.f90

# Link
gfortran mymodule2.o main2.o -o main2.exe

# Execute
./main2.exe

# Clean up
rm *.o
rm *.mod
