# build.sh
# Builds and executes the code

# Compile
# gfortran -c -fdefault-real-8 mymodule.f90 main.f90
gfortran -c -fdefault-real-8 mymodule_if.f90 main.f90

# Link
# gfortran mymodule.o main.o -o main.exe
gfortran mymodule_if.o main.o -o main.exe

# Execute
./main.exe

# Clean up
rm *.o
rm *.mod
