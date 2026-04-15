# build.sh
mkdir -p __build__
# Remove previous files
rm *.exe

# Compilation
gfortran -c -fdefault-real-8 llt.f90
gfortran -c -fdefault-real-8 test_llt.f90
# gfortran -c -fdefault-real-8 test_program.f90

# Linking
gfortran llt.o test_llt.o -o __build__/test_program.exe

# Execution
./__build__/test_program.exe

# Remove unnecessary files
rm *.o
rm *.mod