# build_f2py.sh
# This script builds the Python interface for llt.f90

# Define variables
MOD=llt
FFLAGS="-O2 -fPIC -fdefault-real-8 -fbounds-check -fallow-argument-mismatch"

# Remove previous files
rm *.o
rm *.mod
rm *.so

# Make a copy of the file that maps Fortran types to Python types
cp __build__/f2py_f2cmap .f2py_f2cmap

# Generate the Python module
f2py -c -m ${MOD}_f90 ${MOD}.f90 --f90flags="$FFLAGS" -DF2PY_REPORT_ON_ARRAY_COPY=1

# Clean-up
rm *.o
rm *.mod
rm .f2py_f2cmap
rm *.pyf