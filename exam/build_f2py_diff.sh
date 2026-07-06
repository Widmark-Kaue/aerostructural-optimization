# build_f2py_diff_llt.sh
# This script builds the Python interface for llt.f90
# The differentiated modes are also included here

# Define variables
MOD=asa_module
MOD1=fem_module
MOD2=llt_module
FFLAGS="-O2 -fPIC -fdefault-real-8 -fbounds-check -fallow-argument-mismatch"
CFLAGS="-fPIC"

# Remove previous files
rm *.o
rm *.mod
rm *.so


# Differentiate the codes
rm -rf TapenadeResults_d
rm -rf TapenadeResults_b
python3 send_to_tapenade.py f
python3 send_to_tapenade.py r

# Make a copy of the file that maps Fortran types to Python types
cp f2py_f2cmap .f2py_f2cmap

# BUILD ORIGINAL CODE
f2py -c -m ${MOD} ${MOD1}.f90 ${MOD2}.f90 ${MOD}.f90 --f90flags="$FFLAGS" --f77flags="$FFLAGS" -DF2PY_REPORT_ON_ARRAY_COPY=1

# BUILD FORWARD AD
f2py -c -m ${MOD}_d TapenadeResults_d/${MOD1}_d.f90 TapenadeResults_d/${MOD2}_d.f90 TapenadeResults_d/${MOD}_d.f90 --f90flags="$FFLAGS" --f77flags="$FFLAGS" -DF2PY_REPORT_ON_ARRAY_COPY=1

# BUILD REVERSE AD
f2py -c -m ${MOD}_b ADFirstAidKit/adStack.c ADFirstAidKit/adBuffer.f TapenadeResults_b/${MOD1}_b.f90 TapenadeResults_b/${MOD2}_b.f90 TapenadeResults_b/${MOD}_b.f90 --f90flags="$FFLAGS" --f77flags="$FFLAGS" -DF2PY_REPORT_ON_ARRAY_COPY=1

# Clean-up
rm *.o
rm *.mod
rm .f2py_f2cmap
