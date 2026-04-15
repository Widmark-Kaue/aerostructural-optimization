PROG=fortran_var

rm *.o
rm *.exe


gfortran -c -fdefault-real-8 $PROG.F90

gfortran $PROG.o -o $PROG.exe

 ./$PROG.exe
