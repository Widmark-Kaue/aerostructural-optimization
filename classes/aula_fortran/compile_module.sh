rm *.o
rm *.exe

gfortran -c -fdefault-real-8 mymodule.F90 usemodule.F90

gfortran mymodule.o usemodule.o -o usemodule.exe

./usemodule.exe
