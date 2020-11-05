3 state explicit ADT equation is solved along theta grid using:
1. Fortran with `DLSODA` subroutine from `ODEPACK`
2. Python with interpolation and `LSODA` differential solver from scipy.

Execute everything:
```bash

cd Fortran 
gfortran -c *.f
gfortran adt_angle.f90 *.o
./a.out
cd ..

cd PythonScipyInterp
python3 adt_3s2d.py
cd ..
```
