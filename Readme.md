3 state explicit ADT equation is solved along theta grid in three cases:
1. Fortran with `DLSODA` subroutine from `ODEPACK`
2. Python with interpolation and `LSODA` differential solver from scipy
3. Python with `LSODA` differential solver from scipy and Fortran interpolation routine used in "1" compiled with f2py

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

cd PythonFortInterp
f2py -c interpolfort.f90 -m fInterp
python3 adt_3s2d_f2py.py

```
