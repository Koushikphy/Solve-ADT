# 3 state 2D adt solve

import numpy as np 
from scipy.interpolate import InterpolatedUnivariateSpline as spl
from scipy.integrate import solve_ivp
from numpy import sin,cos,tan

dat = np.loadtxt('./tauth_1d.dat')

grid = dat[:,0]

getTau1, getTau2, getTau3 = spl(grid, dat[:,1]), spl(grid, dat[:,2]), spl(grid, dat[:,3])


# for general 
def res(t, y): # returns value of the differential function at the point

    t12,t13,t23 = getTau1(t),getTau2(t),getTau3(t)
    y12 = -t12 - tan(y[1])*( t13*sin(y[0]) + t23*cos(y[0]))
    y13 = -t13*cos(y[0]) + t23*sin(y[0])
    y23 = -(1.0/cos(y[1]))*(t13*sin(y[0]) + t23*cos(y[0]))

    return y12,y13,y23


sol = solve_ivp(res, [grid[0],grid[-1]], [0,0,0], method='RK45', t_eval=grid, dense_output=True)
res = np.column_stack([sol.t, sol.y.T])

np.savetxt('tmp_th', res, delimiter='\t', fmt='%.8f')