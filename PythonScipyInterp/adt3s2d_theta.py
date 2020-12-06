# 3 state 2D adt solve, preferably use python 3.x
# import sys
import numpy as np 
from numpy import sin,cos,tan
from scipy.integrate import solve_ivp
from scipy.interpolate import RectBivariateSpline as rbs
import mycodes as mc
import os


class adtSolver3S2D():
    def __init__(self, theGrid, phiGrid, tauth, tauph ):
        self.theGrid = theGrid
        self.phiGrid = phiGrid
        self.tauTheFuncs = [rbs(theGrid, phiGrid,tauth[...,i]) for i in [0,1,2]]
        self.tauPhiFuncs = [rbs(theGrid, phiGrid,tauph[...,i]) for i in [0,1,2]]


    def resThe(self,theta,y,phi): # along theta
        tau = [f(theta,phi)[0,0] for f in self.tauTheFuncs]
        # tau = np.clip(tau, 0, np.inf)   # remove any negetive values during interplation
        return self.diffSys(y,tau)


    def resPhi(self,phi,y,theta): # along Phi
        tau = [f(theta,phi)[0,0] for f in self.tauPhiFuncs]
        # tau = np.clip(tau, 0, np.inf)   # remove any negetive values during interplation
        return self.diffSys(y,tau)


    def diffSys(self,y,t):
        # A = A12*A13*A23 #
        y12 = -t[0] - tan(y[1])*( t[1]*sin(y[0]) + t[2]*cos(y[0]))
        y13 = -t[1]*cos(y[0]) + t[2]*sin(y[0])
        y23 = -(1.0/cos(y[1]))*(t[1]*sin(y[0]) + t[2]*cos(y[0]))
        return y12,y13,y23


    def solve(self, thGrid=None, phGrid=None, resFile='angle.dat', initZero=False):
        if thGrid is None: thGrid = self.theGrid
        if phGrid is None: phGrid = self.phiGrid

        if initZero:
            iPhVal = np.full((thGrid.shape[0],3),0)
        else:
            sol    = solve_ivp(self.resThe, [thGrid[0],thGrid[-1]], [0,0,0], method='LSODA', t_eval=thGrid, args=(phGrid[0],))
            iPhVal = sol.y.T                                      #<<<--- initial values for integration along the theta grid

        res  = np.empty((thGrid.shape[0], phGrid.shape[0], 3))
        for i,th in enumerate(thGrid):
            sol = solve_ivp(self.resPhi, [phGrid[0],phGrid[-1]], iPhVal[i], method='LSODA', t_eval=phGrid, args=(th,))
            res[i,...] = sol.y.T

        xx,yy = np.meshgrid(thGrid, phGrid, indexing='ij')
        resToWrite = np.dstack([xx,yy,res])

        with open(resFile,'w') as f:
            for i in resToWrite:
                np.savetxt(f,i ,delimiter='\t',fmt='%.8f')
                f.write('\n')


# # read file
# tauth = np.loadtxt('./TauThetaNACT_new.dat', usecols=(1,2,3,4,5))
# tauph = np.loadtxt('./TauPhiNACT_new.dat',   usecols=(1,2,3,4,5))
# # mirror data
# tauth = mc.mirror(tauth, shape=(66,181,5))
# tauph = mc.mirror(tauph, shape=(66,181,5))
# # shape data
# tauth.shape = (66,361,5)
# tauph.shape = (66,361,5)
# # grid from data
# gridTh = np.deg2rad(tauth[:,0,0])
# gridPh = np.deg2rad(tauth[0,:,1])
# # nacts
# tauth = tauth[:,:,[2,3,4]]
# tauph = tauph[:,:,[2,3,4]]
# # new grid to solve to
# thGrid = np.deg2rad(np.linspace(0,90,91))
# phGrid = np.deg2rad(np.linspace(0,360,361))


# ss = adtSolver3S2D(gridTh, gridPh, tauth, tauth)   #passing dummy `tauph`, as its not required anyway
# ss.solve(thGrid=thGrid, phGrid=phGrid, resFile='angle2.dat')

os.chdir("./RHO_1.0")

# read file
tauth = np.loadtxt('./tauth.dat')
tauph = np.loadtxt('./tauph.dat')
# mirror data
# shape data
tauth.shape = (91,361,5)
tauph.shape = (91,361,5)
# grid from data
gridTh = tauth[:,0,0]
gridPh = tauth[0,:,1]
# nacts
tauth = tauth[:,:,[2,3,4]]
tauph = tauph[:,:,[2,3,4]]


ss = adtSolver3S2D(gridTh, gridPh, tauth, tauph)   #passing dummy `tauph`, as its not required anyway
ss.solve(resFile='angle2.dat')