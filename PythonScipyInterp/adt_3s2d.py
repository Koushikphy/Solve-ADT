# 3 state 2D adt solve, preferably use python 3.x
import sys
import numpy as np 
from numpy import sin,cos,tan
from scipy.integrate import solve_ivp
from scipy.interpolate import RectBivariateSpline as rbs


class adtSolver3S2D():

    def __init__(self, theGrid, phiGrid, tauth, tauph ):
        self.theGrid = theGrid
        self.phiGrid = phiGrid
        self.tauTheFuncs = [rbs(theGrid, phiGrid,tauth[...,i]) for i in [0,1,2]]
        self.tauPhiFuncs = [rbs(theGrid, phiGrid,tauph[...,i]) for i in [0,1,2]]


    def resThe(self,theta,y,phi): # along theta
        tau = [f(theta,phi)[0,0] for f in self.tauTheFuncs]
        return self.diffSys(y,tau)


    def resPhi(self,phi,y,theta): # along Phi
        tau = [f(theta,phi)[0,0] for f in self.tauPhiFuncs]
        return self.diffSys(y,tau)


    def diffSys(self,y,t):
        # A = A12*A13*A23 #
        y12 = -t[0] - tan(y[1])*( t[1]*sin(y[0]) + t[2]*cos(y[0]))
        y13 = -t[1]*cos(y[0]) + t[2]*sin(y[0])
        y23 = -(1.0/cos(y[1]))*(t[1]*sin(y[0]) + t[2]*cos(y[0]))
        return y12,y13,y23


    def solve(self, thGrid=None, phGrid=None, resFile='angle.dat', initZero=False):

        # if `thGrid` and `phGrid` is provided its used for solution
        # if `initZero` is True, then don't do the initial phi integration and start along theta integration with 0, making it symmetric

        if thGrid is None: thGrid = self.theGrid
        if phGrid is None: phGrid = self.phiGrid

        if initZero:
            iThVal  = np.full((phGrid.shape[0],3),0)
        else: # first solve along phi grid for the 1st theta value
            sol     = solve_ivp(self.resPhi, [phGrid[0],phGrid[-1]], [0,0,0], method='LSODA', t_eval=phGrid, args=(thGrid[0],)) #<- initial value [0,0,0]
            iThVal  = sol.y.T  #<<<--- initial values for integration along the theta grid

        res  = np.empty((thGrid.shape[0], phGrid.shape[0], 3))
        for i,phi in enumerate(phGrid):
            sol = solve_ivp(self.resThe, [thGrid[0],thGrid[-1]], iThVal[i], method='LSODA', t_eval=thGrid, args=(phi,))#, rtol=1.0e-8, atol=1.0e-10)
            res[:,i,:] = sol.y.T

        xx,yy = np.meshgrid(thGrid, phGrid, indexing='ij')
        resToWrite = np.dstack([xx,yy,res])

        with open(resFile,'w') as f:
            for i in resToWrite:
                np.savetxt(f,i ,delimiter='\t',fmt='%.8f')
                f.write('\n')




tauth = np.loadtxt('../tauth.dat')
tauph = np.loadtxt('../tauph.dat')

tauth.shape = (46,121,-1)
tauph.shape = (46,121,-1)


# # parse grid and NACTs
gridTh = np.deg2rad(tauth[:,0,0])
gridPh = np.deg2rad(tauth[0,:,1])
tauth = tauth[:,:,[2,3,4]]
tauph = tauph[:,:,[2,3,4]]


thGrid = np.deg2rad(np.linspace(0,90,91))
phGrid = np.deg2rad(np.linspace(0,180,181))

ss = adtSolver3S2D(gridTh, gridPh, tauth, tauph)

ss.solve()                                                     # solve on the original grid
ss.solve(thGrid=thGrid, phGrid=phGrid, resFile='angle_new.dat')  # solve on a different grid
ss.solve(initZero=True, resFile='angle_sym.dat')      # just 1D along theta, no phi integration, useful for making symmetric along phi
