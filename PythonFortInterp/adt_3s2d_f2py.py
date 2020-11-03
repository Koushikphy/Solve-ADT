# 3 state 2D adt solve
import sys
import numpy as np 
from numpy import sin,cos,tan
from scipy.integrate import solve_ivp
# from scipy.interpolate import InterpolatedUnivariateSpline as spl, RectBivariateSpline as rbs
from fInterp import newinterpol as interp 




class adtSolver3S2D():

    def __init__(self, tauPhiFile='../tauph.dat', tauTheFile='../tauth.dat', shape=None):

        tauth = np.loadtxt(tauTheFile)
        tauph = np.loadtxt(tauPhiFile)

        tauth.shape = shape
        tauph.shape = shape
        
        self.theGrid = tauph[:,0,0]
        self.phiGrid = tauph[0,:,1]
        
        self.tauth = tauth[:,:,2:].transpose(2,1,0)
        self.tauph = tauph[:,:,2:].transpose(2,1,0)


    def resThe(self,theta,y,phi): # along theta
        tau = interp(self.tauth, self.theGrid, self.phiGrid, theta, phi) 
        return self.diffSys(y,tau)


    def resPhi(self,phi,y,theta): # along Phi
        tau = interp(self.tauph, self.theGrid, self.phiGrid, theta, phi) 
        return self.diffSys(y,tau)


    def diffSys(self,y,t):
        y12 = -t[0] - tan(y[1])*( t[1]*sin(y[0]) + t[2]*cos(y[0]))
        y13 = -t[1]*cos(y[0]) + t[2]*sin(y[0])
        y23 = -(1.0/cos(y[1]))*(t[1]*sin(y[0]) + t[2]*cos(y[0]))
        return y12,y13,y23


    def solve(self):

        # first solve along phi grid for the 1st theta value
        initVal = [0,0,0]
        rang    = [self.phiGrid[0],self.phiGrid[-1]]
        sol     = solve_ivp(self.resPhi, rang, initVal, method='LSODA', t_eval=self.phiGrid, args=(self.theGrid[0],))
        iThVal  = sol.y.T  #<<<--- initial values for integration along the theta grid

        res  = np.empty((46,121,3))
        rang = [self.theGrid[0],self.theGrid[-1]]
        for i,phi in enumerate(self.phiGrid):
            sol = solve_ivp(self.resThe, rang, iThVal[i], method='LSODA', t_eval=self.theGrid, args=(phi,))#, rtol=1.0e-8, atol=1.0e-10)
            res[:,i,:] = sol.y.T

        xx,yy = np.meshgrid(self.theGrid, self.phiGrid, indexing='ij')
        resToWrite = np.dstack([xx,yy,res])

        with open('angle_fInterp.dat','w') as f:
            for i in resToWrite:
                np.savetxt(f,i ,delimiter='\t',fmt='%.8f')
                f.write('\n')




ss = adtSolver3S2D(shape=(46,121,-1))
ss.solve()
