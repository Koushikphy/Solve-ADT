# Explicit solution for 3S and 4S ADT equation. Also a generalized diabatization function
# preferably use python 3.x
#-----------------------------------------------------------------------------------------

# import sys
import numpy as np 
from numpy import sin,cos,tan
from scipy.integrate import solve_ivp
from scipy.interpolate import RectBivariateSpline as rbs
from scipy.interpolate import InterpolatedUnivariateSpline as spl
from functools import reduce

sec = lambda x : 1.0/cos(x)



class ADT_Base():
    def __init__(self):
        pass 


    def calculate1Dcurve(self, phiGrid, tauph):
        self.phiGrid = phiGrid
        self.diffShape = tauph.shape[-1]
        self.tauPhiFuncs1D = [spl(phiGrid ,tauph[:,i]) for i in  range(self.diffShape)]


    def calculate2Dsurf(self, theGrid, phiGrid, tauth, tauph):
        self.theGrid = theGrid
        self.phiGrid = phiGrid
        self.diffShape = tauph.shape[-1]
        self.tauTheFuncs = [rbs(theGrid , phiGrid,tauth[...,i]) for i in  range(self.diffShape)]
        self.tauPhiFuncs = [rbs(theGrid , phiGrid,tauph[...,i]) for i in  range(self.diffShape)]


    def diffSys(self,y,t):
        # 3 state differential equation
        # A = A12*A13*A23 #
        y12 = -t[0] - tan(y[1])*( t[1]*sin(y[0]) + t[2]*cos(y[0]))
        y13 = -t[1]*cos(y[0]) + t[2]*sin(y[0])
        y23 = -(1.0/cos(y[1]))*(t[1]*sin(y[0]) + t[2]*cos(y[0]))
        return y12,y13,y23


    def diffSys4S(self,y,t):
        # 4 state ADT equations
        y12= -t[0]-sin(y[0])*tan(y[1])*t[1]-cos(y[0])*tan(y[1])*t[2]-sin(y[0])*sec(y[1])*tan(y[3])*t[3]-cos(y[0])*sec(y[1])*tan(y[3])*t[4]

        y13 = -cos(y[0])*t[1]+sin(y[0])*t[2]-cos(y[0])*sin(y[1])*tan(y[3])*t[3]+sin(y[0])*sin(y[1])*tan(y[3])*t[4]-cos(y[1])*tan(y[3])*t[5]

        y23 = -cos(y[1])*( t[1]*sin(y[0])*(sec(y[1]))**2+cos(y[2])*sec(y[3])*(t[5]-t[4]*sin(y[0])*tan(y[1]))*tan(y[4])\
        +sin(y[0])*sec(y[1])*tan(y[1])*t[3]*tan(y[3]) + sin(y[2])*t[3]*sec(y[3])*tan(y[4])\
        +cos(y[0])*(t[2]*(sec(y[1]))**2+tan(y[1])*cos(y[2])*t[3]*sec(y[3])* tan(y[4])+  sec(y[1])*(tan(y[1])*t[4]*tan(y[3]) +sin(y[2])*t[4]*sec(y[3])*tan(y[4]))))

        y14 = -cos(y[0])*cos(y[1])*t[3]+sin(y[0])*cos(y[1])*t[4] +sin(y[1])*t[5]

        y24 = sin(y[2])*(-sin(y[0])*sin(y[1])*t[4]*sec(y[3])+cos(y[1])*t[5]*sec(y[3]))- sin(y[0])*cos(y[2])*t[3]*sec(y[3])\
        +cos(y[0])*(sin(y[1])*sin(y[2])*t[3]*sec(y[3])-cos(y[2])*t[4]*sec(y[3])) 

        y34 = sin(y[0])*(-sin(y[2])*t[3]*sec(y[3])*sec(y[4])+sin(y[1])*cos(y[2])*t[4]*sec(y[3])*sec(y[4]))\
        +cos(y[0])*(-sec(y[4])*sin(y[1])*cos(y[2])*t[3]*sec(y[3])+sin(y[2])*t[4]*sec(y[3]))-cos(y[1])*cos(y[2])*t[5]*sec(y[3])*sec(y[4])

        return y12, y13, y23, y14, y24, y34


    def resPhi1D(self,phi,y):
        tau = [ f(phi) for f in self.tauPhiFuncs1D]
        tau = np.clip(tau, 0, np.inf)   # remove any negetive values during interplation
        return self.diffSys(y,tau)


    def resThe(self,theta,y,phi): # along theta
        tau = [f(theta,phi)[0,0] for f in self.tauTheFuncs]
        tau = np.clip(tau, 0, np.inf)   # remove any negetive values during interplation
        return self.diffSys(y,tau)


    def resPhi(self,phi,y,theta): # along Phi
        tau = [f(theta,phi)[0,0] for f in self.tauPhiFuncs]
        tau = np.clip(tau, 0, np.inf)   # remove any negetive values during interplation
        return self.diffSys(y,tau)


    def solv1D(self, phGrid=None):
        # Solve for 1D
        if phGrid is None: phGrid = self.phiGrid

        sol = solve_ivp(self.resPhi1D, [phGrid[0],phGrid[-1]], np.zeros(self.diffShape), method='LSODA', t_eval=phGrid, dense_output=True)
        return np.column_stack([sol.t, sol.y.T])


    def solveP5(self, thGrid=None, phGrid=None, initZero=False):
        # Path 5:  first along theta grid, then along phi grid
        if thGrid is None: thGrid = self.theGrid
        if phGrid is None: phGrid = self.phiGrid

        if initZero:
            iPhVal = np.full((thGrid.shape[0],self.diffShape),0)
        else:
            sol    = solve_ivp(self.resThe, [thGrid[0],thGrid[-1]], np.zeros(self.diffShape), method='LSODA', t_eval=thGrid, args=(phGrid[0],))
            iPhVal = sol.y.T                                      #<<<--- initial values for integration along the theta grid

        res  = np.empty((thGrid.shape[0], phGrid.shape[0], self.diffShape))
        for i,th in enumerate(thGrid):
            sol = solve_ivp(self.resPhi, [phGrid[0],phGrid[-1]], iPhVal[i], method='LSODA', t_eval=phGrid, args=(th,))
            res[i,...] = sol.y.T

        xx,yy = np.meshgrid(thGrid, phGrid, indexing='ij')
        resToWrite = np.dstack([xx,yy,res])

        resToWrite.shape = (-1,self.diffShape+2)
        return resToWrite


    def solveP1(self, thGrid=None, phGrid=None, initZero=False):
        # Path 1: first along Phi grid, then along Theta grid
        if thGrid is None: thGrid = self.theGrid
        if phGrid is None: phGrid = self.phiGrid

        if initZero:
            iThVal  = np.full((phGrid.shape[0],self.diffShape),0)
        else: # first solve along phi grid for the 1st theta value
            sol     = solve_ivp(self.resPhi, [phGrid[0],phGrid[-1]], np.zeros(self.diffShape), method='LSODA', t_eval=phGrid, args=(thGrid[0],)) 
            iThVal  = sol.y.T  #<<<--- initial values for integration along the theta grid

        res  = np.empty((thGrid.shape[0], phGrid.shape[0], self.diffShape))
        for i,phi in enumerate(phGrid):
            sol = solve_ivp(self.resThe, [thGrid[0],thGrid[-1]], iThVal[i], method='LSODA', t_eval=thGrid, args=(phi,))
            res[:,i,:] = sol.y.T

        xx,yy = np.meshgrid(thGrid, phGrid, indexing='ij')
        resToWrite = np.dstack([xx,yy,res])

        resToWrite.shape = (-1,self.diffShape+2)
        return resToWrite



def ang2diab(enr, ang):
    # A generalized diabatization function from energy and angle
    ns = enr.shape[1]
    ts = ang.shape[1]
    #print(ns,ts)

    assert ts==(ns*(ns-1)/2), "Mismatch in energy and nact"

    def createList():
        # create unique list of diabatic element
        mt = np.arange(ns*ns)
        mt.shape = (ns,ns)
        ls = [mt[i,i] for i in range(ns)] 

        for j in range(ns):
            for i in range(j):
                ls.append(mt[i,j])

        print(ls)
        return ls

    ls = createList()

    def aMat(y): 
        arr = np.zeros((ts, ns, ns))

        k = 0 
        for j in range(ns):
            for i in range(j):
                arr[k,...] = np.eye(ns)
                arr[k,i,j] = sin(y[k])
                arr[k,j,i] = -sin(y[k])
                arr[k,i,i] = arr[k,j,j] = cos(y[k])
                k+=1

        return reduce(np.matmul,arr)

    amat  = np.array([aMat(a) for a in ang])

    #dbDat = np.einsum("ijk,ij,ijl->ikl",amat,enr,amat).reshape(-1,16)[:,[0,4,7,9, 1, 2,5, 3,6,8]]
    dbDat = np.einsum("ijk,ij,ijl->ikl",amat,enr,amat).reshape(-1,ns*ns)[:,ls]
    return dbDat





tauth = np.loadtxt('TauThetaNACT_new.dat', usecols=(1,2,3,4,5))
# prepare the data for ADt
tauth[:,:2] = np.deg2rad(tauth[:,:2])
tauth.shape = (66,361,5)
gridTh = tauth[:,0,0]
gridPh = tauth[0,:,1]
tauth = tauth[...,2:]
# instantiate the ADT object
ss = ADT_Base()
ss.calculate2Dsurf(gridTh, gridPh, tauth, tauth)
# ss = adtSolver3S2DTheta(gridTh, gridPh, tauth, tauth)  # sending dummy tauph, as its not required



nOfGrid = 500 # new grid for ADt
newTh = np.linspace(np.min(gridTh), np.max(gridTh), nOfGrid)
newPh = np.linspace(np.min(gridPh), np.max(gridPh), nOfGrid)

angleDense = ss.solveP1(thGrid=newTh, phGrid=newPh, initZero=True) # do ADT on a new grid


# angle = rectGridInt(angleDense, 0, 1, 91, 361) # interpolate down the angle and save in a file
# angle[:,:2] = np.rad2deg(angle[:,:2])
# writeFile('angle.dat', angle)


# enr = np.loadtxt('./send/Enr_Mrci_new.dat', usecols=(1,2,3,4,5))   # read energy file 
# enr = rectGridInt(enr, 0, 1, 91, 361)                             # interp the enrgy data to the new grid
# grid = angle[:,:2]
# diab = ang2diab(enr[:,2:] ,angle[:,2:]) # do dibatize, returned as 11,22,33,12,13,23

# diab = np.column_stack([grid, diab])  # save
# writeFile('diab.dat', diab)