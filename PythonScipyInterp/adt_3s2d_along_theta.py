# 3 state 2D adt solve, preferably use python 3.x
# import sys
import numpy as np 
from numpy import sin,cos,tan
from scipy.integrate import solve_ivp
from scipy.interpolate import RectBivariateSpline as rbs



def writeFile(file,dat,tc=0,fmt='%.8f'):
    with open(file,'w') as f:
        for i in np.unique(dat[:,tc]):
            np.savetxt(f,dat[dat[:,tc]==i],delimiter='\t', fmt=fmt)
            f.write('\n')
    print("File saved as: "+file)



def rectGridInt(data, grid1Col, grid2Col, newGrid1, newGrid2):
    g1 = np.unique(data[:,grid1Col])
    g2 = np.unique(data[:,grid2Col])
    c1 = g1.shape[0]
    c2 = g2.shape[0]

    ng1 = np.linspace(g1.min(),g1.max(),newGrid1)
    ng2 = np.linspace(g2.min(),g2.max(),newGrid2)
    ng1m, ng2m = np.meshgrid(ng1, ng2, indexing='ij')

    result = []
    for i in range(data.shape[1]):
        if i==grid1Col:
            res = ng1m.reshape(-1)
        elif i==grid2Col:
            res = ng2m.reshape(-1)
        else:
            res = rbs(g1,g2,data[:,i].reshape(c1,c2))(ng1, ng2).reshape(-1)
        result.append(res)
    return np.column_stack(result)




class adtSolver3S2DTheta:
    def __init__(self, theGrid, phiGrid, tauth, tauph):
        self.theGrid = theGrid
        self.phiGrid = phiGrid
        self.tauTheFuncs = [rbs(theGrid , phiGrid,tauth[...,i]) for i in [0,1,2]]
        self.tauPhiFuncs = [rbs(theGrid , phiGrid,tauph[...,i]) for i in [0,1,2]]


    def resThe(self,theta,y,phi): # along theta
        tau = [f(theta,phi)[0,0] for f in self.tauTheFuncs]
        tau = np.clip(tau, 0, np.inf)   # remove any negetive values during interplation
        return self.diffSys(y,tau)


    def resPhi(self,phi,y,theta): # along Phi
        tau = [f(theta,phi)[0,0] for f in self.tauPhiFuncs]
        tau = np.clip(tau, 0, np.inf)   # remove any negetive values during interplation
        return self.diffSys(y,tau)


    def diffSys(self,y,t):
        # A = A12*A13*A23 #
        y12 = -t[0] - tan(y[1])*( t[1]*sin(y[0]) + t[2]*cos(y[0]))
        y13 = -t[1]*cos(y[0]) + t[2]*sin(y[0])
        y23 = -(1.0/cos(y[1]))*(t[1]*sin(y[0]) + t[2]*cos(y[0]))
        return y12,y13,y23


    def solve(self, thGrid=None, phGrid=None, initZero=False):
        if thGrid is None: thGrid = self.theGrid
        if phGrid is None: phGrid = self.phiGrid

        if initZero:
            iThVal  = np.full((phGrid.shape[0],3),0)
        else: # first solve along phi grid for the 1st theta value
            sol     = solve_ivp(self.resPhi, [phGrid[0],phGrid[-1]], [0,0,0], method='LSODA', t_eval=phGrid, args=(thGrid[0],)) 
            iThVal  = sol.y.T  #<<<--- initial values for integration along the theta grid

        res  = np.empty((thGrid.shape[0], phGrid.shape[0], 3))
        for i,phi in enumerate(phGrid):
            sol = solve_ivp(self.resThe, [thGrid[0],thGrid[-1]], iThVal[i], method='LSODA', t_eval=thGrid, args=(phi,))
            res[:,i,:] = sol.y.T

        xx,yy = np.meshgrid(thGrid, phGrid, indexing='ij')
        resToWrite = np.dstack([xx,yy,res])

        resToWrite.shape = (-1,5)
        return resToWrite




def ang2diab(enr, ang):
    def aMat(a,b,c): # 12, 13, 23
        return np.array([
            [cos(a)*cos(b), sin(a)*cos(c)-cos(a)*sin(b)*sin(c), cos(a)*sin(b)*cos(c)+sin(a)*sin(c) ],
            [-sin(a)*cos(b),cos(a)*cos(c)+sin(a)*sin(b)*sin(c), -sin(a)*sin(b)*cos(c)+cos(a)*sin(c)],
            [-sin(b),       -cos(b)*sin(c),                     cos(b)*cos(c)                      ]
        ])
    amat  = np.array([aMat(*a) for a in ang])
    dbDat = np.einsum("ijk,ij,ijl->ikl",amat,enr,amat).reshape(-1,9)[:,[0,4,8,1,2,5]]
    return dbDat






tauth = np.loadtxt('./send/TauThetaNACT_new.dat', usecols=(1,2,3,4,5))
# prepare the data for ADt
tauth[:,:2] = np.deg2rad(tauth[:,:2])
tauth.shape = (66,361,5)
gridTh = tauth[:,0,0]
gridPh = tauth[0,:,1]
tauth = tauth[...,2:]
# instantiate the ADT object
ss = adtSolver3S2DTheta(gridTh, gridPh, tauth, tauth)  # sending dummy tauph, as its not required



nOfGrid = 500 # new grid for ADt
newTh = np.linspace(np.min(gridTh), np.max(gridTh), nOfGrid)
newPh = np.linspace(np.min(gridPh), np.max(gridPh), nOfGrid)

angleDense = ss.solve(thGrid=newTh, phGrid=newPh, initZero=True) # do ADT on a new grid


angle = rectGridInt(angleDense, 0, 1, 91, 361) # interpolate down the angle and save in a file
angle[:,:2] = np.rad2deg(angle[:,:2])
writeFile('angle.dat', angle)


enr = np.loadtxt('./send/Enr_Mrci_new.dat', usecols=(1,2,3,4,5))   # read energy file 
enr = rectGridInt(enr, 0, 1, 91, 361)                             # interp the enrgy data to the new grid
grid = angle[:,:2]
diab = ang2diab(enr[:,2:] ,angle[:,2:]) # do dibatize, returned as 11,22,33,12,13,23

diab = np.column_stack([grid, diab])  # save
writeFile('diab.dat', diab)