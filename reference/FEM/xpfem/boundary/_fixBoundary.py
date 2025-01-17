import numpy as np
from ._Direct import Direct

class fixBoundary(Direct):
    def __init__(self,ndim,method = 'change_1'):
        self.ndim = ndim
        self.method = method
    
    def loadFix(self,filename = 'fixNode.dat'):
        self.fix = np.loadtxt(filename)

    def fixValue(self):
        dof = self.expandNdf(self.fix[:,0:2])
        value = self.fix[:,2]

        super().__init__(dof, value)

    def expandNdf(self,fixNode):
        fixNodeNew = ((fixNode[:, 0] - 1) * self.ndim + fixNode[:, 1]-1).astype('int')
        return fixNodeNew