import numpy as np

class Region:
    def __init__(self, mesh,umat,ndim,mnode,trcounter):
        self.mesh = mesh
        self.umat = umat
        self.ndim = ndim
        self.mnode = mnode
        self.trcounter = trcounter
