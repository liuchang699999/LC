import numpy as np


from ._mesh import Mesh


def readtxt(
        nodeFileName='NLIST.dat', elemFileName='ELIST.dat',dim=None,delimiter=None
):
    from numpy import loadtxt
    node=loadtxt(nodeFileName,delimiter=delimiter)
    elem=loadtxt(elemFileName,delimiter=delimiter)

    node = node[:,1:]

    if dim is None:
        dim = node.shape[1]

    node = node[:,:dim]
    elem = elem[:,2:]

    elem = elem - 1

    mesh = Mesh(node,elem)
    return mesh