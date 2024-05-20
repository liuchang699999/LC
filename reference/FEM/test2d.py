import numpy as np
import xpfem as fem

# read the mesh
dir = 'input\\2D'
mesh=fem.mesh.readtxt(dir+'\\'+'NLIST.dat',dir+'\\'+'ELIST.dat',dim=2)
# Define Sets
set = {'set1':np.arange(0,mesh.sumElem)}

# Define Material
mat1 = fem.LinearElasticPlaneStress(E=2000,nu=0.3,name='Ico')
part = {'set1':mat1} # part 与 set对应
umat = fem.CreatMat(part,set)

ndim, mnode = 2, 4
trcounter = 0

# 有限元区域
region = fem.Region(mesh,umat,ndim,mnode,trcounter)
# 压力边界条件
Numman = fem.transPress(mnode, ndim)
Numman.loadPress(filename=dir+'\\'+'press.dat')
Numman.transpress(mesh.points,mesh.cells)

# 位移边界条件
fix = fem.fixBoundary(ndim,method = 'change_1')
fix.loadFix(filename=dir+'\\'+'fixNode.dat')
fix.fixValue()

Fe_analysis = fem.FEM_Solver_linearElasticity()

Fe_analysis.FEM_steady(region,Numman,fix) # 也可以一步步求解，不用FEM_steady方法

post = fem.Post(ndim,mnode)
post.showDisplacement(region,Fe_analysis.u,direction = 'x',show_deformed = True, scale = 10.0, show_shade = False, opacity = 0.8)
# ux = Fe_analysis.u[0::2]
# post.showScalar(region,ux,name='ux')