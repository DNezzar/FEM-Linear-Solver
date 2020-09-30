#!/usr/bin/env python
# -*- coding: utf8 -*-

import module
import numpy as np
import matplotlib.pyplot as plt

#print(module.mod.solver.__doc__)

#-----------------------INPUT-----------------------

node=np.array([[0,0],[0,3000],[3000,0],[3000,3000],[6000,0],[6000,3000]],dtype=np.float)
elem=np.array([[1,2],[1,3],[2,3],[2,4],[1,4],[3,4],[3,6],[4,5],[4,6],[3,5],[5,6]],dtype=np.int64)
free=np.array([3,4,5,6,7,8,9,11,12],dtype=np.int64)
A=np.float(300)
E=np.float(70000)
force=np.array([0,0,0,-50000,0,0,0,-100000,0,0,0,-50000],dtype=np.float)

#-----------------------SOLVER----------------------
[d,strain,stress]=module.mod.solver(elem,node,free,A,E,force)

#-------------------POST PROCESSING-----------------
scl=10
n1=elem[:,0]
n2=elem[:,1]
x1=node[n1-1,0]
y1=node[n1-1,1]
x2=node[n2-1,0]
y2=node[n2-1,1]
u1=d[2*n1-2]
v1=d[2*n1-1]
u2=d[2*n2-2]
v2=d[2*n2-1]

plt.plot([x1,x2],[y1,y2],'k--')
plt.plot([x1+scl*u1,x2+scl*u2],[y1+scl*v1,y2+scl*v2],'k')
plt.title('Deformed shape of a 2D Truss')
plt.xlabel('x-position')
plt.ylabel('y-position')

plt.show()
