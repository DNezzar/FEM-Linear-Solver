#!/usr/bin/env python
# -*- coding: utf8 -*-

import module

#print(module.mod.linear.__doc__)

A=[[5.0,1.0,-4.0,0.0],[1.0,5.0,0.0,0.0],[-4.0,0.0,5.0,-1.0],[0.0,0.0,-1.0,5.0]]

b=[1.0,0.0,0.0,0.0]

x=module.mod.linear(A,b)

print(x)






