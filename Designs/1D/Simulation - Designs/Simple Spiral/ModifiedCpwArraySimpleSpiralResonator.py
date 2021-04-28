# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 10:43:23 2021

@author: jouanny
"""

import gdspy
import numpy as np
from ModuleResonator import spiral, Cpw

#Importing the base spiral and adding it to a cell.
L = 0.75e-3
w = 500e-9
I = 0.5e-6
D = 70e-6
Centre = (0,-I*1e6/2)

Res = spiral(L,w,I,D,Centre)[0]

List = []

SqSpiral = [{"layer":0, "datatype":3},{"layer":0, "datatype":0}]


for k in range(len(Res)-1):
    if k == 0 or k == range(len(Res)-1)[-1]:
        List.append(gdspy.Rectangle(Res[k], Res[k+1], **SqSpiral[1]))
    else:
        List.append(gdspy.Rectangle(Res[k], Res[k+1], **SqSpiral[0]))
Resonator = gdspy.Cell("Resonator")
Resonator.add(List)
Resonator.flatten()

#Getting sizes+positions of a single resonator
El = np.array([Res[0][0],(np.abs(Res[0][1])-np.abs(Res[1][1]))/2])#Left edge
Er = np.array([Res[2][0],El[1]])#Right edge


#Generating an array from this resonator

N = 3#Number of resonators

Spacing = 5e-6

Array = gdspy.Cell("Array")
array = gdspy.CellArray(Resonator, N, 1, ((Er[0]-El[0]) + Spacing*1e6 ,0))
Array.add(array)
Array.flatten()

#Right edges of the last resonator
ErLast = np.array([Er[0] + (N-1)*((Er[0]-El[0]) + Spacing*1e6),El[1]])

#Importing the CPW

Hr = D
t = 0
Sr = 0

CPW = Cpw(El*1e-6,ErLast*1e-6,Hr,Sr,t=t)

LCpw = gdspy.Polygon(CPW[0])
RCpw = gdspy.Polygon(CPW[1])
TGround = gdspy.Polygon(CPW[2])
BGround = gdspy.Polygon(CPW[3])

cpw = gdspy.Cell("cpw")
cpw.add(LCpw)
cpw.add(RCpw)
cpw.add(TGround)
cpw.add(BGround)
cpw.add(Array)

lib = gdspy.GdsLibrary()

lib.add(cpw)

lib.write_gds("CPW2SingleSimpleresonatorVarySpace-I0_5um.gds")


