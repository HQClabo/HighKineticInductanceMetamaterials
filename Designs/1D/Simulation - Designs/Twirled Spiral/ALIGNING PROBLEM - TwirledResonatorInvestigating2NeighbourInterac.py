# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 10:43:23 2021

@author: jouanny
"""

import gdspy
import numpy as np
from Modules.ModuleResonator import *


List = []

SqSpiral = {"layer":0, "datatype":3}


Centre = (0,0)
L = 0.75e-3#The length of the resonator
W = 500e-9#The width of the wire
I = 10.085e-6#The interspacing, this value here just so we don't have an incomplete branch at the end.

#Res = TwirledSpiral(Centre, L, W, I)

#The part below is just here to generate the spiral and add it to a newly created cell.
carac = {"layer": 0, "datatype": 3}#Characterize what layer, and datatype you want to attribute to a shape you created
Res = TwirledSpiral(Centre, L, W, I, [0,0])#Calling the function for the double spiral
a = gdspy.FlexPath(Res[0], W*1e6, **carac).rotate(np.pi/2, (Res[0][0]))#Creating the polygon for both arms of the spirals.
b = gdspy.FlexPath(Res[1], W*1e6, **carac).rotate(np.pi/2, (Res[0][0]))#I used rotate here as the arms of the spiral weren't in the right position
Resonator=gdspy.Cell("Resonator")#Creation of the cell to host the spiral.
Resonator.add(a)#adding the two arms
Resonator.add(b)

#Getting sizes+positions of a single resonator
El = np.array([Res[0][-1][1]-W*1e6/2, Res[0][0][1] + 1/2*I*1e6 + W*1e6/2])
Er = np.array([Res[1][-1][1] +W*1e6/2 ,Res[0][0][1] + 1/2*I*1e6 + W*1e6/2])


#Generating an array from this resonator
#In this simulation there is a fixed number of resonators, we just vary one resonator to the right to see how the
#resonant behaviour is influenced by the second nearest neighbour.


N = 3#Number of resonators

Spacing = 5e-6

Array = gdspy.Cell("Array")
array = gdspy.CellArray(Resonator, N, 1, ((Er[0]-El[0]) + Spacing*1e6 ,0))
Array.add(array)
Array.flatten()

#Right edges of the last resonator
ErLast = np.array([Er[0] + (N-1)*((Er[0]-El[0]) + Spacing*1e6),El[1]])


#Second array to modulate the spacing of the second resonator
N2 = 2
Spacing2 = 10e-6

Array2 = gdspy.Cell("Array2")
array2 = gdspy.CellArray(Resonator, N2, 1, ((Er[0]-El[0]) + Spacing*1e6 ,0), origin = (ErLast[0] + (Er[0]-El[0])/2 + Spacing2*1e6, ErLast[1]))
Array2.add(array2)
Array2.flatten()

ErTrueLast = np.array([ErLast[0] + (N2)*((Er[0]-El[0]))+ Spacing*1e6 + Spacing2*1e6,El[1]])


#Importing the CPW

Hr = (Res[0][-1][0]-Res[0][-2][0])*1e-6#Height of the capacitor
t = 0
Sr = 0

CPW = Cpw(El*1e-6,ErTrueLast*1e-6,Hr,Sr,t=t)

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
cpw.add(Array2)

lib = gdspy.GdsLibrary()

lib.add(cpw)

lib.write_gds("2ndNeighbour_twirl.gds")


