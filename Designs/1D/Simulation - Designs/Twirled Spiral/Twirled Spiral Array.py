# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 11:16:18 2021

@author: jouanny
"""

#This script is here to test the flattening function with the array function, it might solve the hierarchy problem that there was.

import numpy as np
import gdspy
from Modules.ModuleResonator import *



#Here these are just the parameters I use for the function TwirledSpiral which generates a double spiral.
Centre = (0,0)
L = 0.75e-3#The length of the resonator
W = 500e-9#The width of the wire
I = 10.085e-6#The interspacing, this value here just so we don't have an incomplete branch at the end.

#The part below is just here to generate the spiral and add it to a newly created cell.
carac = {"layer": 0, "datatype": 3}#Characterize what layer, and datatype you want to attribute to a shape you created
ab = TwirledSpiral(Centre, L, W, I, [0,0])#Calling the function for the double spiral
a = gdspy.FlexPath(ab[0], W*1e6, **carac).rotate(np.pi/2, (ab[0][0]))#Creating the polygon for both arms of the spirals.
b = gdspy.FlexPath(ab[1], W*1e6, **carac).rotate(np.pi/2, (ab[0][0]))#I used rotate here as the arms of the spiral weren't in the right position
testcell=gdspy.Cell("testcell")#Creation of the cell to host the spiral.
testcell.add(a)#adding the two arms
testcell.add(b)

#This part here gets the informations from the spiral,
#it's not the most elegant way to do it as if you rotate the spiral once more it messes everything.
#It could be nice to find a way to correct that.
#It works by taking the last point of the polygons and adding the thickness or the interspacing to them
El1 = ab[0][-1][1]*1e-6-W/2, ab[0][0][1]*1e-6 + 1/2*I + W/2
Er1 = ab[1][-1][1]*1e-6 +W/2 ,ab[0][0][1]*1e-6 + 1/2*I + W/2 


#Here, to obtain the size of the resonator I just substract the two sides of the resonator.
Lsize = Er1[0] - El1[0]
#Height of the resonator
Hsize = (ab[0][-1][0]-ab[0][-2][0])*1e-6

Rep = 2#Number of repetitions in unit cell
Srr = 10e-6# spacing between the two resonators/ intracell spacing

#Generating the cell
Cell = gdspy.Cell("Cell")#Create the new cell for the array
CellC = gdspy.CellArray(testcell, Rep, 1, (Lsize*1e6+Srr*1e6,Hsize))
Cell.add(CellC)
U = Cell.flatten()#This flatten command here, is used to transform all cells into polygons. It reduces the risk of hierarchy problem.


#Generating the array
RepArr = 8 #minimum 2
SrrArr = 2*1e-6 #intercell spacing

Arr = gdspy.Cell("Arr")
ArrC = gdspy.CellArray(U, RepArr, 1, ((Lsize*1e6+Srr*1e6)*Rep+SrrArr*1e6-Srr*1e6,Hsize))
Arr.add(ArrC)
V = Arr.flatten()

#El = (0,0)
#Er = (100e-6, 0)
#Hr = 50e-6
Sr = 10e-6#Spacing between the resonator and the capacitor.

Hr = (ab[0][-1][0]-ab[0][-2][0])*1e-6#Height of the capacitor
Er2 = Er1[0] + ((2*Lsize+Srr)+SrrArr)*(RepArr)-Lsize-SrrArr, ab[0][0][1]*1e-6 + 1/2*I + W/2#Right edge of the last resonator ##(Rep-1)*(Lsize+Srr)

#Generating the ground plane + cpw
A = Cpw(El1,Er2,Hr,Sr)
#Generating polygons from the previous function
LCpw = gdspy.Polygon(A[0])
RCpw = gdspy.Polygon(A[1])
TGround = gdspy.Polygon(A[2])
BGround = gdspy.Polygon(A[3])

#We add everything to a new cell
cpw = gdspy.Cell("cpw")
cpw.add(LCpw)
cpw.add(RCpw)
cpw.add(TGround)
cpw.add(BGround)
cpw.add(V)

#We create the library where our cell will be.
lib = gdspy.GdsLibrary()
lib.add(cpw)
lib.write_gds("CPW-16Res-intra10um-inter2um.gds")
