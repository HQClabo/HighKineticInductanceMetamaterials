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


#Dimensions of one resonator (see in vjws journal 26/02/21 for graphical explanation of formula)
Lsize = ab[1][-1][1]-ab[0][-1][1] + W*1e6
#Height of the resonator
Hsize = ab[0][-1][0]-ab[1][-1][0]

Rep = 2#Number of repetitions in unit cell
Srr = 3e-6# spacing between the two resonators/ intracell spacing

#Generating the cell
Cell = gdspy.Cell("Cell")#Create the new cell for the array
CellC = gdspy.CellArray(testcell, Rep, 1, (Lsize+Srr*1e6,Hsize))
Cell.add(CellC)
U = Cell.flatten()#This flatten command here, is used to transform all cells into polygons. It reduces the risk of hierarchy problem.


#Generating the array
RepArr = 4 #minimum 2
SrrArr = 2e-6 #intercell spacing

Arr = gdspy.Cell("Arr")
ArrC = gdspy.CellArray(U, RepArr, 1, ((Lsize+Srr*1e6)*Rep+SrrArr*1e6-Srr*1e6,Hsize))
Arr.add(ArrC)
V = Arr.flatten()

#El = (0,0)
#Er = (100e-6, 0)
#Hr = 50e-6
Sr = 1e-6#Spacing between the resonator and the capacitor.

Hr = Hsize*1e-6 #height of resonator in SI units
Lr = Lsize*1e-6 #length of resonator in SI units
center_first = ab[0][0] #center of first resonator in um
center_last = center_first[0] + (Rep*RepArr-1)*Lsize + RepArr*(Rep-1)*Srr*1e6 + (RepArr-1)*SrrArr*1e6, center_first[1] #center last resonator in um
#Er2 = Er1[0] + ((2*Lsize+Srr)+SrrArr)*(RepArr)-Lsize-SrrArr, ab[0][0][1]*1e-6 + 1/2*I + W/2#Right edge of the last resonator ##(Rep-1)*(Lsize+Srr)

#Generating the ground plane + cpw
A = Cpw_twirled(center_first,center_last,ab[1][-1][1],Hr,Lr,Sr,L, W, I)
#Generating polygons from the previous function
carac2 = {"layer": 0, "datatype": 0}
LCpw = gdspy.Polygon(A[0])
HeadL1 = gdspy.FlexPath(A[1][0], W*1e6, **carac2).rotate(np.pi/2, (A[1][0][0]))
HeadL2 = gdspy.FlexPath(A[1][1], W*1e6, **carac2).rotate(np.pi/2, (A[1][0][0]))
#c = gdspy.FlexPath(A[1][0], W*1e6, **carac2)
#d = gdspy.FlexPath(A[1][1], W*1e6, **carac2)
RCpw = gdspy.Polygon(A[2])
HeadR1 = gdspy.FlexPath(A[3][0], W*1e6, **carac2).rotate(np.pi/2, (A[3][0][0]))
HeadR2 = gdspy.FlexPath(A[3][1], W*1e6, **carac2).rotate(np.pi/2, (A[3][0][0]))
TGround = gdspy.Polygon(A[4])
BGround = gdspy.Polygon(A[5])

#We add everything to a new cell
cpw = gdspy.Cell("cpw")
cpw.add(LCpw)
cpw.add(HeadL1)
cpw.add(HeadL2)
cpw.add(RCpw)
cpw.add(HeadR1)
cpw.add(HeadR2)
cpw.add(TGround)
cpw.add(BGround)
cpw.add(V)
#cpw.add([c,d])

#We create the library where our cell will be.
lib = gdspy.GdsLibrary()
lib.add(cpw)
lib.write_gds("CPW_ghost_8Res_3um_2um.gds")
