# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 15:19:21 2021

@author: vweibel
"""

import gdspy
import ResonatorDesign
import numpy as np

Centre = (0,0)
L = 0.75e-3
W = 500e-9
I = 10.085e-6

carac = {"layer": 0, "datatype": 3}
ab = TwirledSpiral(Centre, L, W, I, [0,0])
a = gdspy.FlexPath(ab[0], W*1e6, **carac).rotate(np.pi/2, (ab[0][0]))
b = gdspy.FlexPath(ab[1], W*1e6, **carac).rotate(np.pi/2, (ab[0][0]))
testcell=gdspy.Cell("testcell")
testcell.add(a)
testcell.add(b)


#Dimensions of one resonator (see in vjws journal 26/02/21 for graphical explanation of formula)
Lsize = ab[1][-1][1]-ab[0][-1][1] + W*1e6
#Height of the resonator
Hsize = ab[0][-1][0]-ab[1][-1][0]

Hr = Hsize*1e-6 #height of resonator in SI units
Lr = Lsize*1e-6 #length of resonator in SI units
center_first = ab[0][0]

Sr = 1e-6

#Generating the ground plane + cpw
A = Cpw_twirled(center_first,center_first,ab[1][-1][1],Hr,Lr,Sr,L, W, I)

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
cpw.add(testcell)

#We create the library where our cell will be.
lib = gdspy.GdsLibrary()
lib.add(cpw)
lib.write_gds("CPW_ghost_1Res_Sr_1um.gds")