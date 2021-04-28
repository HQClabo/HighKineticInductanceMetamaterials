# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 10:40:34 2021

@author: jouanny
"""

import gdspy
from Modules.ModuleResonator import *
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

#Array = gdspy.CellArray(testcell, 3,1,100)



El = ab[0][-1][1]*1e-6-W/2, ab[0][0][1]*1e-6 + 1/2*I + W/2
Er = ab[1][-1][1]*1e-6 +W/2 ,ab[0][0][1]*1e-6 + 1/2*I + W/2 
Hr = (ab[0][-1][0]-ab[0][-2][0])*1e-6 + W

#El = (0,0)
#Er = (100e-6, 0)
#Hr = 50e-6
Sr = 9e-6

A = Cpw_finger(El,Er,Hr,Sr)

LCpw = gdspy.Polygon(A[0])
RCpw = gdspy.Polygon(A[1])
TGround = gdspy.Polygon(A[2])
BGround = gdspy.Polygon(A[3])

cpw = gdspy.Cell("cpw")
cpw.add(LCpw)
cpw.add(RCpw)
cpw.add(TGround)
cpw.add(BGround)
cpw.add(testcell)
#cpw.flatten()

lib = gdspy.GdsLibrary()
lib.add(cpw)
lib.write_gds("CPWFinger1Res_Sr9um.gds")