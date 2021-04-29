# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 10:40:34 2021

@author: jouanny
"""

import gdspy
from ModuleResonator import *
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

El1 = ab[0][-1][1]*1e-6-W/2, ab[0][0][1]*1e-6 + 1/2*I + W/2
Er1 = ab[1][-1][1]*1e-6 +W/2 ,ab[0][0][1]*1e-6 + 1/2*I + W/2 

Srr = 320e-6 #spacing between the two resonators

Centre2 = [-2*(ab[0][-1][1]*1e-6-W/2) + I + W + Srr,0] 
cd = TwirledSpiral(Centre2, L, W, I, [0,0])
c = gdspy.FlexPath(cd[0], W*1e6, **carac).rotate(np.pi/2, (cd[0][0]))
d = gdspy.FlexPath(cd[1], W*1e6, **carac).rotate(np.pi/2, (cd[0][0]))
testcell.add(c)
testcell.add(d)

El2 = cd[1][-1][1]*1e-6-W/2, cd[0][0][1]*1e-6 + 1/2*I + W/2
Er2 = cd[0][-1][0]*1e-6 +W/2,cd[0][0][1]*1e-6 + 1/2*I + W/2


Hr = (ab[0][-1][0]-ab[0][-2][0])*1e-6

#El = (0,0)
#Er = (100e-6, 0)
#Hr = 50e-6
Sr = 1e-6





A = Cpw(El1,Er2,Hr,Sr)

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
lib.write_gds("2Resonator_CPW_Srr320.gds")
