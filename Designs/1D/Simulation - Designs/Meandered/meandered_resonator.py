# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 15:05:19 2021

@author: vweibel
"""

import numpy as np
import gdspy
from ModuleResonator import *        
    
height = 50e-6
bigspacing = 24e-6
smallspacing = 4.5e-6
width=500e-9
intracell = 3e-6
snakedim = [height,bigspacing,smallspacing,width]
carac = {"layer": 0, "datatype": 3}
snake1, center1 = meandered(height,bigspacing,smallspacing,width,start=(0,0))
snake2, center2 = meandered(height,bigspacing,smallspacing,width,start=(snake1[-1][0]+intracell*1e6,snake1[0][1]))
#print(snake1[-1][0],center2,center1)
meander_m = gdspy.FlexPath(snake1,width*1e6,**carac)
meander_w = gdspy.FlexPath(snake2,width*1e6,**carac).rotate(np.pi,center2)

end = snake2[-1]
Sr = 2e-6

A = Cpw_meandered(snake1, snakedim, end, Sr, w = 10e-6, Lc = 215e-6, Lline = 470e-6, T = 50e-6, s = 7.1e-6)
#calculating center of left head for rotation
#Generating polygons from the previous function
LCpw = gdspy.Polygon(A[0])
HeadLeft = gdspy.FlexPath(A[1], width*1e6,**carac).rotate(np.pi,A[2])
RCpw = gdspy.Polygon(A[3])
HeadRight = gdspy.FlexPath(A[4], width*1e6,**carac)
TGround = gdspy.Polygon(A[5])
BGround = gdspy.Polygon(A[6])

#add everything to cell cpw
cpw = gdspy.Cell("cpw")
cpw.add(LCpw)
cpw.add(HeadLeft)
cpw.add(RCpw)
cpw.add(HeadRight)
cpw.add(TGround)
cpw.add(BGround)
cpw.add([meander_m,meander_w])
lib = gdspy.GdsLibrary()
lib.add(cpw)
lib.write_gds('meandered_test_cpw.gds')
