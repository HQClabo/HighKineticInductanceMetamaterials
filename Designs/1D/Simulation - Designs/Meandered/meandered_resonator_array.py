# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 15:05:19 2021

@author: vweibel
"""

import numpy as np
import gdspy
from Modules.ModuleResonator import * 

#parameters for drawing the meander resonator
height = 50e-6 #height of resonator
bigspacing = 24e-6 #spacing between "arms" and "body"
smallspacing = 4.5e-6 #spacing between "body" meanders
width=500e-9 #width of wire
intracell = 3e-6 #intracell spacing
intercell = 1e-6
N = 3 #number of unit cells
snakedim = [height,bigspacing,smallspacing,width] #combine them

#create unit cell with 2 meanders
carac = {"layer": 0, "datatype": 3}
snake1, center1 = meandered(height,bigspacing,smallspacing,width,start=(0,0))
snake2, center2 = meandered(height,bigspacing,smallspacing,width,start=(snake1[-1][0]+intracell*1e6+(width*1e6)/2,snake1[0][1]))
meander_m = gdspy.FlexPath(snake1,width*1e6,**carac)
meander_w = gdspy.FlexPath(snake2,width*1e6,**carac).rotate(np.pi,center2)
unitcell = gdspy.Cell("unit cell")
unitcell.add([meander_m,meander_w])
unitcellsize = snake2[-1][0]-snake1[0][0] + width*1e6
#unitcell = unitcell.flatten()

#determine spacing to pass to array-function = unitcell size + intercell spacing
arrayspacing_x = unitcellsize + intercell*1e6
arrayspacing_y = 0.0

#create array
resonatorarray = gdspy.Cell("resonator array")
manysnakes = gdspy.CellArray(unitcell,N,1,(arrayspacing_x,arrayspacing_y))
resonatorarray.add(manysnakes)
resonatorarray.flatten()

end = N*unitcellsize + (N-1)*intercell*1e6-width*1e6/2, snake1[-1][1]
Sr = 2e-6

A = Cpw_meandered(snake1, snakedim, end, Sr, w = 10e-6, Lc = 215e-6, Lline = 470e-6, T = 50e-6, s = 7.1e-6)
#Generating polygons from the previous function
LCpw = gdspy.Polygon(A[0])
HeadLeft = gdspy.FlexPath(A[1], width*1e6,**carac).rotate(np.pi,A[2])
RCpw = gdspy.Polygon(A[3])
HeadRight = gdspy.FlexPath(A[4], width*1e6,**carac)
TGround = gdspy.Polygon(A[5])
BGround = gdspy.Polygon(A[6])

#add everything to cell cpw
cpw = gdspy.Cell("cpw")
cpw.add([LCpw, HeadLeft,RCpw,HeadRight, TGround, BGround])
cpw.add(resonatorarray)
lib = gdspy.GdsLibrary()
lib.add(cpw)
lib.write_gds('meandered_test_array.gds')
