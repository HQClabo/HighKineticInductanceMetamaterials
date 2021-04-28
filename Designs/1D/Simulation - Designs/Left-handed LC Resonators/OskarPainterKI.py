# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 17:17:32 2021

@author: jouanny
"""

import numpy as np
import gdspy
from Modules.ModuleResonator import *

"""
~~~~~~ parameters to choose ~~~~~~
"""

L = 572.8e-6*2/3             #total length of inductor in SI units
s = 9e-6              #interspacing of inductor in SI units
w = 1e-6              #width of inductor wire
A = 100e-6               #horizontal dimension of capacitor
t = 4e-6                #thickness of capacitor plates
tv = 36e-6              #intracell spacing
tw = 36e-6              #intercell spacing
k = 1/6                 #fraction determining how thick the ground strip between resonators is : k*min(tv,tw)
e = 41e-6             #horizontal dimension of ground patches
f = 48e-6               #vertical dimension of ground patches
r = 58e-6               #vertical spacing of resonator "head" to ground plane
Q = 4                   #number of unit cells
N_ghost = 2             #number of ghosts: For no ghosts, put zero
Sf2r = [0,4.5e-6]       #spacing to feedline in [x,y]-direction
ground_yn = True        #ground in between yes (True) or no (False) 

carac_res = {'layer' : 0, 'datatype' : 3}       #layer + datatype of resonator array
carac_ghost = {'layer' :  1, 'datatype' : 1}    #layer + datatype of ghosts
carac_ground = {'layer' : 0, 'datatype' : 0}    #layer + datatype of ground planes and feedline

"""
~~~~ program structure ~~~~
1) draw unit cell and get relevant coordinates/dimensions (startM/U = 1st coordinate of first meander/condensator of array, stopU = 1st coordinate of last condensator of array,
                                                          center1st = center of first resonator ("down"), center2nd = center of second resonator ("up"))
2) centerQth = center of 2nd resonator in Qth (=last) unitcell
3) draw resonator array (ResArray, only resonators without grounding patches)
4) draw ground array (= grounding patches of meanders & - if desired - ground strips in between resonators) up to the (Q-1)th unitcell
5) add ground to Qth unitcell (why do it like this? - It's to be able to choose spacing between last resonator and ghost)                                                         
6) groundplane_coords contains the coordinates to draw the ground planes and the feedline.
7) get coordinates to draw the ghosts and the grounding of the ghosts. Can tune the spacing of the 1st/last resonator to the ghost by passing also tg = value in SI units)
8) everything stored in the cell variable design and exported as .gds file
"""

unitcell, ground, startM,startU, stopUy, strip_height,center1st, center2nd = unit_cell(L, s, w, A, t, tv,tw,e,f,r,carac = carac_ghost,gamma = k, ground_in_between = ground_yn)
unitcell_size = [2*A*1e6 + tv*1e6 + tw*1e6,0]
stopU = [startU[0]-A*1e6 + Q*unitcell_size[0]-tw*1e6, stopUy]

centerQth = [center1st[0]-A*1e6 + Q*unitcell_size[0] - tw*1e6, center2nd[1]]

design = gdspy.Cell('Resonator Array')

ResArray = gdspy.CellArray(unitcell, Q, 1, unitcell_size)
design.add(ResArray)

GrdArray = gdspy.CellArray(ground, Q-1, 1, unitcell_size)
#add ground to Qth unitcell
#corners of 1st ground patch: down
box_x1, box_y1 = centerQth[0] - tv*1e6 - A*1e6 - e*1e6/2 - w*1e6/2, startM[1]
box_x2, box_y2 = box_x1 + e*1e6, startM[1]-f*1e6
ground_box1 = gdspy.Rectangle([box_x1,box_y1], [box_x2,box_y2])
#corners of 2nd ground patch: up
box_x3, box_y3 = centerQth[0] - e*1e6/2 - w*1e6/2, center1st[1]+r*1e6
box_x4, box_y4 = box_x3 + e*1e6, box_y3 - f*1e6
ground_box2 = gdspy.Rectangle([box_x3,box_y3], [box_x4,box_y4])
groundingQth = gdspy.boolean(ground_box1, ground_box2, 'or')
if ground_yn == True:
    #2nd last ground strip
    strip_x1, strip_y1 = stopU[0] + t*1e6/2 - A*1e6 - 0.5*(tv*1e6 - k*min(tv,tw)*1e6), box_y2
    strip_x2, strip_y2 = strip_x1 - k*min(tv,tw)*1e6, strip_y1 + strip_height
    ground_strip_2nd_last = gdspy.Rectangle([strip_x1,strip_y1], [strip_x2,strip_y2])
    groundingQth = gdspy.boolean(groundingQth, ground_strip_2nd_last, 'or', **carac_ground)

groundplane_coords = waveguide_simulation(Q, unitcell_size, startM, startU, stopU, strip_height,tw, tv, N_ghost, t=t,Sr=Sf2r, Sg=120e-6, T = 60e-6, R=60e-6,f=f,A=A)
#Q, unitcell_size,startM, startU, stopU, strip_height, tw, tv, N_ghost, t = 2e-6, Sr = [1e-6,2e-6], Sg = 60e-6,T=100e-6,D=200e-6,f=24e-6, A = 50e-6
GrdArray = gdspy.boolean(GrdArray, groundingQth, 'or')
ground_plane = gdspy.boolean(groundplane_coords,GrdArray,'or', **carac_ground)

ghosts_coords, ground_ghosts = ghosts(L, s, w, A, t, tv,tw, N_ghost,strip_height,tg = tw,e=e,f=f,r=r,gamma=k,carac = carac_ghost, center_first=center1st,center_last=centerQth,ground_in_between=(ground_yn))
#L,s,w,A,t,tv,tw,N, strip_height, tg = 15e-6, e=20.5e-6, f=24e-6, r=29e-6, gamma=1/4, ground_in_between=True, carac = {'layer' :  1, 'datatype' : 1}, center_first = (0,0), center_last = (0,0)):
ground_plane = gdspy.boolean(ground_plane, ground_ghosts, 'or', **carac_ground)
if N_ghost > 0:
    design.add(ghosts_coords)
    
design.add(ground_plane)

                             
lib = gdspy.GdsLibrary()
lib.add(design)
# lib.write_gds("8-big-LCRes_round.gds")
gdspy.LayoutViewer()

"""
~~~~~ Vincent's original version below ~~~~~~
"""

# def OPKI(L,H,W,l,dh,dl,w, centre=(0,0)):
#     L = L*1e6
#     H = H*1e6
#     W = W*1e6
#     l = l*1e6
#     dh = dh*1e6
#     dl = dl*1e6
#     w = w*1e6
#     centre = centre[0]*1e6, centre[1]*1e6
#     #Ushape, we take the centre at at the bottom/top of the U
#     U=[]
#     x0, y0 = 0, 0
#     U.append([x0, y0])
    
#     x1, y1 = x0, y0 + H
#     U.append([x1, y1])
    
#     x2, y2 = x1 + L, y1
#     U.append([x2, y2])
    
#     x3, y3 = x2, y2 - H
#     U.append([x3, y3])
    
#     U = np.asarray(U)
#     centring = L/2, H
#     U = U - centring
    
#     #ZigZag shape
#     al = 0
#     Z = []
    
#     #1st points
#     Z.append([centre[0] ,centre[1]-W/2])
    
#     #2nd points
#     al += dl
#     a1, b1 = centre[0], centre[1] - dl
#     Z.append([a1, b1])
    
#     #3rd points
    
    
#     return U

# L = 50e-6
# H = 50e-6
# W = 1e-6
# l = 100e-6
# dh = 5e-6
# dl = 20e-6
# w = 50e-9

# Ushape = OPKI(L, H, W, l, dh, dl, w)

# carac = {"layer": 0, "datatype": 3}
# U = gdspy.FlexPath(Ushape, W*1e6, **carac)

# Resonator = gdspy.Cell("Resonator")
# Resonator.add(U)


# lib = gdspy.GdsLibrary()
# lib.add(Resonator)
# lib.write_gds("OKPI1.gds")