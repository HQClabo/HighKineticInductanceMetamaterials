# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 17:17:32 2021

@author: jouanny
"""

import numpy as np
import gdspy
from Modules.ModuleResonator import OPKI_up,unit_cell, unit_cell_GGG, unit_cell_upup,ghosts_U, ghosts_U_GGG, T_feedline_simulation, T_feedline_extended, waveguide_simulation, waveguide_extended_new

def ghost_feedline_upup(Q, unitcell_size,startM, startU, stopU, strip_height, tw, tv, M, U, t = 2e-6, Sg = 10e-6,T=100e-6,R=50e-6,f=24e-6, A = 50e-6,B=60e-6, w=4e-6,cozy = False):
    tw = tw*1e6
    tv = tv*1e6
    t = t*1e6
    Sg = Sg*1e6 #lateral spacing to ground planes
    T = T*1e6 #thickness of ground planes
    R = R*1e6 #ground plane width
    f = f*1e6
    A = A*1e6
    B = B*1e6
    w = w*1e6#start width feedline
    
    M[-1][1] = M[-2][1] + 1/2*(M[-1][1]-M[-2][1])
    
    if cozy == True:
        x3,y3 = startM[0] - 3*A/2 - tw - Sg, startU[1] - B - Sg
        x4,y4 = startM[0] + A/2 + Q*unitcell_size[0] + Sg, startM[1]
    else:
        x3, y3 = startM[0] - 3*A/2 - tw - Sg, startU[1] - B - Sg
        x4, y4 =startM[0] + A/2 + Q*unitcell_size[0] + Sg, startM[1] + f
    
    #(x1,y1) and (x2,y2): corners of outer box, (x3,y3) and (x4,y4): corners of window
    
    x1, y1 = x3 - R, y3 - T
    x2, y2 = x4 + R, y4 + T
    
    width_guide = 11/5 * w

    dyke = (width_guide-w)/2
    
    #draw outer box + window
    window = gdspy.Rectangle([x3,y3], [x4,y4])
    outer_box = gdspy.Rectangle([x1,y1],[x2,y2])
    
    #coordinates/ dimensions for FlexPath: right feedline + guide (etched away)
    Y0 = y3 + 1/2*np.abs(y3-y4)
    points = [[x1,Y0],[startM[0]-3*A/2-tw,Y0]]
    
    #enclosing contact pad/patch to define it/ this big_patch + guide will be etched away
    guide = gdspy.FlexPath([points[0],points[1]],width_guide)
    
    
    #draw the feedline (not etched away)
    feedline = gdspy.FlexPath([points[0],points[1]],w)
    ghost_U = gdspy.FlexPath(U, t).translate(startM[0]-A-tw,0)
    ghost_M = gdspy.FlexPath(M,0.5).translate(startM[0]-A-tw,0)
    
    feedline = gdspy.boolean(feedline,ghost_U,'or')
    feedline = gdspy.boolean(feedline,ghost_M, 'or')

    
    #coordinates for left feedline + guide (FlexPath)
    points2 = [[x2,Y0],[startM[0]-A/2+Q*unitcell_size[0]+A,Y0]]
    
    #enclose contact pad with etched line (right edge)
    guide2 = gdspy.FlexPath([points2[0],points2[1]],width_guide)  
    feedline2 = gdspy.FlexPath([points2[0],points2[1]],w)
    ghost_U2 = gdspy.FlexPath(U, t).translate(startM[0]+Q*unitcell_size[0],0)
    ghost_M2 = gdspy.FlexPath(M,0.5).translate(startM[0]+Q*unitcell_size[0],0)
    
    feedline2 = gdspy.boolean(feedline2,ghost_U2,'or')
    feedline2 = gdspy.boolean(feedline2,ghost_M2, 'or')
    
    
    outer_box = gdspy.boolean(outer_box, guide, 'not')
    outer_box = gdspy.boolean(outer_box, guide2, 'not')
    ground = gdspy.boolean(outer_box,window,'not')
    ground_with_feed = gdspy.boolean(ground,feedline,'or')
    ground_with_feeds = gdspy.boolean(ground_with_feed,feedline2,'or')
    
    return ground_with_feeds


"""
~~~~~~ parameters to choose ~~~~~~
"""
BIG = False                    #use twice the size-parameters
Tfeed = True                 #T-feedline (True) or smooth feedline (False)
GGG = False                 # GGG = separation to ground strip constant instead of ground strip width
ghostfeed = False             #feedline with ghost resonator galvanically connected to it
upup = False

if BIG == True:
    L = 567e-6                          #total length of inductor in SI units
    s = 7e-6                            #interspacing of inductor in SI units
    w = 3e-6                            #width of inductor wire
    Ac = 100e-6                          #horizontal dimension of capacitor
    tc = 9e-6                           #thickness of capacitor plates
    tv = 36e-6                          #intracell spacing
    tw = 36e-6                          #intercell spacing
    tsg = 6e-6
    k = 1/6                             #fraction determining how thick the ground strip between resonators is : k*min(tv,tw)
    ep = 42e-6                          #horizontal dimension of ground patches
    fp = 48e-6                          #vertical dimension of ground patches
    rp = 58e-6                          #vertical spacing of resonator "head" to ground plane
    Sgg = 120e-6                        #spacing last ghost to ground
    Tg = 200e-6                         #thickness of ground planes on each side of the array
    Df = 250e-6                         #from where grounding of feedline widens
    Rg = 400e-6                         #extent of ground planes on each side of the array
    wfT = [360e-6,90e-6,20e-6]          #width feedline for T-geometry
    wfs = [360e-6,90e-6,20e-6,4e-6]     #width feedline for smooth feedline
    Q = 16                               #number of unit cells
    N_ghost = 2                         #number of ghosts: For no ghosts, put zero
    if Tfeed == False:
        Sf2r = [0.0,4.5e-6]             #spacing to feedline in [x,y]-direction
    else:
        Sf2r = [20e-6,0]                #spacing to feedline in [x,y]-direction         #spacing to feedline in [x,y]-direction
else:
    L = 279e-6              #total length of inductor in SI units
    s = 4.5e-6              #interspacing of inductor in SI units
    w = 0.47e-6              #width of inductor wire
    Ac = 50e-6               #horizontal dimension of capacitor
    tc = 2e-6                #thickness of capacitor plates
    tv = 24e-6              #intracell spacing
    tw = 12e-6              #intercell spacing
    k = 1/4                #fraction determining how thick the ground strip between resonators is : k*min(tv,tw)
    ep = 20.5e-6             #horizontal dimension of ground patches
    fp = 24e-6               #vertical dimension of ground patches
    rp = 29e-6               #vertical spacing of resonator "head" to ground plane
    Sgg = 60e-6                        #spacing last ghost to ground
    Tg = 100e-6                         #thickness of ground planes on each side of the array
    Df = 250e-6                         #from where grounding of feedline widens
    Rg = 200e-6                         #extent of ground planes on each side of the array
    wfT = [360e-6,90e-6,10e-6]          #width feedline for T-geometry
    wfs = [360e-6,90e-6,10e-6,2e-6]     #width feedline for smooth feedline
    Q = 10                   #number of unit cells
    N_ghost = 2             #number of ghosts: For no ghosts, put zero
    if Tfeed==True:
        Sf2r = [10e-6,0.0]
    else:
        Sf2r = [0.0,2e-6]       #spacing to feedline in [x,y]-direction 

if Tfeed==True or ghostfeed ==True:
    N_ghost = 0

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

if GGG == True:
    unitcell, ground, startM,startU, stopUy, strip_height,Bc, center1st, center2nd = unit_cell_GGG(L, s, w, Ac, tc, tv,tw,ep,fp,rp,carac = carac_res,ts = tsg, ground_in_between = ground_yn)
elif upup == True:
    unitcell, ground, startM, startU, stopUy, strip_height, Bc, center1st, center2nd = unit_cell_upup(L,s,w,Ac,tc,tv,tw,e=ep,f=fp,r=rp,carac=carac_res,gamma=k,ground_in_between = ground_yn)
else:
    unitcell, ground, startM,startU,stopUy,strip_height, Bc, center1st, center2nd = unit_cell(L,s,w,Ac,tc,tv,tw,e=ep,f=fp,r=rp,carac=carac_res, gamma = k, ground_in_between=ground_yn)
unitcell_size = [2*Ac*1e6 + tv*1e6 + tw*1e6,0]
stopU = [startU[0]-Ac*1e6 + Q*unitcell_size[0]-tw*1e6, stopUy]

centerQth = [center1st[0]-Ac*1e6 + Q*unitcell_size[0] - tw*1e6, center2nd[1]]

design = gdspy.Cell('Resonator Array')

ResArray = gdspy.CellArray(unitcell, Q, 1, unitcell_size)
design.add(ResArray)

if ghostfeed == True:
    GrdArray = gdspy.CellArray(ground, Q, 1, unitcell_size)
else:
    GrdArray = gdspy.CellArray(ground, Q-1, 1, unitcell_size)
#add ground to Qth unitcell
#corners of 1st ground patch: down
box_x1, box_y1 = centerQth[0] - tv*1e6 - Ac*1e6 - ep*1e6/2 - w*1e6/2, startM[1]
box_x2, box_y2 = box_x1 + ep*1e6, startM[1]-fp*1e6
ground_box1 = gdspy.Rectangle([box_x1,box_y1], [box_x2,box_y2])
#corners of 2nd ground patch: up
box_x3, box_y3 = centerQth[0] - ep*1e6/2 - w*1e6/2, center1st[1]+rp*1e6
box_x4, box_y4 = box_x3 + ep*1e6, box_y3 - fp*1e6
ground_box2 = gdspy.Rectangle([box_x3,box_y3], [box_x4,box_y4])
groundingQth = gdspy.boolean(ground_box1, ground_box2, 'or')
if ground_yn == True:
    #2nd last ground strip
    if GGG == True:
         strip_x1, strip_y1 = stopU[0] + tc*1e6/2 - Ac*1e6 - tsg*1e6, box_y2
         strip_x2, strip_y2 = strip_x1 - (tv*1e6 - 2*tsg*1e6), strip_y1 + strip_height    
         if ghostfeed == True:
             strip1_x1,strip1_y1 = startM[0] - Ac/2*1e6 - tsg*1e6, startM[1] - fp*1e6
             strip1_x2,strip1_y2 = strip1_x1 - (tw*1e6 - 2*tsg*1e6), startU[1] - Bc - Sgg*1e6
             gr_strip1 = gdspy.Rectangle([strip1_x1,strip1_y1], [strip1_x2,strip1_y2])   
    else:
        strip_x1, strip_y1 = stopU[0] + tc*1e6/2 - Ac*1e6 - 0.5*(tv*1e6 - k*min(tv,tw)*1e6), box_y2
        strip_x2, strip_y2 = strip_x1 - k*min(tv,tw)*1e6, strip_y1 + strip_height
        if ghostfeed == True:
            strip1_x1,strip1_y1 = startM[0] - Ac/2*1e6 - 0.5*(tw*1e6 - k*min(tv,tw)*1e6), startM[1] +fp*1e6
            strip1_x2,strip1_y2 = strip1_x1 - k*min(tv,tw)*1e6, startU[1] - Bc - Sgg*1e6
            gr_strip1 = gdspy.Rectangle([strip1_x1,strip1_y1], [strip1_x2,strip1_y2])
   
    ground_strip_2nd_last = gdspy.Rectangle([strip_x1,strip_y1], [strip_x2,strip_y2])
    groundingQth = gdspy.boolean(groundingQth, ground_strip_2nd_last, 'or', **carac_ground)


if Tfeed == True:
    groundplane_coords = T_feedline_simulation(Q, unitcell_size, startM, startU, stopU, strip_height, tw, tv,t=tc,Sr=Sf2r,D=Df,Sg = Sgg/2, T = Tg, R = Rg, f=fp,A=Ac,B=Bc*1e-6, w= wfT, cozy = False)
    # groundplane_coords = T_feedline_extended(Q, unitcell_size, startM, startU, stopU, strip_height, tw, tv, N_ghost,t=tc,Sr=Sf2r,D=Df,Sg=Sgg,T=Tg,R=Rg,f=fp,A=Ac,B=Bc*1e-6,laserwriter=True)
elif ghostfeed == True:
    U_g, M_g, Bc = OPKI_up(L,s,w,Ac,tc,compact = False)
    groundplane_coords = ghost_feedline_upup(Q,unitcell_size,startM,startU,stopU,strip_height,tw,tv,M_g,U_g,t=tc,Sg=fp,T=Tg,R=Rg,f=fp,A=Ac,B=Bc*1e-6,w=wfT[2], cozy = True)
else:
    groundplane_coords = waveguide_simulation(Q, unitcell_size, startM, startU, stopU, strip_height, tw, tv, N_ghost, t=tc,Sr=Sf2r,Sg=Sgg,T=Tg,R=Rg,f=fp,A=Ac,w_end=wfs[3],cozy=True)
    # groundplane_coords = waveguide_extended_new(Q, unitcell_size, startM, startU, stopU, strip_height, tw, tv, N_ghost,t=tc,Sr=Sf2r,Sg=Sgg,T=Tg,R=Rg,f=fp,A=Ac,w=wfs)
GrdArray = gdspy.boolean(GrdArray, groundingQth, 'or')
if ghostfeed == True:
    GrdArray = gdspy.boolean(GrdArray,gr_strip1,'or')
ground_plane = gdspy.boolean(groundplane_coords,GrdArray,'or', **carac_ground)

if GGG == True:
    ghosts_coords, ground_ghosts = ghosts_U_GGG(L,s,w,Ac,tc,tv,tw,N_ghost, strip_height, tg = tw, e=ep, f=fp, r=rp, ts=tsg, ground_in_between=ground_yn, carac = carac_ground, center_first = center1st, center_last = centerQth)
else:
    ghosts_coords, ground_ghosts = ghosts_U(L,s,w,Ac,tc,tv,tw,N_ghost,strip_height,tg = tw, e = ep, f=fp, r=rp, gamma=k, ground_in_between=ground_yn, carac = carac_ground, center_first = center1st, center_last = centerQth)


if Tfeed == False and ghostfeed==False:
    ground_plane = gdspy.boolean(ground_plane, ground_ghosts, 'or', **carac_ground)
if N_ghost > 0:
    design.add(ghosts_coords)
    
design.add(ground_plane)

design.flatten()

lib = gdspy.GdsLibrary()
lib.add(design)
gdspy.LayoutViewer()
# lib.write_gds("16LC-Array_trivial_Tfeed.gds")

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