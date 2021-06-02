# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 17:17:32 2021

@author: jouanny
"""

import numpy as np
import gdspy
from Modules.ModuleResonator import *





def T_feedline_extended_negative_upup(Q, unitcell_size,startM, startU, stopU, strip_height, tw, tv, N_ghost, t = 2e-6, Sr = [1e-6,2e-6], Sg = 60e-6,T=100e-6,D=200e-6,R=200e-6,f=24e-6, A = 50e-6,B=60e-6, w=[360e-6,90e-6,10e-6],maskdim = [150e-6,300e-6],laserwriter=False,cozy = False):
    tw = tw*1e6
    tv = tv*1e6
    t = t*1e6
    Sr = Sr[0]*1e6,Sr[1]*1e6 #spacing feedline to resonator [x,y]-direction
    Sg = Sg*1e6 #lateral spacing to ground planes
    T = T*1e6 #thickness of ground planes
    D = D*1e6 #spacing to array region
    R = R*1e6 #ground plane width
    # E = 800
    chip_length = 6200
    f = f*1e6
    A = A*1e6
    B = B*1e6
    w_patch = w[0]*1e6
    w_core = w[1]*1e6
    w_start = w[2]*1e6#start width feedline
    maskdim = maskdim[0]*1e6, maskdim[1]*1e6
    # d_high_prec = 150 #distance from ground s.t. everything beyond is in the low-precision layer
    d_overlap = 1 #overlap between the two precision layers
    
    if cozy == True:
        x3,y3 = startM[0] - A/2 - Sr[0] - 2*w_start, startM[1]-strip_height + f
        x4,y4 = x3 + 4*w_start + 2*Sr[0] + Q*unitcell_size[0] - tw, startM[1] + f
    else:
        x3, y3 = startM[0] - A/2 - Sr[0] - 2*w_start, startM[1] - strip_height + f
        x4, y4 = x3 + 2 * Sr[0] + 4*w_start + Q*unitcell_size[0]-tw, startM[1]
    
    #(x1,y1) and (x2,y2): corners of outer box, (x3,y3) and (x4,y4): corners of window
    
    x1, y1 = x3 - R, startM[1] + T
    x2, y2 = x4 + R, y1 - 2*T - strip_height + 2*f   
    
    width_feedline = [w_patch,w_core,w_start]
    width_guide = [73/45*w for w in width_feedline]
    if laserwriter == False:
        width_guide[1] = 31/30*width_feedline[1]
        width_guide[2] = 11/5 * width_feedline[2]
    else:
        width_guide[1] = 16/15*width_feedline[1]
        width_guide[2] = 11/5 * width_feedline[2]

    bend_radii = [3*w for w in width_guide]
    dyke = [(width_guide[0]-width_feedline[0])/2,(width_guide[1]-width_feedline[1])/2,(width_guide[2]-width_feedline[2])/2]
    
    
    E = 0.5*(chip_length - (x2-x1)) - (width_guide[0]-width_feedline[0])/2
    
    #draw outer box + window
    # big_plane = gdspy.Rectangle([x1,y1], [x2,y2])
    window = gdspy.Rectangle([x3,y3], [x4,y4])
    
    #coordinates/ dimensions for FlexPath: right feedline + guide (etched away)
    Y0 = y3 + 1/2*np.abs(y3-y4)#startM[1] - 1/2*(strip_height )#- 2*f)
    points = [[x1-E,Y0],[x1-7*E/8, Y0],[x1-3*E/4, Y0],[x1,Y0],[x3 - D, Y0],[startU[0] + t/2 - A - Sr[0] - w_start/2, Y0]]
    
    #enclose contact pad with etched line (right edge)
    patch_start = [points[0][0]- dyke[0],points[0][1]]
    # patch_end = [points[2][0]+1,points[2][1]]
    
    #enclosing contact pad/patch to define it/ this big_patch + guide will be etched away
    big_patch = gdspy.FlexPath([patch_start,points[1]], width_guide[0])#,corners=#'circular bend',bend_radius = bend_radii[1])
    big_patch = big_patch.segment(points[2],width_guide[1])
    guide = gdspy.FlexPath([points[2],points[4]],width_guide[1]).segment(points[5],width_guide[2])
    T_bar_guide = gdspy.Rectangle([points[5][0] - w_start/2 - dyke[2], startU[1] - B + dyke[2]], [points[5][0] + w_start/2 + dyke[2], startU[1] - dyke[2]])
    
    guide = gdspy.boolean(guide, T_bar_guide, 'or')
    
    #draw the feedline (not etched away)
    patch = gdspy.FlexPath([points[0],points[1]], width_feedline[0])#,corners='circular bend',bend_radius = bend_radii[1])
    patch = patch.segment(points[2],width_feedline[1])
    feedline = gdspy.FlexPath([points[2],points[4]],width_feedline[1]).segment(points[5],width_feedline[2])#,corners='circular bend',bend_radius=bend_radii[2])
    T_bar = gdspy.Rectangle([points[5][0] - w_start/2, startU[1] - B], [points[5][0] + w_start/2, startU[1]])
    
    feedline = gdspy.boolean(feedline, T_bar, 'or')
    
    pad = gdspy.boolean(big_patch,patch, 'not')
    line = gdspy.boolean(guide,feedline, 'not')
    
    #coordinates for left feedline + guide (FlexPath)
    # Y02 = y2 + 1/2*(2*T + strip_height - 2*f)
    points2 = [[x2+E,Y0],[x2+7*E/8, Y0],[x2+3*E/4, Y0],[x2,Y0],[x4 + D, Y0],[stopU[0] + t/2 + Sr[0] + w_start/2,Y0]]
    
    #enclose contact pad with etched line (right edge)
    patch_start2 = [points2[0][0]+dyke[0],points2[0][1]]
    # patch_end2 = [points2[2][0]-1,points2[2][1]]
    
    big_patch2 = gdspy.FlexPath([patch_start2,points2[1]], width_guide[0])#,corners='circular bend',bend_radius = bend_radii[1])
    big_patch2 = big_patch2.segment(points2[2],width_guide[1])
    guide2 = gdspy.FlexPath([points2[2],points2[4]],width_guide[1]).segment(points2[5],width_guide[2])#,corners='circular bend',bend_radius=bend_radii[2])
    T_bar_guide2 = gdspy.Rectangle([points2[5][0] - w_start/2-dyke[2], stopU[1]+dyke[2]], [points2[5][0] - w_start/2 - dyke[2], stopU[1]-B-dyke[2]])
    
    guide2 = gdspy.boolean(guide2, T_bar_guide2, 'or')
    
    patch2 = gdspy.FlexPath([points2[0],points2[1]], width_feedline[0])#,corners='circular bend',bend_radius = bend_radii[1])
    patch2 = patch2.segment(points2[2],width_feedline[1])   
    feedline2 = gdspy.FlexPath([points2[2],points2[4]],width_feedline[1]).segment(points2[5],width_feedline[2])#,corners='circular bend',bend_radius=bend_radii[2])
    T_bar2 = gdspy.Rectangle([points2[5][0] - w_start/2, stopU[1]], [points2[5][0] + w_start/2, stopU[1]-B])

    feedline2 = gdspy.boolean(feedline2, T_bar2, 'or')

    pad2 = gdspy.boolean(big_patch2, patch2, 'not')
    line2 = gdspy.boolean(guide2,feedline2, 'not')
    
    window_left_guide = gdspy.boolean(window, line, 'or')
    window_both_guides = gdspy.boolean(window_left_guide, line2, 'or')
    window_left_feedline = gdspy.boolean(window_both_guides,feedline ,'not')
    window_both_feedlines = gdspy.boolean(window_left_feedline,feedline2,'not')
    pads = gdspy.boolean(pad,pad2,'or')
    # window_both_feedlines = gdspy.boolean(window_both_feedlines,pads,'or')
    
    mask = gdspy.Rectangle([x1+R-maskdim[0],y1+2*R-maskdim[1]], [x2-R+maskdim[0],y2-2*R+maskdim[1]])
    mask_overlap = gdspy.Rectangle([x1+R - maskdim[0] - d_overlap, y1+2*R-maskdim[1]-d_overlap], [x2-R+maskdim[0] + d_overlap,y2-2*R+maskdim[1] + d_overlap])
    
    return window_both_feedlines,pads, mask, mask_overlap
"""
~~~~~~ parameters to choose ~~~~~~
"""

BIG = False                    #use twice the size-parameters
Tfeed = True                 #T-feedline (True) or smooth feedline (False)
laseryn = False


if BIG == True:
    L = 567e-6                          #total length of inductor in SI units
    s = 7e-6                            #interspacing of inductor in SI units
    w = 3e-6                            #width of inductor wire
    Ac = 100e-6                          #horizontal dimension of capacitor
    tc = 4e-6                           #thickness of capacitor plates
    tv = 12e-6                          #intracell spacing
    tw = 48e-6                          #intercell spacing
    tsg = 6e-6
    k = 1/2                             #fraction determining how thick the ground strip between resonators is : k*min(tv,tw)
    ep = 42e-6                          #horizontal dimension of ground patches
    fp = 48e-6                          #vertical dimension of ground patches
    rp = 58e-6                          #vertical spacing of resonator "head" to ground plane
    Sgg = 120e-6                        #spacing last ghost to ground
    Tg = 200e-6                         #thickness of ground planes on each side of the array
    Df = 500e-6                         #from where grounding of feedline widens
    Rc = 400e-6                         #extent of ground planes on each side of the array
    wfT = [360e-6,90e-6,10e-6]          #width feedline for T-geometry
    wfs = [360e-6,90e-6,20e-6,4e-6]     #width feedline for smooth feedline
    Q = 2                               #number of unit cells
    N_ghost = 2                         #number of ghosts: For no ghosts, put zero
    if Tfeed == False:
        Sf2r = [0.0,4.5e-6]             #spacing to feedline in [x,y]-direction
    else:
        Sf2r = [20e-6,0]                #spacing to feedline in [x,y]-direction
else:
    L = 227e-6              #total length of inductor in SI units
    s = 4.5e-6              #interspacing of inductor in SI units
    w = 0.5e-6              #width of inductor wire
    Ac = 38e-6               #horizontal dimension of capacitor
    tc = 2e-6                #thickness of capacitor plates
    tv = 18e-6              #intracell spacing
    tw = 18e-6              #intercell spacing
    k = 1/6                 #fraction determining how thick the ground strip between resonators is : k*min(tv,tw)
    ep = 20.5e-6             #horizontal dimension of ground patches
    fp = 14e-6#24e-6               #vertical dimension of ground patches
    rp = 29e-6               #vertical spacing of resonator "head" to ground plane
    Sgg = 60e-6                        #spacing last ghost to ground
    Tg = 100e-6                         #thickness of ground planes on each side of the array
    Df = 250e-6                         #from where grounding of feedline widens
    Rc = 200e-6                         #extent of ground planes on each side of the array
    wfT = [360e-6,90e-6,10e-6]          #width feedline for T-geometry
    wfs = [360e-6,90e-6,10e-6,2e-6]     #width feedline for smooth feedline
    Q = 12                   #number of unit cells
    N_ghost = 2             #number of ghosts: For no ghosts, put zero
    if Tfeed==True:
        Sf2r = [10e-6,0.0]
    else:
        Sf2r = [0.0,2e-6]       #spacing to feedline in [x,y]-direction 

if Tfeed==True:
    N_ghost = 0


ground_yn = True        #ground in between yes (True) or no (False) 
carac_lowPres = {'layer' : 2, 'datatype' : 0}       
carac_highPres = {'layer' :  4, 'datatype' : 0}
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

unitcell, ground, startM,startU, stopUy, strip_height, Bc, center1st, center2nd = unit_cell_upup(L, s, w, Ac, tc, tv,tw,ep,fp,rp,gamma = k, ground_in_between = ground_yn)
# unitcell, ground, startM,startU, stopUy, strip_height, Bc, center1st, center2nd = unit_cell_GGG(L, s, w, Ac, tc, tv,tw,ep,fp,rp,ts = tsg, ground_in_between = ground_yn)
unitcell_size = [2*Ac*1e6 + tv*1e6 + tw*1e6,0]
stopU = [startU[0]-Ac*1e6 + Q*unitcell_size[0]-tw*1e6, stopUy]

centerQth = [center1st[0]-Ac*1e6 + Q*unitcell_size[0] - tw*1e6, center2nd[1]]

design = gdspy.Cell('Resonator Array')

ResArray = gdspy.CellArray(unitcell, Q, 1, unitcell_size)
# ResArr_cell = gdspy.Cell('bliblub')
# ResArr_cell.add(ResArray)
# ResArr_cell.flatten()

GrdArray = gdspy.CellArray(ground, Q-1, 1, unitcell_size)
#add ground to Qth unitcell
#corners of 1st ground patch: down
box_x1, box_y1 = centerQth[0] - tv*1e6 - Ac*1e6 - ep*1e6/2 - w*1e6/2, startM[1]
box_x2, box_y2 = box_x1 + ep*1e6, startM[1]+fp*1e6
# box_x2, box_y2 = box_x1 + ep*1e6, startM[1]-fp*1e6
ground_box1 = gdspy.Rectangle([box_x1,box_y1], [box_x2,box_y2])
#corners of 2nd ground patch: up
box_x3, box_y3 = centerQth[0] - ep*1e6/2 - w*1e6/2, startM[1] #center1st[1]+rp*1e6 
box_x4, box_y4 = box_x3 + ep*1e6, box_y3 + fp*1e6#- fp*1e6
ground_box2 = gdspy.Rectangle([box_x3,box_y3], [box_x4,box_y4])
groundingQth = gdspy.boolean(ground_box1, ground_box2, 'or')
if ground_yn == True:
    #2nd last ground strip
    # strip_x1, strip_y1 = stopU[0] + tc*1e6/2 - Ac*1e6 - tsg*1e6, box_y2
    # strip_x2, strip_y2 = strip_x1 - (tv*1e6-2*tsg*1e6), strip_y1 + strip_height
    strip_x1, strip_y1 = stopU[0] + tc*1e6/2 - Ac*1e6 - 0.5*(tv*1e6 - k*min(tv,tw)*1e6), box_y2
    strip_x2, strip_y2 = strip_x1 - k*min(tv,tw)*1e6, strip_y1 - strip_height#+ strip_height
    ground_strip_2nd_last = gdspy.Rectangle([strip_x1,strip_y1], [strip_x2,strip_y2])
    groundingQth = gdspy.boolean(groundingQth, ground_strip_2nd_last, 'or', **carac_ground)


if Tfeed == True:
    groundplane_coords,pads,mask, mask_overlap = T_feedline_extended_negative_upup(Q, unitcell_size,startM, startU, stopU, strip_height, tw, tv, N_ghost, t = tc, Sr = Sf2r, Sg = Sgg,T=Tg,D=Df,R=Rc,f=fp, A = Ac,B=Bc*1e-6, w=wfT,laserwriter=laseryn,maskdim = [100e-6,300e-6])
else:
    groundplane_coords,pads,mask, mask_overlap = waveguide_extended_negative_new(Q,unitcell_size,startM,startU,stopU,strip_height,tw,tv,N_ghost,t=tc, Sr = Sf2r, Sg = Sgg, T = Tg, R = Rc, f = fp, A = Ac, w = wfs,laserwriter=laseryn,maskdim = [240e-6,250e-6])

GrdArray = gdspy.boolean(GrdArray, groundingQth, 'or')
ground_plane = gdspy.boolean(groundplane_coords,GrdArray,'not', **carac_highPres)

if BIG == False:
    ghosts_coords, ground_ghosts = ghosts_U(L, s, w, Ac, tc, tv,tw, N_ghost,strip_height,tg = tw,gamma=k, center_first=center1st,center_last=centerQth,ground_in_between=(ground_yn))
    # ghosts_coords, ground_ghosts = ghosts_U_GGG(L, s, w, Ac, tc, tv,tw, N_ghost,strip_height,tg = tw,ts=tsg, center_first=center1st,center_last=centerQth,ground_in_between=(ground_yn))
else:
    ghosts_coords, ground_ghosts = ghosts_U(L, s, w, Ac, tc, tv,tw, N_ghost,strip_height,tg = tw,e=ep, f=fp, r=rp, gamma=k, ground_in_between=(ground_yn), center_first = center1st, center_last = centerQth)
    # ghosts_coords, ground_ghosts = ghosts_U_GGG(L, s, w, Ac, tc, tv,tw, N_ghost,strip_height,tg = tw,e=ep, f=fp, r=rp, ts=tsg, ground_in_between=(ground_yn), center_first = center1st, center_last = centerQth)


if Tfeed == False:
    ground_plane = gdspy.boolean(ground_plane, ground_ghosts, 'not', **carac_highPres)

negative = gdspy.boolean(ground_plane,ResArray, 'not',**carac_highPres)

layer4 = gdspy.boolean(negative,mask_overlap,'and', layer=4)
layer2_line = gdspy.boolean(negative,mask, 'not')
layer2 = gdspy.boolean(layer2_line,pads, 'or', layer=2)

#just to get the layer properties to the pads
# pads = gdspy.boolean(pads, groundingQth, 'not', **carac_lowPres)


if N_ghost > 0:
    # bla = gdspy.boolean(negative, ghosts_coords, 'not', **carac_ground)
    design.add(ghosts_coords)

design.add(layer4)
design.add(layer2)

lib = gdspy.GdsLibrary()
lib.add(design)
gdspy.LayoutViewer()
lib.write_gds("22LC-Array_normal.gds")
