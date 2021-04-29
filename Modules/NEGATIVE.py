# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 09:28:11 2021

@author: vweibel
"""
import gdspy
import numpy as np

def OPKI_down(L,s,w,A,t, centre=(0,0),return_center=False):
    L = L*1e6 #total length of inductor
    s = s*1e6 #interspacing of turns of inductor
    w = w*1e6 #width of inductor wire
    A = A*1e6 #horizontal dimension of U-capacitor
    t = t*1e6 #thickness of U-capacitor
    centre = centre[0], centre[1]
    
    #calculate horizontal dimension of inductor
    b = 3/5*A - 2*t
    #print('b: ',b)
    
    #calculate number of windings
    N = int((L-2*(s+w))/(b+s))
    # N = round((L-2*(s+w))/(b+s))
    print('N: ',N)
    
    #adapt start & end segment of inductor
    d_prime = 0.5*(L - N*(s+b)) + w
    #print('d_prime: ',d_prime)
    
    #calculate vertical dimension of capacitor
    B = N*(s+w) + d_prime + t
    # B = N*(s+w) + 7/4*d_prime + t
   # print('B: ', B)
    
    #Ushape, we take the centre at at the bottom/top of the U
    U=[]
    x0, y0 = centre[0],centre[1]
    U.append([x0, y0])
    
    x1, y1 = x0, y0 + B - t/2
    U.append([x1, y1])
    
    x2, y2 = x1 + A - t, y1
    U.append([x2, y2])
    
    x3, y3 = x2, y2 - B + t/2
    U.append([x3, y3])
    
    U = np.asarray(U)
    centring = (A-t)/2, B
    U = U - centring
    
    #draw meander from centre of U
    M = []
    v0,w0 = centre[0], centre[1]-t
    M.append([v0,w0])
    v1,w1 = v0, w0 - d_prime + w/2
    M.append([v1,w1])
    
    for i in range(N):
        if i%2 == 0:
            v2, w2 = M[-1][0]  - b/2 + w/2, M[-1][1]
            M.append([v2,w2])
            v3, w3 = v2, w2 - (s+w)
            M.append([v3,w3])
            v4, w4 = v3 + (b/2 - w/2), w3
            M.append([v4,w4])
            # print(i)
        if (i+1)%2==0:
            v5,w5 = M[-1][0] + b/2 - w/2, M[-1][1]
            M.append([v5,w5])
            v6, w6 = v5, w5 - (s+w)
            M.append([v6,w6])
            v7,w7 = v6 - b/2 + w/2, w6
            M.append([v7,w7])
            # print(i)
        if (i+1)==N:
            # print('end turns:', M[-1])
            v8,w8 = M[-1][0], M[-1][1]
            v9, w9 = v8, w8 - d_prime + w/2
            M.append([v9,w9])
            # print('end meander:', M[-1])
        # else:
            # print('not entered')
    
    if return_center == False:
        return U,M
    if return_center == True:
        return U, M, centre

def OPKI_up(L,s,w,A,t, centre=(0,0),return_center = False):
    L = L*1e6 #total length of inductor
    s = s*1e6 #interspacing of turns of inductor
    w = w*1e6 #width of inductor wire
    A = A*1e6 #horizontal dimension of U-capacitor
    t = t*1e6 #thickness of U-capacitor
    # centre = centre[0]*1e6, centre[1]*1e6
    
    #calculate horizontal dimension of inductor
    b = 3/5*A - 2*t
    # print('b: ',b)
    
    #calculate number of windings
    # N = round((L-2*(s+w))/(b+s))
    N = int((L-2*(s+w))/(b+s))
    # print('N: ',N)
    
    #adapt start & end segment of inductor
    d_prime = 0.5*(L - N*(s+b)) + w
    # print('d_prime: ',d_prime)
    
    #calculate vertical dimension of capacitor
    B = N*(s+w) + d_prime + t
    # B = N*(s+w) + 7/4*d_prime + t
    # print('B: ', B)
    
    #Ushape, we take the centre at at the bottom/top of the U
    U=[]
    x0, y0 = centre[0],centre[1]
    U.append([x0, y0])
    
    x1, y1 = x0, y0 - B + t/2
    U.append([x1, y1])
    
    x2, y2 = x1 + A - t, y1
    U.append([x2, y2])
    
    x3, y3 = x2, y2 + B - t/2
    U.append([x3, y3])
    
    U = np.asarray(U)
    centring = (A-t)/2, -B
    U = U - centring
    
    #draw meander from centre of U
    M = []
    v0,w0 = centre[0], centre[1]+t
    M.append([v0,w0])
    v1,w1 = v0, w0 + d_prime - w/2
    M.append([v1,w1])
    
    for i in range(N):
        if i%2 == 0:
            v2, w2 = M[-1][0] - b/2 + w/2, M[-1][1]
            M.append([v2,w2])
            v3, w3 = v2, w2 + (s+w)
            M.append([v3,w3])
            v4, w4 = v3 + (b/2 - w/2), w3
            M.append([v4,w4])
            # print(i)
        if (i+1)%2==0:
            v5,w5 = M[-1][0] + b/2 - w/2, M[-1][1]
            M.append([v5,w5])
            v6, w6 = v5, w5 + (s+w)
            M.append([v6,w6])
            v7,w7 = v6 - b/2 + w/2, w6
            M.append([v7,w7])
            # print(i)
        if (i+1)==N:
            # print('end turns:', M[-1])
            v8,w8 = M[-1][0], M[-1][1]
            v9, w9 = v8, w8 + d_prime - w/2
            M.append([v9,w9])
            # print('end meander:', M[-1])
        # else:
            # print('not entered')
            
    if return_center == False:
        return U,M
    if return_center == True:
        return U, M, centre


"""
------------------------------------------------------------------------------
left-handed LC array: Negative
------------------------------------------------------------------------------
"""

def waveguide_negative(Q, unitcell_size,startM, startU, stopU, strip_height, tw, tv, N_ghost, t = 2e-6, Sr = [1e-6,2e-6], Sg = 60e-6,T=100e-6,R=200e-6,f=24e-6, A = 50e-6):
    tw = tw*1e6
    tv = tv*1e6
    t = t*1e6
    Sr = Sr[0]*1e6,Sr[1]*1e6 #spacing feedline to resonator [x,y]-direction
    Sg = Sg*1e6 #lateral spacing to ground planes
    T = T*1e6 #thickness of ground planes
    R = R*1e6 #width of ground planes
    f = f*1e6
    A = A*1e6
    w_start = 10 #start width feedline
    w_end = 2 #end width feedline
    
    #calculate length ghosts demand
    if (N_ghost-1)%2==0:
        extent_ghosts = tw + tv + N_ghost*A + (N_ghost-1)/2*(tw + tv)
        print('extent ghosts is: ',extent_ghosts, ' (N_ghost odd)')
    else:
        extent_ghosts = tw + tv + N_ghost*A + (N_ghost-1)/2*tw + ((N_ghost-1)/2 - 1)*tv
        print('extent ghosts is: ',extent_ghosts, ' (N_ghost even)')
    
    #(x1,y1) and (x2,y2): corners of outer box, (x3,y3) and (x4,y4): corners of window
    x1, y1 = startM[0] - A/2 - Sg - R - extent_ghosts, startM[1] - f - T
    x2, y2 = x1 + 2*(R + Sg) + Q*unitcell_size[0]-tw + 2*extent_ghosts, y1 + 2*T + strip_height
    x3, y3 = startM[0] - A/2 - Sg - extent_ghosts, startM[1] - f
    x4, y4 = x3 + 2 * Sg + Q*unitcell_size[0]-tw + 2*extent_ghosts, y3 + strip_height
    
    #draw outer box + window
    # big_plane = gdspy.Rectangle([x1,y1], [x2,y2])
    window = gdspy.Rectangle([x3,y3], [x4,y4])
    #subtract window from outer box
    # ground_plane = gdspy.boolean(big_plane,window,'not')
    #draw path guiding the feedline
    points = [[x1, y1 + 2/3*(2*T + strip_height)],[x1 + R/4, y1 + 2/3*(2*T + strip_height)],[x1 + R/2, y1 + 3/5*T],[startU[0] - A + t + Sr[0],y1 + 1/5*T],[startU[0] - A + t + Sr[0],startU[1]-Sr[1]]]
    if Sr[0] > 0:
        points[-2][0] = startU[0] - A + t + Sr[0] + w_end
        points[-1][0] = startU[0] - A + t + Sr[0] + w_end
    width_feedline = [w_start, w_start - 1/4*(w_start - w_end), w_end + 1/4*(w_start-w_end), w_end]
    width_guide = [11/5*w for w in width_feedline]
    bend_radii = [3*w for w in width_guide]
    guide = gdspy.FlexPath([points[0],points[1]],width_guide[0],corners='circular bend',bend_radius=bend_radii[0])
    guide = guide.segment(points[2],width_guide[1]).segment(points[3],width_guide[2]).segment(points[4],width_guide[3]) #
    # ground_plane = gdspy.boolean(ground_plane,guide,'not')
    feedline = gdspy.FlexPath([points[0],points[1]],width_feedline[0],corners='circular bend',bend_radius=bend_radii[0])
    feedline = feedline.segment(points[2],width_feedline[1]).segment(points[3],width_feedline[2]).segment(points[4],width_feedline[3]) #
    
    points2 = [[x2, y2 - 2/3*(2*T + strip_height)],[x2 - R/4, y2 - 2/3*(2*T + strip_height)],[x2 - R/2, y2 - 3/5*T],[stopU[0]+ Sr[0],y2 - 1/5*T],[stopU[0]+ Sr[0],stopU[1]+Sr[1]]]
    if Sr[0] > 0:
        points2[-2][0] = stopU[0]+ Sr[0] + w_end
        points2[-1][0] = stopU[0]+ Sr[0] + w_end
    guide2 = gdspy.FlexPath([points2[0],points2[1]],width_guide[0],corners='circular bend',bend_radius=bend_radii[0])
    guide2 = guide2.segment(points2[2],width_guide[1]).segment(points2[3],width_guide[2]).segment(points2[4],width_guide[3])
    # ground_plane2 = gdspy.boolean(ground_plane,guide2,'not')
    feedline2 = gdspy.FlexPath([points2[0],points2[1]],width_feedline[0],corners='circular bend',bend_radius=bend_radii[0])
    feedline2 = feedline2.segment(points2[2],width_feedline[1]).segment(points2[3],width_feedline[2]).segment(points2[4],width_feedline[3])
    
    
    window_left_guide = gdspy.boolean(window, guide, 'or')
    window_both_guides = gdspy.boolean(window_left_guide, guide2, 'or')
    window_left_feedline = gdspy.boolean(window_both_guides,feedline ,'not')
    window_both_feedlines = gdspy.boolean(window_left_feedline,feedline2,'not')
    
    return window_both_feedlines

def ghosts_negative(L,s,w,A,t,tv,tw,N, strip_height, tg = 15e-6, e=20.5e-6, f=24e-6, r=29e-6, gamma=1/4, ground_in_between=True, carac = {'layer' :  1, 'datatype' : 1}, center_first = (0,0), center_last = (0,0)):
    A = A*1e6
    t = t*1e6
    w = w*1e6
    e = e*1e6
    f = f*1e6
    r = r*1e6
    tv = tv*1e6
    tw = tw*1e6
    tg = tg*1e6
    center_left = center_first
    center_right = center_last
    xw = 0.5*(tw - gamma*min(tw,tv))
    xv = 0.5*(tv - gamma*min(tw,tv))
    xg = 0.5*(tg - gamma*min(tw,tv))
    
    ghostcell = gdspy.Cell('ghosts')
    ground = gdspy.Cell('ghosts ground')
    
    if N == 0 and ground_in_between == True:
        strip_x1_first = center_left[0] - A/2 - xg
        strip_y1_first = center_left[1] + r
        strip_x2_first = strip_x1_first - gamma*min(tw,tv)
        strip_y2_first = strip_y1_first - strip_height
        ground_left = gdspy.Rectangle([strip_x1_first,strip_y1_first], [strip_x2_first,strip_y2_first])
        strip_x1_last = center_last[0] + A/2 + xg
        strip_y1_last = center_last[1] - r
        strip_x2_last = strip_x1_last + gamma*min(tw,tv)
        strip_y2_last = strip_y1_last + strip_height
        ground_right = gdspy.Rectangle([strip_x1_last,strip_y1_last], [strip_x2_last,strip_y2_last])
        ground_edge = gdspy.boolean(ground_left, ground_right, 'or')
        print('N_ghosts is zero.')
        return gdspy.Rectangle([100,100],[200,200]), ground_edge

    for i in range(1,N+1):
        if (i+1)%2 == 0: #odd
            if i == 1:
                center_left = [center_left[0] - (A + tg), center_last[1]]
                center_right = [center_right[0] + (A + tg), center_first[1]]
            else:
                center_left = [center_left[0] - (A+tw), center_last[1]]
                center_right = [center_right[0] + (A+tw), center_first[1]]
            
            Ushape_l, Mshape_l= OPKI_up(L,s,w*1e-6,A*1e-6,t*1e-6,centre=center_left)
            Ushape_r, Mshape_r= OPKI_down(L,s,w*1e-6,A*1e-6,t*1e-6,centre=center_right)
            U_l = gdspy.FlexPath(Ushape_l, t,**carac)
            # M_l = gdspy.FlexPath(Mshape_l,w,**carac)
            # U_r = gdspy.Rectangle(center_left,[center_left[0]+50, center_left[1]+50])
            U_r = gdspy.FlexPath(Ushape_r, t,**carac)
            # M_r = gdspy.FlexPath(Mshape_r,w,**carac)
            ghostcell.add([U_r,U_l])
            #ground patches
            box_x1, box_y1 = Mshape_r[-1][0] - e/2 - w/2, Mshape_r[-1][1]
            box_x2, box_y2 = Mshape_r[-1][0] + w/2 + e/2, Mshape_r[-1][1]-f
            ground_box1 = gdspy.Rectangle([box_x1,box_y1], [box_x2,box_y2])
            box_x3, box_y3 = Mshape_l[-1][0] - e/2 - w/2, Mshape_l[-1][1]
            box_x4, box_y4 = Mshape_l[-1][0] + w/2 + e/2, Mshape_l[-1][1] + f
            ground_box2 = gdspy.Rectangle([box_x3,box_y3], [box_x4,box_y4])
            ground.add([ground_box1,ground_box2])
            
            if i == 1:
                #ground strips
                strip_x1_l = center_left[0] + A/2 + xg
                strip_y1_l = box_y4
                strip_x2_l = strip_x1_l + gamma*min(tw,tv)
                strip_y2_l = strip_y1_l - strip_height
                strip_x1_r = center_right[0] + A/2 + xv
                strip_y1_r = box_y2
                strip_x2_r = strip_x1_r + gamma*min(tw,tv)
                strip_y2_r = strip_y1_r + strip_height
            else:
                #ground strips
                strip_x1_l = center_left[0] + A/2 + xg
                strip_y1_l = box_y4
                strip_x2_l = strip_x1_l + gamma*min(tw,tv)
                strip_y2_l = strip_y1_l - strip_height
                strip_x1_r = center_right[0] + A/2 + xv
                strip_y1_r = box_y2
                strip_x2_r = strip_x1_r + gamma*min(tw,tv)
                strip_y2_r = strip_y1_r + strip_height
            
            if i == N:
                strip_x1_N = center_left[0] - A/2 - xv
                strip_y1_N = box_y4
                strip_x2_N = strip_x1_N - gamma*min(tw,tv)
                strip_y2_N = strip_y1_N - strip_height
                strip_x1_last = center_last[0] + A/2 + xg
                strip_y1_last = center_last[1] - r
                strip_x2_last = strip_x1_last + gamma*min(tw,tv)
                strip_y2_last = strip_y1_last + strip_height
                print('N is odd.')
            
            if ground_in_between == True:
                ground_strip_l = gdspy.Rectangle([strip_x1_l,strip_y1_l], [strip_x2_l,strip_y2_l])
                ground.add(ground_strip_l)
                ground_strip_r = gdspy.Rectangle([strip_x1_r,strip_y1_r], [strip_x2_r,strip_y2_r])
                ground.add(ground_strip_r)
            
            # print('run ',i,'- center_left: ',center_left,' center_right: ',center_right)
        if i%2 == 0: #even
            center_left = [center_left[0] - (A+tv), center_first[1]]
            center_right = [center_right[0] + (A+tv),center_last[1]]
            Ushape_l,Mshape_l = OPKI_down(L, s, w*1e-6,A*1e-6,t*1e-6, centre=center_left)
            Ushape_r,Mshape_r = OPKI_up(L, s, w*1e-6,A*1e-6,t*1e-6, centre = center_right)
            U_l = gdspy.FlexPath(Ushape_l, t,**carac)
            # M_l = gdspy.FlexPath(Mshape_l,w,**carac)
            # U_r = gdspy.Rectangle(center_left,[center_left[0]+50, center_left[1]+50])
            U_r = gdspy.FlexPath(Ushape_r, t,**carac)
            # M_r = gdspy.FlexPath(Mshape_r,w,**carac)
            ghostcell.add([U_r,U_l])
            #ground patches
            box_x1, box_y1 = Mshape_l[-1][0] - e/2 - w/2, Mshape_l[-1][1]
            box_x2, box_y2 = Mshape_l[-1][0] + w/2 + e/2, Mshape_l[-1][1]-f
            ground_box1 = gdspy.Rectangle([box_x1,box_y1], [box_x2,box_y2])
            box_x3, box_y3 = Mshape_r[-1][0] - e/2 - w/2, Mshape_r[-1][1]
            box_x4, box_y4 = Mshape_r[-1][0] + w/2 + e/2, Mshape_r[-1][1] + f
            ground_box2 = gdspy.Rectangle([box_x3,box_y3], [box_x4,box_y4])
            ground.add([ground_box1,ground_box2])
            #ground strips
            strip_x1_l = center_left[0] + A/2 + xv
            strip_y1_l = box_y2
            strip_x2_l = strip_x1_l + gamma*min(tw,tv)
            strip_y2_l = strip_y1_l + strip_height
            strip_x1_r = center_right[0] + A/2 + xw
            strip_y1_r = box_y4
            strip_x2_r = strip_x1_r + gamma*min(tw,tv)
            strip_y2_r = strip_y1_r - strip_height
            
            if i == N:
                strip_x1_N = center_left[0] - A/2 - xv
                strip_y1_N = box_y2
                strip_x2_N = strip_x1_N - gamma*min(tw,tv)
                strip_y2_N = strip_y1_N + strip_height
                strip_x1_last = center_last[0] + A/2 + xg
                strip_y1_last = center_last[1] - r
                strip_x2_last = strip_x1_last + gamma*min(tw,tv)
                strip_y2_last = strip_y1_last + strip_height
                print('N is even.')
            
            if ground_in_between == True:
                ground_strip_l = gdspy.Rectangle([strip_x1_l,strip_y1_l], [strip_x2_l,strip_y2_l])
                ground.add(ground_strip_l)
                ground_strip_r = gdspy.Rectangle([strip_x1_r,strip_y1_r], [strip_x2_r,strip_y2_r])
                ground.add(ground_strip_r)

            
    if ground_in_between == True:
        ground_strip_N = gdspy.Rectangle([strip_x1_N,strip_y1_N], [strip_x2_N,strip_y2_N])
        ground.add(ground_strip_N)
        ground_right = gdspy.Rectangle([strip_x1_last,strip_y1_last], [strip_x2_last,strip_y2_last])
        ground.add(ground_right)
                    
            # print('run ',i,'- center_left: ',center_left,' center_right: ',center_right)
    ghostcell.flatten()     
    return ghostcell,ground

def waveguide_extended_negative(Q, unitcell_size,startM, startU, stopU, strip_height, tw, tv, N_ghost, t = 2e-6, Sr = [1e-6,2e-6], Sg = 60e-6,T=100e-6,R=200e-6,f=24e-6, A = 50e-6):
    tw = tw*1e6
    tv = tv*1e6
    t = t*1e6
    Sr = Sr[0]*1e6,Sr[1]*1e6 #spacing feedline to resonator [x,y]-direction
    Sg = Sg*1e6 #lateral spacing to ground planes
    T = T*1e6 #thickness of ground planes
    R = R*1e6 #width of ground planes
    # E = 800
    chip_length = 6200
    f = f*1e6
    A = A*1e6
    w_patch = 360.0
    w_start = 10.0 #start width feedline
    w_end = 2.0 #end width feedline
    
    
    #calculate length ghosts demand
    if (N_ghost-1)%2==0:
        extent_ghosts = tw + tv + N_ghost*A + (N_ghost-1)/2*(tw + tv)
        print('extent ghosts is: ',extent_ghosts, ' (N_ghost odd)')
    else:
        extent_ghosts = tw + tv + N_ghost*A + (N_ghost-1)/2*tw + ((N_ghost-1)/2 - 1)*tv
        print('extent ghosts is: ',extent_ghosts, ' (N_ghost even)')
    
    #(x1,y1) and (x2,y2): corners of outer box, (x3,y3) and (x4,y4): corners of window
    x1, y1 = startM[0] - A/2 - Sg - R - extent_ghosts, startM[1] - f - T
    x2, y2 = x1 + 2*(R + Sg) + Q*unitcell_size[0]-tw + 2*extent_ghosts, y1 + 2*T + strip_height
    x3, y3 = startM[0] - A/2 - Sg - extent_ghosts, startM[1] - f
    x4, y4 = x3 + 2 * Sg + Q*unitcell_size[0]-tw + 2*extent_ghosts, y3 + strip_height
    
    width_feedline = [w_patch,(w_patch-w_start)/8,w_start, w_start - 1/4*(w_start - w_end), w_end + 1/4*(w_start-w_end), w_end]
    width_guide = [11/5*w for w in width_feedline]
    width_guide[0] = w_patch + 10*w_start
    bend_radii = [3*w for w in width_guide]
    
    
    E = 0.5*(chip_length - (x2-x1)) - (width_guide[0]-width_feedline[0])/2
    
    #draw outer box + window
    big_plane = gdspy.Rectangle([x1,y1], [x2,y2])
    window = gdspy.Rectangle([x3,y3], [x4,y4])
    
    #coordinates/ dimensions for FlexPath: right feedline + guide (etched away)
    points = [[x1-E,y1 + 1/2*(2*T + strip_height)],[x1-3*E/4, y1 + 1/2*(2*T + strip_height)],[x1-E/4, y1 + 1/2*(2*T + strip_height)],[x1-R/3, y1 + 1/2*(2*T + strip_height)],[x1, y1 + 1/2*(2*T + strip_height)],[x1 + R/2, y1 + 3/5*T], [startU[0] - A + t + Sr[0],y1 + 1/5*T],[startU[0] - A + t + Sr[0],startU[1]-Sr[1]]]
    # if Sr[0] > 0:
        # points[-2][0] = startU[0] - A + t + Sr[0] + w_end
        # points[-1][0] = startU[0] - A + t + Sr[0] + w_end
    
    #enclose contact pad with etched line (right edge)
    patch_start = [points[0][0]-(width_guide[0]-width_feedline[0])/2,points[0][1]]
    patch_end = [points[3][0]+1,points[3][1]]
    
    #enclosing contact pad/patch to define it/ this big_patch + guide will be etched away
    big_patch = gdspy.FlexPath([patch_start,points[1]], width_guide[0],corners='circular bend',bend_radius = bend_radii[1])
    big_patch = big_patch.segment(points[2],width_guide[1]).segment(patch_end,width_guide[2])
    guide = gdspy.FlexPath([points[3],points[4]],width_guide[2],corners='circular bend',bend_radius=bend_radii[2])
    guide = guide.segment(points[5],width_guide[3]).segment(points[6],width_guide[4]).segment(points[7],width_guide[5])
    
    #draw the feedline (not etched away)
    patch = gdspy.FlexPath([points[0],points[1]], width_feedline[0],corners='circular bend',bend_radius = bend_radii[1])
    patch = patch.segment(points[2],width_feedline[1]).segment(patch_end,width_feedline[2])
    feedline = gdspy.FlexPath([points[3],points[4]],width_feedline[2],corners='circular bend',bend_radius=bend_radii[2])
    feedline = feedline.segment(points[5],width_feedline[2]).segment(points[6],width_feedline[4]).segment(points[7],width_feedline[5])
    
    pad = gdspy.boolean(big_patch,patch, 'not')
    line = gdspy.boolean(guide,feedline, 'not')
    
    #coordinates for left feedline + guide (FlexPath)
    points2 = [[x2+E,y2 - 1/2*(2*T + strip_height)],[x2+3*E/4, y2 - 1/2*(2*T + strip_height)],[x2+E/4, y2 - 1/2*(2*T + strip_height)],[x2+R/3, y2 - 1/2*(2*T + strip_height)],[x2, y2 - 1/2*(2*T + strip_height)],[x2 - R/2, y2 - 3/5*T],[stopU[0]+ Sr[0],y2 - 1/5*T],[stopU[0]+ Sr[0],stopU[1]+Sr[1]]]
    # if Sr[0] > 0:
        # points2[-2][0] = stopU[0]+ Sr[0] + w_end
        # points2[-1][0] = stopU[0]+ Sr[0] + w_end
    #enclose contact pad with etched line (right edge)
    patch_start2 = [points2[0][0]+(width_guide[0]-width_feedline[0])/2,points2[0][1]]
    patch_end2 = [points2[3][0]-1,points2[3][1]]
    
    big_patch2 = gdspy.FlexPath([patch_start2,points2[1]], width_guide[0],corners='circular bend',bend_radius = bend_radii[1])
    big_patch2 = big_patch2.segment(points2[2],width_guide[1]).segment(patch_end2,width_guide[2])
    guide2 = gdspy.FlexPath([points2[3],points2[4]],width_guide[2],corners='circular bend',bend_radius=bend_radii[2])
    guide2 = guide2.segment(points2[5],width_guide[3]).segment(points2[6],width_guide[4]).segment(points2[7],width_guide[5])
   
    
    patch2 = gdspy.FlexPath([points2[0],points2[1]], width_feedline[0],corners='circular bend',bend_radius = bend_radii[1])
    patch2 = patch2.segment(points2[2],width_feedline[1]).segment(patch_end2,width_feedline[2])    
    feedline2 = gdspy.FlexPath([points2[3],points2[4]],width_feedline[2],corners='circular bend',bend_radius=bend_radii[2])
    feedline2 = feedline2.segment(points2[5],width_feedline[3]).segment(points2[6],width_feedline[4]).segment(points2[7],width_feedline[5])

    pad2 = gdspy.boolean(big_patch2, patch2, 'not')
    line2 = gdspy.boolean(guide2,feedline2, 'not')
    
    window_left_guide = gdspy.boolean(window, line, 'or')
    window_both_guides = gdspy.boolean(window_left_guide, line2, 'or')
    window_left_feedline = gdspy.boolean(window_both_guides,feedline ,'not')
    window_both_feedlines = gdspy.boolean(window_left_feedline,feedline2,'not')
    pads = gdspy.boolean(pad,pad2,'or')
    
    return window_both_feedlines,pads

def waveguide_extended_negative_new(Q, unitcell_size,startM, startU, stopU, strip_height, tw, tv, N_ghost, t = 2e-6, Sr = [1e-6,2e-6], Sg = 60e-6,T=100e-6,R=200e-6,f=24e-6, A = 50e-6):
    tw = tw*1e6
    tv = tv*1e6
    t = t*1e6
    Sr = Sr[0]*1e6,Sr[1]*1e6 #spacing feedline to resonator [x,y]-direction
    Sg = Sg*1e6 #lateral spacing to ground planes
    T = T*1e6 #thickness of ground planes
    R = R*1e6 #width of ground planes
    # E = 800
    chip_length = 6200
    f = f*1e6
    A = A*1e6
    w_patch = 360.0
    w_middle = 90.0
    w_start = 10.0 #start width feedline
    w_end = 2.0 #end width feedline
    
    
    #calculate length ghosts demand
    if (N_ghost-1)%2==0:
        extent_ghosts = tw + tv + N_ghost*A + (N_ghost-1)/2*(tw + tv)
        print('extent ghosts is: ',extent_ghosts, ' (N_ghost odd)')
    else:
        extent_ghosts = tw + tv + N_ghost*A + (N_ghost-1)/2*tw + ((N_ghost-1)/2 - 1)*tv
        print('extent ghosts is: ',extent_ghosts, ' (N_ghost even)')
    
    #(x1,y1) and (x2,y2): corners of outer box, (x3,y3) and (x4,y4): corners of window
    x1, y1 = startM[0] - A/2 - Sg - R - extent_ghosts, startM[1] - f - T
    x2, y2 = x1 + 2*(R + Sg) + Q*unitcell_size[0]-tw + 2*extent_ghosts, y1 + 2*T + strip_height
    x3, y3 = startM[0] - A/2 - Sg - extent_ghosts, startM[1] - f
    x4, y4 = x3 + 2 * Sg + Q*unitcell_size[0]-tw + 2*extent_ghosts, y3 + strip_height
    
    width_feedline = [w_patch,w_middle,w_start, w_start - 1/4*(w_start - w_end), w_end + 1/4*(w_start-w_end), w_end]
    width_guide = [11/5*w for w in width_feedline]
    width_guide[0] = 73/45*width_feedline[0]
    width_guide[1] = 31/30*width_feedline[1]
    bend_radii = [3*w for w in width_guide]
    dyke = [(width_guide[0]-width_feedline[0])/2,(width_guide[1]-width_feedline[1])/2,(width_guide[2]-width_feedline[2])/2]
    
    E = 0.5*(chip_length - (x2-x1)) - (width_guide[0]-width_feedline[0])/2
    
    #draw outer box + window
    big_plane = gdspy.Rectangle([x1,y1], [x2,y2])
    window = gdspy.Rectangle([x3,y3], [x4,y4])
    
    #coordinates/ dimensions for FlexPath: right feedline + guide (etched away)
    Y0 = y1 + 1/2*(2*T + strip_height)
    points = [[x1-E,Y0],[x1-7*E/8, Y0],[x1-3*E/4, Y0],[x1-E/8, Y0],[x1, Y0],[x1 + 3*R/4, y1 + 3/5*T], [startU[0] - A + t + Sr[0],y1 + 1/5*T],[startU[0] - A + t + Sr[0],startU[1]-Sr[1]]]
    
    if Sr[0] > 0:
        points[-2][0] = startU[0] - A + t + Sr[0] + w_end
        points[-1][0] = startU[0] - A + t + Sr[0] + w_end
    
    #enclose contact pad with etched line (right edge)
    patch_start = [points[0][0]-dyke[0],points[0][1]]
    patch_end = [points[2][0]+1,points[2][1]]
    
    #enclosing contact pad/patch to define it/ this big_patch + guide will be etched away
    big_patch = gdspy.FlexPath([patch_start,points[1]], width_guide[0])#,corners=#'circular bend',bend_radius = bend_radii[1])
    big_patch = big_patch.segment(patch_end,width_guide[1])
    guide = gdspy.FlexPath([points[2],points[3]],width_guide[1],corners='circular bend',bend_radius=bend_radii[2])
    guide = guide.segment(points[4], width_guide[2]).segment(points[5],width_feedline[3]).segment(points[6],width_guide[4]).segment(points[7],width_guide[5])
    
    #draw the feedline (not etched away)
    patch = gdspy.FlexPath([points[0],points[1]], width_feedline[0])#,corners='circular bend',bend_radius = bend_radii[1])
    patch = patch.segment(patch_end,width_feedline[1])
    feedline = gdspy.FlexPath([points[2],points[3]],width_feedline[1],corners='circular bend',bend_radius=bend_radii[2])
    feedline = feedline.segment(points[4], width_feedline[2]).segment(points[5],width_feedline[3]).segment(points[6],width_feedline[4]).segment(points[7],width_feedline[5])
    
    pad = gdspy.boolean(big_patch,patch, 'not')
    line = gdspy.boolean(guide,feedline, 'not')
    
    #coordinates for left feedline + guide (FlexPath)
    Y02 = y2 - 1/2*(2*T + strip_height)
    points2 = [[x2+E,Y02],[x2+7*E/8, Y02],[x2+3*E/4, Y02],[x2+E/8, Y02],[x2, Y02],[x2 - 3*R/4, y2 - 3/5*T],[stopU[0] - Sr[0],y2 - 1/5*T],[stopU[0] - Sr[0],stopU[1]+Sr[1]]]
    if Sr[0] > 0:
        points2[-2][0] = stopU[0]+ Sr[0] + w_end
        points2[-1][0] = stopU[0]+ Sr[0] + w_end
    #enclose contact pad with etched line (right edge)
    patch_start2 = [points2[0][0]+dyke[0],points2[0][1]]
    patch_end2 = [points2[2][0]-1,points2[2][1]]
    
    big_patch2 = gdspy.FlexPath([patch_start2,points2[1]], width_guide[0])#,corners=#'circular bend',bend_radius = bend_radii[1])
    big_patch2 = big_patch2.segment(patch_end2,width_guide[1])
    guide2 = gdspy.FlexPath([points2[2],points2[3]],width_guide[1],corners='circular bend',bend_radius=bend_radii[2])
    guide2 = guide2.segment(points2[4], width_guide[2]).segment(points2[5],width_guide[3]).segment(points2[6],width_guide[4]).segment(points2[7],width_guide[5])
   
    
    patch2 = gdspy.FlexPath([points2[0],points2[1]], width_feedline[0])#,corners='circular bend',bend_radius = bend_radii[1])
    patch2 = patch2.segment(patch_end2,width_feedline[1])   
    feedline2 = gdspy.FlexPath([points2[2],points2[3]],width_feedline[1],corners='circular bend',bend_radius=bend_radii[2])
    feedline2 = feedline2.segment(points2[4],width_feedline[2]).segment(points2[5],width_feedline[3]).segment(points2[6],width_feedline[4]).segment(points2[7],width_feedline[5])

    pad2 = gdspy.boolean(big_patch2, patch2, 'not')
    line2 = gdspy.boolean(guide2,feedline2, 'not')
    
    window_left_guide = gdspy.boolean(window, line, 'or')
    window_both_guides = gdspy.boolean(window_left_guide, line2, 'or')
    window_left_feedline = gdspy.boolean(window_both_guides,feedline ,'not')
    window_both_feedlines = gdspy.boolean(window_left_feedline,feedline2,'not')
    pads = gdspy.boolean(pad,pad2,'or')
    
    return window_both_feedlines,pads

def T_feedline_extended_negative(Q, unitcell_size,startM, startU, stopU, strip_height, tw, tv, N_ghost, t = 2e-6, Sr = [1e-6,2e-6], Sg = 60e-6,T=100e-6,D=200e-6,R=200e-6,f=24e-6, A = 50e-6,B=60e-6):
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
    w_patch = 360.0
    w_middle = 90.0
    w_start = 10.0 #start width feedline
    
    
    #(x1,y1) and (x2,y2): corners of outer box, (x3,y3) and (x4,y4): corners of window
    x3, y3 = startM[0] - A/2 - Sr[0] - 2*w_start, startM[1] - f
    x4, y4 = x3 + 2 * Sr[0] + 4*w_start + Q*unitcell_size[0]-tw, y3 + strip_height
    x1, y1 = x3 - R, startM[1] - f - T
    x2, y2 = x4 + R, y1 + 2*T + strip_height
    # x1, y1 = startM[0] - A/2 - Sg - D , startM[1] - f - T
    # x2, y2 = x1 + 2*(D + Sg) + Q*unitcell_size[0]-tw, y1 + 2*T + strip_height
    # x3, y3 = startM[0] - A/2 - Sg, startM[1] - f
    # x4, y4 = x3 + 2 * Sg + Q*unitcell_size[0]-tw, y3 + strip_height
    
    
    width_feedline = [w_patch,w_middle,w_start]
    width_guide = [73/45*w for w in width_feedline]
    width_guide[1] = 47/45*width_feedline[1]
    width_guide[2] = 11/5*width_feedline[2]
    bend_radii = [3*w for w in width_guide]
    dyke = [(width_guide[0]-width_feedline[0])/2,(width_guide[1]-width_feedline[1])/2,(width_guide[2]-width_feedline[2])/2]
    
    
    E = 0.5*(chip_length - (x2-x1)) - (width_guide[0]-width_feedline[0])/2
    
    #draw outer box + window
    # big_plane = gdspy.Rectangle([x1,y1], [x2,y2])
    window = gdspy.Rectangle([x3,y3], [x4,y4])
    
    #coordinates/ dimensions for FlexPath: right feedline + guide (etched away)
    Y0 = y1 + 1/2*(2*T + strip_height)
    points = [[x1-E,Y0],[x1-7*E/8, Y0],[x1-3*E/4, Y0],[x1,Y0],[x3 - D, Y0],[startU[0] + t/2 - A - Sr[0] - w_start/2, Y0]]
    
    #enclose contact pad with etched line (right edge)
    patch_start = [points[0][0]- dyke[0],points[0][1]]
    patch_end = [points[2][0]+1,points[2][1]]
    
    #enclosing contact pad/patch to define it/ this big_patch + guide will be etched away
    big_patch = gdspy.FlexPath([patch_start,points[1]], width_guide[0])#,corners=#'circular bend',bend_radius = bend_radii[1])
    big_patch = big_patch.segment(patch_end,width_guide[1])
    guide = gdspy.FlexPath([points[2],points[4]],width_guide[1]).segment(points[5],width_guide[2])
    T_bar_guide = gdspy.Rectangle([points[5][0] - w_start/2 - dyke[2], startU[1] + B + dyke[2]], [points[5][0] + w_start/2 + dyke[2], startU[1] - dyke[2]])
    
    guide = gdspy.boolean(guide, T_bar_guide, 'or')
    
    #draw the feedline (not etched away)
    patch = gdspy.FlexPath([points[0],points[1]], width_feedline[0])#,corners='circular bend',bend_radius = bend_radii[1])
    patch = patch.segment(patch_end,width_feedline[1])
    feedline = gdspy.FlexPath([points[2],points[4]],width_feedline[1]).segment(points[5],width_feedline[2])#,corners='circular bend',bend_radius=bend_radii[2])
    T_bar = gdspy.Rectangle([points[5][0] - w_start/2, startU[1] + B], [points[5][0] + w_start/2, startU[1]])
    
    feedline = gdspy.boolean(feedline, T_bar, 'or')
    
    pad = gdspy.boolean(big_patch,patch, 'not')
    line = gdspy.boolean(guide,feedline, 'not')
    
    #coordinates for left feedline + guide (FlexPath)
    Y02 = y2 - 1/2*(2*T + strip_height)
    points2 = [[x2+E,Y02],[x2+7*E/8, Y02],[x2+3*E/4, Y02],[x2,Y02],[x4 + D, Y02],[stopU[0] + t/2 + Sr[0] + w_start/2,Y02]]
    
    #enclose contact pad with etched line (right edge)
    patch_start2 = [points2[0][0]+dyke[0],points2[0][1]]
    patch_end2 = [points2[2][0]-1,points2[2][1]]
    
    big_patch2 = gdspy.FlexPath([patch_start2,points2[1]], width_guide[0])#,corners='circular bend',bend_radius = bend_radii[1])
    big_patch2 = big_patch2.segment(patch_end2,width_guide[1])
    guide2 = gdspy.FlexPath([points2[2],points2[4]],width_guide[1]).segment(points2[5],width_guide[2])#,corners='circular bend',bend_radius=bend_radii[2])
    T_bar_guide2 = gdspy.Rectangle([points2[5][0] - w_start/2-dyke[2], stopU[1]+dyke[2]], [points2[5][0] - w_start/2 - dyke[2], stopU[1]-B-dyke[2]])
    
    guide2 = gdspy.boolean(guide2, T_bar_guide2, 'or')
    
    patch2 = gdspy.FlexPath([points2[0],points2[1]], width_feedline[0])#,corners='circular bend',bend_radius = bend_radii[1])
    patch2 = patch2.segment(patch_end2,width_feedline[1])   
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
    
    return window_both_feedlines,pads

"""
~~~~~~~~~~~~ Testing functions ~~~~~~~~
"""


# if __name__ == '__main__':
#     L = 279e-6              #total length of inductor in SI units
#     s = 4.5e-6              #interspacing of inductor in SI units
#     w = 0.5e-6              #width of inductor wire
#     A = 50e-6               #horizontal dimension of capacitor
#     t = 2e-6                #thickness of capacitor plates
#     tv = 12e-6              #intracell spacing
#     tw = 24e-6               #intercell spacing
#     k = 1/3                 #fraction determining how thick the ground strip between resonators is : k*min(tv,tw)
#     e = 20.5e-6             #horizontal dimension of ground patches
#     f = 24e-6               #vertical dimension of ground patches
#     r = 29e-6               #vertical spacing of resonator "head" to ground plane
#     Q = 4                   #number of unit cells
#     N_ghost = 2             #number of ghosts: For no ghosts, put zero
#     Sf2r = [0.0,2e-6]      #spacing to feedline in [x,y]-direction
#     ground_yn = True        #ground in between yes (True) or no (False) 
    
#     unitcell_size = [2*A*1e6 + tv*1e6 + tw*1e6,0]
    
#     test = gdspy.Cell('negative')
    
#     unitcell, ground, startM,startU, stopUy, strip_height,center1st, center2nd = unit_cell(L, s, w, A, t, tv,tw,e,f,r,gamma = k, ground_in_between = ground_yn)
    
#     stopU = [startU[0]-A*1e6 + Q*unitcell_size[0]-tw*1e6, stopUy]
#     centerQth = [center1st[0]-A*1e6 + Q*unitcell_size[0] - tw*1e6, center2nd[1]]
    
#     blib = only_waveguide(startM, startU, strip_height)
#     # blub, blab = waveguide_extended_negative(Q, unitcell_size,startM, startU, stopU, strip_height, tw, tv, N_ghost, Sr = Sf2r, f=24e-6, A = 50e-6)

#     # test.add(blub)
#     # test.add(blab)
#     test.add(blib)
#     lib = gdspy.GdsLibrary()
#     lib.add(test)
#     gdspy.LayoutViewer()
#     lib.write_gds("Only waveguide.gds")