# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 15:12:30 2021

@author: jouanny

Bunch of functions for the topological waveguide/surface project
"""

import numpy as np
import gdspy

"""
--- Twirled spiral: old
""" 

def spiral(L, W, I, D, Centre, SC = 0):
    """
    Generates a list of points to generate multiple Rectangles to form a squared spiral.
    L: Length of the wirestrip \n
    W: Width of the wirestrip \n
    I: Interspacing between wires \n
    D: Diameter
    Centre: List like variable with x and y positions of the centre
    
    """
    #For length, once every two time, count the length in x or y direction
    #Start from the centre
    #Finally once resonator is ready draw coplanar waveguide
    L = L*1e6
    D = D*1e6
    I = I*1e6
    W = W*1e6
    D += I 
    # SET A MINIMAL OF 4 TURNs !!! with the ratio L/D
    points = []
    xc, yc = Centre #Estimating the centre
    x0, y0 = xc - D/2, yc - D/2
    points.append([x0,y0])
    LTemp = 0
    k = 0
    
    xt, yt = points[0]
    x, y = xt + W, yt + D
    points.append([x, y])
    LTemp += points[1][1] - points[0][1]
    
    #Second Rectangle
    x, y = xt + D - I, yt + D - W
    points.append([x, y])
    LTemp += points[2][0] - points[1][0]
    
    #Third Rectangle
    x, y = xt + D - I - W, yt + I
    points.append([x, y])
    LTemp += points[2][1] - points[3][1]
        
    while LTemp < L:
        #put a if with a condition modulo 4 (4k+1)
        #Use break to stop the loop
        #Everything should be ready.
        #increment after if statements
        
        #We build first a basis to build the resonator upon to
        #First Rectangle
        
        
        #Now we can build an "automatic generation"
        
        k += 1
        for i in range(4):
            print('i',i)
            if i == 0: #LEFT
                xt, yt = points[i]
                print(points[i])
                x, y = xt + k*(W + I), yt + k*(W + I)
                if LTemp + (points[-1][0] - x) > L:
                     Lf = L - LTemp
                     xt, yt = points[-1]
                     x, y = xt - Lf, yt + W
                     points.append([x, y])
                     LTemp += Lf
                     break
                else:
                    LTemp += points[-2][0] - points[-1][0]
                    points.append([x ,y])
                #print('FullList',k,points)
                
                print(k+i)
                    
            elif i == 1:#TOP
                xt, yt = points[i]
                x, y = xt + k*(I + W), yt - k*(W + I)
                if LTemp + (y - points[-1][1]) > L:
                    Lf = L - LTemp
                    xt, yt = points[-1]
                    x, y = xt + W, yt + Lf
                    points.append([x, y])
                    LTemp += Lf
                    break
                else:
                    points.append([x ,y])
                    LTemp += points[-1][1] - points[-2][1]
                
            elif i == 2:#RIGHT
                xt, yt = points[i]
                x, y = xt - k*(W + I), yt - k*(I + W)
                if LTemp + (x - points[-1][0]) > L:
                    xt, yt = points[-1]
                    Lf = L - LTemp#Final length
                    x, y = xt + Lf, yt - W
                    points.append([x, y])
                    LTemp += Lf
                    break
                else:
                    points.append([x ,y])
                    LTemp += points[-1][0] - points[-2][0]
                
            else:#BOTTOM
                xt, yt = points[i]
                x, y = xt - k*(I + W), yt + k*(W + I)
                if LTemp + (points[-1][1] - y) > L:
                    xt, yt = points[-1]
                    Lf = L - LTemp#Final length
                    x, y = xt - W, yt - Lf
                    points.append([x, y])
                    LTemp += Lf
                    break
                else:
                    points.append([x ,y])
                    LTemp += points[-2][1] - points[-1][1]
                print(k)
    points[0][1] = points[0][1] + I
            
    return points, LTemp

def LeftWG(Origin, CapaHeight, S,SPad = 150,LProg = 160,LThin = 375,t = 10):
    """
    Program that generates a left coplanar waveguide to couple array of resonators.
    Origin: Center point chosen, I think it would be wise to take the right of the waveguide to deal with the spacing between resonators and waveguide easily.
    CapaHeight: Height of the coupling plate
    S: Spacing between the ground plane and the waveguide.
    """
    
    
    
    CapaHeight = CapaHeight*1e6
    S = S*1e6
    points = []
    
    # SPad = 150
    # LProg = 160
    # LThin = 375
    # t = 10
    
    points.append([-Origin[0]*1e6,-Origin[1]*1e6])
    
    x1, y1 = points[-1][0], points[-1][1] + SPad
    points.append([x1, y1])
    
    x2, y2 = points[-1][0]+SPad, points[-1][1]
    points.append([x2, y2])
    
    x3, y3 = points[-1][0] + LProg, points[-1][1] - SPad/2 + t/2
    points.append([x3, y3])
    
    x4, y4 = points[-1][0] + LThin, points[-1][1]
    points.append([x4, y4])
    
    x5, y5 = points[-1][0], points[-1][1] + (CapaHeight-t)/2
    points.append([x5, y5])
    
    x6, y6 = points[-1][0] + t, points[-1][1]
    points.append([x6, y6])
    
    x7, y7 = points[-1][0], points[-1][1]-CapaHeight
    points.append([x7, y7])
    
    x8, y8 = points[-1][0]-t, points[-1][1]
    points.append([x8, y8])
    
    x9, y9 = points[-1][0], points[-1][1] + (CapaHeight-t)/2
    points.append([x9, y9])
    
    x10, y10 = points[-1][0] - LThin, points[-1][1]
    points.append([x10, y10])
    
    x11, y11 = points[-1][0] - LProg, points[-1][1] - SPad/2 + t/2
    points.append([x11, y11])
    
    x12, y12 = points[0]
    points.append([x12, y12])
    
    #S = 7
    Hprog = 175
    
    points2 = []
    
    a0 ,b0 = x3, y3 + S
    points2.append([a0, b0])
    
    a1, b1 = points2[-1][0] - LProg, points2[-1][1] + Hprog
    points2.append([a1, b1])
    
    a2, b2 = points2[-1][0] - 215, points2[-1][1]
    points2.append([a2, b2])
    
    a3, b3 = points2[-1][0], points2[-1][1] + 300
    points2.append([a3, b3])
    
    a4, b4 = points2[-1][0]+630, points2[-1][1]
    points2.append([a4, b4])
    #print(a4,b4,points2)
    a5, b5 = points2[-1][0], points2[0][1]
    points2.append([a5, b5])
    
    a6, b6 = points2[0][0], points2[0][1]
    points2.append([a6, b6])
    
    
    points3 = []
    
    c0 ,d0 = x10, y10 - S
    points3.append([c0, d0])
    
    c1, d1 = points3[-1][0] - LProg, points3[-1][1] - Hprog
    points3.append([c1, d1])
    
    c2, d2 = points3[-1][0] - 215, points3[-1][1]
    points3.append([c2, d2])
    
    c3, d3 = points3[-1][0], points3[-1][1] - 300
    points3.append([c3, d3])
    
    c4, d4 = points3[-1][0]+630, points3[-1][1]
    points3.append([c4, d4])
    #print(a4,b4,points2)
    c5, d5 = points3[-1][0], points3[0][1]
    points3.append([c5, d5])
    
    c6, d6 = points3[0][0], points3[0][1]
    points3.append([c6, d6])
    
    return np.asarray(points), np.asarray(points2), np.asarray(points3)

def CPWPGround(Centre, H, Cl, Cd, L, t, Sx, Sy, Hg, Lg, lg):
    Centre = Centre[0]*1e6, Centre[0]*1e6
    H = H*1e6
    Cl = Cl*1e6
    Cd = Cd*1e6
    L = L*1e6
    t = t*1e6
    
    Poly = []
    Poly.append([-Centre[0],-Centre[1]])#
    
    x1, y1 = Poly[-1][0], Poly[-1][1] + H
    Poly.append([x1,y1])#
    
    x2, y2 = Poly[-1][0] + Cl, Poly[-1][1]
    Poly.append([x2,y2])#
    
    x3, y3 = Poly[-1][0] + Cd, Poly[0][1] + H/2 + t/2
    Poly.append([x3,y3])#
    
    x4, y4 = Poly[-1][0] + L, Poly[-1][1]
    Poly.append([x4,y4])#
    
    x5, y5 = Poly[-1][0] + Cd, y2
    Poly.append([x5,y5])#
    
    x6, y6 = Poly[-1][0] + Cl, Poly[-1][1]
    Poly.append([x6,y6])#
    
    x7, y7 = Poly[-1][0], Poly[-1][1] - H
    Poly.append([x7,y7])#
    
    x8, y8 = Poly[-1][0] - Cl, Poly[0][1]
    Poly.append([x8,y8])#
    
    x9, y9 = Poly[-1][0] - Cd, Poly[0][1] + H/2 - t/2
    Poly.append([x9,y9])#
    
    x10, y10 = Poly[-1][0] - L, Poly[0][1] + H/2 - t/2
    Poly.append([x10,y10])#
    
    x11, y11 = Poly[0][0] + Cl, Poly[0][1]
    Poly.append([x11,y11])
    
    x12, y12 = Poly[0][0], Poly[0][1]
    Poly.append([x12,y12])
    
    
    
    #Ground
    Ground = []
    a0, b0 = Poly[0][0] - Sx, Poly[0][1] - Sy
    Ground.append([a0, b0])#
    
    a1, b1 = Poly[1][0] - Sx, Poly[1][1] + Sy
    Ground.append([a1, b1])#
    
    a2, b2 = Ground[-1][0] + Sx, Ground[-1][1]
    Ground.append([a2, b2])#
    
    a3, b3 = Poly[3][0] + Sx/2, Poly[3][1] + Sy/2
    Ground.append([a3 ,b3])#
    
    a4, b4 = Poly[4][0] - Sx/2, Poly[4][1] + Sy/2
    Ground.append([a4 ,b4])#
    
    a5, b5 = Poly[5][0] + Sx/2, Poly[5][1] + Sy/2
    Ground.append([a5, b5])
    
    a6, b6 = Poly[6][0] + Sx, Poly[6][1] + Sy
    Ground.append([a6, b6])
    
    a7, b7 = Poly[7][0] - Sx, Poly[7][1] - Sy
    Ground.append([a7, b7])
    
    a8, b8 = Ground[-1][0] - Sx, Ground[-1][1]
    Ground.append([a8, b8])
    
    a9, b9 = Poly[9][0] - Sx/2, Poly[9][1] - Sy/2
    Ground.append([a9, b9])
    
    a10, b10 = Ground[-1][0] - lg/2, Ground[-1][1]
    Ground.append([a10, b10])
    
    a11, b11 = Ground[-1][0], Ground[-1][1] - 460
    Ground.append([a11, b11])
    
    a12, b12 = Ground[-1][0] + lg, Ground[-1][1]
    Ground.append([a12, b12])

    a13, b13 = Ground[-1][0], Ground[-1][1] + Hg
    Ground.append([a13, b13])

    a14, b14 = Ground[-1][0] - Lg, Ground[-1][1]
    Ground.append([a14, b14])

    a15, b15 = Ground[-1][0], Ground[-1][1] - Hg
    Ground.append([a15, b15])
    
    a16, b16 = Ground[-1][0] + lg, Ground[-1][1]
    Ground.append([a16, b16])
    
    a17, b17 = Ground[-1][0], Ground[-1][1] + 460
    Ground.append([a17, b17])
    
    a18, b18 = Ground[-1][0] - lg/2, Ground[-1][1]
    Ground.append([a18, b18])
    
    a19,b19 = Ground[0][0], Ground[0][1]
    Ground.append([a19, b19])
    
    return Poly, Ground

def TwirledSpiral(Centre, L, W, I,Offset = (0,0)):
    """
    Parameters:
        -Centre: Array of two points
        -L: Length of the resonator
        -W: Width of the wire of the resonator
        -I: Interspacing.
    """
    Offset = Offset[0]*1e6, Offset[1]*1e6
    Centre = (Centre[0]*1e6, Centre[1]*1e6)
    L = L*1e6
    W = W*1e6
    I = I*1e6
    LTot = 0
    
    #Code to modify the interspacing to not get edges that are non complete.
    k = 1
    j = 0
    Lt = 0
    Lt += (I+W)/2
    while Lt < L/2:
        for i in range(4):
            if i == 0:
                if Lt + k*(I+W) < L/2:
                    Lt += k*(I+W)
                else:
                    r = L/2 - Lt
                    Lt += r
                    break
                k += 1
            elif i == 1:
                if Lt + k*(I+W) - W/2 < L/2:
                    Lt += k*(I+W)-W/2
                else:
                    r = L/2 - Lt
                    Lt += r
                    break
                k += 1
            elif i == 2:
                if Lt + k*(I+W) < L/2:
                    Lt += k*(I+W)
                else:
                    r = L/2 - Lt
                    Lt += r
                    break 
                k += 1
            else:
                if Lt + k*(I+W) - W/2 < L/2:
                    Lt += k*(I+W) - W/2
                else:
                    r = L/2 - Lt
                    Lt += r
                    break
                k += 1 
        j += 4
    
    # print(k)
    
    # print(L)
    # print(r)
    # Ln = L/2-r
    # print(Ln)
    
    # In = ((L/2-r)/k)-(3/4)*W - (I+W)/2
    # print("In", In)
    
    #It doesn't work.
    
    points = []
    x0, y0 = Centre #LOOK AT HOW WE CAN ESTIMATE LATER.
    points.append([x0,y0])
    
    x1, y1 = x0 + I + W, y0
    points.append([x1, y1])
    i = 0
    k=1
    j=0
    LTemp = 0
    LTemp += (I+W)/2
    while LTemp < L/2:
        for i in range(4):
            if i == 0:
                if LTemp + k*(I+W) < L/2:
                    points.append([points[j+1][0], points[j+1][1]+k*(I+W)])
                    LTemp += k*(I+W)
                    k += 1
                else:
                    l = L/2 - LTemp
                    points.append([points[j+1][0], points[j+1][1]+l])
                    LTemp += l
                    k += 1
                    break
            
            elif i == 1:
                if LTemp + k*(I+W) - W/2 < L/2:
                    points.append([points[j+2][0]-k*(I+W)+W/2, points[j+2][1]])
                    LTemp += k*(I+W)-W/2
                    k += 1
                else:
                    l = L/2 - LTemp
                    points.append([points[j+2][0]-l, points[j+2][1]])
                    LTemp += l
                    k += 1
                    break
                
            elif i == 2:
                if LTemp + k*(I+W) < L/2:
                    points.append([points[j+3][0], points[j+3][1] - k*(I+W)])
                    LTemp += k*(I+W)
                    k += 1
                else:
                    l = L/2 - LTemp
                    points.append([points[j+3][0], points[j+3][1]-l])
                    LTemp += l
                    k += 1
                    break 
                
            else:
                if LTemp + k*(I+W) - W/2 < L/2:
                    points.append([points[j+4][0] + k*(I+W) - W/2, points[j+4][1]])
                    LTemp += k*(I+W) - W/2
                    k += 1
                else:
                    l = L/2 - LTemp
                    points.append([points[j+4][0] + l, points[j+4][1]])
                    LTemp += l
                    k += 1
                    #print(j)
                    break
                
        j += 4
        #2Spiral
    points2 = []
    points2.append([x0+W/2,y0-W/2])
    
    a1, b1 = points2[0][0], points2[0][1] - I - W/2
    points2.append([a1, b1])

    LTot += LTemp
    i = 0
    k=2
    j=0
    LTemp = 0
    LTemp += (I+W)/2 + I + W/2
    while LTemp < L/2:
        for i in range(4):
            if i == 0:
                if LTemp + k*(I+W) - W/2 < L/2:
                    points2.append([points2[j+1][0] + k*(I+W) - W/2, points2[j+1][1]])
                    LTemp += k*(I+W) - W/2
                    k += 1
                else:
                    l = L/2 - LTemp
                    points2.append([points2[j+1][0]+l, points2[j+1][1]])
                    LTemp += l
                    k += 1
                    break
            
            elif i == 1:
                if LTemp + k*(I+W) < L/2:
                    points2.append([points2[j+2][0], points2[j+2][1]+k*(I+W)])
                    LTemp += k*(I+W)
                    k += 1
                else:
                    l = L/2 - LTemp
                    points2.append([points2[j+2][0], points2[j+2][1]+l])
                    LTemp += l
                    k += 1
                    break
                
            elif i == 2:
                if LTemp + k*(I+W) - W/2 < L/2:
                    points2.append([points2[j+3][0] - k*(I+W) + W/2, points2[j+3][1]])
                    LTemp += k*(I+W)-W/2
                    k += 1
                else:
                    l = L/2 - LTemp
                    points2.append([points2[j+3][0] - l, points2[j+3][1]])
                    LTemp += l
                    k += 1
                    break 
                
            else:
                if LTemp - k*(I+W) < L/2:
                    points2.append([points2[j+4][0], points2[j+4][1]-k*(I+W)])
                    LTemp += k*(I+W)
                    k += 1
                else:
                    l = L/2 - LTemp
                    points2.append([points2[j+4][0], points2[j+4][1]-l])
                    LTemp += l
                    k += 1
                    #print(j)
                    break      
        j += 4
    LTot += LTemp
    
    points = np.asarray(points)
    points = points #+ ((points[1][0]-points[0][0])/2)
    points2 = np.asarray(points2)
    points2 = points2 #+ ((points[1][0]-points[0][0])/2)

    
    # EdgesRight = points[-1][0] - points[-2][0], points[-1][1] - points[-2][1]
    # EdgesLeft = points2[-1][0] - points2[-2][0], points2[-1][1] - points2[-2][1]
    
    return points - Offset, points2 - Offset, LTot#, EdgesRight, EdgesLeft


def Cpw(El, Er, Hr, Sr, w = 10e-6, t = 10e-6, Lc = 215e-6, Lline = 470e-6, T = 50e-6, s = 7.1e-6):
    El = El[0]*1e6 ,El[1]*1e6
    Er = Er[0]*1e6, Er[1]*1e6
    Hr = Hr*1e6
    Sr = Sr*1e6
    w = w*1e6
    t = t*1e6
    Lc = Lc*1e6
    Lline = Lline*1e6
    T = T*1e6
    s = s*1e6
    
    #Left CPW
    Lcpw = []
    
    x1, y1 = El[0] - Lline - t - Sr, El[1] + w/2
    Lcpw.append([x1, y1])
    
    x2, y2 = x1 + Lline, y1
    Lcpw.append([x2, y2])
    
    x3, y3 = x2, y2 + (Hr - w)/2
    Lcpw.append([x3, y3])
    
    x4, y4 = x3 + t, y3
    Lcpw.append([x4, y4])
    
    x5, y5 = x4, y4 - Hr
    Lcpw.append([x5, y5])
    
    x6, y6 = x5 - t, y5
    Lcpw.append([x6, y6])
    
    x7, y7 = x6, y6 + (Hr - w)/2
    Lcpw.append([x7, y7])
    
    x8, y8 = x7 - Lline, y7
    Lcpw.append([x8, y8])
    
    x9, y9 = x1, y1
    Lcpw.append([x9, y9])
    
    
    #Right CPW
    Rcpw = []
    
    a1, b1 = Er[0] + Lline + t + Sr, Er[1] + w/2
    Rcpw.append([a1, b1])
    
    a2, b2 = a1 - Lline, b1
    Rcpw.append([a2, b2])
    
    a3, b3 = a2, b2 + (Hr - w)/2
    Rcpw.append([a3, b3])
    
    a4, b4 = a3 - t, b3
    Rcpw.append([a4, b4])
    
    a5, b5 = a4, b4 - Hr
    Rcpw.append([a5, b5])
    
    a6, b6 = a5 + t, b5
    Rcpw.append([a6, b6])
    
    a7, b7 = a6, b6 + (Hr - w)/2
    Rcpw.append([a7, b7])
    
    a8, b8 = a7 + Lline, b7
    Rcpw.append([a8, b8])
    
    a9, b9 = a1, b1
    Rcpw.append([a9, b9])
    
    
    #Top ground
    TGround = []
    
    c1, d1 = x1, y1 + s
    TGround.append([c1, d1])
    
    c2, d2 = c1, d1 + (Hr*10)/4 - s - w/2
    TGround.append([c2, d2])
    
    c3, d3 = a1, d2
    TGround.append([c3, d3])
    
    c4, d4 = c3, d3 - ((Hr*10)/4 - s - w/2)
    TGround.append([c4, d4])
    
    c5, d5 = c4 - Lline + Lc, d4
    TGround.append([c5, d5])
    
    c6, d6 = c5, d5 + ((Hr*10)/4 - s - w/2) - T
    TGround.append([c6, d6])
    
    c7, d7 = c1 + Lline - Lc, d6
    TGround.append([c7, d7])
    
    c8, d8 = c7, d7 - ((Hr*10)/4 - s - w/2) + T
    TGround.append([c8, d8])
    
    c9, d9 = c1, d1
    TGround.append([c9, d9])
    
    
    #Bottom ground
    BGround = []
    
    e1, f1 = x8, y8 - s
    BGround.append([e1, f1])
    
    e2, f2 = e1, f1 - ((Hr*10)/4 - s - w/2)
    BGround.append([e2, f2])
    
    e3, f3 = c3, f2
    BGround.append([e3, f3])
    
    e4, f4 = e3, f3 + ((Hr*10)/4 - s - w/2)
    BGround.append([e4, f4])
    
    e5, f5 = e4 - Lline + Lc, f4
    BGround.append([e5, f5])
    
    e6, f6 = e5, f5 - ((Hr*10)/4 - s - w/2) + T
    BGround.append([e6, f6])
    
    e7, f7 = c1 + Lline - Lc, f6
    BGround.append([e7, f7])
    
    e8, f8 = e7, f7 + ((Hr*10)/4 - s - w/2) - T
    BGround.append([e8, f8])
    
    e9, f9 = e1, f1
    BGround.append([e9, f9])
    
    
    Lcpw = np.asarray(Lcpw)
    Rcpw = np.asarray(Rcpw)
    TGround = np.asarray(TGround)
    BGround = np.asarray(BGround)
    return Lcpw, Rcpw, TGround, BGround

def Cpw_twirled(center1, center2, align_y, Hr, Lr, Sr, L, wg, I, w = 10e-6, Lc = 215e-6, Lline = 470e-6, T = 50e-6, s = 7.1e-6):
    Hr = Hr*1e6
    Sr = Sr*1e6
    Lr = Lr*1e6
    w = w*1e6
    Lc = Lc*1e6
    Lline = Lline*1e6
    T = T*1e6
    s = s*1e6
    I = I*1e6
    L = L*1e6
    wg = wg*1e6
    
   #Left CPW
    Lcpw = []
    
    x1, y1 = center1[0] - Lline - 1.5*Lr - Sr + wg/2, align_y + wg/2 - Hr/2 + w/2
    Lcpw.append([x1, y1])
    
    x2, y2 = x1 + Lline, y1
    Lcpw.append([x2, y2])
    
    x3, y3 = x2, y2 - w
    Lcpw.append([x3, y3])
    
    x4, y4 = x3 - Lline, y3
    Lcpw.append([x4, y4])
    
    x5, y5 = x1, y1
    Lcpw.append([x5, y5])
    
    #calculate start head coordinates x0+W/2,y0-W/2
    startx = (center1[0] - Lr - Sr)*1e-6
    starty = center1[1]*1e-6
    
    
    headl = TwirledSpiral([startx,starty], L*1e-6, wg*1e-6, I*1e-6, [0,0])
    
    
    #Right CPW
    Rcpw = []
    
    a1, b1 = center2[0] + Lline + 1.5*Lr + Sr - wg/2, align_y + wg/2 - Hr/2 + w/2
    Rcpw.append([a1, b1])
    
    a2, b2 = a1 - Lline, b1
    Rcpw.append([a2, b2])
    
    a3, b3 = a2, b2 - w
    Rcpw.append([a3, b3])
    
    a4, b4 = a3 + Lline, b3
    Rcpw.append([a4, b4])
    
    a5, b5 = a1, b1
    Rcpw.append([a5, b5])
    
    #calculate start head coordinates 
    gox = (center2[0] + Lr + Sr )*1e-6 
    goy = center2[1]*1e-6
    
    headr = TwirledSpiral([gox,goy], L*1e-6, wg*1e-6, I*1e-6, [0,0])

    
    #Top ground
    TGround = []
    
    c1, d1 = x1, y1 + s
    TGround.append([c1, d1])
    
    c2, d2 = c1, d1 + (Hr*10)/4 - s - w/2
    TGround.append([c2, d2])
    
    c3, d3 = a1, d2
    TGround.append([c3, d3])
    
    c4, d4 = c3, d1
    TGround.append([c4, d4])
    
    c5, d5 = c4 - Lline + Lc, d4
    TGround.append([c5, d5])
    
    c6, d6 = c5, d3 - T
    TGround.append([c6, d6])
    
    c7, d7 = c1 + Lline - Lc, d6
    TGround.append([c7, d7])
    
    c8, d8 = c7, d1
    TGround.append([c8, d8])
    
    c9, d9 = c1, d1
    TGround.append([c9, d9])
    
    
    #Bottom ground
    BGround = []
    
    e1, f1 = x4, y4 - s
    BGround.append([e1, f1])
    
    e2, f2 = e1, f1 - ((Hr*10)/4 - s - w/2)
    BGround.append([e2, f2])
    
    e3, f3 = c3, f2
    BGround.append([e3, f3])
    
    e4, f4 = e3, f1
    BGround.append([e4, f4])
    
    e5, f5 = e4 - Lline + Lc, f4
    BGround.append([e5, f5])
    
    e6, f6 = e5, f3 + T
    BGround.append([e6, f6])
    
    e7, f7 = e1 + Lline - Lc, f6
    BGround.append([e7, f7])
    
    e8, f8 = e7, f1
    BGround.append([e8, f8])
    
    e9, f9 = e1, f1
    BGround.append([e9, f9])
    
    
    Lcpw = np.asarray(Lcpw)
    Rcpw = np.asarray(Rcpw)
    TGround = np.asarray(TGround)
    BGround = np.asarray(BGround)
    return Lcpw, headl, Rcpw, headr, TGround, BGround

def Cpw_finger(El, Er, Hr, Sr, w = 10e-6, Lc = 215e-6, Lline = 470e-6, T = 50e-6, s = 7.1e-6):
    El = El[0]*1e6 ,El[1]*1e6
    Er = Er[0]*1e6, Er[1]*1e6
    Hr = Hr*1e6
    Sr = Sr*1e6
    w = w*1e6
    Lc = Lc*1e6
    Lline = Lline*1e6
    T = T*1e6
    s = s*1e6
    
     #Left CPW
    Lcpw = []
    
    x1, y1 = El[0] - Lline - Sr, El[1] + w/2
    Lcpw.append([x1, y1])
    
    x2, y2 = x1 + Lline, y1
    Lcpw.append([x2, y2])
    
    x3, y3 = x2, y2 - w
    Lcpw.append([x3, y3])
    
    x4, y4 = x3 - Lline, y3
    Lcpw.append([x4, y4])
    
    x5, y5 = x1, y1
    Lcpw.append([x5, y5])
    
    
    
    #Right CPW
    Rcpw = []
    
    a1, b1 = Er[0] + Lline + Sr, Er[1] + w/2
    Rcpw.append([a1, b1])
    
    a2, b2 = a1 - Lline, b1
    Rcpw.append([a2, b2])
    
    a3, b3 = a2, b2 - w
    Rcpw.append([a3, b3])
    
    a4, b4 = a3 + Lline, b3
    Rcpw.append([a4, b4])
    
    a5, b5 = a1, b1
    Rcpw.append([a5, b5])
      
    
    #Top ground
    TGround = []
    
    c1, d1 = x1, y1 + s
    TGround.append([c1, d1])
    
    c2, d2 = c1, d1 + (Hr*10)/4 - s - w/2
    TGround.append([c2, d2])
    
    c3, d3 = a1, d2
    TGround.append([c3, d3])
    
    c4, d4 = c3, d3 - ((Hr*10)/4 - s - w/2)
    TGround.append([c4, d4])
    
    c5, d5 = c4 - Lline + Lc, d4
    TGround.append([c5, d5])
    
    c6, d6 = c5, d5 + ((Hr*10)/4 - s - w/2) - T
    TGround.append([c6, d6])
    
    c7, d7 = c1 + Lline - Lc, d6
    TGround.append([c7, d7])
    
    c8, d8 = c7, d7 - ((Hr*10)/4 - s - w/2) + T
    TGround.append([c8, d8])
    
    c9, d9 = c1, d1
    TGround.append([c9, d9])
    
    
    #Bottom ground
    BGround = []
    
    e1, f1 = x4, y4 - s
    BGround.append([e1, f1])
    
    e2, f2 = e1, f1 - ((Hr*10)/4 - s - w/2)
    BGround.append([e2, f2])
    
    e3, f3 = c3, f2
    BGround.append([e3, f3])
    
    e4, f4 = e3, f3 + ((Hr*10)/4 - s - w/2)
    BGround.append([e4, f4])
    
    e5, f5 = e4 - Lline + Lc, f4
    BGround.append([e5, f5])
    
    e6, f6 = e5, f5 - ((Hr*10)/4 - s - w/2) + T
    BGround.append([e6, f6])
    
    e7, f7 = c1 + Lline - Lc, f6
    BGround.append([e7, f7])
    
    e8, f8 = e7, f7 + ((Hr*10)/4 - s - w/2) - T
    BGround.append([e8, f8])
    
    e9, f9 = e1, f1
    BGround.append([e9, f9])
    
    
    Lcpw = np.asarray(Lcpw)
    Rcpw = np.asarray(Rcpw)
    TGround = np.asarray(TGround)
    BGround = np.asarray(BGround)
    return Lcpw, Rcpw, TGround, BGround

"""
--- Twirled Spiral: new
"""

def CompTwirlSpir(L, t, centre = (0,0), m = 10): #We can still change the name to something clearer/simpler !
    

    #GDSPY works in um.
    L = L*1e6
    t = t*1e6
    #Ip = I + t #Here we define the spacing, Ip, defined for the path we will define below. It makes syntax easier to understand.
    
    #The part below find a minimal value for I such that I > 10*t
    It = 0
    S = 1
    i = 1
    while True:
        if (L)/(2*S + 1 -2) < (m+1)*t:
            print('The minimal value for the interspacing such that I>10*t is:',It,'um')
            print('aaaa',(It+t)*(2*(S-i) + 1 - 2))
            N = i-1
            break
        else:
            It = (L)/(2*S + 1 -2) - t
            print('It',It)
            Ltemp = (It+t)*(2*S + 1 - 2)
            print('Ltemp',Ltemp)
        i += 1
        S += i
    print(i)
    print(L)
    LB = (It+t)*(2*(S-i)+1-2)
    print('THIS IS THE CALCULATED L',LB)
    #Now we have the minimal for I that we can have, such that we limit the strength of the stray capacitance.
    I = It+t#We now add the thickness of the arms so it looks cleaner in the loops below
    
    #Now we define the structure
    #As we before we work with right-handed and left-handed arms, we define each arm's length to be L/2
    
    Ltest = 0 #Value to check if we have the right length at the end of the loops
    
    #Right-handed part
    Pr = [] #Pr stands for path right
    xr0, yr0 = centre
    Pr.append([xr0,yr0])
    
    xr1, yr1 = xr0 + I/2, yr0
    Pr.append([xr1,yr1])
    
    Ltest += I/2
    
    
    k = 0
    while True:
        k+=1
        if k > (N-1):
            xr_up, yr_up = Pr[-1][0], Pr[-1][1] + (k-1)*I
            Pr.append([xr_up,yr_up])
            Ltest += ((k-1)*I)
            break
        else:
            xr_up, yr_up = Pr[-1][0], Pr[-1][1] + k*I
            Pr.append([xr_up,yr_up])
            Ltest += k*I
        k += 1
        if k > (N-1):
            xr_left, yr_left = Pr[-1][0] - (k-1)*I, Pr[-1][1] 
            Pr.append([xr_left,yr_left])
            Ltest += (k-1)*I
            break
        else:
            xr_left, yr_left = Pr[-1][0] - k*I, Pr[-1][1]
            Pr.append([xr_left,yr_left])
            Ltest += k*I
        k+=1
        if k > (N-1):
            xr_down, yr_down = Pr[-1][0], Pr[-1][1] - (k-1)*I
            Pr.append([xr_down,yr_down])
            Ltest += (k-1)*I
            break
        else:
            xr_down, yr_down = Pr[-1][0], Pr[-1][1] - k*I
            Pr.append([xr_down,yr_down])
            Ltest += k*I
        k+=1
        if k > (N-1):
            xr_right, yr_right = Pr[-1][0] + (k-1)*I, Pr[-1][1]
            Pr.append([xr_right, yr_right])
            Ltest += (k-1)*I
            break
        else:
            xr_right, yr_right = Pr[-1][0] + k*I, Pr[-1][1]
            Pr.append([xr_right, yr_right])
            Ltest += k*I
    Pr = np.asarray(Pr)

        
    #Left-handed part
    Pl = []
    xl0, yl0 = centre
    Pl.append([xl0,yl0])
    
    xl1, yl1 = xl0 - I/2, yl0
    Pl.append([xl1, yl1])
    
    Ltest += I/2
    
    k = 0
    while True:
        k+=1
        if k > (N-1):
            xl_down, yl_down = Pl[-1][0], Pl[-1][1] - (k-1)*I
            Pl.append([xl_down,yl_down])
            Ltest += (k-1)*I
            break
        else:
            xl_down, yl_down = Pl[-1][0], Pl[-1][1] - k*I
            Pl.append([xl_down,yl_down])
            Ltest += k*I
        k += 1
        if k > (N-1):
            xl_right, yl_right = Pl[-1][0] + (k-1)*I, Pl[-1][1]
            Pl.append([xl_right,yl_right])
            Ltest += (k-1)*I
            break
        else:
            xl_right, yl_right = Pl[-1][0] + k*I, Pl[-1][1]
            Pl.append([xl_right,yl_right])
            Ltest += k*I
        k+=1
        if k > (N-1):
            xl_up, yl_up = Pl[-1][0], Pl[-1][1] + (k-1)*I
            Pl.append([xl_up,yl_up])
            Ltest += (k-1)*I
            break
        else:
            xl_up, yl_up = Pl[-1][0], Pl[-1][1] + k*I
            Pl.append([xl_up,yl_up])
            Ltest += k*I
        k+=1
        if k > (N-1):
            xl_left, yl_left = Pl[-1][0] - (k-1)*I, Pl[-1][1]
            Pl.append([xl_left, yl_left])
            Ltest += (k-1)*I
            break
        else:
            xl_left, yl_left = Pl[-1][0] - k*I, Pl[-1][1]
            Pl.append([xl_left, yl_left])
            Ltest += k*I
    Pl = np.asarray(Pl)
    print('Ltest',Ltest)
    print(L-Ltest)
    
    return Pr, Pl, #L-Ltest

def CompTwirlSpirPad(L, t, Lcouple, Pad, centre = (0,0), m = 10): #We can still change the name to something clearer/simpler !
    """
    This function returns three list of points to generate a spiral resonator with an elongated side with a pad to increase the coupling
    Parameters:
        L: The total length of the resonator (um)
        t: The width of the arms of your resonator v
        Lcouple: The length of the couping arm (um)
        Pad: The dimensions of the pad (um)
        centre: The position of the centre of your resonator
        m: The rule on your interspacing such that I > m * t
    Output:
        Pr: List of points for the right arm of the resonator
        Pl: List of points for the left arm of the resonator
        Ptop: List of points to generate the pad at the end of your left arm
        L - Ltest: To test if you obtain the right length at the end.
    
    """
    #GDSPY works in um.
    #This function differ than the other, here we will build a spiral up to the N-1 branch so that we can build the last branch to couple
    #to the feedline and building the two pads
    L = L*1e6
    Lcouple  = Lcouple*1e6
    PT = Pad[0]*1e6, Pad[1]*1e6
    t = t*1e6
    
    
    La = L - Lcouple
    
    
    #Ip = I + t #Here we define the spacing, Ip, defined for the path we will define below. It makes syntax easier to understand.
    
    #The part below find a minimal value for I such that I > 10*t
    It = 0
    S = 1
    i = 1
    while True:
        if (La)/(2*S + i+1) < m*t:
            # print('The minimal value for the interspacing such that I>10*t is:',It,'um')
            # print('aaaa',It*(2*S + 1 -2))
            N = i-1
            break
        else:
            It = (La)/(2*S + i+1)
            # print('It',It)
            Ltemp = It*(2*S + i+1)
            # print('Ltemp',Ltemp)
        i += 1
        S += i
    # print(i)
    # print(L)
    LB = It*(2*(S-i)+1-2)
    # print('THIS IS THE CALCULATED L',LB)
    #Now we have the minimal for I that we can have, such that we limit the strength of the stray capacitance.
    I = It#We now add the thickness of the arms so it looks cleaner in the loops below
    
    #Now we define the structure
    #As we before we work with right-handed and left-handed arms, we define each arm's length to be L/2
    
    Ltest = 0 #Value to check if we have the right length at the end of the loops
    
    #Right-handed part
    Pr = [] #Pr stands for path right
    xr0, yr0 = centre
    Pr.append([xr0,yr0])
    
    xr1, yr1 = xr0 + I/2, yr0
    Pr.append([xr1,yr1])
    
    Ltest += I/2
    
    k = 0
    while True:
        k+=1
        if k > (N):
            xr_up, yr_up = Pr[-1][0], Pr[-1][1] + (k-1)*I
            Pr.append([xr_up,yr_up])
            Ltest += (k-1)*I
            break
        else:
            xr_up, yr_up = Pr[-1][0], Pr[-1][1] + k*I
            Pr.append([xr_up,yr_up])
            Ltest += k*I
        k += 1
        if k > (N):
            xr_left, yr_left = Pr[-1][0] - (k-1)*I, Pr[-1][1]
            Pr.append([xr_left,yr_left])
            Ltest += (k-1)*I
            break
        else:
            xr_left, yr_left = Pr[-1][0] - k*I, Pr[-1][1]
            Pr.append([xr_left,yr_left])
            Ltest += k*I
        k+=1
        if k > (N):
            xr_down, yr_down = Pr[-1][0], Pr[-1][1] - (k-1)*I
            Pr.append([xr_down,yr_down])
            Ltest += (k-1)*I
            break
        else:
            xr_down, yr_down = Pr[-1][0], Pr[-1][1] - k*I
            Pr.append([xr_down,yr_down])
            Ltest += k*I
        k+=1
        if k > (N):
            xr_right, yr_right = Pr[-1][0] + (k-1)*I, Pr[-1][1]
            Pr.append([xr_right, yr_right])
            Ltest += (k-1)*I
            break
        else:
            xr_right, yr_right = Pr[-1][0] + k*I, Pr[-1][1]
            Pr.append([xr_right, yr_right])
            Ltest += k*I
    Pr = np.asarray(Pr)

        
    #Left-handed part
    Ptop = []
    Pl = []
    xl0, yl0 = centre
    Pl.append([xl0,yl0])
    
    xl1, yl1 = xl0 - I/2, yl0
    Pl.append([xl1, yl1])
    
    Ltest += I/2
    
    k = 0
    while True:
        k+=1
        if k > (N):
            # xl_down, yl_down = Pl[-1][0], Pl[-1][1] - (k-1)*I
            # Pl.append([xl_down,yl_down])
            # Ltest += (k-1)*I
            xl_down, yl_down = Pl[-1][0], Pl[-1][1] - Lcouple
            Pl.append([xl_down,yl_down])
            Ltest += Lcouple
            Ptop.append([xl_down-t,yl_down])#Ptop.append([xl_down-t/2,yl_down])
            xl_square, yl_square = Ptop[-1][0] + PT[0], Ptop[-1][1] + PT[1]
            Ptop.append([xl_square, yl_square])
            break
        else:
            xl_down, yl_down = Pl[-1][0], Pl[-1][1] - k*I
            Pl.append([xl_down,yl_down])
            Ltest += k*I
        k += 1
        if k > (N):
            # xl_right, yl_right = Pl[-1][0] + (k-1)*I, Pl[-1][1]
            # Pl.append([xl_right,yl_right])
            # Ltest += (k-1)*I
            xl_right, yl_right = Pl[-1][0] + Lcouple, Pl[-1][1]
            Pl.append([xl_right,yl_right])
            Ltest += Lcouple
            Ptop.append([xl_right,yl_right-t/2])
            xl_square, yl_square = Ptop[-1][0] - PT[0], Ptop[-1][1] + PT[1]
            Ptop.append([xl_square, yl_square])
            break
        else:
            xl_right, yl_right = Pl[-1][0] + k*I, Pl[-1][1]
            Pl.append([xl_right,yl_right])
            Ltest += k*I
        k+=1
        if k > (N):
            # xl_up, yl_up = Pl[-1][0], Pl[-1][1] + (k-1)*I
            # Pl.append([xl_up,yl_up])
            # Ltest += (k-1)*I
            xl_up, yl_up = Pl[-1][0], Pl[-1][1] + Lcouple
            Pl.append([xl_up,yl_up])
            Ltest += Lcouple
            Ptop.append([xl_up+t/2,yl_up])
            xl_square, yl_square = Ptop[-1][0] - PT[0], Ptop[-1][1] - PT[1]
            Ptop.append([xl_square, yl_square])
            break
        else:
            xl_up, yl_up = Pl[-1][0], Pl[-1][1] + k*I
            Pl.append([xl_up,yl_up])
            Ltest += k*I
        k+=1
        if k > (N):
            # xl_left, yl_left = Pl[-1][0] - (k-1)*I, Pl[-1][1]
            # Pl.append([xl_left, yl_left])
            # Ltest += (k-1)*I
            xl_left, yl_left = Pl[-1][0] - Lcouple, Pl[-1][1]
            Pl.append([xl_left, yl_left+t/2])
            Ltest += Lcouple
            Ptop.append([xl_left, yl_left])
            xl_square, yl_square = Ptop[-1][0] + PT[0], Ptop[-1][1] - PT[1]
            Ptop.append([xl_square, yl_square])
            break
        else:
            xl_left, yl_left = Pl[-1][0] - k*I, Pl[-1][1]
            Pl.append([xl_left, yl_left])
            Ltest += k*I
    Pl = np.asarray(Pl)
    # print('Ltest',Ltest)
    # print(L-Ltest)
    
    
    return Pr, Pl, Ptop#, L - Ltest



"""
--- meandered Spirals ------------
"""

def meandered(H,D,m,w,start=(0,0)):
    """

    Parameters
    ----------
    H : height of resonator in SI units
    D : big spacing in SI units
    m : small spacing in SI units
    w : width in SI units
    start : coordinates of bottomleft corner. The default is (0,0).

    Returns
    -------
    coord : array of points [x,y]
        gives the coordinates to draw the meander with FlexPath
    middle : (xc,yc)
        coordinates of center of meander for rotation

    """
    H = H*1e6 #height of resonator
    D = D*1e6 #big spacing
    m = m*1e6 #small spacing
    w = w*1e6 #width of path/wire
    
    
    coord = []
    x0,y0 = start[0],start[1]
    coord.append([x0,y0])

    #this for loop draws the first big "Π"
    for j in range(4):
        # print(j)
        if j == 0:
            x0,y0 = coord[-1]
            coord.append([x0+w/2,y0])
        elif j == 1:
            x1,y1 = coord[-1]
            coord.append([x1,y1+H-w/2])
        elif j == 2:
            x2,y2 = coord[-1]
            coord.append([x2+D+w,y2])
        elif j == 3:
            x3,y3 = coord[-1]
            coord.append([x3,y3-H+w])
            
    #draw 3 x the small meanders        
    for l in range(3):
        # print(l)
    
        #draw 1 small meander ("Π")
        for k in range(4):
            # print('x')
            if k == 0:
                x4,y4 = coord[-1]
                coord.append([x4+w+m,y4])
            elif k == 1:
                x5,y5 = coord[-1]
                coord.append([x5,y5+H-w])
            elif k == 2:
                x6,y6 = coord[-1]
                coord.append([x6+w+m,y6])
            elif k == 3:
                x7,y7 = coord[-1]
                coord.append([x7,y7-H+w])
    x8,y8 = coord[-1]
    coord.append([x8+w+m,y8])
    
    #draw last big meander
    for i in range(3):
        # print(j)
        if i == 0:
            x1,y1 = coord[-1]
            coord.append([x1,y1+H-w])
        elif i == 1:
            x2,y2 = coord[-1]
            coord.append([x2+D+w,y2])
        elif i == 2:
            x3,y3 = coord[-1]
            coord.append([x3,y3-H+w/2])
            
    #discard starting point (not needed to draw) and calculate center coordinates of resonator (used for rotation)            
    del coord[0]
    middle = ((coord[-1][0]-coord[0][0]+w)/2+start[0],H/2+start[1])
    return coord, middle

def Cpw_meandered(snake1, snakedim, end, Sr, w = 10e-6, Lc = 215e-6, Lline = 470e-6, T = 50e-6, s = 7.1e-6):
    """
    

    Parameters
    ----------
    snake1 : array of coordinates of first resonator
    snakedim : array with dimensions of one resonator in SI units
    end : last coordinate of last resonator
    Sr : spacing between feedline and resonator in SI units
    w : width of the feedline in SI units. The default is 10e-6.
    Lc : #spacing between feedline "neck" and ground plane corner. The default is 215e-6.
    Lline : length of feedline without "head". The default is 470e-6.
    T : width of narrow junk of ground planes. The default is 50e-6.
    s : spacing between feedline and ground planes. The default is 7.1e-6.

    Returns
    -------
    Lcpw : coordinates-array to create Polygon for left feedline "arm"
    headl : coordinates-array to create FlexPath for left feedline "head" (meander)
    centerl : tuple of center coordinates of left "head" (used for rotation)
    Rcpw : coordinates-array to create polygon for right feedline "arm"
    headr : coordinates-array to create FlexPath for left feedline "head" (meander)
    TGround : coordinates-array to create polygon as top ground plate
    BGround : coordinates-array to create polygon as bottom ground plate

    """
    snake1 = np.array(snake1) #coordinates of first spiral
    H = snakedim[0]*1e6 #height of resonator
    D = snakedim[1]*1e6
    m = snakedim[2]*1e6
    ws = snakedim[3]*1e6 #width of resonator
    end = end[0], end[1] #last coordinate of last spiral
    
    Sr = Sr*1e6 #spacing feedline to 1st resonator
    w = w*1e6 #width feedline
    Lc = Lc*1e6 #spacing between feedline "neck" and ground plane corner
    Lline = Lline*1e6 #length of feedline without "head"
    T = T*1e6 #width of narrow junk of ground planes
    s = s*1e6 #spacing between feedline and ground planes
    t = snake1[-1][0] - snake1[0][0] + ws #length of feedline "head"
    
    #define left and right edges
    El = snake1[0][0] - ws/2, snake1[0][1] + H/2
    Er = end[0] + ws/2, end[1] + H/2
    
    
    #Left CPW
    Lcpw = []
    
    x1, y1 = El[0] - Lline - t - Sr, El[1] + w/2
    Lcpw.append([x1, y1])
    
    x2, y2 = x1 + Lline, y1
    Lcpw.append([x2, y2])
    
    x3, y3 = x2, y2 - w
    Lcpw.append([x3, y3])
    
    x4, y4 = x3 - Lline, y3
    Lcpw.append([x4, y4])
    
    x5, y5 = x1, y1
    Lcpw.append([x5, y5])
    
    #calculate start head coordinates
    startx = x2
    starty = y3 - H/2 + w/2
    
    headl, centerl = meandered(H*1e-6,D*1e-6,m*1e-6,ws*1e-6,start=(startx,starty))
    
    
    #Right CPW
    Rcpw = []
    
    a1, b1 = Er[0] + Lline + t + Sr, Er[1] + w/2
    Rcpw.append([a1, b1])
    
    a2, b2 = a1 - Lline, b1
    Rcpw.append([a2, b2])
    
    a3, b3 = a2, b2 - w
    Rcpw.append([a3, b3])
    
    a4, b4 = a3 + Lline, b3
    Rcpw.append([a4, b4])
    
    a5, b5 = a1, b1
    Rcpw.append([a5, b5])
    
    #calculate start head coordinates
    gox = end[0] + ws/2 + Sr
    goy = end[1]
    
    headr, trashr = meandered(H*1e-6, D*1e-6, m*1e-6, ws*1e-6, start =(gox,goy))

    
    #Top ground
    TGround = []
    
    c1, d1 = x1, y1 + s
    TGround.append([c1, d1])
    
    c2, d2 = c1, d1 + (H*10)/4 - s - w/2
    TGround.append([c2, d2])
    
    c3, d3 = a1, d2
    TGround.append([c3, d3])
    
    c4, d4 = c3, d1
    TGround.append([c4, d4])
    
    c5, d5 = c4 - Lline + Lc, d4
    TGround.append([c5, d5])
    
    c6, d6 = c5, d3 - T
    TGround.append([c6, d6])
    
    c7, d7 = c1 + Lline - Lc, d6
    TGround.append([c7, d7])
    
    c8, d8 = c7, d1
    TGround.append([c8, d8])
    
    c9, d9 = c1, d1
    TGround.append([c9, d9])
    
    
    #Bottom ground
    BGround = []
    
    e1, f1 = x4, y4 - s
    BGround.append([e1, f1])
    
    e2, f2 = e1, f1 - ((H*10)/4 - s - w/2)
    BGround.append([e2, f2])
    
    e3, f3 = c3, f2
    BGround.append([e3, f3])
    
    e4, f4 = e3, f1
    BGround.append([e4, f4])
    
    e5, f5 = e4 - Lline + Lc, f4
    BGround.append([e5, f5])
    
    e6, f6 = e5, f3 + T
    BGround.append([e6, f6])
    
    e7, f7 = e1 + Lline - Lc, f6
    BGround.append([e7, f7])
    
    e8, f8 = e7, f1
    BGround.append([e8, f8])
    
    e9, f9 = e1, f1
    BGround.append([e9, f9])
    
    
    Lcpw = np.asarray(Lcpw)
    Rcpw = np.asarray(Rcpw)
    TGround = np.asarray(TGround)
    BGround = np.asarray(BGround)
    return Lcpw, headl, centerl, Rcpw, headr, TGround, BGround


"""
------------------------------------------------------------------------------
left-handed LC array
------------------------------------------------------------------------------
"""

def OPKI_down(L,s,w,A,t, centre=(0,0),return_center=False,compact=False):
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
    if compact == False:
        N = int((L-2*(s+w))/(b+s))
    else:
        N = round((L-2*(s+w))/(b+s))
    
    #adapt start & end segment of inductor
    d_prime = 0.5*(L - N*(s+b)) + w
    #print('d_prime: ',d_prime)
    
    #calculate vertical dimension of capacitor
    B = N*(s+w) + d_prime + t
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
        return U,M,B
    if return_center == True:
        return U, M,B, centre

def OPKI_up(L,s,w,A,t, centre=(0,0),return_center = False, compact = False):
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
    if compact == False:
        N = int((L-2*(s+w))/(b+s))
    else:
        N = round((L-2*(s+w))/(b+s))
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
        return U,M,B
    if return_center == True:
        return U, M, B, centre


def grounded_Res(L,s,w,A,t,tv,e=20.5e-6,f=24e-6,r=29e-6,gamma=1/4,ground_in_between=True, center = (0,0)):
    t = t*1e6
    w = w*1e6
    e = e*1e6
    f = f*1e6
    r = r*1e6
    tv = tv*1e6
    
    carac1 = {'layer': 0, 'datatype':3}
    
    Ushape, Mshape, B = OPKI_down(L,s,w*1e-6,A,t*1e-6)
    U = gdspy.FlexPath(Ushape, t, **carac1)
    M = gdspy.FlexPath(Mshape,w,**carac1)
    
    oneRes = gdspy.Cell('1 grounded LC Resonator')
    oneRes.add(U)
    oneRes.add(M)
    
    box_x1, box_y1 = Mshape[-1][0] - e/2 - w/2, Mshape[-1][1]
    box_x2, box_y2 = Mshape[-1][0] + w/2 + e/2, Mshape[-1][1]-f
    
    carac2 = {'layer': 1, 'datatype': 3}
    ground_box = gdspy.Rectangle([box_x1,box_y1], [box_x2,box_y2],**carac2)
    oneRes.add(ground_box)
    
    strip_x1, strip_y1 = Ushape[-1][0] + t/2 + (1-gamma)/2 * tv, box_y2
    strip_x2, strip_y2 = strip_x1 + gamma*tv, Ushape[-2][1] + t/2 + r
    
    if ground_in_between==True:
        ground_strip = gdspy.Rectangle([strip_x1,strip_y1], [strip_x2,strip_y2],**carac2)
        oneRes.add(ground_strip)
        
    return oneRes

def ghosts(L,s,w,A,t,tv,tw,N, strip_height, tg = 15e-6, e=20.5e-6, f=24e-6, r=29e-6, gamma=1/4, ground_in_between=True, carac = {'layer' :  1, 'datatype' : 1}, center_first = (0,0), center_last = (0,0)):
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
            
            Ushape_l, Mshape_l, B= OPKI_up(L,s,w*1e-6,A*1e-6,t*1e-6,centre=center_left)
            Ushape_r, Mshape_r, B= OPKI_down(L,s,w*1e-6,A*1e-6,t*1e-6,centre=center_right)
            U_l = gdspy.FlexPath(Ushape_l, t,**carac)
            M_l = gdspy.FlexPath(Mshape_l,w,**carac)
            U_r = gdspy.FlexPath(Ushape_r, t,**carac)
            M_r = gdspy.FlexPath(Mshape_r,w,**carac)
            # ghostcell.add([U_l,U_r])
            ghostcell.add([U_l,M_l,U_r,M_r])
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
            Ushape_l,Mshape_l, B = OPKI_down(L, s, w*1e-6,A*1e-6,t*1e-6, centre=center_left)
            Ushape_r,Mshape_r, B = OPKI_up(L, s, w*1e-6,A*1e-6,t*1e-6, centre = center_right)
            U_l = gdspy.FlexPath(Ushape_l, t, **carac)
            M_l = gdspy.FlexPath(Mshape_l,w,**carac)
            U_r = gdspy.FlexPath(Ushape_r, t, **carac)
            M_r = gdspy.FlexPath(Mshape_r,w,**carac)
            # ghostcell.add([U_l,U_r])
            ghostcell.add([U_l,M_l,U_r,M_r])
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
    return ghostcell,ground
  
def ghosts_U(L,s,w,A,t,tv,tw,N, strip_height, tg = 15e-6, e=20.5e-6, f=24e-6, r=29e-6, gamma=1/4, ground_in_between=True, carac = {'layer' :  1, 'datatype' : 1}, center_first = (0,0), center_last = (0,0)):
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
            
            Ushape_l, Mshape_l, B= OPKI_up(L,s,w*1e-6,A*1e-6,t*1e-6,centre=center_left)
            Ushape_r, Mshape_r, B= OPKI_down(L,s,w*1e-6,A*1e-6,t*1e-6,centre=center_right)
            U_l = gdspy.FlexPath(Ushape_l, t,**carac)
            # M_l = gdspy.FlexPath(Mshape_l,w,**carac)
            U_r = gdspy.FlexPath(Ushape_r, t,**carac)
            # M_r = gdspy.FlexPath(Mshape_r,w,**carac)
            ghostcell.add([U_l,U_r])
            # ghostcell.add([U_l,M_l,U_r,M_r])
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
            Ushape_l,Mshape_l, B = OPKI_down(L, s, w*1e-6,A*1e-6,t*1e-6, centre=center_left)
            Ushape_r,Mshape_r, B = OPKI_up(L, s, w*1e-6,A*1e-6,t*1e-6, centre = center_right)
            U_l = gdspy.FlexPath(Ushape_l, t, **carac)
            # M_l = gdspy.FlexPath(Mshape_l,w,**carac)
            U_r = gdspy.FlexPath(Ushape_r, t, **carac)
            # M_r = gdspy.FlexPath(Mshape_r,w,**carac)
            ghostcell.add([U_l,U_r])
            # ghostcell.add([U_l,M_l,U_r,M_r])
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
                strip_x1_N = center_left[0] - A/2 - xw
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
    
    return ghostcell,ground

def ghosts_U_GGG(L,s,w,A,t,tv,tw,N, strip_height, tg = 15e-6, e=20.5e-6, f=24e-6, r=29e-6, ts=3e-6, ground_in_between=True, carac = {'layer' :  1, 'datatype' : 1}, center_first = (0,0), center_last = (0,0)):
    A = A*1e6
    t = t*1e6
    w = w*1e6
    e = e*1e6
    f = f*1e6
    r = r*1e6
    tv = tv*1e6
    tw = tw*1e6
    tg = tg*1e6
    ts = ts*1e6
    center_left = center_first
    center_right = center_last
    xw = tw - 2*ts
    xv = tv - 2*ts
    xg = tg - 2*ts
    
    if xw < 0 or xv < 0 or xg < 0:
        print('ERROR. Ground strip width negative.')
    
    ghostcell = gdspy.Cell('ghosts')
    ground = gdspy.Cell('ghosts ground')
    
    if N == 0 and ground_in_between == True:
        strip_x1_first = center_left[0] - A/2 - ts
        strip_y1_first = center_left[1] + r
        strip_x2_first = strip_x1_first - xg
        strip_y2_first = strip_y1_first - strip_height
        ground_left = gdspy.Rectangle([strip_x1_first,strip_y1_first], [strip_x2_first,strip_y2_first])
        strip_x1_last = center_last[0] + A/2 + ts
        strip_y1_last = center_last[1] - r
        strip_x2_last = strip_x1_last + xg
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
            
            Ushape_l, Mshape_l,B = OPKI_up(L,s,w*1e-6,A*1e-6,t*1e-6,centre=center_left)
            Ushape_r, Mshape_r, B= OPKI_down(L,s,w*1e-6,A*1e-6,t*1e-6,centre=center_right)
            U_l = gdspy.FlexPath(Ushape_l, t,**carac)
            # M_l = gdspy.FlexPath(Mshape_l,w,**carac)
            U_r = gdspy.FlexPath(Ushape_r, t,**carac)
            # M_r = gdspy.FlexPath(Mshape_r,w,**carac)
            ghostcell.add([U_l,U_r])
            # ghostcell.add([U_l,M_l,U_r,M_r])
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
                strip_x1_l = center_left[0] + A/2 + ts
                strip_y1_l = box_y4
                strip_x2_l = strip_x1_l + xg
                strip_y2_l = strip_y1_l - strip_height
                strip_x1_r = center_right[0] + A/2 + ts
                strip_y1_r = box_y2
                strip_x2_r = strip_x1_r + xv
                strip_y2_r = strip_y1_r + strip_height
            else:
                #ground strips
                strip_x1_l = center_left[0] + A/2 + ts
                strip_y1_l = box_y4
                strip_x2_l = strip_x1_l + xv
                strip_y2_l = strip_y1_l - strip_height
                strip_x1_r = center_right[0] + A/2 + ts
                strip_y1_r = box_y2
                strip_x2_r = strip_x1_r + xv
                strip_y2_r = strip_y1_r + strip_height
            
            if i == N:
                strip_x1_N = center_left[0] - A/2 - ts
                strip_y1_N = box_y4
                strip_x2_N = strip_x1_N - xv
                strip_y2_N = strip_y1_N - strip_height
                strip_x1_last = center_last[0] + A/2 + ts
                strip_y1_last = center_last[1] - r
                strip_x2_last = strip_x1_last + xg
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
            Ushape_l,Mshape_l, B = OPKI_down(L, s, w*1e-6,A*1e-6,t*1e-6, centre=center_left)
            Ushape_r,Mshape_r,B = OPKI_up(L, s, w*1e-6,A*1e-6,t*1e-6, centre = center_right)
            U_l = gdspy.FlexPath(Ushape_l, t, **carac)
            # M_l = gdspy.FlexPath(Mshape_l,w,**carac)
            U_r = gdspy.FlexPath(Ushape_r, t, **carac)
            # M_r = gdspy.FlexPath(Mshape_r,w,**carac)
            ghostcell.add([U_l,U_r])
            # ghostcell.add([U_l,M_l,U_r,M_r])
            #ground patches
            box_x1, box_y1 = Mshape_l[-1][0] - e/2 - w/2, Mshape_l[-1][1]
            box_x2, box_y2 = Mshape_l[-1][0] + w/2 + e/2, Mshape_l[-1][1]-f
            ground_box1 = gdspy.Rectangle([box_x1,box_y1], [box_x2,box_y2])
            box_x3, box_y3 = Mshape_r[-1][0] - e/2 - w/2, Mshape_r[-1][1]
            box_x4, box_y4 = Mshape_r[-1][0] + w/2 + e/2, Mshape_r[-1][1] + f
            ground_box2 = gdspy.Rectangle([box_x3,box_y3], [box_x4,box_y4])
            ground.add([ground_box1,ground_box2])
            #ground strips
            strip_x1_l = center_left[0] + A/2 + ts
            strip_y1_l = box_y2
            strip_x2_l = strip_x1_l + xv
            strip_y2_l = strip_y1_l + strip_height
            strip_x1_r = center_right[0] + A/2 + ts
            strip_y1_r = box_y4
            strip_x2_r = strip_x1_r + xw
            strip_y2_r = strip_y1_r - strip_height
            
            if i == N:
                strip_x1_N = center_left[0] - A/2 - ts
                strip_y1_N = box_y2
                strip_x2_N = strip_x1_N - xw
                strip_y2_N = strip_y1_N + strip_height
                strip_x1_last = center_last[0] + A/2 + ts
                strip_y1_last = center_last[1] - r
                strip_x2_last = strip_x1_last + xg
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
    
    return ghostcell,ground
    
def unit_cell(L,s,w,A,t,tv,tw,e=20.5e-6,f=24e-6,r=29e-6,carac = {'layer' : 0, 'datatype' : 3},gamma=1/4,ground_in_between=True,compactRes=False):
    A = A*1e6
    t = t*1e6
    w = w*1e6
    e = e*1e6
    f = f*1e6
    r = r*1e6
    tv = tv*1e6
    tw = tw*1e6
    xv = 0.5*(tv - gamma*min(tv,tw))
    xw = 0.5*(tw - gamma*min(tv,tw))
    # print(xv,xw)
    

    if compactRes == True:
        #get coordinates for first resonator
        Ushape_d, Mshape_d, B, centre_d = OPKI_down(L,s,w*1e-6,A*1e-6,t*1e-6,return_center=(True), compact=True)
        
        #corners of 1st ground patch: down
        box_x1, box_y1 = Mshape_d[-1][0] - e/2 - w/2, Mshape_d[-1][1]
        box_x2, box_y2 = Mshape_d[-1][0] + w/2 + e/2, Mshape_d[-1][1]-f
        ground_box1 = gdspy.Rectangle([box_x1,box_y1], [box_x2,box_y2])
        
        #center of first resonator (used to draw the second one)
        centre_up = (Ushape_d[-1][0] + t/2 + tv + A/2,box_y2 + r)
        #get coordinates for second resonator
        Ushape_u, Mshape_u,B = OPKI_up(L,s,w*1e-6,A*1e-6,t*1e-6,centre=centre_up, compact = True)
    else:
        #get coordinates for first resonator
        Ushape_d, Mshape_d, B, centre_d = OPKI_down(L,s,w*1e-6,A*1e-6,t*1e-6,return_center=(True))
        
        #corners of 1st ground patch: down
        box_x1, box_y1 = Mshape_d[-1][0] - e/2 - w/2, Mshape_d[-1][1]
        box_x2, box_y2 = Mshape_d[-1][0] + w/2 + e/2, Mshape_d[-1][1]-f
        ground_box1 = gdspy.Rectangle([box_x1,box_y1], [box_x2,box_y2])
        
        #center of first resonator (used to draw the second one)
        centre_up = (Ushape_d[-1][0] + t/2 + tv + A/2,box_y2 + r)
        #get coordinates for second resonator
        Ushape_u, Mshape_u,B = OPKI_up(L,s,w*1e-6,A*1e-6,t*1e-6,centre=centre_up)
   
    #draw the two resonator
    U_down = gdspy.FlexPath(Ushape_d, t, **carac)
    M_down = gdspy.FlexPath(Mshape_d,w,**carac)  
    U_up = gdspy.FlexPath(Ushape_u, t, **carac)
    M_up = gdspy.FlexPath(Mshape_u, w, **carac)    
    
    #corners of 2nd ground patch: up
    box_x3, box_y3 = Mshape_u[-1][0] - e/2 - w/2, Mshape_u[-1][1]
    box_x4, box_y4 = Mshape_u[-1][0] + w/2 + e/2, Mshape_u[-1][1] + f
    ground_box2 = gdspy.Rectangle([box_x3,box_y3], [box_x4,box_y4])
    
    #corners of 1st ground strip in between resonators
    strip_x1, strip_y1 = Ushape_d[-1][0] + t/2 + xv, box_y2
    strip_x2, strip_y2 = strip_x1 + gamma*min(tv,tw), Ushape_d[-2][1] + t/2 + r
    #corners of 2nd ground strip in between resonators
    strip_x3, strip_y3 = Ushape_u[-1][0] + t/2 + xw, box_y4
    strip_x4, strip_y4 = strip_x3 + gamma*min(tv,tw), Ushape_u[-2][1] - t/2 - r
    
    
    if tv < tw:
        print('trivial')
    if tv == tw:
        print('normal')
    if tv > tw:
        print('topological')
    
    unitcell = gdspy.Cell('unit cell')
    ground = gdspy.Cell('ground')
    unitcell.add(U_down)
    unitcell.add(M_down)
    ground.add(ground_box1)
    unitcell.add(U_up)
    unitcell.add(M_up)
    ground.add(ground_box2)
    
    if ground_in_between==True:
        ground_strip1 = gdspy.Rectangle([strip_x1,strip_y1], [strip_x2,strip_y2])
        ground.add(ground_strip1)
        ground_strip2 = gdspy.Rectangle([strip_x3,strip_y3], [strip_x4,strip_y4])
        ground.add(ground_strip2)
        
    startM = Mshape_d[-1]
    startU = Ushape_d[-1]
    stopUy = Ushape_u[-1][1]
    strip_height = np.abs(strip_y2 - strip_y1)

    return unitcell, ground, startM, startU, stopUy, strip_height, B, centre_d, centre_up

def unit_cell_GGG(L,s,w,A,t,tv,tw,e=20.5e-6,f=24e-6,r=29e-6,carac = {'layer' : 0, 'datatype' : 3},ts=3e-6,ground_in_between=True,compactRes=False):
    A = A*1e6
    t = t*1e6
    w = w*1e6
    e = e*1e6
    f = f*1e6
    r = r*1e6
    tv = tv*1e6
    tw = tw*1e6
    ts = ts*1e6
    xv = tv - 2*ts
    xw = tw-2*ts
    if xv < 0 or xw < 0:
        print('ERROR: Cannot draw ground strips with those parameters, negative width!')
        return gdspy.Rectangle([-20,-20],[20,20]), gdspy.Rectangle([25,-20],[30,20]), [0,0], [0,0], [0,0], 40, 40, [0,0], [0,0]
    # print(xv,xw)
    

    if compactRes == True:
        #get coordinates for first resonator
        Ushape_d, Mshape_d, B, centre_d = OPKI_down(L,s,w*1e-6,A*1e-6,t*1e-6,return_center=(True), compact=True)
        
        #corners of 1st ground patch: down
        box_x1, box_y1 = Mshape_d[-1][0] - e/2 - w/2, Mshape_d[-1][1]
        box_x2, box_y2 = Mshape_d[-1][0] + w/2 + e/2, Mshape_d[-1][1]-f
        ground_box1 = gdspy.Rectangle([box_x1,box_y1], [box_x2,box_y2])
        
        #center of first resonator (used to draw the second one)
        centre_up = (Ushape_d[-1][0] + t/2 + tv + A/2,box_y2 + r)
        #get coordinates for second resonator
        Ushape_u, Mshape_u, B = OPKI_up(L,s,w*1e-6,A*1e-6,t*1e-6,centre=centre_up, compact = True)
    else:
        #get coordinates for first resonator
        Ushape_d, Mshape_d, B, centre_d = OPKI_down(L,s,w*1e-6,A*1e-6,t*1e-6,return_center=(True))
        
        #corners of 1st ground patch: down
        box_x1, box_y1 = Mshape_d[-1][0] - e/2 - w/2, Mshape_d[-1][1]
        box_x2, box_y2 = Mshape_d[-1][0] + w/2 + e/2, Mshape_d[-1][1]-f
        ground_box1 = gdspy.Rectangle([box_x1,box_y1], [box_x2,box_y2])
        
        #center of first resonator (used to draw the second one)
        centre_up = (Ushape_d[-1][0] + t/2 + tv + A/2,box_y2 + r)
        #get coordinates for second resonator
        Ushape_u, Mshape_u,B = OPKI_up(L,s,w*1e-6,A*1e-6,t*1e-6,centre=centre_up)
   
    #draw the two resonator
    U_down = gdspy.FlexPath(Ushape_d, t, **carac)
    M_down = gdspy.FlexPath(Mshape_d,w,**carac)  
    U_up = gdspy.FlexPath(Ushape_u, t, **carac)
    M_up = gdspy.FlexPath(Mshape_u, w, **carac)    
    
    #corners of 2nd ground patch: up
    box_x3, box_y3 = Mshape_u[-1][0] - e/2 - w/2, Mshape_u[-1][1]
    box_x4, box_y4 = Mshape_u[-1][0] + w/2 + e/2, Mshape_u[-1][1] + f
    ground_box2 = gdspy.Rectangle([box_x3,box_y3], [box_x4,box_y4])
    
    #corners of 1st ground strip in between resonators
    strip_x1, strip_y1 = Ushape_d[-1][0] + t/2 + ts, box_y2
    strip_x2, strip_y2 = strip_x1 + xv, Ushape_d[-2][1] + t/2 + r
    #corners of 2nd ground strip in between resonators
    strip_x3, strip_y3 = Ushape_u[-1][0] + t/2 + ts, box_y4
    strip_x4, strip_y4 = strip_x3 + xw, Ushape_u[-2][1] - t/2 - r
    
    
    if tv < tw:
        print('trivial')
    if tv == tw:
        print('normal')
    if tv > tw:
        print('topological')
    
    unitcell = gdspy.Cell('unit cell')
    ground = gdspy.Cell('ground')
    unitcell.add(U_down)
    unitcell.add(M_down)
    ground.add(ground_box1)
    unitcell.add(U_up)
    unitcell.add(M_up)
    ground.add(ground_box2)
    
    if ground_in_between==True:
        ground_strip1 = gdspy.Rectangle([strip_x1,strip_y1], [strip_x2,strip_y2])
        ground.add(ground_strip1)
        ground_strip2 = gdspy.Rectangle([strip_x3,strip_y3], [strip_x4,strip_y4])
        ground.add(ground_strip2)
        
    startM = Mshape_d[-1]
    startU = Ushape_d[-1]
    stopUy = Ushape_u[-1][1]
    strip_height = np.abs(strip_y2 - strip_y1)

    return unitcell, ground, startM, startU, stopUy, strip_height, B, centre_d, centre_up

def unit_cell_upup(L,s,w,A,t,tv,tw,e=20.5e-6,f=24e-6,r=29e-6,carac = {'layer' : 0, 'datatype' : 3},gamma=1/4,ground_in_between=True,compactRes=False):
    A = A*1e6
    t = t*1e6
    w = w*1e6
    e = e*1e6
    f = f*1e6
    r = r*1e6
    tv = tv*1e6
    tw = tw*1e6
    xv = 0.5*(tv - gamma*min(tv,tw))
    xw = 0.5*(tw - gamma*min(tv,tw))
    # print(xv,xw)
    

    if compactRes == True:
        #get coordinates for first resonator
        Ushape_1, Mshape_1, B, centre_1 = OPKI_up(L,s,w*1e-6,A*1e-6,t*1e-6,return_center=(True), compact=True)
        
        #corners of 1st ground patch: up
        box_x1, box_y1 = Mshape_1[-1][0] - e/2 - w/2, Mshape_1[-1][1]
        box_x2, box_y2 = Mshape_1[-1][0] + w/2 + e/2, Mshape_1[-1][1]+f
        ground_box1 = gdspy.Rectangle([box_x1,box_y1], [box_x2,box_y2])
        
        #center of resonator (used to draw the second one)
        centre_2 = (Ushape_1[-1][0] + t/2 + tv + A/2,centre_1[1])
        #get coordinates for second resonator
        Ushape_2, Mshape_2, B = OPKI_up(L,s,w*1e-6,A*1e-6,t*1e-6,centre=centre_2, compact = True)
    else:
        #get coordinates for first resonator
        Ushape_1, Mshape_1, B, centre_1 = OPKI_up(L,s,w*1e-6,A*1e-6,t*1e-6,return_center=(True))
        
        #corners of 1st ground patch: down
        box_x1, box_y1 = Mshape_1[-1][0] - e/2 - w/2, Mshape_1[-1][1]
        box_x2, box_y2 = Mshape_1[-1][0] + w/2 + e/2, Mshape_1[-1][1]+f
        ground_box1 = gdspy.Rectangle([box_x1,box_y1], [box_x2,box_y2])
        
        #center of resonator (used to draw the second one)
        centre_2 = (Ushape_1[-1][0] + t/2 + tv + A/2,centre_1[1])
        #get coordinates for second resonator
        Ushape_2, Mshape_2 = OPKI_up(L,s,w*1e-6,A*1e-6,t*1e-6,centre=centre_2)
   
    #draw the two resonator
    U_1 = gdspy.FlexPath(Ushape_1, t, **carac)
    M_1 = gdspy.FlexPath(Mshape_1,w,**carac)  
    U_2 = gdspy.FlexPath(Ushape_2, t, **carac)
    M_2 = gdspy.FlexPath(Mshape_2, w, **carac)    
    
    #corners of 2nd ground patch: up
    box_x3, box_y3 = Mshape_2[-1][0] - e/2 - w/2, Mshape_2[-1][1]
    box_x4, box_y4 = Mshape_2[-1][0] + w/2 + e/2, Mshape_2[-1][1] + f
    ground_box2 = gdspy.Rectangle([box_x3,box_y3], [box_x4,box_y4])
    
    #corners of 1st ground strip in between resonators
    strip_x1, strip_y1 = Ushape_1[-1][0] + t/2 + xv, box_y2
    strip_x2, strip_y2 = strip_x1 + gamma*min(tv,tw), Ushape_1[-2][1] - t/2 - r
    #corners of 2nd ground strip in between resonators
    strip_x3, strip_y3 = Ushape_2[-1][0] + t/2 + xw, box_y4
    strip_x4, strip_y4 = strip_x3 + gamma*min(tv,tw), Ushape_2[-2][1] - t/2 - r
    
    
    if tv < tw:
        print('trivial')
    if tv == tw:
        print('normal')
    if tv > tw:
        print('topological')
    
    unitcell = gdspy.Cell('unit cell')
    ground = gdspy.Cell('ground')
    unitcell.add(U_1)
    unitcell.add(M_1)
    ground.add(ground_box1)
    unitcell.add(U_2)
    unitcell.add(M_2)
    ground.add(ground_box2)
    
    if ground_in_between==True:
        ground_strip1 = gdspy.Rectangle([strip_x1,strip_y1], [strip_x2,strip_y2])
        ground.add(ground_strip1)
        ground_strip2 = gdspy.Rectangle([strip_x3,strip_y3], [strip_x4,strip_y4])
        ground.add(ground_strip2)
        
    startM = Mshape_1[-1]
    startU = Ushape_1[-1]
    stopUy = Ushape_2[-1][1]
    strip_height = np.abs(strip_y2 - strip_y1)

    return unitcell, ground, startM, startU, stopUy, strip_height, B, centre_1, centre_2


def waveguide(Q, unitcell_size,startM, startU, stopU, strip_height, tw, tv, N_ghost, t = 2e-6, Sr = [1e-6,2e-6], Sg = 60e-6,T=100e-6,R=200e-6,f=24e-6, A = 50e-6):
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
        print('extent ghosts is: ',extent_ghosts, ' (rN_ghost odd)')
    else:
        extent_ghosts = tw + N_ghost*A + N_ghost/2*tw + N_ghost/2*tv
        print('extent ghosts is: ',extent_ghosts, ' (N_ghost even)')
    
    #(x1,y1) and (x2,y2): corners of outer box, (x3,y3) and (x4,y4): corners of window
    x1, y1 = startM[0] - A/2 - Sg - R - extent_ghosts, startM[1] - f - T
    x2, y2 = x1 + 2*(R + Sg) + Q*unitcell_size[0]-tw + 2*extent_ghosts, y1 + 2*T + strip_height
    x3, y3 = startM[0] - A/2 - Sg - extent_ghosts, startM[1] - f
    x4, y4 = x3 + 2 * Sg + Q*unitcell_size[0]-tw + 2*extent_ghosts, y3 + strip_height
    
    #draw outer box + window
    big_plane = gdspy.Rectangle([x1,y1], [x2,y2])
    window = gdspy.Rectangle([x3,y3], [x4,y4])
    #subtract window from outer box
    ground_plane = gdspy.boolean(big_plane,window,'not')
    #draw path guiding the feedline
    points = [[x1 , y1 + 2/3*(2*T + strip_height)],[x1 + R/4, y1 + 2/3*(2*T + strip_height)],[x1 + R/2, y1 + 3/5*T],[startU[0] - A + t + Sr[0],y1 + 1/5*T],[startU[0] - A + t + Sr[0],startU[1]-Sr[1]]]
    if Sr[0] > 0:
        points[-2][0] = startU[0] - A + t + Sr[0] + w_end
        points[-1][0] = startU[0] - A + t + Sr[0] + w_end
    width_feedline = [w_start, w_start - 1/4*(w_start - w_end), w_end + 1/4*(w_start-w_end), w_end]
    width_guide = [11/5*w for w in width_feedline]
    bend_radii = [3*w for w in width_guide]
    guide = gdspy.FlexPath([points[0],points[1]],width_guide[0],corners='circular bend',bend_radius=bend_radii[0])
    guide = guide.segment(points[2],width_guide[1]).segment(points[3],width_guide[2]).segment(points[4],width_guide[3]) #
    ground_plane = gdspy.boolean(ground_plane,guide,'not')
    feedline = gdspy.FlexPath([points[0],points[1]],width_feedline[0],corners='circular bend',bend_radius=bend_radii[0])
    feedline = feedline.segment(points[2],width_feedline[1]).segment(points[3],width_feedline[2]).segment(points[4],width_feedline[3]) #
    
    points2 = [[x2, y2 - 2/3*(2*T + strip_height)],[x2 - R/4, y2 - 2/3*(2*T + strip_height)],[x2 - R/2, y2 - 3/5*T],[stopU[0] - Sr[0],y2 - 1/5*T],[stopU[0] - Sr[0],stopU[1]+Sr[1]]]
    if Sr[0] > 0:
        points2[-2][0] = stopU[0] - Sr[0] - w_end
        points2[-1][0] = stopU[0] - Sr[0] - w_end
    guide2 = gdspy.FlexPath([points2[0],points2[1]],width_guide[0],corners='circular bend',bend_radius=bend_radii[0])
    guide2 = guide2.segment(points2[2],width_guide[1]).segment(points2[3],width_guide[2]).segment(points2[4],width_guide[3])
    ground_plane2 = gdspy.boolean(ground_plane,guide2,'not')
    feedline2 = gdspy.FlexPath([points2[0],points2[1]],width_feedline[0],corners='circular bend',bend_radius=bend_radii[0])
    feedline2 = feedline2.segment(points2[2],width_feedline[1]).segment(points2[3],width_feedline[2]).segment(points2[4],width_feedline[3])
    
    ground_plane_w_feed = gdspy.boolean(ground_plane2,feedline,'or')
    ground_plane_w_feed = gdspy.boolean(ground_plane_w_feed, feedline2, 'or')
    
    # print(points2)
    # print(width_feedline)
    
    return ground_plane_w_feed#,points2

def waveguide_extended(Q, unitcell_size,startM, startU, stopU, strip_height, tw, tv, N_ghost, t = 2e-6, Sr = [1e-6,2e-6], Sg = 60e-6,T=100e-6,R=200e-6,f=24e-6, A = 50e-6):
    tw = tw*1e6
    tv = tv*1e6
    t = t*1e6
    Sr = Sr[0]*1e6,Sr[1]*1e6 #spacing feedline to resonator [x,y]-direction
    Sg = Sg*1e6 #lateral spacing to ground planes
    T = T*1e6 #thickness of ground planes
    R = R*1e6 #width of ground planes
    E = 800
    f = f*1e6
    A = A*1e6
    w_patch = 360.0
    w_start = 10.0 #start width feedline
    w_end = 2.0 #end width feedline
    
    #calculate length ghosts demand
    if (N_ghost-1)%2==0:
        extent_ghosts = tw + tv + N_ghost*A + (N_ghost-1)/2*(tw + tv)
        print('extent ghosts is: ',extent_ghosts, ' (rN_ghost odd)')
    else:
        extent_ghosts = tw + N_ghost*A + N_ghost/2*tw + N_ghost/2*tv
        print('extent ghosts is: ',extent_ghosts, ' (N_ghost even)')
        
    #(x1,y1) and (x2,y2): corners of outer box, (x3,y3) and (x4,y4): corners of window
    x1, y1 = startM[0] - A/2 - Sg - R - extent_ghosts, startM[1] - f - T
    x2, y2 = x1 + 2*(R + Sg) + Q*unitcell_size[0]-tw + 2*extent_ghosts, y1 + 2*T + strip_height
    x3, y3 = startM[0] - A/2 - Sg - extent_ghosts, startM[1] - f
    x4, y4 = x3 + 2 * Sg + Q*unitcell_size[0]-tw + 2*extent_ghosts, y3 + strip_height
    
    
    #draw outer box + window
    big_plane = gdspy.Rectangle([x1,y1], [x2,y2])
    window = gdspy.Rectangle([x3,y3], [x4,y4])
    #subtract window from outer box
    ground_plane = gdspy.boolean(big_plane,window,'not')
    #coordinates/dimensions of right feedline + "feedline-window"
    points = [[x1-E,y1 + 1/2*(2*T + strip_height)],[x1-3*E/4, y1 + 1/2*(2*T + strip_height)],[x1-E/4, y1 + 1/2*(2*T + strip_height)],[x1-R/3, y1 + 1/2*(2*T + strip_height)],[x1, y1 + 1/2*(2*T + strip_height)],[x1 + R/4, y1 + 1/2*(2*T + strip_height)],[x1 + R/2, y1 + 3/5*T], [startU[0] - A + t + Sr[0],y1 + 1/5*T],[startU[0] - A + t + Sr[0],startU[1]-Sr[1]]]
    if Sr[0] > 0:
        points[-2][0] = startU[0] - A + t + Sr[0] + w_end
        points[-1][0] = startU[0] - A + t + Sr[0] + w_end
    width_feedline = [w_patch,(w_patch-w_start)/8,w_start, w_start - 1/4*(w_start - w_end), w_end + 1/4*(w_start-w_end), w_end]
    width_guide = [11/5*w for w in width_feedline]
    width_guide[0] = w_patch + 10*w_start
    bend_radii = [3*w for w in width_guide]
    #such that patch is enclosed in "feedline-window"
    patch_start = [points[0][0]-(width_guide[0]-width_feedline[0])/2,points[0][1]]
    
    #draw the super big ground box (contact pads included)
    x5, y5 = x1 - E, y1 + 1/2*(2*T + strip_height) - 3*width_guide[0]/4
    x6, y6 = x2 + E, y2 - 1/2*(2*T + strip_height) + 3*width_guide[0]/4
    super_big_box = gdspy.Rectangle([x5,y5],[x6,y6])
    
    ground_plane = gdspy.boolean(ground_plane, super_big_box, 'or')
    
    #draw the right "feedline-window"
    big_patch = gdspy.FlexPath([patch_start,points[1]], width_guide[0],corners='circular bend',bend_radius = bend_radii[1])
    big_patch = big_patch.segment(points[2],width_guide[1]).segment(points[3],width_guide[2])
    guide = gdspy.FlexPath([points[3],points[4]],width_guide[2],corners='circular bend',bend_radius=bend_radii[2])
    guide = guide.segment(points[6],width_guide[3]).segment(points[7],width_guide[4]).segment(points[8],width_guide[5])
    
    guide_join = gdspy.boolean(big_patch,guide,'or')
    
    ground_plane = gdspy.boolean(ground_plane,guide_join,'not')
    
    patch = gdspy.FlexPath([points[0],points[1]], width_feedline[0],corners='circular bend',bend_radius = bend_radii[1])
    patch = patch.segment(points[2],width_feedline[1]).segment(points[3],width_feedline[2])
    feedline = gdspy.FlexPath([points[3],points[4]],width_feedline[2],corners='circular bend',bend_radius=bend_radii[2])
    feedline = feedline.segment(points[6],width_feedline[2]).segment(points[7],width_feedline[4]).segment(points[8],width_feedline[5])
    
    feedline_join = gdspy.boolean(patch,feedline,'or')
    
    points2 = [[x2+E,y2 - 1/2*(2*T + strip_height)],[x2+3*E/4, y2 - 1/2*(2*T + strip_height)],[x2+E/4, y2 - 1/2*(2*T + strip_height)],[x2+R/3, y2 - 1/2*(2*T + strip_height)],[x2, y2 - 1/2*(2*T + strip_height)],[x2 - R/4, y2 - 1/2*(2*T + strip_height)],[x2 - R/2, y2 - 3/5*T],[stopU[0] - Sr[0],y2 - 1/5*T],[stopU[0] - Sr[0],stopU[1]+Sr[1]]]
    if Sr[0] > 0:
        points2[-2][0] = stopU[0] - Sr[0] - w_end
        points2[-1][0] = stopU[0] - Sr[0] - w_end
    patch_start2 = [points2[0][0]+(width_guide[0]-width_feedline[0])/2,points2[0][1]]
    print('points: ' , points)
    
    big_patch2 = gdspy.FlexPath([patch_start2,points2[1]], width_guide[0],corners='circular bend',bend_radius = bend_radii[1])
    big_patch2 = big_patch2.segment(points2[2],width_guide[1]).segment(points2[3],width_guide[2])
    guide2 = gdspy.FlexPath([points2[3],points2[4]],width_guide[2],corners='circular bend',bend_radius=bend_radii[2])
    guide2 = guide2.segment(points2[6],width_guide[3]).segment(points2[7],width_guide[4]).segment(points2[8],width_guide[5])
    #.segment(points2[5],width_guide[2])
    guide_join2 = gdspy.boolean(big_patch2,guide2,'or')
    
    ground_plane2 = gdspy.boolean(ground_plane,guide_join2,'not')
    
    patch2 = gdspy.FlexPath([points2[0],points2[1]], width_feedline[0],corners='circular bend',bend_radius = bend_radii[1])
    patch2 = patch2.segment(points2[2],width_feedline[1]).segment(points2[3],width_feedline[2])
    
    feedline2 = gdspy.FlexPath([points2[3],points2[4]],width_feedline[2],corners='circular bend',bend_radius=bend_radii[2])
    feedline2 = feedline2.segment(points2[6],width_feedline[3]).segment(points2[7],width_feedline[4]).segment(points2[8],width_feedline[5])
    #.segment(points2[5],width_feedline[2])
    feedline_join2 = gdspy.boolean(patch2,feedline2,'or')
    
    
    
    ground_plane_w_feed = gdspy.boolean(ground_plane2,feedline_join,'or')
    ground_plane_w_feed = gdspy.boolean(ground_plane_w_feed, feedline_join2, 'or')
    
    
    
    return ground_plane_w_feed

def waveguide_simulation(Q, unitcell_size,startM, startU, stopU, strip_height, tw, tv, N_ghost, t = 2e-6, Sr = [1e-6,2e-6], Sg = 60e-6,T=100e-6,R=200e-6,f=24e-6, A = 50e-6,w_end = 2e-6,cozy = False):
    tw = tw*1e6
    tv = tv*1e6
    t = t*1e6
    Sr = Sr[0]*1e6,Sr[1]*1e6 #spacing feedline to resonator [x,y]-direction
    Sg = Sg*1e6 #lateral spacing to ground planes
    T = T*1e6 #thickness of ground planes
    R = R*1e6 #width of ground planes
    f = f*1e6
    A = A*1e6
    w_end = w_end*1e6 #end width feedline
    
     #calculate length ghosts demand
    if (N_ghost-1)%2==0:
        extent_ghosts = tw + tv + N_ghost*A + (N_ghost-1)/2*(tw + tv)
        print('extent ghosts is: ',extent_ghosts, ' (rN_ghost odd)')
    else:
        extent_ghosts = tw + N_ghost*A + N_ghost/2*tw + N_ghost/2*tv
        print('extent ghosts is: ',extent_ghosts, ' (N_ghost even)')
    
    
    if cozy == False:
        #(x1,y1) and (x2,y2): corners of outer box, (x3,y3) and (x4,y4): corners of window
        x1, y1 = startM[0] - A/2 - Sg - R - extent_ghosts, startM[1] - f - T
        x2, y2 = x1 + 2*(R + Sg) + Q*unitcell_size[0]-tw + 2*extent_ghosts, y1 + 2*T + strip_height
        x3, y3 = startM[0] - A/2 - Sg - extent_ghosts, startM[1] - f
        x4, y4 = x3 + 2 * Sg + Q*unitcell_size[0]-tw + 2*extent_ghosts, y3 + strip_height
    else:
        #(x1,y1) and (x2,y2): corners of outer box, (x3,y3) and (x4,y4): corners of window
        x1, y1 = startM[0] - A/2 - Sg - R - extent_ghosts, startM[1] - T
        x2, y2 = x1 + 2*(R + Sg) + Q*unitcell_size[0]-tw + 2*extent_ghosts, y1 + 2*T + strip_height - 2*f
        x3, y3 = startM[0] - A/2 - Sg - extent_ghosts, startM[1]
        x4, y4 = x3 + 2 * Sg + Q*unitcell_size[0]-tw + 2*extent_ghosts, y3 + strip_height - 2*f
    
    #draw outer box + window
    big_plane = gdspy.Rectangle([x1,y1], [x2,y2])
    window = gdspy.Rectangle([x3,y3], [x4,y4])
    #subtract window from outer box
    ground_plane = gdspy.boolean(big_plane,window,'not')
    #draw path guiding the feedline
    points = [[startU[0] - A + t + Sr[0],y1],[startU[0] - A + t + Sr[0],startU[1]-Sr[1]]]
    if Sr[0] > 0:
        points[-2][0] = startU[0] - A + t + Sr[0] + w_end
        points[-1][0] = startU[0] - A + t + Sr[0] + w_end
    guide = gdspy.FlexPath(points,11/5*w_end)
    ground_plane = gdspy.boolean(ground_plane,guide,'not')
    feedline = gdspy.FlexPath(points,w_end)
    
    points2 = [[stopU[0] - Sr[0],y2],[stopU[0] - Sr[0],stopU[1]+Sr[1]]]
    if Sr[0] > 0:
        points2[-2][0] = stopU[0] - Sr[0] - w_end
        points2[-1][0] = stopU[0] - Sr[0] - w_end
    guide2 = gdspy.FlexPath(points2,11/5*w_end)
    ground_plane2 = gdspy.boolean(ground_plane,guide2,'not')
    feedline2 = gdspy.FlexPath(points2,w_end)
    
    ground_plane_w_feed = gdspy.boolean(ground_plane2,feedline,'or')
    ground_plane_w_feed = gdspy.boolean(ground_plane_w_feed, feedline2, 'or')
    
    
    return ground_plane_w_feed

def only_waveguide(startM, startU, strip_height, t = 2e-6, Sg = 60e-6,T=100e-6,R=200e-6,f=24e-6, A = 50e-6):
    t = t*1e6
    Sg = Sg*1e6 #lateral spacing to ground planes
    T = T*1e6 #thickness of ground planes
    R = R*1e6 #width of ground planes
    f = f*1e6
    A = A*1e6
    w_start = 10 #start width feedline
    w_end = 2 #end width feedline
    
    #(x1,y1) and (x2,y2): corners of outer box, (x3,y3) and (x4,y4): corners of window
    x1, y1 = startM[0] - A/2 - Sg - R, startM[1] - f - T
    x2, y2 = x1 + 2*(R + Sg), y1 + 2*T + strip_height
    # x3, y3 = startM[0] - A/2 - Sg, startM[1] - f
    # x4, y4 = x3 + 2 * Sg, y3 + strip_height
    
    #draw outer box
    big_plane = gdspy.Rectangle([x1,y1], [x2,y2])

    #draw path guiding the feedline
    points = [[x1 , y1 + 2/3*(2*T + strip_height)],[x1 + R/4, y1 + 2/3*(2*T + strip_height)],[x1 + R/2, y1 + 3/5*T],[startU[0] - A + t,y1 + 1/5*T],[startU[0] - A + t,startU[1]],[startU[0] - A + t,y2 - 1/5*T],[x2 - R/2, y2 - 3/5*T],[x2 - R/4, y2 - 2/3*(2*T + strip_height)],[x2, y2 - 2/3*(2*T + strip_height)]]
    width_feedline = [w_start, w_start - 1/4*(w_start - w_end), w_end + 1/4*(w_start-w_end), w_end]
    width_guide = [11/5*w for w in width_feedline]
    bend_radii = [3*w for w in width_guide]
    
    guide = gdspy.FlexPath([points[0],points[1]],width_guide[0],corners='circular bend',bend_radius=bend_radii[0])
    guide = guide.segment(points[2],width_guide[1]).segment(points[3],width_guide[2]).segment(points[4],width_guide[3]) #
    guide = guide.segment(points[5],width_guide[3]).segment(points[6],width_guide[2]).segment(points[7],width_guide[1]).segment(points[8],width_guide[0])
   
    feedline = gdspy.FlexPath([points[0],points[1]],width_feedline[0],corners='circular bend',bend_radius=bend_radii[0])
    feedline = feedline.segment(points[2],width_feedline[1]).segment(points[3],width_feedline[2]).segment(points[4],width_feedline[3]) #
    feedline = feedline.segment(points[5],width_feedline[3]).segment(points[6],width_feedline[2]).segment(points[7],width_feedline[1]).segment(points[8],width_feedline[0])
    
    ground_plane = gdspy.boolean(big_plane,guide,'not')         

    ground_plane_w_feed = gdspy.boolean(ground_plane,feedline,'or')
    
    # print(points2)
    # print(width_feedline)
    
    return ground_plane_w_feed#,points2

def T_feedline_extended(Q, unitcell_size,startM, startU, stopU, strip_height, tw, tv, N_ghost, t = 2e-6, Sr = [1e-6,2e-6], D = 500e-6, Sg = 60e-6,T=100e-6,R=200e-6,f=24e-6, A = 50e-6, B = 60e-6,laserwriter=False):
    tw = tw*1e6
    tv = tv*1e6
    t = t*1e6
    Sr = Sr[0]*1e6,Sr[1]*1e6 #spacing feedline to resonator [x,y]-direction
    D = D*1e6  #spacing to array region
    Sg = Sg*1e6 #lateral spacing to ground planes
    T = T*1e6 #thickness of ground planes
    R = R*1e6 #width of ground planes
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
    x1, y1 = x3 - R, y3 - T
    x2, y2 = x4 + R, y1 + 2*T + strip_height
    
    width_feedline = [w_patch,w_middle,w_start]
    width_guide = [73/45*w for w in width_feedline]
    if laserwriter == False:
        width_guide[1] = 31/30*width_feedline[1]
        width_guide[2] = 11/5 * width_feedline[2]
    else:
        width_guide[1] = 16/15*width_feedline[1]
        width_guide[2] = 11/5 * width_feedline[2]
        
    dyke = [(width_guide[0]-width_feedline[0])/2,(width_guide[1]-width_feedline[1])/2,(width_guide[2]-width_feedline[2])/2]
    
    
    E = 0.5*(chip_length - (x2-x1)) - (width_guide[0]-width_feedline[0])/2
    
    #draw outer box + window
    # big_plane = gdspy.Rectangle([x1,y1], [x2,y2])
    window = gdspy.Rectangle([x3,y3], [x4,y4])
    
    
    #draw the super big ground box (contact pads included)
    Y0 = y1 + 1/2*(2*T + strip_height)
    x5, y5 = x1 - E, Y0 - 3*width_guide[0]/4
    x6, y6 = x2 + E, Y0 + 3*width_guide[0]/4
    super_big_box = gdspy.Rectangle([x5,y5],[x6,y6])
    #subtract window from outer box
    ground_plane = gdspy.boolean(super_big_box, window, 'not')
    
    #coordinates/ dimensions for FlexPath: left feedline + guide (etched away)
    points = [[x1-E,Y0],[x1-7*E/8, Y0],[x1-3*E/4, Y0],[x1,Y0],[x3 - D, Y0],[startU[0] + t/2 - A - Sr[0] - w_start/2, Y0]]
    
    #enclose contact pad with etched line (left edge)
    patch_start = [points[0][0]- dyke[0],points[0][1]]
    
    #enclosing contact pad/patch to define it/ guide will be etched away
    guide = gdspy.FlexPath([patch_start, points[1]], width_guide[0]).segment(points[2],width_guide[1]).segment(points[4],width_guide[1]).segment(points[5],width_guide[2])
    T_bar_guide = gdspy.Rectangle([points[5][0] - w_start/2 - dyke[2], startU[1] + B + dyke[2]], [points[5][0] + w_start/2 + dyke[2], startU[1] - dyke[2]])
    
    guide_T = gdspy.boolean(guide, T_bar_guide, 'or')
    ground_plane = gdspy.boolean(ground_plane,guide_T,'not')
    
    #draw the left feedline (not etched away)
    feedline = gdspy.FlexPath([points[0],points[1]], width_feedline[0]).segment(points[2],width_feedline[1]).segment(points[4],width_feedline[1]).segment(points[5],width_feedline[2])
    T_bar = gdspy.Rectangle([points[5][0] - w_start/2, startU[1] + B], [points[5][0] + w_start/2, startU[1]])
    
    feedline = gdspy.boolean(feedline, T_bar, 'or')
    
    #coordinates for right feedline + guide (FlexPath)
    Y02 = y2 - 1/2*(2*T + strip_height)
    points2 = [[x2+E,Y02],[x2+7*E/8, Y02],[x2+3*E/4, Y02],[x2,Y02],[x4 + D, Y02],[stopU[0] + t/2 + Sr[0] + w_start/2,Y02]]
    
    #enclose contact pad with etched line (right edge)
    patch_start2 = [points2[0][0]+dyke[0],points2[0][1]]
    
    #draw right guide that will be etched away (leaving the feedline inside the ground plane)
    guide2 = gdspy.FlexPath([patch_start2, points2[1]], width_guide[0]).segment(points2[2],width_guide[1]).segment(points2[4],width_guide[1]).segment(points2[5],width_guide[2])
    T_bar_guide2 = gdspy.Rectangle([points2[5][0] - w_start/2-dyke[2], stopU[1]+dyke[2]], [points2[5][0] - w_start/2 - dyke[2], stopU[1]-B-dyke[2]])
    
    guide_T2 = gdspy.boolean(guide2, T_bar_guide2, 'or')
    ground_plane2 = gdspy.boolean(ground_plane,guide_T2,'not')
    
    #draw the right feedline
    feedline2 = gdspy.FlexPath([points2[0],points2[1]],width_feedline[0]).segment(points2[2],width_feedline[1]).segment(points2[4],width_feedline[1]).segment(points2[5],width_feedline[2])
    T_bar2 = gdspy.Rectangle([points2[5][0] - w_start/2, stopU[1]], [points2[5][0] + w_start/2, stopU[1]-B])

    feedline2 = gdspy.boolean(feedline2, T_bar2, 'or')
    
    ground_plane_w_feed = gdspy.boolean(ground_plane2,feedline,'or')
    ground_plane_w_feed = gdspy.boolean(ground_plane_w_feed, feedline2, 'or')
    
    return ground_plane_w_feed

def T_feedline_simulation(Q, unitcell_size,startM, startU, stopU, strip_height, tw, tv, t = 2e-6, Sr = [1e-6,2e-6], D = 500e-6, Sg = 60e-6,T=100e-6,R=200e-6,f=24e-6, A = 50e-6, B = 60e-6, cozy = True, w = [360e-6,90e-6,20e-6]):
    tw = tw*1e6
    tv = tv*1e6
    t = t*1e6
    Sr = Sr[0]*1e6,Sr[1]*1e6 #spacing feedline to resonator [x,y]-direction
    D = D*1e6  #spacing to array region
    Sg = Sg*1e6 #lateral spacing to ground planes
    T = T*1e6 #thickness of ground planes
    R = R*1e6 #width of ground planes
    f = f*1e6
    A = A*1e6
    B = B*1e6
    w_patch = w[0]*1e6
    w_middle = w[1]*1e6
    w_start = w[2]*1e6 #start width feedline
    
    
    if cozy == False:
        #(x1,y1) and (x2,y2): corners of outer box, (x3,y3) and (x4,y4): corners of window
        x1, y1 = startM[0] - A/2 - Sg - R, startM[1] - f - T
        x2, y2 = x1 + 2*(R + Sg) + Q*unitcell_size[0]-tw, y1 + 2*T + strip_height
        x3, y3 = startM[0] - A/2 - Sg, startM[1] - f
        x4, y4 = x3 + 2 * Sg + Q*unitcell_size[0]-tw, y3 + strip_height
    else:
        #(x1,y1) and (x2,y2): corners of outer box, (x3,y3) and (x4,y4): corners of window
        x1, y1 = startM[0] - A/2 - Sg - R, startM[1] - T
        x2, y2 = x1 + 2*(R + Sg) + Q*unitcell_size[0]-tw, y1 + 2*T + strip_height - 2*f
        x3, y3 = startM[0] - A/2 - Sg, startM[1]
        x4, y4 = x3 + 2 * Sg + Q*unitcell_size[0]-tw, y3 + strip_height - 2*f
    
    
    width_feedline = [w_patch,w_middle,w_start]
    width_guide = [73/45*w for w in width_feedline]
    width_guide[1] = 31/30*width_feedline[1]
    width_guide[2] = 11/5 * width_feedline[2]
    dyke = [(width_guide[0]-width_feedline[0])/2,(width_guide[1]-width_feedline[1])/2,(width_guide[2]-width_feedline[2])/2]

    
    #draw outer box + window
    big_plane = gdspy.Rectangle([x1,y1], [x2,y2])
    window = gdspy.Rectangle([x3,y3], [x4,y4])
    
    #subtract window from outer box
    ground_plane = gdspy.boolean(big_plane, window, 'not')
    
    #coordinates/ dimensions for FlexPath: left feedline + guide (etched away)
    Y0 = y1 + 1/2*(2*T + strip_height)
    points = [[x1,Y0],[startU[0] + t/2 - A - Sr[0] - w_start/2, Y0]]
    
    #enclosing contact pad/patch to define it/ guide will be etched away
    guide = gdspy.FlexPath([points[0], points[1]], width_guide[2])
    # T_bar_guide = gdspy.Rectangle([points[1][0] - w_start/2 - dyke[2], startU[1] + B + dyke[2]], [points[1][0] + w_start/2 + dyke[2], startU[1] - dyke[2]])
    
    # guide_T = gdspy.boolean(guide, T_bar_guide, 'or')
    ground_plane = gdspy.boolean(ground_plane,guide,'not')
    
    #draw the left feedline (not etched away)
    feedline = gdspy.FlexPath([points[0],points[1]], width_feedline[2])
    T_bar = gdspy.Rectangle([points[1][0] - w_start/2, startU[1] + B], [points[1][0] + w_start/2, startU[1]])
    
    feedline = gdspy.boolean(feedline, T_bar, 'or')
    
    #coordinates for right feedline + guide (FlexPath)
    Y02 = y2 - 1/2*(2*T + strip_height)
    points2 = [[x2,Y02],[stopU[0] + t/2 + Sr[0] + w_start/2,Y02]]
    
    
    #draw right guide that will be etched away (leaving the feedline inside the ground plane)
    guide2 = gdspy.FlexPath([points2[0], points2[1]], width_guide[2])
    # T_bar_guide2 = gdspy.Rectangle([points2[1][0] - w_start/2-dyke[2], stopU[1]+dyke[2]], [points2[1][0] - w_start/2 - dyke[2], stopU[1]-B-dyke[2]])
    
    # guide_T2 = gdspy.boolean(guide2, T_bar_guide2, 'or')
    ground_plane2 = gdspy.boolean(ground_plane,guide2,'not')
    
    #draw the right feedline
    feedline2 = gdspy.FlexPath([points2[0],points2[1]],width_feedline[2])
    T_bar2 = gdspy.Rectangle([points2[1][0] - w_start/2, stopU[1]], [points2[1][0] + w_start/2, stopU[1]-B])

    feedline2 = gdspy.boolean(feedline2, T_bar2, 'or')
    
    ground_plane_w_feed = gdspy.boolean(ground_plane2,feedline,'or')
    ground_plane_w_feed = gdspy.boolean(ground_plane_w_feed, feedline2, 'or')
    
    return ground_plane_w_feed

def waveguide_extended_new(Q, unitcell_size,startM, startU, stopU, strip_height, tw, tv, N_ghost, t = 2e-6, Sr = [1e-6,2e-6], Sg = 60e-6,T=100e-6,R=200e-6,f=24e-6, A = 50e-6, w = [360e-6,90e-6,10e-6,2e-6],cozy=False):
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
    w_patch = w[0]*1e6
    w_core = w[1]*1e6
    w_start = w[2]*1e6 #start width feedline
    w_end = w[3]*1e6 #end width feedline
    
    
    #calculate length ghosts demand
    if (N_ghost-1)%2==0:
        extent_ghosts = tw + tv + N_ghost*A + (N_ghost-1)/2*(tw + tv)
        print('extent ghosts is: ',extent_ghosts, ' (rN_ghost odd)')
    else:
        extent_ghosts = tw + N_ghost*A + N_ghost/2*tw + N_ghost/2*tv
        print('extent ghosts is: ',extent_ghosts, ' (N_ghost even)')
    
    if cozy == False:
        #(x1,y1) and (x2,y2): corners of outer box, (x3,y3) and (x4,y4): corners of window
        x1, y1 = startM[0] - A/2 - Sg - R - extent_ghosts, startM[1] - f - T
        x2, y2 = x1 + 2*(R + Sg) + Q*unitcell_size[0]-tw + 2*extent_ghosts, y1 + 2*T + strip_height
        x3, y3 = startM[0] - A/2 - Sg - extent_ghosts, startM[1] - f
        x4, y4 = x3 + 2 * Sg + Q*unitcell_size[0]-tw + 2*extent_ghosts, y3 + strip_height
    else:
        #(x1,y1) and (x2,y2): corners of outer box, (x3,y3) and (x4,y4): corners of window
        x1, y1 = startM[0] - A/2 - Sg - R - extent_ghosts, startM[1] - T
        x2, y2 = x1 + 2*(R + Sg) + Q*unitcell_size[0]-tw + 2*extent_ghosts, y1 + 2*T + strip_height - 2*f
        x3, y3 = startM[0] - A/2 - Sg - extent_ghosts, startM[1]
        x4, y4 = x3 + 2 * Sg + Q*unitcell_size[0]-tw + 2*extent_ghosts, y3 + strip_height - 2*f
    
    
    width_feedline = [w_patch,w_core,w_start, w_start - 1/4*(w_start - w_end), w_end + 1/4*(w_start-w_end), w_end]
    width_guide = [11/5*w for w in width_feedline]
    width_guide[0] = 73/45*width_feedline[0]
    width_guide[1] = 31/30*width_feedline[1]
    bend_radii = [3*w for w in width_guide]
    dyke = [(width_guide[0]-width_feedline[0])/2,(width_guide[1]-width_feedline[1])/2,(width_guide[2]-width_feedline[2])/2]
    
    E = 0.5*(chip_length - (x2-x1)) - (width_guide[0]-width_feedline[0])/2
    
    #draw the super big ground box (contact pads included)
    Y0 = y1 + 1/2*(2*T + strip_height)
    x5, y5 = x1 - E, Y0 - 3*width_guide[0]/4
    x6, y6 = x2 + E, Y0 + 3*width_guide[0]/4
    super_big_box = gdspy.Rectangle([x5,y5],[x6,y6])
    
    #draw outer box + window
    big_plane = gdspy.Rectangle([x1,y1], [x2,y2])
    window = gdspy.Rectangle([x3,y3], [x4,y4])
    
    #subtract window from outer box
    ground_plane = gdspy.boolean(super_big_box, window, 'not')
    
    #coordinates/ dimensions for FlexPath: right feedline + guide (etched away)
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
    guide = guide.segment(points[4], width_guide[2]).segment(points[5],width_guide[3]).segment(points[6],width_guide[4]).segment(points[7],width_guide[5])
    
    valley = gdspy.boolean(big_patch, guide, 'or')
    
    #draw the feedline (not etched away)
    patch = gdspy.FlexPath([points[0],points[1]], width_feedline[0])#,corners='circular bend',bend_radius = bend_radii[1])
    patch = patch.segment(patch_end,width_feedline[1])
    feedline = gdspy.FlexPath([points[2],points[3]],width_feedline[1],corners='circular bend',bend_radius=bend_radii[2])
    feedline = feedline.segment(points[4], width_feedline[2]).segment(points[5],width_feedline[3]).segment(points[6],width_feedline[4]).segment(points[7],width_feedline[5])
    
    line = gdspy.boolean(patch,feedline,'or')
    
    ground_plane_valley = gdspy.boolean(ground_plane,valley, 'not')
    ground_plane_line = gdspy.boolean(ground_plane_valley,line, 'or')
    
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
   
    valley2 = gdspy.boolean(big_patch2, guide2, 'or')
    
    patch2 = gdspy.FlexPath([points2[0],points2[1]], width_feedline[0])#,corners='circular bend',bend_radius = bend_radii[1])
    patch2 = patch2.segment(patch_end2,width_feedline[1])   
    feedline2 = gdspy.FlexPath([points2[2],points2[3]],width_feedline[1],corners='circular bend',bend_radius=bend_radii[2])
    feedline2 = feedline2.segment(points2[4],width_feedline[2]).segment(points2[5],width_feedline[3]).segment(points2[6],width_feedline[4]).segment(points2[7],width_feedline[5])

    line2 = gdspy.boolean(patch2,feedline2,'or')

    ground_plane_valley2 = gdspy.boolean(ground_plane_line,valley2, 'not')
    ground_plane_line2 = gdspy.boolean(ground_plane_valley2,line2, 'or')
    
    # window_left_guide = gdspy.boolean(window, line, 'or')
    # window_both_guides = gdspy.boolean(window_left_guide, line2, 'or')
    # window_left_feedline = gdspy.boolean(window_both_guides,feedline ,'not')
    # window_both_feedlines = gdspy.boolean(window_left_feedline,feedline2,'not')
    # pads = gdspy.boolean(pad,pad2,'or')
    
    # mask = gdspy.Rectangle([x1+R-150,y1-150], [x2-R+150,y2+150])
    # mask_overlap = gdspy.Rectangle([x1+R - 151, y1-151], [x2-R+151,y2+151])
    
    return ground_plane_line2#,pads,mask,mask_overlap

"""
------------------------------------------------------------------------------
left-handed LC array: Etch mask
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
        print('extent ghosts is: ',extent_ghosts, ' (rN_ghost odd)')
    else:
        extent_ghosts = tw + N_ghost*A + N_ghost/2*tw + N_ghost/2*tv
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
        print('extent ghosts is: ',extent_ghosts, ' (rN_ghost odd)')
    else:
        extent_ghosts = tw + N_ghost*A + N_ghost/2*tw + N_ghost/2*tv
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
    if Sr[0] > 0:
        points[-2][0] = startU[0] - A + t + Sr[0] + w_end
        points[-1][0] = startU[0] - A + t + Sr[0] + w_end
    
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
    if Sr[0] > 0:
        points2[-2][0] = stopU[0]+ Sr[0] + w_end
        points2[-1][0] = stopU[0]+ Sr[0] + w_end
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

def waveguide_extended_negative_new(Q, unitcell_size,startM, startU, stopU, strip_height, tw, tv, N_ghost, t = 2e-6, Sr = [1e-6,2e-6], Sg = 60e-6,T=100e-6,R=200e-6,f=24e-6, A = 50e-6, w = [360e-6,90e-6,10e-6,2e-6],maskdim = [150e-6,300e-6],cozy = False,laserwriter=False):
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
    w_patch = w[0]*1e6
    w_core = w[1]*1e6
    w_start = w[2]*1e6 #start width feedline
    w_end = w[3]*1e6 #end width feedline
    maskdim = maskdim[0]*1e6, maskdim[1]*1e6
    d_overlap = 5
    
    #calculate length ghosts demand
    if (N_ghost-1)%2==0:
        extent_ghosts = tw + tv + N_ghost*A + (N_ghost-1)/2*(tw + tv)
        print('extent ghosts is: ',extent_ghosts, ' (rN_ghost odd)')
    else:
        extent_ghosts = tw + N_ghost*A + N_ghost/2*tw + N_ghost/2*tv
        print('extent ghosts is: ',extent_ghosts, ' (N_ghost even)')
    
    
    if cozy == False:
        #(x1,y1) and (x2,y2): corners of outer box, (x3,y3) and (x4,y4): corners of window
        x1, y1 = startM[0] - A/2 - Sg - R - extent_ghosts, startM[1] - f - T
        x2, y2 = x1 + 2*(R + Sg) + Q*unitcell_size[0]-tw + 2*extent_ghosts, y1 + 2*T + strip_height
        x3, y3 = startM[0] - A/2 - Sg - extent_ghosts, startM[1] - f
        x4, y4 = x3 + 2 * Sg + Q*unitcell_size[0]-tw + 2*extent_ghosts, y3 + strip_height
    else:
        #(x1,y1) and (x2,y2): corners of outer box, (x3,y3) and (x4,y4): corners of window
        x1, y1 = startM[0] - A/2 - Sg - R - extent_ghosts, startM[1] - T
        x2, y2 = x1 + 2*(R + Sg) + Q*unitcell_size[0]-tw + 2*extent_ghosts, y1 + 2*T + strip_height - 2*f
        x3, y3 = startM[0] - A/2 - Sg - extent_ghosts, startM[1]
        x4, y4 = x3 + 2 * Sg + Q*unitcell_size[0]-tw + 2*extent_ghosts, y3 + strip_height - 2*f
    
    width_feedline = [w_patch,w_core,w_start, w_start - 1/4*(w_start - w_end), w_end + 1/4*(w_start-w_end), w_end]
    
    if laserwriter == True:
        width_guide = [w+6 for w in width_feedline]
        width_guide[0] = 73/45*width_feedline[0]
    else:
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
    # patch_end = [points[2][0]+1,points[2][1]]
    
    #enclosing contact pad/patch to define it/ this big_patch + guide will be etched away
    big_patch = gdspy.FlexPath([patch_start,points[1]], width_guide[0])#,corners=#'circular bend',bend_radius = bend_radii[1])
    big_patch = big_patch.segment(points[2],width_guide[1])
    guide = gdspy.FlexPath([points[2],points[3]],width_guide[1],corners='circular bend',bend_radius=bend_radii[2])
    guide = guide.segment(points[4], width_guide[2]).segment(points[5],width_guide[3]).segment(points[6],width_guide[4]).segment(points[7],width_guide[5])
    
    #draw the feedline (not etched away)
    patch = gdspy.FlexPath([points[0],points[1]], width_feedline[0])#,corners='circular bend',bend_radius = bend_radii[1])
    patch = patch.segment(points[2],width_feedline[1])
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
    # patch_end2 = [points2[2][0]-1,points2[2][1]]
    
    big_patch2 = gdspy.FlexPath([patch_start2,points2[1]], width_guide[0])#,corners=#'circular bend',bend_radius = bend_radii[1])
    big_patch2 = big_patch2.segment(points2[2],width_guide[1])
    guide2 = gdspy.FlexPath([points2[2],points2[3]],width_guide[1],corners='circular bend',bend_radius=bend_radii[2])
    guide2 = guide2.segment(points2[4], width_guide[2]).segment(points2[5],width_guide[3]).segment(points2[6],width_guide[4]).segment(points2[7],width_guide[5])
   
    
    patch2 = gdspy.FlexPath([points2[0],points2[1]], width_feedline[0])#,corners='circular bend',bend_radius = bend_radii[1])
    patch2 = patch2.segment(points2[2],width_feedline[1])   
    feedline2 = gdspy.FlexPath([points2[2],points2[3]],width_feedline[1],corners='circular bend',bend_radius=bend_radii[2])
    feedline2 = feedline2.segment(points2[4],width_feedline[2]).segment(points2[5],width_feedline[3]).segment(points2[6],width_feedline[4]).segment(points2[7],width_feedline[5])

    pad2 = gdspy.boolean(big_patch2, patch2, 'not')
    line2 = gdspy.boolean(guide2,feedline2, 'not')
    
    window_left_guide = gdspy.boolean(window, line, 'or')
    window_both_guides = gdspy.boolean(window_left_guide, line2, 'or')
    window_left_feedline = gdspy.boolean(window_both_guides,feedline ,'not')
    window_both_feedlines = gdspy.boolean(window_left_feedline,feedline2,'not')
    pads = gdspy.boolean(pad,pad2,'or')
    
    # mask = gdspy.Rectangle([x1+R-maskdim[0],y1+2*R-maskdim[1]], [x2-R+maskdim[0],y2-2*R+maskdim[1]])
    # mask_overlap = gdspy.Rectangle([x1+R - maskdim[0] - d_overlap, y1+2*R-maskdim[1]-d_overlap], [x2-R+maskdim[0] + d_overlap,y2-2*R+maskdim[1] + d_overlap])
    mask = gdspy.Rectangle([x1-maskdim[0],y1-maskdim[1]], [x2+maskdim[0],y2+maskdim[1]])
    mask_overlap = gdspy.Rectangle([x1 - maskdim[0] - d_overlap, y1-maskdim[1]-d_overlap], [x2+maskdim[0] + d_overlap,y2+maskdim[1] + d_overlap])
    
    return window_both_feedlines,pads,mask,mask_overlap

def T_feedline_extended_negative(Q, unitcell_size,startM, startU, stopU, strip_height, tw, tv, N_ghost, t = 2e-6, Sr = [1e-6,2e-6], Sg = 60e-6,T=100e-6,D=200e-6,R=200e-6,f=24e-6, A = 50e-6,B=60e-6, w=[360e-6,90e-6,10e-6],maskdim = [150e-6,300e-6],laserwriter=False):
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
    
    #(x1,y1) and (x2,y2): corners of outer box, (x3,y3) and (x4,y4): corners of window
    x3, y3 = startM[0] - A/2 - Sr[0] - 2*w_start, startM[1] - f
    x4, y4 = x3 + 2 * Sr[0] + 4*w_start + Q*unitcell_size[0]-tw, y3 + strip_height
    x1, y1 = x3 - R, startM[1] - f - T
    x2, y2 = x4 + R, y1 + 2*T + strip_height    
    
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
    Y0 = y1 + 1/2*(2*T + strip_height)
    points = [[x1-E,Y0],[x1-7*E/8, Y0],[x1-3*E/4, Y0],[x1,Y0],[x3 - D, Y0],[startU[0] + t/2 - A - Sr[0] - w_start/2, Y0]]
    
    #enclose contact pad with etched line (right edge)
    patch_start = [points[0][0]- dyke[0],points[0][1]]
    # patch_end = [points[2][0]+1,points[2][1]]
    
    #enclosing contact pad/patch to define it/ this big_patch + guide will be etched away
    big_patch = gdspy.FlexPath([patch_start,points[1]], width_guide[0])#,corners=#'circular bend',bend_radius = bend_radii[1])
    big_patch = big_patch.segment(points[2],width_guide[1])
    guide = gdspy.FlexPath([points[2],points[4]],width_guide[1]).segment(points[5],width_guide[2])
    T_bar_guide = gdspy.Rectangle([points[5][0] - w_start/2 - dyke[2], startU[1] + B + dyke[2]], [points[5][0] + w_start/2 + dyke[2], startU[1] - dyke[2]])
    
    guide = gdspy.boolean(guide, T_bar_guide, 'or')
    
    #draw the feedline (not etched away)
    patch = gdspy.FlexPath([points[0],points[1]], width_feedline[0])#,corners='circular bend',bend_radius = bend_radii[1])
    patch = patch.segment(points[2],width_feedline[1])
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
------------------------------------------------------------------------------
hanged waveguide geometry (positive + negative)
------------------------------------------------------------------------------
"""
def hanged_waveguide(L,w,N,dim_wells,first_object,c=30e-6,d=70e-6,P=360e-6,Q=530e-6):
    """
    Parameters
    ----------
    L : float
        LENGTH OF TRANSMISSION LINE WITHOUT PADS. (SI units)
    w : float
        WIDTH OF TRANSMISSION LINE. (SI units)
    N : integer
        NUMBER OF WELLS ALONG TRANSMISSION LINE.
    dim_wells : N-dimensional list of tuples
        DIMENSION OF WELLS in [x,y] direction
    c : float, optional
        SPACING TRANSMISSION LINE TO GROUND PLANE. The default is 30e-6. (SI units)
    d : float, optional
        SPACING PADS TO GROUND PLANE. The default is 70e-6. (SI units)
    P : float, optional
        HEIGHT OF PADS. The default is 360e-6. (SI units)
    Q : float, optional
        WIDTH OF PADS. The default is 530e-6. (SI units)
    first_object : list of 1 tuple and 1 number
        GIVE TOP-COORDINATES (i.e. x-coordinate of center, y-coordinate = top edge of object; in um) and spacing to feedline

    Returns
    -------
    waveguide : TYPE
        DESCRIPTION.
    wells : TYPE
        DESCRIPTION.

    """
    import gdspy
    import ctypes
    import numpy as np
    
    L = L*1e6
    w = w*1e6
    A = [a[1] for a in dim_wells]
    B = [b[0] for b in dim_wells]
    c = c*1e6
    d = d*1e6
    P = P*1e6
    Q = Q*1e6
    R = max(A) + P
    
    
    if len(A)%N > 0:
        print('Wrong amount of specified well depths.')
        ctypes.windll.user32.MessageBoxW(0,'Wrong amount of specified well depths.','hanged_waveguide',1)
        return
    if len(B)%N > 0:
        print('Wrong amount of specified well widths!')
        ctypes.windll.user32.MessageBoxW(0,'Wrong amount of specified well widths!','hanged_waveguide',1)
        return
    
    #calculate alpha (spacing between wells)
    alpha = 1/(N+1)*(L - sum(B))
    if alpha < 0:
        print('Please choose smaller well widths or a smaller number of them, cannot be spaced properly!')
        ctypes.windll.user32.MessageBoxW(0,'Please choose smaller well widths or a smaller number of them, cannot be spaced properly!','hanged_waveguide',1)
        return
    
    
    #draw big ground box
    origin = first_object[0][0]- B[0]/2 - alpha - 2*Q - d, first_object[0][1] + first_object[1]*1e6 + w/2 - P/2 - d - R
    grdbox_width = 2*d + 4*Q + L
    grdbox_height = 2*R + 2*d + P
    grdbox = gdspy.Rectangle(origin, [origin[0]+grdbox_width, origin[1]+grdbox_height])
    
    
    
    #draw feedline
    p1,q1 = origin[0]+Q+d, first_object[0][1] + first_object[1]*1e6 + w/2 
    p2,q2 = p1 + Q/2, q1
    p3,q3 = p1 + Q, q1
    p4,q4 = p3 + L, q1 
    p5,q5 = p4 + Q/2, q1
    p6,q6 = p4 + Q, q1
    
    feedline = gdspy.FlexPath([[p1,q1],[p2,q2]],P)
    feedline = feedline.segment([p3,q3],w).segment([p4,q4],w).segment([p5,q5],P).segment([p6,q6],P)
    guide = gdspy.FlexPath([[p1-d,q1],[p2,q2]],P+2*d)
    guide = guide.segment([p3,q3],w+2*c).segment([p4,q4],w+2*c).segment([p5,q5],P+2*d).segment([p6+d,q6],P+2*d)
    
    #draw the wells
    x11,y11 = origin[0] + d + 2*Q + alpha, origin[1] + R + d + P/2 - w/2 - c
    x12,y12 = x11 + B[0], y11 - A[0]
    
    wells = gdspy.Cell('wells')
    
    well_1 = gdspy.Rectangle([x11,y11],[x12,y12])
    wells.add(well_1)
    
    center_wells_x = [first_object[0][0]]
    
    w2 = x12
    x1, y1 = 0.0, 0.0
    x2, y2 = 0.0, 0.0
    
    for i in range(1,len(A)):
        x1,y1 = w2 + alpha, y11
        x2,y2 = x1 + B[i], y1 - A[i]
        well_i = gdspy.Rectangle([x1,y1], [x2,y2])
        center_i_x = x1 + B[i]/2
        center_wells_x.append(center_i_x)
        wells.add(well_i)
        w2 = x2
    
    
    guide = gdspy.boolean(guide,wells,'or')
    waveguide = gdspy.boolean(grdbox, guide, 'not')
    waveguide = gdspy.boolean(waveguide, feedline,'or')
    
    # bottom_transline_edge = q1 - w/2
    
    shift_x = [center_wells_x[0]]
    
    for i in range(1,len(A)):
        dxh= np.abs(center_wells_x[i] - center_wells_x[0])
        shift_x.append(dxh)
    
    return waveguide, shift_x

def hanged_waveguide_negative(L,w,N,dim_wells,first_object,c=30e-6,d=70e-6,P=360e-6,Q=530e-6):
    """
    Parameters
    ----------
    L : float
        LENGTH OF TRANSMISSION LINE WITHOUT PADS. (SI units)
    w : float
        WIDTH OF TRANSMISSION LINE. (SI units)
    N : integer
        NUMBER OF WELLS ALONG TRANSMISSION LINE.
    dim_wells : N-dimensional list of tuples
        DIMENSION OF WELLS in [x,y] direction
    c : float, optional
        SPACING TRANSMISSION LINE TO GROUND PLANE. The default is 30e-6. (SI units)
    d : float, optional
        SPACING PADS TO GROUND PLANE. The default is 70e-6. (SI units)
    P : float, optional
        HEIGHT OF PADS. The default is 360e-6. (SI units)
    Q : float, optional
        WIDTH OF PADS. The default is 530e-6. (SI units)
    first_object : list of 1 tuple and 1 number
        GIVE TOP-COORDINATES (i.e. x-coordinate of center, y-coordinate = top edge of object; in um) and spacing to feedline

    Returns
    -------
    waveguide : TYPE
        DESCRIPTION.
    wells : TYPE
        DESCRIPTION.

    """
    import gdspy
    import ctypes
    import numpy as np
    
    L = L*1e6
    w = w*1e6
    A = [a[1] for a in dim_wells]
    B = [b[0] for b in dim_wells]
    c = c*1e6
    d = d*1e6
    P = P*1e6
    Q = Q*1e6
    R = max(A) + P
    
    
    if len(A)%N > 0:
        print('Wrong amount of specified well depths.')
        ctypes.windll.user32.MessageBoxW(0,'Wrong amount of specified well depths.','hanged_waveguide',1)
        return
    if len(B)%N > 0:
        print('Wrong amount of specified well widths!')
        ctypes.windll.user32.MessageBoxW(0,'Wrong amount of specified well widths!','hanged_waveguide',1)
        return
    
    #calculate alpha (spacing between wells)
    alpha = 1/(N+1)*(L - sum(B))
    if alpha < 0:
        print('Please choose smaller well widths or a smaller number of them, cannot be spaced properly!')
        ctypes.windll.user32.MessageBoxW(0,'Please choose smaller well widths or a smaller number of them, cannot be spaced properly!','hanged_waveguide',1)
        return
    
    
    #draw big ground box
    origin = first_object[0][0]- B[0]/2 - alpha - 2*Q - d, first_object[0][1] + first_object[1]*1e6 + w/2 - P/2 - d - R
    grdbox_width = 2*d + 4*Q + L
    grdbox_height = 2*R + 2*d + P
    grdbox = gdspy.Rectangle(origin, [origin[0]+grdbox_width, origin[1]+grdbox_height])
    
    
    
    #draw feedline
    p1,q1 = origin[0]+Q+d, first_object[0][1] + first_object[1]*1e6 + w/2 
    p2,q2 = p1 + Q/2, q1
    p3,q3 = p1 + Q, q1
    p4,q4 = p3 + L, q1 
    p5,q5 = p4 + Q/2, q1
    p6,q6 = p4 + Q, q1
    
    feedline = gdspy.FlexPath([[p1,q1],[p2,q2]],P)
    feedline = feedline.segment([p3,q3],w).segment([p4,q4],w).segment([p5,q5],P).segment([p6,q6],P)
    guide = gdspy.FlexPath([[p1-d,q1],[p2,q2]],P+2*d)
    guide = guide.segment([p3,q3],w+2*c).segment([p4,q4],w+2*c).segment([p5,q5],P+2*d).segment([p6+d,q6],P+2*d)
    
    #draw the wells
    x11,y11 = origin[0] + d + 2*Q + alpha, origin[1] + R + d + P/2 - w/2 - c
    x12,y12 = x11 + B[0], y11 - A[0]
    
    wells = gdspy.Cell('wells')
    
    well_1 = gdspy.Rectangle([x11,y11],[x12,y12])
    wells.add(well_1)
    
    center_wells_x = [first_object[0][0]]
    
    w2 = x12
    x1, y1 = 0.0, 0.0
    x2, y2 = 0.0, 0.0
    
    for i in range(1,len(A)):
        x1,y1 = w2 + alpha, y11
        x2,y2 = x1 + B[i], y1 - A[i]
        well_i = gdspy.Rectangle([x1,y1], [x2,y2])
        center_i_x = x1 + B[i]/2
        center_wells_x.append(center_i_x)
        wells.add(well_i)
        w2 = x2
    
    
    guide = gdspy.boolean(guide,wells,'or')
    waveguide = gdspy.boolean(guide,feedline, 'not')
    
    shift_x = [center_wells_x[0]]
    
    for i in range(1,len(A)):
        dxh= np.abs(center_wells_x[i] - center_wells_x[0])
        shift_x.append(dxh)
    
    return waveguide, shift_x



"""
~~~~~~~~~~~~ Testing functions ~~~~~~~~
"""


if __name__ == '__main__':
    L = 279e-6              #total length of inductor in SI units
    s = 4.5e-6              #interspacing of inductor in SI units
    w = 0.5e-6              #width of inductor wire
    Ac = 50e-6               #horizontal dimension of capacitor
    tc = 2e-6                #thickness of capacitor plates
    tv = 12e-6              #intracell spacing
    tw = 24e-6               #intercell spacing
    k = 1/3                 #fraction determining how thick the ground strip between resonators is : k*min(tv,tw)
    ep = 20.5e-6             #horizontal dimension of ground patches
    fp = 24e-6               #vertical dimension of ground patches
    rp = 29e-6               #vertical spacing of resonator "head" to ground plane
    Q = 4                   #number of unit cells
    N_ghost = 2             #number of ghosts: For no ghosts, put zero
    Sf2r = [10e-6,0]      #spacing to feedline in [x,y]-direction
    ground_yn = True        #ground in between yes (True) or no (False) 
    
    unitcell_size = [2*Ac*1e6 + tv*1e6 + tw*1e6,0]
    
    test = gdspy.Cell('negative')
    test_new = gdspy.Cell('negative new')

    
    unitcell, ground, startM,startU, stopUy, strip_height, Bc, center1st, center2nd = unit_cell_GGG(L, s, w, Ac, tc, tv,tw,ep,fp,rp,ground_in_between = ground_yn)
    
    ResArray = gdspy.CellArray(unitcell, Q, 1, unitcell_size)
    test.add(ResArray)
    
    stopU = [startU[0]-Ac*1e6 + Q*unitcell_size[0]-tw*1e6, stopUy]
    centerQth = [center1st[0]-Ac*1e6 + Q*unitcell_size[0] - tw*1e6, center2nd[1]]
    
    # blib = only_waveguide(startM, startU, strip_height)
    # blub = waveguide_extended_new(Q, unitcell_size,startM, startU, stopU, strip_height, tw, tv, N_ghost, t = tc, Sr = Sf2r,f=fp, A = Ac)
    blob = waveguide_extended_negative_new(Q, unitcell_size, startM, startU, stopU, strip_height, tw, tv, N_ghost, Sr=Sf2r,laserwriter=True,maskdim=[300e-6,250e-6])

    blub = T_feedline_extended_negative(Q,unitcell_size,startM,startU,stopU,strip_height,tw,tv,N_ghost,t=tc,Sr=Sf2r,f=fp,A=Ac,B=Bc*1e-6,laserwriter=True,maskdim=[300e-6,300e-6])

    # blub = T_feedline_simulation(Q, unitcell_size,startM, startU, stopU, strip_height, tw, tv, t = tc, Sr = Sf2r,f=fp, A = Ac, B = Bc*1e-6, R=10e-6)

    # blob, blab = ghosts_U_GGG(L,s,w,Ac,tc,tv,tw,N_ghost, strip_height, tg = tw, e=ep, f=fp, r=rp, ground_in_between=ground_yn, center_first = center1st, center_last = centerQth)

    test.add(blub)
    # test.add(blab)
    test_new.add(blob)
    # test.add(bleb)
    lib = gdspy.GdsLibrary()
    lib.add(test)
    lib.add(test_new)
    gdspy.LayoutViewer()
    # lib.write_gds("Only waveguide.gds")