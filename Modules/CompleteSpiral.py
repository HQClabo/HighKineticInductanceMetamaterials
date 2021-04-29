# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 09:05:28 2021

@author: jouanny
"""

import gdspy
import numpy as np
import matplotlib.pyplot as plt

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


# L = 0.75e-3
# t = 500e-9
# Lcouple = 100e-6
# Pad = 20e-6, 20e-6
# S = CompTwirlSpirPad(L, t, Lcouple, Pad)


# Right = gdspy.FlexPath(S[0],t*1e6, ends = 'extended')
# Left = gdspy.FlexPath(S[1], t*1e6)
# PadTop = gdspy.Rectangle(S[2][0],S[2][1])

# Spiral = gdspy.Cell('Spiral')
# Spiral.add(Left)
# Spiral.add(Right)
# Spiral.add(PadTop)


# lib = gdspy.GdsLibrary()

# lib.add(Spiral)

# lib.write_gds('SpiralPad.gds')
    
# L = 2e-3
# t = 250e-9

# Spiral = CompTwirlSpir(L, t)

# CellSpiral = gdspy.Cell('RightArm')
# RightArm = gdspy.FlexPath(Spiral[0],t*1e6, ends = 'extended')
# LeftArm = gdspy.FlexPath(Spiral[1],t*1e6, ends = 'extended')
# CellSpiral.add(RightArm)
# CellSpiral.add(LeftArm)


# lib = gdspy.GdsLibrary()
# lib.add(CellSpiral)

# lib.write_gds("NewSpiral.gds")

# Llist = np.linspace(200e-6, 2000e-6, 1000)
# Llist = np.asarray(Llist)
# Q = []
# for k in Llist:
#     Q.append(CompTwirlSpir(k, t)[2])


# plt.title('Error on the length')    
# plt.plot(Llist*1e6, Q)
# plt.xlabel(r'L ($\mu m$)')
# plt.ylabel(r'$\varepsilon_L (\mu m)$')
# plt.grid()
# plt.show()

