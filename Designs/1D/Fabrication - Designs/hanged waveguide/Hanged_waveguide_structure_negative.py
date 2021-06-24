# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 12:46:01 2021

@author: vweibel
"""
import gdspy
from Modules.ModuleResonator import CompTwirlSpirPad, CompTwirlSpir, OPKI_down, hanged_waveguide_negative
import numpy as np

def hanged_waveguide_negative(L,w,N,dim_wells,first_object,c=30e-6,d=70e-6,P=360e-6,Q=530e-6,overlap=1e-6):
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
    overlap = overlap*1e6
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
    mask = gdspy.Cell('mask')
    mask_overlap = gdspy.Cell('mask overlap')
    
    k = 1/20
    
    well_1 = gdspy.Rectangle([x11,y11],[x12,y12])
    mask_1 = gdspy.Rectangle([x11-alpha*k,y11+c],[x12+alpha*k,y12-alpha*k])
    mask_overlap_1 = gdspy.Rectangle([x11-alpha*k-overlap,y11+c],[x12+alpha*k+overlap,y12-alpha*k-overlap])
    wells.add(well_1)
    mask.add(mask_1)
    mask_overlap.add(mask_overlap_1)
    
    center_wells_x = [first_object[0][0]]
    
    w2 = x12
    x1, y1 = 0.0, 0.0
    x2, y2 = 0.0, 0.0
    
    for i in range(1,len(A)):
        x1,y1 = w2 + alpha, y11
        x2,y2 = x1 + B[i], y1 - A[i]
        well_i = gdspy.Rectangle([x1,y1], [x2,y2])
        mask_i = gdspy.Rectangle([x1-alpha*k,y1+c],[x2+alpha*k,y2-alpha*k])
        mask_overlap_i = gdspy.Rectangle([x1-alpha*k-overlap,y1+c], [x2+alpha*k+overlap,y2-alpha*k-overlap])
        center_i_x = x1 + B[i]/2
        center_wells_x.append(center_i_x)
        wells.add(well_i)
        mask.add(mask_i)
        mask_overlap.add(mask_overlap_i)
        w2 = x2
    
    
    guide = gdspy.boolean(guide,wells,'or')
    waveguide = gdspy.boolean(guide,feedline, 'not')
    
    shift_x = [center_wells_x[0]]
    
    for i in range(1,len(A)):
        dxh= np.abs(center_wells_x[i] - center_wells_x[0])
        shift_x.append(dxh)
    
    return waveguide, shift_x, mask, mask_overlap


'''
specifiy relevant parameters
'''
n = 4                                 #how many objects to put along the transmission line
alpha = 2                             #alpha gives the ratio between (well width):(object width)
beta = 1.2                            #beta is the ration between (well depth):(object height)
L_feed = 4.79e-3                      #length of the transmission line (without pads)
w_feed = 2e-4                         #width of transmission line
c_feed = 4e-6                        #separation of transmission line to ground
L_spir = [250e-6,279e-6,380e-6,500e-6]       #length of spiral
s_M = 4.5e-6                          #interspacing of inductor
t_spir = 0.5e-6                       #width of spiral
L_couple =100e-6                      #length of spiral parallel to transmission line, ending in square capacitor
dim_pad1 = np.array([25e-6,25e-6])    #dimension of capacitor pad of spiral
Ac = 50e-6                            # condensator width of LC-Resonator
tc = 2e-6                             # condensator thickness of LC-Resonator
center1 = [0.0,0.0]                   #center of first object, used as a reference point for everything else
N_array = 1                           #number of unit cells in spiral array
M_uc = 2                              #how many objects in one unit cell
tv = 5e-6                             #intracell spacing of spiral array
tw = 10e-6                            #intercell spacing of spiral array
spacings_to_feed = 60e-6               #spacing of each object to the feedline

spiral_yn = False
array_in_2nd = False

'''
Next, each object is drawn with its center at the reference point (i.e. origin of the first object).
Later, they will be translated to their right position. Why do like this? To get the dimensions of each object to
being able to draw the hanged waveguide around them.
Naming logic: every variable used for the i-th object has the number i in it.
The spirals (/ and pads) need to be rotated sometimes around reference point (=center1) after generation.
Whenever a center needs to be specified, use the reference point (=center1).

OBJECT IS A SINGLE SPIRAL WITH/WITHOUT PAD AT THE END
1) call function to calculate the coordinates of the spiral (/ and the pad) (CompTwirlSpir/CompTwirlSpirPad)
2) draw the spiral (/ and the pad)
3) create cell and add the spiral (/ and the pad)
4) get dimensions of the object by calling .get_bounding_box() on the cell -> gives the corners of the rectangle bounding the cell. From those, obtaining size of object is straightforward.
5) add dimensions of well to the well_size list (append [size_in_x*alpha, size_in_y*beta])
ONLY FOR FIRST OBJECT:
6) calculate the reference point to draw the hanged waveguide. It has the x-coordinate of its origin (, i.e. center1[0])
    and the y-coordinate of the topmost point of the object (i.e. center1[1] + size1[1]/2)

OBJECT IS AN ARRAY
1) draw the object of which the array consists by calling its function (here: spiral -> CompTwirlSpir) and actually drawing it
2) add this object to a cell and get the size of the cell to know the dimensions of one object
3) create unit cell with gdspy.CellArray(cell_to_repeat, 1, M_uc, [spacing in x, -(spacing in y)]).
    Since we want a vertical array, we choose 1 column and M_uc rows and spacing in x = 0.
    The spacing in y is given by the size of one object plus the intracell spacing. Add a negative sign to the spacing -(spacing in y), such that array is drawn from top to bottom.
    It is more convienient since then, the origin of the array lies at the center of the topmost spiral.
    Also, obtain size of unitcell by again calling .get_bounding_box() on the unitcell cell (badum tzz)
4) flatten the unitcell-CellArray with .flatten() s.t. it will be displayed nicely in KLayout/Sonnet
5) create the array from the unitcell using the same logic as for the unit cell:
    gdspy.CellArray(unit_cell, 1, N_array, [0.0, -(spacing in y)]) with
    spacing in y = size of unit cell + intercell spacing
    And again, determine size of array and flatten it.
6) add dimensions of well to the well_size list (append [size_in_x*alpha, size_in_y*beta])
'''

#first object: n = 1
if spiral_yn == True:
    Twirl_right1, Twirl_left1,pad1_corners = CompTwirlSpirPad(L_spir[0], t_spir, L_couple, dim_pad1)
    Spiral1 = gdspy.Cell('Spiral 1')
    RightArm1 = gdspy.FlexPath(Twirl_right1,t_spir*1e6, ends = 'extended').rotate(np.pi/2,center=center1)
    LeftArm1 = gdspy.FlexPath(Twirl_left1,t_spir*1e6).rotate(np.pi/2,center=center1)
    pad1 = gdspy.Rectangle(pad1_corners[0],pad1_corners[1]).rotate(np.pi/2,center=center1)
    Spiral1.add([RightArm1,LeftArm1,pad1])
    corners1 = Spiral1.get_bounding_box()
    size1 = [np.abs(corners1[1][0]-corners1[0][0]),np.abs(corners1[1][1]-corners1[0][1])]
    well_size = [[size1[0]*alpha,size1[1]*beta]]
    top1 = [corners1[0][0]+size1[0]/2, center1[1]+size1[1]/2]
else: #draw a LC-Resonator
    U_shape1, M_shape1, B1 = OPKI_down(L_spir[0],s_M, t_spir, Ac, tc,centre=center1,compact=False)
    U1 = gdspy.FlexPath(U_shape1, tc*1e6)
    M1 = gdspy.FlexPath(M_shape1, t_spir*1e6)
    LCRes1 = gdspy.Cell('LC-Resonator 1')
    LCRes1.add([U1,M1])
    corners1 = LCRes1.get_bounding_box()
    size1 = [np.abs(corners1[1][0]-corners1[0][0]),np.abs(corners1[1][1]-corners1[0][1])]
    well_size = [[size1[0]*alpha,size1[1]+spacings_to_feed*1e6-c_feed*1e6]]
    top1 = center1



#second object: n = 2
if spiral_yn == True and array_in_2nd == False: #2nd object is a Twirled Spiral with pad
    L_spir2 = 1e-3
    Twirl_right2, Twirl_left2,pad2_corners = CompTwirlSpirPad(L_spir2, t_spir, L_couple*1.26, dim_pad1)
    Spiral2 = gdspy.Cell('Spiral 2')
    RightArm2 = gdspy.FlexPath(Twirl_right2,t_spir*1e6, ends = 'extended').rotate(-np.pi/2,center=center1)
    LeftArm2 = gdspy.FlexPath(Twirl_left2,t_spir*1e6).rotate(-np.pi/2,center=center1)
    pad2 = gdspy.Rectangle(pad2_corners[0],pad2_corners[1]).rotate(-np.pi/2,center=center1)
    Spiral2.add([RightArm2,LeftArm2,pad2])
    corners2 = Spiral2.get_bounding_box()
    size2 = [np.abs(corners2[1][0]-corners2[0][0]),np.abs(corners2[1][1]-corners2[0][1])]
    dx2 = top1[0] - (corners2[0][0] + size2[0]/2)
    dy2 = top1[1] - (corners2[0][1] + size2[1]/2)
    RightArm2.translate(dx2, dy2)
    LeftArm2.translate(dx2, dy2)
    pad2.translate(dx2, dy2)
    well_size.append([size2[0]*alpha, size2[1]*beta])
elif spiral_yn == True and array_in_2nd == True: #2nd object is a twirled spiral array
    Twirl_right2, Twirl_left2 = CompTwirlSpir(L_spir[1], t_spir)
    Spiral2 = gdspy.Cell('Spiral 2')
    RightArm2 = gdspy.FlexPath(Twirl_right2,t_spir*1e6, ends = 'extended')#.rotate(np.pi/2,center=center1)
    LeftArm2 = gdspy.FlexPath(Twirl_left2,t_spir*1e6, ends = 'extended')#.rotate(np.pi/2,center=center1)
    Spiral2.add([RightArm2,LeftArm2])
    cornerssmall2 = Spiral2.get_bounding_box()
    sizeone2 =[np.abs(cornerssmall2[1][0]-cornerssmall2[0][0]), np.abs(cornerssmall2[1][1]-cornerssmall2[0][1])]
    dx2 = top1[0] - (cornerssmall2[0][0] + sizeone2[0]/2)
    dy2 = top1[1] - (cornerssmall2[0][1] + sizeone2[1]/2)
    RightArm2.translate(dx2, dy2)
    LeftArm2.translate(dx2, dy2)
    unitcell2 = gdspy.Cell('unitcell 2')
    Unitcell2_arr = gdspy.CellArray(Spiral2, 1, M_uc, [0,-(sizeone2[1]+tv*1e6)])  #minus sign, so that CellArray fills up from top to bottom -> origin lies at center of first spiral of unit cell
    uc_corners2 = Unitcell2_arr.get_bounding_box()
    unitcell_size2 = [np.abs(uc_corners2[1][0]-uc_corners2[0][0]),np.abs(uc_corners2[1][1]-uc_corners2[0][1])]
    unitcell2.add(Unitcell2_arr)
    unitcell2.flatten()
    Spiral2arr = gdspy.Cell('Spiral 2 Array')
    SpiralArray = gdspy.CellArray(unitcell2, 1, N_array, [0,-(unitcell_size2[1]+tw*1e6)]) #minus sign, so that CellArray fills up from top to bottom -> origin lies at center of spiral on the top of the array (the one below the transmission line)
    Spiral2arr.add(SpiralArray)
    Spiral2arr.flatten()
    corners2 = SpiralArray.get_bounding_box()
    size2 =[np.abs(corners2[1][0]-corners2[0][0]), np.abs(corners2[1][1]-corners2[0][1])]
    well_size.append([size2[0]*alpha,size2[1]*beta])
elif spiral_yn == False and array_in_2nd == False: #2nd object is a LC-Resonator
    U_shape2, M_shape2, B2 = OPKI_down(L_spir[1],s_M, t_spir, Ac, tc,centre=center1,compact=False)
    U2 = gdspy.FlexPath(U_shape2, tc*1e6)
    M2 = gdspy.FlexPath(M_shape2, t_spir*1e6)
    LCRes2 = gdspy.Cell('LC-Resonator 2')
    LCRes2.add([U2,M2])
    corners2 = LCRes2.get_bounding_box()
    size2 = [np.abs(corners2[1][0]-corners2[0][0]),np.abs(corners2[1][1]-corners2[0][1])]
    well_size.append([size2[0]*alpha,size2[1]+spacings_to_feed*1e6-c_feed*1e6])

#third object: n = 3
if spiral_yn == True: #3rd object is a spiral
    Twirl_right3, Twirl_left3, pad3_corners = CompTwirlSpirPad(L_spir[2], t_spir, L_couple, dim_pad1)
    Spiral3 = gdspy.Cell('Spiral 3')
    RightArm3 = gdspy.FlexPath(Twirl_right3,t_spir*1e6, ends = 'extended').rotate(np.pi,center=center1)
    LeftArm3 = gdspy.FlexPath(Twirl_left3,t_spir*1e6).rotate(np.pi,center=center1)
    pad3 = gdspy.Rectangle(pad3_corners[0],pad3_corners[1]).rotate(np.pi,center=center1)
    Spiral3.add([RightArm3,LeftArm3,pad3])
    corners3 = Spiral3.get_bounding_box()
    size3 = [np.abs(corners3[1][0]-corners3[0][0]),np.abs(corners3[1][1]-corners3[0][1])]
    dx3 = top1[0] - (corners3[0][0] + size3[0]/2)
    dy3 = top1[1] - (corners3[0][1] + size3[1]/2)
    RightArm3.translate(dx3, dy3)
    LeftArm3.translate(dx3, dy3)
    pad3.translate(dx3, dy3)
    well_size.append([size3[0]*alpha,size3[1]*beta])
else: #3rd object is a LC-Resonator
    U_shape3, M_shape3, B3 = OPKI_down(L_spir[2],s_M, t_spir, Ac, tc,centre=center1,compact=False)
    U3 = gdspy.FlexPath(U_shape3, tc*1e6)
    M3 = gdspy.FlexPath(M_shape3, t_spir*1e6)
    LCRes3 = gdspy.Cell('LC-Resonator 3')
    LCRes3.add([U3,M3])
    corners3 = LCRes3.get_bounding_box()
    size3 = [np.abs(corners3[1][0]-corners3[0][0]),np.abs(corners3[1][1]-corners3[0][1])]
    well_size.append([size3[0]*alpha,size3[1]+spacings_to_feed*1e6-c_feed*1e6])
    
#fourth object: n = 4
if spiral_yn == True:
    Twirl_right4, Twirl_left4, pad4_corners = CompTwirlSpirPad(L_spir[3], t_spir, L_couple, dim_pad1)
    Spiral4 = gdspy.Cell('Spiral 4')
    RightArm4 = gdspy.FlexPath(Twirl_right4,t_spir*1e6, ends = 'extended').rotate(np.pi,center=center1)
    LeftArm4 = gdspy.FlexPath(Twirl_left4,t_spir*1e6).rotate(np.pi,center=center1)
    pad4 = gdspy.Rectangle(pad4_corners[0],pad4_corners[1]).rotate(np.pi,center=center1)
    Spiral4.add([RightArm4,LeftArm4,pad4])
    corners4 = Spiral4.get_bounding_box()
    size4 = [np.abs(corners4[1][0]-corners4[0][0]),np.abs(corners4[1][1]-corners4[0][1])]
    dx4 = top1[0] - (corners4[0][0] + size4[0]/2)
    dy4 = top1[1] - (corners4[0][1] + size4[1]/2)
    RightArm4.translate(dx4, dy4)
    LeftArm4.translate(dx4, dy4)
    pad4.translate(dx4, dy4)
    well_size.append([size4[0]*alpha,size4[1]*beta])
else:
    U_shape4, M_shape4, B4 = OPKI_down(L_spir[3],s_M, t_spir, Ac, tc,centre=center1,compact=False)
    U4 = gdspy.FlexPath(U_shape4, tc*1e6)
    M4 = gdspy.FlexPath(M_shape4, t_spir*1e6)
    LCRes4 = gdspy.Cell('LC-Resonator 4')
    LCRes4.add([U4,M4])
    corners4 = LCRes4.get_bounding_box()
    size4 = [np.abs(corners4[1][0]-corners4[0][0]),np.abs(corners4[1][1]-corners4[0][1])]
    well_size.append([size4[0]*alpha,size4[1]+spacings_to_feed*1e6-c_feed*1e6])

'''
After creation of every object to be hanged on the transmission line (no resonators were harmed in production of this programm),
the waveguide needs to be actually be drawn.
The function accepts the following arguments:
hw.hanged_waveguide(L_feed, w_feed,n, well_size, carac_first_object),
where well_size is a list with n tuples as entries and carac_first_object is a list with 2 entries:
    [top1, spacings_to_feed[0]], giving the reference point to draw the waveguide (top1) and the spacing of the first object to the transmission line.
The function outputs a PolygonSet and a list of shifts (need later to translate the objects to the wells).
Add the PolygonSet to the cell with the design ('Hanged Waveguide').
'''


waveguide, shift,mask, mask_overlap = hanged_waveguide_negative(L_feed,w_feed,n,well_size,[top1, spacings_to_feed],c=c_feed)

structure = gdspy.Cell('Hanged Waveguide')
# structure.add(waveguide)

'''
SHIFT FOR THE FIRST SPIRAL
no shift needed, it already sits at the perfect spot as everything was drawn around it.
Its cell can directly be added to the structure.

SHIFT FOR i-TH SPIRAL
The spiral has its origin at its center -> shift the center in x by shift[i],
shift in y by half of the size of the spiral (- sizei[1]/2).
Note that for spirals, each arm needs to be translated seperately, as only gdspy.FlexPath, but not gdspy.Cell
have the attribute .translate(dx,dy).

SHIFT FOR j-TH SPIRAL ARRAY
The array has its center at the center of the first, i.e. the topmost spiral. So the shift is the same as
for just one spiral. Shift in x: shift[j], shift in y: - sizeonej[1]/2,
where sizeonej refers to the size of one spiral of the array
Take care to translate the CellArray, not the cell where the CellArray was put. Only gdspy.CellArray has the attribute .translate(dx,dy)

SHIFT FOR i-th SINGLE LC-RESONATOR
Since the LC-Resonators are drawn such that their center lies at their middle on top, they are already all perfectly
aligned to each other horizontally, so no need to shift in y.
Shift in x: shift[i]

'''
if spiral_yn == True:
    negative = gdspy.boolean(waveguide,RightArm1,'not')
    negative = gdspy.boolean(negative,LeftArm1,'not')
    negative = gdspy.boolean(negative,pad1,'not')
    obj3r = RightArm3.translate(shift[2],-size3[1]/2)
    obj3l = LeftArm3.translate(shift[2],-size3[1]/2)
    obj3p = pad3.translate(shift[2],-size3[1]/2)
    negative = gdspy.boolean(negative,obj3r,'not')
    negative = gdspy.boolean(negative,obj3l,'not')
    negative = gdspy.boolean(negative,obj3p,'not')
    if array_in_2nd == True:
        obj2 = SpiralArray.translate(shift[1],-sizeone2[1]/2)
        negative = gdspy.boolean(negative,obj2, 'not')
    else:
        obj2r = RightArm2.translate(shift[1],-size2[1]/2)
        obj2l = LeftArm2.translate(shift[1],-size2[1]/2)
        obj2p = pad2.translate(shift[1],-size2[1]/2)
        negative = gdspy.boolean(negative,obj2r,'not')
        negative = gdspy.boolean(negative,obj2l,'not')
        negative = gdspy.boolean(negative,obj2p,'not')
else:
    negative = gdspy.boolean(waveguide,U1,'not')
    negative = gdspy.boolean(negative,M1, 'not')
    obj2u = U2.translate(shift[1], 0)
    obj2m = M2.translate(shift[1], 0)
    obj3u = U3.translate(shift[2], 0)
    obj3m = M3.translate(shift[2], 0)
    obj4u = U4.translate(shift[3], 0)
    obj4m = M4.translate(shift[3], 0)
    negative = gdspy.boolean(negative,obj2u,'not')
    negative = gdspy.boolean(negative,obj2m,'not')
    negative = gdspy.boolean(negative,obj3u,'not')
    negative = gdspy.boolean(negative,obj3m,'not')
    negative = gdspy.boolean(negative,obj4u,'not')
    negative = gdspy.boolean(negative,obj4m,'not')
    
'''
Add to the cell with the design and flatten it just to be sure.
'''

# structure.add(negative)
# structure.flatten()

layer2 = gdspy.boolean(negative,mask,'not',layer=2)
layer4 = gdspy.boolean(negative,mask_overlap,'and',layer=4)

structure.add(layer2)
structure.add(layer4)

lib = gdspy.GdsLibrary()
lib.add(structure)
lib.write_gds('hanged_L250_L279_L380_L500_s_60um.gds')