# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 12:46:01 2021

@author: vweibel
"""
import gdspy
import Modules.ModuleResonator as hw
# import Modules.CompTwirlSpirPad as cspad
import numpy as np

'''
specifiy relevant parameters
'''

n = 3                                 #how many objects to put along the transmission line
alpha = 2                             #alpha gives the ratio between (well width):(object width)
beta = 1.2                            #beta is the ration between (well depth):(object height)
L_feed = 4.79e-3                      #length of the transmission line (without pads)
w_feed = 2e-4                         #width of transmission line
L_spir = [750e-6,650e-6,600e-6]                      #length of spiral
t_spir = 500e-9                         #width of spiral
L_couple =100e-6                     #length of spiral parallel to transmission line, ending in square capacitor
dim_pad1 = np.array([25e-6,25e-6])    #dimension of capacitor pad of spiral
center1 = [0.0,0.0]                   #center of first object, used as a reference point for everything else
N_array = 1                           #number of unit cells in spiral array
M_uc = 2                              #how many objects in one unit cell
tv = 5e-6                            #intracell spacing of spiral array
tw = 10e-6                            #intercell spacing of spiral array
spacings_to_feed = 5e-6               #spacing of each object to the feedline

'''
Next, each object is drawn with its center at the reference point (i.e. origin of the first object).
Later, they will be translated to their right position. Why do like this? To get the dimensions of each object to
being able to draw the hanged waveguide around them.
Naming logic: every variable used for the i-th object has the number i in it.
The spirals (/ and pads) need to be rotated by pi/2 around reference point (=center1) after generation.
Whenever a center needs to be specified, use the reference point (=center1).

OBJECT IS A SINGLE SPIRAL WITH/WITHOUT PAD AT THE END
1) call function to calculate the coordinates of the spiral (/ and the pad) (cs.CompTwirlSpir/cs.CompTwirlSpirPad)
2) draw the spiral (/ and the pad)
3) create cell and add the spiral (/ and the pad)
4) get dimensions of the object by calling .get_bounding_box() on the cell -> gives the corners of the rectangle bounding the cell. From those, obtaining size of object is straightforward.
5) add dimensions of well to the well_size list (append [size_in_x*alpha, size_in_y*beta])
ONLY FOR FIRST OBJECT:
6) calculate the reference point to draw the hanged waveguide. It has the x-coordinate of its origin (, i.e. center1[0])
    and the y-coordinate of the topmost point of the object (i.e. center1[1] + size1[1]/2)

OBJECT IS AN ARRAY
1) draw the object of which the array consists by calling its function (here: spiral -> cs.CompTwirlSpir) and actually drawing it
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
Twirl_right1, Twirl_left1,pad1_corners = hw.CompTwirlSpirPad(L_spir[0], t_spir, L_couple, dim_pad1)
Spiral1 = gdspy.Cell('Spiral 1')
RightArm1 = gdspy.FlexPath(Twirl_right1,t_spir*1e6, ends = 'extended').rotate(np.pi/2,center=center1)
LeftArm1 = gdspy.FlexPath(Twirl_left1,t_spir*1e6).rotate(np.pi/2,center=center1)
pad1 = gdspy.Rectangle(pad1_corners[0],pad1_corners[1]).rotate(np.pi/2,center=center1)
Spiral1.add([RightArm1,LeftArm1,pad1])
corners1 = Spiral1.get_bounding_box()
size1 = [np.abs(corners1[1][0]-corners1[0][0]),np.abs(corners1[1][1]-corners1[0][1])]
well_size = [[size1[0]*alpha,size1[1]*beta]]
top1 = [corners1[0][0]+size1[0]/2, center1[1]+size1[1]/2]

#second object: n = 2
# Twirl_right2, Twirl_left2,pad2_corners = hw.CompTwirlSpirPad(L_spir*1.5, t_spir, L_couple*1.2, dim_pad1)
# Spiral2 = gdspy.Cell('Spiral 2')
# RightArm2 = gdspy.FlexPath(Twirl_right2,t_spir*1e6, ends = 'extended').rotate(np.pi,center=center1)
# LeftArm2 = gdspy.FlexPath(Twirl_left2,t_spir*1e6).rotate(np.pi,center=center1)
# pad2 = gdspy.Rectangle(pad2_corners[0],pad2_corners[1]).rotate(np.pi,center=center1)
# Spiral2.add([RightArm2,LeftArm2,pad2])
# corners2 = Spiral2.get_bounding_box()
# size2 = [np.abs(corners2[1][0]-corners2[0][0]),np.abs(corners2[1][1]-corners2[0][1])]
# dx2 = top1[0] - (corners2[0][0] + size2[0]/2)
# dy2 = top1[1] - (corners2[0][1] + size2[1]/2)
# RightArm2.translate(dx2, dy2)
# LeftArm2.translate(dx2, dy2)
# pad2.translate(dx2, dy2)
# well_size.append([size2[0]*alpha, size2[1]*beta])
Twirl_right2, Twirl_left2 = hw.CompTwirlSpir(L_spir[1], t_spir)
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


#third object: n = 3
Twirl_right3, Twirl_left3, pad3_corners = hw.CompTwirlSpirPad(L_spir[2], t_spir, L_couple, dim_pad1)
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


waveguide, shift = hw.hanged_waveguide_negative(L_feed,w_feed,n,well_size,[top1, spacings_to_feed])

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
'''
obj2 = SpiralArray.translate(shift[1],-sizeone2[1]/2)#-spacings_to_feed[1]*1e6-size1[1]/2)#-size2[1])
# obj2r = RightArm2.translate(shift[1],-size2[1]/2)#-spacings_to_feed[2]*1e6-size3[1]/2)
# obj2l = LeftArm2.translate(shift[1],-size2[1]/2)#-spacings_to_feed[2]*1e6-size3[1]/2)
# obj2p = pad2.translate(shift[1],-size2[1]/2)# -spacings_to_feed[2]*1e6-size3[1]/2)
obj3r = RightArm3.translate(shift[2],-size3[1]/2)#-spacings_to_feed[2]*1e6-size3[1]/2)
obj3l = LeftArm3.translate(shift[2],-size3[1]/2)#-spacings_to_feed[2]*1e6-size3[1]/2)
obj3p = pad3.translate(shift[2],-size3[1]/2)# -spacings_to_feed[2]*1e6-size3[1]/2)

'''
Add the translated objects to the cell with the design and flatten it just to be sure.
'''

negative = gdspy.boolean(waveguide,obj3r,'not')
negative = gdspy.boolean(negative,obj3l,'not')
negative = gdspy.boolean(negative,obj3p,'not')
# negative = gdspy.boolean(negative,obj2r,'not')
# negative = gdspy.boolean(negative,obj2l,'not')
# negative = gdspy.boolean(negative,obj2p,'not')
negative = gdspy.boolean(negative,RightArm1,'not')
negative = gdspy.boolean(negative,LeftArm1,'not')
negative = gdspy.boolean(negative,pad1,'not')
negative = gdspy.boolean(negative, obj2, 'not')

structure.add(negative)
# structure.add(Spiral1)
# structure.add([obj2r,obj2l,obj2p,obj3r,obj3l,obj3p])
structure.flatten()

lib = gdspy.GdsLibrary()
lib.add(structure)
# lib.write_gds('3PadRes-L790um-3_2L-2L.gds')