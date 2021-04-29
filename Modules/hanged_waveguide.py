# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 14:51:59 2021

@author: vweibel
"""
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


# if __name__ == '__main__':
#     import gdspy
#     import CompleteSpiral
#     import numpy as np
    
#     n = 3
#     L_guide = 4.79e-3
#     w_guide = 1e-4
#     L_spir = 790e-6
#     t_spir = 500e-9
#     center_spir = [0.0,0.0]
    
#     Twirl_right, Twirl_left = CompleteSpiral.CompTwirlSpir(L_spir, t_spir,centre = center_spir)
#     Height_spir = np.abs(Twirl_left[-1][0]-Twirl_right[-1][0]) + t_spir*1e6
#     Width_spir = np.abs(Twirl_right[-2][1] - Twirl_right[-1][1]) + t_spir*1e6
    
#     size = [[3*Width_spir,2*Height_spir]]*n
#     well_widths = [3*Width_spir*1e-6]*n
#     well_depths = [2*Height_spir*1e-6]*n
#     spacings = [2e-6]*n
    
#     first_spir = [[center_spir[0],center_spir[1]+Height_spir/2],spacings[0]] 
   
#     whuy = hanged_waveguide_negative(L_guide,w_guide,n,size,first_spir) 
   
#     CellSpiral = gdspy.Cell('RightArm')
    
#     RightArm = gdspy.FlexPath(Twirl_right,t_spir*1e6, ends = 'extended').rotate(np.pi/2,center=center_spir)
#     LeftArm = gdspy.FlexPath(Twirl_left,t_spir*1e6, ends = 'extended').rotate(np.pi/2,center=center_spir)
#     CellSpiral.add(RightArm)
#     CellSpiral.add(LeftArm)
    
    
#     # for i in range(1,n):
#     #     center_spir = whuy[1][i], whuy[2] - spacings[i]*1e6- Height_spir/2
#     #     Twirl_right, Twirl_left = CompleteSpiral.CompTwirlSpir(L_spir, t_spir,centre=center_spir)
#     #     RightArm = gdspy.FlexPath(Twirl_right,t_spir*1e6, ends = 'extended').rotate(np.pi/2,center=center_spir)
#     #     LeftArm = gdspy.FlexPath(Twirl_left,t_spir*1e6, ends = 'extended').rotate(np.pi/2,center=center_spir)
#     #     CellSpiral.add(RightArm)
#     #     CellSpiral.add(LeftArm)

#     CellSpiral.flatten()
    
#     test = gdspy.Cell('Test')
    
#     test.add(whuy[0])
#     test.add(CellSpiral)
#     test.flatten()
#     lib = gdspy.GdsLibrary()
#     lib.add(test)
#     lib.write_gds('testing.gds')
#     # gdspy.LayoutViewer()



    