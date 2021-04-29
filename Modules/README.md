#Overview
The package Modules contains all functions needed to draw our resonator shapes + waveguides with gdspy.
Currently, they are all included in one file named ModuleResonator.py.
Most generally, include
```python
import Modules.ModuleResonator
```
in the preamble and to call a function,
```python
Modules.ModuleResonator.spiral(**args)
```
The module is divided according to the different resonator shapes and their corresponding waveguide geometry. They follow this order:
1. Twirled spiral: old
2. Twirled spiral: new
3. Meandered resonator
4. Left-handed LC-Resonator metamaterial (design + etch mask)
5. Hanged waveguide geometry (design + etch mask)

#4. Left-handed LC-Resonator metamaterial
##Etch mask
###waveguide_negative(Q,unitcell_size,startM,startU,stopU,strip_height,tw,tv,N_ghost,t=2e-6, Sr = [1e-6,2e-6],Sg = 60e-6, T = 100e-6, R = 200e-6, f = 24e-6, A = 50e-6)
_returns window_both_feedlines (type: PolygonSet)_
Draws the etch mask for a waveguide with a bended feedline without contact pads.
![waveguide_negative][https://github.com/HQClabo/HighKineticInductanceMetamaterials/blob/main/Modules/pictures/waveguide_negative.png]

###ghosts_negative(L,s,w,A,t,tv,tw,N, strip_height, tg = 15e-6, e=20.5e-6, f=24e-6, r=29e-6, gamma=1/4, ground_in_between=True, carac = {'layer' :  1, 'datatype' : 1}, center_first = (0,0), center_last = (0,0))
_returns ghostcell, ground (type: Cell, Cell)_
Draws the ghosts and their ground strips and ground pads. For N = 0 (no ghosts), ghostcell is just a rectangle (do not add to cell). However, add ground to cell, it contains outermost ground strips.
![ghosts_negative][https://github.com/HQClabo/HighKineticInductanceMetamaterials/blob/main/Modules/pictures/ghosts_negative.png]

###waveguide_extended_negative(Q, unitcell_size,startM, startU, stopU, strip_height, tw, tv, N_ghost, t = 2e-6, Sr = [1e-6,2e-6], Sg = 60e-6,T=100e-6,R=200e-6,f=24e-6, A = 50e-6)
_returns window_both_feedlines, pads (type: PolygonSet, PolygonSet)_
Draws the waveguide with contact pads continously narrowed down to feedline end (not 50 Ohm).
![waveguide_extended_negative][https://github.com/HQClabo/HighKineticInductanceMetamaterials/blob/main/Modules/pictures/waveguide_extended_negative.png]

###waveguide_extended_negative_new(Q, unitcell_size,startM, startU, stopU, strip_height, tw, tv, N_ghost, t = 2e-6, Sr = [1e-6,2e-6], Sg = 60e-6,T=100e-6,R=200e-6,f=24e-6, A = 50e-6, w = [360e-6,90e-6,10e-6,2e-6])
_returns window_both_feedlines,pads,mask,mask_overlap (type: PolygonSet, PolygonSet, Rectangle, Rectangle)_
New feedline geometry to have a 50 Ohm-line in the middle.
![waveguide_extended_negative_new][https://github.com/HQClabo/HighKineticInductanceMetamaterials/blob/main/Modules/pictures/waveguide_extended_negative_new.png]

###T_feedline_extended_negative(Q, unitcell_size,startM, startU, stopU, strip_height, tw, tv, N_ghost, t = 2e-6, Sr = [1e-6,2e-6], Sg = 60e-6,T=100e-6,D=200e-6,R=200e-6,f=24e-6, A = 50e-6,B=60e-6, w=[360e-6,90e-6,10e-6])
_returns window_both_feedlines,pads,mask,mask_overlap_
Draws the waveguide with contact pads and a T-feedline to the resonator. In the middle matched to 50 $\Omega$.
![T_feedline_extended_negative][https://github.com/HQClabo/HighKineticInductanceMetamaterials/blob/main/Modules/pictures/T_feedline_extended_negative.png]