# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 08:59:35 2021

@author: vweibel
"""

import scipy.linalg as slin
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio

def lorentzian(x,x0,gamma):
    return 1/(np.pi*gamma)*(gamma**2/((x-x0)**2 + gamma**2))

def S21_mag_single_peak(f,f0,k,S_max):
    """
    f : frequency array
    f0 : resonance frequency
    k : kappa/2pi
    """
    S21 = k/(k-1j*(f-f0))
    return S_max*k/np.sqrt(k**2 + (f-f0)**2)

def SSH_Hamiltonian(N,w0,v,w,s,delta):
    H = np.zeros((N,N))
    for i in range(N):
        H[i,i] = w0
        if i == 1:
            H[i,i-1] = v
            H[i-1,i] = v
        if i >= 2 and i%2 == 0:
            H[i,i-1] = w
            H[i-1,i] = w
        if i >= 3 and (i+1)%2 == 0:
            H[i,i-1] = v
            H[i-1,i] = v
        if i >= 2:
            H[i,i-2] = s
            H[i-2,i] = s
    H[0,0] = w0 + delta
    H[N-1,N-1] = w0 + delta
    return H


w0 = 4.86e9
ratio = np.linspace(0,5,num=500)
v = 100e6
w = v*2 #cw_array = v*ratio
s = 0e6
N = 16
second_neighbour = True
span = 50e6
k = 1e6
Smax = 0.008
delta = np.linspace(-500e6,500e6,num=51)


# path = r"C:\Users\vweibel\switchdrive\Private\Masterthesis\Measurements\2021-06\Data\16LC-ghost-Resonator\16LCghost_topolo\VNA\06-22-2021\10.34.45_powersweep_in_transmission"

# file = sio.loadmat(path+r'\16LCghost_topolo_powersweep_in_transmission.mat')

# FREQ = file['Frequency'][0]
# POW = file['Power'][0]
# S21_MAG = np.abs(file['ComplexResponse']) #convert to linear scale

lambdas, eigvecs = slin.eig(SSH_Hamiltonian(N, w0, v, w, 0, 0))

# plt.plot(eigvecs[2])

H_array = []
w_array = []

for d in delta:
    H_array.append(SSH_Hamiltonian(N, w0, v, w, s, d))
    w_array.append(np.real(slin.eigvals(H_array[-1])))

plt.rcdefaults()
plt.rcParams["font.sans-serif"] = "Yu Gothic"
plt.rcParams["font.size"] = 16
plt.figure(dpi=200)
plt.plot(w_array, delta,'ko',markersize=0.9)
plt.xlabel(r'Eigenvalues')
plt.ylabel(r'$\delta$')
plt.show()

# for cw in cw_array:
#     H_array.append(SSH_Hamiltonian(N,w0,v,cw,0,0))
#     w_array.append(np.real(slin.eigvals(H_array[-1])))

# plt.plot(ratio,w_array,'ko',markersize=0.9)

# plt.plot(w_array[25],'o')

# w = [[(l-span),(l+span)] for l in lambdas]

# t = [Smax for l in lambdas]


# for j in range(len(lambdas)):
#     x = np.linspace(w[j][0],w[j][1],num=2000)
#     plt.plot(x,S21_mag_single_peak(x,lambdas[j],k,Smax),'b')
#     plt.plot(lambdas[j], t[j],'ro')
#     plt.ylim(0,Smax)
# # plt.plot(FREQ,S21_MAG[6],label=str(POW[6])+' dBm', color='r')
# plt.legend()