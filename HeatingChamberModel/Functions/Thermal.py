import numpy as np

"""
Description:
Thermal analysis of the design. The methods are used to determine conductive constants [W/K], used to solve the heating functions.

Author:
Pim Stohr

Version:
1.0

Critical parameters:
L		Length of the chamber
r1		Internal radius 
r2		External radius of shell
A		Area
t		thickness of the connection
k		Thermal conductivity coefficient
"""

#Conductive constant for a 1-D case
def cond_1D(k,A,t):
  #Fourier law applied to a 1D case
  constant	= k*A/t                       #[W/K]
  return constant

#Conductive constant for a cylindrical shell
def cond_shell(k,L,r1,r2):
  #Fourier law applied to a cylindrical shell
  constant	= 2*k*np.pi*L/np.log(r2/r1)   #[W/K]
  return constant

#Combining conductive constants
def cond_comb(C1,C2):
  #Combining the conductive constants
  constant	= C1*C2/(C1+C2)                #[W/K]
  return constant