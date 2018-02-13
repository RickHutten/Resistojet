import numpy as np
from scipy import interpolate as ip
from matplotlib.mlab import griddata as gd

"""
Description:
This file is used to determine propellant properties at enthalpy and pressures. The file uses two functions, the temperatures function T(h,p) and phase function a(h,p)

Author: 
Pim Stohr

Version: 
1.0

Critical parameters:
p		Pressure of the propellant
T		Temperature of the propellant
h		Enthalpy (relative) of the propellant 
a		Binary phase indicator (0 is liquid, 1 is gas)
prop 	Propellant under consideration (string)
"""


def properties(prop):
    if prop == "N2":
        M=28        #[g/mol] Molar mass
        gamma=1.40   #[-] Specific heat ratio
    elif prop == "Ar":
        M=40         #[g/mol] Molar mass
        gamma=1.66  #[-] Specific heat ratio
    elif prop == "H2O":
        M=18        #[g/mol] Molar mass
        gamma=1.33   #[-] Specific heat ratio
    elif prop == "NH3":
        M=17        #[g/mol] Molar mass
        gamma=1.32   #[-] Specific heat ratio
    else:
        M=18
        gamma=1.33
    R=8.314462175/M*1000

    return(R,gamma,)