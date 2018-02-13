import numpy as np
from Functions import Power_budget as PB
from Functions import Prop_heating as heat
from scipy import interpolate as ip
from scipy import optimize as op
import matplotlib.pyplot as plt

#Validation data
PPI     =   20          #[inch-1] Pore density of medium
eps     =   0.9         #[-] porosity of medium
p_in    =   3.5e5       #[Pa] Input pressure
T_in    =   10+273.15   #[K] Input temperature
T_s     =   210+273.15  #[K] Structure temperature
L_c     =   0.150       #[m] Length of chamber
D       =   26e-3       #[m] Diameter of chamber
C_D     =   0.7         #[-] Discharge coefficient
A_t     =               #[m2] Throat area
prop    =   "N2"        #[-] Propellant


#Loading enthalpy data
data_name   = '_data_'+prop+'.npy'                      #Name of data
data_folder = 'Functions/Data/'                         #Name of folder
p_lu        = np.linspace(1e5,10e5,10)                  #Pressure bounds
h_lu        = 1e3*np.load(data_folder+'h'+data_name)    #[J/kg] (200x1 vector)
T_data      = np.load(data_folder+'T'+data_name)        #[K] (lookup table)

#Making T_spline
T_spline    = ip.interp2d(p_lu,h_lu,T_data.T)
func        = lambda h : T_spline(p_in,h)-T_in  #Function to be solved

#Solving for the enthalpy
h_guess     =   100                             #[J/kg]
h_in        =   op.fsolve(func,h_guess)         #[J/kg]

#Model
[m_flow,h,T,p,rho,phase,mu]=PB.varying_mflow(C_D,A_t,D,h_in,p_in,T_s,L_c,prop,eps,PPI)





