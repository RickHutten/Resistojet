import numpy as np
from Functions import Power_budget as pb
from Functions import Nozzle as nz
from scipy import interpolate as ip
from scipy import optimize as op

import matplotlib.pyplot as plt

#Dimensions of VLM
W_int   =   3e-3
H_int   =   0.15e-3

#Input
eps     =   0.67         #[-] Porosity of the medium
A_s     =   2/W_int +2/H_int #[m-1] Specific heater area
D_strut =   63.2e-6     #[m] Diameter of the struts
p_in    =   5e5         #[Pa] Input pressure
T_in    =   20+273.15   #[K] Input temperature
T_s     =   550  #[K] Structure temperature
L_c     =   9e-3       #[m] Length of the channel
A_c     =   W_int*H_int #[m2] Channel area
D       =   4/A_s       #[m] Characteristic diameter
C_D     =   0.7         #[-] Discharge coefficient
A_t     =   0.15e-3*30.1e-6 /2#[m2] Nozzle throat area
A_e     =   0.5*0.15e-3*500e-6            #[m2] Nozzle exit area
prop    =   "H2O"        #[-] Propellant

print A_t
print D
print A_c
print A_s
print D_strut*1e3
#Loading enthalpy data
data_name   = '_data_'+prop+'.npy'
data_folder = 'Functions/Data/'
p_lu        = np.linspace(1e5,10e5,10)  #Pressure bounds
h_lu        = 1e3*np.load(data_folder+'h'+data_name)  #[J/kg] (200x1 vector)
T_data      = np.load(data_folder+'T'+data_name)  #[K] (lookup table)

#Making T_spline
T_spline    =   ip.interp2d(p_lu,h_lu,T_data.T)
func        =   lambda h : T_spline(p_in,h)-T_in  #Function to be solved

#Solving for input enthalpy
h_guess     =   h_lu[0]
h_in        =   op.fsolve(func,h_guess)         #[J/kg]

#Heating chamber model
[m_flow,h,T,p,rho,phase,mu] =  pb.HCM_temp( C_D, A_t, D,
                                                 h_in, p_in, T_s,
                                                 L_c, prop, eps, A_s=A_s,
                                                 A_c=A_c, D_strut=D_strut)

print m_flow
#Plotting results
x=np.linspace(0,1,1e3+1)
plt.figure()
plt.subplot(221)
plt.plot(x,T,'r')
plt.grid()
plt.title('Temperature')
plt.subplot(222)
plt.plot(x,(p)/p_in,'b')
plt.grid()
plt.title('Pressure')
#plt.axis([0,1,0.9,1 ])
plt.subplot(223)
plt.plot(x,rho,'c')
plt.grid()
plt.title('Density')
plt.subplot(224)
plt.plot(x,mu*1e6,'k')
plt.grid()
plt.title('Viscosity')

#Ideal rocket theory
[F,Isp] =   nz.IRT(p[-1],T[-1],m_flow,A_t,A_e,prop,v_eff=0.85)

print F,Isp

DATA1=np.load('VLM_study.npy')
DATA2 =np.concatenate((DATA1,p.reshape(-1,1),T.reshape(-1,1)),axis=0)
#print DATA.shape
#print p.shape,T.shape

#np.save('VLM_study.npy',DATA2)

plt.show()
