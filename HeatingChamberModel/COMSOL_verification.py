import numpy as np
from Functions import Power_budget as pb
from Functions import Nozzle as nz
from Functions import Thermal as th
from Functions import Prop_heating as pp
from scipy import interpolate as ip
from scipy import optimize as op
import matplotlib.pyplot as plt
import time as dt

#Verification data
PPI=40          #[inch-1]
eps=0.95        #[-]
p_in=5e5        #[Pa]
T_in=20+273.15  #[K]
T_s=120+273.15  #[K]
L_c=0.150       #[m]
D=26e-3         #[m]
mflow=1e-3      #[kg/s]
Ac=np.pi/4*(D**2) #[m2]
prop="N2"

#Loading enthalpy data
data_name   = '_data_'+prop+'.npy'                      #Name of data
data_folder = 'Functions/Data/'                         #Name of folder
p_lu        = np.linspace(1e5,10e5,10)                  #Pressure bounds
h_lu        = 1e3*np.load(data_folder+'h'+data_name)    #[J/kg] (200x1 vector)
T_data      = np.load(data_folder+'T'+data_name)        #[K] (lookup table)

#Making T_spline
T_spline    = ip.interp2d(p_lu,h_lu,T_data.T)
func        = lambda h : T_in-T_spline(p_in,h)          #Function to be solved

#Solving for the enthalpy
h_guess     =   h_lu[0]                                     #[J/kg]
h_in        =   op.fsolve(func,h_guess)                 #[J/kg]

[h,T,p,rho,phase,mu]=pp.heating(mflow,D,h_in,p_in,T_s,L_c,prop,eps,PPI)

x=np.linspace(0,1,num=1001)
plt.figure(1)
plt.subplot(221)
plt.plot(x,T,'r')
#plt.axis([0,1,400,420])
plt.title('Temperature')
plt.grid()
plt.subplot(222)
plt.plot(x,p/p_in,'b')
#plt.axis([0,1,0.5,1.1])
plt.title('Pressure drop')
plt.grid()
plt.subplot(223)
plt.plot(x,mflow*D/mu/Ac,'c')
#plt.axis([0,1,1e5,5e5])
plt.grid()
plt.title('Re')
plt.subplot(224)
plt.plot(x,mu*1e5,'g')
#plt.axis([0,1,500,2000])
plt.title('mu')
plt.grid()
print p[-1],T[-1]

#Comsol data
p_com=np.asarray([5,4.999935,4.999866,4.999793,4.999717,4.999638,4.999556,4.999471])
T_com=np.asarray([293.15,303.34,312.38,321.24,329.20,336.36,342.59,348.15])
x_com=np.linspace(0,1,8)

plt.figure()
plt.plot(x,T,'r',label='Python model')
plt.plot(x_com,T_com,'or',label='COMSOL model')
plt.grid()
plt.xlabel('Chamber length [-]')
plt.ylabel('Temperature [K]')
plt.legend(loc=2)

plt.figure()
plt.plot(x,(5e5-p)/5,'b',label='Python model')
plt.plot(x_com,(5-p_com)*1e5,'ob',label='COMSOL model')
plt.grid()
plt.xlabel('Chamber length [-]')
plt.ylabel('Pressure drop [Pa]')
plt.legend(loc=2)

plt.show()