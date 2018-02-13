import numpy as np
from Functions import Power_budget as PB
from Functions import Prop_heating as pp
from scipy import interpolate as ip
from scipy import optimize as op
import matplotlib.pyplot as plt

#Validation data
PPI=40          #[inch-1]
eps=0.9         #[-]
p_in=3.5e5        #[Pa]
T_s=11+273.15   #[K]
L_c=0.150       #[m]
D=26e-3         #[m]
Ac=np.pi/4*D**2 #[m2]
prop="H2O"
h_guess=1050e3


#Determining input enthalpy
data_name   = '_data_'+prop+'.npy'
data_folder = 'Functions/Data/'
p_lu        = np.linspace(1e5,10e5,10)  #Pressure bounds
h_lu        = 1e3*np.load(data_folder+'h'+data_name)  #[J/kg] (200x1 vector)
phase_data  = np.load(data_folder+'phase'+data_name)  #[-] (lookup table)
mu_data     = 1e-6*np.load(data_folder+'mu'+data_name)#[Pa s] (lookup table)

mu_spline     = ip.interp2d(p_lu,h_lu,mu_data.T)
phase_spline  = ip.interp2d(p_lu,h_lu,phase_data.T)


#Input Case
v_x   =   [0.18,0.30,0.53]
if PPI==20:
    v_x   =   [0.13,0.22,0.28,0.40,0.51,0.63,0.73]
m =   106

#Initialize
i=0

#Reserving space
dp_dx=np.zeros_like(v_x)
h_val=np.zeros_like(v_x)

for v in v_x:

    #Solver
    func= lambda h : phase_spline(p_in,h)-v

    #Mass flow
    mflow=Ac*m                          #[kg/s]
    h_in=   op.fsolve(func,h_guess)     #[J/kg]
    if (h_in < h_guess+1e3) & (h_in > h_guess-1e3):
        print 'Oops'

    #Model
    [h,T,p,rho,phase,mu] = pp.heating(mflow,D,h_in,p_in,T_s,L_c,prop,eps,PPI)

    #Pressure drop over the channel
    dp_dx[i]=   (p[0]-p[-1])/L_c   #Pressure drop

    #Heating over the channel
    Pgas    =   (h[0]-h[-1])*mflow     #Power towards the gas
    kgas    =   Pgas/(T[0]-T_s)  #Power per Kelvin difference
    h_val[i]    =   kgas/L_c/(np.pi*D)           #Power per Kelvin per tube surface area

    #Next loop
    i=i+1

    print i


print('Pressure drop over the section: ')
print dp_dx
print('Heat transfer over the section: ')
print h_val

print('Start temp')
print T[0]

fname="Validation_v_"+str(PPI)+".txt"
DATA=np.concatenate((h_val,dp_dx))
np.savetxt(fname,(h_val,dp_dx),delimiter="\t")

x=np.linspace(0,1,num=1001)
plt.figure(1)
plt.subplot(221)
plt.plot(x,T,'r')
#plt.axis([0,1,400,420])
plt.title('Temperature')
plt.grid()
plt.subplot(222)
plt.plot(x,p/p_in,'b')
plt.axis([0,1,0.5,1.1])
plt.title('Pressure drop')
plt.grid()
plt.subplot(223)
plt.plot(x,mflow*D/mu/Ac,'c')
plt.axis([0,1,1e5,5e5])
plt.grid()
plt.title('Re')
plt.subplot(224)
plt.plot(x,h/1e3,'g')
plt.axis([0,1,500,2000])
plt.title('Enthalpy')
plt.grid()

plt.show()
