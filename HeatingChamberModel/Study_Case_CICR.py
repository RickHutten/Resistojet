import numpy as np
from Functions import Power_budget as pb
from Functions import Nozzle as nz
from Functions import Thermal as th
from scipy import interpolate as ip
from scipy import optimize as op
import matplotlib.pyplot as plt
import time as dt

#Validation data
PPI     =   50                                          #[inch-1] Pore density of medium
eps     =   0.952                                         #[-] porosity of medium
p_in    =   5e5                                       #[Pa] Input pressure
T_in    =   15+273.15                                   #[K] Input temperature
P_in    =   np.linspace(0,40,41)                                          #[W] Input power
L_c     =   35e-3                                       #[m] Length of chamber
D_ex    =   27e-3                                       #[m] Diameter of the exterior shell
A_ex    =   np.pi*D_ex*L_c+np.pi/4*D_ex **2             #[m2] Exterior chamber area
e_ex    =   0.07                                        #[-] Emissivity of the exterior
D       =   15e-3                                       #[m] Diameter of chamber
C_D     =   0.7                                         #[-] Discharge coefficient
A_t     =   np.pi/4*(0.15e-3)**2                        #[m2] Throat area
A_e     =   np.pi/4*(3e-3)**2                           #[m2] Nozzle exit area
prop    =   "N2"                                        #[-] Propellant

#Thermal analysis of design
r1_steel= 11.5e-3                                       #[m] Inner radius of first shell
r1_PI   = D/2                                           #[m] Inner radius of second shell
r2_steel= D_ex/2                                        #[m] External radius of first shell
k_steel = 16                                            #[W/K/m] conductive constant of steel
k_PI    = 0.16                                          #[W/K/m] conductive constant of Polyimide
k_PE    = 0.40                                          #[W/K/m] conductive constant of Polyethylene
A_m     = 2e-4                                          #[m2] Surface area of the mounting
t_m     = 10e-3                                         #[m] Thickness of the mounting

#Determining conductive constants
k_m     = th.cond_1D(k_PE,A_m,t_m)                      #[W/K] Conductive heat constant of the mounting
k_c1    = th.cond_shell(k_steel,L_c,r2_steel,r1_steel)  #[W/K] Conductive heat constant of first chamber shell
k_c2    = th.cond_shell(k_PI,L_c,r1_PI,r1_steel)        #[W/K] Conductive heat constant of second chamber shell
k_c     = th.cond_comb(k_c1,k_c2)                       #[W/K] Conductive heat constant of the chamber
print "KM", k_m
print A_t
print D
print L_c
print A_ex
print k_m,k_c
print k_c2

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

#Save routine
i  = 0
DATA= np.zeros((41,7))

for P in P_in:
    #Model
    [m_flow,T_s,h,T,p,rho,phase,mu] =    pb.HCM( C_D, A_t, D,
                                             h_in, p_in, P,
                                             k_m, k_c, A_ex, e_ex,
                                             L_c, prop, eps, PPI=PPI)
    eff = (h[-1]-h[0])*m_flow/P
    print T_s,eff
    #Ideal rocket theory
    [F,Isp] =   nz.IRT(p[-1],T[-1],m_flow,A_t,A_e,prop) #[N][s] Vacuum thrust and specific impulse

    #Saving data
    DATA[i,0]   = P
    DATA[i,1]   = F
    DATA[i,2]   = Isp
    DATA[i,3]   = eff
    DATA[i,4]   = T_s
    DATA[i,5]   = m_flow
    DATA[i,6]   = phase[-1]

    #Next loop
    i           = i+1



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
# plt.axis([0,1,0.9,1 ])
plt.subplot(223)
plt.plot(x,phase,'c')
plt.grid()
plt.title('Vapor Quality')
plt.subplot(224)
plt.plot(x,h/1e3,'k')
plt.grid()
plt.title('Enthalpy')

# pin=str(int(p_in/1e5))
# #Saving data on HDD
# fname   = '-CICR_study-' + prop + '-' + pin +'bar.npy'
# stamp   =  dt.strftime('%Y%m%d%H%M%S')
# np.save(stamp+fname,DATA)

plt.show()
