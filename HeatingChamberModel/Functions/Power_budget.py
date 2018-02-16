import numpy as np
from Functions import Prop_heating as heating
from Inputs import Propellants as pp
from scipy import optimize as op

"""
Description:
This file controls the function to solve the power budget for the thruster. In the function, the link is made between the thermal model and the cfd module.

Author:
Pim Stohr

Version:
0.1

Critical parameters:
P_in		Input power
m_flow	    Mass flow
p_in		Inlet pressure
prop		Propellant
k_m		    Conductive heat loss parameter for the mounting (W/K)
k_c		    Conductive heat loss parameter for the chamber (W/K)
A_ex		Exterior area of the chamber
e_ex		Emmissivity of the chamber exterior
"""

#Constant over all functions
T_inf=288.15        #[K] Exterior temperature
sigma=5.67037321e-8 #[W/K4/m2] Stefan constant

#Function to solve the heating problem with varying mass flow (and constant temperature)
def HCM_temp(C_D,A_t,D_char,h_in,p_in,T_struc,L_cham,prop,eps,PPI=0,A_s=0,A_c=0,D_strut=0):
  #Finding propellant properties
  [R,gam]=pp.properties(prop)
  print R

  #Vandekerckhove constant
  VdK = np.sqrt(gam*(2/(gam+1))**((gam+1)/(gam-1)))

  #Initial conditions
  m_flow  = C_D*A_t*p_in*VdK/np.sqrt(T_struc*R)  #[kg/s] Assuming perfect heating
  error   = 1       #[-] Initial error to start loop

  #Tolerance for break
  tol     = 1e-3    #[-]

  #Counter to correct overflow
  N_max   = 1e2     #[-] Counter max value
  N       = 0       #[-] Counter initial value

  while (error > tol)&(N<N_max):
      #Simulating flow through the chamber
      [h,T,p,rho,phase,mu]=heating.heating(m_flow,D_char,h_in,p_in,T_struc,L_cham,prop,eps,PPI,A_s,A_c,D_strut)

      #Updating the massflow
      m_old   = m_flow
      m_flow  = C_D*A_t*p[-1]*VdK/np.sqrt(T[-1]*R)
      error   = np.abs((m_flow-m_old)/m_old)

      #Fail safe to correct when spitting occurs
      if phase[-1] < 1:
        N = N_max
        print('Loop broken due to spitting')

      #Increasing the counter
      N = N + 1

      #Providing feedback
      print 'Finding steady state...'

  #Returning output
  return (m_flow,h,T,p,rho,phase,mu)

#Function to solve the heating problem at varying structure temperatures (and constant mass flow)
def HCM_mflow(m_flow,D_char,h_in,p_in,P_in,k_m,k_c,A_ex,e_ex,L_cham,prop,eps,PPI=0,A_s=0,A_c=0,D_strut=0):
    #Initial conditions
    eff     =   0.5         #[-] Heating efficiency initial guess
    error   =   1           #[-] Initial error to start loop
    tol     =   0.1         #[-] Tolerance for break
    P_loss  =   (1-eff)*P_in#[W] Initial guess for power lost

    T_struc =   Temp_in(k_m,k_c,A_ex,e_ex,P_loss)

    #Counter to correct overflow
    N_max   = 50.           #[-] Counter max value
    N       = 0             #[-] Counter initial value

    while (error > tol)&(N<N_max):
        #Simulating flow through the chamber
        [h,T,p,rho,phase,mu]=heating.heating(m_flow,D_char,h_in,p_in,T_struc,L_cham,prop,eps,PPI,A_s,A_c,D_strut)

        #Updating the temperature of the structure
        P_gas   =   (h[-1]-h[0])*m_flow     #[W] Power towards propellant
        eff     =   P_gas/P_in              #[-] Efficiency value
        P_loss  =   (1-eff)*P_in            #[W] Heat lost

        #Updating the temperature
        T_old   =   T_struc                 #Saving old value

        if P_loss < 0:  # Not possible to
            P_loss = 0

        # Calculate new structure temperature with smoothing, needed for low temperatures
        T_struc =   T_old * (N/N_max) + (1-N/N_max)*Temp_in(k_m,k_c,A_ex,e_ex,P_loss)
        # T_struc = Temp_in(k_m,k_c,A_ex,e_ex,P_loss)
        print T_struc
        error   =   np.abs((T_struc-T_old))   # Determining absolute error

        #Updating overflow parameter
        N       =   N+1

    return (eff, P_in, T_struc,h,T,p)

#Function to solve the heating problem at varying structure temperatures and varying mass flow)
def HCM(C_D,A_t,D_char,h_in,p_in,P_in,k_m,k_c,A_ex,e_ex,L_cham,prop,eps,PPI=0,A_s=0,A_c=0,D_strut=0):
    #Finding propellant properties
    [R,gam]=pp.properties(prop)

    #Vandekerckhove constant
    VdK = np.sqrt(gam*(2/(gam+1))**((gam+1)/(gam-1)))

    #Initial conditions
    eff     =   0.50        #[-] Heating efficiency initial guess
    error   =   1           #[-] Initial error to start loop
    tol     =   1e-3        #[-] Tolerance for break
    P_loss  =   (1-eff)*P_in#[W] Initial guess for power lost

    T_struc =   Temp_in(k_m,k_c,A_ex,e_ex,P_loss)
    m_flow  =   C_D*A_t*p_in*VdK/np.sqrt(T_struc*R)  #[kg/s] Assuming perfect heating

    #Counter to correct overflow
    N_max   = 1e2           #[-] Counter max value
    N       = 0             #[-] Counter initial value

    while (error > tol)&(N<N_max):
        #Simulating flow through the chamber
        [h,T,p,rho,phase,mu]=heating.heating(m_flow,D_char,h_in,p_in,T_struc,L_cham,prop,eps,PPI,A_s,A_c,D_strut)

        #Updating the massflow
        m_old   = m_flow                                #[kg/s] Saving old value
        m_flow  = C_D*A_t*p[-1]*VdK/np.sqrt(T[-1]*R)        #[kg/s] New massflow
        error2  = np.abs((m_flow-m_old)/m_old)          #[-] Determining error

        #Updating the efficiency of the thruster
        P_gas   =   (h[-1]-h[0])*m_flow     #[W] Power towards propellant
        eff     =   P_gas/P_in              #[-] Efficiency value
        P_loss_new  =   (P_in-P_gas)            #[W] Heat lost
        P_loss  = 0.8*P_loss + 0.2*P_loss_new #[W] Updating Ploss, with a smoothing

        #Updating the temperature
        T_old   =   T_struc                             #[K] Saving old value
        T_struc =   Temp_in(k_m,k_c,A_ex,e_ex,P_loss)   #[K] New structure temperature
        error1  =   np.abs((T_struc-T_old)/T_old)       #[-] Determining error

        #Steady state is reached when both values converge
        error   = max(error1,error2)

        #Updating overflow parameter
        N       =   N+1

        #Providing feedback
        print 'Finding steady state...'

        #Fail safe to correct when spitting occurs
        if phase[-1] < 1-1e-5: #Small deviation due to machine errors
            N = N_max
            print('Loop broken due to spitting')
            print('T at nozzle: ',T[-1])
            print('Input power: ',P_in)


    return (m_flow,T_struc,h,T,p,rho,phase,mu)

#Function to solve the interior temperature at a power level
def Temp_in(k_m,k_c,A_ex,e_ex,P_loss):
    #Saving parameters in named
    C1  =   k_m             #[W/K] Conductive heat loss of exterior
    C2  =   sigma*A_ex*e_ex #[W/K4] Radiative heat loss of exterior
    C3  =   k_c             #[W/K] Conductive heat loss of interior

    #Determing exterior temperature
    T_ex    =   Temp_ext(C1,C2,P_loss)  #[K] Exterior temperature

    #Determing interior temperature
    T_in   =   P_loss/C3 + T_ex

    #Returing the interior temperature
    return T_in

#Function to solve the exterior temperature at a power level
def Temp_ext(C1,C2,P_loss):
  #Startpoint of the solver
  T_guess   = 200+273.15

  #Result
  if C1 == 0:
      T_ex = ((C2*T_inf**4 + P_loss) / C2)**(1/4.)
  else:
      # Function to solve
      func = lambda T: C1 * (T - T_inf) + C2 * (T ** 4 - T_inf ** 4) - P_loss
      T_ex  = op.fsolve(func,T_guess)

  #Returning the required value
  return T_ex