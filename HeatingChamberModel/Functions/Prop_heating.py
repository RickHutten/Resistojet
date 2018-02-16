import numpy as np
from scipy import interpolate as ip
from scipy import special as sp
import os.path as pt

import ChemicalProperties.reader as reader

"""
Description:
This package models the heating of propellant in the chamber. Supported functions are Darcy's law for pressure drop and convective heating.
The file uses dedicated lookup tables per propellant. These are constructed in the subfolder.

Author: Pim Stohr
Version: 1.0

Critical parameters:
m_flow 	    Mass flow through the chamber
A_flow	    Flow area through the chamber
T_in		Temperature at the inlet
p_in		Pressure at the inlet
rho_in	    Density at the inlet (required if liquid)
D_char      Characteristic diameter of flow through the chamber
T_struc	    Temperature of the structure (constant value)
A_spec	    Specific area (surface area per volume)
L_cham	    Length of the chamber
prop		Type of propellant
perm		Permeability value used in Darcy's law, if no input is given, Darcy's law is ommitted.
alpha       Fraction of flow in gaseous phase
Nu          Nusselt number
"""

# Constants over the functions
C_2f  = 1#2.75                      # Correction parameters for heat transfer in two-phase flow
N     = 1e3                       # Amount of grid points
error = 5e-2                      # Region in which the flow is deemed two-phased

# Relevant number for lookup data
p_lu  = np.linspace(1e5,10e5,10)  # Pressure bounds
# The temperature bounds are for reference only
T_max = 1275                      # Temperature bounds
T_min = 275                       # Temperature bounds
dT    = 5                         # Step in temperature
# The enthalpy is less well defined due to its variation with pressure

# Function which determines the temperature, velocity and pressure in the chamber
def heating(m_flow,D_char,h_in,p_in,T_struc,L_cham,prop,eps,PPI=0,A_s=0,A_c=0,D_strut=0):

  # Flow area in axisymmetric case
  A_flow=eps*np.pi/4*D_char**2    # [m2] Flow area in the chamber
  # print m_flow/A_flow
  # Metal foam properties
  if PPI > 0:
    (perm,A_spec,D_p)=specs(eps,PPI)# [m2] permeability and [m-1]specific area
  else:
    A_spec  = A_s                                 # Input specific area
    perm    = eps**3 /((1-eps)**2) *D_strut**2/180     # Kozeny Carman permeability
    D_p     = 0                                   # Unused variable used to select Nusselt and press drop

   # Overwriting the specific area when input is given
  if A_s > 0:
    A_spec  = A_s

  # Overriding the flow area if there is an input
  if A_c > 0:
    A_flow  = A_c*eps

  # Constants over the iteration
  dx=L_cham/N                     # [m] width of a cell

  # print A_spec,perm,D_p

  # Loading the relevant lookup table
  data_name   = '_data_'+prop+'.npy'
  here        = pt.dirname(pt.abspath(__file__))
  data_folder = pt.join(here,'Data')

  T_data    = np.load(pt.join(data_folder,'T'+data_name))      # [K] (lookup table)
  phase_data= np.load(pt.join(data_folder,'phase'+data_name))  # [-] (lookup table)
  mu_data   = 1e-6*np.load(pt.join(data_folder,'mu'+data_name))# [Pa s] (lookup table)
  k_data    = np.load(pt.join(data_folder,'k'+data_name))      # [W/mK] (lookup table)
  rho_data  = np.load(pt.join(data_folder,'rho'+data_name))    # [kg/m3] (lookup table)
  cp_data   = np.load(pt.join(data_folder,'cp'+data_name))    # [J/kg/K] (lookup table)
  h_lu      = 1e3*np.load(pt.join(data_folder,'h'+data_name))  # [J/kg] (vector)

  # Generating interpolated functions
  T_spline      = ip.interp2d(p_lu,h_lu,T_data.T)
  phase_spline  = ip.interp2d(p_lu,h_lu,phase_data.T)
  mu_spline     = ip.interp2d(p_lu,h_lu,mu_data.T)
  k_spline      = ip.interp2d(p_lu,h_lu,k_data.T)
  rho_spline    = ip.interp2d(p_lu,h_lu,rho_data.T)
  cp_spline     = ip.interp2d(p_lu,h_lu,cp_data.T)

  # Reserving space for outer oriented variables
  T     = np.zeros(int(N+1))
  p     = np.zeros(int(N+1))
  rho   = np.zeros(int(N+1))
  k     = np.zeros(int(N+1))
  phase = np.zeros(int(N+1))
  h     = np.zeros(int(N+1))
  mu    = np.zeros(int(N+1))
  cp    = np.zeros(int(N+1))
  ReDp  = np.zeros(int(N+1))
  Pr    = np.zeros(int(N+1))

  # Inner oriented variables
  h_con = np.zeros(int(N))

  # Initial conditions
  h[0]    = h_in                          # [J/kg]
  p[0]    = p_in                          # [Pa]
  T[0]    = T_spline(p[0],h[0])           # [K]
  k[0]    = k_spline(p[0],h[0])           # [W/m/K]
  mu[0]   = mu_spline(p[0],h[0])          # [Pa s]
  cp[0]   = cp_spline(p[0],h[0])          # [J/kg/K]
  phase[0]= phase_spline(p[0],h[0])       # [-] mass fraction
  rho[0]  = rho_spline(p[0],h[0])         # [kg/m3]
  V       = m_flow/A_flow/rho[0]          # [m/s]
  ReDp[0] = rho[0]*V*D_p/mu[0]            # [-]
  Pr[0]   = cp[0]*mu[0]/k[0]              # [-]


  # Flow-wise solving of the problem
  for i in range(0,int(N)):

    # Checking phase of the flow
    if (phase[i]>0) & (phase[i]<(1-error)):
      C_H     = C_2f                          # Using correction factor of Chen in two-phase flow
    else:
      C_H     = 1                             # Restoring outside two-phase flow

    # Determining which nusselt relation is applicable
    if PPI > 0:
      PeDp    = ReDp[i]*Pr[i]
      Nu      = 2.1234+0.0029*(PeDp/(1-eps))**1.1853   # [-] Low regime estimation
      #Using the high regime equation above transition
      if PeDp > 50:
        Nu = 14.9706*(PeDp/ (1-eps))**(0.2368)# [-] High regime estimation
    else:
      Nu      = 12                            # [-] constant laminar flow

    #  Calculating derivatives at start of cell
    h_con[i]  = Nu/D_char*C_H*A_spec*k[i]       # [W/K/m3]
    qh        = h_con[i]*(T_struc-T[i])         # [W/m3]
    dh        = qh/rho[i] *dx/V                 # [J/kg]

    # Calculate pressure drop
    dp        = darcy(V,perm,mu[i],eps,rho[i],D_char,D_p)# [Pa]
    #print "Velocity:", V,  "Rho:", rho[i], T[i], p[i], h[i]
    # Calculating values at the end of cell
    h[i+1]    = h[i]+dh                         # [J/kg]
    p[i+1]    = p[i]+dp*dx                      # [Pa]

    # Interpolated values from lookup tables
    T[i+1]    = T_spline(p[i+1],h[i+1],)        # [K] Temperature of mixture
    k[i+1]    = k_spline(p[i+1],h[i+1])         # [W/m/K] Conductivity of mixture
    mu[i+1]   = mu_spline(p[i+1],h[i+1])        # [Pa s] Viscosity of mixture
    phase[i+1]= phase_spline(p[i+1],h[i+1])     # [-] mass-fraction: 1 is gas 0 is liquid
    rho[i+1]  = rho_spline(p[i+1],h[i+1])       # [kg/m3] Density of mixture
    #rho[i+1]  = reader.get_density_pressure_temperature(p[i+1], T[i+1])
    cp[i+1]   = cp_spline(p[i+1],h[i+1])        # [J/kg/K] Specific heat of the mixture

    # Determining the speed of the flow
    V         = m_flow/rho[i+1]/A_flow          # Assuming continuum flow [m/s]

    # Determining characteristic numbers
    ReDp[i+1] = m_flow*D_p/mu[i+1]/A_flow       # [-]
    Pr[i+1]   = cp[i+1]*mu[i+1]/k[i+1]          # [-]
    #print rho[i]

  #Returning the output
  return (h,T,p,rho,phase,mu)
  

# Function to determine the pressure drop over porous media through Brinkman's extended Darcy's law
def darcy(V,perm,mu,eps,rho,D,D_p):

  # No extension
  if D_p == 0:
    fp =0


  else:
    #Forchheimer extension
    df_dp   = 1.18*np.sqrt((1-eps)/3/np.pi) /(1-np.exp(-(1-eps)/0.04))
    F       = 0.000212*(1-eps)**(-0.132) * df_dp**-1.63                      #Calmidi
    fp      = rho*F/np.sqrt(perm)*V**2

    # Brinkmann extension
    # Da      = 4*perm/D**2                                                     # [-] Darcy number
    # Da_eps  = np.sqrt(Da/eps)                                                 # [-] Intermediate step for convenience
    # fp      = sp.jv(0,Da_eps)/(2*Da_eps*sp.jv(1,Da_eps) -sp.jv(0,Da_eps))     # [-]

  # Pressure drop
  dp  =   -mu/perm *V - fp

  return dp

# Return the permeability and specific area of the porous medium
def specs(e,PPI):
    # Determining the tortuosity
    C1=(1-np.exp(-(1-e)/0.04))
    C2=(1.18*np.sqrt((1-e)/3/np.pi) /C1)**2
    echi=np.pi*(1-C2)/4                 # Tortuousity (Bhattacharya)
    chi=e/echi

    # Determining cell sizes
    dc=0.0254/PPI                       # CURC size  (definition)
    dp=(3-chi)/2 *dc                    # Pore size  (Du Plessis)
    df=dp*np.sqrt(C2)                   # Fibre size (Calmidi)

        # Specific area of metal foam
    Asf=3*np.pi*df/(0.59*dp)**2 *C1     # Specific area  (Calmidi and Mahajan)

    # Permeability from Calmidi
    perm=0.00073*(1-e)**(-0.224) *(df/dp)**(-1.11) * dp**2  # [m2]

    # Equivalent spherical diameter
    Sv=4*e/dp/(1-e)
    Dp=6/Sv

    # Returning output
    return (perm,Asf,Dp)

