import numpy as np
from Inputs import Propellants as pp
from scipy import optimize as op

"""
Description:
This module holds function that are used when calculating with a set nozzle. The functions are used to calculate a nozzle constant as well as give an estimate of the discharge coefficient

Author:
Pim Stohr

Version:
1.0

Critical variables:
p_c		Pressure in chamber
T_c		Temperature in chamber
m_flow	Mass flow through the chamber
prop	Propellant that is used
gam		Specific heat ratio
R		Universal gas constant
VdK		Vanderkerckhove constant
A_t		Throat area
A_e		Exit area
"""

#Function to determine propulsion propeties from chamber properties
def IRT(p_c,T_c,m_flow,A_t,A_e,prop,v_eff=1):
	#Importing propellant properties
	[R,gam]=pp.properties(prop)

	#Pressure ratio in the exit
	p_r	=	p_r_solve(A_e/A_t,gam)				#[-] Pressure ratio between exit and chamber
	p_e	=	p_r*p_c								#[Pa] Pressure at the exit of nozzle

	#Limit velocity
	U_l	=	np.sqrt(2*gam/(gam-1)*R*T_c)		#[m/s] Limit velocity
	U_e	=	U_l*np.sqrt(1-p_r**((gam-1)/gam))	#[m/s] Exit velocity

	#Accounting for nozzle efficiency
	U_e =	U_e*v_eff

	#Vaccuum thrust
	F	=	U_e*m_flow+A_e*p_e					#[N] Vacuum thrust

	#Corresponding specific impulse
	Isp	=	F/m_flow/9.81						#[s] Vacuum specific impulse

	#Returning output
	return (F,Isp)

#Function to solve the pressure ratio for area ratio
def p_r_solve(A_r,g):
	#Vandekerckhove constant
  	VdK = np.sqrt(g*(2/(g+1))**((g+1)/(g-1)))

	#Constants for convenience
	C1	=	2/(g-1)
	C2	=	2/g
	C3	=	(g-1)/g

	#Setting up function to solve
	func	=	lambda p : VdK/np.sqrt(C1*p**C2 *(1-p**C3)) - A_r	#Area ratio function

	#Solving for p_r
	p_guess	=	1/A_r/1e3							#[Pa] Initial guess for the solver
	p_r		=	op.fsolve(func,p_guess)		#[-] Solved for the root value

	#Returning output
	return p_r