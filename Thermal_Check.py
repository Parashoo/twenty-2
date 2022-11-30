import numpy as np
import math

#Different input if material changes
max_temp_lug=500
min_temp_lug=350
max_temp_wall=282
min_temp_wall=225
alpha_c_lug=23.6E-6 #Thermal Expansion Coefficient of Lug material

#Values do not change with iteration:
assembly_temp=288 #Assumed temperature during assembly.
alpha_c_wall=23.6E-6    #Thermal Expansion Coefficient of Wall Material (aluminum)

#Values can be assigned from earlier scripts
phi=    #Force Ratio, 4.10
alpha_b=23.6E-6    #Thermal Expansion Coefficient of Fastener Material, 4.10
E_b=    #Young's Modulus of the Fastener, 4.10
A_sm=   #Area on which the Shear force acts

delta_T_max_lug=max_temp_lug-assembly_temp
delta_T_min_lug=min_temp_lug-assembly_temp
delta_T_max_wall=max_temp_wall-assembly_temp
delta_T_min_wall=min_temp_wall-assembly_temp

def induced_load(alpha_c,alpha_b,delta_T,E_b,A_sm,phi):
    F_delta_T=(alpha_c-alpha_b)*delta_T*E_b*A_sm*(1-phi)
    return(F_delta_T)

F_deltaT_max_lug=induced_load(alpha_c_lug,alpha_b,delta_T_max_lug,E_b,A_sm,phi)
F_deltaT_min_lug=induced_load(alpha_c_lug,alpha_b,delta_T_min_lug,E_b,A_sm,phi)
F_deltaT_max_wall=induced_load(alpha_c_wall,alpha_b,delta_T_max_wall,E_b,A_sm,phi)
F_deltaT_min_wall=induced_load(alpha_c_wall,alpha_b,delta_T_min_wall,E_b,A_sm,phi)
