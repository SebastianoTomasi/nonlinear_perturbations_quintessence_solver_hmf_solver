# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 15:32:54 2023

@author: Sebastiano Tomasi
"""
import numpy as np
from numpy import sqrt

import scipy as sp
import scipy.interpolate

import os
# Get the absolute path of the script Python file
script_dir = os.path.dirname(os.path.abspath(__file__))
# Change the current working directory
os.chdir(script_dir)

import sys
sys.path.append("../data_modules")
from data_classes import friedmann_sol
import cosmological_functions as cosm_func
import simulation_parameters as params

sys.path.append("../utility_modules")
import numerical_methods as mynm
import plotting_functions as mypl

#%% Global variables

asterisks_lenght=int(60)

"""Create the class object to store the results"""
result=friedmann_sol()


#%%

def solve(time_domain=False,does_print_info=False):
    """This function solve for the background equations, given the parameters in the simulation_parameters module.
        input:
            - time_domain if True, computes all the bakground quantities also as function of time. 
            - doese_print_info if True prints the age of the universe and the deceleration parameter."""
            
    """Compute the dark density evolution"""
    dark_density_evolution_numerical_a=cosm_func.solve_dark_density_evolution_numerical_a()
    scale_parameter_values=dark_density_evolution_numerical_a[0]
    dark_density_evolution_a=sp.interpolate.interp1d(scale_parameter_values, dark_density_evolution_numerical_a[1],
                                        fill_value="extrapolate", assume_sorted=True)
    
    
    rescaled_hubble_function_a=cosm_func.create_rescaled_hubble_function_a(dark_density_evolution_a)
    
    """Compute also the de_eos numerically in order to plot it."""
    de_eos_numerical_a=np.array([scale_parameter_values,params.de_eos_a(scale_parameter_values)])
    
    
    """Effective eos"""
    result.effective_eos_numerical_a=[scale_parameter_values,
                                      (params.de_eos_a(scale_parameter_values)*dark_density_evolution_numerical_a[1]+  \
                                       cosm_func.rad_density_evolution_a(scale_parameter_values)/3)/    \
                                          (cosm_func.matter_density_evolution_a(scale_parameter_values)+   \
                                         dark_density_evolution_numerical_a[1]+   \
                                        cosm_func.rad_density_evolution_a(scale_parameter_values))]
    
    """Approximate dark density evolution using a step transition."""
    aux=[]
    for i in scale_parameter_values:
        aux.append(cosm_func.step_dark_density_evolution_a(i))
    appox_dark_density_evolution_numerical_a=np.array([scale_parameter_values,aux])
    
    """Rescaled hubble function H/H0, we need it to compute the density parameters"""
    rescaled_hubble_functions_numerical_a=np.array([scale_parameter_values,
                                          rescaled_hubble_function_a(scale_parameter_values)])
    
    dark_density_parameter_numerical_a=np.array([scale_parameter_values,
              dark_density_evolution_a(scale_parameter_values)/rescaled_hubble_functions_numerical_a[1]**2])
    matter_density_parameter_numerical_a=np.array([scale_parameter_values,
              cosm_func.matter_density_evolution_a(scale_parameter_values)/rescaled_hubble_functions_numerical_a[1]**2])
    rad_density_numerical_a=np.array([scale_parameter_values,
              cosm_func.rad_density_evolution_a(scale_parameter_values)/rescaled_hubble_functions_numerical_a[1]**2])
    
    """Saving result in the Friedmann results class"""
    result.dark_density_evolution_numerical_a=dark_density_evolution_numerical_a
    result.appox_dark_density_evolution_numerical_a=appox_dark_density_evolution_numerical_a
    result.rescaled_hubble_functions_numerical_a=rescaled_hubble_functions_numerical_a
    result.dark_density_parameter_numerical_a=dark_density_parameter_numerical_a
    result.matter_density_parameter_numerical_a=matter_density_parameter_numerical_a
    result.rad_density_numerical_a=rad_density_numerical_a
    result.de_eos_numerical_a=de_eos_numerical_a
    
    if time_domain:
        """Integrate the friedmann equation"""
        friedmann_equation_integrand=cosm_func.create_log_friedmann_equation_integrand(dark_density_evolution_a)
        hubble_constant_times_t = mynm.integrate(f=friedmann_equation_integrand,
                                                  a=np.log(params.a_min), b=np.log(params.a_max),
                                                  atol=params.friedmann_atol,rtol=params.friedmann_rtol,max_step=params.friedmann_max_stepsize)
        hubble_constant_times_t[0]=np.exp(hubble_constant_times_t[0])#It is y=ln(a)
        scale_parameter_values=hubble_constant_times_t[0]
    
        """Compute t(a)==time_a"""
        time=hubble_constant_times_t[1]/params.hubble_constant
        universe_age = time[-1]#In giga yeras
        
        """Invert t(a) to obtain a(t). In more detail we have a(t) s.t a(t0)=1 
        where t_0 is the age of the universe"""
        scale_parameter_numerical_t=np.array([time,scale_parameter_values])
        
        """Compute the hubble function H(t)=\dot{a}/a and the second derivative of a"""
        scale_parameter_t=sp.interpolate.interp1d(time, scale_parameter_values,
                                          fill_value="extrapolate", assume_sorted=True,kind="quadratic")
        h=mynm.rms_neigbour_distance(time)#Compute average spacing between time points
        new_lenght=int(2/h)#Use h to obtain the new time array lenght
        time=np.linspace(time[0],time[-1],new_lenght)#Build a new time array with equally spaced points
        scale_parameter_values=scale_parameter_t(time)#Compute the corresponding scale parameter values
        
        

        
        scale_parameter_derivative_numerical_t=mynm.derivate(scale_parameter_t,time[0],time[-1], new_lenght)
        
        hubble_function_numerical_t=np.array([time,scale_parameter_derivative_numerical_t[1]/scale_parameter_values])
    
        """We can now calculate the deceleration parameter today: q_0"""
        scale_parameter_derivative_t=sp.interpolate.interp1d(scale_parameter_derivative_numerical_t[0], 
                                                            scale_parameter_derivative_numerical_t[1],
                                            fill_value="extrapolate", assume_sorted=True,kind="quadratic")
        scale_parameter_2derivative_numerical_t=mynm.derivate(scale_parameter_derivative_t,time[0],time[-1], new_lenght)
        
        dec_param_now=-scale_parameter_2derivative_numerical_t[1][-1]/params.hubble_constant**2

        """Save results"""
        result.scale_parameter_numerical_t=scale_parameter_numerical_t
        result.hubble_function_numerical_t=hubble_function_numerical_t
        result.universe_age=universe_age
        result.deceleration_parameter=dec_param_now
        
        """Compare the age of the universe of the model to LCDM"""
        if does_print_info:
            if params.omega_matter_now==1:
                theoretical_universe_age=2/(3*params.hubble_constant)
            else:
                theoretical_universe_age= cosm_func.time_a_LCDM(1)
            print("*"*asterisks_lenght)
            print("MODEL: t0=",round(universe_age,5)," Gy\t\tLCDM: t0=",round(theoretical_universe_age,5)," Gy")
            print("*"*asterisks_lenght)
            print("MODEL: q0 = ",round(dec_param_now,3),"""\t\t\tLCDM: q0 = """,round(params.omega_matter_now/2-params.omega_dark_now,3))
            print("*"*asterisks_lenght)
    return result





# backgr_dat=np.loadtxt("D:/Synced/OneDrive/python_code/data/class/de_eos_1/w_i/cs2_equal_0/0_background.dat",comments="#")
mine=solve(False,True)
#%%
from numpy import pi,log,sqrt,exp,tanh


omega0 = params.omega_matter_now  # Example value, replace with actual value
omegaL = params.omega_dark_now  # Example value, replace with actual value
gamma = params.trans_steepness  # Example value, replace with actual value
wf = params.w_f  # Example value, replace with actual value
wi = params.w_i  # Example value, replace with actual value
zt = params.trans_z  # Example value, replace with actual value

def anal_h(a):
    term1 = omega0 / pow(a, 3)
    term2_1 = pow(a, -3 * (1 + (wf + wi) / 2))
    term2_2 = omegaL
    term2_3 = (pow(a, -gamma) * pow(1 + zt, -gamma) + pow(a, gamma) * pow(1 + zt, gamma)) / (pow(1 + zt, -gamma) + pow(1 + zt, gamma))
    term2_4 = pow(term2_3, (3 * (-wf + wi)) / (2 * gamma))
    term2 = term2_1 * term2_2 * term2_4
    
    result = sqrt(term1 + term2)
    return result

    

def anal_ahhp(a):
    term1 = 1.0 / (2 * pow(a, 3))
    term2 = -omega0
    term3 = pow(a, -3 * wi)
    term4 = omegaL
    term5 = (1 + pow(1 + zt, 2 * gamma)) / (1 + pow(a * (1 + zt), 2 * gamma))
    term6 = pow(term5, (3 * (wf - wi)) / (2 * gamma))
    term7 = 1 + wi + (1 + wf) * pow(a * (1 + zt), 2 * gamma)
    term8 = 1 + pow(a * (1 + zt), 2 * gamma)
    
    result = (term1 * 3 * (term2 - (term3 * term4 * term6 * term7) / term8))
    return result


num_h=mine.rescaled_hubble_functions_numerical_a

anal=[[],[]]

anal[0]=num_h[0]
anal[1]=anal_h(num_h[0])


mypl.plot([num_h,anal],legend=["num","anal"],
           xscale="log",
           yscale="log"
          )
# mypl.plot(num_h,legend=["num"],
#            xscale="log",
#            yscale="log"
#           )

# mypl.plot(anal,legend=["anal"],
#            xscale="log",
#            yscale="log"
#           )


num_h[1]=np.abs(num_h[0]*num_h[1]*mynm.Nderivate(num_h)[1])

anal=[[],[]]

anal[0]=num_h[0]
anal[1]=np.abs(anal_ahhp(num_h[0]))


mypl.plot([num_h,anal],legend=["num1","anal1"],
            xscale="log",
            yscale="log",
            dotted=True,
          )
mypl.plot(num_h,legend=["num"],
            xscale="log",
            yscale="log"
          )

mypl.plot(anal,legend=["anal"],
            xscale="log",
            yscale="log"
          )

# a=1/(backgr_dat[:,0]+1)
# rho_fld=backgr_dat[:,11]
# compare=[a,rho_fld]
# compare=[a,rho_fld/rho_fld[-1]]


# my_a=mine.dark_density_evolution_numerical_a[0]
# my_rho_fld=mine.dark_density_evolution_numerical_a[1]/mine.dark_density_evolution_numerical_a[1][-1]

# # my_a=mine.de_eos_numerical_a[0]
# # my_rho_fld=mine.de_eos_numerical_a[1]


# my_rho_fld=[my_a,my_rho_fld]

# mypl.plot([compare,my_rho_fld],
#            yscale="log",
#            xscale="log",
#           legend=["class","mine"])
# mypl.plot(compare)



