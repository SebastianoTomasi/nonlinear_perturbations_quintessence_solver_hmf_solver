# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 17:29:04 2023

@author: Sebastiano Tomasi
"""

import os
# Get the absolute path of the script Python file
script_dir = os.path.dirname(os.path.abspath(__file__))
# Change the current working directory
os.chdir(script_dir)

import scipy as sp
# import scipy.interpolate
import numpy as np

import mass_function_aux as aux
import sys
sys.path.append("../utility_modules")
import import_export as myie
import numerical_methods as mynm
import cosmological_functions as cosm_func
import plotting_functions as mypl

sys.path.append("../data_modules")
import simulation_parameters as params




if  params.this_run_specifier_1=="perturbed_de":
    this_run_specifier_4="cs2_equal_0"
elif params.this_run_specifier_1=="unperturbed_de":
    this_run_specifier_4="cs2_equal_1"

this_run_specifier_0= "nonlinear_perturbations"



hmfs=[]
legend=[]
for i in range(1,5):#Cycle on the models
    this_run_specifier_2="de_eos_"+str(i)
    print(i)
    for redshift_index in range(1,len(params.z)+1):#cycle on the redshifts
        save_data_to_path="../data/mass_function/perturbed_de/"+this_run_specifier_2+"/"+params.this_run_specifier_3+"/"
        save_data_to_path_i=save_data_to_path+str(redshift_index-1)+"_"+params.fitting_func.lower()+"_hmf"
        legend.append("$w_"+str(i)+"$\t z="+str(params.z[redshift_index-1]))
        data=myie.import_from_txt_multicolumn(save_data_to_path_i)
        
        hmfs.append(data[0])
        
        
#%%
import plotting_functions as mypl

# zoomed_xlim=(2.06*1e12,2.09*1e12)
# zoomed_ylim=(0.00028,0.000290)
# zoomed_position=(0.09, 0.47, #x_pos,y_pos
#                   0.38,#lenght 
#                   0.26)#height

zoomed_xlim=(0.93*1e14,0.947e14)
zoomed_ylim=(1e-9,1.15e-9)
zoomed_position=(0.06, 0.48, #x_pos,y_pos
                  0.38,#lenght 
                  0.26)#height

mypl.plot(hmfs,
            xscale="log",
            yscale="log",
            legend=legend,
            xlabel=["Mass  $[M_\odot h^{-1}]$"],
            ylabel=[r"$\frac{dn}{d\ln (M)}$   $[h^3$ ${\rm Mpc}^{-3}]$"],
            title=None,
            ncol=3,
            zoomed=True, 
            zoomed_xlim=zoomed_xlim, 
            zoomed_ylim=zoomed_ylim, 
            zoomed_position=zoomed_position,
            save=True,
            name="hmf_models_comparison"
            )
















