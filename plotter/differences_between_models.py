# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 13:03:57 2023

@author: Sebastiano Tomasi
"""

import scipy as sp
import numpy as np
import sys

import os
# Get the absolute path of the script Python file
script_dir = os.path.dirname(os.path.abspath(__file__))
# Change the current working directory
os.chdir(script_dir)

import matplotlib.pyplot as pl


sys.path.append("../utility_modules")
import lcdm_model as lcdm
import import_export as myie
import cosmological_functions as cosm_func
import plotting_functions as mypl

save = True
format_=".pdf"#".pdf"

default_save_plot_dir="../data/defoult_plot_folder/"

this_run_specifier_1 = "perturbed_de"
this_run_specifier_3 = "w_i"


delta_m=[]
linear_matter_density_contrasts_at_collapse_z=[]
for this_run_specifier_2 in ["de_eos_1","de_eos_2","de_eos_3","de_eos_4"]:
    """Import the LCDM to compare evrithing else aganist."""
    lcdm_linear_matter_density_contrasts_at_collapse_z = myie.import_from_txt_multicolumn(
        "../data/nonlinear_perturbations/LCDM/delta_c")
    
        
    get_data_from_path = "../data/nonlinear_perturbations/" + \
        this_run_specifier_1+"/"+this_run_specifier_2+"/"+this_run_specifier_3+"/"
    
    try:
        # if this_run_specifier_2=="de_eos_1":
        linear_matter_density_contrasts_at_collapse_z.append(myie.import_from_txt_multicolumn(
            get_data_from_path+"delta_c")[0])
        delta_m.append(myie.import_from_txt_multicolumn(get_data_from_path+"growth_factor")[0])
        # else:
        #     linear_matter_density_contrasts_at_collapse_z.append(myie.import_from_txt_multicolumn(
        #         get_data_from_path+"delta_c"))
        #     delta_m.append(myie.import_from_txt_multicolumn(get_data_from_path+"growth_factor"))
            
    except FileNotFoundError:
        print(f"Error: file not found. Check if the data is available at {get_data_from_path}")
    

    
    legend=["$w_1$","$w_2$","$w_3$","$w_4$"]
    legend.append("$\Lambda$CDM")
    
    
    



"""Append the lcdm model curves"""
linear_matter_density_contrasts_at_collapse_z.append(
    lcdm_linear_matter_density_contrasts_at_collapse_z[0])

"""Use the lcdm_model to do comparisons"""
LCDM = lcdm.LCDM()
linear_perturbations_lcdm_numerical_a = LCDM.compute_perturbations()
lcdm_scale_parameter = linear_perturbations_lcdm_numerical_a[0]
linear_perturbations_lcdm_a = sp.interpolate.interp1d(lcdm_scale_parameter,
                                                      linear_perturbations_lcdm_numerical_a[1],
                                                          fill_value="extrapolate", assume_sorted=False)
#%%

mypl.plot(delta_m, r"$a$", r"$\delta_{\rm{m}}$",
            title="",
            legend=legend,
            func_to_compare=linear_perturbations_lcdm_a,
            save=save,
            # save_plot_dir=save_plot_to_path, 
            format_=format_,
            dotted=False,
            zoomed=True, 
            zoomed_xlim=(0.932,0.96),
            # zoomed_xticks=(0.9,1),
            zoomed_ylim=(0.7,0.71),
            # ùzoomed_yticks=(0.92,0.96), 
            zoomed_position=(0.5, 0.1, #x_pos,y_pos
                              0.45,#lenght 
                              0.4),#height
            name=default_save_plot_dir+"density_contrast_matter_models_a",
            # xlim=(0.9,1),ylim=(0.7,0.85)
            )


        
mypl.plot(linear_matter_density_contrasts_at_collapse_z, legend=legend,
            # title="Linearly extrapolated matter density contrast at collapse.",
            title="",
            xlabel=r"$z_{\rm{c}}$", ylabel=r"$\delta_{\rm{c}}$",
            zoomed=True, 
            zoomed_xlim=(8,8.5),
            # zoomed_xticks=(0.9,1),
            zoomed_ylim=(1.685785,1.68586),
            # ùzoomed_yticks=(0.92,0.96), 
            zoomed_position=(0.2, 0.1, #x_pos,y_pos
                              0.52,#lenght 
                              0.6),#height
            save=save,
            name=default_save_plot_dir+"delta_models_c",
            # save_plot_dir=save_plot_to_path,
            format_=format_)