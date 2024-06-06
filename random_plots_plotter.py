# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 21:49:58 2023

@author: sebas
"""

import numpy as np
from numpy import sin
from numpy import cos
from numpy import log
from numpy import exp
from numpy import e
from numpy import pi
from numpy import sqrt
from numpy import sinh
from numpy import cosh
from numpy import arcsinh
from numpy import arccosh

import scipy as sp

import sys
sys.path.append("../data_modules")
sys.path.append("../utility_modules")
sys.path.append("../friedmann_solver")
sys.path.append("../perturbation_solver")

import matplotlib.pyplot as pl

import quintessence_solver 
import numerical_methods as mynm
from simulation_parameters import *
import friedmann_solver as fr_sol

import plotting_functions as mypl
import dark_energy_eos as eos
import lcdm_model
#%%

save=False
omega_rad_now=round(0.2473/(67.4)**2,5)

a=np.linspace(0,3/2*pi,int(1e5))

def omega_m_lcdm_a1(a):
    omega_m_lcdm_a=omega_matter_now/(omega_matter_now+omega_dark_now*a**3)
    return omega_m_lcdm_a

def omega_l_lcdm_a1(a):
    omega_l_lcdm_a=omega_dark_now*a**3/(omega_matter_now+omega_dark_now*a**3)
    return omega_l_lcdm_a 


def omega_m_lcdm_a2(a):
    omega_m_lcdm_a=a*omega_matter_now/(omega_matter_now*a+omega_dark_now*a**4+omega_rad_now)
    return omega_m_lcdm_a

def omega_l_lcdm_a2(a):
    omega_l_lcdm_a=omega_dark_now*a**4/(omega_matter_now*a+omega_dark_now*a**4+omega_rad_now)
    return omega_l_lcdm_a 


def omega_r_lcdm_a(a):
    omega_r_lcdm_a=omega_rad_now/(omega_matter_now*a+omega_dark_now*a**4+omega_rad_now)
    return omega_r_lcdm_a

def zeta(alfa):
    zeta=9/2*(alfa-sin(alfa))**2/(1-cos(alfa))**3
    return zeta

fun1=zeta
fun2=omega_l_lcdm_a1

fun3=omega_m_lcdm_a2
fun4=omega_l_lcdm_a2
fun5=omega_r_lcdm_a

pl1=np.array([a,fun1(a)])
pl2=np.array([a,fun2(a)])
pl3=np.array([a,fun1(a)+fun2(a)])

pl4=np.array([a,fun3(a)])
pl5=np.array([a,fun4(a)])
pl6=np.array([a,fun5(a)])
pl7=np.array([a,fun3(a)+fun4(a)+fun5(a)])

mypl.m_plot(pl1,
          
          xlim=None,
          ylim=None,
          xscale="linear",
          yscale="linear",
          xlabel=r'$\alpha$',
          ylabel="$ \zeta$",
          legend=None,
          title="Matter overdensity in EdS universe",
          name="zeta_overdensity",
          dotted=False,
          y_ticks=[1,5.552,146.85],
          y_ticklables=[1,5.552,146.85],
          x_ticks=[0,pi,3*pi/2],
          x_ticklables=["0",r"$\pi$",r"$\frac{3\pi}{2}$"],
          save=True)


















