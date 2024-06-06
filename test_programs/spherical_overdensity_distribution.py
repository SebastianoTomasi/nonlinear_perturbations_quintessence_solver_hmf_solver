# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 11:45:03 2023

@author: Sebastiano
"""

import numpy as np
from numpy import sqrt
from numpy import pi
from numpy import e

import matplotlib.pyplot as plt
import matplotlib

N=50
L=1
spatial_dimensions=2
V=L**spatial_dimensions

R_medio=0.01
sigma_R=R_medio*0.5

def center_probability_density(N,V):
    """Returns the probability that a spherical overdensity is contained in
    the volume dV if the universe has volume V and N overdensities in it."""
    return N/V

def radius_probability_density(Ri,R_medio,sigma_R):
    """Returns the probability density that a spherical overdensity has
    radius Ri"""
    return 1/(sigma_R*sqrt(2*pi))*e**(-0.5*((Ri-R_medio)/sigma_R)**2)

def generate_radius():
    return np.random.normal(loc=R_medio,scale=sigma_R,size=None)

def generate_center():
    return np.random.uniform(low=0,high=L,size=spatial_dimensions)

def generate_overdensities():
    overdensities=[]
    for i in range(N):
        overdensities.append([generate_center(),generate_radius()])
    return overdensities

overdensities=generate_overdensities()


# Note that the patches won't be added to the axes, instead a collection will
patches = [plt.Circle(center, size) for center, size in overdensities]

fig, ax = plt.subplots(figsize=(10,15),dpi=256)


coll = matplotlib.collections.PatchCollection(patches, facecolors='black')
ax.add_collection(coll)
ax.set_aspect(aspect=1)

ax.margins(0.01)
plt.xlim((0,L))
plt.ylim((0,L))
plt.show()
