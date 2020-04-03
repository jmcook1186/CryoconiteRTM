"""
Driver script for CryoconiteRTM

This is where all user-defined vairables are set and output values are plotted

AUTHOR: JOSEPH COOK, April 2020
www.tothepoles.co.uk
ww.github.com/jmcook1186

"""

import numpy as np
import math
import matplotlib.pyplot as plt
from ControlFuncs import ControlFuncs


########################
# 1 DEFINE HOLE GEOMETRY
########################

hole_d = 7
hole_w = 4
hole_water_d = 3
hole_ar = hole_d/hole_w
point = hole_w/2 # horizontal distance from LH wall to desired location on hole floor
cryoconite_albedo = np.ones(470)*0.2 #constant albedo across wavelength for now
WL = np.arange(0.3,5,0.01)

Mie = True
GeometricOptics = False

####################
## 2. CONFIGURE RTM 
####################

SZA = 20 # SZA is calculated as degrees from the vertical (zenith) - i.e. 0 is illumination from directly overhead
incoming = np.genfromtxt ('/home/joe/Code/CryoconiteRTM/Data/mlw_sfc_flx_frc_clr.txt', delimiter=",") # path to I* file

#############################################
## 3. SET PHYSICAL PROPERTIES OF THE ICE/SNOW
#############################################

rho_snw = [600] # density of each layer (unit = kg m-3)
# if using Mie optical properties, set spherical grain radius
rds_snw = [1500]
# if using GeometricOptics, set grain side_length and depth
side_length = [15000] 
depth = [15000]


#############################################################
# END OF USER INPUT (i.e. leave all remaining code unchanged)
#############################################################

dir_energy_absorbed_by_cryoconite, diffuse_energy_absorbed_by_cryoconite, total_incoming_energy, \
    dir_energy_at_hole_floor, diffuse_energy_at_hole_floor = ControlFuncs.CalculateFluxes(\
    hole_d, hole_w, hole_water_d, point, cryoconite_albedo, WL, Mie, GeometricOptics, SZA, incoming,\
    rho_snw, rds_snw, side_length, depth)


total_energy_absorbed_by_cryoconite = diffuse_energy_absorbed_by_cryoconite + dir_energy_absorbed_by_cryoconite

BB_output = np.sum(total_energy_absorbed_by_cryoconite, axis=1)


plt.figure()
plt.plot(WL,total_energy_absorbed_by_cryoconite.T,label='total energy absorbed by CC (w/m2)')
plt.plot(WL,incoming,label='total incoming energy (w/m2)')
plt.xlim(0.3,5),plt.xlabel('Wavelength (microns)')
plt.ylabel('Energy (W/m2)')
plt.legend(loc='best')

plt.show()

print("BROADBAND ENERGY ABSORBED =", BB_output)
print("TOTAL ENERGY in WM2 = ", np.sum(incoming))
print("WIDTH/DEPTH RATIO = ", hole_w/hole_d)
print("FLOOR/SURFACE I* RATIO = ", BB_output/np.sum(incoming))