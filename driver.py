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
from TwoStreamFuncs import TwoStreamFuncs


########################
# 1 DEFINE HOLE GEOMETRY
########################

hole_d = 40
hole_w = 40
hole_water_d = 30
hole_ar = hole_d/hole_w
point = hole_w/2 # horizontal distance from LH wall to desired location on hole floor
cryoconite_albedo = np.ones(470)*0.2 #constant albedo across wavelength for now
WL = np.arange(0.3,5,0.01)

####################
## 2. CONFIGURE RTM 
####################

solzen = 1
density = [600]
grain_rds = [500]
layer_type = [1]
dz = [hole_d/100] # cm to m
algae = 0
incoming_i = 4
DIRECT = True

# density,grain_rds,layer_type,dz,algae,solzen, incoming_i, DIRECT
params = TwoStreamFuncs.generate_ice_physical_params(density,grain_rds,layer_type,dz,algae,solzen,incoming_i,DIRECT)
incoming = TwoStreamFuncs.generate_incoming_irradiance(params)


#############################################################
# END OF USER INPUT (i.e. leave all remaining code unchanged)
#############################################################

dir_energy_absorbed_by_cryoconite, diffuse_energy_absorbed_by_cryoconite, total_incoming_energy, \
    dir_energy_at_hole_floor, diffuse_energy_at_hole_floor = ControlFuncs.CalculateFluxes(\
    hole_d, hole_w, hole_water_d, point, cryoconite_albedo, WL, params)

total_energy_absorbed_by_cryoconite = diffuse_energy_absorbed_by_cryoconite + dir_energy_absorbed_by_cryoconite

BB_output = np.sum(total_energy_absorbed_by_cryoconite)


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