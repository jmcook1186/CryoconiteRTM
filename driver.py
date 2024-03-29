"""
Driver script for CryoconiteRTM

This is where all user-defined vairables are set and output values are plotted

AUTHOR: JOSEPH COOK, April 2020
www.tothepoles.co.uk
www.github.com/jmcook1186

"""

import numpy as np
import math
import matplotlib.pyplot as plt
from ControlFuncs import ControlFuncs
from TwoStreamFuncs import TwoStreamFuncs


########################
# 1 DEFINE HOLE GEOMETRY
########################

hole_d = 10
hole_w = 50
hole_water_d = 5
hole_ar = hole_d/hole_w
# point = hole_w/2 # horizontal distance from LH wall to desired location on hole floor
cryoconite_albedo = np.ones(470)*0.2 #constant albedo across wavelength for now
WL = np.arange(0.3,5,0.01)

####################
## 2. CONFIGURE RTM 
####################

solzen = 15
density = [700]
grain_rds = [700]
layer_type = [1]
dz = [hole_d/100] # cm to m
algae = 0
incoming_i = 4
DIRECT = True
tolerance = 1e-10 #how close to zero doe the flux need to get before we stop iterating internal reflections?

# create named tuple containing snicar input params
params = TwoStreamFuncs.generate_ice_physical_params(density,grain_rds,layer_type,dz,algae,solzen,incoming_i,DIRECT)
incoming = TwoStreamFuncs.generate_incoming_irradiance(params)


BB_output_by_point = []
total_energy_absorbed_by_cryoconite_by_point = []



#############################################################
# END OF USER INPUT (i.e. leave all remaining code unchanged)
#############################################################

# Validation function will raise errors if input data is invalid
ControlFuncs.Validate_Input_Data(hole_d, hole_w, hole_water_d, solzen)

for point in np.arange(0, hole_w, 1):

    # function calls
    dir_energy_absorbed_by_cryoconite, diffuse_energy_absorbed_by_cryoconite, total_incoming_energy, \
        dir_energy_at_hole_floor, diffuse_energy_at_hole_floor = ControlFuncs.CalculateFluxes(\
        hole_d, hole_w, hole_water_d, point, cryoconite_albedo, WL, params, tolerance)


    total_energy_absorbed_by_cryoconite = diffuse_energy_absorbed_by_cryoconite + dir_energy_absorbed_by_cryoconite
    total_energy_absorbed_by_cryoconite_by_point.append(total_energy_absorbed_by_cryoconite)

    BB_output = np.sum(total_energy_absorbed_by_cryoconite)
    BB_output_by_point.append(BB_output)

total_energy_absorbed_by_cryoconite_by_point = np.array(total_energy_absorbed_by_cryoconite_by_point)
total_energy_absorbed_by_cryoconite_by_point = np.mean(total_energy_absorbed_by_cryoconite_by_point,axis=1)
BB_output = np.mean(BB_output_by_point)


# plots and printing
plt.figure()
plt.plot(WL,total_energy_absorbed_by_cryoconite.T,label='total energy absorbed by CC (w/m2)')
plt.plot(WL,incoming,label='total incoming energy (w/m2)')
plt.xlim(0.3,5),plt.xlabel('Wavelength (microns)')
plt.ylabel('Energy (W/m2)')
plt.legend(loc='best')

plt.savefig('/home/joe/Code/CryoconiteRTM/Out.jpg')

print("BROADBAND ENERGY ABSORBED =", np.round(BB_output,3))
print("TOTAL INCOMING ENERGY in WM2 = ", np.round(np.sum(incoming),3))
print("WIDTH/DEPTH RATIO = ", np.round(hole_w/hole_d,3))
print("FLOOR/SURFACE I* RATIO = ", np.round(BB_output/np.sum(incoming),3))
