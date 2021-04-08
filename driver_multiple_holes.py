"""
Driver script for CryoconiteRTM

This is where all user-defined vairables are set and output values are plotted

AUTHOR: JOSEPH COOK, April 2020
www.tothepoles.co.uk
www.github.com/jmcook1186

TODOs

1) calculate albedo
2) add functions for 2D surface ( % coverage by holes, dimensions, total albedo)

"""

import numpy as np
import math
import matplotlib.pyplot as plt
from ControlFuncs import ControlFuncs
from TwoStreamFuncs import TwoStreamFuncs

########################
# 1 DEFINE HOLE GEOMETRY
########################

hole_d = [10, 15, 20, 20] # depth of each hole
hole_w = [20, 30, 40, 50] # width of each hole
hole_water_d = [5, 10, 10, 20] #water depth in each hole
n_holes = [2, 1, 2, 5] # how many of each hole?
hole_areas= np.pi*((np.array(hole_w)/200)**2) # in m
study_area = 1 # total study area in m2

if np.sum(hole_areas>study_area):
    raise ValueError (f"Cryconite area = {np.sum(hole_areas)} Total study area is less than total cryocontie area")

cryoconite_albedo = np.ones(470)*0.2 #constant albedo across wavelength for now
WL = np.arange(0.3,5,0.01)

####################
## 2. CONFIGURE RTM 
####################

solzen = 25
density = [700]
grain_rds = [700]
layer_type = [1]
dz = [np.max(hole_d)/100] # cm to m
algae = 0
incoming_i = 4
DIRECT = True
tolerance = 1e-10 #how close to zero doe the flux need to get before we stop iterating internal reflections?

# create named tuple containing snicar input params
params = TwoStreamFuncs.generate_ice_physical_params(density,grain_rds,layer_type,dz,algae,solzen,incoming_i,DIRECT)
incoming = TwoStreamFuncs.generate_incoming_irradiance(params)

BB_out_all_holes = np.zeros(shape=(len(hole_w),1))
absolute_energy_absorbed_all_holes = np.zeros(shape=(len(hole_w),1))
#############################################################
# END OF USER INPUT (i.e. leave all remaining code unchanged)
#############################################################

# Validation function will raise errors if input data is invalid
total_energy_absorbed_per_hole = []

for hole in np.arange(0,len(hole_w),1):

    BB_output_by_point = []
    total_energy_absorbed_by_cryoconite_by_point = []

    hole_depth = hole_d[hole]
    hole_width = hole_w[hole]
    hole_water_depth = hole_water_d[hole]
    hole_area = np.pi*((hole_width/200)**2) #in m

    ControlFuncs.Validate_Input_Data(hole_depth, hole_width, hole_water_depth, solzen)

    for point in np.arange(0, hole_width, 1):

        # function calls
        dir_energy_absorbed_by_cryoconite, diffuse_energy_absorbed_by_cryoconite, total_incoming_energy, \
            dir_energy_at_hole_floor, diffuse_energy_at_hole_floor = ControlFuncs.CalculateFluxes(\
            hole_depth, hole_width, hole_water_depth, point, cryoconite_albedo, WL, params, tolerance)

        total_energy_absorbed_by_cryoconite = diffuse_energy_absorbed_by_cryoconite + dir_energy_absorbed_by_cryoconite
        total_energy_absorbed_by_cryoconite_by_point.append(total_energy_absorbed_by_cryoconite)

        BB_output = np.sum(total_energy_absorbed_by_cryoconite)
        BB_output_by_point.append(BB_output)
    

    total_energy_absorbed_per_hole.append(np.sum(np.array(total_energy_absorbed_by_cryoconite)))


## BROADBAND ONLY ATM
total_energy_absorbed_per_hole = np.array(total_energy_absorbed_per_hole)*hole_areas
total_energy_absorbed_per_hole = total_energy_absorbed_per_hole*n_holes
total_absorbed_all_holes = np.sum(total_energy_absorbed_per_hole)


print(total_energy_absorbed_per_hole)
print(f"total area covered by cryoconite holes = {np.sum(hole_areas*n_holes)}")
print(f"For a {study_area} m2 area of ice with {np.sum(n_holes)} cryoconite holes:\n")
print(f"The total irradiance flux is {np.sum(incoming)*study_area}")
print(np.round(total_absorbed_all_holes,3)," Watts are absorbed by the cryoconite hole")
pc_abs = (np.sum(total_absorbed_all_holes)/(np.sum(incoming)*study_area))*100
print(f"% of incoming absorbed = {pc_abs}")