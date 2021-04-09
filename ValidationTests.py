

import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
from ControlFuncs import ControlFuncs
from TwoStreamFuncs import TwoStreamFuncs
########################
# IMPORT FIELD DATA
#######################

FieldDF = pd.read_csv('./TestData/FieldMeasurements.csv')
depths = FieldDF['HoleDepth(mm)']/10
widths = FieldDF['HoleWidth(mm)']/10
FieldRatios = FieldDF['Ratio']
modelRatio = np.zeros(len(depths))
sensor_height = 5 # sensor height in cm

for i in range(len(depths)):

    hole_d = depths[i] - sensor_height
    hole_w = widths[i]
    hole_water_d = np.round(hole_d * 0.7,0)

    #print("D: ", hole_d, "W: ", hole_w, "WatD: ", hole_water_d)

    hole_ar = int(hole_d/hole_w)
    point = hole_w/2 # horizontal distance from LH wall to desired location on hole floor
    cryoconite_albedo = np.ones(470)*0.2 #constant albedo across wavelength for now
    WL = np.arange(0.3,5,0.01)

    #############################################
    ## 3. SET PHYSICAL PROPERTIES OF THE ICE/SNOW
    #############################################

    solzen = 25
    density = [850]
    grain_rds = [850]
    layer_type = [1]
    dz = [hole_d/100] # cm to m
    algae = 0
    incoming_i = 4
    DIRECT = True
    n_internal_reflections = 100

    # create named tuple containing snicar input params
    params = TwoStreamFuncs.generate_ice_physical_params(density,grain_rds,layer_type,dz,algae,solzen,incoming_i,DIRECT)
    incoming = TwoStreamFuncs.generate_incoming_irradiance(params)

    #############################################################
    # END OF USER INPUT (i.e. leave all remaining code unchanged)
    #############################################################

    dir_energy_absorbed_by_cryoconite, diffuse_energy_absorbed_by_cryoconite, total_incoming_energy,\
    dir_energy_at_hole_floor, diffuse_energy_at_hole_floor = ControlFuncs.CalculateFluxes(\
        hole_d, hole_w, hole_water_d, point, cryoconite_albedo, WL, params, n_internal_reflections)

    total_energy_at_hole_floor = diffuse_energy_at_hole_floor + dir_energy_at_hole_floor

    BB_output = np.sum(total_energy_at_hole_floor[0:40])

    modelRatio[i] = BB_output/np.sum(incoming)


error = modelRatio - FieldRatios
norm_error = abs(error)

print("mean error = ", np.mean(norm_error))
print("STD error = ", np.std(norm_error))
plt.scatter(range(len(error)),norm_error)
plt.title("Absolute error: field measurements v simulations for hole floor irradiance")
plt.ylim(0,0.5)
plt.ylabel("Absolute error"), plt.xlabel("Measurement ID")
plt.savefig('/home/joe/Code/CryoconiteRTM/Assets/ValidationTests.jpg')