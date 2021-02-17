"""
Driver script for validation testing CryoconiteRTM

This validation is achieved by comparing predicted energy at the hole floor to field measurements
made by lowering pyranometers into cryoconite holes in the field.

AUTHOR: JOSEPH COOK, April 2020
www.tothepoles.co.uk
ww.github.com/jmcook1186

"""

import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
from ControlFuncs import ControlFuncs

########################
# IMPORT FIELD DATA
#######################

FieldDF = pd.read_csv('./TestData/FieldMeasurements.csv')
depths = FieldDF['HoleDepth(mm)']/10
widths = FieldDF['HoleWidth(mm)']/10
FieldRatios = FieldDF['Ratio']
modelRatio = np.zeros(len(depths))
sensor_height = 4 # sensor height in cm
incoming = np.genfromtxt ('/home/joe/Code/CryoconiteRTM/Data/mlw_sfc_flx_frc_clr.txt', delimiter=",") # path to I* file


for i in range(len(depths)):

    hole_d = depths[i] - sensor_height
    hole_w = widths[i]
    hole_water_d = np.round(hole_d * 0.5,0)

    #print("D: ", hole_d, "W: ", hole_w, "WatD: ", hole_water_d)

    hole_ar = int(hole_d/hole_w)
    point = hole_w/2 # horizontal distance from LH wall to desired location on hole floor
    cryoconite_albedo = np.ones(470)*0.2 #constant albedo across wavelength for now
    WL = np.arange(0.3,5,0.01)
    Mie = False
    GeometricOptics = True
    SZA = 20 # SZA is calculated as degrees from the vertical (zenith) - i.e. 0 is illumination from directly overhead

    #############################################
    ## 3. SET PHYSICAL PROPERTIES OF THE ICE/SNOW
    #############################################

    # These values are set to realistic values from the literature
    # for SW Greenland

    rho_snw = [600] # density of each layer (unit = kg m-3)
    # if using Mie optical properties, set spherical grain radius
    rds_snw = [2000]
    # if using GeometricOptics, set grain side_length and depth
    side_length = [4000] 
    depth = [4000]

    #############################################################
    # END OF USER INPUT (i.e. leave all remaining code unchanged)
    #############################################################

    dir_energy_absorbed_by_cryoconite, diffuse_energy_absorbed_by_cryoconite, total_incoming_energy,\
    dir_energy_at_hole_floor, diffuse_energy_at_hole_floor = ControlFuncs.CalculateFluxes(
        hole_d, hole_w, hole_water_d, point, cryoconite_albedo, WL, Mie, GeometricOptics, SZA, incoming,\
         rho_snw, rds_snw, side_length, depth)

    total_energy_at_hole_floor = diffuse_energy_at_hole_floor + dir_energy_at_hole_floor

    BB_output = np.sum(total_energy_at_hole_floor[:,0:40], axis=1)

    modelRatio[i] = BB_output/np.sum(incoming)


error = modelRatio - FieldRatios
norm_error = abs(error)

print("mean error = ", np.mean(norm_error))
print("STD error = ", np.std(norm_error))
plt.scatter(range(len(error)),norm_error)
plt.title("Absolute error: field measurements v simulations for hole floor irradiance")
plt.ylim(0,0.5)
plt.ylabel("Absolute error"), plt.xlabel("Measurement ID")
plt.show()