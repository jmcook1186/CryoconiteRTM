

import numpy as np
import math
import matplotlib.pyplot as plt
from ControlFuncs import ControlFuncs
from TwoStreamFuncs import TwoStreamFuncs

########################
# 1 DEFINE HOLE GEOMETRY
########################

hole_d = [0.2, 0.2] # depth of each hole (m)
hole_w = [0.5, 0.05] # width of each hole (m)
hole_water_d = [0.10, 0.10] # water depth in each hole (m)
n_holes = [3, 2] # how many of each hole (indexes match with dims)?
hole_areas= np.pi*((np.array(hole_w)/2)**2) # in m^2
total_cryoconite_area = np.sum(n_holes*hole_areas)
study_area = 1 # total study area in m^2

#constant albedo across wavelength for now
# can be udpated with measured spectrum later
cryoconite_albedo = np.ones(470)*0.2 
WL = np.arange(0.3,5,0.01)

####################
## 2. CONFIGURE RTM 
####################

solzen = 15
density = [500]
grain_rds = [500]
layer_type = [1]
dz = [0.2] # cm to m
algae = 0
incoming_i = 4
DIRECT = True
tolerance = 1e-10 #how close to zero doe the flux need to get before we stop iterating internal reflections?

# create named tuple containing snicar input params
params = TwoStreamFuncs.generate_ice_physical_params(density,grain_rds,layer_type,dz,algae,solzen,incoming_i,DIRECT)
incoming = TwoStreamFuncs.generate_incoming_irradiance(params)

BB_out_all_holes = np.zeros(shape=(len(hole_w),len(WL)))
absolute_energy_absorbed_all_holes = np.zeros(shape=(len(hole_w),1))


#############################################################
# END OF USER INPUT (i.e. leave all remaining code unchanged)
#############################################################

# Validation function will raise errors if input data is invalid
total_energy_absorbed_per_hole = []
spectral_energy_absorbed_per_hole= np.zeros(shape=(len(hole_d), len(WL)))
spectral_energy_escaping_per_hole= np.zeros(shape=(len(hole_d), len(WL)))
energy_escaping_internal_reflections_per_hole = np.zeros(shape=(len(hole_d), len(WL)))


for hole in np.arange(0,len(hole_w),1):

    BB_output_by_point = []
    total_energy_absorbed_by_cryoconite_by_point = []

    hole_depth = hole_d[hole]
    hole_width = hole_w[hole]
    hole_water_depth = hole_water_d[hole]
    hole_area = np.pi*((hole_width/200)**2) #in m

    ControlFuncs.Validate_Input_Data(hole_depth, hole_width, hole_water_depth, solzen)

    energy_escaping_internal_reflections_per_point=np.zeros(shape=(int(hole_width*100),len(WL)))

    for point in np.arange(0, int(hole_width*100), 1):

        # function calls
        energy_escaping_internal_reflections,dir_energy_absorbed_by_cryoconite, diffuse_energy_absorbed_by_cryoconite, total_incoming_energy, \
            dir_energy_at_hole_floor, diffuse_energy_at_hole_floor, F_top_pls, reflected_from_water_surface = ControlFuncs.CalculateFluxes(\
            hole_depth, hole_width, hole_water_depth, point, cryoconite_albedo, WL, params, tolerance)

        total_energy_absorbed_by_cryoconite = diffuse_energy_absorbed_by_cryoconite + dir_energy_absorbed_by_cryoconite
        total_energy_absorbed_by_cryoconite_by_point.append(total_energy_absorbed_by_cryoconite)

        BB_output = np.sum(total_energy_absorbed_by_cryoconite)
        BB_output_by_point.append(BB_output)
        energy_escaping_internal_reflections_per_point[point,:] = energy_escaping_internal_reflections
    
    energy_escaping_internal_reflections_per_hole[hole,:] = energy_escaping_internal_reflections_per_point.mean(axis=0)
    spectral_energy_absorbed_per_hole[hole,:] = np.array(total_energy_absorbed_by_cryoconite_by_point).sum(axis=0)
    total_energy_absorbed_per_hole.append(np.sum(np.array(total_energy_absorbed_by_cryoconite)))

print(energy_escaping_internal_reflections_per_hole.shape)

#spectral_energy in all holes
for i in range(len(n_holes)):

    energy_escaping_internal_reflections_all_holes = energy_escaping_internal_reflections_per_hole[i,:]*n_holes[i]
print(energy_escaping_internal_reflections_all_holes.shape)
    



up1 = np.array(F_top_pls[10:])*(study_area-total_cryoconite_area)
up2 = np.array(reflected_from_water_surface)*total_cryoconite_area
up3 = np.array(energy_escaping_internal_reflections_all_holes)*total_cryoconite_area

print(up1.shape,up2.shape,up3.shape)

new_up = up1+up2+up3

albedo = new_up/incoming
albedo[albedo<0]=0.0001

plt.plot(new_up, color='b', marker='x', label='up total')
plt.plot(incoming,color='r', label='incoming')
plt.plot(F_top_pls[10:]*(study_area-total_cryoconite_area),color='g', label='up from non-conite areas')
plt.plot(reflected_from_water_surface*total_cryoconite_area, color='k', label='reflected from water surf')
plt.plot(energy_escaping_internal_reflections_all_holes*total_cryoconite_area, color='k',linestyle='--', label='escaping internal refl')

plt.legend(loc='best')
# plt.plot(new_up),plt.plot(incoming)
# plt.plot(F_top_pls[10:]),plt.plot(spectral_energy_absorbed_all_holes), plt.plot(reflected_from_water_surface,color='r')
# plt.plot(incoming, color='k'),plt.plot(spectral_energy_absorbed_all_holes,color='b'),plt.plot(F_top_pls[10:],color='g')

plt.savefig('test.jpg')

plt.figure()
plt.plot(albedo),plt.ylim(0,1),plt.savefig('albedo.jpg')

total_energy_absorbed_per_hole = np.array(total_energy_absorbed_per_hole)*hole_areas
total_energy_absorbed_per_hole = total_energy_absorbed_per_hole*n_holes
total_absorbed_all_holes = np.sum(total_energy_absorbed_per_hole)


print(total_energy_absorbed_per_hole)
print(f"total area covered by cryoconite holes = {np.sum(hole_areas*n_holes)}")
print(f"For a {study_area} m2 area of ice with {np.sum(n_holes)} cryoconite holes:\n")
print(f"The total irradiance flux is {np.sum(incoming)*study_area}")
print(np.round(total_absorbed_all_holes,3)," Watts are absorbed by the cryoconite holes")
pc_abs = (np.sum(total_absorbed_all_holes)/(np.sum(incoming)*study_area))*100
print(f"% of incoming absorbed = {pc_abs}")