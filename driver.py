import numpy as np
import math
import matplotlib.pyplot as plt
from SpecReflFuncs import specFuncs



# TODO: work out why SZA=45 breaks 
# TODO: make spectral!
# TODO: test hole geometries
# TODO: consider the use of floor_strike_d to examine energy distribution across floor 


######################
# Define Hole geometry
######################

hole_d = 50
hole_w = 20
hole_water_d = 30
hole_ar = hole_d/hole_w
point = hole_w/2 # horizontal distance from LH wall to desired location on hole floor
cryoconite_albedo = np.ones(470)*0.2 #constant albedo across wavelength for now
WL = np.arange(0.3,5,0.01)

# define incoming irradiance
SZA = 70
theta = 90-SZA # calculated from SZA
# define incoming irradiance in Wm-2
incoming = np.genfromtxt ('/home/joe/Code/CryoconiteRTM/Data/mlw_sfc_flx_frc_clr.txt', delimiter=",")

# define n and k for air (array of ones)
nAir = np.ones(shape=(470))
kAir = np.zeros(shape=(470))

# import spectral refractive index for water
nWat = np.genfromtxt('/home/joe/Code/CryoconiteRTM/Data/water_n.csv', delimiter=",")
kWat = np.genfromtxt ('/home/joe/Code/CryoconiteRTM/Data/water_k.csv', delimiter=",")
kWat = kWat[0:-1:10] #every 10th element to match resolution of SNICAR

nIce = np.genfromtxt('/home/joe/Code/CryoconiteRTM/Data/ice_n.csv', delimiter=",")
kIce = np.genfromtxt ('/home/joe/Code/CryoconiteRTM/Data/ice_k.csv', delimiter=",")


# call critical angle function to determine whether the direct beam reaches the hole floor at "point"
ang_crit, hyp = specFuncs.critical_angle(theta, hole_d, hole_w, point) # calculate critical angle and path lengh (hyp)

t_theta = specFuncs.trans_angle(SZA,nAir,nWat) # calculate adjusted solar elevation angle after direct beam refracted at air-water boundary
t_theta_mean = np.mean(t_theta)

print("critical angle = ", np.round(ang_crit,2))
print("mean transmitted angle = ", np.round(t_theta_mean,2))

# calculate losses expected at each type of transition (air/water, water/ice)
R_airtowat = []
R_wattoice = []
energy_at_hole_floor = []
upbeam = []
incoming_new = np.copy(incoming)


for i in range(len(nAir)):
    
    result1 = specFuncs.fresnel(nAir[i],nWat[i],kAir[i],kWat[i],theta)
    R_airtowat.append(result1)
    result2 = specFuncs.fresnel(nWat[i],nIce[i],kWat[i],kIce[i],t_theta[i])
    R_wattoice.append(result2)

    # direct beam only hits point on hole floor when the refracted illumination angle exceeds the critical angle
    if t_theta[i] > ang_crit:

        energy_at_hole_floor.append((incoming[i] * (1-R_airtowat[i])))

    else:
        # the following function calculates the number of reflections between the hole walls before and after entering the water
        # and prior to the beam striking the hole floor

        n_air_reflections, n_wat_reflections, total_reflections, floor_strike_d = specFuncs.test_multiple_reflections(theta, t_theta[i], hole_d, hole_w, hole_water_d)

        # REDUCE POWER IN INCOMING BEAM AT EACH REFLECTION

        for j in range(int(n_air_reflections)):

            incoming_new[i] = incoming_new[i]*R_airtowat[i] # only the reflected energy remains available for further reflections and possible
            # interaction with the cryoconite layer

        for j in range(int(n_wat_reflections)):

            incoming_new[i] = incoming_new[i]*R_wattoice[i]

        energy_at_hole_floor.append(incoming_new)



energy_absorbed_by_cryoconite = energy_at_hole_floor * cryoconite_albedo
upwelling = energy_at_hole_floor - energy_absorbed_by_cryoconite
total_absorbed = np.sum(energy_absorbed_by_cryoconite)
total_upwelling = np.sum(upwelling)

print("Energy received by cryoconite = {}".format(np.round(total_absorbed,3)))
print("Energy returned to upwelling field = {}".format(np.round(total_upwelling,3)))

print("Total downwelling energy",np.sum(incoming))
print("Total energy absorbed by cryoconite",total_absorbed)

print("Proportion of incoming absorbed by cryoconite = {}%".format(np.round((total_absorbed/np.sum(incoming))*100,2)))