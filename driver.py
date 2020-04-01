import numpy as np
import math
import matplotlib.pyplot as plt
from SpecReflFuncs import specFuncs


# TODO: Find and debug edge cases (e.g. small water depths can lead to -ve #reflections)
# TODO: consider the use of floor_strike_d to examine energy distribution across floor 


######################
# Define Hole geometry
######################

hole_d = 30
hole_w = 50
hole_water_d = 20
hole_ar = hole_d/hole_w
point = 40 # horizontal distance from LH wall to desired location on hole floor
cryoconite_albedo = np.ones(470)*0.2 #constant albedo across wavelength for now
WL = np.arange(0.3,5,0.01)

# define incoming irradiance
SZA = 30 # SZA is calculated as degrees from the vertical (zenith) - i.e. 0 is illumination from directly overhead
theta = 90-SZA # calculated from SZA
# define incoming irradiance in Wm-2
incoming = np.genfromtxt ('/home/joe/Code/CryoconiteRTM/Data/mlw_sfc_flx_frc_clr.txt', delimiter=",")

# define n and k for air (array of ones)
nAir = np.ones(shape=(470))
kAir = np.zeros(shape=(470))+0.00000001

# import spectral refractive index for water
nWat = np.genfromtxt('/home/joe/Code/CryoconiteRTM/Data/water_n.csv', delimiter=",")
kWat = np.genfromtxt ('/home/joe/Code/CryoconiteRTM/Data/water_k.csv', delimiter=",")
kWat = kWat[0:-1:10] #every 10th element to match resolution of SNICAR

nIce = np.genfromtxt('/home/joe/Code/CryoconiteRTM/Data/ice_n.csv', delimiter=",")
nIce[nIce<1.0] = 1.0 # prevent math domain error - this is a negligible adjustment to a few wavelengths
kIce = np.genfromtxt ('/home/joe/Code/CryoconiteRTM/Data/ice_k.csv', delimiter=",")

# call critical angle function to determine whether the direct beam reaches the hole floor at "point"
ang_crit, hyp = specFuncs.critical_angle(theta, hole_d, hole_w, point) # calculate critical angle and path lengh (hyp)

# if there is no water, no refraction of incoming beam occurs so t_theta = theta
if hole_water_d == 0:
    t_theta = theta
    t_theta_mean = theta
else:
    t_theta = specFuncs.trans_angle(theta,nAir,nWat) # calculate adjusted solar elevation angle after direct beam refracted at air-water boundary
    t_theta_mean = np.mean(t_theta) # calculate a mean across wavelengths

print("critical angle = ", np.round(ang_crit,2))
print("mean transmitted angle = ", np.round(t_theta_mean,2))


# set up empty lists
R_airtowat = []
R_wattoice = []
energy_at_hole_floor = []
upbeam = []
incoming_new = np.copy(incoming)
tot_refl = []
wat_refl = []
air_refl = []


for i in range(len(WL)):
    
    # calculate losses expected at each type of transition (air/water, water/ice)
    result1 = specFuncs.fresnel(nAir[i],nWat[i],kAir[i],kWat[i],theta)
    R_airtowat.append(result1)

    result2 = specFuncs.fresnel(nWat[i],nIce[i],kWat[i],kIce[i],t_theta[i])
    R_wattoice.append(result2)

    # direct beam only hits point on hole floor when the refracted illumination angle exceeds the critical angle
    if t_theta[i] > ang_crit:
        
        print("DIRECT BEAM HITS HOLE FLOOR")

        incoming_new[i] = incoming_new[i]*(1-R_airtowat[i])

        energy_at_hole_floor.append(incoming_new[i]) # subract energy lost in reflection from water surface

        air_refl.append(0)
        wat_refl.append(0)
        tot_refl.append(0)

    else:
        # the following function calculates the number of reflections between the hole walls before and after entering the water
        # and prior to the beam striking the hole floor

        n_air_reflections, n_wat_reflections, total_reflections, floor_strike_d = specFuncs.test_multiple_reflections(theta, t_theta[i], hole_d, hole_w, hole_water_d, nAir, nWat, verbose=False)
 
        tot_refl.append(total_reflections)
        wat_refl.append(n_wat_reflections)
        air_refl.append(n_air_reflections)

        # REDUCE POWER IN INCOMING BEAM AT EACH REFLECTION
        for j in range(int(n_air_reflections)+1):
            # add one to account for the specular reflection occurring when beam hits water surface (not in air_reflections or wat_reflections)
            incoming_new[i] = incoming_new[i]*R_airtowat[i] # only the reflected energy remains available for further reflections and possible
            # interaction with the cryoconite layer

        for j in range(int(n_wat_reflections)):
            
            incoming_new[i] = incoming_new[i]*R_wattoice[i]

        energy_at_hole_floor.append(incoming_new[i])



total_incoming_energy = np.sum(incoming)

energy_at_hole_floor = np.array(energy_at_hole_floor)
total_energy_at_hole_floor = np.sum(energy_at_hole_floor)

energy_absorbed_by_cryoconite = energy_at_hole_floor * (1-cryoconite_albedo)
total_energy_absorbed_by_cryoconite = np.sum(energy_absorbed_by_cryoconite)
print("total incoming energy = ",np.round(total_incoming_energy,5))

print("total_energy_absorbed_by_cryoconite = ", np.round(total_energy_absorbed_by_cryoconite,2))
print("total energy at hole floor = ", np.round(total_energy_at_hole_floor,2))

print("proportion of incoming energy reaching hole floor = ", np.round(total_energy_at_hole_floor/total_incoming_energy,5))
print("proportion of incoming energy absorbed = ", np.round(total_energy_absorbed_by_cryoconite/total_incoming_energy,5))

plt.plot(WL,energy_at_hole_floor.T,label='energy at hole floor (w/m2)')
plt.plot(WL,energy_absorbed_by_cryoconite.T,label='energy absorbed by CC (w/m2)')
plt.plot(WL,incoming,label='total incoming energy (w/m2)')
plt.xlim(0.3,5),plt.xlabel('Wavelength (microns)')
plt.ylabel('Energy (W/m2)')
plt.legend(loc='best')
plt.show()

plt.figure()
plt.xlim(0.3,5),plt.xlabel('Wavelength (micron')
plt.ylabel('Number of reflections')
plt.plot(WL,tot_refl,label='total reflections')
plt.plot(WL,air_refl, label='reflections above water surface')
plt.plot(WL,wat_refl, label='reflections below water surface')
plt.legend(loc='best')
plt.show()
