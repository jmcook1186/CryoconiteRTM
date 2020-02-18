import numpy as np
import math
import matplotlib.pyplot as plt
from SpecReflFuncs import specFuncs


###############
# Define Hole geometry
###############
incoming = 1.0

hole_d = 20
hole_w = 10
hole_water_d = 10
hole_ar = hole_d/hole_w
point = hole_w/2 # horizontal distance from LH wall to desired location on hole floor
cryoconite_albedo = 0.2

nWat = 1.3333333
nIce = 1.31
nAir = 1.0

kWat = 0.00001
kIce = 0.00002
kAir = 0.0000001

energy =  []

SZA = 45

theta = 90-SZA

# call critical angle function to determine whether the direct beam reaches the hole floor at "point"
ang_crit, hyp = specFuncs.critical_angle(theta, hole_d, hole_w, point) # calculate critical angle and path lengh (hyp)
t_theta = specFuncs.trans_angle(SZA,nAir,nWat) # calculate adjusted solar elevation angle after direct beam refracted at air-water boundary

print("critical angle = ", ang_crit)
print("t_theta = ", t_theta)

# direct beam only hits point on hole floor when the refracted illumination angle exceeds the critical angle

if t_theta > ang_crit:

    print("\nDIRECT ILLUMINATION: PROCEEDING WITH FRESNEL CALCULATION AT AIR/WATER BOUNDARY\n")

    R_airtowat = specFuncs.fresnel(nAir,nWat,kAir,kWat,theta)

    print("DIRECT BEAM REFLECTED = ", incoming*(R_airtowat))
    downbeam = incoming * (1-R_airtowat)
    upbeam = downbeam*cryoconite_albedo

    print("Energy received by cryoconite = {}".format(downbeam))
    print("Energy returned to upwelling field = {}".format(upbeam))


else:
    # the following function calculates the number of reflections between the hole walls before and after entering the water
    # and prior to the beam striking the hole floor
    n_air_reflections, n_wat_reflections, total_reflections = specFuncs.test_multiple_reflections(theta, t_theta, hole_d, hole_w, hole_water_d)

    # calculate losses at each of air/water and water/ice
    R_airtowat = specFuncs.fresnel(nAir,nWat,kAir,kWat,theta)
    R_wattoice = specFuncs.fresnel(nWat,nIce,kWat,kIce,t_theta)

    print("R air/wat = ",R_airtowat)
    print("R wat/ice = ", R_wattoice)

    for i in range(int(n_air_reflections)):
        incoming = incoming*R_airtowat # only the reflected energy remains available for further reflections and possible
        # interaction with the cryoconite layer

    for i in range(int(n_wat_reflections)):
        incoming = incoming*R_wattoice

    downbeam = incoming

print("beam power at cryoconite layer: ",downbeam)

