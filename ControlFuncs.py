"""
Class ControlFuncs contains functions that call direct and diffuse flux calculations and returns the
energy at the hole floor. This includes calls to SpecReflFuncs and SNICAR. Hard coded variables used
to configure the radiative transfer model and beam reflection losses are defined in here.

Functions in this class include:

1) CalculateFluxes
    Sets variable values and makes calls to external functions to calculate energy flux at cryoconite sediment layer

AUTHOR: JOSEPH COOK, April 2020
www.tothepoles.co.uk
ww.github.com/jmcook1186
"""



class ControlFuncs:

    def __init__(self):


        return

    def Validate_Input_Data(hole_d, hole_w, hole_water_d, solzen):

        if hole_water_d > hole_d:
            raise ValueError("ERROR: The water is deeper than the cryoconite hole")
        else:
            pass

        if (solzen ==0) or (solzen > 85):
            raise ValueError("ERROR: Please adjust solar zenith to be 1 - 85 degrees")
        else:
            pass

        if (hole_w == 0) or (hole_d ==0):
            raise ValueError("ERROR: Either the hole width or depth is set to zero")
        else:
            pass

        return 


    def CalculateFluxes(hole_d, hole_w, hole_water_d, point, cryoconite_albedo, WL, params, n_internal_reflections):

        import numpy as np
        import math
        from SpecReflFuncs import specFuncs
        from TwoStreamFuncs import TwoStreamFuncs

        #############################################
        # HARD CODED AND DERIVED VARIABLE DEFINITIONS
        #############################################

        dz = [hole_d/100] # thickness of each vertical layer (unit = m)
        R_sfc = np.mean(cryoconite_albedo) # reflectance of underlying surface - set across all wavelengths
        theta = 90-params.solzen # calculated from SZA
        nAir = np.ones(shape=(470))+0.0003 # define n and k for air (array of ones)
        kAir = np.zeros(shape=(470))+0.00000001

        # import spectral refractive index for water
        nWat = np.genfromtxt('/home/joe/Code/CryoconiteRTM/Data/water_n.csv', delimiter=",")
        kWat = np.genfromtxt ('/home/joe/Code/CryoconiteRTM/Data/water_k.csv', delimiter=",")
        kWat = kWat[0:-1:10] #every 10th element to match resolution of SNICAR

        nIce = np.genfromtxt('/home/joe/Code/CryoconiteRTM/Data/ice_n.csv', delimiter=",")
        nIce[nIce<1.0] = 1.0 # prevent math domain error - this is a negligible adjustment to a few wavelengths
        kIce = np.genfromtxt ('/home/joe/Code/CryoconiteRTM/Data/ice_k.csv', delimiter=",")

        incoming = TwoStreamFuncs.generate_incoming_irradiance(params)

        # set up empty lists
        R_airtowat = []
        R_wattoice = []
        dir_energy_at_hole_floor = []
        upbeam = []
        incoming_new = np.copy(incoming)

        ####################################
        # CALCULATE TRANSPORT OF DIRECT BEAM
        ####################################

        # 1) Illumination geometry

        # if there is no water, no refraction of incoming beam occurs so t_theta = theta
        if hole_water_d == 0:

            t_theta = theta
            t_theta_mean = theta

        else: # calculate adjusted solar elevation angle after direct beam refracted at air-water boundary
            t_theta = specFuncs.trans_angle(theta,nAir,nWat) 

        # call critical angle function to determine whether the direct beam reaches the hole floor at "point"
        ang_crit = specFuncs.critical_angle(theta, hole_d, hole_w, point) # calculate critical angle

        # 2) For each wavelength, calculate losses at medium boundaries and apply for n interactions

        for i in range(len(WL)):
            
            # calculate losses expected at each type of transition (air/water, water/ice)
            result1 = specFuncs.fresnel(nAir[i],nWat[i],kAir[i],kWat[i],theta)
            R_airtowat.append(result1)

            result2 = specFuncs.fresnel(nWat[i],nIce[i],kWat[i],kIce[i],t_theta[i])
            R_wattoice.append(result2)

            if t_theta[i] > ang_crit:

                # direct beam only hits point on hole floor when the refracted illumination angle 
                # exceeds the critical angle,
                
                if hole_water_d > 0: 
                    # if there is water, some energy is lost when beam enters from air
                    incoming_new[i] = incoming_new[i]*(1-R_airtowat[i])

                # the beam energy at the hole floor is just the total incoming after
                #  air/water fresnel loss
                dir_energy_at_hole_floor.append(incoming_new[i]) 

                # there are no reflections and the beam does not hit the hole wall
                beamHitsWall = False
                SurfStrike_d = hole_w/2
                n_wat_reflections = 0


            else:
                # if the beam does not directly illuminate the point on the floor, 
                # it may still reach the floor after multiple reflections betwen the 
                # hole walls. The following function calculates the number of reflections
                # between the hole walls before and after entering the water and prior 
                # to the beam striking the hole floor

                n_air_reflections, n_wat_reflections, total_reflections, SurfStrike_d,\
                beamHitsWall = specFuncs.test_multiple_reflections(
                    theta, t_theta[i], hole_d, hole_w, hole_water_d, nAir, nWat, verbose=False)
        
                # Use calculated # reflections to remove energy from the beam
                for j in range(int(n_air_reflections)+1):
                    # add one to account for the specular reflection occurring when beam 
                    # hits water surface (not in air_reflections or wat_reflections)
                    incoming_new[i] = incoming_new[i]*R_airtowat[i] 

                for j in range(int(n_wat_reflections)):
                    
                    incoming_new[i] = incoming_new[i]*R_wattoice[i]
                
                # this gives the amount of energy available at the hole floor
                dir_energy_at_hole_floor.append(incoming_new[i]) 
                
        
        # 3) calculate absorptive losses due to transport through water
        # First calculate path length in water, then calculate loss
        # using path length and absorption coefficient
        PathLengthInWat = specFuncs.CalculatePathLength(hole_water_d, hole_w, beamHitsWall, t_theta[i],\
        SurfStrike_d, n_wat_reflections, ang_crit)        
        dir_energy_at_hole_floor = specFuncs.AttenuateBeam(PathLengthInWat, kWat, dir_energy_at_hole_floor, WL)


        ####################################################
        ## CALCULATE DIFFUSE ENERGY FLUX REACHING HOLE FLOOR
        ####################################################

        albedo, BBA, F_btm_net = TwoStreamFuncs.call_snicar(params)
        albedo = albedo[10:]
        F_btm_net = F_btm_net[10:]
        
        ###############################################
        # CALCULATE ENERGY ABSORBED AT CRYOCONITE LAYER
        ###############################################

        dir_energy_absorbed_by_cryoconite = dir_energy_at_hole_floor * (1-cryoconite_albedo)
        diffuse_energy_at_hole_floor = F_btm_net

        diffuse_energy_absorbed_by_cryoconite = diffuse_energy_at_hole_floor * (1-cryoconite_albedo)

        total_incoming_energy = np.sum(incoming)

        ####################################
        ## ACCOUNT FOR INTERNAL REFLECTIONS
        ####################################

        energy_escaping_internal_reflections = []
        energy_lost_internal_reflections = []
        
        # loop through wavelengths
        for i in range(len(WL)): 
            escaped, loss, cryoconite_abs = specFuncs.internal_reflection(hole_water_d, cryoconite_albedo, WL[i], nAir[i], kAir[i], nWat[i], kWat[i], n_internal_reflections,\
                dir_energy_at_hole_floor[i], diffuse_energy_at_hole_floor[i])

            energy_escaping_internal_reflections.append(escaped)
            energy_lost_internal_reflections.append(loss)


        diffuse_energy_absorbed_by_cryoconite = diffuse_energy_absorbed_by_cryoconite + cryoconite_abs
        
        return dir_energy_absorbed_by_cryoconite, diffuse_energy_absorbed_by_cryoconite, total_incoming_energy,\
        dir_energy_at_hole_floor, diffuse_energy_at_hole_floor
