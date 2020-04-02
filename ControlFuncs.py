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

    def CalculateFluxes(hole_d, hole_w, hole_water_d, point, cryoconite_albedo, WL, Mie, GeometricOptics, SZA, incoming, 
    rho_snw, rds_snw, side_length, depth):

        import numpy as np
        import math
        from snicar8d_mie import snicar8d_mie
        from snicar8d_GO import snicar8d_GO
        from SpecReflFuncs import specFuncs

        #############################################
        # HARD CODED AND DERIVED VARIABLE DEFINITIONS
        #############################################


        dz = [hole_d/100] # thickness of each vertical layer (unit = m)
        nbr_lyr = len(dz)  # number of snow layers
        R_sfc = np.mean(cryoconite_albedo) # reflectance of undrlying surface - set across all wavelengths
        theta = 90-SZA # calculated from SZA
        nAir = np.ones(shape=(470)) # define n and k for air (array of ones)
        kAir = np.zeros(shape=(470))+0.00000001

        # import spectral refractive index for water
        nWat = np.genfromtxt('/home/joe/Code/CryoconiteRTM/Data/water_n.csv', delimiter=",")
        kWat = np.genfromtxt ('/home/joe/Code/CryoconiteRTM/Data/water_k.csv', delimiter=",")
        kWat = kWat[0:-1:10] #every 10th element to match resolution of SNICAR

        nIce = np.genfromtxt('/home/joe/Code/CryoconiteRTM/Data/ice_n.csv', delimiter=",")
        nIce[nIce<1.0] = 1.0 # prevent math domain error - this is a negligible adjustment to a few wavelengths
        kIce = np.genfromtxt ('/home/joe/Code/CryoconiteRTM/Data/ice_k.csv', delimiter=",")

        DIRECT   = 0        # 1= Direct-beam incident flux, 0= Diffuse incident flux
        APRX_TYP = 1        # 1= Eddington, 2= Quadrature, 3= Hemispheric Mean
        DELTA    = 1        # 1= Apply Delta approximation, 0= No delta
        coszen   = math.cos(SZA*(np.pi/180))    # if DIRECT give cosine of solar zenith angle 

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

        nbr_aer = 16 # Define total number of different LAPs/aerosols in model

        # set filename stubs
        stb1 = 'algae_geom_' # %name stub 1
        stb2 = '.nc'  # file extension
        wrkdir2 = '/home/joe/Code/BioSNICAR_GO_PY/Data/Algal_Optical_Props/' # working directory
        snw_stb1 = 'snw_alg_' # name stub for snow algae

        # CHOOSE DIMENSIONS OF GLACIER ALGAE 1
        algae_r = 1 # algae radius
        algae_l = 1 # algae length
        glacier_algae1 = str(wrkdir2+stb1+str(algae_r)+'_'+str(algae_l)+stb2) # create filename string

        # CHOOSE DIMENSIONS OF GLACIER ALGAE 2
        algae2_r = 1 # algae radius
        algae2_l = 1 # algae length
        glacier_algae2 = str(wrkdir2+stb1+str(algae2_r)+'_'+str(algae2_l)+stb2) # create filename string

        # CHOOSE SNOW ALGAE DIAMETER
        snw_algae_r = 1 # snow algae diameter
        snw_alg = str(wrkdir2+snw_stb1+str(snw_algae_r)+stb2) # create filename string


        mss_cnc_soot1 = [0]    # uncoated black carbon
        mss_cnc_soot2 = [0]    # coated black carbon
        mss_cnc_dust1 = [0]    # global average dust 1
        mss_cnc_dust2 = [0]    # global average dust 2
        mss_cnc_dust3 = [0]    # global average dust 3
        mss_cnc_dust4 = [0]    # global average dust 4
        mss_cnc_ash1 = [0]    # volcanic ash species 1
        mss_cnc_GRISdust1 = [0]    # GRIS dust 1 (Cook et al. 2019 "mean")
        mss_cnc_GRISdust2 = [0]    # GRIS dust 2 (Cook et al. 2019 HIGH)
        mss_cnc_GRISdust3 = [0]    # GRIS dust 3 (Cook et al. 2019 LOW)
        mss_cnc_GRISdustP1 = [0]  # GRIS dust 1 (Polashenki2015: low hematite)
        mss_cnc_GRISdustP2 = [0]  # GRIS dust 1 (Polashenki2015: median hematite)
        mss_cnc_GRISdustP3 = [0]  # GRIS dust 1 (Polashenki2015: median hematite)
        mss_cnc_snw_alg = [0]    # Snow Algae (spherical, C nivalis)
        mss_cnc_glacier_algae1 = [0]    # glacier algae type1
        mss_cnc_glacier_algae2 = [0]    # glacier algae type2


        # SET FILE NAMES CONTAINING OPTICAL PARAMETERS FOR ALL IMPURITIES:

        FILE_soot1  = 'mie_sot_ChC90_dns_1317.nc'
        FILE_soot2  = 'miecot_slfsot_ChC90_dns_1317.nc'
        FILE_dust1  = 'aer_dst_bln_20060904_01.nc'
        FILE_dust2  = 'aer_dst_bln_20060904_02.nc'
        FILE_dust3  = 'aer_dst_bln_20060904_03.nc'
        FILE_dust4  = 'aer_dst_bln_20060904_04.nc'
        FILE_ash1  = 'volc_ash_mtsthelens_20081011.nc'
        FILE_GRISdust1 = 'dust_greenland_Cook_CENTRAL_20190911.nc'
        FILE_GRISdust2 = 'dust_greenland_Cook_HIGH_20190911.nc'
        FILE_GRISdust3 = 'dust_greenland_Cook_LOW_20190911.nc'
        FILE_GRISdustP1 = 'dust_greenland_L_20150308.nc'
        FILE_GRISdustP2 = 'dust_greenland_C_20150308.nc'
        FILE_GRISdustP3 = 'dust_greenland_H_20150308.nc'
        FILE_snw_alg  = snw_alg # snow algae (c nivalis)
        FILE_glacier_algae1 = glacier_algae1 # Glacier algae
        FILE_glacier_algae2 = glacier_algae2 # Glacier algae


        # call snicar functions

        if Mie == True and GeometricOptics == True:
            print("ERROR: BOTH MIE AND GO SELECTED: PLEASE CHOOSE ONE")

        elif Mie == True:
            [wvl, albedo, BBA, BBAVIS, BBANIR, abs_slr, abs_slr_tot, abs_vis_tot, heat_rt, total_insolation, F_btm_net] = snicar8d_mie(DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz, rho_snw, rds_snw, nbr_lyr, nbr_aer, mss_cnc_soot1,
            mss_cnc_soot2, mss_cnc_dust1, mss_cnc_dust2, mss_cnc_dust3, mss_cnc_dust4, mss_cnc_ash1, mss_cnc_GRISdust1, 
            mss_cnc_GRISdust2, mss_cnc_GRISdust3, mss_cnc_GRISdustP1, mss_cnc_GRISdustP2, mss_cnc_GRISdustP3, 
            mss_cnc_snw_alg, mss_cnc_glacier_algae1, mss_cnc_glacier_algae2, FILE_soot1, FILE_soot2, FILE_dust1, FILE_dust2,
            FILE_dust3, FILE_dust4, FILE_ash1, FILE_GRISdust1, FILE_GRISdust2, FILE_GRISdust3, FILE_GRISdustP1, FILE_GRISdustP2, 
            FILE_GRISdustP3, FILE_snw_alg, FILE_glacier_algae1, FILE_glacier_algae2)

        elif GeometricOptics == True:

            [wvl, albedo, BBA, BBAVIS, BBANIR, abs_slr, abs_slr_tot, abs_vis_tot, heat_rt, total_insolation, F_btm_net] = snicar8d_GO(DIRECT, APRX_TYP, DELTA, coszen, R_sfc, dz, rho_snw, side_length, depth, nbr_lyr, nbr_aer, mss_cnc_soot1,
            mss_cnc_soot2, mss_cnc_dust1, mss_cnc_dust2, mss_cnc_dust3, mss_cnc_dust4, mss_cnc_ash1, mss_cnc_GRISdust1, 
            mss_cnc_GRISdust2, mss_cnc_GRISdust3, mss_cnc_GRISdustP1, mss_cnc_GRISdustP2, mss_cnc_GRISdustP3, 
            mss_cnc_snw_alg, mss_cnc_glacier_algae1, mss_cnc_glacier_algae2, FILE_soot1, FILE_soot2, FILE_dust1, FILE_dust2,
            FILE_dust3, FILE_dust4, FILE_ash1, FILE_GRISdust1, FILE_GRISdust2, FILE_GRISdust3, FILE_GRISdustP1, FILE_GRISdustP2, 
            FILE_GRISdustP3, FILE_snw_alg, FILE_glacier_algae1, FILE_glacier_algae2)

        else:
            print("NEITHER MIE NOR GO SELECTED: PLEASE CHOOSE ONE")

        
        ########################################
        # CALCULATE AND DISPLAY OUTPUT VARIABLES
        ########################################

        dir_energy_absorbed_by_cryoconite = dir_energy_at_hole_floor * (1-cryoconite_albedo)

        diffuse_energy_absorbed_by_cryoconite = F_btm_net * (1-cryoconite_albedo)

        total_incoming_energy = np.sum(incoming)
        
        return dir_energy_absorbed_by_cryoconite, diffuse_energy_absorbed_by_cryoconite, total_incoming_energy