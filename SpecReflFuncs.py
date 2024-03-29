"""

Class specFuncs contains functions for specular reflection, fresnel and snells law
functions that enable calculation of losses and transfers from the direct beam.

Functions in this class include:

1) critical angle
    Calculates minimum transmitted elevation angle required for given point on hole floor to be directly illuminated

2)trans_angle
    Calculates angle of beam after it has been refracted at the air/water boundary

3) test_multiple_reflections
    Calculates the total nmber of reflections occurring above and below the water surface for beams that hit the hole walls
    Contains embedded methods:
    a) DoesBeamHitWall
        Checks whether the beam strikes the hole wall before the water surface
    b) ReflectionsInAir
        Calculates the number of times the beam reflects between the hole walls before hitting the water surface
    c) UpdateSurfStrike
        Calculates the distance from the sunward wall that the beam strikes the water surface, and also sets the boolean
        "REVERSE" to true if the beam has reflected an odd number of times and is therefore incident from the sun-facing side.
    d) ReflectionsInWater
        Calculates how many times the beam reflects between the walls under the water before hitting the hole floor

4) fresnel
    Calculates the magnitude of energy losses expected at each material boundary (air/water, water/ice)

AUTHOR: JOSEPH COOK, April 2020
www.tothepoles.co.uk
ww.github.com/jmcook1186

"""

class specFuncs:

    def __init__(self):


        return


    def critical_angle(theta, hole_d, hole_w, point):

        """
        function takes hole dimensions and calculates angle between point on hole floor and 
        edge of hole aperture. This is the critical angle that must be exceeded for a direct beam 
        to reach that point on the hole floor. Code compares the critical angle to the
        solar elevation angle (90 - SZA) and returns a Boolean where 1 = direct illumination,
        0 = no direct illumination.
        
        """
        
        import math
        import numpy as np
        
        # calculate angle between "point" on hole floor and hole aperture
        
        d = hole_w - point # horizontal distance along hole floor from point to hole wall

        # use pythagorus' theorem to find length of hypotenuse, then inverse sin (opp/hyp) 
        # to find missing angle, ang_x

        hyp = math.sqrt(d**2+hole_d**2)
        ang_crit = math.asin(hole_d/hyp)
        ang_crit = ang_crit * 180/(np.pi) # convert rads to degrees
       
        return ang_crit


    def trans_angle(theta, nAir, nWat):
        
        """
        function calculates the adjusted angle of the direct beam accounting for refraction at the
        air/water interface.

        This is necessary because refraction changes the value of theta required to directly
        illuminate the hole floor

        """
        import math
        import numpy as np

        # calculate adjusted angle after incoming direct beam has entered water
        # Snells Law: sin(transmitted_angle) = (nAir*sin(incident angle)) /  nWat
        # the relevant angle to use for the incoming beam is the SZA not the elevation angle because
        # Snell's law expects angles measured form vertical normal
        # the predicted angle (out_ang) is the angle between the transmitted beam and the vertical normal
        # beneath the water surface. The transmitted beam elevation angle (t_theta) is 
        # equal to 90 - out_ang.

        theta_rad = theta * (np.pi/180)
                
        t_theta = np.zeros(len(nAir))

        for i in range(len(nAir)):

            t_theta_temp = math.asin((nAir[i] * math.sin(theta_rad))/nWat[i])

            t_theta[i] = t_theta_temp * 180/np.pi

        return  t_theta



    def test_multiple_reflections(theta, t_theta, hole_d, hole_w, hole_water_d, nAir, nWat, verbose = False):
        
        """
        function tests whether a single reflection from the hole wall causes the beam to hit the hole floor.
        If not, how many reflections between the walls occur before the beam reaches the hole floor?
        
        """

        import numpy as np
        import math

        SZA = 90-theta

        def DoesBeamHitWall(theta, hole_d, hole_water_d, hole_w):
            """
            Tests whether the beam hits the wall ABOVE the water surface
            """
            ang = theta
            ang_rad = ang*(np.pi/180) # radians

            SurfToWat = hole_d - hole_water_d
            SurfStrike_d = SurfToWat/(math.tan(ang_rad))

            if SurfStrike_d >= hole_w:
                
                BeamHitsWall = True

            else:
                BeamHitsWall = False

                            
            return BeamHitsWall, SurfStrike_d, SurfToWat


        def ReflectionsInAir(hole_w, theta, SurfToWat, BeamHitsWall = None):
            
            """
            calculates number of times beam reflects from walls before hitting water surface and provides 
            height above surface of final beam strike on wall before beam hits water

            """

            if BeamHitsWall:

                ang = theta*(np.pi/180)
                beam_d_air = hole_w*(math.tan(ang))

                n_air_reflections = np.floor(SurfToWat/beam_d_air)
                residual_d = SurfToWat - (n_air_reflections*beam_d_air)
            
            else:

                beam_d_air = SurfToWat
                n_air_reflections = 0
                residual_d = 0
            
            return n_air_reflections, residual_d, beam_d_air


        def UpdateSurfStrike_d(n_air_reflections, theta, SurfStrike_d, residual_d):
            
            """
            if there are reflections above the water surface, the surface strike position needs to be updated. In this initial
            calculation, the beam overshooting the hole boundaries is the diagnostic tool used to determine whether the beam 
            hits the wall rather than the water. Therefore, when that occurs, at least one reflection must be taken into account
            and used to reposition the beam on the water surface after n reflections.

            REVERSE is a boolean that flags whether the beam has changed direction due to an uneven number of reflections, if
            FALSE the beam is still arriving at the hole floor from the sunward side, if False the beam has reflected and is 
            arriving from the non-sunwards side (this changes the floor geometry calculations) 


            NB SurfStrike_D is always expressed as distance from sunwards wall
            """            

            if n_air_reflections == 0:

                SurfStrike_d = SurfStrike_d

                REVERSE = False
            
            elif n_air_reflections %2 != 0:

                ang_top = 90-theta
                ang_bot = 180-(ang_top+90)
                opp = residual_d
                adj = opp/math.tan((ang_bot*(np.pi/180)))
                SurfStrike_d = hole_w - adj
                REVERSE = True

            else:
                ang = theta
                ang = ang*(np.pi/180)
                opp = residual_d
                adj = opp/math.tan(ang)
                SurfStrike_d = adj
                REVERSE = False

            return SurfStrike_d, REVERSE

        
        def ReflectionsInWater(hole_water_d, SZA, theta, t_theta, hole_w, SurfStrike_d, nAir, nWat, REVERSE):
           
            base = hole_water_d / math.tan(t_theta*(np.pi/180))
            n_wat_reflections = 0
            reflect = False

            if REVERSE:

                if base <= SurfStrike_d:

                    n_wat_reflections = 0
                    beam_d_wat = 0

                    reflect = False
                    
                else:

                    reflect = True
                    # depth gained by first subsurface beam
                    beam_d_wat = math.tan(t_theta*(np.pi/180)) * SurfStrike_d
                    
                    while reflect == True:

                        ang_top = 90-t_theta

                        depth_gained = hole_w / math.tan(ang_top * np.pi/180) 

                        beam_d_wat += depth_gained

                        if beam_d_wat >= hole_water_d:
                            
                            reflect = False

                        n_wat_reflections += 1
                               
            
            else:

                reflect = True

                if base <= hole_w - SurfStrike_d:

                    n_wat_reflections = 0
                    beam_d_wat = 0
                    reflect = False


                else:

                    beam_d_wat = math.tan(t_theta*(np.pi/180)) * (hole_w - SurfStrike_d)
                    
                    while reflect == True:

                        ang_top = 90-t_theta
                        depth_gained = hole_w / math.tan(ang_top * (np.pi/180)) 
                        beam_d_wat += depth_gained
                       

                        if beam_d_wat >= hole_water_d:
                            
                            reflect = False

                        n_wat_reflections += 1
            

            return n_wat_reflections, beam_d_wat


       ############################################################


        beamHitsWall, SurfStrike_d, SurfToWat = DoesBeamHitWall(theta, hole_d, hole_water_d, hole_w)
    
        n_air_reflections, residual_d, beam_d_air = ReflectionsInAir(hole_w, theta,SurfToWat, beamHitsWall)

        SurfStrike_d, REVERSE = UpdateSurfStrike_d(n_air_reflections, theta, SurfStrike_d, residual_d)

        n_wat_reflections, beam_d_wat = ReflectionsInWater(hole_water_d, SZA, theta, t_theta, hole_w, SurfStrike_d, nAir, nWat, REVERSE)

        total_reflections = n_air_reflections + n_wat_reflections


        if verbose:

            print("during each reflection the beam descends by {} cm".format(np.round(beam_d_air,2)))
            print("there are {} reflections above the water surface".format(n_air_reflections))
            print("the beam reaches the water surface after {} reflections".format(n_air_reflections))
            print("the beam hits the water {}cm from the nearest wall".format(np.round(SurfStrike_d,2)))
            print("hole width",hole_w)
            print("the beam hits the vertical wall {}cm below the surface".format(np.round(beam_d_wat,2)))
            print("in each subsurface reflection the beam descends by {} cm".format(np.round(beam_d_wat,2)))
            print("the total number of subsurface reflections is {}".format(n_wat_reflections))

                    
        return n_air_reflections, n_wat_reflections, total_reflections, SurfStrike_d, beamHitsWall



    def CalculatePathLength(hole_water_d, hole_w, BeamHitsWall, t_theta, SurfStrike_d, n_wat_reflections, ang_crit):
        
        """
        
        Function calculates the total path length travelled by beam in water

        """
        import math
        import numpy as np

        if t_theta > ang_crit:
            
            if hole_water_d > 0:
                
                ang_top = 90 - t_theta
                beam_d_wat = hole_water_d
                base = beam_d_wat * math.tan(ang_top*(np.pi/180))
                PathLengthInWat = math.sqrt(beam_d_wat**2 + base **2)

            else:

                PathLengthInWat = 0
        
        else:

            if hole_water_d == 0:


                PathLengthInWat = 0
            
            else:

                if BeamHitsWall:
                    # start with first subsurface path (surface strike poition to wall) 

                    beam_d_wat1 = math.tan(t_theta*(np.pi/180)) * SurfStrike_d

                    hyp1 = math.sqrt((beam_d_wat1**2)+(SurfStrike_d**2))

                    # subsequent paths traverse the entire hole width, for the number of reflections
                    # defined in n_watreflections
                    ang_top = 90 - t_theta
                    
                    beam_d_wat2 = hole_w / math.tan(ang_top*(np.pi/180))

                    hyp2 = math.sqrt((beam_d_wat2**2)+(hole_w)**2)

                    # final path is from the last wall reflection to the floor, which may not be 
                    # traverse the full hole width

                    residual_d = hole_water_d - (beam_d_wat1 + (beam_d_wat2 * n_wat_reflections))

                    base = residual_d* math.tan(ang_top*(np.pi/180))

                    hyp3 = math.sqrt((residual_d**2)+(base**2))

                    # total path length is sum of all hypotenuses

                    PathLengthInWat = hyp1 + (hyp2 * n_wat_reflections) + hyp3

                else:

                    ang_top = 90 - t_theta
                    base = hole_water_d * math.tan(ang_top*(np.pi/180))
                    PathLengthInWat = math.sqrt(base**2 + hole_water_d**2)

        return PathLengthInWat


    def AttenuateBeam(PathLengthInWat, kWat, dir_energy_at_hole_floor, WL):
        """
        Function uses imaginary refractive index to calculate absorption coefficient (1/m). The
        product of absorption coefficient and path length (m) is dimensionless attenuation factor used
        to attenuate the beam due to absorption by water

        """
        import numpy as np

        abs_coeff = np.zeros(len(WL))
        
        if PathLengthInWat != 0:

            for i in range(len(WL)):

                abs_coeff[i] = 4*np.pi*kWat[i] / WL[i]

            norm_abs_coeff = abs_coeff * (PathLengthInWat) # multiply abs coeff (/m) by path length in m
            
            energyLoss = dir_energy_at_hole_floor * norm_abs_coeff


        else:

            energyLoss = np.zeros(len(WL))

        
        dir_energy_at_hole_floor = dir_energy_at_hole_floor - energyLoss
        dir_energy_at_hole_floor[dir_energy_at_hole_floor < 0] =0

        return dir_energy_at_hole_floor


    def fresnel(n1, n2, k1, k2, theta):
        
        """
        function computes fresnel reflection coefficients given the difference in real and imaginary 
        refractive indices between two media (n, k) and the illumination angle, theta

        http://www.oceanopticsbook.info/view/surfaces/the_level_sea_surface

        """
        
        import math
        import numpy as np

        if theta == 90:
            print("\nIncident angle is 90 degrees")
            Rf = (((n1-1)**2)+k1**2) / (((n2-1)**2)+k2**2)
        
        else:

            theta = math.radians(theta)
            theta_2 = math.asin(1/n2*(math.sin(theta)))
            
            Rf = 0.5*(
                ((math.sin(theta-theta_2)/math.sin(theta+theta_2))**2) +
                ((math.tan(theta-theta_2)/math.tan(theta+theta_2))**2))

        return Rf


    def internal_reflection(hole_water_d, cryoconite_albedo, WL, nAir, kAir, nWat, kWat, tolerance,\
        dir_energy_at_hole_floor, diffuse_energy_at_hole_floor):

        """
        calculates losses due to upwelling flux being reflected back down into the water column to then
        be partially absorbed by the cyrocontie layer, then repeat n times with diminishing flux

        Calculates total escaped energy after n intenal reflections too 
        
        """
    
        import numpy as np
             

        # define energy arriving at hole floor at first iteration
        energy_arriving_at_floor = dir_energy_at_hole_floor + diffuse_energy_at_hole_floor

        # since we assume the energy is diffuse after interating with cryoconite layer
        # we effectively start thinking of the ystem in a two-stream ay rather
        # than tracing rays. Therefore, path length is the hole water depth.
        path_length = hole_water_d # in cm

        Rf = [] # empty list for R values at range of angles


        # calculate reflected portion of upwelling energy using Fresnel equation
        # loop through angles then average for diffuse flux
        # suppress /0 warning for this calculation (this is expected, and later corrected) 
        np.seterr(divide='ignore', invalid='ignore')

        for theta in np.arange(1,89,1): 
            theta_rad = np.radians(theta)
            Rf.append(((nWat * np.sqrt( 1- (((nWat/nAir)*np.sin(theta_rad))**2)) - nAir*np.cos(theta_rad))/ (nWat * np.sqrt( 1- (((nWat/nAir)*np.sin(theta_rad))**2)) * nAir*np.cos(theta_rad)))**2)
        
        Rf = np.array(Rf)
        Rf[np.isnan(Rf)] = 0.9999999999   # nans are angles > Brewster's angle == total internal reflection
        Rf[Rf>1] = 0.99999   # correct for any invalid values
        diffuse_Rf = np.mean(Rf) # diffuse Fresnel reflection
        
        # turn warnings back on
        np.seterr(divide='warn', invalid='warn')

        # energy upwelling after absorption by cryoconite layer
        upwelling_energy = energy_arriving_at_floor * (1-cryoconite_albedo)
        
        # initiate cumulative vars for total loss and total absorbed by cryoconite
        loss = 0
        cryoconite_abs = 0

        # calculate absorption coefficient of water column
        abs_coeff = 4*np.pi*kWat / WL
        norm_abs_coeff = abs_coeff * (path_length) # multiply abs coeff (/m) by path length in m
        
        # begin loop, one iteration for one internal reflection
        # i.e. each n tracks light from immediately above cryoconite layer t o water surface
        # then back down to the cryoconite layer    
        while np.sum(upwelling_energy) > tolerance:

            loss_upwards = upwelling_energy*norm_abs_coeff # loss as beam passes upwards through water column
            loss_at_boundary = (upwelling_energy - loss_upwards) * (1-diffuse_Rf) #fresnel loss at boundary
            down_flux = upwelling_energy - loss_upwards - loss_at_boundary 
            loss_downwards = down_flux * norm_abs_coeff # loss as beam passes downwards through water column
            loss_at_cryoconite = (down_flux-loss_downwards) * (1-cryoconite_albedo) # loss due to abs by cryoconite

            # total losses in iteration
            total_iteration_loss = loss_upwards + loss_at_boundary + loss_downwards + loss_at_cryoconite

            # update cumulative loss and abs at cryoconite 
            loss += total_iteration_loss
            cryoconite_abs += loss_at_cryoconite

            # update upwelling flux after one complete iteration
            upwelling_energy = upwelling_energy - total_iteration_loss

        # total energy escaping to atmosphere (i.e. contributing to surface albedo) 
        # is original upwelling flux minus total losses after n iterations
        escaped = energy_arriving_at_floor - loss

        return escaped, loss, cryoconite_abs
