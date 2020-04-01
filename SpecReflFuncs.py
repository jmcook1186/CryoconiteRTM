"""
Class specFuncs contains functions for specular reflection and cryoconite hole geometry calculations

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
       
        return ang_crit, hyp


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
            """dis
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
                        
                        n_wat_reflections += 1

                        if beam_d_wat >= hole_water_d:
                            
                            reflect = False

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
                        n_wat_reflections += 1

                        if beam_d_wat >= hole_water_d:
                            
                            reflect = False


            floor_strike_d = 0


            return n_wat_reflections, beam_d_wat, floor_strike_d



        ############################################################


        beamHitsWall, SurfStrike_d, SurfToWat = DoesBeamHitWall(theta, hole_d, hole_water_d, hole_w)
    
        n_air_reflections, residual_d, beam_d_air = ReflectionsInAir(hole_w, theta,SurfToWat, beamHitsWall)

        SurfStrike_d, REVERSE = UpdateSurfStrike_d(n_air_reflections, theta, SurfStrike_d, residual_d)

        n_wat_reflections, beam_d_wat, floor_strike_d = ReflectionsInWater(hole_water_d, SZA, theta, t_theta, hole_w, SurfStrike_d, nAir, nWat, REVERSE)

        total_reflections = n_air_reflections + n_wat_reflections


        if verbose:
            print("during each reflection the beam descends by {} cm".format(np.round(beam_d_air,2)))
            print("there are {} reflections above the water surface".format(n_air_reflections))
            print("the beam reaches the water surface after {} reflections".format(n_air_reflections))
            print("the beam hits the water {}cm from the nearest wall".format(np.round(SurfStrike_d,2)))
            print("hole width",hole_w)
            print("SurfStrike_d = ",SurfStrike_d)
            print("the beam hits the vertical wall {}cm below the surface".format(np.round(beam_d_wat,2)))
            print("in each subsurface reflection the beam descends by {} cm".format(np.round(beam_d_wat,2)))
            print("the total number of subsurface reflections is {}".format(n_wat_reflections))

                    
        return n_air_reflections, n_wat_reflections, total_reflections, floor_strike_d



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

