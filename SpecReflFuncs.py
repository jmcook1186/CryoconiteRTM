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


    def trans_angle(SZA, nAir, nWat):
        
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
        t_theta = []

        for i in range(len(nAir)):
            
            sin_ang = nAir[i] * math.sin(SZA/180*np.pi) / nWat[i]
            out_ang = math.asin(sin_ang)*180/np.pi
            t_theta_temp = 90 - out_ang # 180 - (out_ang+90)

            check = t_theta_temp + out_ang

            if abs(90-check) > 0.0001:
            
                print("\nSNELL ANGLE CALCULATION: TEST FAILED")
        
            else:
                pass

            t_theta.append(t_theta_temp)

        return t_theta


    def test_multiple_reflections(theta, t_theta, hole_d, hole_w, hole_water_d, verbose = False):
        """
        function tests whether a single reflection from the hole wall causes the beam to hit the hole floor.
        If not, how many reflections between the walls occur before the beam reaches the hole floor?
        """

        import numpy as np
        import math


        # 1) use incoming angle (theta) and the vertical distance between the hole aperture and the water surface (adj)
        # to calculate the horizontal distance from the sunward-wall to the beam-strike point (opp) using trigonometric relation:
        # tan(angle) = opp/adj  
        
        ang1 = theta
        ang1_rad = ang1*(np.pi/180) # radians

        SurfToWat = hole_d - hole_water_d
        opp1 = SurfToWat*(math.tan(ang1_rad))

        # opp1 in this case represents the distance from the sunwards wall that the beam strikes the water surface
        # for very oblique angles or narrow holes this might be greater than the hole width - in this case
        # the beam hits the vertical wall above the water surface and at least one reflection before
        # reaching the water surface.

        # Each reflection from the vertical walls will be at angle theta and the beam depth in the hole
        # will increase by a constant amount during each reflection until the water depth is reached. We
        # need to know the number of reflections occurring above the water surface and the horizontal 
        # position of the beam when it eventually strikes the water. Because of the law of reflection, the angle
        # theta will always be the angle of incidence regardless of the number of reflections occurring.

        if opp1 >= hole_w: #i.e. if the beam strikes the wall above the water

            adj1 = hole_w # take the hole width as a known length adjacent to angle 180-(theta+90) 
            ang1 = 180-(theta+90)
            ang1_rad = ang1*(np.pi/180) # radians
            opp1 = adj1*(math.tan(ang1_rad))  # opp = tan(angle)*adj to find vertical distance gained in single reflection
            beam_d_air = opp1   # re-alias for clarity

            if verbose:
                print("during each reflection the beam descends by {} cm".format(np.round(beam_d_air,2)))
            
            n_air_reflections = np.floor((SurfToWat)/(beam_d_air)) # number of reflections in air before hitting water surface
            
            residual_d = SurfToWat - (n_air_reflections*beam_d_air) # SurfToWat is hole_d - water_d, n_refl*beam_d is total
            # depth acounted for by complete reflections, residual is depth gained between final reflection and
            # striking water surface
            
            # calculate horizontal distance along water surface of beam strike after reflections - i.e. use residual_d as
            # known vertical distance and angle theta as known angle to calculate horizontal distance

            ang_top = 90-theta
            ang_bot = 180-(ang_top+90)
            opp = residual_d
            adj = opp/math.tan((ang_bot*(np.pi/180)))  
            horizontal_strike_d = adj

            if verbose:
                print("there are {} reflections above the water surface".format(n_air_reflections))
                print("the beam reaches the water surface after {} reflections".format(n_air_reflections))
                print("the beam hits the water {}cm from the nearest wall".format(np.round(horizontal_strike_d,2)))

            # now subtract horizontal_strike_d from hole width to get the width of the triangle formed by the subsurface
            # beam. We have two angles in the subsurface triangle - t_theta and the right angle. Therefore we
            # have one side and all three angles once the final angle has been calculated (angles add to 180).
            # We can therefore calculate the length of the missing side (opp) using tan_ang = opp/adj
            # opp2 is the depth beneath the water surface where the beam strikes the wall

            adj2 = hole_w - horizontal_strike_d
            ang2 = 180-(t_theta+90)
            ang2_rad = ang2*(np.pi/180)
            opp2 = adj2*(math.tan(ang2_rad))

            beam_d_wat = opp2
            
            if verbose:
                print("the beam hits the vertical wall {}cm below the surface".format(np.round(beam_d_wat,2)))
                print("in each subsurface reflection the beam descends by {} cm".format(np.round(beam_d_wat,2)))

            n_wat_reflections = np.floor(hole_water_d/(beam_d_wat)) # nmber of reflections in air before hitting water surface

            if verbose:
                print("the total number of subsurface reflections is {}".format(n_wat_reflections))

            residual_d = SurfToWat - (n_wat_reflections*beam_d_wat) #adj1 is hole_d - water_d, n_refl*beam_d is total
            # depth acounted for by complete reflections, residual is depth gained between final reflection and
            # striking water surface
            #
            ang_top = 90-t_theta
            ang_bot = 180-(ang_top+90)
            opp = residual_d
            adj = opp/math.tan((ang_bot*(np.pi/180)))  
            floor_strike_d = adj

        else:
            
            # if no reflections above water surface, jump straight to calculating suburface depth for beam striking hole wall
            
            adj2 = hole_w - opp1
            ang2 = 180-(t_theta+90)
            ang2_rad = ang2*(np.pi/180)
            opp2 = adj2*(math.tan(ang2_rad))

            beam_d = opp2

            n_air_reflections = 0

            if verbose:
                print("there are no reflections above the water surface")
                print("the beam hits the vertical wall {}cm below the surface".format(np.round(beam_d,2)))

            # after the beam hits the wall it reflects back again at angle theta and gains depth
            # if the cumulative depth is less than the hole depth it reflects again until the hole depth is reached
            # The total number of complete reflections (n_reflections) is the number of loss terms to include. The
            # number of reflections multiplied by the depth gained per reflection, subtracted from the hole depth gives
            # the residual depth - the height of the triangle formed by the beam on the final reflection and the hole floor.
            # trigonometry allows us to solve for the length of the horizontal, giving the position on the hole floor that
            # the beam strikes.

            if verbose:
                print("in each subsurface reflection the beam descends by {} cm".format(np.round(beam_d,2)))
            
            n_wat_reflections = np.floor(hole_water_d/beam_d) # number of reflections in air before hitting water surface
    
            if verbose:
                print("the total number of subsurface reflections is {}".format(n_wat_reflections))

            
            residual_d = hole_d - (n_wat_reflections*beam_d) #adj1 is hole_d - water_d, n_refl*beam_d is total
            # depth acounted for by complete reflections, residual is depth gained between final reflection and
            # striking water surface
            
            ang_top = 90-t_theta
            ang_bot = 180-(ang_top+90)
            opp = residual_d
            adj = opp/math.tan((ang_bot*(np.pi/180)))  
            floor_strike_d = adj

        total_reflections = n_air_reflections + n_wat_reflections
                    
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



import math

math.asin(1/1.33*(math.sin(10)))