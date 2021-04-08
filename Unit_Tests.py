
import numpy as np
import math
import matplotlib.pyplot as plt
from SpecReflFuncs import specFuncs


# SET UP SOME CONSTANTS NEEDED TO RUN FUNCS
WL = np.arange(0.3,5,0.01)

# define n and k for air (array of ones)
nAir = np.ones(shape=(470))
kAir = np.zeros(shape=(470))

# import spectral refractive index for water
nWat = np.genfromtxt('/home/joe/Code/CryoconiteRTM/Data/water_n.csv', delimiter=",")
kWat = np.genfromtxt ('/home/joe/Code/CryoconiteRTM/Data/water_k.csv', delimiter=",")
kWat = kWat[0:-1:10] #every 10th element to match resolution of SNICAR

nIce = np.genfromtxt('/home/joe/Code/CryoconiteRTM/Data/ice_n.csv', delimiter=",")
kIce = np.genfromtxt ('/home/joe/Code/CryoconiteRTM/Data/ice_k.csv', delimiter=",")


# CHECK FRESNEL FUNCTION
#run function and compare against benchmark data

def check_fresnel(nAir,nWat,kAir,kWat,WL,plot_figs=True):

    """
    Checks the Fresnel reflection calculations aainst benchmark data

    """

    FresnelTargetData = np.genfromtxt('/home/joe/Code/CryoconiteRTM/TestData/FresnelTargetData.csv')
    test = np.zeros(shape=(FresnelTargetData.shape))
    passes = 0
    fails = 0


    for i in np.arange(1,90,1):
        for j in range(len(WL)):
            SZA = i
            R = specFuncs.fresnel(nAir[j],nWat[j],kAir[j],kWat[j],i)
            test[i,j] = R
    
    

    for i in range(FresnelTargetData.shape[0]):

        for j in range(FresnelTargetData.shape[1]):
            
            if test[i,j] - FresnelTargetData[i,j] > 0.00001:
                fails +=1
                
            else:
                passes +=1
    
    if plot_figs:
        plt.plot(test)
        plt.title('Fresnel reflection: each curve is a different wavelength')
        plt.ylabel('Fresnel Rf (proportion of incident reflected)')
        plt.xlabel('Solar zenith angle (0 = overhead, 90 = horizon')
        plt.show()

    return print("\n*** {}/{} FRESNEL TESTS PASSED WITH TOLERANCE 1E-5***\n".format(passes,passes+fails))



def check_trans_angle(nAir, nWat):

    """

    Checks that the transmitted angle s correctly by plotting for manual verification

    """
    
    WL = np.arange(0.3,5,0.01)
    outList = []
    thetaList = []

    for theta in np.arange(0,90,3):

        t_theta = specFuncs.trans_angle(theta,nAir,nWat)
    
        plt.plot(WL,t_theta, label='Theta = {}'.format(theta))
    
    plt.title('Transmitted angle (each curve is a different incoming SZA)')
    plt.ylabel('Transmitted angle (elevation from horizontal)')
    plt.xlim(0.3,5)
    plt.xlabel('Wavelength (microns)')
    plt.legend(bbox_to_anchor=(1, 1), loc='upper left', ncol=2)
    plt.tight_layout()
    plt.show()

    return


def check_multiple_reflections(nAir,nWat):
    
    """
    checks multiple reflections between hole walls by plotting for manual verification

    """

    n_air_reflections_list = []
    n_wat_reflections_list = []
    total_reflections_list = []
    floor_strike_d_list = []
    t_theta_list = []
    
    for theta in np.arange(1,90,1):

        hole_d = 60
        hole_w = 30
        hole_water_d = 40
        point = hole_w/2

        t_theta = specFuncs.trans_angle(theta,nAir,nWat)

        print("\nelevation angle = ",theta," transmitted angle (400 nm) =  ", t_theta[20])

        n_air_reflections, n_wat_reflections, total_reflections, SurfStrike_d, \
        beamHitsWall = specFuncs.test_multiple_reflections(
            theta, t_theta[0], hole_d, hole_w, hole_water_d, nAir, nWat, verbose = False)

        n_air_reflections_list.append(n_air_reflections)
        n_wat_reflections_list.append(n_wat_reflections)
        total_reflections_list.append(total_reflections)
        
    plt.plot(np.arange(1,90,1), n_air_reflections_list,label='# air refls')
    plt.plot(np.arange(1,90,1), n_wat_reflections_list, label='# wat refls')
    plt.title('Multiple reflections for hole depth 60 cm, hole width 30cm, hole water depth 40 cm')
    plt.xlabel('Solar Elevation Angle (degrees)')
    plt.ylabel('Number of reflections')
    #plt.plot(total_reflections_list, label='# total refls')
    plt.legend(loc='best')

    plt.show()

    return


def check_internal_reflections(nAir, kAir, nWat, kWat, WL):

    """
    checks internal reflections are correctly calculated using several assertions
    
    """

    from SpecReflFuncs import specFuncs

    hole_water_d = 20
    cryoconite_albedo = 0.2
    tolerance =1e-10
    out = []

    # reference values for dir & diff energy arriving at hole floor and diff energy absirbed by cconite
    dir_energy_at_hole_floor = np.genfromtxt('/home/joe/Code/CryoconiteRTM/TestData/dir_energy_at_floor_unit_test.csv')
    diffuse_energy_at_hole_floor = np.genfromtxt('/home/joe/Code/CryoconiteRTM/TestData/diffuse_energy_at_floor_unit_test.csv')
    diffuse_energy_absorbed_by_cryoconite = np.genfromtxt('/home/joe/Code/CryoconiteRTM/TestData/diffuse_energy_absorbed_by_cryoconite.csv')

    total_energy_arriving_at_floor = np.sum([diffuse_energy_at_hole_floor,dir_energy_at_hole_floor])

    # test changing water depth
    for hole_water_d in np.arange(1,20,1): # only go to 10 because that was used as hole-water_d to generate flux files

        energy_escaping_internal_reflections = []
        energy_lost_internal_reflections = []

        for i in range(len(WL)): 

            escaped, loss, cryoconite_abs = specFuncs.internal_reflection(hole_water_d, cryoconite_albedo, WL[i], nAir[i], kAir[i], nWat[i], kWat[i], tolerance,\
                dir_energy_at_hole_floor[i], diffuse_energy_at_hole_floor[i])

            energy_escaping_internal_reflections.append(escaped)
            energy_lost_internal_reflections.append(loss)

        energy_lost= np.sum(energy_escaping_internal_reflections) + np.sum(energy_lost_internal_reflections)
    
        # CHECK 1) all energy arriving at hole floor must either escape or be absorbed, if not, conservation of energy violation
        assert abs(energy_lost - total_energy_arriving_at_floor) < tolerance, \
            f"hole_water_d= {hole_water_d}, energy_lost ={energy_lost}: CONSERVATION OF ENERGY ERROR"

    out.append(energy_lost)

    # CHECK 2) as water height increases, energy loss should increase (by small amount) - i.e. out should be monononically increasing
    assert all(x<y for x,y in zip(out, out[1:])), "increasing water depth leads to less energy loss"

    # NOTE: these tests do not provide full coverage - to fully test requires
    # iteration through hole depths, hole widths, solar zeniths...

    # if we make it through the assertions we can print a success statement
    print("*** Internal reflection unit tests passed successfully ***")

    return

# WHICH FUNCTIONS TO TEST?
# check_fresnel(nAir,nWat,kAir,kWat,WL,plot_figs=True)
# check_multiple_reflections(nAir,nWat)
# check_trans_angle(nAir, nWat)
check_internal_reflections(nAir,kAir,nWat,kWat,WL)
