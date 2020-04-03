
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



# WHICH FUNCTIONS TO TEST?
# check_fresnel(nAir,nWat,kAir,kWat,WL,plot_figs=True)
# check_multiple_reflections(nAir,nWat)
check_trans_angle(nAir, nWat)

