# CryoconiteRTM

Codes to calculate absorption of energy in a cryoconite hole under different hole geometry and illumination conditions. This model includes the effect of a collimated direct beam which can either directly illuminate the hole floor or reach the hole floor after reflecting one or more times from the hole walls. The energy in the beam is
subject to losses due to Fresnel reflection at the air/water interface and multiple reflections between the hole walls.

In the current form, this model works on the assumption that the conditions at a point in the centre of the hole floor is representative of entire hole floor (i.e. if centre is not directly illuminated none of the floor is directly illuminated, if centre of hole is directly illuminated so is the rest of the floor).

# How to Use

## Environment
This software runs best in the provided environment BioSNICAR_py. This can be created from the provided .yaml file using 

```
conda create -f BioSNICAR_py.yaml

```

The user must also have the repository BioSNICAR_GO_PY downloaded onto their
computer. There are currently some hard-coded paths to that repository in this
code that I will try to remove later in the development.


## Running the model
The model is run from "driver.py" from the terminal. In-script annpotations show clearly the user-defined variables whose values can be changed. I have moved all derived and hard-coded variable setting to external scripts. Simply set the hole geometry and illumination conditions as required and run the script. The output is a plot showing the total incoming irradiance and the spectral energy absorbed by the cryoconite layer at the hole floor.

# Background

## Theory
This model calculates the energy absorbed by layers of cryoconite at the base of cryoconote holes of known geometry and under known illumination conditions. This is achieved in two stages: 1) calculation of illumination by a collimated beam assumed to arrive from a point source at zenith angle SZA. This can be direct illumination or illumination after one or more reflections from the hole walls. 2) Energy arriving from all directions due to multiple scattering in the ice around the cryoconite hole. This is calculated by calculating using a two-stream model based on SNICAR. The energy at hole depth d in a simulated homogenous layer of ice is added to the collimated beam flux. The total energy (collimated beam + diffuse) is multiplied by the cryoconite albedo to calculate the energy absorbed by the cryoconite sediment at the hole floor.



Schematic of simple no-reflection case 

![HoleGeom1](/Assets/cryoconite_geometry1.jpg)



Schematic of more comples multiple reflection case 

![HoleGeom2](/Assets/cryoconite_geometry2.jpg)




The hole is considered as a 2D idealised parallel and vertical-sided and flat floored shape with a continuous and uniform layer of cryoconite covering the floor, overlain by a column of water. The direct beam might illuminate the hole floor directly. To determine whether the floor is directly illuminated, first there is a test to determine whether the beam hits the water surface or the hole wall. If the beam does not hit the water surface, it may hit the hole wall above the water surface. If the beam hits the wall, the beam is considered to create a right triangle with angle theta and base of length equal to the hole width. The depth gained by the beam when it reflects is the vertical side of that triangle that can be solved for using basic trigonometry. If that depth is less than the dostance between ice surface and water surface, another reflection occurs. This can be repeated until the beam reaches the water surface. Another trigonometric calculation provides the distance from the sunwards wall that the beam hits the water surface. Losses occur at each reflection between air and ice and at the entry to the water, according to Fresnel reflection. When the beam enters the water it is refracted, so the incoming angle is adjusted to a new transmitted angle.

If the beam does hit the water surface without hitting the hole wall, there is a loss due to fresnel reflection, then the beam is refracted according to Snell's law. 

Once the beam enters the water, it has a new angle of incidence due to refraction, and a known position on the water surface that the beam entered the water column. There is then a chance that the beam illuminates the hole floor without reflecting from any other surfaces. In this case, the only other loss is due to absorption in the water column, according to the path length and the absorption coefficient of water. However, the beam may hit a wall before hitting the floor. In this case, there is an additional Fresnel loss each time the beam reflects from the hole wall.

Eventually, after n reflections, the beam strikes the hole floor. The amount of energy absorbed by cryoconite is the remaining energy after the incoming energy at the surface has been reduced at each reflective loss and by absorption in the water (taking into account the path length after multiple reflections), multiplied by 1 - albedo of the cryoconite.

There is also energy arriving from the surrounding ice, since light incident around the cryoconite hole that does not enter the hole scatteres multiple times and penetrates the ice to some depth, likely coming into contact with the cryoconite layer where is has a strng chance of being absorbed. This is accounted for using a to-stream radiative transfer calculation where the ice is assumed to be a homogenous layer with depth equal to the hole depth. The ice configuration is user-defined, by default there are no surface impurities but his can easily be updated if necessary. The incoming irradiance is identical to that used for the collimated beam. The actinic flux at depth = hole_depth is then added to the energy field and actinic flux * 1-cryoconite albedo is added to the absorbed energy.


# Testing
Two types of testing have been undertaken - a) unit testing, where the target data for specific functions are gathered from known theory, and b) validation testing, where the model predictions are compared against empirical field measurements. 

## Unit testing
The key functions tested are the Fresnel reflection calculations, the Snell's Law calculations and the multiple reflection functions. Each of these have well-known theory underpinning them. In the case of the fresnel calculations at an air/ice boundary, there is a known relationship between illumination angle and the magnitude of fresnel reflection with a curved shape that varies slightly with wavelength. In the Snell's law calculations, the relationship between incident and transmitted angles are also easily calculated using an external script, again, varying with wavelength. For the multiple reflections, the test is really just with logic - the questions I asked of the data were: does the number of multiple reflections increase when the solar angle is more oblique? Does the number of multiple reflections increase when the hole becomes deeper and/or narrower? Does the ratio of in-air to in-water reflections change sensibly as the in-hole water depth increases and decreases? Are the patterns of in-air reflections and in-water reflections consistent?

## Representative Test Outputs:
### Fresnel Reflection Functions

![FresnelTests](/Assets/FresnelTests.jpg)

### Transmitted Angle Function

![TransAngle](/Assets/TransAngleTests.jpg)

### Multiple Reflection Function

![MultipleReflectionTests](/Assets/MultipleReflectionsTests.jpg)


## Validation tests
The validation tests compared predicted values against field measured values. The field measurements are PAR-pyranometers that were simultaneously positioned on cryoconite hole floor and on the ice surface immediately adjacent to the hole edge. The hole width and depth was measured for each hole. The sensor was 4cm tall, meaning the pyranometer aperture sat 4cm above the hole floor, so 4cm was subtracted from the hole depth when defining the simulated hole depth. The values for SZA were estimated from the measurement time and location, the ice grain size and density were taken from literature values for nearby sites. The values used were SZA = 20, densoity = 600, grain radius = 2000, hole water depth = 0.5 * hole depth. The ratio between the surface and hole floor irradiance is used as the target data to be simulated by the model. The absolute error between measured and modelled irradiance is plotted below. Since the field measurements were made with pyranometers, the simulated irradiance was limited to the wavelength range 400-700 nm. Under these conditions, the mean absolute error was 0.10 +/- 0.08 (1SD).

![ValidationTests](/Assets/ValidationTests.jpg)


# TODOs and next steps

## Field validation
The current validation has been achieved by comparison with field measurements made with a PAR-sensitive pyranometer lowered into cryoconite holes of known width/depth during a field expedition to SW Greenland in 2010.
These data were not collected specifically for the purpose of validating this model, so ther erae some metadata missing that have had to be infilled by inference, e.g. representative SZA estimated from measurement time and location, density and grain sizes for surrounding ice estimated from literature values. This model was written in 2020, and the covid19 situation has made the prospects of dedicated field work for validating this model extremely unlikely, but in the future a dedicated field campaign for model validation would be very useful.

## Spatial mapping
Currently the model assumes the conditions at the centre of the hole are representative of the conditions across the hole floor. A major future development would be to distribute te energy realistically across the hole floor. This might be achievable by using a sequence 1-N representing evenly spaced points across the diameter of the hole in the principle solar plane, rather than a single point in the centre of the hole, and determining the critical angle required for direct illumination for each position. This would allow a proportion of directly illuminated positions to be calculated. Then, the multiple reflections of the direct beam could be calculated for all points on the sun-facing hole wall, with the output being not only the energy in the beam when it hits the floor, but the position along the transect that the beam hits. Then, the magnitude of the beam could be added to the energy in that position. The diffuse flux wcould be assumed to be uniform across the transect. Then, there is some challenge to integrate the values in the transect around 360 degrees of the hole.


# Permissions

This code is provided with absolutely no warranty of any kind. I am still actively developing this repository and advise against using it for any downstream purpose. Collaboration requests to joc102@aber.ac.uk or @tothepoles (Instagram/Twitter).