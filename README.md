# CryoconiteRTM

## Introduction
This software uses a simple ray-tracing method to calculate the absorption of indicent solar radiation in a cryoconite of given geometry. Two modes are avalable - the simplest is "hole-mode" in which the absorption of energy within a single cryoconite hole is calculated. The user provides the depth, water depth and width (diameter) of a cryoconite hole along with the illumination conditions and surrounding ice density.

The second mode is "surface-mode". Here, the absorption (soon, albedo) of a patch of ice is calculated. The user defines the physical properties of the ice surface and populates it with any number of individual cryoconite holes each with their own dimensions. The total energy absorbed by all the cryoconite holes in the area is calculated.

In either case, the model includes the effect of:

1) a collimated direct beam which can either directly illuminate the hole floor or reach the hole floor after reflecting one or more times from the hole walls. The energy in the beam is subject to losses due to Fresnel reflection at the air/water interface and multiple reflections between the hole walls. 
2) A diffuse flux that illuminates the cryoconite hole floor from above and below after passing through the surrounding ice
3) The incident energy reaching, but not absorbed by, the cryoconite layer becomes a diffuse upwards flux that can reflect back down from the underside of the water surface multiple times.

The total energy absorbed by the hole is the mean across n points across the hole floor, each of which has a different illumination geometry determined by the solar angle and the aspect ratio of the hole (ratio of width to depth). 


## How to Use

### Environment
This software runs best in the provided environment BioSNICAR_py. This can be created from the provided .yaml file using 

```
conda create -f BioSNICAR_py.yaml

```

The user must also have the repository BioSNICAR_GO_PY downloaded onto their
computer. There are currently some hard-coded paths to that repository in this
code that I will try to remove later in the development.


### Running the model
The model is run from "driver.py" from the terminal. In-script annpotations show clearly the user-defined variables whose values can be changed. I have moved all derived and hard-coded variable setting to external scripts. Simply set the hole geometry and illumination conditions as required and run the script. The output is a plot showing the total incoming irradiance and the spectral energy absorbed by the cryoconite layer at the hole floor.

## Background

### Theory
This model calculates the energy absorbed by layers of cryoconite at the base of cryoconote holes of known geometry and under known illumination conditions. This is achieved in two stages: 1) calculation of illumination by a collimated beam assumed to arrive from a point source at zenith angle SZA. This can be direct illumination or illumination after one or more reflections from the hole walls. 2) Energy arriving from all directions due to multiple scattering in the ice around the cryoconite hole. This is calculated by calculating using a two-stream model based on SNICAR. The energy at hole depth d in a simulated homogenous layer of ice is added to the collimated beam flux. The total energy (collimated beam + diffuse) is multiplied by the cryoconite albedo to calculate the energy absorbed by the cryoconite sediment at the hole floor.



Schematic of simple no-reflection case 

![HoleGeom1](/Assets/cryoconite_geometry1.jpg)



Schematic of more comples multiple reflection case 

![HoleGeom2](/Assets/cryoconite_geometry2.jpg)


The hole is considered as a 2D idealised parallel and vertical-sided and flat floored shape with a continuous and uniform layer of cryoconite covering the floor, overlain by a column of water and assumed to be azimuthally-symmetrical. The direct beam might illuminate the hole floor directly. To determine whether the floor is directly illuminated, first there is a test to determine whether the beam hits the water surface or the hole wall. If the beam does not hit the water surface, it may hit the hole wall above the water surface. If the beam hits the wall, the beam is considered to create a right triangle with angle theta and base of length equal to the hole width. The depth gained by the beam when it reflects is the vertical side of that triangle that can be solved for using basic trigonometry. If that depth is less than the dostance between ice surface and water surface, another reflection occurs. This can be repeated until the beam reaches the water surface. Another trigonometric calculation provides the distance from the sunwards wall that the beam hits the water surface. Losses occur at each reflection between air and ice and at the entry to the water, according to Fresnel reflection. When the beam enters the water it is refracted, so the incoming angle is adjusted to a new transmitted angle.

If the beam does hit the water surface without hitting the hole wall, there is a loss due to fresnel reflection, then the beam is refracted according to Snell's law. 

Once the beam enters the water, it has a new angle of incidence due to refraction, and a known position on the water surface that the beam entered the water column. There is then a chance that the beam illuminates the hole floor without reflecting from any other surfaces. In this case, the only other loss is due to absorption in the water column, according to the path length and the absorption coefficient of water. However, the beam may hit a wall before hitting the floor. In this case, there is an additional Fresnel loss each time the beam reflects from the hole wall.

Eventually, after n reflections, the beam strikes the hole floor. The amount of energy absorbed by cryoconite is the remaining energy after the incoming energy at the surface has been reduced at each reflective loss and by absorption in the water (taking into account the path length after multiple reflections), multiplied by 1 - albedo of the cryoconite.

There is also energy arriving from the surrounding ice, since light incident around the cryoconite hole that does not enter the hole scatteres multiple times and penetrates the ice to some depth, likely coming into contact with the cryoconite layer where is has a strong chance of being absorbed. This is accounted for using a two-stream radiative transfer calculation where the ice is assumed to be a homogenous layer with depth equal to the hole depth. The ice configuration is user-defined. The incoming irradiance is identical to that used for the collimated beam. The actinic flux at depth = hole_depth is then added to the energy field and actinic flux * 1-cryoconite albedo is added to the absorbed energy.

The energy reaching the cryoconite but not absorbed is then considered to be an upwelling diffuse flux. This then interacts with the underside od the water surface, some of which will be transmitted ("escaping" the hole) and some of which will be reflected back down into the hole as a downwards flux. This happens a infinite number of times with diminishing flux in each iteration. In the model the user defines a threshold (default 1e-10) for the upwelling flux - when the upwelling flux is smaller than the threshold the iteration for internal reflections ends.

The sum of all three fluxes gives the total "escaped" and total "absorbed". These can be expressed in Watts, Watts/m2 or as a percentage of the incoming flux.


## Example output

This plot shows the spectral energy absorbed at the hole floor compared to the incoming irradiance for a single hole.

![ExOutput](/Assets/Out.jpg)

## Testing

Two types of testing have been undertaken - a) unit testing, where the target data for specific functions are gathered from known theory, and b) validation testing, where the model predictions are compared against empirical field measurements. 

### Unit testing
The key functions tested are the Fresnel reflection calculations, the Snell's Law calculations and the multiple reflection functions. Each of these have well-known theory underpinning them. In the case of the fresnel calculations at an air/ice boundary, there is a known relationship between illumination angle and the magnitude of fresnel reflection with a curved shape that varies slightly with wavelength. In the Snell's law calculations, the relationship between incident and transmitted angles are also easily calculated using an external script, again, varying with wavelength. For the multiple reflections, the test is really just with logic - the questions I asked of the data were: does the number of multiple reflections increase when the solar angle is more oblique? Does the number of multiple reflections increase when the hole becomes deeper and/or narrower? Does the ratio of in-air to in-water reflections change sensibly as the in-hole water depth increases and decreases? Are the patterns of in-air reflections and in-water reflections consistent?

### Representative Test Outputs:
#### Fresnel Reflection Functions

![FresnelTests](/Assets/FresnelTests.jpg)

#### Transmitted Angle Function

![TransAngle](/Assets/TransAngleTests.jpg)

#### Multiple Reflection Function

![MultipleReflectionTests](/Assets/MultipleReflectionsTests.jpg)


### Validation tests
The validation tests compared predicted values against field measured values. The field measurements are PAR-pyranometers that were simultaneously positioned on cryoconite hole floor and on the ice surface immediately adjacent to the hole edge. The hole width and depth was measured for each hole. The sensor was 4cm tall, meaning the pyranometer aperture sat 4cm above the hole floor, so 4cm was subtracted from the hole depth when defining the simulated hole depth. The values for SZA were estimated from the measurement time and location, the ice grain size and density were taken from literature values for nearby sites. The values used were SZA = 20, densoity = 700, bubble radius = 700, hole water depth = 0.7 * hole depth. The ratio between the surface and hole floor irradiance is used as the target data to be simulated by the model. The absolute error between measured and modelled irradiance is plotted below. Since the field measurements were made with pyranometers, the simulated irradiance was limited to the wavelength range 400-700 nm. Under these conditions, the mean absolute error was 0.10 +/- 0.08 (1SD).

![ValidationTests](/Assets/ValidationTests.jpg)


## TODOs and next steps

### Albedo

At the moment the model calculates energy absorption by the various surface components. The next step will be to use this information to calculate spectral and broadband albed for a surface. It might also be possible to use UAV images as input data for cryoconite hole coverage.

### Field validation

The current validation has been achieved by comparison with field measurements made with a PAR-sensitive pyranometer lowered into cryoconite holes of known width/depth during a field expedition to SW Greenland in 2010.
These data were not collected specifically for the purpose of validating this model, so ther erae some metadata missing that have had to be infilled by inference, e.g. representative SZA estimated from measurement time and location, density and grain sizes for surrounding ice estimated from literature values. This model was written in 2020, and the covid19 situation has made the prospects of dedicated field work for validating this model extremely unlikely, but in the future a dedicated field campaign for model validation would be very useful.


## Permissions

This code is provided with absolutely no warranty of any kind. I am still actively developing this repository and advise against using it for any downstream purpose. Collaboration requests to jcook@envs.au.dk or @tothepoles (Instagram/Twitter).