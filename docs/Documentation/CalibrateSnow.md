# VIC Model Calibration Using Snow Parameters

1.  MAX_SNOW_TEMP - Maximum air temperature at which snow can fall. Usually set to 1.5 C. This value is set in the global control file.
2.  MIN_SNOW_TEMP - Minimum air temperature at which rain can fall. Usually set to -0.5 C This value is set in the global control file.
3.  SNOW_ROUGH - Snow surface roughness. Values used are in the range of 5 mm to 20 cm. This value have effect on sublimation, sensibel/latent heat fluxes. A higher value gives more sublimation, higher sensible heat fluxes. Increased snow surface roughness can be used to force the snowmelt to come alitle earlyer in the melt periods. It should be noticed that the snow surface roughness is not only the snowgrainsize in the top layer, but is a factor that describes how much irregularityes there are to make turbulent eddies to transport energy.
4.  SNOW ALBEDO - Very important for calculation of the snowmelt is the initial snow albedo and how the decay in albedo is calculated as the snow ages and metamorphosis changes the snow structure. All the snow albedo parameters are set in the file snow.h.

    snow.h default values:
    ```
    #define NEW_SNOW_ALB      0.85
    #define SNOW_ALB_ACCUM_A  0.94
    #define SNOW_ALB_ACCUM_B  0.58
    #define SNOW_ALB_THAW_A   0.82
    #define SNOW_ALB_THAW_B   0.46
    ```

    Initial albedo of new snow is by default set to 0.85\. All the equations describing albedo follow Bras 1990 page 262.

    ![](Bernt/calibration/Image4.gif)

    Where _a_ is new snow albedo, _b_ and _c_ are parameters, _t_ is number of days since last snowfall.

    ![](Bernt/calibration/Image5.gif)

    for the accumulation season

    ![](Bernt/calibration/Image6.gif)
    
    for the melt season

    New snow albedo should be in the range 0.75 to 0.95. There should be no spatial variability in the new snow albedo. These parameters should be used in a consistance way for the whole basin. Do the changes in snow.h and recompile before running the model again.
