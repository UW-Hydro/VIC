# Notes on VIC Modes of Operation

A (large) number of options exist for tailoring VIC simulations to your application. Most of these are set in the [global parameter file](GlobalParam.md) (run-time options), while others are set in the file vicNl_def.h and require re-compiling VIC to take effect (compile-time options). The most commonly used options are FULL_ENERGY (which determines whether VIC will iterate to find the surface temperature that balances the energy budget) and FROZEN_SOIL (which determines whether VIC simulates the phase change of water in the soil). The notes below show some comparisons of the effects of the different options on model performance, in terms of both results and speed.

## *   [Modes of Operation](TechnicalNotes/Modes.md)

## *   [Time Step Comparisons](TechnicalNotes/Timestep.md)

## *   [Frozen Soil Time Step Comparisons](TechnicalNotes/Timestepfrozsoil.md)

## *   [Ground Heat Flux Comparisons](TechnicalNotes/GroundHeatFlux.md)

## *   [VIC/GOES Radiation Evaluation](http://ftp.hydro.washington.edu/pub/niklas/VIC_radiation_evaluation_010703.pdf)

## *   [Stability Correction for Aerodynamic Resistance](TechnicalNotes/Stability.md)

## Archived Notes

## *   [Solar Radiation](TechnicalNotes/Radiation.md)
