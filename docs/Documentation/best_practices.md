# Introduction
The VIC model is a computation tool that we use to answer hydrologic sciences question. As is the case with any tool, there are proper and improper ways to use it. This page highlights a few "best practices" that will help VIC users avoid common pitfalls. For more information on the motivation and background for these best practices, see Wilson et al. (2014).

## Use the computer to record history
VIC includes tools for saving runtime log files. VIC writes status updates, warnings, and errors to the user during runtime. This can be done by setting the `LOG_DIR` variable in the global parameter file ([example](./Drivers/Classic/GlobalParam.md)). These files are often useful for debugging problems with VIC simulations and serve as a way to document completed VIC runs. We recommend keeping these log files after each simulation.

## Output files should have informative names
Most of the VIC drivers allow users to specify input and output filenames. For example, the Classic Driver uses the `OUTFILE {prefix}` syntax in the global parameter file, where `{prefix}` represents an the first part of an arbitrary output filename. Although VIC defaults to prefix values of `fluxes`, `snow`, or `snowbands`, these are not particularly informative when they are used repeatedly. We suggest using prefixes that indicate the important parameters of the VIC simulation. For example a prefix of `vic_5.0.0_erai_3h_fullengy_snow` may represent snow variables from a VIC version 5.0.0 simulation using 3-hourly Era-Interim forcings with `FULL_ENERGY == TRUE`. There are obviously many other ways to specify simulation metadata in the input and output files.

## Use Git for version control
We strongly recommend using the Git version control software when modifying VIC. Using Git will allows users to save and document changes made to the model source code. Highlights of what Git offers a VIC user/developer:

  - document changes to the VIC source
  - seamlessly switch between model versions
  - share changes/improvements to the model (e.g. Github)
  - backup model versions (e.g. Github)

For more information on using Git with VIC, see our [Git tutorial](../Development/working-with-git.md).

#### References

1. Wilson, G., Aruliah, D.A., Brown, C.T., Hong, N.P.C., Davis, M., Guy, R.T., Haddock, S.H., Huff, K.D., Mitchell, I.M., Plumbley, M.D. and Waugh, B., 2014. Best practices for scientific computing. _PLoS Biol_, *12(1)*, p.e1001745, [10.1371/journal.pbio.1001745](http://dx.doi.org/10.1371/journal.pbio.1001745).
