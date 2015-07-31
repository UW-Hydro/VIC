# Variable Infiltration Capacity (VIC) Macroscale Hydrologic Model

[![Build Status](https://travis-ci.org/UW-Hydro/VIC.png?branch=develop)](https://travis-ci.org/UW-Hydro/VIC) [![VIC Users Listserve](https://img.shields.io/badge/VIC%20Users%20Listserve-Active-blue.svg)](https://mailman.u.washington.edu/mailman/listinfo/vic_users) [![Join the chat at https://gitter.im/UW-Hydro/VIC](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/UW-Hydro/VIC?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge) [![GitHub license](https://img.shields.io/badge/license-GPLv2-blue.svg)](https://raw.githubusercontent.com/UW-Hydro/VIC/master/LICENSE.txt) [![Documentation Status](https://readthedocs.org/projects/vic/badge/?version=latest)](http://vic.readthedocs.org/en/latest/) [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.22307.svg)](http://dx.doi.org/10.5281/zenodo.22307)

VIC ([Liang et al., 1994](Documentation/References.md)) is a macroscale hydrologic model that solves full water and energy balances, originally developed by Xu Liang at the University of Washington.  VIC is a research model and in its various forms it has been applied to most of the major river basins around the world, as well as being applied [globally](links.md).  The VIC model is distributed under the [GNU GPL v2.0](http://www.gnu.org/licenses/gpl-2.0.html) license. If you make use of this model, please acknowledge the appropriate references listed on the [references](Documentation/References.md) page.

Development and maintenance of the current official version of the VIC model is conducted at the [University of Washington](http://www.washington.edu), [Department of Civil and Environmental Engineering](http://www.ce.washington.edu), under the direction of [Bart Nijssen](http://www.hydro.washington.edu/~nijssen/).  Every new application addresses new problems and conditions which the model may not currently be able to handle, and as such the model is always under development. The VIC model is an open source develment project, contributions from outside the [University of Washington Land Surface Hydrology Group](http://www.hydro.washington.edu/) are welcome.  Archived, current, beta, and development versions of the model are available via the group's [GitHub repository](https://github.com/UW-Hydro).

-----

# News
** VIC source code is now open source **

VIC is now an open-source model.  We encourage any users who have modified (or would like to modify) VIC, either to fix bugs or develop new features, to contact us and coordinate development work with us.  The VIC model source code is archived in Git and is publicly available through [GitHub](https://github.com).  More details are available on the [Model Development](Development/ModelDevelopment.md) page.  Instructions on how to use Git and GitHub for VIC development are available [here](Development/working-with-git.md).

** VIC administrator email account is no longer monitored **

The VIC administrator email account is no longer actively monitored.  Any emails to the VIC administrator will be automatically forwarded to the vic_users email list.  Therefore, if you have a question about running VIC, please just email the vic_users list.

** 0.5-degree Global VIC parameters and meteorological forcings available, April 5, 2013 **

Visit [the "Datasets" page](Datasets/Datasets.md) for more details.

----

# Release History

** VIC 4.2 Released, November 19, 2014 **

VIC release 4.2 contains the following new features relative to 4.1.2:

*   Partial vegetation cover*   Non-climatological time-varying vegetation parameters
*   Simulation of photosynthesis
*   Simulation of soil carbon storage and fluxes
*   Updates to frozen soil options
*   Removal of some obsolete options and general clean-up of the sample global parameter file*   Numerous bug fixes
For more details, see the description of the [current release](Development/CurrentVersion.md).

** VIC 4.1.2 Released, Feb. 8, 2012 **

-----
