# Variable Infiltration Capacity (VIC) Model

[![Build Status](https://travis-ci.org/UW-Hydro/VIC.png?branch=develop)](https://travis-ci.org/UW-Hydro/VIC) [![Join the chat at https://gitter.im/UW-Hydro/VIC](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/UW-Hydro/VIC?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

VIC Users Listserve: [![VIC Users Listserve](https://img.shields.io/badge/VIC%20Users%20Listserve-Active-blue.svg)](https://mailman.u.washington.edu/mailman/listinfo/vic_users)

Developers Gitter Room: [![Join the chat at https://gitter.im/UW-Hydro/VIC](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/UW-Hydro/VIC?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

License: [![GitHub license](https://img.shields.io/badge/license-GPLv2-blue.svg)](https://raw.githubusercontent.com/UW-Hydro/VIC/master/LICENSE.txt)

VIC Documentation: [![Documentation Status](https://readthedocs.org/projects/vic/badge/?version=latest)](http://vic.readthedocs.org/en/latest/)

Current Release DOI: [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.22307.svg)](http://dx.doi.org/10.5281/zenodo.22307)

------

This repository serves as the public source code repository of the **Variable Infiltration Capacity Model**, better known as **VIC**. The [official VIC website](http://vic.readthedocs.org) it still maintained at the [University of Washington Computational Hydrology Group](http://uw-hydro.github.io/) headed by [Bart Nijssen](http://uw-hydro.github.io/current_member/bart_nijssen/).

The Variable Infiltration Capacity (VIC) macroscale hydrological model (MHM) has been developed over the last two decades at the [University of Washington](http://uw-hydro.github.io/) and [Princeton University](http://hydrology.princeton.edu) in collaboration with a large number of other researchers around the globe. A skeletal first version of the VIC model was  introduced to the community by [Wood et al. [1992]](http://dx.doi.org/10.1029/91JD01786) and a greatly expanded version, from which current variations evolved, is described by [Liang et al. [1994]](http://dx.doi.org/10.1029/94jd00483). As compared to other MHMs, VIC’s distinguishing hydrological features are its representation of subgrid variability in soil storage capacity as a spatial probability distribution to which surface runoff is related, and its parameterization of base flow, which occurs from a lower soil moisture zone as a nonlinear recession. Movement of moisture between the soil layers is modeled as gravity drainage, with the unsaturated hydraulic conductivity a function of the degree of saturation of the soil. Spatial variability in soil properties, at scales smaller than the grid scale, is represented statistically, without assigning infiltration parameters to specific subgrid locations. Over time, many additional features and representations of physical processes have been added to the model. VIC has been used in a large number of regional and continental scale (even global) hydrological studies. A selection of VIC applications by the University of Washington Computational Hydrology Group can be found on the [VIC references page](http://vic.readthedocs.org/en/latest/about/references.md) and the [publications page](http://uw-hydro.github.io/publications/).

Development and maintenance of the current official version of the VIC model is conducted at the University of Washington, Department of Civil and Environmental Engineering, under the direction of Dennis P. Lettenmaier. By placing the original source code archive on GitHub, we hope to encourage a more collaborative development environment. A tutorial on how to use the VIC git repository and how to contribute your changes to the VIC model can be found on the [github VIC wiki](https://github.com/UW-Hydro/VIC/wiki). The most stable version of the model is in the `master` branch, while beta versions of releases under development can be obtained from the `development` branch of this repository.

VIC is a research model developed by graduate students, post-docs and research scientists over a long period of time (since the early 1990s). Every new VIC application addresses new problems and conditions which the model may not currently be able to handle. As a result, the model is always under development. Because of the incremental nature of this development, not all sections of the code are equally mature and not every combination of model options has been exhaustively tested or is guaranteed to work. While you are more than welcome to use VIC in your own research endeavors, ***the model code comes with no guarantees, expressed or implied, as to suitability, completeness, accuracy, and whatever other claim you would like to make***. In addition, the model has no graphical user interface, nor does it include a large set of analysis tools, since most people want to use their own set of tools.

While we would like to hear about your particular application (especially a copy of any published paper), we cannot give you individual support in setting up and running the model. The [VIC website](http://vic.readthedocs.org) includes reasonably complete instructions on how to run the model, as well as the opportunity to sign up for the VIC Users Email List.

If you make use of this model, please acknowledge the appropriate references listed on the [references page](http://vic.readthedocs.org/en/latest/Documentation/References/).

If you want to contact us directly about the VIC model, you can send
e-mail to `vic_users at u.washington.edu`. Please include "VIC Model" in the subject line.
