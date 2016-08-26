# Contributing

VIC is Open Source software. This means that the code is made available for free, but also that development, maintenance and support need to be community efforts. Our rationale for moving the VIC model development to the open source community is that we want:

- to encourage other researchers and developers to contribute to VIC development, and
- to facilitate transparent development and use of the model.

## Support

There is no official support for the VIC model, other than the VIC documentation, the VIC source code archive and the description of the model in the literature. Any additional support relies on volunteer efforts by the  VIC development community, which means that no one is being paid to provide VIC support. The following resources are available:

- [VIC Source code repository](https://github.com/UW-Hydro/VIC) : Source code distribution, coordination of model development, bug fixes, and releases.
- [VIC documentation](http://vic.readthedocs.org) at readthedocs.org.

We expect that the user comes prepared with some understanding of the model and scientific computing. As such, these items are specifically not supported by the VIC development community:

- Building and running the VIC model on platforms other than LINUX, UNIX, and OSX.
- Using LINUX, UNIX, or OSX operating systems.
- Development of project specific features.
- Configuring individual model applications.

## Getting In-touch

The VIC user and developer community is quite large and is spread all around the world.  Depending on the question or issue your facing, there are a number of avenues to get in-touch:

VIC Users Listserve: [![VIC Users Listserv](https://img.shields.io/badge/VIC%20Users%20Listserve-Active-blue.svg)](https://mailman.u.washington.edu/mailman/listinfo/vic_users)

Developers Gitter Room: [![Join the chat at https://gitter.im/UW-Hydro/VIC](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/UW-Hydro/VIC?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

VIC Github Issues [![GitHub issues](https://img.shields.io/github/issues/UW-Hydro/VIC.svg)](https://github.com/UW-Hydro/VIC/issues)

## Submitting Issues
#### Submitting Bug Reports

If you think you have found a bug in VIC, please file an issue on Github [here](https://github.com/UW-Hydro/VIC/issues). Please include the following information in your bug report:

- Version of VIC that you are using (e.g. VIC.4.2.b)
- Name and version of the C compiler you are using
- Operating system
- A description of relevant model settings
- A summary of the bug or error message you are getting

If you can provide more information that is great. If you know how to run the model in a debugger, you may be able to pinpoint where the problem occurs.

#### Proposing New Features

The VIC model is under active development.  If you would like to propose a new feature, driver, or extension to the VIC model, please raise a Github issue [here](https://github.com/UW-Hydro/VIC/issues). Also, because VIC is an open source model with no official support for development, be prepared to contribute to the implementation of your feature request. Also note that features that are only of interest to you are unlikely to be implemented in the main source code repo (although you are of course free to modify the code in any way you see fit).

## Contributing to VIC
#### Git Workflow
We have developed some documentation to help you get started if you are new to Git but want to contribute to VIC:

- [Working with Git](http://vic.readthedocs.org/en/latest/Development/working-with-git/)
- [Git Workflow](http://vic.readthedocs.org/en/latest/Development/git-workflow/)

#### Style Guide
To be added soon.  For now, please remember to comment your code well and [use `uncrustify` to format your code properly](https://github.com/UW-Hydro/VIC/tree/master/tools/code_format).
