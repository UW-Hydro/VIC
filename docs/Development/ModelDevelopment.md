# Open Source

We encourage any users who have modified (or would like to modify) VIC, either to fix bugs or develop new features, to contact us and coordinate development work with us. The VIC model source code is archived in Git and is publicly available through [GitHub](https://github.com). To access the source code, visit [GitHub](https://github.com), create an account, and visit [github.com/UW-Hydro/VIC](https://github.com/UW-Hydro/VIC). The routing model will soon be added to the archive as well.

VIC is an [open source](open-source-philosophy.md) development model and is released under the terms of the [GNU General Public License v2 (GPL).](http://www.gnu.org/licenses/old-licenses/gpl-2.0.html)

Instructions for using Git and GitHub to access the VIC code and contribute changes are here:

*   **[Working with VIC in the git version control system](working-with-git.md)**

# Model Development

Information about the VIC codebase:

- [Current Versions of VIC and Routing Model](ReleaseNotes.md)
- [Known issues](https://github.com/UW-Hydro/VIC/issues)
- [Contributing](Contributing.md)

All development activity is coordinated via the [VIC github page](https://github.com/UW-Hydro/VIC), where you can also find all archived, current, beta, and development versions of the model.

# Release Naming Convention

Beginning with the release of VIC 4.2, a release naming convention of VIC.MAJOR.MINOR.PATCH has been formally adopted. This naming system was introduced by [Semantic Versioning](http://semver.org/spec/v2.0.0.html) and defined as:

1.  MAJOR version: when changes are included that are not backward compatible,
2.  MINOR version: when features and new functionality are added in a backwards-compatible manner, and
3.  PATCH version: when backwards-compatible bug fixes are made.

# Release Roadmap

Regardless of the release type (major, minor, patch), the VIC model release process follows this roadmap.

1. Identification of features and fixes to be included in the release.  This is coordinated through the [VIC github page](https://github.com/UW-Hydro/VIC).
1. Features and fixes are added to the appropriate branch (`develop` or `hotfix`) via Github pull requests.
1. The model is tested using the [VIC test framework](Testing.md).  Each of the following tests should be run and verified for any release.
    - Unit tests
    - System tests
    - Example tests
    - Science tests
    - Release tests
1. Once all tests have been run and the release is ready, a few minor changes need to be made to the source code prior to a release:
    - Update the version string for classic and image drivers (`vic_version.h`)
    - Update the version string for the python driver (`setup.py`)
    - Update the documentation.  At a minimum this means updating [ReleaseNotes.md]('ReleaseNotes.md').
1.  Make the release!
    - Merge the `release` branch into the `master` branch.
    - Create the tag. This is best done on the command line as `git tag -f -a VIC.5.1.2 -m 'release tag of VIC 5.1.2, bug fix update for this example`.
    - Publish the release on [Github](https://github.com/UW-Hydro/VIC/releases).
    - Update the VIC documentation and _Readme_ files with the releases Zenodo DOI (e.g. [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.35303.svg)](http://dx.doi.org/10.5281/zenodo.35303)).
    - Update the `develop` branch with the changes from the `master` branch.
    - Reset the version strings in the `develop` branch.
