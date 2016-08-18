# Development Workflow
The basic workflow described here follows a workflow originally [outlined](http://nvie.com/posts/a-successful-git-branching-model/) by Vincent Driessen.  The workflow is built around the Git version control system.  A basic description of the branching strategy and release management used by our research group is presented here.  We use a central truth repository (https://github.com/UW-Hydro/VIC) that contains our main branches.

## Main Branches
Both of the main branches are published on the Github page and are controlled only by those within the admin group.  The repository is organized into several branches for different purposes.  In general, any new development will want to start with a branch from the develop branch.

**1. master**– The master branch represents the official release of the code.  This branch is updated by new releases from the develop/release branches and by hotfixes.

**2. develop** – The develop branch represents the bleeding edge of the code.  Because, new releases, hotfixes, and new features update the development, we recommend that all new development begins by branching from the develop branch.
## Feature Branches
For most developers, branching off the develop branch is the best choice when developing new features.  This allows for easy merging of these new features into the main source code and facilitates seamless rebasing of long running feature development with the current develop branch.   We merge completed feature branches into the develop branch for inclusion in the next release.   A developer may have many feature branches, which may or may not be intended for inclusion in the main source code.  Utilizing the frictionless context switching supported by Git, the development of features using this strategy provides developers with a clean and simple method of switching between features and the main branches during the development process.

## Support Branches
Often times VIC development is driven by projects that require very specific modifications to the code that would not be appropriate for inclusion in a major release.  The use of support branches allows for the continued development of the trunk while "supporting" project specific versions of the VIC code.  Instead of completely tossing these changes, we put the project's version of VIC in a support branch and continue developing.  Support branches are essentially branches that are not expected to merge back into the development branch.

## Admin Branches
Although anyone could create these branches, they are designed for the preparation or hotfix of an updated master branch and therefore should only be used by members of the admin group.

**1. release** – The release branch supports the preparation of a new release.  It includes any changes needed for public release or minor bug fixes.

**2. hotfix** -  The hotfix branch facilitates mid-release bug fixes of the master branch.  The key point of the hotfix branch is that it does not incorporate any new features from the develop branch, rather it is a branch off the master that addresses a specific issue or set of issues.  When the hotfix is applied, the development branch is updated to reflect the hotfix changes.

## Naming Conventions
* Master branch – `master`
* Develop branch – `develop`
* Feature branch – `feature/{feature_name}`
* Hotfix branch – `hotfix/{hotfix_name}`
* Release branch – `release/{release_name}`
* Support branch – `support/VIC.{base_release_number}.{feature_branch_name}`
* Release name – `VIC.{major.minor.patch}`
* Support release name – `VIC.{base_release_number}.{feature_branch_name}.{##}`

# User Permissions
Using Github to host the central or truth repository of our models allows us to easily control contributor permissions.  Currently we split permission levels into 3 levels, Owners, Model Admins, and Developers.

1. Owners have full access to all repositories and have admin rights to the organization.

2. Model Admins have full access to specific repositories.  They may push, pull, or make administrative changes to those repositories associated with their model.

3. Developers have pull only access to any of the publicly listed repositories under the UW-hydro name.  If a developer would like a feature branch merged into the main repository, a pull request must be submitted and a Model Admin may merge it in.
