# Using Xcode as a development environment for VIC

## Project setup

The easiest way to use Xcode is to set up a project that uses an external build system (the Makefile provided with VIC).

The following works for the develop branch in VIC 4.1.2 and with Xcode 5.1.

* Open Xcode

* Select `File --> New --> Project`

* Choose `OS X | Other | External Build System`

* Fill out the details on the page. Make sure that the `Build Tool` field reflects the path for the `make` executable. The default `/usr/bin/make` works fine, but if you want to use another variety of make you can specify the details here. For `Product Name` you can pick anything you like, for example `vic_xcode`, which is what will be used in the rest of this readme.

* Specify the path for the Xcode project. The easiest would be outside your git tree, unless you want to store the project file in git.

You now have an Xcode project, but need to add the files and point to the VIC makefile.

* Select `File --> Add Files to "vic_xcode"`

* Select the VIC directories that contain the source code (note you do not want to add all the `*.o` files, so do a `make clean` first.

    * Do not select `Copy items into destination group's folder (if needed)`.

    * You can either use `Create groups for any added folders` or `Create folder references for any added folders` (check the Xcode help for the difference).

    * If you have multiple targets, make sure to select the `vic_xcode` target in the `Add to targets` list.

    * Select `Add`.

* Select `View --> Navigators --> Show Project Navigator`

* Click on the project name in the left hand side bar or select the target if you have multiple targets.

* In the window that shows the `External Build Tool Configuration` add the directory where `make` should be run (the directory with the top level makefile).

## Build

* You can now just run build, for example by clicking the big play button in the upper left corner.

* Because you are using an external build system, you can control the build settings as you do on the command-line, i.e. by editing the makefile and re-building. For example, Xcode by default will use the Apple LLVM compiler. If you want to use another one, the easiest way is just to specify the full path in the Makefile (`CC= `).

## Notes

* If you want to add the Xcode project to git, beware that there is user data in the project directory. You may want to add most of that to your `.gitignore` file.



