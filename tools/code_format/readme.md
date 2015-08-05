*VIC Code Format Conventions*

**Variable Names**

 1. Use self-describing variable names, even if they have a large number of characters. For example, `canopy_evap` is preferable to `ce`

 1. Local Variables: `lowercase_separated_by_underscores`

 1. Global Variables:  `CAPITALIZED_WITH_UNDERSCORES`

 1. Define every single variable, including units, with each variable on a separate line. That is
    ```C
    double  canopy_latent;  /**< latent heat flux from the canopy (W/m^2) */
    ```

1. Be explicit about integer types.  For example use `unsigned int` rather than simply `unsigned`.

**Commenting**

 1. Include a comment block at the start of each function, enum, and structure, providing a brief description, including references where applicable.

 1. Do not commit code with large sections of code commented out. The version control system should be used for tracking different versions. There is no need to do it manually by commenting out code, which is much more likely to lead to confusion.

 1. Comment often. Strive for a comment on every line of code. It takes much less time to "*comment as you go*" than try and comment afterwards.

 1. When there is potential for confusion, define units for the LHS of the assignment.

**Indenting**

Use four spaces to indent

 * subroutines within a module
 * if statements (if - else)
 * for and while loops

**Functions**

 1. Organize functions into modules (even if the module contains just one subroutine). 

    ```C
    /**************************************************************************
     * @brief    Determines from the air temperature what fraction of incoming
     *           precipitation is frozen and unfrozen (snow and rain).
     *************************************************************************/
    double
    calc_rainonly(double air_temp,
                  double prec,
                  double MAX_SNOW_TEMP,
                  double MIN_RAIN_TEMP)
    {
        double rainonly;

        rainonly = 0.;
        if (MAX_SNOW_TEMP <= MIN_RAIN_TEMP) {
            log_err("MAX_SNOW_TEMP must be greater then MIN_RAIN_TEMP");
        }
        if (air_temp < MAX_SNOW_TEMP && air_temp > MIN_RAIN_TEMP) {
            rainonly = (air_temp - MIN_RAIN_TEMP) /
                       (MAX_SNOW_TEMP - MIN_RAIN_TEMP) * prec;
        }
        else if (air_temp >= MAX_SNOW_TEMP) {
            rainonly = prec;
        }

        return(rainonly);
    }

    ```

 **Licensing**

 1. VIC is distributed under the GPL-V2 license. This means that all code modifications must be shared with the community. All code files should include a header that contains the following (this should be copied verbatim other than the text in `{}`).

        The Variable Infiltration Capacity (VIC) macroscale hydrological model
        Copyright (C) {} The Land Surface Hydrology Group, Department of Civil
        and Environmental Engineering, University of Washington.

        The VIC model is free software; you can redistribute it and/or
        modify it under the terms of the GNU General Public License
        as published by the Free Software Foundation; either version 2
        of the License, or (at your option) any later version.

        This program is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        GNU General Public License for more details.

        You should have received a copy of the GNU General Public License along with
        this program; if not, write to the Free Software Foundation, Inc.,
        51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

**Source Code Beautifier**

1.  VIC source code is formated using the [Uncrustify](http://uncrustify.sourceforge.net/) source code beautifier.  Uncrustify should be run prior to committing changes to the Git repository.

2.  To format a VIC source code file run:

    ```bash
    uncrustify -c uncrustify_VIC_c.cfg [files ...]
    ```

    or with replacement and without backup files:

    ```bash
    uncrustify -c uncrustify_VIC_c.cfg --replace --no-backup [files ...]
    ```

    finally, to format all files in the VIC repository:
    ```bash
    ./run_uncrustify.bash
    ```

**Miscellaneous**

1.  When typecasting variables, add a single space after typecast before the variable name, e.g. `new_var = (int) my_var;`

1.  Use integer and float printf format specifiers for all integer and floating point types, regarless of signedness or precision.
