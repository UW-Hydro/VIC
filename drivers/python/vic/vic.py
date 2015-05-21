#!/usr/bin/env python
"""
  @section DESCRIPTION

  Front end for VIC Python driver

  @section LICENSE

  The Variable Infiltration Capacity (VIC) macroscale hydrological model
  Copyright (C) 2014 The Land Surface Hydrology Group, Department of Civil
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
"""

import argparse
import vic


# -------------------------------------------------------------------- #
def main():
    """
    Get the script and path to the config_file
    """
    # Parse arguments
    parser = argparse.ArgumentParser(description='Python interface for the '
                                                 'Variable Infiltration '
                                                 'Capacity Model (VIC)')
    parser.add_argument("global_config_file", type=str,
                        help="Input global configuration file",
                        default=None)
    parser.add_argument("--version", action='store_true',
                        help="Return RVIC version string")

    args = parser.parse_args()

    if args.version:
        print(vic.version.short_version)
        return

    if args.config_file is not None:
        vic.init()

        vic.run()

        vic.final()
    return
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
if __name__ == "__main__":
    main()
# -------------------------------------------------------------------- #
