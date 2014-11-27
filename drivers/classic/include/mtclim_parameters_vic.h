/******************************************************************************
 * @section DESCRIPTION
 *
 * MTCLIM Parameters header file.
 *
 * @section LICENSE
 *
 * The Variable Infiltration Capacity (VIC) macroscale hydrological model
 * Copyright (C) 2014 The Land Surface Hydrology Group, Department of Civil
 * and Environmental Engineering, University of Washington.
 *
 * The VIC model is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 *****************************************************************************/

#ifndef MTCLIM_PARAMETERS_VIC_H
#define MTCLIM_PARAMETERS_VIC_H

/* parameters for the Tair algorithm */
#define TDAYCOEF     0.45  /* (dim) daylight air temperature coefficient (dim) */

/* parameters for the snowpack algorithm */
#define SNOW_TCRIT   -6.0  /* (deg C) critical temperature for snowmelt */
#define SNOW_TRATE  0.042  /* (cm/degC/day) snowmelt rate */

/* parameters for the radiation algorithm */
#define TBASE       0.870  /* (dim) max inst. trans., 0m, nadir, dry atm */
#define ABASE     -6.1e-5  /* (1/Pa) vapor pressure effect on transmittance */
#define C             1.5  /* (dim) radiation parameter */
#define B0          0.031  /* (dim) radiation parameter */
#define B1          0.201  /* (dim) radiation parameter */
#define B2          0.185  /* (dim) radiation parameter */
#define RAIN_SCALAR  0.75  /* (dim) correction to trans. for rain day */
#define DIF_ALB       0.6  /* (dim) diffuse albedo for horizon correction */
#define SC_INT       1.32  /* (MJ/m2/day) snow correction intercept */
#define SC_SLOPE    0.096  /* (MJ/m2/day/cm) snow correction slope */

#endif

