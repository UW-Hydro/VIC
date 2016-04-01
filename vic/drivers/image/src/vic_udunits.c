/******************************************************************************
 * @section DESCRIPTION
 *
 * Save model state.
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

#include <vic_driver_shared_image.h>

/******************************************************************************
 * @brief    udunits2 support.
 *****************************************************************************/
cv_converter*
udunits_conversion(const char*const from,
                   const char*const to)
{
    ut_system   *unitSystem;
    ut_unit     *ut_from;
    ut_unit     *ut_to;
    cv_converter*converter;

    // Read udunits unitSystem
    ut_set_error_message_handler(ut_ignore);
    unitSystem = ut_read_xml(NULL);

    ut_from = ut_parse(unitSystem, from, UT_ASCII);
    ut_to = ut_parse(unitSystem, to, UT_ASCII);

    // Check
    if (!ut_are_convertible(ut_from, ut_to)) {
        log_warn(
            "Cannot convert, or does not recognize the unit %s, assuming %s",
            from, to);
        ut_from = ut_to;
    }

    converter = ut_get_converter(ut_from, ut_to);
    return(converter);
}
