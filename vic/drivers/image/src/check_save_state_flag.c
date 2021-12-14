/******************************************************************************
 * @section DESCRIPTION
 *
 * Function to check whether model state should be saved for the current
 * time step
 *****************************************************************************/

#include <vic_driver_image.h>

/******************************************************************************
 * @brief   Function to check whether model state should be saved for the
 *          current time step
 *****************************************************************************/
bool
check_save_state_flag(size_t      current,
                      dmy_struct *dmy_offset)
{
    extern global_param_struct global_param;
    extern dmy_struct         *dmy;

    double                     offset;
    double                     time_num;

    // Advance dmy by one timestep because dmy is the "timestep-beginning"
    // timestamp, but we want to check whether the end of the current
    // time step is the user-specified output state time
    offset = global_param.dt / (double) SEC_PER_DAY;
    time_num = date2num(global_param.time_origin_num, &dmy[current], 0,
                        global_param.calendar, TIME_UNITS_DAYS);
    time_num += offset;
    num2date(global_param.time_origin_num, time_num, 0,
             global_param.calendar, TIME_UNITS_DAYS,
             dmy_offset);

    // Check if the end of the current time step is equal to the state output
    // timestep specified by user
    if (dmy_offset->year == global_param.stateyear &&
        dmy_offset->month == global_param.statemonth &&
        dmy_offset->day == global_param.stateday &&
        dmy_offset->dayseconds == global_param.statesec) {
        return true;
    }
    else {
        return false;
    }
}
