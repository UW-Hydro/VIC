/******************************************************************************
 * @section DESCRIPTION
 *
 * Routines to calculate and store model runtime timing.
 *****************************************************************************/

#include <vic_driver_shared_all.h>

/******************************************************************************
 * @brief    Get wall time
 *****************************************************************************/
double
get_wall_time()
{
    struct timeval time;
    if (gettimeofday(&time, NULL)) {
        log_err("Unable to get time of day")
    }
    return (double) time.tv_sec + (double) time.tv_usec * 0.000001;
}

/******************************************************************************
 * @brief    Get CPU time
 *****************************************************************************/
double
get_cpu_time()
{
    return (double) clock() / CLOCKS_PER_SEC;
}

/******************************************************************************
 * @brief    Initialize timer values
 *****************************************************************************/
void
timer_init(timer_struct *t)
{
    t->start_wall = 0;
    t->start_cpu = 0;

    t->delta_wall = 0;
    t->delta_cpu = 0;
}

/******************************************************************************
 * @brief    Start timer
 *****************************************************************************/
void
timer_start(timer_struct *t)
{
    timer_init(t);

    t->start_wall = get_wall_time();
    t->start_cpu = get_cpu_time();
}

/******************************************************************************
 * @brief    Stop timer
 *****************************************************************************/
void
timer_stop(timer_struct *t)
{
    t->stop_wall = get_wall_time();
    t->stop_cpu = get_cpu_time();

    t->delta_wall += t->stop_wall - t->start_wall;
    t->delta_cpu += t->stop_cpu - t->start_cpu;
}

/******************************************************************************
 * @brief    Continue timer without resetting counters
 *****************************************************************************/
void
timer_continue(timer_struct *t)
{
    t->start_wall = get_wall_time();
    t->start_cpu = get_cpu_time();
}
