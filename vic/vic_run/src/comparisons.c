/******************************************************************************
 * @section DESCRIPTION
 *
 * floating point comparison utilities
 *****************************************************************************/

#include <vic_def.h>

/******************************************************************************
 * @brief    returns false if two floats are not equal up to desired tolerance
 *****************************************************************************/
bool
assert_close_float(float x,
                   float y,
                   float rtol,
                   float abs_tol)
{
    if (fabs(x - y) <= abs_tol + rtol * fabs(y)) {
        return true;
    }
    return false;
}

/******************************************************************************
 * @brief    returns false if two doubles are not equal up to desired tolerance
 *****************************************************************************/
bool
assert_close_double(double x,
                    double y,
                    double rtol,
                    double abs_tol)
{
    if (fabs(x - y) <= abs_tol + rtol * fabs(y)) {
        return true;
    }
    return false;
}
