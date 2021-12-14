/******************************************************************************
 * @section DESCRIPTION
 *
 * Function to shift columns or rows one position (used in convolution)
 *****************************************************************************/

#include <rout.h>

/******************************************************************************
 * @brief   Function to shift columns or rows one position
 *****************************************************************************/
void
cshift(double *data,
       int     nx,
       int     ny,
       int     axis,
       int     direction)
{
    int    x, y;
    double b;

    if (axis == 0 && direction == 1) {
        for (y = 0; y != ny; y++) {
            b = *(data + y);
            for (x = 0; x != nx - 1; x++) {
                *(data + y + ny * x) = *(data + y + ny * (x + 1));
            }
            *(data + y + ny * x) = b;
        }
    }

    if (axis == 0 && direction == -1) {
        for (y = 0; y != ny; y++) {
            b = *(data + y + ny * (nx - 1));
            for (x = nx - 1; x >= 0; x--) {
                *(data + y + ny * (x + 1)) = *(data + y + ny * x);
            }
            *(data + y) = b;
        }
    }

    if (axis == 1 && direction == 1) {
        for (x = 0; x < nx; x++) {
            b = *(data + x * ny);
            for (y = 0; y != ny; y++) {
                *(data + y + ny * x) = *(data + y + 1 + ny * x);
            }
            *(data + y - 1 + ny * x) = b;
        }
    }

    if (axis == 1 && direction == -1) {
        for (x = 0; x < nx; x++) {
            b = *(data + ny - 1 + ny * x);
            for (y = ny - 2; y >= 0; y--) {
                *(data + y + 1 + ny * x) = *(data + y + ny * x);
            }
            *(data + x * ny) = b;
        }
    }
}
