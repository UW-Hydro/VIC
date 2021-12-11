/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine redistributes soil properties based on the thermal solutions
 * found for the current time step.
 *****************************************************************************/

#include <vic_run.h>

/******************************************************************************
 * @brief    This subroutine allocates memory for a 2-dimensional double array
 *****************************************************************************/
void
malloc_2d_double(size_t   *shape,
                 double ***array)
{
    size_t i;

    *array = malloc(shape[0] * sizeof(*(*array)));
    check_alloc_status(*array, "Memory allocation error in.");
    for (i = 0; i < shape[0]; i++) {
        (*array)[i] = malloc(shape[1] * sizeof(*((*array)[i])));
        check_alloc_status((*array)[i], "Memory allocation error in.");
    }
}

/******************************************************************************
 * @brief    This subroutine allocates memory for a 3-dimensional double array
 *****************************************************************************/
void
malloc_3d_double(size_t    *shape,
                 double ****array)
{
    size_t i;
    size_t j;

    *array = malloc(shape[0] * sizeof(*(*array)));
    check_alloc_status(*array, "Memory allocation error.");

    for (i = 0; i < shape[0]; i++) {
        (*array)[i] = malloc(shape[1] * sizeof(*((*array)[i])));
        check_alloc_status((*array)[i], "Memory allocation error.");
        for (j = 0; j < shape[1]; j++) {
            (*array)[i][j] = malloc(shape[2] * sizeof(*((*array)[i][j])));
            check_alloc_status((*array)[i][j], "Memory allocation error.");
        }
    }
}

/******************************************************************************
 * @brief    This subroutine frees memory for a 2-dimensional double array
 *****************************************************************************/
void
free_2d_double(size_t  *shape,
               double **array)
{
    size_t i;

    for (i = 0; i < shape[0]; i++) {
        free(array[i]);
    }
    free(array);
}

/******************************************************************************
 * @brief    This subroutine frees memory for a 3-dimensional double array
 *****************************************************************************/
void
free_3d_double(size_t   *shape,
               double ***array)
{
    size_t i;
    size_t j;

    for (i = 0; i < shape[0]; i++) {
        for (j = 0; j < shape[1]; j++) {
            free(array[i][j]);
        }
        free(array[i]);
    }
    free(array);
}
