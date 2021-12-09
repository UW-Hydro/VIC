/******************************************************************************
 * @section DESCRIPTION
 *
 * Allocate memory for VIC structures.
 *****************************************************************************/

#include <vic_driver_cesm.h>

/******************************************************************************
 * @brief    Allocate memory for VIC structures.
 *****************************************************************************/
void
vic_cesm_alloc(void)
{
    extern x2l_data_struct *x2l_vic;
    extern l2x_data_struct *l2x_vic;
    extern domain_struct    local_domain;

    debug("In vic_cesm_alloc");

    // allocate memory for x2l_vic structure
    x2l_vic = malloc(local_domain.ncells_active * sizeof(*x2l_vic));
    check_alloc_status(x2l_vic, "Memory allocation error.");
    // initialize x2l data
    initialize_x2l_data();

    // allocate memory for l2x_vic structure
    l2x_vic = malloc(local_domain.ncells_active * sizeof(*l2x_vic));
    check_alloc_status(l2x_vic, "Memory allocation error.");
    // initialize l2x data
    initialize_l2x_data();

    // allocate the rest of the image mode structures
    vic_alloc();
}
