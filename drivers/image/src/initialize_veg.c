#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_image.h>

void
initialize_veg(veg_var_struct **veg_var,
               size_t           Nveg)

/**********************************************************************
   initialize_veg		Dag Lohmann	 January 1996

   This routine initailizes the vegetation variable array.
**********************************************************************/
{
    extern option_struct options;

    size_t               i, j, k;

    for (i = 0; i < Nveg; i++) {
        for (j = 0; j < options.SNOW_BAND; j++) {
            veg_var[i][j].Wdew = 0.0;
            veg_var[i][j].throughfall = 0.0;
            veg_var[i][j].LAI = 0.0;
            veg_var[i][j].Wdmax = 0.0;
            veg_var[i][j].vegcover = 0.0;
            veg_var[i][j].albedo = 0.0;
        }
        if (options.CARBON) {
            for (j = 0; j < options.SNOW_BAND; j++) {
                veg_var[i][j].aPAR = 0.0;
                veg_var[i][j].Ci = 0.0;
                veg_var[i][j].rc = 0.0;
                veg_var[i][j].GPP = 0.0;
                veg_var[i][j].Rphoto = 0.0;
                veg_var[i][j].Rdark = 0.0;
                veg_var[i][j].Rmaint = 0.0;
                veg_var[i][j].Rgrowth = 0.0;
                veg_var[i][j].Raut = 0.0;
                veg_var[i][j].NPP = 0.0;
                veg_var[i][j].AnnualNPP = 0.0;
                veg_var[i][j].AnnualNPPPrev = 0.0;
                for (k = 0; k < options.Ncanopy; k++) {
                    veg_var[i][j].NscaleFactor[k] = 0.0;
                    veg_var[i][j].aPARLayer[k] = 0.0;
                    veg_var[i][j].CiLayer[k] = 0.0;
                    veg_var[i][j].rsLayer[k] = 0.0;
                }
            }
        }
    }
}
