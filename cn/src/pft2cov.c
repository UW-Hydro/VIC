#include <vic_def.h>
#include <vic_run.h>

double pft2cov(double pft_vals[21], int veg_class)

/****************************************************************************
                                                                           
  pft2cov:  converts UMD vegetation cover values from CLM PFT values.

****************************************************************************/
{

  double cov_val;

  switch(veg_class)
    {
    case 0: cov_val = pft_vals[2];
      break;
    case 1: cov_val = pft_vals[5];
      break;
    case 2: cov_val = pft_vals[3];
      break;
    case 3: cov_val = pft_vals[8];
      break;
    case 4: cov_val = 0.5 * pft_vals[2] + 0.5 * pft_vals[8];
      break;
    case 5: cov_val = 0.4 * pft_vals[2] + 0.4 * pft_vals[8] + \
	0.1 * pft_vals[9] + 0.1 * pft_vals[10];
      break;
    case 6: cov_val = 0.1 * pft_vals[14] + 0.3 * pft_vals[13] + \
	0.2 * pft_vals[12] + 0.2 * pft_vals[2] + 0.2 * pft_vals[8];
      break;
    case 7: cov_val = 0.3 * pft_vals[9] + 0.3 * pft_vals[10] + \
	0.2 * pft_vals[11] + 0.1 * pft_vals[2] + 0.1 * pft_vals[8];
      break;
    case 8: cov_val = 0.2 * pft_vals[9] + 0.2 * pft_vals[10] + \
	0.2 * pft_vals[11] + 0.05 * pft_vals[14] + 0.15 * pft_vals[13] + \
	0.1 * pft_vals[12] + 0.1 * pft_vals[0];
      break;
    case 9: cov_val = 0.1 * pft_vals[14] + 0.3 * pft_vals[13] + \
	0.2 * pft_vals[12] + 0.4 * pft_vals[0];
      break;
    case 10: cov_val = pft_vals[17];
      break;
    default: cov_val = pft_vals[0];
      break;
    }

  return(cov_val);

}

