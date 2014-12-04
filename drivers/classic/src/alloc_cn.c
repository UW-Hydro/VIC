#include <stdio.h>
#include <stdlib.h>
#include <vic_def.h>
#define MAX_PFT 21

static char vcid[] = "$Id$";

/****************************************************************************/
/*			       alloc_cn()                                */
/****************************************************************************/
void alloc_cn(int Nbands, int Nnode, cn_data_struct **cn)
/*******************************************************************
  alloc_cn       Created by Michael Brunke

  Allocates data structure containing CN quantities in CN data
  structure.

  Modifications:

*******************************************************************/
{

  int iband, iveg, k;

  *cn = (cn_data_struct *) calloc(Nbands, sizeof(cn_data_struct));
  /*  if (*cn == NULL)
      { */
      for(iband = 0; iband < Nbands; iband++)
	{
	  for(iveg = 0; iveg < MAX_PFT; iveg++)
	    {
	      (*cn)[iband].LAI[iveg] = 0.0;
	      (*cn)[iband].dormant_flag[iveg] = 0.0;
	      (*cn)[iband].days_active[iveg] = 0.0;
	      (*cn)[iband].onset_flag[iveg] = 0.0;
	      (*cn)[iband].onset_counter[iveg] = 0.0;
	      (*cn)[iband].onset_gddflag[iveg] = 0.0;
	      (*cn)[iband].onset_fdd[iveg] = 0.0;
	      (*cn)[iband].onset_gdd[iveg] = 0.0;
	      (*cn)[iband].onset_swi[iveg] = 0.0;
	      (*cn)[iband].offset_flag[iveg] = 0.0;
	      (*cn)[iband].offset_counter[iveg] = 0.0;
	      (*cn)[iband].offset_fdd[iveg] = 0.0;
	      (*cn)[iband].offset_swi[iveg] = 0.0;
	      (*cn)[iband].lgsf[iveg] = 0.0;
	      (*cn)[iband].bglfr[iveg] = 0.0;
	      (*cn)[iband].bgtr[iveg] = 0.0;
	      (*cn)[iband].dayl[iveg] = 0.0;
	      (*cn)[iband].prev_dayl[iveg] = 0.0;
	      (*cn)[iband].annavg_t2m[iveg] = 0.0;
	      (*cn)[iband].tempavg_t2m[iveg] = 0.0;
	      (*cn)[iband].gpp2[iveg] = 0.0;
	      (*cn)[iband].availc[iveg] = 0.0;
	      (*cn)[iband].xsmrpool_recover[iveg] = 0.0;
	      (*cn)[iband].alloc_pnow[iveg] = 0.0;
	      (*cn)[iband].c_allometry[iveg] = 0.0;
	      (*cn)[iband].n_allometry[iveg] = 0.0;
	      (*cn)[iband].tempsum_potential_gpp[iveg] = 0.0;
	      (*cn)[iband].annsum_potential_gpp[iveg] = 0.0;
	      (*cn)[iband].tempmax_retransn[iveg] = 0.0;
	      (*cn)[iband].annmax_retransn[iveg] = 0.0;
	      (*cn)[iband].avail_retransn[iveg] = 0.0;
	      (*cn)[iband].plant_nalloc[iveg] = 0.0;
	      (*cn)[iband].plant_calloc[iveg] = 0.0;
	      (*cn)[iband].excess_cflux[iveg] = 0.0;
	      (*cn)[iband].downreg[iveg] = 0.0;
	      (*cn)[iband].prev_leafc_to_litter[iveg] = 0.0;
	      (*cn)[iband].prev_frootc_to_litter[iveg] = 0.0;
	      (*cn)[iband].tempsum_npp[iveg] = 0.0;
	      (*cn)[iband].annsum_npp[iveg] = 0.0;
	      (*cn)[iband].gpp[iveg] = 0.0;
	      (*cn)[iband].npp[iveg] = 0.0;
	      (*cn)[iband].ar[iveg] = 0.0;
	      (*cn)[iband].leafc[iveg] = 0.0;
	      (*cn)[iband].leafc_storage[iveg] = 0.0;
	      (*cn)[iband].leafc_xfer[iveg] = 0.0;
	      (*cn)[iband].frootc[iveg] = 0.0;
	      (*cn)[iband].frootc_storage[iveg] = 0.0;
	      (*cn)[iband].frootc_xfer[iveg] = 0.0;
	      (*cn)[iband].livestemc[iveg] = 0.0;
	      (*cn)[iband].livestemc_storage[iveg] = 0.0;
	      (*cn)[iband].livestemc_xfer[iveg] = 0.0;
	      (*cn)[iband].deadstemc[iveg] = 0.0;
	      (*cn)[iband].deadstemc_storage[iveg] = 0.0;
	      (*cn)[iband].deadstemc_xfer[iveg] = 0.0;
	      (*cn)[iband].livecrootc[iveg] = 0.0;
	      (*cn)[iband].livecrootc_storage[iveg] = 0.0;
	      (*cn)[iband].livecrootc_xfer[iveg] = 0.0;
	      (*cn)[iband].deadcrootc[iveg] = 0.0;
	      (*cn)[iband].deadcrootc_storage[iveg] = 0.0;
	      (*cn)[iband].deadcrootc_xfer[iveg] = 0.0;
	      (*cn)[iband].gresp_storage[iveg] = 0.0;
	      (*cn)[iband].gresp_xfer[iveg] = 0.0;
	      (*cn)[iband].cpool[iveg] = 0.0;
	      (*cn)[iband].xsmrpool[iveg] = 0.0;
	      (*cn)[iband].pft_ctrunc[iveg] = 0.0;
	      (*cn)[iband].totvegc[iveg] = 0.0;
	      (*cn)[iband].woodc[iveg] = 0.0;
	      (*cn)[iband].leafn[iveg] = 0.0;
	      (*cn)[iband].leafn_storage[iveg] = 0.0;
	      (*cn)[iband].leafn_xfer[iveg] = 0.0;
	      (*cn)[iband].frootn[iveg] = 0.0;
	      (*cn)[iband].frootn_storage[iveg] = 0.0;
	      (*cn)[iband].frootn_xfer[iveg] = 0.0;
	      (*cn)[iband].livestemn[iveg] = 0.0;
	      (*cn)[iband].livestemn_storage[iveg] = 0.0;
	      (*cn)[iband].livestemn_xfer[iveg] = 0.0;
	      (*cn)[iband].deadstemn[iveg] = 0.0;
	      (*cn)[iband].deadstemn_storage[iveg] = 0.0;
	      (*cn)[iband].deadstemn_xfer[iveg] = 0.0;
	      (*cn)[iband].livecrootn[iveg] = 0.0;
	      (*cn)[iband].livecrootn_storage[iveg] = 0.0;
	      (*cn)[iband].livecrootn_xfer[iveg] = 0.0;
	      (*cn)[iband].deadcrootn[iveg] = 0.0;
	      (*cn)[iband].deadcrootn_storage[iveg] = 0.0;
	      (*cn)[iband].deadcrootn_xfer[iveg] = 0.0;
	      (*cn)[iband].retransn[iveg] = 0.0;
	      (*cn)[iband].npool[iveg] = 0.0;
	      (*cn)[iband].pft_ntrunc[iveg] = 0.0;
	    }

	  (*cn)[iband].decl = 0.0;
	  (*cn)[iband].fpi = 0.0;
	  (*cn)[iband].fpg = 0.0;
	  (*cn)[iband].annsum_counter = 0.0;
	  (*cn)[iband].cannsum_npp = 0.0;
	  (*cn)[iband].cannavg_t2m = 0.0;
	  for(k = 0; k < Nnode; k++)
	    (*cn)[iband].watfc[k] = 0.0;
	  (*cn)[iband].me = 0.0;
	  (*cn)[iband].fire_prob = 0.0;
	  (*cn)[iband].mean_fire_prob = 0.0;
	  (*cn)[iband].fireseasonl = 0.0;
	  (*cn)[iband].ann_farea_burned = 0.0;
	  (*cn)[iband].cwdc = 0.0;
	  (*cn)[iband].litr1c = 0.0;
	  (*cn)[iband].litr2c = 0.0;
	  (*cn)[iband].litr3c = 0.0;
	  (*cn)[iband].soil1c = 0.0;
	  (*cn)[iband].soil2c = 0.0;
	  (*cn)[iband].soil3c = 0.0;
	  (*cn)[iband].soil4c = 0.0;
	  (*cn)[iband].seedc = 0.0;
	  (*cn)[iband].col_ctrunc = 0.0;
	  (*cn)[iband].totlitc = 0.0;
          (*cn)[iband].totsomc = 0.0;
	  (*cn)[iband].totcolc = 0.0;
	  (*cn)[iband].prod10c = 0.0;
	  (*cn)[iband].prod100c = 0.0;
	  (*cn)[iband].cwdn = 0.0;
	  (*cn)[iband].litr1n = 0.0;
	  (*cn)[iband].litr2n = 0.0;
	  (*cn)[iband].litr3n = 0.0;
	  (*cn)[iband].soil1n = 0.0;
	  (*cn)[iband].soil2n = 0.0;
	  (*cn)[iband].soil3n = 0.0;
	  (*cn)[iband].soil4n = 0.0;
	  (*cn)[iband].sminn = 0.0;
	  (*cn)[iband].seedn = 0.0;
	  (*cn)[iband].col_ntrunc = 0.0;
	  (*cn)[iband].totcoln = 0.0;
	  (*cn)[iband].prod10n = 0.0;
	  (*cn)[iband].prod100n = 0.0;
	  (*cn)[iband].hr = 0.0;
	  (*cn)[iband].nee = 0.0;
	  (*cn)[iband].nep = 0.0;

	  }

      /*    } */

}

/****************************************************************************/
/*	      		  free_cn()                                         */
/****************************************************************************/
void free_cn(cn_data_struct **cn)
/***************************************************************************
  Modifications:
***************************************************************************/
{

  if (*cn == NULL)
    return;

  free(*cn);
}
