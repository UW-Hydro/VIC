#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

#if OUTPUT_FORCE
void write_forcing_file(atmos_data_struct *atmos,
			int                nrecs,
			int                dt,
			outfiles_struct   *outfiles) 
/**********************************************************************
  write_forcing_file          Keith Cherkauer           July 19, 2000

  This routine writes the complete forcing data files for use in 
  future simulations.

  xx-xx-01 Modified to output pressures, which are handled internally
           in kPa, as Pa for backward compatability.              KAC
  2005-Mar-24 Added support for ALMA variables.			TJB

**********************************************************************/
{
  extern option_struct options;

  int                 rec, j;
  short int          *tmp_siptr;
  unsigned short int *tmp_usiptr;
  float              *tmp_fptr;
  atmos_data_alma_struct *atmos_alma;

  atmos_alma = (atmos_data_alma_struct*)calloc(1,sizeof(atmos_data_alma_struct));

  // Create forcing files
  for ( rec = 0; rec < nrecs; rec++ ) {
    for ( j = 0; j < NF; j++ ) {
      if (options.ALMA_OUTPUT) {
        // ALMA Standard
        conv_force_vic2alma(atmos, rec, j, dt, atmos_alma); 
        if (options.BINARY_OUTPUT) {
          /* Write a binary forcing file */
	  tmp_usiptr = (unsigned short int *)calloc(1,sizeof(unsigned short int));
	  /* SWdown * 50 */
	  tmp_usiptr[0] = (unsigned short int)(atmos_alma->SWdown*50.);
	  /* LWdown * 80 */
	  tmp_usiptr[1] = (unsigned short int)(atmos_alma->LWdown*80.);
	  /* Tair * 80 */
	  tmp_usiptr[2] = (unsigned short int)(atmos_alma->Tair*80.);
	  /* Qair * 10000000 */
	  tmp_usiptr[3] = (unsigned short int)(atmos_alma->Qair*10000000.);
	  /* Psurf * 0.01 */
	  tmp_usiptr[4] = (unsigned short int)(atmos_alma->Psurf*0.01);
	  /* Rainf * 400000000 */
	  tmp_usiptr[5] = (unsigned short int)(atmos_alma->Rainf*400000000.);
	  /* Snowf * 400000000 */
	  tmp_usiptr[6] = (unsigned short int)(atmos_alma->Snowf*400000000.);
	  /* Wind * 100 */
	  tmp_usiptr[7] = (unsigned short int)(atmos_alma->Wind*100.);
	  fwrite(tmp_usiptr,8,sizeof(unsigned short int),outfiles->fluxes);   

	  free((char *)tmp_usiptr);
        }
        else {
          /* Write an ASCII forcing file */
          fprintf(outfiles->fluxes,"%.7e \t%.7e \t%.7e \t%.7e \t%.7e \t%.7e \t%.7e \t%.7e\n", 
            atmos_alma->SWdown, atmos_alma->LWdown, atmos_alma->Tair,
            atmos_alma->Qair, atmos_alma->Psurf, atmos_alma->Rainf,
            atmos_alma->Snowf, atmos_alma->Wind);
        }
      }
      else {
        // VIC Standard
        if (options.BINARY_OUTPUT) {
  	  /* Write a binary forcing file */
  	  tmp_siptr = (short int *)calloc(1,sizeof(short int));
	  tmp_usiptr = (unsigned short int *)calloc(1,sizeof(unsigned short int));
	  /* precipitation * 40 */
	  tmp_usiptr[0] = (unsigned short int)(atmos[rec].prec[j]*40.);
	  fwrite(tmp_usiptr,1,sizeof(unsigned short int),outfiles->fluxes);
	  /* air temperature * 100 */
	  tmp_siptr[0] = (short int)(atmos[rec].air_temp[j]*100.);
	  fwrite(tmp_siptr,1,sizeof(short int),outfiles->fluxes);
	  /* shortwave * 50 */
	  tmp_usiptr[0] = (unsigned short int)(atmos[rec].shortwave[j]*50.);
	  fwrite(tmp_usiptr,1,sizeof(unsigned short int),outfiles->fluxes);
	  /* longwave * 80 */
	  tmp_usiptr[0] = (unsigned short int)(atmos[rec].longwave[j]*80.);
	  fwrite(tmp_usiptr,1,sizeof(unsigned short int),outfiles->fluxes);
	  /* density * 100 */
	  tmp_usiptr[0] = (unsigned short int)(atmos[rec].density[j]*100.);
	  fwrite(tmp_usiptr,1,sizeof(unsigned short int),outfiles->fluxes);
	  /* pressure * 0.1 */
	  tmp_usiptr[0] = (unsigned short int)(atmos[rec].pressure[j]*0.1);
	  fwrite(tmp_usiptr,1,sizeof(unsigned short int),outfiles->fluxes);
	  /* vp * 1 */
	  tmp_siptr[0] = (short int)(atmos[rec].vp[j]*1.);
	  fwrite(tmp_siptr,1,sizeof(short int),outfiles->fluxes);
	  /* wind * 100 */
	  tmp_usiptr[0] = (unsigned short int)(atmos[rec].wind[j]*100.);
	  fwrite(tmp_usiptr,1,sizeof(unsigned short int),outfiles->fluxes);   

	  free((char *)tmp_siptr);
	  free((char *)tmp_usiptr);
        }
        else {
	  /* Write an ASCII forcing file */
	  fprintf(outfiles->fluxes,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", 
		atmos[rec].prec[j], atmos[rec].air_temp[j], 
		atmos[rec].shortwave[j], atmos[rec].longwave[j], 
		atmos[rec].density[j], atmos[rec].pressure[j]/1000., 
		atmos[rec].vp[j]/1000., atmos[rec].wind[j]);
	}
      }
    }
  }

  free ((char *)atmos_alma);

}
#endif /* OUTPUT_FORCE */
