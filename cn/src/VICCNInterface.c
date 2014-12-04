#include <stdio.h>
#include <stdlib.h>
#include <vic_def.h>
#include <vic_run.h>
#include <string.h>
#define MAX_DIST 2
#define MAX_PFT 20
#define NUMRAD 2

extern void vic2clmtype_(int *nlevgrnd, int *rec, int *nrec, \
			 int *adspinup, int *yr, int *mo, int *day, \
			 int *secs, double *jday, int *yrnxt, double *jdynxt, \
			 double *dt, double *lat, double *lon, \
			 int *begg, int *endg, int *begc,      \
			 int *endc, int *begp, int *endp, \
			 int *num_soilc, int *num_soilp,	      \
			 double *psfc, double *Tair, double *vp, double *vpd, \
			 double *lwrad, double *swrad, double swrd[], \
			 double swri[], double alb[], \
			 double z[][MAX_NODES+2], double dz[][MAX_NODES+2], \
			 double baseflow[],	      \
			 double moist[][MAX_NODES+2], \
			 double ice[][MAX_NODES+2], \
			 double tsoisno[][MAX_NODES+2], \
			 double t2m[], double tveg[][MAX_PFT+1], \
			 double snowdep[], double fwet[][MAX_PFT+1], 	\
			 double rootf[][MAX_PFT][MAX_NODES],	\
			 double satpsi[][MAX_BANDS][MAX_NODES],		\
			 double soipsi[][MAX_BANDS][MAX_NODES], \
			 double coeff[][MAX_BANDS][MAX_NODES], \
			 double rveg[][MAX_PFT+1], double *zo, double *zos, \
			 double zov[][MAX_PFT+1], double displ[][MAX_PFT+1], \
			 double LAI[][MAX_PFT+1], double soilcfast[], \
			 double soilcmid[], double soilcslo1[], \
			 double soilcslo2[], double litrlabc[], \
			 double litrcellc[], double litrligc[], \
			 double cwdc[], double leafc[][MAX_PFT+1], \
			 double frootc[][MAX_PFT+1], \
			 double livestemc[][MAX_PFT+1], \
			 double deadstemc[][MAX_PFT+1], \
			 double livecrootc[][MAX_PFT+1], \
			 double deadcrootc[][MAX_PFT+1], \
			 double woodc[][MAX_PFT+1], double soilnfast[], \
			 double soilnmid[], double soilnslo1[], \
			 double soilnslo2[], double soilminn[], \
			 double litrlabn[], double litrcelln[], \
			 double litrlign[], double cwdn[], \
			 double leafn[][MAX_PFT+1], \
			 double frootn[][MAX_PFT+1],	\
			 double livestemn[][MAX_PFT+1], \
			 double deadstemn[][MAX_PFT+1], \
			 double livecrootn[][MAX_PFT+1], \
			 double deadcrootn[][MAX_PFT+1], \
			 double totvegc[][MAX_PFT+1], \
			 double totlitc[], double totsomc[], \
			 double gpp2[][MAX_PFT+1], double gpp[][MAX_PFT+1], \
			 double npp[][MAX_PFT+1], double darkr[][MAX_PFT+1], double mr[][MAX_PFT+1], \
			 double gr[][MAX_PFT+1], double ar[][MAX_PFT+1], \
			 double hr[], double lithr[], \
			 double nee[], double nep[],	    \
			 double dormant_flag[][MAX_PFT+1],	  \
			 double days_active[][MAX_PFT+1], \
			 double onset_flag[][MAX_PFT+1], \
			 double onset_counter[][MAX_PFT+1], \
			 double onset_gddflag[][MAX_PFT+1], \
			 double onset_fdd[][MAX_PFT+1], \
			 double onset_gdd[][MAX_PFT+1], \
			 double onset_swi[][MAX_PFT+1], \
			 double offset_flag[][MAX_PFT+1], \
			 double offset_counter[][MAX_PFT+1], \
			 double offset_fdd[][MAX_PFT+1], \
			 double offset_swi[][MAX_PFT+1], \
			 double lgsf[][MAX_PFT+1], double bglfr[][MAX_PFT+1], \
			 double bgtr[][MAX_PFT+1], double dayl[][MAX_PFT+1], \
			 double prev_dayl[][MAX_PFT+1], \
			 double annavg_t2m[][MAX_PFT+1], \
			 double tempavg_t2m[][MAX_PFT+1], \
			 double availc[][MAX_PFT+1], \
			 double xsmrpool_recover[][MAX_PFT+1], \
			 double alloc_pnow[][MAX_PFT+1], \
			 double c_allometry[][MAX_PFT+1], \
			 double n_allometry[][MAX_PFT+1], \
			 double plant_ndemand[][MAX_PFT+1], \
			 double tempsum_potential_gpp[][MAX_PFT+1], \
			 double annsum_potential_gpp[][MAX_PFT+1], \
			 double tempmax_retransn[][MAX_PFT+1], \
			 double annmax_retransn[][MAX_PFT+1], \
			 double avail_retransn[][MAX_PFT+1], \
			 double plant_nalloc[][MAX_PFT+1], \
			 double plant_calloc[][MAX_PFT+1], \
			 double excess_cflux[][MAX_PFT+1], \
			 double downreg[][MAX_PFT+1], \
			 double prev_leafc_to_litter[][MAX_PFT+1], \
			 double prev_frootc_to_litter[][MAX_PFT+1], \
			 double tempsum_npp[][MAX_PFT+1], \
			 double annsum_npp[][MAX_PFT+1], \
			 double leafc_storage[][MAX_PFT+1], \
			 double leafc_xfer[][MAX_PFT+1], \
			 double frootc_storage[][MAX_PFT+1], \
			 double frootc_xfer[][MAX_PFT+1], \
			 double livestemc_storage[][MAX_PFT+1], \
			 double livestemc_xfer[][MAX_PFT+1], \
			 double deadstemc_storage[][MAX_PFT+1], \
			 double deadstemc_xfer[][MAX_PFT+1], \
			 double livecrootc_storage[][MAX_PFT+1], \
			 double livecrootc_xfer[][MAX_PFT+1], \
			 double deadcrootc_storage[][MAX_PFT+1], \
			 double deadcrootc_xfer[][MAX_PFT+1], \
			 double gresp_storage[][MAX_PFT+1], \
			 double gresp_xfer[][MAX_PFT+1], \
			 double cpool[][MAX_PFT+1], \
			 double xsmrpool[][MAX_PFT+1], \
			 double pft_ctrunc[][MAX_PFT+1], \
			 double leafn_storage[][MAX_PFT+1], \
			 double leafn_xfer[][MAX_PFT+1], \
			 double frootn_storage[][MAX_PFT+1], \
			 double frootn_xfer[][MAX_PFT+1], \
			 double livestemn_storage[][MAX_PFT+1], \
			 double livestemn_xfer[][MAX_PFT+1], \
			 double deadstemn_storage[][MAX_PFT+1], \
			 double deadstemn_xfer[][MAX_PFT+1], \
			 double livecrootn_storage[][MAX_PFT+1], \
			 double livecrootn_xfer[][MAX_PFT+1], \
			 double deadcrootn_storage[][MAX_PFT+1], \
			 double deadcrootn_xfer[][MAX_PFT+1], \
			 double retransn[][MAX_PFT+1], \
			 double npool[][MAX_PFT+1], \
			 double pft_ntrunc[][MAX_PFT+1], \
			 double decl[], double fpi[], double fpg[], \
			 double annsum_counter[], double cannsum_npp[], \
			 double cannavg_t2m[], double watfc[][MAX_NODES], \
			 double me[], double fire_prob[], \
			 double mean_fire_prob[],		      \
			 double fireseasonl[], double farea_burned[], \
			 double ann_farea_burned[], double seedc[], \
			 double col_ctrunc[], double totcolc[], \
			 double prod10c[], double prod100c[], \
			 double seedn[], double col_ntrunc[], \
			 double totcoln[], double prod10n[], \
			 double prod100n[], double litfall[][MAX_PFT+1], \
			 double fpsn[][MAX_PFT+1], double ci[][MAX_PFT+1], \
			 double rc[][MAX_PFT+1], double apar[][MAX_PFT+1], \
			 int *init_state);

int VICCNInterface(int                 rec,
                   int                 Nveg,
		   int                 Npft,
		   dmy_struct          *dmy, 
		   global_param_struct *global_param,
                   atmos_data_struct   *atmos,
		   all_vars_struct    *all_vars,
		   soil_con_struct     *soil_con,
		   veg_con_struct      *veg_con,
                   veg_lib_struct      *veg_lib,
                   cn_data_struct      *cn)

/**********************************************************************
	VICCNInterface	Michael Brunke		July 23, 2014

Prepares VIC data for use in CN.

Modifications:
  24-Jul-2014 Replaced dist by bands.                                   MAB
  15-Oct-2014 Added veg_var and veg_hist structures as input.           MAB

**********************************************************************/

{

  extern option_struct options;

  double z[MAX_BANDS][MAX_NODES+2];             /* Layer depth (m) */
  double dz[MAX_BANDS][MAX_NODES+2];            /* Layer thickness (m) */
  double baseflow[MAX_BANDS];                   /* Baseflow */
  double moist[MAX_BANDS][MAX_NODES+2];          /* Soil moisture */
  double ice[MAX_BANDS][MAX_NODES+2];            /* Soil ice lens */
  double t_soisno[MAX_BANDS][MAX_NODES+2];     /* Snow/soil temperature */
  double t2m[MAX_PFT];                         /* 2-m air temperature */
  double snowdep[MAX_BANDS];                    /* Snow depth */
  double fwet[MAX_BANDS][MAX_PFT+1];            /* Wet veg. fraction */
  double rootfr[MAX_BANDS][MAX_PFT][MAX_NODES];/* Root fraction in ea. layer */
  double psisat[2][MAX_BANDS][MAX_NODES];      /* Saturated matric potential */
  double soipsi[2][MAX_BANDS][MAX_NODES];         /* Matric potential */
  double coeff[2][MAX_BANDS][MAX_NODES];          /* Clapp-Hornberger coeff. */
  double rveg[MAX_BANDS][MAX_PFT+1];           /* Leaf resistance */
  double *rcan;                               /* Canopy resistance */
  double LAI[MAX_BANDS][MAX_PFT+1];             /* Leaf area index */
  double soilcfast[MAX_BANDS];                  /* Fast soil C pool */
  double soilcmid[MAX_BANDS];                   /* Medium soil C pool */
  double soilcslo1[MAX_BANDS];                  /* Slow soil C pool */
  double soilcslo2[MAX_BANDS];                  /* Slowest soil C pool */
  double litrlabc[MAX_BANDS];                   /* Litter labile C pool */
  double litrcellc[MAX_BANDS];                  /* Litter cellulose C pool */
  double litrligc[MAX_BANDS];                   /* Litter lignin C pool */
  double cwdc[MAX_BANDS];                       /* Coarse woody debris C pool */
  double leafc[MAX_BANDS][MAX_PFT+1];           /* Leaf C pool */
  double frootc[MAX_BANDS][MAX_PFT+1];          /* Fine root C pool */
  double livestemc[MAX_BANDS][MAX_PFT+1];       /* Live stem C pool */
  double deadstemc[MAX_BANDS][MAX_PFT+1];       /* Dead stem C pool */
  double livecrootc[MAX_BANDS][MAX_PFT+1];      /* Live coarse root C pool */
  double deadcrootc[MAX_BANDS][MAX_PFT+1];      /* Dead coarse root C pool */
  double woodc[MAX_BANDS][MAX_PFT+1];           /* Wood C pool */
  double totvegc[MAX_BANDS][MAX_PFT+1];         /* Total vegetation C */
  double totlitc[MAX_BANDS];                    /* Total litter C */
  double totsomc[MAX_BANDS];                    /* Total SOM C */
  double soilnfast[MAX_BANDS];                  /* Fast soil N pool */
  double soilnmid[MAX_BANDS];                   /* Medium soil N pool */
  double soilnslo1[MAX_BANDS];                  /* Slow soil N pool */
  double soilnslo2[MAX_BANDS];                  /* Slowest soil N pool */
  double soilminn[MAX_BANDS];                   /* Mineral soil N pool */
  double litrlabn[MAX_BANDS];                   /* Litter labile N pool */
  double litrcelln[MAX_BANDS];                  /* Litter cellulose N pool */
  double litrlign[MAX_BANDS];                   /* Litter lignin N pool */
  double cwdn[MAX_BANDS];                       /* Coarse woody debris N pool */
  double leafn[MAX_BANDS][MAX_PFT+1];           /* Leaf N pool */
  double frootn[MAX_BANDS][MAX_PFT+1];          /* Fine root N pool */
  double livestemn[MAX_BANDS][MAX_PFT+1];       /* Live stem N pool */
  double deadstemn[MAX_BANDS][MAX_PFT+1];       /* Dead stem N pool */
  double livecrootn[MAX_BANDS][MAX_PFT+1];      /* Live coarse root N pool */
  double deadcrootn[MAX_BANDS][MAX_PFT+1];      /* Dead coarse root N pool */
  double totlitn[MAX_BANDS];                    /* Total litter N */
  double gpp2[MAX_BANDS][MAX_PFT+1];            /* GPP before downregulation */
  double gpp[MAX_BANDS][MAX_PFT+1];             /* Gross primary production */
  double npp[MAX_BANDS][MAX_PFT+1];             /* Net primary production */
  double darkr[MAX_BANDS][MAX_PFT+1];           /* Leaf maint respiration */
  double mr[MAX_BANDS][MAX_PFT+1];              /* Maintenance respiration */
  double gr[MAX_BANDS][MAX_PFT+1];              /* Growth respiration */
  double ar[MAX_BANDS][MAX_PFT+1];              /* Autotrophic respiration */
  double hr[MAX_BANDS];                         /* Heterotrophic respiration */
  double lithr[MAX_BANDS];                      /* Litter hetero respiration */
  double nee[MAX_BANDS];                        /* Net ecosystem exchange */
  double nep[MAX_BANDS];                        /* Net ecosystem production */

  double dormant_flag[MAX_BANDS][MAX_PFT+1];    /* Dormancy flag */
  double days_active[MAX_BANDS][MAX_PFT+1];     /* # days since last dormancy */
  double onset_flag[MAX_BANDS][MAX_PFT+1];      /* Onset flag */
  double onset_counter[MAX_BANDS][MAX_PFT+1];   /* Onset days counter */
  double onset_gddflag[MAX_BANDS][MAX_PFT+1];  /* Onset flag for grow deg sum */
  double onset_fdd[MAX_BANDS][MAX_PFT+1];       /* Onset freeze days counter */
  double onset_gdd[MAX_BANDS][MAX_PFT+1];       /* Onset growing deg days */
  double onset_swi[MAX_BANDS][MAX_PFT+1];       /* Onset soil water index */
  double offset_flag[MAX_BANDS][MAX_PFT+1];      /* Offset flag */
  double offset_counter[MAX_BANDS][MAX_PFT+1];   /* Offset days counter */
  double offset_fdd[MAX_BANDS][MAX_PFT+1];       /* Offset freeze days counter */
  double offset_swi[MAX_BANDS][MAX_PFT+1];       /* Offset soil water index */
  double lgsf[MAX_BANDS][MAX_PFT+1];            /* Long growing season factor */
  double bglfr[MAX_BANDS][MAX_PFT+1];           /* Background litterfall rate */
  double bgtr[MAX_BANDS][MAX_PFT+1];            /* Background transfer growth */
  double dayl[MAX_BANDS][MAX_PFT+1];            /* Daylength (s) */
  double prev_dayl[MAX_BANDS][MAX_PFT+1];       /* Previous daylength (s) */
  double annavg_t2m[MAX_BANDS][MAX_PFT+1];      /* Annual avg 2m air temp. */
  double tempavg_t2m[MAX_BANDS][MAX_PFT+1];     /* Temp. avg. 2m air temp. */
  double availc[MAX_BANDS][MAX_PFT+1];          /* C flux avail for alloc */
  double xsmrpool_recover[MAX_BANDS][MAX_PFT+1];/* C flx assigned to recovery */
  double alloc_pnow[MAX_BANDS][MAX_PFT+1];    /* Fract of alloc as new growth */
  double c_allometry[MAX_BANDS][MAX_PFT+1];     /* C allocation index */
  double n_allometry[MAX_BANDS][MAX_PFT+1];     /* N allocation index */
  double plant_ndemand[MAX_BANDS][MAX_PFT+1];   /* N flux to support GPP */
  double tempsum_potential_gpp[MAX_BANDS][MAX_PFT+1];/* Temp ann sum of pot GPP */
  double annsum_potential_gpp[MAX_BANDS][MAX_PFT+1];/* Ann sum of pot GPP */
  double tempmax_retransn[MAX_BANDS][MAX_PFT+1];/* Temp ann max retrans N */
  double annmax_retransn[MAX_BANDS][MAX_PFT+1]; /* Ann max retrans N pool */
  double avail_retransn[MAX_BANDS][MAX_PFT+1];  /* N flux from retrans pool */
  double plant_nalloc[MAX_BANDS][MAX_PFT+1];    /* Total allocated N flux */
  double plant_calloc[MAX_BANDS][MAX_PFT+1];    /* Total allocated C flux */
  double excess_cflux[MAX_BANDS][MAX_PFT+1];    /* C flux not allocated */
  double downreg[MAX_BANDS][MAX_PFT+1];        /* Fract reduct GPP from N lim */
  double prev_leafc_to_litter[MAX_BANDS][MAX_PFT+1];/* Prev leaf C litterfall */
  double prev_frootc_to_litter[MAX_BANDS][MAX_PFT+1];/* Prev froot C litterfall */
  double tempsum_npp[MAX_BANDS][MAX_PFT+1];     /* Temp ann sum of NPP */
  double annsum_npp[MAX_BANDS][MAX_PFT+1];      /* Annual sum of NPP */
  double leafc_storage[MAX_BANDS][MAX_PFT+1];   /* Leaf C storage */
  double leafc_xfer[MAX_BANDS][MAX_PFT+1];      /* Leaf C transfer */
  double frootc_storage[MAX_BANDS][MAX_PFT+1];  /* Fine root C storage */
  double frootc_xfer[MAX_BANDS][MAX_PFT+1];     /* Fine root C transfer */
  double livestemc_storage[MAX_BANDS][MAX_PFT+1];/* Live stem C storage */
  double livestemc_xfer[MAX_BANDS][MAX_PFT+1];  /* Live stem C transfer */
  double deadstemc_storage[MAX_BANDS][MAX_PFT+1]; /* Dead stem C storage */
  double deadstemc_xfer[MAX_BANDS][MAX_PFT+1];  /* Dead stem C transfer */
  double livecrootc_storage[MAX_BANDS][MAX_PFT+1];/* Live coarse root C storage */
  double livecrootc_xfer[MAX_BANDS][MAX_PFT+1]; /* Live coarse root C transfer */
  double deadcrootc_storage[MAX_BANDS][MAX_PFT+1];/* Dead coarse root C storage */
  double deadcrootc_xfer[MAX_BANDS][MAX_PFT+1]; /* Dead coarse root C transfer */
  double gresp_storage[MAX_BANDS][MAX_PFT+1];   /* Growth respiration storage */
  double gresp_xfer[MAX_BANDS][MAX_PFT+1];      /* Growth respiration transfer */
  double cpool[MAX_BANDS][MAX_PFT+1];           /* Temp photosynthate C pool */
  double xsmrpool[MAX_BANDS][MAX_PFT+1];        /* Abstract C pool */
  double pft_ctrunc[MAX_BANDS][MAX_PFT+1];      /* PFT sink for C truncation */
  double leafn_storage[MAX_BANDS][MAX_PFT+1];   /* Leaf N storage */
  double leafn_xfer[MAX_BANDS][MAX_PFT+1];      /* Leaf N storage */
  double frootn_storage[MAX_BANDS][MAX_PFT+1];  /* Fine root N storage */
  double frootn_xfer[MAX_BANDS][MAX_PFT+1];     /* Fine root N transfer */
  double livestemn_storage[MAX_BANDS][MAX_PFT+1]; /* Live stem N storage */
  double livestemn_xfer[MAX_BANDS][MAX_PFT+1];  /* Live stem N transfer */
  double deadstemn_storage[MAX_BANDS][MAX_PFT+1]; /* Dead stem N storage */
  double deadstemn_xfer[MAX_BANDS][MAX_PFT+1];  /* Dead stem N transfer */
  double livecrootn_storage[MAX_BANDS][MAX_PFT+1]; /* Live coarse root N storage */
  double livecrootn_xfer[MAX_BANDS][MAX_PFT+1]; /* Live coarse root N transfer */
  double deadcrootn_storage[MAX_BANDS][MAX_PFT+1]; /* Dead coarse root N storage */
  double deadcrootn_xfer[MAX_BANDS][MAX_PFT+1]; /* Dead coarse root N transfer */
  double retransn[MAX_BANDS][MAX_PFT+1];        /* Retranslocated N */
  double npool[MAX_BANDS][MAX_PFT+1];           /* Temp. photosynthate N pool */
  double pft_ntrunc[MAX_BANDS][MAX_PFT+1];      /* PFT sink for N truncation */
  double decl[MAX_BANDS];                       /* Solar declination angle */
  double fpi[MAX_BANDS];                        /* Fract. pot. immobilization */
  double fpg[MAX_BANDS];                        /* Fract. potential GPP */
  double annsum_counter[MAX_BANDS];             /* Secs since last ann accum. turnover */
  double cannsum_npp[MAX_BANDS];                /* Annual sum of NPP */
  double cannavg_t2m[MAX_BANDS];                /* Annual avg. of 2-m air temp. */
  double watfc[MAX_BANDS][MAX_NODES];           /* Volumetric soil moisture at field capcity */
  double me[MAX_BANDS];                         /* Moisture of extinction */
  double fire_prob[MAX_BANDS];                  /* Daily fire probability */
  double mean_fire_prob[MAX_BANDS];             /* E-fold mean daily fire prob. */
  double fireseasonl[MAX_BANDS];                /* Ann. fire season length */
  double farea_burned[MAX_BANDS];               /* Timestep fract. area burned */
  double ann_farea_burned[MAX_BANDS];           /* Ann. tot. fract. area burned */
  double seedc[MAX_BANDS];                      /* Col-lev pool for seeding new PFTs */
  double col_ctrunc[MAX_BANDS];                 /* Col-lev sink for C trunc */
  double totcolc[MAX_BANDS];                    /* Total column C */
  double prod10c[MAX_BANDS];                    /* Wood product C pool, 10-yr lifespan */
  double prod100c[MAX_BANDS];                   /* Wood product C pool, 100-yr lifespan */
  double seedn[MAX_BANDS];                      /* Col-lev pool for seeding new PFTs */
  double col_ntrunc[MAX_BANDS];                 /* Col-lev sink for N trunc */
  double totcoln[MAX_BANDS];                    /* Total column N */
  double prod10n[MAX_BANDS];                    /* Wood product N pool, 10-yr lifespan */
  double prod100n[MAX_BANDS];                   /* Wood product N pool, 100-yr lifespan */
  double litfall[MAX_BANDS][MAX_PFT+1];         /* Litterfall */
  double fpsn[MAX_BANDS][MAX_PFT+1];             /* Photosynthesis */
  double ci[MAX_BANDS][MAX_PFT+1];               /* Intracellular C02 */
  double rc[MAX_BANDS][MAX_PFT+1];               /* Canopy stomatal resistance */
  double apar[MAX_BANDS][MAX_PFT+1];            /* absorbed PAR */

  double Tair;                                 /* Air temperature */
  double *Tcan;                                /* Canopy temperature */
  double Tveg[MAX_BANDS][MAX_PFT+1];            /* Leaf temperature */
  double psfc;                                 /* Surface pressure */
  double vp;                                   /* Atmos. vapor pressure */
  double vpd;                                  /* Vapor pressure deficit */
  double lwrad;                                /* LW radiation */
  double swrad;                                /* SW radiation */
  double swrd[NUMRAD];                         /* direct SW radiation */
  double swri[NUMRAD];                         /* diffuse SW radiation */
  double alb[MAX_BANDS];                        /* surface albedo */
  double zo, zos, zov[MAX_BANDS][MAX_PFT+1];      /* roughness lengths */
  double displ[MAX_BANDS][MAX_PFT+1];             /* displacement height */
  double lat;
  double lon;

  int i;
  int iveg;
  int band;
  int Nbands;
  int lidx, nidx;
  int Nlayers;
  int Nnodes;
  int Nrecs;
  int adspinup;
  int lbg;
  int ubg;
  int lbc;
  int ubc;
  int lbp;
  int ubp;
  int num_soilc;
  int num_soilp;
  int init_state;
  int year, month, day, hour, secs, yrnxt;
  double jday, jdynxt;
  int err = 0;
  int veg_class;
  double factor;
  double precip;
  double sucsat[2][MAX_LAYERS];      /* Saturated matric potential in layers */
  double soisuc[2][MAX_LAYERS];      /* Soil matric potential */
  double theta[MAX_LAYERS];          /* Volumetric water content */
  double watsat[MAX_LAYERS];         /* Saturated volumetric water content */
  double psand;                      /* % sand */
  double pclay;                      /* % clay */
  double bsw[2][MAX_LAYERS];         /* Clapp-Hornberger coefficient */
  double ksat;                       /* Saturated water conductivity */
  double watfrac[MAX_LAYERS];        /* Water content at field capacity */
  double fsat;                       /* Saturation fraction */
  double Lsum;                        /* cumulative depth of moisture layer */
  double fsum[MAX_PFT + 1];
  char PAST_BOTTOM;

  cell_data_struct **cell;
  energy_bal_struct **energy;
  snow_data_struct  **snow;
  veg_var_struct **veg;

  double dt;
  double b;

  /* Set local pointers */
  cell = all_vars->cell;
  energy = all_vars->energy;
  snow = all_vars->snow;
  veg = all_vars->veg_var;

  if(options.CARBON == CN_ADECOMP)
    adspinup = 1;
  else
    adspinup = 0;
  if(options.INIT_STATE)
    init_state = 1;
  else
    init_state = 0;

  /* Set number of snow bands */
  Nbands = options.SNOW_BAND;

  /* Set number of soil layers */
  Nlayers = options.Nlayer;

  /* Set number of nodes */
  Nnodes = options.Nnode;

  /* Set number of records */
  Nrecs = global_param->nrecs;

  /* Convert VIC time step to seconds */
  dt = (double) global_param->dt * 3600.0;

  /* Get current date from dmy */
  year = dmy[rec].year;
  month = dmy[rec].month;
  day = dmy[rec].day;
  hour = dmy[rec].hour;
  secs = hour * 3600.0;
  jday = dmy[rec].day_in_year;

  /* Get next date from dmy, MAB 10/8/13 */
  if(rec + 1 < Nrecs)
    {
    yrnxt = dmy[rec+1].year;
    jdynxt = dmy[rec+1].day_in_year;
    }
  else
    {
    yrnxt = -1;
    jdynxt = -1.0;
    }

  /* Get latitude/longitude, MAB 8/27/13 */
  lat = (double) soil_con->lat;
  lon = (double) soil_con->lng;

  /* Get soil/snow roughness, MAB 10/11/13 */
  zo = soil_con->rough;
  zos = soil_con->snow_rough;

  /* Assign atmos structure variables to local variables, MAB 8/19/13 */

  Tair = atmos->air_temp[NR];
  vp = atmos->vp[NR];
  vpd = atmos->vpd[NR];
  psfc = atmos->pressure[NR];
  lwrad = atmos->longwave[NR];
  swrad = atmos->shortwave[NR];
  precip = atmos->prec[NR];

  /* Calculation of direct and diffuse SW radiation elements, MAB 11/26/13 */
  /* Based off of CLM's data atmosphere model.  THIS SHOULD NOT BE ACTIVE  */
  /* WHEN COUPLED TO WRF!  These should be provided by WRF directly if     */
  /* possible. */

  swrd[0] = 0.28 * swrad;
  swrd[1] = 0.31 * swrad;
  swri[0] = 0.24 * swrad;
  swri[1] = 0.17 * swrad;

  factor = 1.0;
  if(lat > -60.0 && lat < -50.0)
    factor = 1.0 - (lat + 60.0) * (0.05 / 10.0);
  else if(lat >= -50.0 && lat <= 30.0)
    factor = 0.95;
  else if(lat > 30.0 && lat < 40.0)
    factor = 1.0 - (40.0 - lat) * (0.05 / 10.0);

  swrd[0] *= factor;
  swri[0] *= factor;
  swrd[1] *= factor;
  swri[1] *= factor;

  /* Layer properties, MAB 10/10/13 */

  for(band = 0; band < Nbands; band++)
    {
      for(nidx = 0; nidx < Nnodes; nidx++)
	{
	  z[band][nidx+2] = soil_con->Zsum_node[nidx];
	  dz[band][nidx+2] = soil_con->dz_node[nidx];
	}
    }

  /* Average baseflow and soil moisture over all bands for each 
     vegetated type */

  for(band = 0; band < Nbands; band++)
    {

      baseflow[band] = 0.0;
      for(nidx = 0; nidx < Nnodes + 2; nidx++)
	{
	  moist[band][nidx] = 0.0;
          ice[band][nidx] = 0.0;
	}

      for(iveg = 0; iveg <= Nveg; iveg++)
	{

	  if(iveg < Nveg)
	    {

		  baseflow[band] += cell[iveg][band].baseflow;

		  moist[band][0] += snow[iveg][band].surf_water;
		  moist[band][1] += snow[iveg][band].pack_water;
                  ice[band][1] += snow[iveg][band].swq;
		  for(nidx = 0; nidx < Nnodes; nidx++)
		    {
		      moist[band][nidx+2] += energy[iveg][band].moist[nidx];
		      ice[band][nidx+2] += (energy[iveg][band].ice[nidx] * \
			snow[iveg][band].density);
		    }

	    }

	}

      baseflow[band] /= (Nveg - 1);

      for(nidx = 0; nidx < Nnodes + 2; nidx++)
	{
	  moist[band][nidx] /= (Nveg - 1);
          ice[band][nidx] /= (Nveg - 1);
	  if(nidx < 2)
	    {
	    moist[band][nidx] *= RHO_W;
	    }
	  else
	    {
	    moist[band][nidx] *= (soil_con->dz_node[nidx-2] * RHO_W);
	    ice[band][nidx] *= (soil_con->dz_node[nidx-2] * ice_density);
	    }
	}

    }

  /* Combine snow and soil temperatures into one array */

  for(band = 0; band < Nbands; band++)
    {

      alb[band] = 0.0;
      for(nidx = 0; nidx < Nnodes; nidx++)
	{
	  t_soisno[band][nidx] = 0.0;
	}

      for(iveg = 0; iveg < Nveg; iveg++)
	{
	      /* Add snow and pack temperatures */

	      t_soisno[band][0] += snow[iveg][band].surf_temp;
	      t_soisno[band][1] += snow[iveg][band].pack_temp;

              /* Add node temperatures */
	      for(nidx = 0; nidx < Nnodes; nidx++)
		{
		  t_soisno[band][nidx + 2] += energy[iveg][band].T[nidx];
		}

	      alb[band] += energy[iveg][band].AlbedoUnder;

	}

      /* Average temperatures */

      for(nidx = 0; nidx < Nnodes + 2; nidx++)
	{
	  t_soisno[band][nidx] /= Nveg;
	}

      /* Average albedo */
      alb[band] /= Nveg;

    }

  /* Average 2-m air temperatures over snow sub-steps */

  for(iveg = 0; iveg < Npft; iveg++)
    {
      t2m[iveg] = atmos->air_temp[NR];
    }

  /* Average canopy temperature and resistance */

  Tcan = (double *) calloc(Nveg, sizeof(double));
  rcan = (double *) calloc(Nveg, sizeof(double));
  for(iveg = 0; iveg < Nveg; iveg++)
    {
      Tcan[iveg] = 0.0;
      rcan[iveg] = 0.0;
    }

  for(band = 0; band < Nbands; band++)
    {

      for(iveg = 0; iveg <= MAX_PFT; iveg++)
	{
	  Tveg[band][iveg] = 0.0;
	  rveg[band][iveg] = 0.0;
          zov[band][iveg] = 0.0;
          displ[band][iveg] = 0.0;
	  fsum[iveg] = 0.0;
	}

      for(iveg = 0; iveg < Nveg; iveg++)
	{
	  Tcan[iveg] += energy[iveg][band].Tcanopy;
	  rcan[iveg] += cell[iveg][band].aero_resist[1];

	  switch(veg_con[iveg].veg_class)
	    {
	    case 0: Tveg[band][2] += (Tcan[iveg] * veg_con[iveg].Cv);
	      rveg[band][2] += (rcan[iveg] * veg_con[iveg].Cv);
              zov[band][2] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv);
              displ[band][2] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] * \
		 veg_con[iveg].Cv);
	      fsum[2] += veg_con[iveg].Cv;
	      break;
	    case 1: Tveg[band][5] += (Tcan[iveg] * veg_con[iveg].Cv);
	      rveg[band][5] += (rcan[iveg] * veg_con[iveg].Cv);
              zov[band][5] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv);
              displ[band][5] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv);
	      fsum[5] += veg_con[iveg].Cv;
	      break;
	    case 2: Tveg[band][3] += (Tcan[iveg] * veg_con[iveg].Cv);
	      rveg[band][3] += (rcan[iveg] * veg_con[iveg].Cv);
              zov[band][3] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv);
              displ[band][3] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv);
	      fsum[3] += veg_con[iveg].Cv;
	      break;
	    case 3: Tveg[band][8] += (Tcan[iveg] * veg_con[iveg].Cv);
	      rveg[band][8] += (rcan[iveg] * veg_con[iveg].Cv);
              zov[band][8] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv);
              displ[band][8] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv);
	      fsum[8] += veg_con[iveg].Cv;
	      break;
	    case 4: Tveg[band][2] += (Tcan[iveg] * veg_con[iveg].Cv * 0.5);
	      Tveg[band][8] += (Tcan[iveg] * veg_con[iveg].Cv * 0.5);
	      rveg[band][2] += (rcan[iveg] * veg_con[iveg].Cv * 0.5);
	      rveg[band][8] += (rcan[iveg] * veg_con[iveg].Cv * 0.5);
              zov[band][8] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.5);
              zov[band][2] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.5);
              displ[band][8] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.5);
              displ[band][2] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.5);
	      fsum[2] += veg_con[iveg].Cv * 0.5;
	      fsum[8] += veg_con[iveg].Cv * 0.5;
	      break;
	    case 5: Tveg[band][2] += (Tcan[iveg] * veg_con[iveg].Cv * 0.4);
	      Tveg[band][8] += (Tcan[iveg] * veg_con[iveg].Cv * 0.4);
	      Tveg[band][9] += (Tcan[iveg] * veg_con[iveg].Cv * 0.1);
	      Tveg[band][10] += (Tcan[iveg] * veg_con[iveg].Cv * 0.1);
	      rveg[band][2] += (rcan[iveg] * veg_con[iveg].Cv * 0.4);
	      rveg[band][8] += (rcan[iveg] * veg_con[iveg].Cv * 0.4);
	      rveg[band][9] += (rcan[iveg] * veg_con[iveg].Cv * 0.1);
	      rveg[band][10] += (rcan[iveg] * veg_con[iveg].Cv * 0.1);
              zov[band][2] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.4);
              zov[band][8] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.4);
              zov[band][9] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.1);
              zov[band][10] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.1);
              displ[band][2] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.4);
              displ[band][8] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.4);
              displ[band][9] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.1);
              displ[band][10] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.1);
	      fsum[2] += veg_con[iveg].Cv * 0.4;
	      fsum[8] += veg_con[iveg].Cv * 0.4;
	      fsum[9] += veg_con[iveg].Cv * 0.1;
	      fsum[10] += veg_con[iveg].Cv * 0.1;
	      break;
	    case 6: Tveg[band][14] += (Tcan[iveg] * veg_con[iveg].Cv * 0.1);
	      Tveg[band][13] += (Tcan[iveg] * veg_con[iveg].Cv * 0.3);
	      Tveg[band][12] += (Tcan[iveg] * veg_con[iveg].Cv * 0.2);
	      Tveg[band][2] += (Tcan[iveg] * veg_con[iveg].Cv * 0.2);
	      Tveg[band][8] += (Tcan[iveg] * veg_con[iveg].Cv * 0.2);
	      rveg[band][14] += (rcan[iveg] * veg_con[iveg].Cv * 0.1);
	      rveg[band][13] += (rcan[iveg] * veg_con[iveg].Cv * 0.3);
	      rveg[band][12] += (rcan[iveg] * veg_con[iveg].Cv * 0.2);
	      rveg[band][2] += (rcan[iveg] * veg_con[iveg].Cv * 0.2);
	      rveg[band][8] += (rcan[iveg] * veg_con[iveg].Cv * 0.2);
              zov[band][14] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.1);
              zov[band][13] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.3);
              zov[band][12] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.2);
              zov[band][2] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.2);
              zov[band][8] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.2);
              displ[band][14] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.1);
              displ[band][13] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.3);
              displ[band][12] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.2);
              displ[band][2] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.2);
              displ[band][8] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.2);
	      fsum[14] += veg_con[iveg].Cv * 0.1;
	      fsum[13] += veg_con[iveg].Cv * 0.3;
	      fsum[12] += veg_con[iveg].Cv * 0.2;
	      fsum[2] += veg_con[iveg].Cv * 0.2;
	      fsum[8] += veg_con[iveg].Cv * 0.2;
	      break;
	    case 7: Tveg[band][9] += (Tcan[iveg] * veg_con[iveg].Cv * 0.3);
	      Tveg[band][10] += (Tcan[iveg] * veg_con[iveg].Cv * 0.3);
	      Tveg[band][11] += (Tcan[iveg] * veg_con[iveg].Cv * 0.2);
	      Tveg[band][2] += (Tcan[iveg] * veg_con[iveg].Cv * 0.1);
	      Tveg[band][8] += (Tcan[iveg] * veg_con[iveg].Cv * 0.1);
	      rveg[band][9] += (rcan[iveg] * veg_con[iveg].Cv * 0.3);
	      rveg[band][10] += (rcan[iveg] * veg_con[iveg].Cv * 0.3);
	      rveg[band][11] += (rcan[iveg] * veg_con[iveg].Cv * 0.2);
	      rveg[band][2] += (rcan[iveg] * veg_con[iveg].Cv * 0.1);
	      rveg[band][8] += (rcan[iveg] * veg_con[iveg].Cv * 0.1);
              zov[band][9] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.3);
              zov[band][10] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.3);
              zov[band][11] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.2);
              zov[band][8] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.1);
              zov[band][2] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.1);
              displ[band][9] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.3);
              displ[band][10] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.3);
              displ[band][11] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.2);
              displ[band][8] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.1);
              displ[band][2] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.1);
	      fsum[9] += veg_con[iveg].Cv * 0.3;
	      fsum[10] += veg_con[iveg].Cv * 0.3;
	      fsum[11] += veg_con[iveg].Cv * 0.2;
	      fsum[8] += veg_con[iveg].Cv * 0.1;
	      fsum[2] += veg_con[iveg].Cv * 0.1;
	      break;
	    case 8: Tveg[band][9] += (Tcan[iveg] * veg_con[iveg].Cv * 0.2);
	      Tveg[band][10] += (Tcan[iveg] * veg_con[iveg].Cv * 0.2);
	      Tveg[band][11] += (Tcan[iveg] * veg_con[iveg].Cv * 0.2);
	      Tveg[band][14] += (Tcan[iveg] * veg_con[iveg].Cv * 0.05);
	      Tveg[band][13] += (Tcan[iveg] * veg_con[iveg].Cv * 0.15);
	      Tveg[band][12] += (Tcan[iveg] * veg_con[iveg].Cv * 0.1);
	      rveg[band][9] += (rcan[iveg] * veg_con[iveg].Cv * 0.2);
	      rveg[band][10] += (rcan[iveg] * veg_con[iveg].Cv * 0.2);
	      rveg[band][11] += (rcan[iveg] * veg_con[iveg].Cv * 0.2);
	      rveg[band][14] += (rcan[iveg] * veg_con[iveg].Cv * 0.05);
	      rveg[band][13] += (rcan[iveg] * veg_con[iveg].Cv * 0.15);
	      rveg[band][12] += (rcan[iveg] * veg_con[iveg].Cv * 0.1);
              zov[band][9] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.2);
              zov[band][10] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.2);
              zov[band][11] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.2);
              zov[band][14] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.05);
              zov[band][13] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.15);
              zov[band][12] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.1);
              displ[band][9] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.2);
              displ[band][10] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.2);
              displ[band][11] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.2);
              displ[band][14] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.05);
              displ[band][13] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.15);
              displ[band][12] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.1);
	      fsum[9] += veg_con[iveg].Cv * 0.2;
	      fsum[10] += veg_con[iveg].Cv * 0.2;
	      fsum[11] += veg_con[iveg].Cv * 0.2;
	      fsum[14] += veg_con[iveg].Cv * 0.05;
	      fsum[13] += veg_con[iveg].Cv * 0.15;
	      fsum[12] += veg_con[iveg].Cv * 0.1;
	      break;
	    case 9: Tveg[band][14] += (Tcan[iveg] * veg_con[iveg].Cv * 0.1);
	      Tveg[band][13] += (Tcan[iveg] * veg_con[iveg].Cv * 0.3);
	      Tveg[band][12] += (Tcan[iveg] * veg_con[iveg].Cv * 0.2);
	      rveg[band][14] += (rcan[iveg] * veg_con[iveg].Cv * 0.1);
	      rveg[band][13] += (rcan[iveg] * veg_con[iveg].Cv * 0.3);
	      rveg[band][12] += (rcan[iveg] * veg_con[iveg].Cv * 0.2);
              zov[band][14] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.1);
              zov[band][13] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.3);
              zov[band][12] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.2);
              displ[band][14] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.1);
              displ[band][13] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.3);
              displ[band][12] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.2);
	      fsum[14] += veg_con[iveg].Cv * 0.1;
	      fsum[13] += veg_con[iveg].Cv * 0.3;
	      fsum[12] += veg_con[iveg].Cv * 0.2;
	      break;
	    case 10: Tveg[band][17] += (Tcan[iveg] * veg_con[iveg].Cv);
	      rveg[band][17] += (rcan[iveg] * veg_con[iveg].Cv);
              zov[band][17] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv);
              displ[band][17] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv);
	      fsum[17] += veg_con[iveg].Cv;
	      break;
	    }

	  /* if(rec == 0 && nspinup == 0) */
	  /* if(rec == 0)

	    {

		  switch(veg_con[iveg].veg_class)
		    {
		    case 0: LAI[band][2] += (veg[iveg][band].LAI * \
					 veg_con[iveg].Cv);
		      break;
		    case 1: LAI[band][5] += (veg[iveg][band].LAI * \
					  veg_con[iveg].Cv);
		      break;
		    case 2: LAI[band][3] += (veg[iveg][band].LAI * \
					  veg_con[iveg].Cv);
		      break;
		    case 3: LAI[band][8] += (veg[iveg][band].LAI * \
					  veg_con[iveg].Cv);
		      break;
		    case 4: LAI[band][2] += (veg[iveg][band].LAI * \
					  veg_con[iveg].Cv * 0.5);
		      LAI[band][8] += (veg[iveg][band].LAI *	\
				    veg_con[iveg].Cv * 0.5);
		      break;
		    case 5: LAI[band][2] += (veg[iveg][band].LAI * \
					  veg_con[iveg].Cv * 0.4);
		      LAI[band][8] += (veg[iveg][band].LAI *	\
				    veg_con[iveg].Cv * 0.4);
		      LAI[band][9] += (veg[iveg][band].LAI *	\
				    veg_con[iveg].Cv * 0.1);
		      LAI[band][10] += (veg[iveg][band].LAI *	\
				     veg_con[iveg].Cv * 0.1);
		      break;
		    case 6: LAI[band][14] += (veg[iveg][band].LAI * \
					   veg_con[iveg].Cv * 0.1);
		      LAI[band][13] += (veg[iveg][band].LAI *	\
				     veg_con[iveg].Cv * 0.3);
		      LAI[band][12] += (veg[iveg][band].LAI *	\
				     veg_con[iveg].Cv * 0.2);
		      LAI[band][2] += (veg[iveg][band].LAI *	\
				    veg_con[iveg].Cv * 0.2);
		      LAI[band][8] += (veg[iveg][band].LAI *	\
				    veg_con[iveg].Cv * 0.2);
		      break;
		    case 7: LAI[band][9] += (veg[iveg][band].LAI * \
					  veg_con[iveg].Cv * 0.3);
		      LAI[band][10] += (veg[iveg][band].LAI *	\
				     veg_con[iveg].Cv * 0.3);
		      LAI[band][11] += (veg[iveg][band].LAI *	\
				     veg_con[iveg].Cv * 0.2);
		      LAI[band][2] += (veg[iveg][band].LAI *	\
				    veg_con[iveg].Cv * 0.1);
		      LAI[band][8] += (veg[iveg][band].LAI *	\
				    veg_con[iveg].Cv * 0.1);
		      break;
		    case 8: LAI[band][9] += (veg[iveg][band].LAI * \
					  veg_con[iveg].Cv * 0.2);
		      LAI[band][10] += (veg[iveg][band].LAI *	\
				     veg_con[iveg].Cv * 0.2);
		      LAI[band][11] += (veg[iveg][band].LAI *	\
				     veg_con[iveg].Cv * 0.2);
		      LAI[band][14] += (veg[iveg][band].LAI * 
				     veg_con[iveg].Cv * 0.05);
		      LAI[band][13] += (veg[iveg][band].LAI *	\
				     veg_con[iveg].Cv * 0.15);
		      LAI[band][12] += (veg[iveg][band].LAI *	\
				     veg_con[iveg].Cv * 0.1);
		      break;
		    case 9: LAI[band][14] += (veg[iveg][band].LAI * \
					   veg_con[iveg].Cv * 0.1);
		      LAI[band][13] += (veg[iveg][band].LAI *	\
				     veg_con[iveg].Cv * 0.3);
		      LAI[band][12] += (veg[iveg][band].LAI *	\
				     veg_con[iveg].Cv * 0.2);
		      break;
		    case 10: LAI[band][17] += (veg[iveg][band].LAI * \
					    veg_con[iveg].Cv);
		      break;
		    }

		    } */

	}

      fwet[band][0] = 0.0;
      decl[band] = cn[band].decl;
      fpi[band] = cn[band].fpi;
      fpg[band] = cn[band].fpg;
      annsum_counter[band] = cn[band].annsum_counter;
      cannsum_npp[band] = cn[band].cannsum_npp;
      cannavg_t2m[band] = cn[band].cannavg_t2m;
      for(nidx = 0; nidx < Nnodes; nidx++)
	{
	  watfc[band][nidx] = cn[band].watfc[nidx];
	}
      me[band] = cn[band].me;
      fire_prob[band] = cn[band].fire_prob;
      mean_fire_prob[band] = cn[band].mean_fire_prob;
      fireseasonl[band] = cn[band].fireseasonl;
      farea_burned[band] = cn[band].farea_burned;
      ann_farea_burned[band] = cn[band].ann_farea_burned;
      cwdc[band] = cn[band].cwdc;
      litrlabc[band] = cn[band].litr1c;
      litrcellc[band] = cn[band].litr2c;
      litrligc[band] = cn[band].litr3c;
      soilcfast[band] = cn[band].soil1c;
      soilcmid[band] = cn[band].soil2c;
      soilcslo1[band] = cn[band].soil3c;
      soilcslo2[band] = cn[band].soil4c;
      seedc[band] = cn[band].seedc;
      col_ctrunc[band] = cn[band].col_ctrunc;
      totlitc[band] = cn[band].totlitc;
      totcolc[band] = cn[band].totcolc;
      prod10c[band] = cn[band].prod10c;
      prod100c[band] = cn[band].prod100c;
      totsomc[band] = 0.0;
      cwdn[band] = cn[band].cwdn;
      litrlabn[band] = cn[band].litr1n;
      litrcelln[band] = cn[band].litr2n;
      litrlign[band] = cn[band].litr3n;
      soilnfast[band] = cn[band].soil1n;
      soilnmid[band] = cn[band].soil2n;
      soilnslo1[band] = cn[band].soil3n;
      soilnslo2[band] = cn[band].soil4n;
      soilminn[band] = cn[band].sminn;
      seedn[band] = cn[band].seedn;
      col_ntrunc[band] = cn[band].col_ntrunc;
      totcoln[band] = cn[band].totcoln;
      prod10n[band] = cn[band].prod10n;
      prod100n[band] = cn[band].prod100n;
      hr[band] = cn[band].hr;
      nee[band] = cn[band].nee;
      nep[band] = cn[band].nep;
      /*      if(rec == 0 && nspinup == 0  && options.CN_SPINUP && init_state == 0) */
      if(rec == 0 && init_state == 0)
	{
	  LAI[band][0] = 0.0;
	  dormant_flag[band][0] = 0.0;
	  days_active[band][0] = 0.0;
	  onset_flag[band][0] = 0.0;
	  onset_counter[band][0] = 0.0;
	  onset_gddflag[band][0] = 0.0;
	  onset_fdd[band][0] = 0.0;
	  onset_gdd[band][0] = 0.0;
	  onset_swi[band][0] = 0.0;
	  offset_flag[band][0] = 0.0;
          offset_counter[band][0] = 0.0;
	  offset_fdd[band][0] = 0.0;
	  offset_swi[band][0] = 0.0;
	  lgsf[band][0] = 0.0;
	  bglfr[band][0] = 0.0;
	  bgtr[band][0] = 0.0;
          dayl[band][0] = 0.0;
	  prev_dayl[band][0] = 0.0;
          annavg_t2m[band][0] = 0.0;
	  tempavg_t2m[band][0] = 0.0;
	  gpp[band][0] = 0.0;
	  availc[band][0] = 0.0;
	  xsmrpool_recover[band][0] = 0.0;
	  alloc_pnow[band][0] = 0.0;
	  c_allometry[band][0] = 0.0;
	  n_allometry[band][0] = 0.0;
	  plant_ndemand[band][0] = 0.0;
	  tempsum_potential_gpp[band][0] = 0.0;
	  annsum_potential_gpp[band][0] = 0.0;
	  tempmax_retransn[band][0] = 0.0;
	  annmax_retransn[band][0] = 0.0;
	  avail_retransn[band][0] = 0.0;
	  plant_nalloc[band][0] = 0.0;
	  plant_calloc[band][0] = 0.0;
	  excess_cflux[band][0] = 0.0;
	  downreg[band][0] = 0.0;
	  prev_leafc_to_litter[band][0] = 0.0;
	  prev_frootc_to_litter[band][0] = 0.0;
	  tempsum_npp[band][0] = 0.0;
	  annsum_npp[band][0] = 0.0;
	  gpp2[band][0] = 0.0;
	  gpp[band][0] = 0.0;
	  npp[band][0] = 0.0;
	  ar[band][0] = 0.0;
	  leafc[band][0] = 0.0;
	  leafc_storage[band][0] = 0.0;
	  leafc_xfer[band][0] = 0.0;
	  frootc[band][0] = 0.0;
	  frootc_storage[band][0] = 0.0;
	  frootc_xfer[band][0] = 0.0;
	  livestemc[band][0] = 0.0;
	  livestemc_xfer[band][0] = 0.0;
	  deadstemc[band][0] = 0.0;
	  deadstemc_storage[band][0] = 0.0;
	  deadstemc_xfer[band][0] = 0.0;
	  livecrootc[band][0] = 0.0;
	  livecrootc_storage[band][0] = 0.0;
	  livecrootc_xfer[band][0] = 0.0;
	  deadcrootc[band][0] = 0.0;
	  deadcrootc_storage[band][0] = 0.0;
	  deadcrootc_xfer[band][0] = 0.0;
	  gresp_storage[band][0] = 0.0;
	  gresp_xfer[band][0] = 0.0;
	  cpool[band][0] = 0.0;
	  xsmrpool[band][0] = 0.0;
	  pft_ctrunc[band][0] = 0.0;
	  totvegc[band][0] = 0.0;
	  leafn[band][0] = 0.0;
	  leafn_storage[band][0] = 0.0;
	  leafn_xfer[band][0] = 0.0;
	  frootn[band][0] = 0.0;
	  frootn_storage[band][0] = 0.0;
	  frootn_xfer[band][0] = 0.0;
	  livestemn[band][0] = 0.0;
	  livestemn_storage[band][0] = 0.0;
	  livestemn_xfer[band][0] = 0.0;
	  deadstemn[band][0] = 0.0;
	  deadstemn_storage[band][0] = 0.0;
	  deadstemn_xfer[band][0] = 0.0;
	  livecrootn[band][0] = 0.0;
	  livecrootn_storage[band][0] = 0.0;
	  livecrootn_xfer[band][0] = 0.0;
	  retransn[band][0] = 0.0;
	  npool[band][0] = 0.0;
	  pft_ntrunc[band][0] = 0.0;
	}
      else
	{
	  LAI[band][0] = cn[band].LAI[0];
	  dormant_flag[band][0] = cn[band].dormant_flag[0];
	  days_active[band][0] = cn[band].days_active[0];
	  onset_flag[band][0] = cn[band].onset_flag[0];
	  onset_counter[band][0] = cn[band].onset_counter[0];
	  onset_gddflag[band][0] = cn[band].onset_gddflag[0];
	  onset_fdd[band][0] = cn[band].onset_fdd[0];
	  onset_gdd[band][0] = cn[band].onset_gdd[0];
	  onset_swi[band][0] = cn[band].onset_swi[0];
	  offset_flag[band][0] = cn[band].offset_flag[0];
          offset_counter[band][0] = cn[band].offset_counter[0];
	  offset_fdd[band][0] = cn[band].offset_fdd[0];
	  offset_swi[band][0] = cn[band].offset_swi[0];
	  lgsf[band][0] = cn[band].lgsf[0];
	  bglfr[band][0] = cn[band].bglfr[0];
	  bgtr[band][0] = cn[band].bgtr[0];
          dayl[band][0] = cn[band].dayl[0];
	  prev_dayl[band][0] = cn[band].prev_dayl[0];
          annavg_t2m[band][0] = cn[band].annavg_t2m[0];
	  tempavg_t2m[band][0] = cn[band].tempavg_t2m[0];
	  gpp[band][0] = cn[band].gpp[0];
	  availc[band][0] = cn[band].availc[0];
	  xsmrpool_recover[band][0] = cn[band].xsmrpool_recover[0];
	  alloc_pnow[band][0] = cn[band].alloc_pnow[0];
	  c_allometry[band][0] = cn[band].c_allometry[0];
	  n_allometry[band][0] = cn[band].n_allometry[0];
	  plant_ndemand[band][0] = cn[band].plant_ndemand[0];
	  tempsum_potential_gpp[band][0] = cn[band].tempsum_potential_gpp[0];
	  annsum_potential_gpp[band][0] = cn[band].annsum_potential_gpp[0];
	  tempmax_retransn[band][0] = cn[band].tempmax_retransn[0];
	  annmax_retransn[band][0] = cn[band].annmax_retransn[0];
	  avail_retransn[band][0] = cn[band].avail_retransn[0];
	  plant_nalloc[band][0] = cn[band].plant_nalloc[0];
	  plant_calloc[band][0] = cn[band].plant_calloc[0];
	  excess_cflux[band][0] = cn[band].excess_cflux[0];
	  downreg[band][0] = cn[band].downreg[0];
	  prev_leafc_to_litter[band][0] = cn[band].prev_leafc_to_litter[0];
	  prev_frootc_to_litter[band][0] = cn[band].prev_frootc_to_litter[0];
	  tempsum_npp[band][0] = cn[band].tempsum_npp[0];
	  annsum_npp[band][0] = cn[band].annsum_npp[0];
	  gpp2[band][0] = cn[band].gpp2[0];
	  gpp[band][0] = cn[band].gpp[0];
	  npp[band][0] = cn[band].npp[0];
	  ar[band][0] = cn[band].ar[0];
	  leafc[band][0] = cn[band].leafc[0];
	  leafc_storage[band][0] = cn[band].leafc_storage[0];
	  leafc_xfer[band][0] = cn[band].leafc_xfer[0];
	  frootc[band][0] = cn[band].frootc[0];
	  frootc_storage[band][0] = cn[band].frootc_storage[0];
	  frootc_xfer[band][0] = cn[band].frootc_xfer[0];
	  livestemc[band][0] = cn[band].livestemc[0];
	  livestemc_xfer[band][0] = cn[band].livestemc_xfer[0];
	  deadstemc[band][0] = cn[band].deadstemc[0];
	  deadstemc_storage[band][0] = cn[band].deadstemc_storage[0];
	  deadstemc_xfer[band][0] = cn[band].deadstemc_xfer[0];
	  livecrootc[band][0] = cn[band].livecrootc[0];
	  livecrootc_storage[band][0] = cn[band].livecrootc_storage[0];
	  livecrootc_xfer[band][0] = cn[band].livecrootc_xfer[0];
	  deadcrootc[band][0] = cn[band].deadcrootc[0];
	  deadcrootc_storage[band][0] = cn[band].deadcrootc_storage[0];
	  deadcrootc_xfer[band][0] = cn[band].deadcrootc_xfer[0];
	  gresp_storage[band][0] = cn[band].gresp_storage[0];
	  gresp_xfer[band][0] = cn[band].gresp_xfer[0];
	  cpool[band][0] = cn[band].cpool[0];
	  xsmrpool[band][0] = cn[band].xsmrpool[0];
	  pft_ctrunc[band][0] = cn[band].pft_ctrunc[0];
	  totvegc[band][0] = cn[band].totvegc[0];
	  leafn[band][0] = cn[band].leafn[0];
	  leafn_storage[band][0] = cn[band].leafn_storage[0];
	  leafn_xfer[band][0] = cn[band].leafn_xfer[0];
	  frootn[band][0] = cn[band].frootn[0];
	  frootn_storage[band][0] = cn[band].frootn_storage[0];
	  frootn_xfer[band][0] = cn[band].frootn_xfer[0];
	  livestemn[band][0] = cn[band].livestemn[0];
	  livestemn_storage[band][0] = cn[band].livestemn_storage[0];
	  livestemn_xfer[band][0] = cn[band].livestemn_xfer[0];
	  deadstemn[band][0] = cn[band].deadstemn[0];
	  deadstemn_storage[band][0] = cn[band].deadstemn_storage[0];
	  deadstemn_xfer[band][0] = cn[band].deadstemn_xfer[0];
	  livecrootn[band][0] = cn[band].livecrootn[0];
	  livecrootn_storage[band][0] = cn[band].livecrootn_storage[0];
	  livecrootn_xfer[band][0] = cn[band].livecrootn_xfer[0];
	  retransn[band][0] = cn[band].retransn[0];
	  npool[band][0] = cn[band].npool[0];
	  pft_ntrunc[band][0] = cn[band].pft_ntrunc[0];
	}
      for(iveg = 1; iveg <= MAX_PFT; iveg++)
	{
	  /*	  if(rec == 0 && nspinup == 0  && options.CN_SPINUP && init_state == 0) */
	  if(rec == 0 && init_state == 0)
	    {
	      LAI[band][iveg] = 0.0;
	      dormant_flag[band][iveg] = 0.0;
	      days_active[band][iveg] = 0.0;
	      onset_flag[band][iveg] = 0.0;
	      onset_counter[band][iveg] = 0.0;
	      onset_gddflag[band][iveg] = 0.0;
	      onset_fdd[band][iveg] = 0.0;
	      onset_gdd[band][iveg] = 0.0;
	      onset_swi[band][iveg] = 0.0;
	      offset_flag[band][iveg] = 0.0;
	      offset_counter[band][iveg] = 0.0;
	      offset_fdd[band][iveg] = 0.0;
	      offset_swi[band][iveg] = 0.0;
	      lgsf[band][iveg] = 0.0;
	      bglfr[band][iveg] = 0.0;
	      bgtr[band][iveg] = 0.0;
	      dayl[band][iveg] = 0.0;
	      prev_dayl[band][iveg] = 0.0;
	      annavg_t2m[band][iveg] = 0.0;
	      tempavg_t2m[band][iveg] = 0.0;
	      gpp[band][iveg] = 0.0;
	      availc[band][iveg] = 0.0;
	      xsmrpool_recover[band][iveg] = 0.0;
	      alloc_pnow[band][iveg] = 0.0;
	      c_allometry[band][iveg] = 0.0;
	      n_allometry[band][iveg] = 0.0;
	      plant_ndemand[band][iveg] = 0.0;
	      tempsum_potential_gpp[band][iveg] = 0.0;
	      annsum_potential_gpp[band][iveg] = 0.0;
	      tempmax_retransn[band][iveg] = 0.0;
	      annmax_retransn[band][iveg] = 0.0;
	      avail_retransn[band][iveg] = 0.0;
	      plant_nalloc[band][iveg] = 0.0;
	      plant_calloc[band][iveg] = 0.0;
	      excess_cflux[band][iveg] = 0.0;
	      downreg[band][iveg] = 0.0;
	      prev_leafc_to_litter[band][iveg] = 0.0;
	      prev_frootc_to_litter[band][iveg] = 0.0;
	      tempsum_npp[band][iveg] = 0.0;
	      annsum_npp[band][iveg] = 0.0;
	      gpp2[band][iveg] = 0.0;
	      gpp[band][iveg] = 0.0;
	      npp[band][iveg] = 0.0;
	      ar[band][iveg] = 0.0;
	      leafc[band][iveg] = 0.0;
	      leafc_storage[band][iveg] = 0.0;
	      leafc_xfer[band][iveg] = 0.0;
	      frootc[band][iveg] = 0.0;
	      frootc_storage[band][iveg] = 0.0;
	      frootc_xfer[band][iveg] = 0.0;
	      livestemc[band][iveg] = 0.0;
	      livestemc_xfer[band][iveg] = 0.0;
	      deadstemc[band][iveg] = 0.0;
	      deadstemc_storage[band][iveg] = 0.0;
	      deadstemc_xfer[band][iveg] = 0.0;
	      livecrootc[band][iveg] = 0.0;
	      livecrootc_storage[band][iveg] = 0.0;
	      livecrootc_xfer[band][iveg] = 0.0;
	      deadcrootc[band][iveg] = 0.0;
	      deadcrootc_storage[band][iveg] = 0.0;
	      deadcrootc_xfer[band][iveg] = 0.0;
	      gresp_storage[band][iveg] = 0.0;
	      gresp_xfer[band][iveg] = 0.0;
	      cpool[band][iveg] = 0.0;
	      xsmrpool[band][iveg] = 0.0;
	      pft_ctrunc[band][iveg] = 0.0;
	      totvegc[band][iveg] = 0.0;
	      leafn[band][iveg] = 0.0;
	      leafn_storage[band][iveg] = 0.0;
	      leafn_xfer[band][iveg] = 0.0;
	      frootn[band][iveg] = 0.0;
	      frootn_storage[band][iveg] = 0.0;
	      frootn_xfer[band][iveg] = 0.0;
	      livestemn[band][iveg] = 0.0;
	      livestemn_storage[band][iveg] = 0.0;
	      livestemn_xfer[band][iveg] = 0.0;
	      deadstemn[band][iveg] = 0.0;
	      deadstemn_storage[band][iveg] = 0.0;
	      deadstemn_xfer[band][iveg] = 0.0;
	      livecrootn[band][iveg] = 0.0;
	      livecrootn_storage[band][iveg] = 0.0;
	      livecrootn_xfer[band][iveg] = 0.0;
	      retransn[band][iveg] = 0.0;
	      npool[band][iveg] = 0.0;
	      pft_ntrunc[band][iveg] = 0.0;
	    }
	  else
	    {
	      LAI[band][iveg] = cn[band].LAI[iveg];
	      days_active[band][iveg] = cn[band].days_active[iveg];
	      onset_flag[band][iveg] = cn[band].onset_flag[iveg];
	      onset_counter[band][iveg] = cn[band].onset_counter[iveg];
	      onset_gddflag[band][iveg] = cn[band].onset_gddflag[iveg];
	      onset_fdd[band][iveg] = cn[band].onset_fdd[iveg];
	      onset_gdd[band][iveg] = cn[band].onset_gdd[iveg];
	      onset_swi[band][iveg] = cn[band].onset_swi[iveg];
	      offset_flag[band][iveg] = cn[band].offset_flag[iveg];
	      offset_counter[band][iveg] = cn[band].offset_counter[iveg];
	      offset_fdd[band][iveg] = cn[band].offset_fdd[iveg];
	      offset_swi[band][iveg] = cn[band].offset_swi[iveg];
	      lgsf[band][iveg] = cn[band].lgsf[iveg];
	      bglfr[band][iveg] = cn[band].bglfr[iveg];
	      bgtr[band][iveg] = cn[band].bgtr[iveg];
	      dayl[band][iveg] = cn[band].dayl[iveg];
	      prev_dayl[band][iveg] = cn[band].prev_dayl[iveg];
	      annavg_t2m[band][iveg] = cn[band].annavg_t2m[iveg];
	      tempavg_t2m[band][iveg] = cn[band].tempavg_t2m[iveg];
	      gpp[band][iveg] = cn[band].gpp[iveg];
	      availc[band][iveg] = cn[band].availc[iveg];
	      xsmrpool_recover[band][iveg] = cn[band].xsmrpool_recover[iveg];
	      alloc_pnow[band][iveg] = cn[band].alloc_pnow[iveg];
	      c_allometry[band][iveg] = cn[band].c_allometry[iveg];
	      n_allometry[band][iveg] = cn[band].n_allometry[iveg];
	      plant_ndemand[band][iveg] = cn[band].plant_ndemand[iveg];
	      tempsum_potential_gpp[band][iveg] = \
		cn[band].tempsum_potential_gpp[iveg];
	      annsum_potential_gpp[band][iveg] = \
		cn[band].annsum_potential_gpp[iveg];
	      tempmax_retransn[band][iveg] = cn[band].tempmax_retransn[iveg];
	      annmax_retransn[band][iveg] = cn[band].annmax_retransn[iveg];
	      avail_retransn[band][iveg] = cn[band].avail_retransn[iveg];
	      plant_nalloc[band][iveg] = cn[band].plant_nalloc[iveg];
	      plant_calloc[band][iveg] = cn[band].plant_calloc[iveg];
	      excess_cflux[band][iveg] = cn[band].excess_cflux[iveg];
	      downreg[band][iveg] = cn[band].downreg[iveg];
	      prev_leafc_to_litter[band][iveg] = \
		cn[band].prev_leafc_to_litter[iveg];
	      prev_frootc_to_litter[band][iveg] = \
		cn[band].prev_frootc_to_litter[iveg];
	      tempsum_npp[band][iveg] = cn[band].tempsum_npp[iveg];
	      annsum_npp[band][iveg] = cn[band].annsum_npp[iveg];
	      gpp2[band][iveg] = cn[band].gpp2[iveg];
	      gpp[band][iveg] = cn[band].gpp[iveg];
	      npp[band][iveg] = cn[band].npp[iveg];
	      ar[band][iveg] = cn[band].ar[iveg];
	      leafc[band][iveg] = cn[band].leafc[iveg];
	      leafc_storage[band][iveg] = cn[band].leafc_storage[iveg];
	      leafc_xfer[band][iveg] = cn[band].leafc_xfer[iveg];
	      frootc[band][iveg] = cn[band].frootc[iveg];
	      frootc_storage[band][iveg] = cn[band].frootc_storage[iveg];
	      frootc_xfer[band][iveg] = cn[band].frootc_xfer[iveg];
	      livestemc[band][iveg] = cn[band].livestemc[iveg];
	      livestemc_xfer[band][iveg] = cn[band].livestemc_xfer[iveg];
	      deadstemc[band][iveg] = cn[band].deadstemc[iveg];
	      deadstemc_storage[band][iveg] = cn[band].deadstemc_storage[iveg];
	      deadstemc_xfer[band][iveg] = cn[band].deadstemc_xfer[iveg];
	      livecrootc[band][iveg] = cn[band].livecrootc[iveg];
	      livecrootc_storage[band][iveg] = \
		cn[band].livecrootc_storage[iveg];
	      livecrootc_xfer[band][iveg] = cn[band].livecrootc_xfer[iveg];
	      deadcrootc[band][iveg] = cn[band].deadcrootc[iveg];
	      deadcrootc_storage[band][iveg] = \
		cn[band].deadcrootc_storage[iveg];
	      deadcrootc_xfer[band][iveg] = cn[band].deadcrootc_xfer[iveg];
	      gresp_storage[band][iveg] = cn[band].gresp_storage[iveg];
	      gresp_xfer[band][iveg] = cn[band].gresp_xfer[iveg];
	      cpool[band][iveg] = cn[band].cpool[iveg];
	      xsmrpool[band][iveg] = cn[band].xsmrpool[iveg];
	      pft_ctrunc[band][iveg] = cn[band].pft_ctrunc[iveg];
	      totvegc[band][iveg] = cn[band].totvegc[iveg];
	      leafn[band][iveg] = cn[band].leafn[iveg];
	      leafn_storage[band][iveg] = cn[band].leafn_storage[iveg];
	      leafn_xfer[band][iveg] = cn[band].leafn_xfer[iveg];
	      frootn[band][iveg] = cn[band].frootn[iveg];
	      frootn_storage[band][iveg] = cn[band].frootn_storage[iveg];
	      frootn_xfer[band][iveg] = cn[band].frootn_xfer[iveg];
	      livestemn[band][iveg] = cn[band].livestemn[iveg];
	      livestemn_storage[band][iveg] = cn[band].livestemn_storage[iveg];
	      livestemn_xfer[band][iveg] = cn[band].livestemn_xfer[iveg];
	      deadstemn[band][iveg] = cn[band].deadstemn[iveg];
	      deadstemn_storage[band][iveg] = cn[band].deadstemn_storage[iveg];
	      deadstemn_xfer[band][iveg] = cn[band].deadstemn_xfer[iveg];
	      livecrootn[band][iveg] = cn[band].livecrootn[iveg];
	      livecrootn_storage[band][iveg] = \
		cn[band].livecrootn_storage[iveg];
	      livecrootn_xfer[band][iveg] = cn[band].livecrootn_xfer[iveg];
	      retransn[band][iveg] = cn[band].retransn[iveg];
	      npool[band][iveg] = cn[band].npool[iveg];
	      pft_ntrunc[band][iveg] = cn[band].pft_ntrunc[iveg];
	    }

	  if(fsum[iveg] != 0.0)
	    {
	      Tveg[band][iveg] /= fsum[iveg];
	      rveg[band][iveg] /= fsum[iveg];
	      zov[band][iveg] /= fsum[iveg];
	      displ[band][iveg] /= fsum[iveg];
	    }

	    if(precip > 0.0)
	      fwet[band][iveg] = 1.0;
	    else
	      fwet[band][iveg] = 0.0;
	}

    }

  for(iveg = 0; iveg < Nveg; iveg++)
    {
      Tcan[iveg] /= Nbands;
      rcan[iveg] /= Nbands;
    }

  /* Snow depth */
  for(band = 0; band < Nbands; band++)
    {

      snowdep[band] = 0.0;

      for(iveg = 0; iveg < Nveg; iveg++)
	{
	      snowdep[band] += snow[iveg][band].depth;
	}

      snowdep[band] /= Nveg;

      if(snowdep[band] > 1.e-3)
	{
	  dz[band][0] = 1.0e-3;
	  dz[band][1] = snowdep[band] - 1.0e-3;

	  z[band][1] = -0.5 * (snowdep[band] - 1.0e-3);
	  z[band][0] = -5.0e-4 - dz[band][1];
	}
      else
	{
	  z[band][0] = 0.0;
	  z[band][1] = 0.0;

	  dz[band][0] = 0.0;
	  dz[band][1] = 0.0;
	}

    }

  /* Distribute root fraction across snow bands */
  for(band = 0; band < Nbands; band++)
    {
      for(iveg = 0; iveg < Npft; iveg++)
	{

	  Lsum = 0.0;
	  lidx = 0;
	  PAST_BOTTOM = FALSE;

	  for(nidx = 0; nidx < Nnodes; nidx++)
	    {

	      if(soil_con->Zsum_node[nidx] == Lsum + soil_con->depth[lidx] && \
		 nidx != 0 && lidx != Nlayers - 1)
		  rootfr[band][iveg][nidx] = (veg_con->root[lidx] + \
					      veg_con->root[lidx+1]) / 2.0;
	      else
		rootfr[band][iveg][nidx] = veg_con->root[lidx];

	      if(soil_con->Zsum_node[nidx] > Lsum + soil_con->depth[lidx] && \
		 !PAST_BOTTOM)
		{
		  Lsum += soil_con->depth[lidx];
		  lidx++;
		  if(lidx == Nlayers)
		    {
		      PAST_BOTTOM = TRUE;
		      lidx = Nlayers - 1;
		    }
		}
	    }
	}

  /* Place saturated matric potential into array and calculate soil matric 
     potential for each layer (Need to convert from cm to MPa) */

      for(lidx = 0; lidx < Nlayers; lidx++)
	{
	  theta[lidx] = 0.0;
	}

      for(iveg = 0; iveg < Nveg; iveg++)
	{
	  for(lidx = 0; lidx < Nlayers; lidx++)
	    {
	      theta[lidx] += cell[iveg][band].layer[lidx].moist;
	    }
	}

      for(lidx = 0; lidx < Nlayers; lidx++)
	{
	  theta[lidx] /= (Nveg - 1);
	  theta[lidx] /= (soil_con->depth[lidx] * 1000.0);
	}

      for(lidx = 0; lidx < Nlayers; lidx++)
	{

	  psand = soil_con->quartz[lidx] * soil_con->bulk_dens_min[lidx] / \
	    soil_con->soil_dens_min[lidx] * 100.0;
          pclay = (1.0 - soil_con->quartz[lidx] - soil_con->porosity[lidx]) * \
	    soil_con->bulk_dens_min[lidx] / soil_con->soil_dens_min[lidx] * \
	    100.0;
	  if(pclay < 0.0)
	    pclay = 0.0;

	  bsw[0][lidx] = 2.91 + 0.159 * pclay;
          sucsat[0][lidx] = 10.0 * (10.0 * (1.88 - 0.0131 * psand));
	  watsat[lidx] = 0.489 - 0.00126 * psand;
          ksat = 0.0070556 * pow(10., -0.884 + 0.0153 * psand);
	  fsat = theta[band] / soil_con->porosity[lidx];
	  if(fsat < 0.001)
	    fsat = 0.001;
	  soisuc[0][lidx] = sucsat[0][lidx] * pow(fsat, bsw[0][lidx]);
	  if(soisuc[0][lidx] > 15.0)
	    soisuc[0][lidx] = 15.0;
	  if(soisuc[0][lidx] < 0.0)
	    soisuc[0][lidx] = 0.0;
	  watfrac[lidx] = watsat[lidx] * pow(0.1 / ksat / 86400., 1.0 / \
					     (2.0 * bsw[0][lidx] + 3.0));

	  bsw[1][lidx] = -(3.10 + 0.157 * pclay - 0.003 * psand);
          sucsat[1][lidx] = -1.0 * (exp((1.54 - 0.0095 * psand + \
					       0.0063 *	(100.0	- psand \
							 - pclay )) *	\
				     log(10.0)) * 9.8e-5);
	  watsat[lidx] = (50.5 - 0.142 * psand - 0.037 * pclay) / 100.0;
	  fsat = theta[lidx] / soil_con->porosity[lidx];
	  if(fsat < 0.001)
	    fsat = 0.001;
	  soisuc[1][lidx] = sucsat[1][lidx] * pow(fsat, bsw[1][lidx]);
	  if(soisuc[1][lidx] < -15.0)
	    soisuc[1][lidx] = -15.0;
	  if(soisuc[1][lidx] > 0.0)
	    soisuc[1][lidx] = 0.0;

	}

      Lsum = 0.0;
      lidx = 0;
      PAST_BOTTOM = FALSE;

      for(nidx = 0; nidx < Nnodes; nidx++)
	{

	  for(band = 0; band < Nbands; band++)
	    {

	      for(i = 0; i < 2; i++)

		{

		  if(soil_con->Zsum_node[nidx] == Lsum + \
		     soil_con->depth[lidx] && nidx != 0 && lidx != Nlayers - 1)
		    {
		      coeff[i][band][nidx] = (bsw[i][lidx] + bsw[i][lidx+1]) \
			/ 2.0;
		      psisat[i][band][nidx] = (sucsat[i][lidx] +	\
					   sucsat[i][lidx+1]) / 2.0;
		      soipsi[i][band][nidx] = (soisuc[i][lidx] +	\
					   soisuc[i][lidx+1]) / 2.0;
		    }
		  else
		    {
		      coeff[i][band][nidx] = bsw[i][lidx];
		      psisat[i][band][nidx] = sucsat[i][lidx];
		      soipsi[i][band][nidx] = soisuc[i][lidx];
		    }

		}

	      if(soil_con->Zsum_node[nidx] == Lsum + soil_con->depth[lidx] \
		 && nidx != 0 && lidx != Nlayers - 1)
		{
		  watfc[band][nidx] = (watfrac[lidx] + watfrac[lidx+1]) / 2.0;
		}
	      else
		{
		  watfc[band][nidx] = watfrac[lidx];
		}

	      if(soil_con->Zsum_node[nidx] > Lsum + soil_con->depth[lidx] && \
		 !PAST_BOTTOM)
		{
		  Lsum += soil_con->depth[lidx];
		  lidx++;
		  if(lidx == Nlayers)
		    {
		      PAST_BOTTOM = TRUE;
		      lidx = Nlayers - 1;
		    }
		}

	    }

	}
    }

  /* This will need to be adapted when implemented into RASM or when VIC is */
  /* run regionally. */

  lbg = 1;
  ubg = 1;
  lbc = 1;
  ubc = Nbands;
  lbp = 1;
  ubp = Npft;
  num_soilc = Nbands;
  num_soilp = Npft;

  vic2clmtype_(&Nnodes, &rec, &Nrecs, &adspinup, &year, &month, \
	       &day, &secs, &jday, &yrnxt, &jdynxt, &dt, &lat, &lon, &lbg, \
	       &ubg, &lbc, &ubc, &lbp, &ubp, &num_soilc, &num_soilp, &psfc, \
	       &Tair, &vp, &vpd, &lwrad, &swrad, swrd, swri, alb, z, dz, \
	       baseflow, moist, ice, t_soisno, t2m, Tveg, snowdep, fwet, \
	       rootfr, psisat, soipsi, coeff, rveg, &zo, &zos, zov, displ, \
	       LAI, soilcfast, soilcmid, soilcslo1, soilcslo2, litrlabc, \
	       litrcellc, litrligc, cwdc, leafc, frootc, livestemc, \
	       deadstemc, livecrootc, deadcrootc, woodc, soilnfast, soilnmid, \
	       soilnslo1, soilnslo2, soilminn, litrlabn, litrcelln, litrlign, \
	       cwdn, leafn, frootn, livestemn, deadstemn, livecrootn, \
	       deadcrootn, totvegc, totlitc, totsomc, gpp2, gpp, npp, darkr, mr, gr, ar, hr, \
	       lithr, nee, nep, dormant_flag, days_active, onset_flag, onset_counter, \
	       onset_gddflag, onset_fdd, onset_gdd, onset_swi, offset_flag, \
	       offset_counter, offset_fdd, offset_swi, lgsf, bglfr,	\
	       bgtr, dayl, prev_dayl, annavg_t2m, tempavg_t2m, availc, \
	       xsmrpool_recover, alloc_pnow, c_allometry, n_allometry, \
	       plant_ndemand, tempsum_potential_gpp, annsum_potential_gpp, \
	       tempmax_retransn, annmax_retransn, avail_retransn, \
	       plant_nalloc, plant_calloc, excess_cflux, downreg, \
	       prev_leafc_to_litter, prev_frootc_to_litter, tempsum_npp, \
	       annsum_npp, leafc_storage, leafc_xfer, frootc_storage, \
	       frootc_xfer, livestemc_storage, livestemc_xfer, \
	       deadstemc_storage, deadstemc_xfer, livecrootc_storage, \
	       livecrootc_xfer, deadcrootc_storage, deadcrootc_xfer, \
	       gresp_storage, gresp_xfer, cpool, xsmrpool, pft_ctrunc, \
	       leafn_storage, leafn_xfer, frootn_storage, frootn_xfer, \
	       livestemn_storage, livestemn_xfer, deadstemn_storage, \
	       deadstemn_xfer, livecrootn_storage, livecrootn_xfer, \
	       deadcrootn_storage, deadcrootn_xfer, retransn, npool, \
	       pft_ntrunc, decl, fpi, fpg, annsum_counter, cannsum_npp, \
	       cannavg_t2m, watfc, me, fire_prob, mean_fire_prob, \
	       fireseasonl, farea_burned, ann_farea_burned, seedc, \
	       col_ctrunc, totcolc, prod10c, prod100c, seedn, \
	       col_ntrunc, totcoln, prod10n, prod100n, litfall, fpsn, ci, rc, \
	       apar, &init_state);

  for(band = 0; band < Nbands; band++)
    {
      for(iveg = 0; iveg < MAX_PFT; iveg++)
	{
	  cn[band].LAI[iveg] = LAI[band][iveg];
	  cn[band].dormant_flag[iveg] = dormant_flag[band][iveg];
	  cn[band].days_active[iveg] = days_active[band][iveg];
	  cn[band].onset_flag[iveg] = onset_flag[band][iveg];
	  cn[band].onset_counter[iveg] = onset_counter[band][iveg];
	  cn[band].onset_gddflag[iveg] = onset_gddflag[band][iveg];
	  cn[band].onset_fdd[iveg] = onset_fdd[band][iveg];
	  cn[band].onset_gdd[iveg] = onset_gdd[band][iveg];
	  cn[band].onset_swi[iveg] = onset_swi[band][iveg];
	  cn[band].offset_flag[iveg] = offset_flag[band][iveg];
	  cn[band].offset_counter[iveg] = offset_counter[band][iveg];
	  cn[band].offset_fdd[iveg] = offset_fdd[band][iveg];
	  cn[band].offset_swi[iveg] = offset_swi[band][iveg];
	  cn[band].lgsf[iveg] = lgsf[band][iveg];
	  cn[band].bglfr[iveg] = bglfr[band][iveg];
	  cn[band].bgtr[iveg] = bgtr[band][iveg];
	  cn[band].dayl[iveg] = dayl[band][iveg];
	  cn[band].prev_dayl[iveg] = prev_dayl[band][iveg];
	  cn[band].annavg_t2m[iveg] = annavg_t2m[band][iveg];
	  cn[band].tempavg_t2m[iveg] = tempavg_t2m[band][iveg];
	  cn[band].availc[iveg] = availc[band][iveg];
	  cn[band].xsmrpool_recover[iveg] = xsmrpool_recover[band][iveg];
	  cn[band].alloc_pnow[iveg] = alloc_pnow[band][iveg];
	  cn[band].c_allometry[iveg] = c_allometry[band][iveg];
	  cn[band].n_allometry[iveg] = n_allometry[band][iveg];
	  cn[band].plant_ndemand[iveg] = plant_ndemand[band][iveg];
	  cn[band].tempsum_potential_gpp[iveg] = tempsum_potential_gpp[band][iveg];
	  cn[band].annsum_potential_gpp[iveg] = annsum_potential_gpp[band][iveg];
	  cn[band].tempmax_retransn[iveg] = tempmax_retransn[band][iveg];
	  cn[band].annmax_retransn[iveg] = annmax_retransn[band][iveg];
	  cn[band].avail_retransn[iveg] = avail_retransn[band][iveg];
	  cn[band].plant_nalloc[iveg] = plant_nalloc[band][iveg];
	  cn[band].plant_calloc[iveg] = plant_calloc[band][iveg];
	  cn[band].excess_cflux[iveg] = excess_cflux[band][iveg];
	  cn[band].downreg[iveg] = downreg[band][iveg];
	  cn[band].prev_leafc_to_litter[iveg] = prev_leafc_to_litter[band][iveg];
	  cn[band].prev_frootc_to_litter[iveg] = prev_frootc_to_litter[band][iveg];
	  cn[band].tempsum_npp[iveg] = tempsum_npp[band][iveg];
	  cn[band].annsum_npp[iveg] = annsum_npp[band][iveg];
	  cn[band].gpp2[iveg] = gpp2[band][iveg];
	  cn[band].gpp[iveg] = gpp[band][iveg];
	  cn[band].npp[iveg] = npp[band][iveg];
	  cn[band].ar[iveg] = ar[band][iveg];
	  cn[band].leafc[iveg] = leafc[band][iveg];
	  cn[band].leafc_storage[iveg] = leafc_storage[band][iveg];
	  cn[band].leafc_xfer[iveg] = leafc_xfer[band][iveg];
	  cn[band].frootc[iveg] = frootc[band][iveg];
	  cn[band].frootc_storage[iveg] = frootc_storage[band][iveg];
	  cn[band].frootc_xfer[iveg] = frootc_xfer[band][iveg];
	  cn[band].livestemc[iveg] = livestemc[band][iveg];
	  cn[band].livestemc_storage[iveg] = livestemc_storage[band][iveg];
	  cn[band].livestemc_xfer[iveg] = livestemc_xfer[band][iveg];
	  cn[band].deadstemc[iveg] = deadstemc[band][iveg];
	  cn[band].deadstemc_storage[iveg] = deadstemc_storage[band][iveg];
	  cn[band].deadstemc_xfer[iveg] = deadstemc_xfer[band][iveg];
	  cn[band].livecrootc[iveg] = livecrootc[band][iveg];
	  cn[band].livecrootc_storage[iveg] = livecrootc_storage[band][iveg];
	  cn[band].livecrootc_xfer[iveg] = livecrootc_xfer[band][iveg];
	  cn[band].deadcrootc[iveg] = deadcrootc[band][iveg];
	  cn[band].deadcrootc_storage[iveg] = deadcrootc_storage[band][iveg];
	  cn[band].deadcrootc_xfer[iveg] = deadcrootc_xfer[band][iveg];
	  cn[band].gresp_storage[iveg] = gresp_storage[band][iveg];
	  cn[band].gresp_xfer[iveg] = gresp_xfer[band][iveg];
	  cn[band].cpool[iveg] = cpool[band][iveg];
	  cn[band].xsmrpool[iveg] = xsmrpool[band][iveg];
	  cn[band].pft_ctrunc[iveg] = pft_ctrunc[band][iveg];
	  cn[band].totvegc[iveg] = totvegc[band][iveg];
	  cn[band].woodc[iveg] = woodc[band][iveg];
	  cn[band].leafn[iveg] = leafn[band][iveg];
	  cn[band].leafn_storage[iveg] = leafn_storage[band][iveg];
	  cn[band].leafn_xfer[iveg] = leafn_xfer[band][iveg];
	  cn[band].frootn[iveg] = frootn[band][iveg];
	  cn[band].frootn_storage[iveg] = frootn_storage[band][iveg];
	  cn[band].frootn_xfer[iveg] = frootn_xfer[band][iveg];
	  cn[band].livestemn[iveg] = livestemn[band][iveg];
	  cn[band].livestemn_storage[iveg] = livestemn_storage[band][iveg];
	  cn[band].livestemn_xfer[iveg] = livestemn_xfer[band][iveg];
	  cn[band].deadstemn[iveg] = deadstemn[band][iveg];
	  cn[band].deadstemn_storage[iveg] = deadstemn_storage[band][iveg];
	  cn[band].deadstemn_xfer[iveg] = deadstemn_xfer[band][iveg];
	  cn[band].livecrootn[iveg] = livecrootn[band][iveg];
	  cn[band].livecrootn_storage[iveg] = livecrootn_storage[band][iveg];
	  cn[band].livecrootn_xfer[iveg] = livecrootn_xfer[band][iveg];
	  cn[band].deadcrootn[iveg] = deadcrootn[band][iveg];
	  cn[band].deadcrootn_storage[iveg] = deadcrootn_storage[band][iveg];
	  cn[band].deadcrootn_xfer[iveg] = deadcrootn_xfer[band][iveg];
	  cn[band].retransn[iveg] = retransn[band][iveg];
	  cn[band].npool[iveg] = npool[band][iveg];
	  cn[band].pft_ntrunc[iveg] = pft_ntrunc[band][iveg];
	}
      cn[band].decl = decl[band];
      cn[band].fpi = fpi[band];
      cn[band].fpg = fpg[band];
      cn[band].annsum_counter = annsum_counter[band];
      cn[band].cannsum_npp = cannsum_npp[band];
      cn[band].cannavg_t2m = cannavg_t2m[band];
      for(nidx = 0; nidx < Nnodes; nidx++)
	{
	  cn[band].watfc[nidx] = watfc[band][nidx];
	}
      cn[band].me = me[band];
      cn[band].fire_prob = fire_prob[band];
      cn[band].mean_fire_prob = mean_fire_prob[band];
      cn[band].fireseasonl = fireseasonl[band];
      cn[band].ann_farea_burned = ann_farea_burned[band];
      cn[band].cwdc = cwdc[band];
      cn[band].litr1c = litrlabc[band];
      cn[band].litr2c = litrcellc[band];
      cn[band].litr3c = litrligc[band];
      cn[band].soil1c = soilcfast[band];
      cn[band].soil2c = soilcmid[band];
      cn[band].soil3c = soilcslo1[band];
      cn[band].soil4c = soilcslo2[band];
      cn[band].seedc = seedc[band];
      cn[band].col_ctrunc = col_ctrunc[band];
      cn[band].totlitc = totlitc[band];
      cn[band].totsomc = totsomc[band];
      cn[band].totcolc = totcolc[band];
      cn[band].prod10c = prod10c[band];
      cn[band].prod100c = prod100c[band];
      cn[band].cwdn = cwdn[band];
      cn[band].litr1n = litrlabn[band];
      cn[band].litr2n = litrcelln[band];
      cn[band].litr3n = litrlign[band];
      cn[band].soil1n = soilnfast[band];
      cn[band].soil2n = soilnmid[band];
      cn[band].soil3n = soilnslo1[band];
      cn[band].soil4n = soilnslo2[band];
      cn[band].sminn = soilminn[band];
      cn[band].col_ntrunc = col_ntrunc[band];
      cn[band].seedn = seedn[band];
      cn[band].totcoln = totcoln[band];
      cn[band].prod10n = prod10n[band];
      cn[band].prod100n = prod100n[band];
      cn[band].hr = hr[band];
      cn[band].nee = nee[band];
      cn[band].nep = nep[band];

      for(iveg = 0; iveg < Nveg; iveg++)

	{

	  cell[iveg][band].CLitter = cn[band].totlitc;
	  cell[iveg][band].CInter = cn[band].soil2c;
	  cell[iveg][band].CSlow = cn[band].soil3c;
	  cell[iveg][band].RhLitter = lithr[band];
	  cell[iveg][band].RhTot = cn[band].hr;

	  veg[iveg][band].LAI = pft2cov(cn[band].LAI, veg_con[iveg].veg_class);
	  veg[iveg][band].aPAR = pft2cov(apar[band], \
					 veg_con[iveg].veg_class);
	  veg[iveg][band].Ci = pft2cov(ci[band], veg_con[iveg].veg_class);
	  veg[iveg][band].rc = pft2cov(rc[band], veg_con[iveg].veg_class);
	  veg[iveg][band].GPP = pft2cov(cn[band].gpp, veg_con[iveg].veg_class);
	  veg[iveg][band].Rphoto = pft2cov(fpsn[band], \
					   veg_con[iveg].veg_class);
	  veg[iveg][band].Rdark = pft2cov(darkr[band], \
					   veg_con[iveg].veg_class);
	  veg[iveg][band].Rmaint = pft2cov(mr[band], \
					   veg_con[iveg].veg_class);
	  veg[iveg][band].Rgrowth = pft2cov(gr[band], \
					    veg_con[iveg].veg_class);
	  veg[iveg][band].Raut = pft2cov(cn[band].ar, veg_con[iveg].veg_class);
	  veg[iveg][band].NPP = pft2cov(cn[band].npp, veg_con[iveg].veg_class);
	  veg[iveg][band].Litterfall = pft2cov(litfall[band], \
					       veg_con[iveg].veg_class);
	  veg[iveg][band].AnnualNPP = pft2cov(cn[band].annsum_npp, \
					      veg_con[iveg].veg_class);

	}

    }

  return(err);

}
