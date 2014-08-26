/*
 * SUMMARY:      SnowInterception.c - simulates snow interception and release
 * USAGE:        
 *
 * AUTHOR:       Brian Connelly and Pascal Storck
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       pstorck@u.washington.edu
 * ORIG-DATE:    29-Aug-1996 at 13:42:17
 * LAST-MOD: Fri Apr 25 14:58:10 2003 by Keith Cherkauer <cherkaue@u.washington.edu>
 * DESCRIPTION:  Calculates the interception and subsequent release of
 *               by the forest canopy using an energy balance approach
 * DESCRIP-END.
 * FUNCTIONS:    SnowInterception()
 * COMMENTS:     Modified for use with VIC-NL code by Keith Cherkauer
 *               on 4-9-98   
 */

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

/*****************************************************************************
  Function name: SnowInterception()

  Purpose      : Calculate snow interception and release by the canopy

  Comments     : Only the top canopy layer is taken into account for snow
                 interception.  Snow interception by lower canopy is
                 disregarded.  Rain water CAN be intercepted by lower canopy
                 layers (similar to InterceptionStorage()).
                 Of course:  NO vegetation -> NO interception

  Modifications:
  06-98 included maximum structural loading to prevent the model
        from loading the canopy with more snow than it can handle
        structurally.                                             PXS
  09-98 aerodynamic resistances in the canopy when snow has been  
        intercepted is increased by a factor of 10:  include REF
        Journal of Hydrology, 1998                            KAC, GO'D
  11-00 energy balance components are now returned to the VIC model
        so that they can be reported as energy balance components
        in model output.                                           KAC
  xx-xx-02 Modified to handle new variables required to close the
        canopy energy balance.                                     KAC
  2-20-03 Added check of intercepted before mass balance error 
        calculation.  Since interception quantity can be greater
        than Wd_max when snow is present, the routine could return
        dew values higher than maximum in a time step when intercepted
        snow melted.  If there is no other snow, and distributed 
        precipitation is active this could cause the model to crash
        in redistribute_during_storm as it will be unable to conserve
        water when too much water is in the canopy.                 KAC
  04-Jun-04 Fixed typo in line 730.  Changed SPATIAL_FRoST to
	    SPATIAL_FROST.					TJB
  04-Jun-04 Added descriptive error message to beginning of screen dump
	    in error_print_canopy_energy_bal.			TJB
  21-Sep-04 Added ErrorString to store error messages from
	    root_brent.						TJB
  28-Sep-04 Added Ra_used to store the aerodynamic resistance used in
	    flux calculations.					TJB
  2007-Apr-11 Modified to handle grid cell errors by returning to the
	      main subroutine, rather than ending the simulation.	KAC
  2007-Aug-27 Modified to drop canopy snow if it is especially thin, which
	      should improve the numeric stability of the canopy energy
	      balance solution.						KAC via TJB
  2007-Aug-31 Checked root_brent return value against -998 rather than
	      -9998.							JCA
  2009-May-22 Added TFALLBACK value to options.CONTINUEONERROR.  This
	      allows simulation to continue when energy balance fails
	      to converge by using previous T value.			TJB
  2009-Jun-19 Added T flag to indicate whether TFALLBACK occurred.	TJB
  2009-Sep-19 Added T fbcount to count TFALLBACK occurrences.		TJB
  2009-Dec-11 Replaced "assert" statements with "if" statements.	TJB
  2010-Apr-26 Replaced individual forcing variables with atmos_data
	      in argument list.						TJB
  2013-Jul-25 Added photosynthesis terms.				TJB
  2013-Dec-27 Moved SPATIAL_FROST to options_struct.			TJB
  2014-Mar-28 Removed DIST_PRCP option.					TJB
  2014-May-05 Added logic to handle LAI = 0.				TJB
*****************************************************************************/
int snow_intercept(double  Dt,
		   double  F,  
		   double  LAI, 
		   double  Le, 
		   double  LongOverIn, // incominf LW from sky
		   double  LongUnderOut, // incoming LW from understroy
		   double  MaxInt, // maximum interception capacity
		   double  ShortOverIn,  // incoming SW to overstory
		   double  ShortUnderIn,  // incoming SW to understory
		   double  Tcanopy, // canopy air temperature
		   double  bare_albedo, // albedo of snow-free ground
		   double *AdvectedEnergy,
		   double *AlbedoOver, // overstory albedo
		   double *IntRain, // intercepted rain
		   double *IntSnow, // intercepted snow
		   double *LatentHeat, // latent heat from overstory
		   double *LatentHeatSub, // sublimation heat from overstory
		   double *LongOverOut, // longwave emitted by canopy
		   double *MeltEnergy, 
		   double *NetLongOver,
		   double *NetShortOver,
		   double *Ra,
		   double *Ra_used,
		   double *RainFall,
		   double *SensibleHeat,
		   double *SnowFall, 
		   double *Tfoliage, 
		   char   *Tfoliage_fbflag, 
		   int    *Tfoliage_fbcount, 
		   double *TempIntStorage, 
		   double *VaporMassFlux,
		   double *Wind,   
		   double *displacement,
		   double *ref_height,
		   double *roughness,
		   float  *root, 
		   int     UnderStory,
		   int     band, 
		   int     hour, 
		   int     iveg, 
		   int     month, 
		   int     rec,
		   int     hidx,
		   int     veg_class,
		   double *CanopLayerBnd,
		   double *dryFrac,
		   atmos_data_struct *atmos,
		   layer_data_struct *layer,
		   soil_con_struct   *soil_con,
		   veg_var_struct    *veg_var)
{

  extern option_struct options;

  /* double AdvectedEnergy; */         /* Energy advected by the rain (W/m2) */
  double BlownSnow;              /* Depth of snow blown of the canopy (m) */
  double DeltaSnowInt;           /* Change in the physical swe of snow
				    interceped on the branches. (m) */
  double Drip;                   /* Amount of drip from intercepted snow as a
				    result of snowmelt (m) */
  double ExcessSnowMelt;         /* Snowmelt in excess of the water holding
				    capacity of the tree (m) */
  double EsSnow;                 /* saturated vapor pressure in the snow pack
                                   (Pa)  */
  double InitialSnowInt;         /* Initial intercepted snow (m) */ 
  double InitialWaterInt;        /* Initial intercepted water (snow and rain)
                                    (m) */ 
  double IntRainOrg;
  double Ls;                     /* Latent heat of sublimation (J/(kg K) */
  double MassBalanceError;       /* Mass blalnce to make sure no water is
				    being destroyed/created (m) */
  double MaxWaterInt;            /* Water interception capacity (m) */  
  double MaxSnowInt;             /* Snow interception capacity (m) */
  double NetRadiation;
  double PotSnowMelt;            /* Potential snow melt (m) */
  double RainThroughFall;        /* Amount of rain reaching to the ground (m)
				  */ 
  double RefreezeEnergy;         /* Energy available for refreezing or melt */
  double ReleasedMass;           /* Amount of mass release of intercepted snow
				    (m) */ 
  /* double SensibleHeat; */           /* Sensible heat flux (W/m2) */
  double SnowThroughFall;        /* Amount of snow reaching to the ground (m)
				  */ 
  double Tmp;                    /* Temporary variable */

  double Imax1;                  /* maxium water intecept regardless of temp */
  double IntRainFract;           /* Fraction of intercpeted water which is 
				    liquid */
  double IntSnowFract;           /* Fraction of intercepted water which is 
				    solid */
  double Overload;               /* temp variable to calculated structural 
				    overloading */
  double Qnet;                   /* temporary storage of energy balance 
				    error */
  double Tupper;
  double Tlower;
  double Evap;
  double OldTfoliage;

  double  AirDens;
  double  EactAir;
  double  Press; // atmospheric pressure
  double  Vpd; // vapor pressure defficit
  double  shortwave; //
  double  Catm; //

  char ErrorString[MAXSTRING];

  AirDens   = atmos->density[hidx];
  EactAir   = atmos->vp[hidx];
  Press     = atmos->pressure[hidx];
  Vpd       = atmos->vpd[hidx];
  shortwave = atmos->shortwave[hidx];
  Catm      = atmos->Catm[hidx];

  /* Initialize Tfoliage_fbflag */
  *Tfoliage_fbflag = 0;

  /* Convert Units from VIC (mm -> m) */
  *RainFall /= 1000.;
  *SnowFall /= 1000.;
  *IntRain  /= 1000.;
  MaxInt    /= 1000.;
  IntRainOrg = *IntRain;

  /* Initialize Drip, H2O balance, and mass release variables. */
  
  InitialWaterInt = *IntSnow + *IntRain;
  
  *IntSnow /= F;
  *IntRain /= F;
  
  InitialSnowInt = *IntSnow;
  
  Drip = 0.0;
  ReleasedMass = 0.0;

  OldTfoliage = *Tfoliage;

  /* Determine the maximum snow interception water equivalent.           
     Kobayashi, D., 1986, Snow Accumulation on a Narrow Board,           
     Cold Regions Science and Technology, (13), pp. 239-245.           
     Figure 4. */  
  
  Imax1 = 4.0* LAI_SNOW_MULTIPLIER * LAI;

  if ((*Tfoliage) < -1.0 && (*Tfoliage) > -3.0)
    MaxSnowInt = ((*Tfoliage)*3.0/2.0) + (11.0/2.0);
  else if ((*Tfoliage) > -1.0) 
    MaxSnowInt = 4.0;  
  else
    MaxSnowInt = 1.0;
  
  /* therefore LAI_ratio decreases as temp decreases */
  
  MaxSnowInt *= LAI_SNOW_MULTIPLIER * LAI;
  
  /* Calculate snow interception. */  
  
  if (MaxSnowInt > 0) {
    DeltaSnowInt = (1-*IntSnow/MaxSnowInt) * *SnowFall; 
    if (DeltaSnowInt + *IntSnow > MaxSnowInt) 
      DeltaSnowInt = MaxSnowInt - *IntSnow;
    if (DeltaSnowInt < 0.0)  
      DeltaSnowInt = 0.0;
  }
  else {
    DeltaSnowInt = 0.0;
  }
  
  /* Reduce the amount of intercepted snow if windy and cold.         
     Ringyo Shikenjo Tokyo, #54, 1952.                                
     Bulletin of the Govt. Forest Exp. Station,                       
     Govt. Forest Exp. Station, Meguro, Tokyo, Japan.                 
     FORSTX 634.9072 R475r #54.                                       
     Page 146, Figure 10.                                               
     
     Reduce the amount of intercepted snow if snowing, windy, and     
     cold (< -3 to -5 C).                                             
     Schmidt and Troendle 1992 western snow conference paper. */  
  
  if ((*Tfoliage) < -3.0 && DeltaSnowInt > 0.0 && Wind[1] > 1.0) {
    BlownSnow = (0.2 * Wind[1] - 0.2) * DeltaSnowInt;
    if (BlownSnow >= DeltaSnowInt) 
      BlownSnow = DeltaSnowInt;
    DeltaSnowInt -= BlownSnow;
  }
  
  /* now update snowfall and total accumulated intercepted snow amounts */

  if (*IntSnow +  DeltaSnowInt > Imax1) DeltaSnowInt =0.0; 
  
  /* pixel depth    */ 
  SnowThroughFall = (*SnowFall - DeltaSnowInt) * F + (*SnowFall) * (1 - F);
  
  // Snow in canopy too thin for EB calculations; let it fall through
  if ( *SnowFall == 0 && *IntSnow < MIN_SWQ_EB_THRES ) {
    SnowThroughFall += *IntSnow;
    DeltaSnowInt -= *IntSnow;
  }

  /* physical depth */
  *IntSnow += DeltaSnowInt;
  if (*IntSnow < SMALL)
    *IntSnow = 0.0;
  
  /* Calculate amount of rain intercepted on branches and stored in
     intercepted snow. */  
  
  /* physical depth */
  MaxWaterInt = LIQUID_WATER_CAPACITY * (*IntSnow) + MaxInt;
  
  if ((*IntRain + *RainFall) <= MaxWaterInt) {
    /* physical depth */
    *IntRain += *RainFall;
    /* pixel depth */
    RainThroughFall = *RainFall * (1 - F);      
  }
  else {
    /* pixel depth */
    RainThroughFall = (*IntRain + *RainFall - MaxWaterInt) * F + 
      (*RainFall * (1 - F));
    /* physical depth */
    *IntRain = MaxWaterInt;
  }

  // Liquid water in canopy too thin for EB calculations; let it fall through
  if ( *RainFall == 0 && *IntRain < MIN_SWQ_EB_THRES ) {
    RainThroughFall += *IntRain;
    *IntRain = 0.0;
  }

  /* at this point we have calculated the amount of snowfall intercepted and
     the amount of rainfall intercepted.  These values have been 
     appropriately subtracted from SnowFall and RainFall to determine 
     SnowThroughfall and RainThroughfall.  However, we can end up with the 
     condition that the total intercepted rain plus intercepted snow is 
     greater than the maximum bearing capacity of the tree regardless of air 
     temp (Imax1).  The following routine will adjust *IntRain and *IntSnow 
     by triggering mass release due to overloading.  Of course since *IntRain
     and *IntSnow are mixed, we need to slough them of as fixed fractions  */

  if (*IntRain + *IntSnow > Imax1) { /*then trigger structural unloading*/
    Overload = (*IntSnow + *IntRain) - Imax1;
    IntRainFract= *IntRain/(*IntRain + *IntSnow);
    IntSnowFract = *IntSnow/(*IntRain + *IntSnow);
    *IntRain = *IntRain - Overload*IntRainFract;
    *IntSnow = *IntSnow - Overload*IntSnowFract;
    RainThroughFall = RainThroughFall + (Overload * IntRainFract) * F;
    SnowThroughFall = SnowThroughFall + (Overload * IntSnowFract) * F;
  }

  // If we've lost all intercepted moisture, we've essentially lost the thermal
  // mass of the canopy and Tfoliage should be equal to Tcanopy
  if (*IntRain + *IntSnow < SMALL) {
    *Tfoliage = Tcanopy;
  }

  /* Calculate the net radiation at the canopy surface, using the canopy 
     temperature.  The outgoing longwave is subtracted twice, because the 
     canopy radiates in two directions */

  Tupper = Tlower = MISSING;

  if ( *IntSnow > 0 || *SnowFall > 0 ) {
    /* Snow present or accumulating in the canopy */

    *AlbedoOver = NEW_SNOW_ALB; // albedo of intercepted snow in canopy
    *NetShortOver = (1. - *AlbedoOver) * ShortOverIn; // net SW in canopy

    Qnet = solve_canopy_energy_bal(0., band, month, rec, Dt, 
				   soil_con->elevation, soil_con->max_moist,
				   soil_con->Wcr, soil_con->Wpwp, 
				   soil_con->depth, 
				   soil_con->frost_fract, 
				   AirDens, EactAir, Press, Le, 
				   Tcanopy, Vpd,
				   shortwave, Catm, dryFrac,
				   &Evap, Ra, Ra_used,
				   *RainFall, Wind, UnderStory, iveg, 
				   veg_class, displacement, ref_height, 
				   roughness, root, CanopLayerBnd, 
				   IntRainOrg, *IntSnow, 
				   IntRain, layer, veg_var, 
				   LongOverIn, LongUnderOut, 
				   *NetShortOver, 
				   AdvectedEnergy, 
				   LatentHeat, LatentHeatSub, 
				   LongOverOut, NetLongOver, &NetRadiation, 
				   &RefreezeEnergy, SensibleHeat, 
				   VaporMassFlux);

    if ( Qnet != 0 ) {
      /* Intercepted snow not melting - need to find temperature */
      Tupper = 0;
      if ( (*Tfoliage) <= 0.) Tlower = (*Tfoliage) - SNOW_DT;
      else Tlower = -SNOW_DT;
    }
    else *Tfoliage = 0.;  // intercepted snow is melting


  }
  else {
    /* No snow in canopy */
    *AlbedoOver = bare_albedo;
    *NetShortOver = (1. - *AlbedoOver) * ShortOverIn; // net SW in canopy
    Qnet = -9999;
    Tupper = (*Tfoliage) + SNOW_DT;
    Tlower = (*Tfoliage) - SNOW_DT;

  }

  if ( Tupper != MISSING && Tlower != MISSING ) {

    *Tfoliage = root_brent(Tlower, Tupper, ErrorString, func_canopy_energy_bal,  band, 
			   month, rec, Dt, soil_con->elevation, soil_con->max_moist, 
			   soil_con->Wcr, soil_con->Wpwp, soil_con->depth, 
			   soil_con->frost_fract, 
			   AirDens, EactAir, Press, Le, 
			   Tcanopy, Vpd, shortwave, Catm, dryFrac,
			   &Evap, Ra, Ra_used,
			   *RainFall, Wind, UnderStory, iveg, 
			   veg_class, displacement, ref_height, 
			   roughness, root, CanopLayerBnd, IntRainOrg, *IntSnow, 
			   IntRain, layer, veg_var, 
			   LongOverIn, LongUnderOut, 
			   *NetShortOver, 
			   AdvectedEnergy, 
			   LatentHeat, LatentHeatSub, 
			   LongOverOut, NetLongOver, &NetRadiation, 
			   &RefreezeEnergy, SensibleHeat, 
			   VaporMassFlux);
    
    if ( *Tfoliage <= -998 ) {
      if (options.TFALLBACK) {
        *Tfoliage = OldTfoliage;
        *Tfoliage_fbflag = 1;
        (*Tfoliage_fbcount)++;
      }
      else { 
        Qnet = error_calc_canopy_energy_bal(*Tfoliage, band, month, rec, Dt, 
					    soil_con->elevation, soil_con->max_moist, 
					    soil_con->Wcr, soil_con->Wpwp, 
					    soil_con->depth, 
					    soil_con->frost_fract, 
					    AirDens, EactAir, Press, Le, 
					    Tcanopy, Vpd, shortwave, Catm, dryFrac,
					    &Evap, Ra, Ra_used,
					    *RainFall, Wind, UnderStory, iveg, 
					    veg_class, displacement, ref_height, 
					    roughness, root, CanopLayerBnd, IntRainOrg, *IntSnow, 
					    IntRain, layer, veg_var, 
					    LongOverIn, LongUnderOut, *NetShortOver, 
					    AdvectedEnergy, 
					    LatentHeat, LatentHeatSub, 
					    LongOverOut, NetLongOver, &NetRadiation, 
					    &RefreezeEnergy, SensibleHeat, 
					    VaporMassFlux, ErrorString);
        return( ERROR );
      }
    }
    
    Qnet = solve_canopy_energy_bal(*Tfoliage, band, month, rec, Dt, 
				   soil_con->elevation, soil_con->max_moist, soil_con->Wcr, 
				   soil_con->Wpwp, soil_con->depth, 
				   soil_con->frost_fract, 
				   AirDens, EactAir, Press, Le, 
				   Tcanopy, Vpd, shortwave, Catm, dryFrac,
				   &Evap, Ra, Ra_used,
				   *RainFall, Wind, UnderStory, iveg, 
				   veg_class, displacement, ref_height, 
				   roughness, root, CanopLayerBnd,
				   IntRainOrg, *IntSnow, 
				   IntRain, layer, veg_var, 
				   LongOverIn, LongUnderOut, 
				   *NetShortOver, 
				   AdvectedEnergy, 
				   LatentHeat, LatentHeatSub, 
				   LongOverOut, NetLongOver, &NetRadiation, 
				   &RefreezeEnergy, SensibleHeat, 
				   VaporMassFlux);

  }

  if ( *IntSnow <= 0 )
    RainThroughFall = veg_var->throughfall / 1000.;

  RefreezeEnergy *= Dt;

  /* if RefreezeEnergy is positive it means energy is available to melt the
     intercepted snow in the canopy.  If it is negative, it means that 
     intercepted water will be refrozen */
  
  /* Update maximum water interception storage */
  
  MaxWaterInt = LIQUID_WATER_CAPACITY * (*IntSnow) + MaxInt;

  /* Convert the vapor mass flux from a flux to a depth per interval */
  *VaporMassFlux *= Dt;
  
  if ( *Tfoliage == 0 ) {

    if (-(*VaporMassFlux) > *IntRain) {
      *VaporMassFlux = -(*IntRain);
      *IntRain = 0.;
    }
    else
      *IntRain += *VaporMassFlux;

    if ( RefreezeEnergy < 0 ) {
      /* intercepted snow is ripe, melt can occur */
      PotSnowMelt = min((-RefreezeEnergy/Lf/RHO_W), *IntSnow);
      *MeltEnergy -= (Lf * PotSnowMelt * RHO_W) / (Dt);
    }
    else {
      /* snow temperature is below freezing, no melt occurs */
      PotSnowMelt = 0;
      *MeltEnergy -= (Lf * PotSnowMelt * RHO_W) / (Dt);
    }
    
    if ((*IntRain + PotSnowMelt) <= MaxWaterInt) {

      *IntSnow -= PotSnowMelt;
      *IntRain += PotSnowMelt;
      PotSnowMelt = 0.0;

    }
    
    else {

      ExcessSnowMelt = PotSnowMelt + *IntRain - MaxWaterInt;
      
      *IntSnow -= MaxWaterInt - (*IntRain);
      *IntRain = MaxWaterInt;
      if (*IntSnow < 0.0) 
        *IntSnow = 0.0;
      
      if (SnowThroughFall > 0.0 && 
	  InitialSnowInt <= MIN_INTERCEPTION_STORAGE) {
        /* Water in excess of MaxWaterInt has been generated.  If it is 
           snowing and there was little intercepted snow at the beginning 
	   of the time step ( <= MIN_INTERCEPTION_STORAGE), then allow the 
	   snow to melt as it is intercepted */
        Drip += ExcessSnowMelt; 
        *IntSnow -= ExcessSnowMelt;
        if (*IntSnow < 0.0) 
          *IntSnow = 0.0;
      }
      else 
      /* Else, SnowThroughFall = 0.0 or SnowThroughFall > 0.0 and there is a 
         substantial amount of intercepted snow at the beginning of the time 
         step ( > MIN_INTERCEPTION_STORAGE).  Snow melt may generate mass 
         release. */
        *TempIntStorage += ExcessSnowMelt;
      
      MassRelease(IntSnow, TempIntStorage, &ReleasedMass, &Drip);
    }
    
    /* If intercepted snow has melted, add the water it held to drip */
    
    MaxWaterInt = LIQUID_WATER_CAPACITY * (*IntSnow) + MaxInt;
    if (*IntRain > MaxWaterInt) {
      Drip += *IntRain - MaxWaterInt;
      *IntRain = MaxWaterInt;
    }
  }  
  
  else /* else (RefreezeEnergy <= 0.0) */ {
    
    /* Reset *TempIntStorage to 0.0 when energy balance is negative */
    
    *TempIntStorage = 0.0;
    
    /* Refreeze as much surface water as you can */
    
    if (-RefreezeEnergy > - (*IntRain) * Lf) {
      *IntSnow += fabs(RefreezeEnergy) / Lf;
      *IntRain -= fabs(RefreezeEnergy) / Lf;

      *MeltEnergy += (fabs(RefreezeEnergy) * RHO_W) / (Dt);

      RefreezeEnergy = 0.0;
    }

    else {
      
      /* All of the water in the surface layer has been frozen. */
      
      *IntSnow += *IntRain;
     
      /* Added on April 8 as a test */
      /*       RefreezeEnergy += *IntRain*Lf; */
      /*       *VaporMassFlux = MAX(*VaporMassFlux,  */
      /*                            RefreezeEnergy/(Ls * RHO_W)); */
      
      /* Energy released by freezing of intercepted water is added to the 
         MeltEnergy */

      *MeltEnergy += (Lf * *IntRain * RHO_W) / (Dt);
      *IntRain = 0.0;
      
    } 
    
    if (-(*VaporMassFlux) > *IntSnow) {
      *VaporMassFlux = -(*IntSnow);
      *IntSnow = 0.0;
    }
    else
      *IntSnow += *VaporMassFlux;
  } 
  
  *IntSnow *= F;
  *IntRain *= F;
  *MeltEnergy *= F;
  *VaporMassFlux *= F;
  Drip           *= F;
  ReleasedMass   *= F;
  
  if ( *IntSnow == 0 && *IntRain > MaxInt ) {
    // if snow has melted, make sure canopy is not retaining extra water
    RainThroughFall += *IntRain - MaxInt;
    *IntRain = MaxInt;
  }

  /* Calculate intercepted H2O balance. */  
  
  MassBalanceError = (InitialWaterInt - (*IntSnow + *IntRain))
                   + (*SnowFall + *RainFall) - (SnowThroughFall
                   + RainThroughFall + Drip + ReleasedMass)
                   + *VaporMassFlux;

  *RainFall = RainThroughFall + Drip;
  *SnowFall = SnowThroughFall + ReleasedMass;

  /* Convert Units to VIC (m -> mm) */
  *VaporMassFlux *= -1.;
  *RainFall *= 1000.;
  *SnowFall *= 1000.;
  *IntRain  *= 1000.;

  /*** FIX THIS ***/
  *MeltEnergy = RefreezeEnergy / Dt;

  return( 0 );

}

double solve_canopy_energy_bal(double Tfoliage, ...)
{
  va_list  ap;
  double   Qnet;

  va_start(ap, Tfoliage);
  Qnet = func_canopy_energy_bal(Tfoliage, ap);
  va_end(ap);

  return Qnet;
}

double error_calc_canopy_energy_bal(double Tfoliage, ...)
{
  va_list  ap;
  double   Qnet;

  va_start(ap, Tfoliage);
  Qnet = error_print_canopy_energy_bal(Tfoliage, ap);
  va_end(ap);

  return Qnet;
}

double error_print_canopy_energy_bal(double Tfoliage, va_list ap)
{  

  extern option_struct options;

  /* General Model Parameters */
  int     band;
  int     month;
  int     rec;

  double  delta_t;
  double  elevation;

  double *Wmax;
  double *Wcr;
  double *Wpwp;
  double *depth;
  double *frost_fract;

  /* Atmopheric Condition and Forcings */
  double  AirDens;
  double  EactAir;
  double  Press;
  double  Le;
  double  Tcanopy;
  double  Vpd;
  double  shortwave;
  double  Catm;
  double  *dryFrac;

  double *Evap;
  double *Ra;
  double *Ra_used;
  double Rainfall;
  double *Wind;

  /* Vegetation Terms */
  int     UnderStory;
  int     iveg;
  int     veg_class;

  double *displacement;
  double *ref_height;
  double *roughness;

  float  *root;
  double *CanopLayerBnd;

  /* Water Flux Terms */
  double  IntRain;
  double  IntSnow;

  double *Wdew;

  layer_data_struct *layer;
  veg_var_struct    *veg_var;

  /* Energy Flux Terms */
  double  LongOverIn;
  double  LongUnderOut;
  double  NetShortOver;

  double *AdvectedEnergy;
  double *LatentHeat;
  double *LatentHeatSub;
  double *LongOverOut;
  double *NetLongOver;
  double *NetRadiation;
  double *RefreezeEnergy;
  double *SensibleHeat;
  double *VaporMassFlux;

  char *ErrorString;
  int   cidx;

  /** Read variables from variable length argument list **/

  /* General Model Parameters */
  band    = (int) va_arg(ap, int);
  month   = (int) va_arg(ap, int);
  rec     = (int) va_arg(ap, int);

  delta_t   = (double) va_arg(ap, double);
  elevation = (double) va_arg(ap, double);

  Wmax  = (double *) va_arg(ap, double *);
  Wcr   = (double *) va_arg(ap, double *);
  Wpwp  = (double *) va_arg(ap, double *);
  depth = (double *) va_arg(ap, double *);
  frost_fract = (double *) va_arg(ap, double *);

  /* Atmopheric Condition and Forcings */
  AirDens = (double) va_arg(ap, double);
  EactAir = (double) va_arg(ap, double);
  Press   = (double) va_arg(ap, double);
  Le      = (double) va_arg(ap, double);
  Tcanopy = (double) va_arg(ap, double);
  Vpd     = (double) va_arg(ap, double);
  shortwave=(double) va_arg(ap, double);
  Catm    = (double) va_arg(ap, double);
  dryFrac = (double *) va_arg(ap, double *);

  Evap     = (double *) va_arg(ap, double *);
  Ra       = (double *) va_arg(ap, double *);
  Ra_used  = (double *) va_arg(ap, double *);
  Rainfall = (double) va_arg(ap, double);
  Wind     = (double *) va_arg(ap, double *);

  /* Vegetation Terms */
  UnderStory = (int) va_arg(ap, int);
  iveg       = (int) va_arg(ap, int);
  veg_class  = (int) va_arg(ap, int);

  displacement = (double *) va_arg(ap, double *);
  ref_height   = (double *) va_arg(ap, double *);
  roughness    = (double *) va_arg(ap, double *);

  root = (float *) va_arg(ap, float *);
  CanopLayerBnd = (double *) va_arg(ap, double *);

  /* Water Flux Terms */
  IntRain = (double) va_arg(ap, double);
  IntSnow = (double) va_arg(ap, double);

  Wdew    = (double *) va_arg(ap, double *);

  layer   = (layer_data_struct *) va_arg(ap, layer_data_struct *);
  veg_var = (veg_var_struct *) va_arg(ap, veg_var_struct *);

  /* Energy Flux Terms */
  LongOverIn       = (double) va_arg(ap, double);
  LongUnderOut     = (double) va_arg(ap, double);
  NetShortOver     = (double) va_arg(ap, double);

  AdvectedEnergy     = (double *) va_arg(ap, double *);
  LatentHeat         = (double *) va_arg(ap, double *);
  LatentHeatSub      = (double *) va_arg(ap, double *);
  LongOverOut        = (double *) va_arg(ap, double *);
  NetLongOver        = (double *) va_arg(ap, double *);
  NetRadiation       = (double *) va_arg(ap, double *);
  RefreezeEnergy     = (double *) va_arg(ap, double *);
  SensibleHeat       = (double *) va_arg(ap, double *);
  VaporMassFlux      = (double *) va_arg(ap, double *);
  ErrorString        = (char *) va_arg(ap, char *);

  /** Print variable info */
  fprintf(stderr, "%s", ErrorString);
  fprintf(stderr, "ERROR: snow_intercept failed to converge to a solution in root_brent.  Variable values will be dumped to the screen, check for invalid values.\n");

  /* General Model Parameters */
  printf("band = %i\n", band);
  printf("month = %i\n",     month);
  printf("rec = %i\n", rec);
  
  printf("delta_t = %f\n", delta_t);
  printf("elevation = %f\n",  elevation);
  
  printf("*Wmax = %f\n", *Wmax);
  printf("*Wcr = %f\n", *Wcr);
  printf("*Wpwp = %f\n", *Wpwp);
  printf("*depth = %f\n", *depth);
#if SPATIAL_FRoST
  printf(" = %f\n", *frost_fract);
#endif
  
  /* Atmopheric Condition and Forcings */
  printf("AirDens = %f\n",  AirDens);
  printf("EactAir = %f\n",  EactAir);
  printf("Press = %f\n",  Press);
  printf("Le = %f\n",  Le);
  printf("Ra = [%f, %f]\n",  Ra[1], Ra[UnderStory]);
  printf("Ra_used = %f\n",  *Ra_used);
  printf("Tcanopy = %f\n",  Tcanopy);
  printf("Vpd = %f\n",  Vpd);
  printf("shortwave = %f\n",  shortwave);
  printf("Catm = %f\n",  Catm);
  printf("dryFrac = %f\n",  *dryFrac);

  printf("Evap = %f\n", *Evap);
  printf("Rainfall = %f\n", Rainfall);
  printf("Wind = [%f, %f]\n",  Wind[1], Wind[UnderStory]);

  /* Vegetation Terms */
  printf("UnderStory = %i\n",     UnderStory);
  printf("iveg = %i\n",     iveg);
  printf("veg_class = %i\n",     veg_class);

  printf("displacement = [%f, %f]\n",  displacement[1], displacement[UnderStory]);
  printf("ref_height = [%f, %f]\n",  ref_height[1], ref_height[UnderStory]);
  printf("roughness = [%f, %f]\n",  roughness[1], roughness[UnderStory]);

  printf("root = %f\n", *root);

  if (options.CARBON) {
    printf("CanopLayerBnd =");
    for (cidx=0; cidx<options.Ncanopy; cidx++) {
      printf(" %f", CanopLayerBnd[cidx]);
    }
    printf("\n");
  }

  /* Water Flux Terms */
  printf("IntRain = %f\n",  IntRain);
  printf("IntSnow = %f\n",  IntSnow);

  printf("Wdew = %f\n", *Wdew);

  write_layer(layer, iveg, options.Nlayer, frost_fract, depth);
  write_vegvar(&(veg_var[0]),iveg);

  /* Energy Flux Terms */
  fprintf(stderr, "LongOverIn = %f\n",  LongOverIn);
  fprintf(stderr, "LongUnderOut = %f\n",  LongUnderOut);
  fprintf(stderr, "NetShortOver = %f\n",  NetShortOver);

  fprintf(stderr, "*AdvectedEnergy = %f\n", *AdvectedEnergy);
  fprintf(stderr, "*LatentHeat = %f\n", *LatentHeat);
  fprintf(stderr, "*LatentHeatSub = %f\n", *LatentHeatSub);
  fprintf(stderr, "*LongOverOut = %f\n", *LongOverOut);
  fprintf(stderr, "*NetLongOver = %f\n", *NetLongOver);
  fprintf(stderr, "*NetRadiation = %f\n", *NetRadiation);
  fprintf(stderr, "*RefreezeEnergy = %f\n", *RefreezeEnergy);
  fprintf(stderr, "*SensibleHeat = %f\n", *SensibleHeat);
  fprintf(stderr, "*VaporMassFlux = %f\n", *VaporMassFlux);

  /* call error handling routine */
  fprintf(stderr,"**********\n**********\nFinished dumping snow_intercept variables.\nTry increasing SNOW_DT to get model to complete cell.\nThen check output for instabilities.\n**********\n**********\n");

  return( ERROR );

}
