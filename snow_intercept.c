/*
 * SUMMARY:      SnowInterception.c - simulates snow interception and release
 * USAGE:        
 *
 * AUTHOR:       Brian Connelly and Pascal Storck
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       pstorck@u.washington.edu
 * ORIG-DATE:    29-Aug-1996 at 13:42:17
 * LAST-MOD: Thu Sep  3 15:59:59 1998 by VIC Administrator <vicadmin@u.washington.edu>
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

/*****************************************************************************
  Function name: SnowInterception()

  Purpose      : Calculate snow interception and release by the canopy

  Required     :
    int Dt                 - Model timestep (hours)
    double F                - Fractional coverage
    double LAI              - Leaf Area Index
    double MaxInt           - Maximum rainfall interception storage (m)
    double BaseRa           - Aerodynamic resistance (uncorrected for
                             stability) (s/m)
    double AirDens          - Density of air (kg/m3)
    double EactAir          - Actual vapor pressure of air (Pa) 
    double Lv               - Latent heat of vaporization (J/kg3)
    PIXRAD *LocalRad       - Components of radiation balance for current pixel
                             (W/m2) 
    double Press            - Air pressure (Pa)
    double Tair             - Air temperature (C) 
    double Vpd	           - Vapor pressure deficit (Pa) 
    double Wind             - Wind speed (m/s)
    double *RainFall        - Amount of rain (m)
    double *Snowfall        - Amount of snow (m)
    double *IntRain         - Intercepted rain (m) 
    double *IntSnow         - Snow water equivalent of intercepted snow (m)
    double *TempIntStorage  - Temporary storage for snowmelt and rainfall
                             involved in mass release calculations (m)
    double *VaporMassFlux   - Vapor mass flux to/from intercepted snow
                             (m/timestep)
    double *Tcanopy         - Canopy temperature (C)
    double *MeltEnergy      - Energy used in heating and melting of the snow 
                             (W/m2)

  Returns      : none

  Modifies     :
    double *RainFall        - Amount of rain (m)
    double *Snowfall        - Amount of snow (m)
    double *IntRain         - Intercepted rain (m) 
    double *IntSnow         - Snow water equivalent of intercepted snow (m)
    double *TempIntStorage  - Temporary storage for snowmelt and rainfall
                             involved in mass release calculations (m)
    double *VaporMassFlux   - Vapor mass flux to/from intercepted snow
                             (m/timestep)  
    double *Tcanopy         - Canopy temperature (C)

  Comments     : Only the top canopy layer is taken into account for snow
                 interception.  Snow interception by lower canopy is
                 disregarded.  Rain water CAN be intercepted by lower canopy
                 layers (similar to InterceptionStorage()).
                 Of course:  NO vegetation -> NO interception

  Modifications:
  06-98 included maximum structural loading to prevent the model
        from loading the canopy with more snow than it can handle
        structurally.                                             PXS

*****************************************************************************/
void snow_intercept(double Dt, double F,  double LAI, 
		    double MaxInt, double Ra, double AirDens,
		    double EactAir, double Lv, double Shortwave,
		    double Longwave, double Press, double Tair, 
		    double Vpd, double Wind,  double *RainFall,
		    double *SnowFall, double *IntRain, double *IntSnow,
		    double *TempIntStorage, double *VaporMassFlux,
		    double *Tcanopy, double *MeltEnergy, int month, int rec,
		    int hour)
{
  extern option_struct options;
  FILE *ftmp;

  const char *Routine = "SnowInterception";
  double AdvectedEnergy;         /* Energy advected by the rain (W/m2) */
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
  double LatentHeat;             /* Latent heat flux (W/m2) */
  double LongOut;                /* Longwave radiation emitted by canopy 
                                   (W/m2) */
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
  double SensibleHeat;           /* Sensible heat flux (W/m2) */
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

  /* Convert Units from VIC (mm -> m) */
  *RainFall /= 1000.;
  *SnowFall /= 1000.;
  *IntRain  /= 1000.;
  MaxInt    /= 1000.;

  /* Initialize Drip, H2O balance, and mass release variables. */
  
  InitialWaterInt = *IntSnow + *IntRain;
  
  *IntSnow /= F;
  *IntRain /= F;
  
  InitialSnowInt = *IntSnow;
  
  Drip = 0.0;
  ReleasedMass = 0.0;
  
  /* Determine the maximum snow interception water equivalent.           
     Kobayashi, D., 1986, Snow Accumulation on a Narrow Board,           
     Cold Regions Science and Technology, (13), pp. 239-245.           
     Figure 4. */  
  
  Imax1 = 4.0* LAI_SNOW_MULTIPLIER * LAI;

  if (Tair < -1.0 && Tair > -3.0)
    MaxSnowInt = (Tair*3.0/2.0) + (11.0/2.0);
  else if (Tair > -1.0) 
    MaxSnowInt = 4.0;  
  else
    MaxSnowInt = 1.0;
  
  /* therefore LAI_ratio decreases as temp decreases */
  
  MaxSnowInt *= LAI_SNOW_MULTIPLIER * LAI;
  
  /* Calculate snow interception. */  
  
  DeltaSnowInt = (1-*IntSnow/MaxSnowInt) * *SnowFall; 
  if (DeltaSnowInt + *IntSnow > MaxSnowInt) 
    DeltaSnowInt = MaxSnowInt - *IntSnow;
  if (DeltaSnowInt < 0.0)  
    DeltaSnowInt = 0.0;
  
  /* Reduce the amount of intercepted snow if windy and cold.         
     Ringyo Shikenjo Tokyo, #54, 1952.                                
     Bulletin of the Govt. Forest Exp. Station,                       
     Govt. Forest Exp. Station, Meguro, Tokyo, Japan.                 
     FORSTX 634.9072 R475r #54.                                       
     Page 146, Figure 10.                                               
     
     Reduce the amount of intercepted snow if snowing, windy, and     
     cold (< -3 to -5 C).                                             
     Schmidt and Troendle 1992 western snow conference paper. */  
  
  if (Tair < -3.0 && DeltaSnowInt > 0.0 && Wind > 1.0) {
    BlownSnow = (0.2 * Wind - 0.2) * DeltaSnowInt;
    if (BlownSnow >= DeltaSnowInt) 
      BlownSnow = DeltaSnowInt;
    DeltaSnowInt -= BlownSnow;
  }
  
  /* now update snowfall and total accumulated intercepted snow amounts */

  if (*IntSnow +  DeltaSnowInt > Imax1) DeltaSnowInt =0.0; 
  
  /* pixel depth    */ 
  SnowThroughFall = (*SnowFall - DeltaSnowInt) * F + (*SnowFall) * (1 - F);
  
  /* physical depth */
  *IntSnow += DeltaSnowInt;
  
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
    RainThroughFall = RainThroughFall + (Overload*IntRainFract)*F;
    SnowThroughFall = SnowThroughFall + (Overload*IntSnowFract)*F;
  }
  
  /* The canopy temperature is assumed to be equal to the air temperature if 
     the air temperature is below 0C, otherwise the canopy temperature is 
     equal to 0C */
  
  if (Tair > 0.)
    *Tcanopy = 0.;
  else
    *Tcanopy = Tair;

  /* Calculate the net radiation at the canopy surface, using the canopy 
     temperature.  The outgoing longwave is subtracted twice, because the 
     canopy radiates in two directions */

  Tmp = *Tcanopy + 273.15;
  LongOut = STEFAN_B * (Tmp * Tmp * Tmp * Tmp);
  NetRadiation = (1.-NEW_SNOW_ALB)*Shortwave + Longwave - 2 * F * LongOut;
  NetRadiation /= F;

  /* Calculate the vapor mass flux between the canopy and the surrounding 
     air mass */
  
  EsSnow = svp(*Tcanopy); 
  if (*Tcanopy < 0.0)
    EsSnow *= 1.0 + .00972 * *Tcanopy + .000042 
      * pow((double)*Tcanopy,(double)2.0);
  *VaporMassFlux = AirDens * (0.622/Press) * (EactAir - EsSnow) / Ra; 
  *VaporMassFlux /= RHO_W; 

/*****
printf("
  AirTemp = %f    EactAir %f   EsSnow %f   Ra %f   AirDens %f   Press %f VaporMassFlux %f \n",
  Tair,           EactAir,     EsSnow,     Ra,      AirDens,     Press); 
*****/

  if (Vpd == 0.0 && *VaporMassFlux < 0.0)
    *VaporMassFlux = 0.0;
  
  ftmp = fopen("canopy_intercept.out","a");
  if(!options.FULL_ENERGY)
    fprintf(ftmp,"%f %f %f %f %f %f %f %f %f %f\n",
	    (float)rec+(float)hour/24., Tair, EactAir, EsSnow, Ra, 
	    AirDens, Press, Shortwave, Longwave, 
	    *VaporMassFlux * Dt * SECPHOUR); 
  else
    fprintf(ftmp,"%f %f %f %f %f %f %f %f %f %f\n",
	    (float)rec/24., Tair, EactAir, EsSnow, Ra, AirDens, 
	    Press, Shortwave, Longwave, *VaporMassFlux * Dt * SECPHOUR); 
  fclose(ftmp);

  /* Calculate the latent heat flux */

  Ls = (677. - 0.07 * *Tcanopy) * 4.1868 * 1000.;
  LatentHeat = Ls * *VaporMassFlux * RHO_W;

  /* Calculate the sensible heat flux */

  SensibleHeat = AirDens * Cp * (Tair - *Tcanopy)/Ra;
  
  /* Calculate the advected energy */

  AdvectedEnergy = (4186.8 * Tair * *RainFall)/(Dt * SECPHOUR);

  /* Calculate the amount of energy available for refreezing */

  RefreezeEnergy = SensibleHeat + LatentHeat + NetRadiation + AdvectedEnergy;




/*****
fprintf(stderr,"%i\t%i\t%.4lg\t%.4lg\t%.4lg\t%.4lg\t%.4lg\t%.4lg\t%.4lg\t%.4lg\t",
	rec,month,LAI,MaxInt,Ra,AirDens,EactAir,Lv,Shortwave,Longwave);
fprintf(stderr,"%.4lg\t%.4lg\t%.4lg\t%.4lg\t%.4lg\t%.4lg\t%.4lg\t%.4lg\t",
	Press,Tair,Vpd,Wind,*RainFall,*SnowFall,*IntRain,*IntSnow);
fprintf(stderr,"%.4lg\t%.4lg\t%.4lg\t%.4lg\t%.4lg\t",
	*TempIntStorage,*VaporMassFlux,*Tcanopy,RefreezeEnergy,NetRadiation);
fprintf(stderr,"%.4lg\t%.4lg\t%.4lg\n",
	LatentHeat,SensibleHeat,AdvectedEnergy);
*****/


  RefreezeEnergy *= Dt * SECPHOUR;

  /* if RefreezeEnergy is positive it means energy is available to melt the
     intercepted snow in the canopy.  If it is negative, it means that 
     intercepted water will be refrozen */
  
  /* Update maximum water interception storage */
  
  MaxWaterInt = LIQUID_WATER_CAPACITY * (*IntSnow) + MaxInt;

  /* Convert the vapor mass flux from a flux to a depth per interval */
  *VaporMassFlux *= Dt * SECPHOUR;
  
  if (RefreezeEnergy > 0.0) {

    if (-(*VaporMassFlux) > *IntRain) {
      *VaporMassFlux = -(*IntRain);
      *IntRain = 0.;
    }
    else
      *IntRain += *VaporMassFlux;

    PotSnowMelt = min((RefreezeEnergy/Lf/RHO_W), *IntSnow);

    *MeltEnergy -= (Lf * PotSnowMelt * RHO_W) / (Dt *SECPHOUR);
    
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
    
    if (RefreezeEnergy > - (*IntRain) * Lf) {
      *IntSnow += fabs(RefreezeEnergy) / Lf;
      *IntRain -= fabs(RefreezeEnergy) / Lf;

      *MeltEnergy += (fabs(RefreezeEnergy) * RHO_W) / (Dt *SECPHOUR);

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

      *MeltEnergy += (Lf * *IntRain * RHO_W) / (Dt *SECPHOUR);
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

}
