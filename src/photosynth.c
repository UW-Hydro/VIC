/********************************************************************************
  filename  : photosynth.c (imported from BETHY model)
  purpose   : Calculate photosynthesis, based on Farquhar (C3) and Collatz (C4)
              formulations
  interface : - input :
			- Ctype:          photosynthesis pathway (PHOTO_C3 or PHOTO_C4)
			- MaxCarboxRate:  maximum carboxlyation rate at 25 deg C       (mol(CO2)/m2 leaf area s)
			- MaxETransport:  maximum electron transport rate at 25 deg C  (mol(CO2)/m2 leaf area s) (C3 plants)
			- CO2Specificity: CO2 specificity at 25 deg C                  (mol(CO2)/m2 leaf area s) (C4 plants)
			- NscaleFactor:   nitrogen scaling factor at max carbox rate (Vm) and max electron transport rate (Jm)
			- Tfoliage:       vegetation temperature                       (deg C)
			- PIRRIN:         total irradiance at the surface              (mol photons/m2s)
			- aPAR:           absorbed photosynthetically active radiation (mol photons/m2 leaf area s)
			- Psurf:          near-surface atmospheric pressure            (Pa)
			- Catm:           CO2 mixing ratio in the atmosphere           (mol(CO2)/mol(air))
			- mode:           'ci': take Ci as an input, and compute photosynthesis and rs
					  'rs': take rs as an input, and compute photosynthesis and Ci

              - input or output, depending on mode:
			- rs:             stomatal resistance (per leaf area) (s/m)               (mode = 'rs': input; mode = 'ci': output)
			- Ci:             leaf-internal CO2 mixing ratio      (mol(CO2)/mol(air)) (mode = 'ci': input; mode = 'rs': output)

              - output:
			- Rdark:          'dark' respiration                  (mol(CO2)/m2 leaf area s)
			- Rphoto:         photorespiration                    (mol(CO2)/m2 leaf area s)
			- Agross:         gross assimilation (photosynthesis) (mol(CO2)/m2 leaf area s)
  programmer: Ted Bohn
  date      : October 20, 2006
  changes   :
  references: 
********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vicNl.h>

static char vcid[] = "$Id: $";

void photosynth(char    Ctype,
                double  MaxCarboxRate,
                double  MaxETransport,
                double  CO2Specificity,
                double  NscaleFactor,
                double  Tfoliage, 
                double  PIRRIN, 
                double  aPAR, 
                double  Psurf, 
                double  Catm, 
                char   *mode, 
                double *rs, 
                double *Ci, 
                double *Rdark, 
                double *Rphoto, 
                double *Agross) 
{
  extern option_struct options;
  double T;
  double T1;
  double T0;
  double Vcmax;
  double KC;
  double KO;
  double gamma;
  double Jmax;
  double K;
  double JE;
  double JC;
  double J0;
  double J1;
  double K1;
  double K2;
  double W1;
  double W2;
  double r0;
  double B;
  double C;
  double tmp;

  T1 = 25 + KELVIN;
  T  = Tfoliage + KELVIN;      // Canopy or Vegetation Temperature in Kelvin
  T0 = T - T1;                 // T relative to 25 degree Celsius, means T - 25

  /********************************************************************************
  ! Vcmax and Jmax are not only temperature dependent but also differ inside the canopy.
  ! This is due to the fact that the plant distributes its Nitrogen content and
  ! therefore Rubisco content so that, the place with the most incoming light got the
  ! most Rubisco. Therefore, it is assumed that the Rubisco content falls
  ! exponentially inside the canopy. This is reflected directly in the values of Vcmax
  ! and Jmax at 25 Celsius (Vcmax * nscl),  Knorr (107/108)
  ********************************************************************************/
  Vcmax = MaxCarboxRate * NscaleFactor * exp(EV*(T0/T1)/(Rgas*T));

  /********************************************************************************
  ! Determine temperature-dependent rates, compensation point, and 'dark' respiration
  ********************************************************************************/
  if (Ctype == PHOTO_C3) {

    /********************************************************************************
    ! C3 Plants
    ********************************************************************************/
    /********************************************************************************
    ! Temperature-dependent rates
    !
    ! Rate (with Aktivationenergy) vegetation temperature dependence is
    !  k = k(25C) * exp((Veg Temp - 25) * aktivation energy
    !                    / 298.16 * Rgas * (Veg Temp + KELVIN))
    ! => k = k0 * exp( T0 * E / T1 / Rgas / T ),  WHERE Rgas is the gas constant (8.314)
    ! This holds for OX-Oxygen partial pressure, KC-Michaelis-Menten constant for CO2,
    !    KO-Michaelis-Menten constant for O2, Vcmax-carboxylation capacity,
    !    Rdark-Dark respiration, K-PEPcase CO2 specivity
    !    Knorr (106)
    ********************************************************************************/
    KC = KC0 * exp(EC*(T0/T1)/(Rgas*T));
    KO = KO0 * exp(EO*(T0/T1)/(Rgas*T));
    /********************************************************************************
    ! CO2 compensation point without leaf respiration, gamma* is assumed to be linearly
    ! dependent on vegetation temperature, gamma* = 1.7 * TC (if gamma* in microMol/Mol)
    ! Here, gamma in Mol/Mol,       Knorr (105)
    ********************************************************************************/
    gamma = 1.7E-6 * Tfoliage;
    if (gamma < 0) gamma = 0;
    /********************************************************************************
    ! The temperature dependence of the electron transport capacity follows
    ! Farquhar(1988) with a linear temperature dependence according to the vegetation
    ! temperature
    !  J = J(25C) * TC / 25 WHERE J(25) = J0 * NscaleFactor
    ! minMaxETransport=1E-12
    ********************************************************************************/
    Jmax = MaxETransport * NscaleFactor * Tfoliage/25.;
    if (Jmax < minMaxETrans) Jmax = minMaxETrans;
    if ( Jmax > minMaxETrans)
       J1 = ALC3 * aPAR * Jmax / sqrt(Jmax*Jmax + (ALC3*aPAR)*(ALC3*aPAR));
    else
       J1 = 0.;

    /********************************************************************************
    !  Compute 'dark' respiration
    ! Following Farquhar et al. (1980), the dark respiration at 25C is proportional
    ! to Vcmax at 25C, therefore Rdark = const * Vcmax, but the temperature dependence
    ! goes with ER (for respiration) and not with EV (for Vcmax)
    ********************************************************************************/
    *Rdark = FRDC3 * MaxCarboxRate * NscaleFactor * exp(ER*(T0/T1)/(Rgas*T)) * hiTinhib(Tfoliage) * darkinhib(PIRRIN);


  }
  else if (Ctype == PHOTO_C4) {

    /********************************************************************************
    ! C4 Plants
    ********************************************************************************/
    /********************************************************************************
    ! Temperature-dependent rates
    !
    ! For C4 plants the Farquhar equations are replaced by the set of equations of
    !  Collatz et al. 1992:
    !  A = min{JC, JE} - Rdark
    !  JC = k * CI
    !  JE = 1/2/Theta *[Vcmax + Ji - sqrt((Vcmax+Ji)^2 - 4*Theta*Vcmax*Ji)]      with
    !  Ji = alphai * Ipar / Epar with Ji=aPAR in Mol(Photons)
    !        Knorr (114a-d)
    !  alphai is the integrated quantum efficiency for C4 plants (ALC4 = 0.04,
    !    compared to the efficiency of C3 plants, ALC3 = 0.28)
    !  Theta is the curve PARAMETER (0.83) which makes the change between
    !   Vcmax and K limitation smooth
    !  K is the PECase CO2 specifity instead of the electron transport capacity
    !   within C3 plants
    !  Ci is the stomatal CO2 concentration = Cimin + (Ci0 - Cimin)* GC/GC0 with
    !    Cimin = 0.3 and 0.15 Catm respectivly (Catm is the CO2 mixing ratio)
    !
    ! The factor 1E3 comes that Jmax for C3 is in microMol and K is in milliMol,
    !   which is not considered in INITVEGDATA
    ! K scales of course with EK
    ********************************************************************************/
    K = CO2Specificity * 1.E3 * NscaleFactor * exp(EK*(T0/T1)/(Rgas*T));
    /********************************************************************************
    !  Compute 'dark' respiration
    !  same as C3, just the 25 degree Celsius proportional factor is different
    !    0.011 for C3,  0.0042 for C4
    ********************************************************************************/
    *Rdark = FRDC4 * MaxCarboxRate * NscaleFactor * exp(ER*(T0/T1)/(Rgas*T)) * hiTinhib(Tfoliage) * darkinhib(PIRRIN);

  } // End computation of T-dependent rates and dark respiration

  if (!strcasecmp(mode,"ci")) {

    /********************************************************************************
    ! If Ci given, compute gross photosynthesis components at given leaf-internal CO2
    ********************************************************************************/
    if (Ctype == PHOTO_C3) {

      /********************************************************************************
      ! C3 Plants
      ********************************************************************************/
      /********************************************************************************
      !  The assimilation follows the Farquhar (1980) formulation for C3 plants
      !  A = min{JC, JE} - Rdark
      !  JC = Vcmax * (Ci - gamma) / (Ci + KC * (1 + OX/KO))
      !  JE = J * (Ci - gamma) / 4 / (Ci + 2 * gamma)      with
      !   J = alpha * I * Jmax / sqrt(Jmax^2 + alpha^2 * I^2) with I=aPAR in Mol(Photons)
      !        Knorr (102a-c, 103)
      !  Here J = J1 and A is the gross photosynthesis (Agross), i.e. still including the
      !          respiratory part Rdark
      ********************************************************************************/
      JE = J1 * ((*Ci) - gamma) / 4. / ((*Ci) + 2. * gamma);
      JC = Vcmax * ((*Ci) - gamma) / ( (*Ci) + KC * (1. + OX / KO) );

    }
    else if (Ctype == PHOTO_C4) {

      /********************************************************************************
      ! C4 Plants
      ********************************************************************************/
      /********************************************************************************
      !  JE = 1/2/Theta *[Vcmax + Ji - sqrt((Vcmax+Ji)^2 - 4*Theta*Vcmax*Ji)]
      !    Ji = ALC4 * aPAR
      !  J0 is the sum of the first two terms in JE
      ********************************************************************************/
      J0 = (ALC4 * aPAR + Vcmax) /  2. / THETA;
      /********************************************************************************
      !  last 2 terms:  with J0^2 = 1/4/Theta^2*(Vcmax+Ji)^2
      !       sqrt(1/4/Theta^2)*sqrt((Vcmax+Ji)^2 - 4*Theta*Vcmax*Ji))
      !   = sqrt (J0^2 - Vcmax*Ji/Theta)
      ********************************************************************************/
      JE = J0 - sqrt (J0*J0 - Vcmax * ALC4 * aPAR / THETA);
      /********************************************************************************
      !         see above
      ********************************************************************************/
      JC = K * (*Ci);

//fprintf(stderr,"other mode\nJE %f JC %f\n",JE,JC);
    }

  } // End computation of gross photosynthesis components at given Ci
  else {

    /********************************************************************************
    ! If rs given, compute gross photosynthesis components at given stomatal resistance
    ********************************************************************************/
    if (Ctype == PHOTO_C3) {

      /********************************************************************************
      ! C3 Plants
      ********************************************************************************/
      /********************************************************************************
      !  Remember:
      !  A = min{JC, JE} - Rdark
      !  JC = Vcmax * (Ci - gamma) / (Ci + KC * (1 + OX/KO))
      !  JE = J * (Ci - gamma) / 4 / (Ci + 2 * gamma)      with
      !   J = alpha * I * Jmax / sqrt(Jmax^2 + alpha^2 * I^2) with I=aPAR in Mol(Photons)
      !        Knorr (102a-c, 103)
      ! J = J1
      ********************************************************************************/
      
      /********************************************************************************
      !         Helping friends K1, W1, W2, K2
      ********************************************************************************/
      K1 = 2. * gamma;
      W1 = J1 / 4.;
      W2 = Vcmax;
      K2 = KC * (1. + OX / KO);
      /********************************************************************************
      ! A = gs / 1.6 * (Catm - Ci) * Psurf / Rgas / T
      ! <=> Ci = Catm - 1.6 * Rgas * T / Psurf / gs * A = Catm - A / G0
      ! Let rs = 1/gs, where gs = stomatal conductance
      ! and r0 = 1/G0
      ! So Ci = Catm - 1.6*(Rgas*T/Psurf)*rs * A = Catm - A * r0
      ********************************************************************************/
      r0 = (*rs) * 1.6*Rgas*T/Psurf;
      /********************************************************************************
      ! A = min{JC, JE} - Rdark
      ! => A = JC - Rdark OR A = JE - Rdark
      ! Set this (A =) in Ci formula above
      ! Set Ci in
      !  JE = J * (Ci - gamma) / 4 / (Ci + 2 * gamma)
      ! => quadratic formula in JE
      ! 0 = JE^2 -(Rdark+J/4+(Catm+2*gamma)/r0)*JE +J/(4*r0)*(Catm-gamma)+J/4*Rdark
      ********************************************************************************/
      B = (*Rdark) + W1 + (Catm + K1)/r0;
      C = W1*(Catm - gamma)/r0 + W1*(*Rdark);
      /********************************************************************************
      ! with 0 = x^2 + bx + c
      !      x1/2 = -b/2 +/- sqrt(b^2/4 - c)
      !  take '-' as minimum value of formula
      ********************************************************************************/
      tmp = B*B/4 - C;
      if (tmp < 0) tmp = 0;
      JE = B/2. - sqrt(tmp);
//fprintf(stderr,"Rdark %.10f W1 %.10f Catm %.10f K1 %.10f rs %.10f r0 %.10f gamma %.10f B %.10f C %.10f JE %.10f\n",*Rdark,W1,Catm,K1,*rs,r0,gamma,B,C,JE);
      /********************************************************************************
      ! Set Ci in
      !  JC = Vcmax * (Ci - gamma) / (Ci + KC * (1 + OX/KO))
      ! WRITE JC = Vcmax * (Ci - gamma) / (Ci + K2)
      ! => quadratic formula in JC
      ! 0 = JC^2 -(Rdark+Vcmax+(Catm+K2)/r0)*JC +Vcmax/r0*(Catm-gamma)+Rdark*Vcmax
      ********************************************************************************/
      B = (*Rdark) + W2 + (Catm + K2)/r0;
      C = W2*(Catm - gamma)/r0 + W2*(*Rdark);
      tmp = B*B/4 - C;
      if (tmp < 0) tmp = 0;
      JC = B/2. - sqrt(tmp);
//fprintf(stderr,"Rdark %.10f W2 %.10f Catm %.10f K2 %.10f rs %.10f r0 %.10f gamma %.10f B %.10f C %.10f JC %.10f\n",*Rdark,W2,Catm,K2,*rs,r0,gamma,B,C,JC);

    }
    else if (Ctype == PHOTO_C4) {

      /********************************************************************************
      ! C4 Plants
      ********************************************************************************/
      /********************************************************************************
      ! Recall:
      !  Collatz et al. 1992:
      !  A = min{JC, JE} - Rdark
      !  JC = k * Ci
      !  JE = 1/2/Theta *[Vcmax + Ji - sqrt((Vcmax+Ji)^2 - 4*Theta*Vcmax*Ji)]      with
      !   Ji = alphai * Ipar / Epar with Ji=aPAR in Mol(Photons)           and
      !   aPAR = Ipar / Epar;  ALC4=alphai; J0=1/2/Theta *(Vcmax + Ji);
      !   Ci = Catm - 1.6 * Rgas * T / Psurf / gs * A = Catm - A / G0                and
      !   => A = JC - Rdark OR A = JE - Rdark
      ! Let rs = 1/gs, where gs = stomatal conductance
      ! and r0 = 1/G0
      ! So Ci = Catm - 1.6*(Rgas*T/Psurf)*rs * A = Catm - A * r0
      ********************************************************************************/
      r0 = (*rs) * 1.6*Rgas*T/Psurf;
      /********************************************************************************
      !  J0=1/2/Theta *(Vcmax + Ji) = (alphai * aPAR + Vcmax) / 2 / Theta
      ********************************************************************************/
      J0 = (ALC4 * aPAR + Vcmax) /  2. / THETA;
      /********************************************************************************
      !  JE = J0 - sqrt( J0^2 - Vcmax*alphai*aPAR/Theta)
      ********************************************************************************/
      JE = J0 - sqrt (J0*J0 - Vcmax * ALC4 * aPAR / THETA);
      /********************************************************************************
      !  JC = (Catm/r0 + Rdark) / (1 + 1/(K*r0))
      ********************************************************************************/
      JC = (Catm/r0 + (*Rdark)) / (1. + 1/(K*r0));

    }

  } // End computation of gross photosynthesis components at given rs

  /********************************************************************************
  ! Compute gross assimilation (photosynthesis)
  ********************************************************************************/
  if (JE < JC) {
    /* light limitation */
    *Agross = JE*hiTinhib(Tfoliage);
  }
  else {
    /* CO2 limitation */
    *Agross = JC*hiTinhib(Tfoliage);
  }

  /********************************************************************************
  ! If rs given, compute leaf-internal CO2 concentration
  ********************************************************************************/
  if (!strcasecmp(mode,"rs")) {

    /********************************************************************************
    ! A = gs / 1.6 * (Catm - Ci) * p / Rgas / T
    ! <=> Ci = Catm - 1.6 * Rgas * T / p / gs * A = Catm - A / G0
    ! Let rs = 1/gs, where gs = stomatal conductance
    ! and r0 = 1/G0
    ! So Ci = Catm - 1.6*(Rgas*T/Psurf)*rs * A = Catm - A * r0
    !   with A = Net assimilation = NPP
    ! (Catm is the CO2 mixing ratio)
    ********************************************************************************/
    if (r0 > 1.e6) r0 = 1.e6;
    *Ci = Catm - ((*Agross)-(*Rdark))*r0;
    if (*Ci < 0) *Ci = 0;

  }
//fprintf(stderr,"JE %.10f JC %.10f Agross %.10f Ci %.10f\n",JE*1e12,JC*1e12,*Agross*1e12,*Ci);

  /********************************************************************************
  ! Compute photorespiration
  ********************************************************************************/
  if (Ctype == PHOTO_C3) {

    /********************************************************************************
     ! Photorespiration for C3 plants: Carboxylation controlled assimilation
     ! JC = Assimilation - Photorespiration
     ! JC = Vcmax * (Ci - gamma) / (Ci + K2)
     ! Photorespiration = Vcmax * gamma / (Ci + K2)
    ********************************************************************************/
     *Rphoto = Vcmax * gamma / ( (*Ci) + KC * ( 1. + OX / KO ) ) * hiTinhib(Tfoliage);

  }
  else {

    /********************************************************************************
     ! Photorespiration is 0 for C4 plants
    ********************************************************************************/
    *Rphoto = 0.;

  }

  /********************************************************************************
  ! If ci given, compute stomatal resistance
  ********************************************************************************/
  if (!strcasecmp(mode,"ci")) {

    /********************************************************************************
    ! Diffusion equation Flux = (Catm - CI) / resistence, rs
    !   conductance gs = 1 / rs  =>  Flux = (Catm-CI) * gs
    !   (Catm ... CO2mixingRatio)
    !   Flux of CO2 is Agross * amount, Assimilation rate * amount
    !   Agross is here, Gross Assimilation, though Agross-Rdark = (net) Assimilation rate
    !   the amount comes from the ideal gas equation pV=nRgasT => n/V = p / RgasT
    !   the stomatal conductance for CO2 is less { the conductance of H2O by
    !   the factor of 1.6: gs(CO2) = gs(H2O) / 1.6, due to its lower mobiblity due
    !   to its higher mass
    !   => A (net) = gs/1.6 * (Catm-CI) * p/RgasT
    !   => gs = A(net)*1.6*RgasT/p/(Catm-CI)
    ********************************************************************************/
    if ((*Agross)-(*Rdark) < SMALL) {
      *rs = HUGE_RESIST;
    }
    else {
      *rs = 0.625*(Catm-(*Ci))/((*Agross)-(*Rdark))*(Psurf/(Rgas*T));
    }
    if (*rs > HUGE_RESIST) *rs = HUGE_RESIST;

  }

}

/********************************************************************************
! High temperature inhibition
********************************************************************************/
/********************************************************************************
! FUNCTION which inhibits Assimilation and Respiration at temperatures above
! 55 Celsius from Collatz et al., Physiological and environmental regulation
! of stomatal conductance, photosynthesis and transpiration: a model that
! includes a laminar boundary layer,
! Agricultural and Forest Meteorology, 54, pp. 107-136, 1991
********************************************************************************/

double hiTinhib(double T) {

  double hiTinhib;

  // T = Vegetation temperature in degrees Celsius
  hiTinhib = 1. / ( 1. + exp( 1.3 * ( T - 55. ) ) );

  return hiTinhib;

}

/********************************************************************************
! Dark inhibition
********************************************************************************/
/********************************************************************************
! FUNCTION which inhibits Dark-Respiration in light
! after Brooks and Farquhar, Effect of temperature on the CO2/O2 specifity on RuBisCO
! and the rate of respiration in the light, Planta 165, 397-406, 1985
!
! It is fitted to inhibit the dark-respiration to 50% of its uninhibited value
! up from 50 umol/m^2s.
********************************************************************************/

double darkinhib(double IRR) {

  double darkinhib;

  // IRR = Total irridiance at the surface in mol/m^2s
  if ( IRR < 0 ) {
     darkinhib = 0.;
  }
  else {
     darkinhib = 0.5 + 0.5*exp(- IRR * 1.e6 / 10. );
  }

  return darkinhib;

}

