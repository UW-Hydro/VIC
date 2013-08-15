/********************************************************************************
  filename  : faparl.c (imported from BETHY model)
  purpose   : Compute each canopy layer's leaf area index (LAIlayer) and
              absorbed photsynthetic active radiation (aPAR) from total LAI
              and direct irradiation.
  interface : - input :
			- CanopLayerBnd:         Array of upper boundaries of canopy layers,
                                       expressed as fractions of total LAI
			- LAItotal:    Total canopy LAI (m2 leaf area / m2 ground area)
			- AlbSoilPAR:  Soil albedo in PAR wavelength range
			- CosZen:      Cosine of solar zenith angle
			- Fdir:        Fraction of incoming solar radiation that
                                       is direct (as opposed to diffuse)
              - output:
                        - LAIlayer:    LAI per canopy layer (m2 leaf area / m2 ground area)
                        - aPAR:        The absorbed PAR (W / m2 leaf area)
                                       per canopy layer, normalized for PAR = 1 W

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

void faparl(double  *CanopLayerBnd,
            double   LAItotal,
            double   AlbSoilPAR,
            double   CosZen,
            double   Fdir,
            double  *LAIlayer,
            double  *aPAR)
{
  extern option_struct options;
  double         FC;
  int            cidx;
  double         ZH;
  double         ZP1;
  double         ZP0;
  double         K0;
  double         X0;
  double         X1;
  double         X2;
  double         Q0;
  double         Q1;
  double         F;
  double         B0;
  double         B1;
  double         B4;
  double         EKL0;
  double         EHL0;
  double         EKL;
  double         EHL;

  /*---------------------------------------------------------------------------
  ! Compute fractional vegetation cover per PFT
  ! NOTE: not the same as the area fraction of the vegetation tile;
  ! here we mean the amount of "clumping" of the vegetation, or
  ! how much bare ground is visible between the plants.
  ! This will be used to relate fluxes per leaf area to fluxes per ground area.
  !---------------------------------------------------------------------------*/
  if (LAItotal < LaiLimit)
    FC = LAItotal / LaiLimit * FcMax;
  else
    FC = FcMax;
  if (FC < FcMin)
    FC = FcMin;

  /*---------------------------------------------------------------------------
  ! Initialize ABSORBED PAR PER LEAF AREA
  ! Partition LAItotal proportionally among the canopy layers
  !---------------------------------------------------------------------------*/
  for (cidx = 0; cidx < options.Ncanopy; cidx++) {
    aPAR[cidx] = 0.;
    if (cidx > 0)
      LAIlayer[cidx] = LAItotal * (CanopLayerBnd[cidx] - CanopLayerBnd[cidx-1]);
    else
      LAIlayer[cidx] = LAItotal * CanopLayerBnd[cidx];
    if (LAIlayer[cidx] < LaiMin)
      LAIlayer[cidx] = LaiMin;
  }

  /*---------------------------------------------------------------------------
  !     COMPUTE ABSORBED PAR PER LEAF AREA
  !---------------------------------------------------------------------------*/
  if (CosZen >= ZenithMinPar) {

    /*-------------------------------------------------------------------
    ! The Absorbed Par per leaf Area which is used later for the Net Assimilation,
    ! is calculated via the two stream approximation of Sellers (1985):
    !  muq * d(Rdown)/dl + (1-(1-b)*omega)*Rdown - omega*b*Rup   = omega*muq*k*(1-b0)*R
    ! -muq * d(Rup)/dl   + (1-(1-b)*omega)*Rup   - omega*b*Rdown = omega*muq*k*b0*R
    !  with
    !   Rdown - downwards diffusive flux
    !   Rup   - upwards diffusive Flux
    !   R     - direct flux, with R(0) = dPAR * RPAR with RPAR incoming irridiance in PAR
    !           and R = R0 * EXP(-kl) with R0=R(0) - exponential extinction law of
    !           Monsi and Racki (1953)
    !           or similar Lambert-Beer's law
    !   b     - forward scatter fraction of diffusive flux
    !   b0    - forward scatter fraction of direct flux
    !   k     - extinction coefficient for the direct flux
    !   muq = int(1/K(mu))dmu|_0^1 - the integral of 1/K over the downward hemisphere
    !   omega - the single leaf albedo
    !
    !  The general solutions are (kl,hl=k*l,h*l):
    !  Rup   = q2*R0*EXP(-kl) + p1*B1*EXP(hl) + p2*B2*EXP(-hl)
    !  Rdown =-q1*R0*EXP(-kl) +    B1*EXP(hl) +    B2*EXP(-hl)
    !   with
    !    h  = sqrt( (1-(1-b)*omega)^2/muq^2 - omega^2*b^2/muq^2 )
    !    p1 = ( (1-(1-b)*omega) + muq * h )/omega/b
    !    p2 = ( (1-(1-b)*omega) - muq * h )/omega/b
    !   -q1 = (omega^2*b*muq*k* (1-b0) + omega*    b0  *muq*k*(1-(1-b)*omega - muq*k))/
    !         ((1-(1-b)*omega)^2-muq^2*k^2-omega^2*b^2)
    !    q2 = (omega^2*b*muq*k*    b0  + omega* (1-b0) *muq*k*(1-(1-b)*omega + muq*k))/
    !         ((1-(1-b)*omega)^2-muq^2*k^2-omega^2*b^2)
    !    B1/B2 from boundary conditions
    !-------------------------------------------------------------------
    !  Make two assumptions:
    !  1) the distribution of leaf angles is isotropic
    !  2) the leaf reflectivity and transmissivity are equal (the sum = omega)
    !  => b=0.5, b0=0.5, k=0.5/mu with mu=cos(theta) solar zenith angle => muq=1
    !
    !  => k  = 1/2/mu
    !     h  = sqrt( 1 - omega )
    !     p1 = ( 1-omega/2 + h )/omega/2
    !     p2 = ( 1-omega/2 - h )/omega/2
    !   ! p2 = 1 / p1 !
    !     q1 = ( (1 + 2*mu)*omega/2 )/(1-4*mu^2*(1-omega)) = ( k*(k + 1)*omega/2 )/
    !                                                        (k^2-1-omega)
    !     q2 = ( (1 - 2*mu)*omega/2 )/(1-4*mu^2*(1-omega)) = ( k*(k - 1)*omega/2 )/
    !                                                        (k^2-1-omega)
    !
    ! Determine B1 and B2 from the boundary conditions:
    !  1) Rdown(0) equals the incoming diffuse radiation
    !     Rdown(0) = (1-dPAR)*RPAR
    !  => Rdown(0) + R(0) = (1-dPAR)*RPAR + dPAR*RPAR = RPAR as total incoming PAR
    !  2) the reflection at the lower boundary of the canopy is given by the soil
    !     reflectance
    !     Rup(LAI) = AlbSoilPAR * (R(LAI) + Rdown(LAI))
    !  Here: faparl gets AlbSoilPAR as Variable AlbSoilPAR, LAI is the total canopy LAI, LAItotal
    !
    !  => B1 = + ( eta*R0 - (Rd+q1*R0) * gamma1 )/(gamma1 - gamma2)
    !     B2 = - ( eta*R0 - (Rd+q1*R0) * gamma2 )/(gamma1 - gamma2)
    !   with
    !     eta    = AlbSoilPAR * (1-q1)-q2) * EXP(-k*LAI)
    !     gamma1 = ( p1 - AlbSoilPAR) * EXP( + h*LAI)
    !     gamma2 = ( p2 - AlbSoilPAR) * EXP( - h*LAI)
    !     Rd     = Rdown(0) = (1-dPAR)*RPAR
    !------------------------------------------------------------------
    ! THAT IS THE COMPLETE SOLUTION OF THE TWO STREAM APPROXIMATION UNDER THE BOUNDARY
    ! CONDITIONS AND ASSUMPTIONS MENTIONED ABOVE !!!!!!!!!!!!!!!!!!
    !
    ! Therefore, the absorbed Radiation inside the canopy is:
    !   aPAR = -d/dl ( R(l) + Rdown - Rup)
    !        = (1-q1-q2)*k*R0*EXP(-kl) - (1-p1)*h*B1*EXP(hl) + (1-p2)*h*B2*EXP(-hl)
    ! But the absorbed PAR per canopy layer in discrete steps is:
    !   aPAR = 1/(D(cidx)-D(cidx-1)) * int(-d/dl(R(l) + Rdown - Rup))dl|_D(cidx)^D(cidx-1)
    !            = (R(D(cidx-1)) + Rdown(D(cidx-1)) - Rup(D(cidx-1)) - R(D(cidx)) + Rdown(D(cidx))
    !              - Rup(D(cidx)) / ((D(cidx)-D(cidx-1))
    !  and  R(l)+Rdown-Rup = (1-q1-q2)*  R0*EXP(-kl) + (1-p1)*  B1*EXP(hl) + (1-p2) *
    !                        B2*EXP(-hl)
    !------------------------------------------------------------------
    ! The clumping of the vegetation is taken into account, defining LAIc = LAI / fc
    ! as an effective LAI.
    ! Taken this into account, l = l/fc but the solutions stay as they are because
    ! of the differentiations are take d/dl according to the NEW l=l/fc
    ! Only aPAR has to be multiplied with fc at the END because D(cidx)-D(cidx-1) is still
    ! the old l
    !-----------------------------------------------------------------*/
    /*---------------------------------
    !  h = sqrt( 1 - omega )
    !---------------------------------*/
    ZH = sqrt(1. - OMEGA);
    /*---------------------------------
    !  p1 = ( 1-omega/2 + h )/omega/
    !---------------------------------*/
    ZP1 = (1. - OMEGA / 2. + ZH) / OMEGA * 2.;
    /*---------------------------------------
    ! p2 = ( 1-omega/2 - h )/omega/2 = 1 / p1
    !---------------------------------------*/
    ZP0 = 1. / ZP1;
    /*---------------------------------
    ! k = 0.5/mu
    !---------------------------------*/
    K0 = 0.5 / CosZen;
    if (K0 == ZH) K0 = K0 + 1E-12;
    if (K0 == -ZH) K0 = K0 + 1E-12;
    /*-------------------------------------------------------------
    ! denominator of q1 and q2
    !-------------------------------------------------------------*/
    X0 = (1. - 4. * CosZen*CosZen * ZH*ZH);
    /*-------------------------------------------------------------
    ! q1 = ( (1 + 2*mu)*omega/2 )/(1-4*mu^2*(1-omega))
    !    = ( k*(k + 1)*omega/2 )/(k^2-1-omega)
    !-------------------------------------------------------------*/
    Q1 = ((1. + 2. * CosZen) * OMEGA / 2.) / X0;
    /*-------------------------------------------------------------
    ! q2 = ( (1 - 2*mu)*omega/2 )/(1-4*mu^2*(1-omega))
    !    = ( k*(k - 1)*omega/2 )/(k^2-1-omega)
    !-------------------------------------------------------------*/
    Q0 = ((1. - 2. * CosZen) * OMEGA / 2.) / X0;
    /*----------------------------------------------------------
    ! EXP(-k*LAI/fc)
    !----------------------------------------------------------*/
    EKL = exp(-K0 / FC * LAItotal);
    /*-----------------------------------------------------------
    ! EXP(-h*LAI/fc)
    !----------------------------------------------------------*/
    EHL = exp(-ZH / FC * LAItotal);
    /*----------------------------------------------------------
    ! gamma1 = ( p1 - AlbSoilPAR) * EXP( + h*LAI)
    !----------------------------------------------------------*/
    X1 = (ZP1 - AlbSoilPAR) / EHL;
    /*----------------------------------------------------------
    ! gamma2 = ( p2 - AlbSoilPAR) * EXP( - h*LAI)
    !----------------------------------------------------------*/
    X0 = (ZP0 - AlbSoilPAR) * EHL;
    /*----------------------------------------------------------
    ! eta = AlbSoilPAR * (1-q1)-q2) * EXP(-k*LAI)
    !----------------------------------------------------------*/
    X2 = (AlbSoilPAR * (1. - Q1) - Q0) * EKL;
    /*-----------------------------------------------------------
    ! F = 1 - dPAR + dPAR * q1
    ! => F*RPAR = Rd + q1*R0
    ! i.e. calculation takes RPAR=1
    !----------------------------------------------------------*/
    F = 1. - Fdir + Q1 * Fdir;
    /*------------------------------------------------------------
    ! B1*RPAR = B1, B2*RPAR = B2, B4*RPAR = R(0) + Rdown(0) - Rup(0)
    !  B1 = + ( eta*R0 - (Rd+q1*R0) * gamma1 )/(gamma1 - gamma2)
    !-----------------------------------------------------------*/
    B1 = (X2 * Fdir - F * X0) / (X1 - X0);
    /*----------------------------------------------------------
    !  B2 = - ( eta*R0 - (Rd+q1*R0) * gamma2 )/(gamma1 - gamma2)
    !----------------------------------------------------------*/
    B0 = (X2 * Fdir - F * X1) / (X0 - X1);
    /*----------------------------------------------------------
    !  R(l)+Rdown-Rup = (1-q1-q2)* R0*EXP(-kl) + (1-p1)* B1*EXP(hl) + (1-p2)* B2*EXP(-hl)
    !---------------------------------------------------------------------------*/
    B4 = (1. - Q0 - Q1) * Fdir + (1. - ZP1) * B1 + (1. - ZP0) * B0;


    /*-----------------------------------------------------------
    ! Loop over all canopy layers except the final layer
    !-----------------------------------------------------------*/
    for (cidx=0; cidx<options.Ncanopy-1; cidx++) {
//      /*---------------------------------------------------------------------------
//      ! p2=1/p1
//      !---------------------------------------------------------------------------*/
//      ZP0 = 1. / ZP1;
      /*---------------------------------------------------------------------------
      ! EXP(-k*l/fc)
      ! with l=LAI*CanopLayerBnd(K), i.e. l is element of [0,LAI], i.e. 0,LAI/3,2LAI/3,LAI
      !---------------------------------------------------------------------------*/
      EKL0 =  exp(-K0 / FC * CanopLayerBnd[cidx] * LAItotal);
      /*---------------------------------------------------------------------------
      ! EXP(-h*l/fc)
      !---------------------------------------------------------------------------*/
      EHL0 = exp(-ZH / FC * CanopLayerBnd[cidx] * LAItotal);
      /*---------------------------------------------------------------------------
      ! R(D(cidx))+ Rdown(D(cidx))- Rup(D(cidx))=
      !      (1-q1-q2)*R0*EXP(-kl)+(1-p1)*B1*EXP(hl)+(1-p2)*B2*EXP(-hl)
      !  i.e. X0*RPAR = above
      !---------------------------------------------------------------------------*/
      X0 = (1. - Q0 - Q1) * EKL0 * Fdir + (1. - ZP1) * B1 / EHL0 + (1. - ZP0) * B0 * EHL0;
      /*---------------------------------------------------------------
      ! aPAR(cidx) = (R(D(cidx-1))+Rdown(D(cidx-1))-Rup(D(cidx-1)) - R(D(cidx))+
      ! Rdown(D(cidx))-Rup(D(cidx)) / ((D(cidx)-D(cidx-1))
      ! Here aPAR only the nominator; the division is made outside faparl
      !--------------------------------------------------------------*/
      aPAR[cidx] = B4 - X0;
      /*--------------------------------------------------------------
      ! R(D(cidx-1))+Rdown(D(cidx-1))-Rup(D(cidx-1) in next step
      !--------------------------------------------------------------*/
      B4 = X0;
    }   // # canopy layers - 1


    /*--------------------------------------------------------------
    ! Now the same for the last layer
    !--------------------------------------------------------------*/
    /*------------------------------------------------------------
    ! R(D(NL))+ Rdown(D(NL))- Rup(D(NL))=
    !      (1-q1-q2)*R0*EXP(-kLAI)+(1-p1)*B1*EXP(hLAI)+(1-p2)*B2*EXP(-hLAI)
    !  i.e. X0*RPAR = above for the lower boundary
    !-----------------------------------------------------------*/
    X0 = (1. - Q0 - Q1) * EKL * Fdir + (1. - ZP1) * B1 / EHL + (1. - ZP0) * B0 * EHL;
    /*-----------------------------------------------------------------
    ! aPAR(NL) = (R(D(NL-1))+Rdown(D(NL-1))-Rup(D(NL-1)) - R(D(NL))+
    ! Rdown(D(NL))-Rup(D(NL)) / ((D(NL)-D(NL-1))
    ! Here aPAR only the denominator; the division is made outside faparl
    !-----------------------------------------------------------------*/
    aPAR[options.Ncanopy-1] = B4 - X0;


    /*-----------------------------------------------------------
    ! Multiplication of aPAR with fc to allow for vegetation "clumping"
    !----------------------------------------------------------*/
    for (cidx=0; cidx<options.Ncanopy; cidx++) {
      aPAR[cidx] = aPAR[cidx]*FC;
    }


  }  // CosZen >= ZenithMinPar

}

