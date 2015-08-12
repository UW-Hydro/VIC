# Solar Radiation

For the past few years, the VIC hydrologic model has used a parameterization of incoming shortwave radiation based on the daily range in temperature, as described by Bristow and Campbell (1984). As we have implemented it, this formulation has been essentially 'calibrated' for the Pacific Northwest. The move to continental, and even global, modeling raises questions about how well the daily temperature range indexes solar radiation over large regions. Thornton and Running (1999) proposed a revised algorithm, with sensitivity to a wider climate range when compared with 40 synoptic stations in the continental U.S. This note describes a comparison between two variations of the Bristow and Campbell (1984) algorithm, with two variations of the Thornton and Running (1999) algorithm for climatic extremes in the tropics and (sub) Arctic.

## Bristow and Campbell Algorithm

Total transmittance T<sub>t</sub> is defined by Bristow and Campbell (1984) as the ratio of daily measured irradiance to daily extraterrestrial insolation. The algorithm for total transmittance, T<sub>t</sub> as a function of ![](radiation_images/Image32.gif) has the form:

![](radiation_images/Image35.gif)is the monthly mean daily temperature range (Bristow and Campbell 1984). The coefficient A, which represents the total transmittance on a cloudless day, can be further parameterized as a function of elevation (John Kimball, personal communication), as A=A<sub>msl</sub>+A<sub>lapse</sub>\*elevation.

The performance of two different implementations of the Bristow and Campbell algortihm was compared:

## VIC Parameters

For historic VIC applications, the parameters A<sub>msl</sub>, A<sub>lapse,</sub> B and C were held constant at 0.6, 0.0000295 /meter, 0.003 and 2.4\. These values were optimized for application in the Columbia River Basin, located in the Pacific Northwest (Kimball, personal comunication).

## B&C Parameters:

For implementation of the original Bristow and Campbell algorithm, A<sub>msl</sub>, A<sub>lapse,</sub> and C were held constant at 0.715, 0.0000295 /meter and 2.4\. B was calculated according to equation 2\. In addition, daily ![](radiation_images/Image32.gif)was reduced by 25% on rainy days and on the day before rain occurred if![](radiation_images/Image32.gif)changed by more than 2<sup>o</sup> C.

## Thornton and Running Algorithm

Thornton and Running have altered the basic form of the Bristow and Campbell algorithm in order to provide a minimum bound on transmittance:

They have also expanded the original model to include the effects of pressure, solar zenith angle and vapor pressure on the clear-sky transmittance (A):

where ![](radiation_images/Image39.gif)is the instantaneous transmittance at the reference elevation, at nadir, for a dry atmosphere, m<sub>theta</sub>is the optical air mass for a given zenith angle, theta, P is surface air pressure, e is vapor pressure (Pa) and alpha is a slope parameter (Pa<sup>-1</sup>) describing the influence of e on ![](radiation_images/Image43.gif). B has also been expanded as a three-parameter exponential decay curve:

## MTCLIM parameters:

The revised solar radiation algorithm has been incorporated into MTCLIM version 4.2\. Due to the dependence on vapor pressure, this algorithm requires input of dewpoint temperature. If T<sub>dew</sub> is not available, it is estimated based on the method of Kimball et al. (1997). Since the Kimball method also requires radiation input, a solution is found by iteration.

Again, the performance of two implementations of the algorithm was compared. One implementation with observed T<sub>dew</sub> and one implementation using the iterative procedure estimating both T<sub>dew</sub> and radiation.

MTCLIM 4.2 was used unmodified, with ![](radiation_images/Image45.gif)and![](radiation_images/Image46.gif)

In general, the MTCLIM formulation provides a tighter fit for both the high-latitude and low-latitude sites, although the improvement appears to be greater in the tropics. There is not a large difference in results when T<sub>dew</sub> is not provided as input. Based on this analysis, the VIC hydrologic model will be updated to include MTCLIM 4.2 for the generation of incident solar radiation.

## References

Bristow, K.L. and G.S. Campbell (1984). On the relationship between incoming solar radiation and daily maximum and minimum temperature, _Agricultural and Forest Meteorology_ 31, 159-166.

Kimball, J.S., S.W. Running and R. Nemani (1997). An improved method for estimating surface humidity from daily minimum temperature, _Agricultural and Forest Meteorology_ 85, 87-98.

Thornton, P.E. and S.W. Running (1999). An improved algorithm for estimating incident daily solar radiation from measurements of temperature, humidity and precipitation, _Agricultural and Forest Meteorology_ 93, 211-228.
