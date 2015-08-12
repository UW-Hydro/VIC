# VIC Model Calibration Using Soil Parameters

There are five soil parameters that comprise the set of variables to rely on for model calibration. More detailed descriptions of each parameter are included in the [VIC Soil Parameter File documentation](SoilParam.md). These parameters, their typical ranges, and a short description of the effect of each one on simulated hydrographs follows.

1.  Ds - [>0 to 1] This is the fraction of Ds<sub>max</sub> where non-linear (rapidly increasing) baseflow begins. With a higher value of Ds, the baseflow will be higher at lower water content in lowest soil layer.
2.  Ds<sub>max</sub> - [>0 to ~30, depends on hydraulic conductivity] This is the maximum baseflow that can occur from the lowest soil layer (in mm/day).
3.  Ws - [>0 to 1] This is the fraction of the maximum soil moisture (of the lowest soil layer) where non-linear baseflow occurs. This is analogous to Ds. A higher value of Ws will raise the water content required for rapidly increasing, non-linear baseflow, which will tend to delay runoff peaks.

4.  b<sub>inf</sub> - [>0 to ~0.4] This parameter defines the shape of the Variable Infiltration Capacity curve. It describes the amount of available infiltration capacity as a function of relative saturated gridcell area. A higher value of b<sub>inf</sub> gives lower infiltration and yields higher surface runoff.
5.  Soil Depth (of each layer) - [typically 0.1 to 1.5 meters] Soil depth effects many model variables. In general, for runoff considerations, thicker soil depths slow down (baseflow dominated) seasonal peak flows and increase the loss due to evapotranspiration.

Before calibration, it is important to understand the dynamics of how baseflow is calculated. Baseflow at low soil moisture content is calculated from

![](Bernt/calibration/Image7.gif).

In Figure 1 four situations of water content in lower layer vs baseflow is plotted. Case A showes a standard situation where Ds = 0.100, Dsmax = 10.0 mm/day, Ws = 0.60 and maximum water content in lower l ayer (W2max) is 176 mm. We can clearly see that the baseflow becomes non linear when the moisture content in lower layer exceeds ca 105 mm. Before this the baseflow is between 0 and 1 mm/timestep. At higher water content in lower layer the baseflow rapidly increases towards maximum baseflow at 10 mm/t imestep.

If a basin have a typical snow accumulation with a following dominant snowmelt period, the baseflow parameters will be very sensitive. Baseflow in January will be low because of the low water content in lower layer. When snowmelt starts water infiltrates and starts to fill the soil layers. At a given point (Ws) the water content in lower layer will be higher than the point where non linear baseflow occurs and baseflow will increase rapidly.

![](Bernt/calibration/Image8.gif)

Figure 1: Baseflow Dynamics

In case B the Ds have changed to 0.300\. The effect of this is that baseflow will be higher at lower watercontent in lower layer, the bottom layer will fill slower in the spring and nonlinear baseflow occurs later in the season. In case C the Ws is changed to 0.85\. The result is that water content in bottom layer have to exceed 155 mm before non linear baseflow occurs, thus the spring peak will be delayed. The last case is D and this sh owes parameters used for the Waneta subbasin in the Columbia River Basin. We needed to slow down the spring snowmelt peak and increase the evapotranspiration and therefor we made soils thick, Ds litle and Ws at a moderate level (Ws = 0.7). If you want to play with the baseflow parameters there is a spreadsheet available that can be used to do calculations and plots for different situations.

After each model run it is a good thing to plot up the soil moisture content for a few gridcells and see how it works. Also plot the timeseries of t he SWE for some years. Make sure that the different parts of the water balance is right. Evaporation, runoff and precip have values that are close to w hat you would expect. There is a program available to calculate the average basin timeseries for precip, evap, runoff, baseflow, moister in all layers. Do a timeseries plot of this data and see how the soil moisture content varies over the season.
