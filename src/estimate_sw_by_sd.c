/********************************************************************************
  filename  : estimate_sw_by_sd.c
  purpose   : Estimate daily short wave radiation by sunshine duration.
              通过日照时数估算日短波辐射量

  interface : - input :
                输入参数
			- sunshine duration (h)
			  日照时数
			- latitude (degree)
			  纬度
			- day sequence (1 - 365/366, 1 at Jan. 1st)
			  日序

              - output: - daily short wave radiation (W/m^2)
                输出：日短波辐射量
  programmer: Sibada, SCUT
  date      : April 27, 2016
  references: National Standard of People’s Republic of China GB/T 20481-2006.
              The grade of meteorological drought. Beijing: China Standard Press, 2006:11-12
              [中华人民共和国国家标准 GB/T 20481-2006 气象干旱等级. 北京: 中国标准出版社, 2006:11-12]
  changes   :

********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>

double estimate_sw_by_sd(double sunshine_duration,
                         double latitude,
                         int day_sequence)
{
    double Rs;                  /* daily short wave radiation (W/m^2)
                                            日短波辐射量 */
    double Ra;                  /* daily extraterrestrial Radiation (W/m^2)
                                            日地球外辐射 */
    double as;                  /* regression constant as
                                            回归常数 */
    double bs;                  /* regression constant bs
                                            回归常数 */
    double dr;                  /* mean distance between sun and earth
                                            日地平均距离 */
    double omegas;              /* sunrise hour angle
                                            日出时角 */
    double delta;               /* solar declination
                                            太阳磁偏角 */
    int J;                   /* day sequence
                                            日序 */
    double N;                   /* theoretical sunshine duration
                                            理论日照时数 */
    double n;                   /* sunshine duration
                                            日照时数 */
    J   =   day_sequence;

    dr  =   1 + 0.033 * cos(0.01721421 * J);  /* dr = 1 + 0.033*cos(2*PI*J/365) */
    delta   =   0.408 * sin(0.01721421 * J - 1.39);
    omegas  =   acos( -tan(latitude) * tan(delta) );

    N   =   7.63943727 * omegas;    /* N = 24/PI * omegas */
    n   =   sunshine_duration;

    as = 0.25;
    bs = 0.50;

    Ra  =   435.0235 * dr * (omegas * sin(latitude)*sin(delta) + cos(latitude)*cos(delta)*sin(omegas));
                                /* Ra = Gsc / PI *[omegas * sin (latitude)*sin(delta) + cos(latitude)*cos(delta)*sin(omegas)
                                 * Gsc - solar constance 太阳常数, Gsc = 0.0820MJ/m^2/min */
    Rs  =   (as + bs * n / N) * Ra;

    return Rs;
}
