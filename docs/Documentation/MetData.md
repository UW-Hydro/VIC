# Raw Daily Meteorological Data for the U.S.

This page has links to pages with state-by-state raw data for:

*   Precipitation
*   Maximum Daily Temperature
*   Minimum Daily Temperature

The raw data is in a format to be used by the various preprocessing scripts (e.g., _preproc_<parameter>.scr_) for creating meteorological forcing files for the VIC model. The first set of data that needs to be processed is the daily data through 1997\. The scripts for processing this raw data create a list of stations that have been filtered according to the following criteria:

*   Minimum of 10 years of record
*   Minimum of 50% coverage (i.e., >50% of days in period of record have valid data)
*   Station has a recorded elevation

After using this data, daily data for more recent years can be appended to the data, also using programs described in the section on file preparation for VIC. The more recent data contains one file for each state which contains Precipitation, Tmax, and Tmin. This more recent data contains data for all active stations in each state, and does not include any filtering. This is the reason for relying on the station lists prepared for the pre-1998 data for the later data also.

The raw data can be downloaded on a state-by-state basis here:

For Station Inception through 12/31/1997 (data and station information files for each state, which must be used together):

*   [Daily Precipitation Data](http://www.hydro.washington.edu/Lettenmaier/Data/prcp_thru_97.html)
*   [Daily Maximum Temperature Data](http://www.hydro.washington.edu/Lettenmaier/Data/tmax_thru_97.html)
*   [Minimum Temperature Data](http://www.hydro.washington.edu/Lettenmaier/Data/tmin_thru_97.html)

Daily data for more recent years::

*   [Daily Data from 1998](http://www.hydro.washington.edu/Lettenmaier/Data/data_1998.html)
*   [Daily Data from 1999](http://www.hydro.washington.edu/Lettenmaier/Data/data_1999.html)
*   [Daily Data through July, 2000](http://www.hydro.washington.edu/Lettenmaier/Data/data_2000.html)
*   [Daily Data, Aug-Dec, 2000](http://www.hydro.washington.edu/Lettenmaier/Data/data_2000b.html)
