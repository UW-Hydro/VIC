#!/usr/bin/python

#----------------------------------------------------
# Program to convert VIC fluxes files to NetCDF file
# will ask the user wich variable he wants to export
# and also for wich years. Assumes there is data
# for the entire time period, from 1-jan to 31-dec
# SET UP FOR DAILY TIME STEP. FLUX FILE SHOUD NOT
# CONTAIN HOUR RECORD!!
#----------------------------------------------------

#------------------------------------------------
# Writen by Daniel de Castro Victoria
# dvictori@cena.usp.br or daniel.victoria@gmail.com
# Needs python libraries Numeric and Scientific
# 03-dec-2004
#-------------------------------------------------

import os, sys, string
# handle dates...
import datetime
# NetCDF and Numeric
from Scientific.IO.NetCDF import *
from Numeric import *

# checking user input
if len(sys.argv) != 2:
    print "Wrong user input"
    print "Convert VIC fluxes files to NetCDF"
    print "usage flux2cdf.py <vic flux dir>"
    print "VIC FLUX DIR SHOULD CONTAIN TRAILING /"
    sys.exit()

if sys.argv[1][-1] != "/":
    print "VIC FLUX DIR SHOULD CONTAIN TRAILING /"
    print "fixing it for you..."
    sys.argv[1] = sys.argv[1] + "/"

print "IMPORTANT: "+sys.argv[1]+" SHOULD CONTAIN ONLY FLUXES FILES!!!"

# building file list and sorted lat lon list
file_list = os.listdir(sys.argv[1])

lat_t = []
lon_t = []
lat = []
lon = []

for f in file_list:
    lat_t.append(float(string.split(f, "_")[1]))
    lon_t.append(float(string.split(f, "_")[2]))

for i in lat_t:
    if i not in lat:
        lat.append(i)

for i in lon_t:
    if i not in lon:
        lon.append(i)


# putting in order. Lat should be from top to botom
# lon from left to rigth
lon.sort()
lat.sort()
lat.reverse()

del(lat_t)
del(lon_t)

#determining the parameter to use
print "Choose output parameter"
print "1 - Precipitation"
print "2 - Evapotranspiration"
print "3 - Runoff"
print "4 - Base flow"
print "5 - Interception"
print "6 - Soil moisture"
varini = input('Choose output (1 a 6)>')

#getting the collumn right
if varini < 6:
    var = varini + 2
elif varini == 6:        #more than one soil layer...
    camada = input('which soil layer?>')
    var = varini + 1 + camada

#set name of out_file. Named after parameter choice
if var == 3:
    var_txt = "ppt"
    var_name = "Precipitation"
elif var == 4:
    var_txt = "evap"
    var_name = "Evapotranspiration"
elif var == 5:
    var_txt = "runoff"
    var_name = "Runoff"
elif var == 6:
    var_txt = "base"
    var_name = "Baseflow"
elif var == 7:
    var_txt = "intercep"
    var_name = "Interception"
else:
    var_txt = "soil_"+str(camada)
    var_name = "Soil moisture, layer %i", camada

# for what date?
start_year = input("Enter start year:")
end_year = input("End year:")

inidate = datetime.date(start_year,01,01)
enddate = datetime.date(end_year,12,31)

days = enddate.toordinal() - inidate.toordinal()+1

print "Go grab a coffe, this could take a while..."

#
# create array containig all data
# This is going to be huge. Create an array with -9999 (NoData)
# Then populate the array by reading each flux file
#

all_data = zeros([days,len(lat),len(lon)], Float)-9999

c = len(file_list)

# for each file in list
for f in file_list:
    # get lat & lon and it's index
    latitude = float(string.split(f, sep="_")[1])
    longitude = float(string.split(f, sep="_")[2])
    lat_id = lat.index(latitude)
    lon_id = lon.index(longitude)

    print "%i files to write." % c
    c = c -1
    
    infile = open(sys.argv[1]+f, "r")
    lixo = infile.readlines()
    infile.close()
    dado = []

    for l in lixo:
        if int(string.split(l, sep="\t")[0]) in range(inidate.year, enddate.year+1):
            dado.append(float(string.split(l, sep="\t")[var]))
        # putting data inside array.
        # Since data has lat & lon fixed uses dimension [:,lat_index,lon_index]

    all_data[:,lat_id,lon_id] = dado
    


#
# writing NetCDF
#

ncfile = NetCDFFile(var_txt+".nc", "w")

ncfile.Conventions = "COARDS"
ncfile.history = "Created using flux2cdf.py. " + datetime.date.today().isoformat()
ncfile.production = "VIC output"

ncfile.start_date = inidate.isoformat()
ncfile.end_date = enddate.isoformat()

#create dimensions
ncfile.createDimension("X", len(lon))
ncfile.createDimension("Y", len(lat))
ncfile.createDimension("T", days)

#create variables
latvar = ncfile.createVariable("Y", Float, ("Y",))
latvar.long_name = "Latitude"
latvar.units = "degrees_north"
latvar[:] = lat

lonvar = ncfile.createVariable("X", Float, ("X",))
lonvar.long_name = "Longitude"
lonvar.units = "degrees_east"
lonvar[:] = lon

timevar = ncfile.createVariable("T", Float, ("T",))
timevar.long_name = "Time"
timevar.units = "days since " + inidate.isoformat()
timevar[:] = range(0, days)

data_var = ncfile.createVariable(var_txt, Float, ("T","Y","X"))
data_var.long_name = var_name+" calculated by VIC"
data_var.missing_value = -9999.0
data_var.units = "milimeters"
data_var[:] = all_data

ncfile.close()
