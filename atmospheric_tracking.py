#!/usr/bin/env python3

from scipy.io.netcdf import *
from netCDF4 import Dataset
from numpy import *
import matplotlib.pyplot as plt
import sys

import netCDF4 as nc
from datetime import datetime
today = datetime.today()

# Provide the month number on the command line, otherwise take July
try:
    month=int(sys.argv[1])
except:
    month=07

fn=Dataset('data_n/Monthly_Avg_EvapLoss.nc') # Evaporation from reservoirs
lonas=fn.variables['lon'][:]
latas=fn.variables['lat'][:]
evaps=fn.variables['totEvM'+(str(month).lstrip('0'))][:]

print(evaps.shape,sum(evaps))

f=Dataset('data_in/utrack_climatology_0.5_'+str(month).zfill(2)+'.nc')

lats=arange(90,-90,-0.5)
lons=arange(0,360,0.5)

forward_footprint_tot = zeros(shape=(len(lats), len(lons)))

def get_closest_index(lats,lat):
        import operator
        lat_index, min_value = min(enumerate(abs(lats-lat)), key=operator.itemgetter(1))
        return lat_index

def get_footprints(latitude,longitude,evap):
    # Determine the closest indices for the location
    latidx=get_closest_index(lats,latitude)
    lonidx=get_closest_index(lons,longitude)

    # Determine the forward footprint, where the ET from the input location will subsequently rain out.
    fp=f.variables['moisture_flow'][latidx,lonidx]
    fp=fp*-0.1
    fp=e**fp
    fp[fp==1]=0
    forward_fp=fp*evap

    return forward_fp

for i,lona in enumerate(lonas):
    for j,lata in enumerate(latas):
        evap=evaps[j,i]
        if evap > 0:
            forward_footprint =get_footprints(lata,lona,evap)
            forward_footprint[forward_footprint.mask]=0.
            forward_footprint.mask=False
            forward_footprint=forward_footprint/sum(forward_footprint)*evap
            forward_footprint_tot=forward_footprint_tot+forward_footprint

lons2=arange(-180,180,0.5) #shift longitude
forward_footprint2=forward_footprint_tot.copy()
forward_footprint2[:,0:360]=forward_footprint_tot[:,360:720]
forward_footprint2[:,360:720]=forward_footprint_tot[:,0:360]

f = nc.Dataset('output.nc','w', format='NETCDF4')
f.createDimension('longitude', len(lons2))
f.createDimension('latitude', len(lats))
longitude = f.createVariable('longitude', 'f4', 'longitude')
latitude = f.createVariable('latitude', 'f4', 'latitude')
ff = f.createVariable('forward_footprint', 'f4', ('latitude', 'longitude'), fill_value=-999.)
longitude[:] = lons2
latitude[:] = lats
ff[:,:] = forward_footprint2

#Add global attributes
f.description = "evapotraspiration footprint"
f.history = "created " + today.strftime("%d/%m/%y")
f.Conventions = 'CF-1.6'

# add attributes to dimension varables
longitude.units = 'degrees_east'
longitude.long_name = 'longitude'

latitude.units = 'degrees_north'
latitude.long_name = 'latitude'

f.close()

