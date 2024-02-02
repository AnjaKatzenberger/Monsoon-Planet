#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#-------------------------------------------------------------
### GENERAL CHARACTERISTICS OF THE STANDARD RUN - part2 ###
# Anja Katzeberger, anja.katzenberger@pik-potsdam.de
#-------------------------------------------------------------

# calculates monthly mean of precipitation spatial distribution (10 years average) and surface temperature and precip meridional distribution
# Creates figure with spatial rainfall distribution (mm/day) and zonal distribution of precip surface temperate (°C) for the 12 months

# The directories of the input files as well as directories for saving the figures must be adapted. 

import netCDF4 as nc
import statistics as stat
import matplotlib.pyplot as plt
import numpy as np
import iris.plot as iplt
import pylab as pl
from scipy.signal import argrelextrema
import os
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
import iris 
import cartopy.crs as ccrs
import iris.plot as iplt
import iris.quickplot as qplt
from iris.time import PartialDateTime
import iris.analysis
import numpy as np


lat_low = 10
width = 50
var = "precip_tsurf"
ppm = 280
startyear = 91

month_list = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']

filename = "/p/projects/climber3/anjaka/POEM/work/slab/B_data/MonsoonPlanet_slab_50_output.nc"
save_folder = '/home/anjaka/spatial_distribution_' + var + '_' + str(lat_low) + "N_" +  str(width) + '_' + str(ppm)+'ppm_months'    

#filename = "/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/A_data/MonsoonPlanet_slab_50.nc"
#save_folder = '/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Neu'   

if not os.path.exists(save_folder):
    os.makedirs(save_folder)


#%%--------------------------
### Loading Iris cubes
#---------------------------
cubes = iris.load(filename,['precip','t_surf'])
 
precip = cubes[0]
tsurf = cubes[1] 

# converting units from kg m-2 day-1 to mm day-1 (pr) and from Kelvin to Celsius (tsurf)
precip.convert_units('kg m-2 day-1')
precip.rename('precipitation_rate')
precip.units = 'mm day-1'

tsurf.convert_units('Celsius')
    
#%%-------------------------
### Averaging data
#------------------------
# Cutting first years until equilibrium is reached
# Averaging over all remaining 10 Januarys, Februarys,..

cutoff_date = PartialDateTime(year= startyear - 10, month=12, day=31)

precip_start = {}
precip_av = {}
precip_sum_lat = {}
precip_sum = {}
precip_av_lon = {}
precip_av_lat = {}
for i in range(1, 13):
    precip_start[i] = precip.extract(iris.Constraint(time=lambda cell: cell.point.month == i and cutoff_date < cell.point))
    precip_av[i] = precip_start[i].collapsed('time', iris.analysis.MEAN)
    precip_av_lat[i] = precip_av[i].collapsed('longitude', iris.analysis.MEAN) 

tsurf_start = {}
tsurf_av = {}
tsurf_av_lat = {}
for i in range(1, 13):
    tsurf_start[i] = tsurf.extract(iris.Constraint(time=lambda cell: cell.point.month == i and cutoff_date < cell.point))
    tsurf_av[i] = tsurf_start[i].collapsed('time', iris.analysis.MEAN)
    tsurf_av_lat[i] = tsurf_av[i].collapsed('longitude', iris.analysis.MEAN)



#%%--------------------------
# Displaying 2D map
#--------------------------
# Load a Cynthia Brewer palette.
brewer_cmap = mpl_cm.get_cmap("brewer_Blues_09")

## left panel: precip distribution
for m in range(1,13):
    fig = plt.figure(figsize=(10,5))
    ax1 = fig.add_subplot(121)
    pos1 = ax1.get_position()
    ax1.set_position([0.03,pos1.y0,pos1.width,pos1.height])
    ax1.set_ylim(-90,90)
    ax1 = qplt.contourf(precip_av[m], np.arange(0,10,1), cmap=brewer_cmap, extend = "max")
    ax1 = plt.gca().set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree())
    ax1 = plt.gca().set_xticks([-120,-60,0,60,120], crs=ccrs.PlateCarree())
    plt.axhline(lat_low, linewidth = 0.7, color = 'black')
    plt.axhline(lat_low+width, linewidth = 0.7, color = 'black')
    plt.axhline(0, linewidth = 0.7, color = 'black',linestyle = '--')
    ax1 = plt.title('')
         
### right panel: meridional distribution    
    ax2 = fig.add_subplot(122)
    ax3 = ax2.twiny()
    pos2 = ax2.get_position()
    ax2.set_position([pos1.x0+pos1.width+0.005,0.34,pos1.width/4,0.459*pos1.height])
    ax2.set_xlim([0,30])
    ax2.set_xticks((0,15,30))
    ax2.tick_params(axis='x', colors='#015482')
    ax2.set_xlabel('Precipitation [mm/day]', color = '#015482')
    ax2 = iplt.plot(precip_av_lat[m],precip.coord('latitude'), color = '#015482',axes = ax2)

    ax3 = plt.xlim([-45,55])
    ax3 = plt.xticks((-30,0,30))
    ax3 = plt.xlabel('Surface temperature [°C]', color = 'orange')
    ax3 = iplt.plot(tsurf_av_lat[m],tsurf.coord('latitude'), color = 'orange')

        
    # creates horizonal line at tsurf max in tsurf zonal distribution
   # tsurf_max = tsurf_av_lat[m].collapsed('latitude', iris.analysis.MAX)
   # plt.axhline(tsurf_av_lat[m].coord('latitude').points[tsurf_av_lat[m].data == tsurf_max.data],linestyle = 'dashed', color = "grey")
    
    plt.axhline(lat_low, linewidth = 0.7, color = 'black')
    plt.axhline(lat_low + width, linewidth = 0.7, color = 'black')
    plt.axhline(0, linewidth = 0.7, color = 'black',linestyle = '--')
    plt.yticks(())
    plt.ylim([-90,90])
                                                                                                                                  
    dirsave = save_folder + '/spatial_distribution_' + str(m) + '.pdf'    
    plt.savefig(dirsave, bbox_inches='tight')
