#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  3 14:08:51 2021

@author: anjakatzenberger
"""

#--------------------------------------------------------------------###
###               Data Analysis - Monsoon planet                     ###
###                     Energy Budget                                ###
###    Anja Katzeberger, anja.katzenberger@pik-potsdam.de            ###
#--------------------------------------------------------------------###

# calculates monthly sensible heat flux (var: shflx )
# and latent heat release (precipitation over land * L)
# and radiation (sum and difference of individual radiation components)
# and convergence (inferred)

import os
import matplotlib
import matplotlib.pyplot as plt
import iris 
from iris.time import PartialDateTime
import iris.analysis
import numpy as np


lat_low = 10
width = 50
ppm = 280
month_list = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
data_year = 91


#%%---------------------------
### Directories: input and output
#---------------------------
filename = "/p/projects/climber3/anjaka/POEM/work/slab/MonsunPlanet_stripe_10N_output.nc"

filename_waterplanet = "/p/projects/climber3/anjaka/POEM/work/slab/Aquaplanet_output.nc"

save_folder = '/home/anjaka/energy_budget_' + str(ppm) + 'ppm'    
if not os.path.exists(save_folder):
    os.makedirs(save_folder)


#%%---------------------------
### Loading Iris cubes
#---------------------------

### Striped Planet
cubes = iris.load(filename,['precip','shflx','netrad_toa', 'swdn_sfc', 'swup_sfc', 'lwdn_sfc','lwup_sfc'])
cubes_w = iris.load(filename_waterplanet,['precip','shflx','netrad_toa', 'swdn_sfc', 'swup_sfc', 'lwdn_sfc','lwup_sfc'])
print(cubes)
print(cubes_w)

precip = cubes[0]
sh = cubes[1] # sensible heat flux

netrad_toa = cubes[2] # net radiation (sw+lw)
swdn_s = cubes[3]
swup_s = cubes[4]
lwdn_s = cubes[5]
lwup_s = cubes[6]


### Waterplanet
precip_w = cubes_w[0]
sh_w = cubes_w[1] # sensible heat flux

netrad_toa_w = cubes_w[2] # net radiation (sw+lw)
swdn_s_w = cubes_w[3]
swup_s_w = cubes_w[4]
lwdn_s_w = cubes_w[5]
lwup_s_w = cubes_w[6]


# important NOT to convert for latent heat unit
#precip.convert_units('kg m-2 day-1')
#precip.rename('precipitation_rate')
#precip.units = 'mm day-1'
    
#%%-------------------------
### Preprocessing data
#------------------------
# Cutting first years until equilibrium is reached
# Averaging over all remaining 10 Januarys, Februarys,..
# Extracting the region of interest (10-30°N)

cutoff_date = PartialDateTime(year = data_year - 10, month=12, day=31)
land_stripe = iris.Constraint(latitude = lambda v: lat_low <= v <= 30)

### Striped Planet
precip_land = precip.extract(land_stripe)
sh_land = sh.extract(land_stripe)
netrad_toa_land = netrad_toa.extract(land_stripe)
swdn_s_land = swdn_s.extract(land_stripe)
swup_s_land = swup_s.extract(land_stripe)
lwdn_s_land = lwdn_s.extract(land_stripe)
lwup_s_land = lwup_s.extract(land_stripe)

### Water Planet
precip_w_land = precip_w.extract(land_stripe)
sh_w_land = sh_w.extract(land_stripe)
netrad_toa_w_land = netrad_toa_w.extract(land_stripe)
swdn_s_w_land = swdn_s_w.extract(land_stripe)
swup_s_w_land = swup_s_w.extract(land_stripe)
lwdn_s_w_land = lwdn_s_w.extract(land_stripe)
lwup_s_w_land = lwup_s_w.extract(land_stripe)



#%%
#------------------------
###Processing data
#------------------------

### Striped Planet
sh_start = {}
sh_av = {}
sh_av_mean= {}
sh_mean = []
for i in range(1, 13):
    sh_start[i] = sh_land.extract(iris.Constraint(time=lambda cell: cell.point.month == i and cutoff_date < cell.point))
    sh_av[i] = sh_start[i].collapsed('time', iris.analysis.MEAN)
    sh_av[i].coord('latitude').guess_bounds()
    sh_av[i].coord('longitude').guess_bounds()
    aws = iris.analysis.cartography.area_weights(sh_av[i])
    sh_av_mean[i] = sh_av[i].collapsed(['latitude','longitude'], iris.analysis.MEAN, weights=aws)
    sh_mean.append(sh_av_mean[i].data)

precip_start = {}
precip_av = {}
precip_av_mean = {}
precip_mean = []
for i in range(1, 13):
    precip_start[i] = precip_land.extract(iris.Constraint(time=lambda cell: cell.point.month == i and cutoff_date < cell.point))
    precip_av[i] = precip_start[i].collapsed('time', iris.analysis.MEAN)
    precip_av[i].coord('latitude').guess_bounds()
    precip_av[i].coord('longitude').guess_bounds()    
    aws = iris.analysis.cartography.area_weights(precip_av[i])
    precip_av_mean[i] = precip_av[i].collapsed(['latitude','longitude'], iris.analysis.MEAN, weights=aws)
    precip_mean.append(precip_av_mean[i].data)

netrad_toa_start = {}
netrad_toa_av = {}
netrad_toa_av_mean = {}
netrad_toa_mean = []
for i in range(1, 13):
    netrad_toa_start[i] = netrad_toa_land.extract(iris.Constraint(time=lambda cell: cell.point.month == i and cutoff_date < cell.point))
    netrad_toa_av[i] = netrad_toa_start[i].collapsed('time', iris.analysis.MEAN)
    netrad_toa_av[i].coord('latitude').guess_bounds()
    netrad_toa_av[i].coord('longitude').guess_bounds()    
    aws = iris.analysis.cartography.area_weights(netrad_toa_av[i])
    netrad_toa_av_mean[i] = netrad_toa_av[i].collapsed(['latitude','longitude'], iris.analysis.MEAN, weights=aws)
    netrad_toa_mean.append(netrad_toa_av_mean[i].data)

swdn_s_start = {}
swdn_s_av = {}
swdn_s_av_mean = {}
swdn_s_mean = []
for i in range(1, 13):
    swdn_s_start[i] = swdn_s_land.extract(iris.Constraint(time=lambda cell: cell.point.month == i and cutoff_date < cell.point))
    swdn_s_av[i] = swdn_s_start[i].collapsed('time', iris.analysis.MEAN)
    swdn_s_av[i].coord('latitude').guess_bounds()
    swdn_s_av[i].coord('longitude').guess_bounds()    
    aws = iris.analysis.cartography.area_weights(swdn_s_av[i])
    swdn_s_av_mean[i] = swdn_s_av[i].collapsed(['latitude','longitude'], iris.analysis.MEAN, weights=aws)
    swdn_s_mean.append(swdn_s_av_mean[i].data)

swup_s_start = {}
swup_s_av = {}
swup_s_av_mean = {}
swup_s_mean = []
for i in range(1, 13):
    swup_s_start[i] = swup_s_land.extract(iris.Constraint(time=lambda cell: cell.point.month == i and cutoff_date < cell.point))
    swup_s_av[i] = swup_s_start[i].collapsed('time', iris.analysis.MEAN)
    swup_s_av[i].coord('latitude').guess_bounds()
    swup_s_av[i].coord('longitude').guess_bounds()    
    aws = iris.analysis.cartography.area_weights(swup_s_av[i])
    swup_s_av_mean[i] = swup_s_av[i].collapsed(['latitude','longitude'], iris.analysis.MEAN, weights=aws)
    swup_s_mean.append(swup_s_av_mean[i].data)
    
lwdn_s_start = {}
lwdn_s_av = {}
lwdn_s_av_mean = {}
lwdn_s_mean = []
for i in range(1, 13):
    lwdn_s_start[i] = lwdn_s_land.extract(iris.Constraint(time=lambda cell: cell.point.month == i and cutoff_date < cell.point))
    lwdn_s_av[i] = lwdn_s_start[i].collapsed('time', iris.analysis.MEAN)
    lwdn_s_av[i].coord('latitude').guess_bounds()
    lwdn_s_av[i].coord('longitude').guess_bounds()    
    aws = iris.analysis.cartography.area_weights(lwdn_s_av[i])
    lwdn_s_av_mean[i] = lwdn_s_av[i].collapsed(['latitude','longitude'], iris.analysis.MEAN, weights=aws)
    lwdn_s_mean.append(lwdn_s_av_mean[i].data)

lwup_s_start = {}
lwup_s_av = {}
lwup_s_av_mean = {}
lwup_s_mean = []
for i in range(1, 13):
    lwup_s_start[i] = lwup_s_land.extract(iris.Constraint(time=lambda cell: cell.point.month == i and cutoff_date < cell.point))
    lwup_s_av[i] = lwup_s_start[i].collapsed('time', iris.analysis.MEAN)
    lwup_s_av[i].coord('latitude').guess_bounds()
    lwup_s_av[i].coord('longitude').guess_bounds()    
    aws = iris.analysis.cartography.area_weights(lwup_s_av[i])
    lwup_s_av_mean[i] = lwup_s_av[i].collapsed(['latitude','longitude'], iris.analysis.MEAN, weights=aws)
    lwup_s_mean.append(lwup_s_av_mean[i].data)


### Water Planet
cutoff_date = PartialDateTime(year = 41 - 10, month=12, day=31)

sh_w_start = {}
sh_w_av = {}
sh_w_av_mean= {}
sh_w_mean = []
for i in range(1, 13):
    sh_w_start[i] = sh_w_land.extract(iris.Constraint(time=lambda cell: cell.point.month == i and cutoff_date < cell.point))
    sh_w_av[i] = sh_w_start[i].collapsed('time', iris.analysis.MEAN)
    sh_w_av[i].coord('latitude').guess_bounds()
    sh_w_av[i].coord('longitude').guess_bounds()
    aws = iris.analysis.cartography.area_weights(sh_w_av[i])
    sh_w_av_mean[i] = sh_w_av[i].collapsed(['latitude','longitude'], iris.analysis.MEAN, weights=aws)
    sh_w_mean.append(sh_w_av_mean[i].data)

precip_w_start = {}
precip_w_av = {}
precip_w_av_mean= {}
precip_w_mean = []
for i in range(1, 13):
    precip_w_start[i] = precip_w_land.extract(iris.Constraint(time=lambda cell: cell.point.month == i and cutoff_date < cell.point))
    precip_w_av[i] = precip_w_start[i].collapsed('time', iris.analysis.MEAN)
    precip_w_av[i].coord('latitude').guess_bounds()
    precip_w_av[i].coord('longitude').guess_bounds()
    aws = iris.analysis.cartography.area_weights(precip_w_av[i])
    precip_w_av_mean[i] = precip_w_av[i].collapsed(['latitude','longitude'], iris.analysis.MEAN, weights=aws)
    precip_w_mean.append(precip_w_av_mean[i].data)

netrad_toa_w_start = {}
netrad_toa_w_av = {}
netrad_toa_w_av_mean= {}
netrad_toa_w_mean = []
for i in range(1, 13):
    netrad_toa_w_start[i] = netrad_toa_w_land.extract(iris.Constraint(time=lambda cell: cell.point.month == i and cutoff_date < cell.point))
    netrad_toa_w_av[i] = netrad_toa_w_start[i].collapsed('time', iris.analysis.MEAN)
    netrad_toa_w_av[i].coord('latitude').guess_bounds()
    netrad_toa_w_av[i].coord('longitude').guess_bounds()
    aws = iris.analysis.cartography.area_weights(netrad_toa_w_av[i])
    netrad_toa_w_av_mean[i] = netrad_toa_w_av[i].collapsed(['latitude','longitude'], iris.analysis.MEAN, weights=aws)
    netrad_toa_w_mean.append(netrad_toa_w_av_mean[i].data)


swdn_s_w_start = {}
swdn_s_w_av = {}
swdn_s_w_av_mean= {}
swdn_s_w_mean = []
for i in range(1, 13):
    swdn_s_w_start[i] = swdn_s_w_land.extract(iris.Constraint(time=lambda cell: cell.point.month == i and cutoff_date < cell.point))
    swdn_s_w_av[i] = swdn_s_w_start[i].collapsed('time', iris.analysis.MEAN)
    swdn_s_w_av[i].coord('latitude').guess_bounds()
    swdn_s_w_av[i].coord('longitude').guess_bounds()
    aws = iris.analysis.cartography.area_weights(swdn_s_w_av[i])
    swdn_s_w_av_mean[i] = swdn_s_w_av[i].collapsed(['latitude','longitude'], iris.analysis.MEAN, weights=aws)
    swdn_s_w_mean.append(swdn_s_w_av_mean[i].data)

swup_s_w_start = {}
swup_s_w_av = {}
swup_s_w_av_mean= {}
swup_s_w_mean = []
for i in range(1, 13):
    swup_s_w_start[i] = swup_s_w_land.extract(iris.Constraint(time=lambda cell: cell.point.month == i and cutoff_date < cell.point))
    swup_s_w_av[i] = swup_s_w_start[i].collapsed('time', iris.analysis.MEAN)
    swup_s_w_av[i].coord('latitude').guess_bounds()
    swup_s_w_av[i].coord('longitude').guess_bounds()
    aws = iris.analysis.cartography.area_weights(swup_s_w_av[i])
    swup_s_w_av_mean[i] = swup_s_w_av[i].collapsed(['latitude','longitude'], iris.analysis.MEAN, weights=aws)
    swup_s_w_mean.append(swup_s_w_av_mean[i].data)

lwdn_s_w_start = {}
lwdn_s_w_av = {}
lwdn_s_w_av_mean= {}
lwdn_s_w_mean = []
for i in range(1, 13):
    lwdn_s_w_start[i] = lwdn_s_w_land.extract(iris.Constraint(time=lambda cell: cell.point.month == i and cutoff_date < cell.point))
    lwdn_s_w_av[i] = lwdn_s_w_start[i].collapsed('time', iris.analysis.MEAN)
    lwdn_s_w_av[i].coord('latitude').guess_bounds()
    lwdn_s_w_av[i].coord('longitude').guess_bounds()
    aws = iris.analysis.cartography.area_weights(lwdn_s_w_av[i])
    lwdn_s_w_av_mean[i] = lwdn_s_w_av[i].collapsed(['latitude','longitude'], iris.analysis.MEAN, weights=aws)
    lwdn_s_w_mean.append(lwdn_s_w_av_mean[i].data)
    
lwup_s_w_start = {}
lwup_s_w_av = {}
lwup_s_w_av_mean= {}
lwup_s_w_mean = []
for i in range(1, 13):
    lwup_s_w_start[i] = lwup_s_w_land.extract(iris.Constraint(time=lambda cell: cell.point.month == i and cutoff_date < cell.point))
    lwup_s_w_av[i] = lwup_s_w_start[i].collapsed('time', iris.analysis.MEAN)
    lwup_s_w_av[i].coord('latitude').guess_bounds()
    lwup_s_w_av[i].coord('longitude').guess_bounds()
    aws = iris.analysis.cartography.area_weights(lwup_s_w_av[i])
    lwup_s_w_av_mean[i] = lwup_s_w_av[i].collapsed(['latitude','longitude'], iris.analysis.MEAN, weights=aws)
    lwup_s_w_mean.append(lwup_s_w_av_mean[i].data)

#%%
### Calculations 
    
# 1.) calculate latent heat release by condensation
    
L = 2.6 * 1000000 #latent heat of condensation (J/kg) 
# Monsoon Planet
latent_heat = [x * L for x in precip_mean]
# Aquaplanet
latent_heat_w = [x * L for x in precip_w_mean]



# 2.) Radiation calculation

# MonsoonPlanet
netrad_s = []
radiation = []
for i in range(0,12):
    print(i)
    netrad_s.append(swdn_s_mean[i] - swup_s_mean[i] + lwdn_s_mean[i] - lwup_s_mean[i])
    radiation.append(netrad_toa_mean[i] - netrad_s[i])

#Aquaplanet
netrad_s_w = []
radiation_w = []

for i in range(0,12):
    netrad_s_w.append(swdn_s_w_mean[i] - swup_s_w_mean[i] + lwdn_s_w_mean[i] - lwup_s_w_mean[i])
    radiation_w.append(netrad_toa_w_mean[i] - netrad_s_w[i])


# 3.) Convergence
    
convergence = []
for i in range(0,12):
    convergence.append(-(sh_mean[i]+radiation[i]+latent_heat[i]))

convergence_w = []
for i in range(0,12):
    convergence_w.append(-(sh_w_mean[i]+radiation_w[i]+latent_heat_w[i]))


# 4.) Sensible heat flux
# directly available in dataset



#--------------------------
# Displaying 2D map
#--------------------------

plt.figure(figsize=(10,5))
plt.xlim(1,12)
plt.ylim([-180,250])
plt.xlabel('Month')
plt.ylabel('Energy flux [$W/m^2$]')

plt.plot(np.arange(1,13,1), sh_mean, color = "red", label = "sensible") #sensible
plt.plot(np.arange(1,13,1),latent_heat, color = "blue", label = "latent") # latent
plt.plot(np.arange(1,13,1), radiation, color = "black", label = "radiative") #radiative
plt.plot(np.arange(1,13,1),convergence, color = "darkorange", label = "covergence")

plt.plot(np.arange(1,13,1), sh_w_mean, color = "red",  linestyle = "dashed") #sensible
plt.plot(np.arange(1,13,1),latent_heat_w, color = "blue",  linestyle = "dashed") # latent
plt.plot(np.arange(1,13,1), radiation_w, color = "black",  linestyle = "dashed") #radiative
plt.plot(np.arange(1,13,1),convergence_w, color = "darkorange", linestyle = "dashed")

plt.subplots_adjust(right=0.65)
plt.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0, frameon = False)
plt.axhline(0, color = "black", linewidth = 0.8)
 
dirsave = save_folder + '/budget_Monsoonplanet_plus_waterplanet.pdf'    
plt.savefig(dirsave, bbox_inches='tight')
plt.show()
plt.close()
    
    
