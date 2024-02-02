#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#-------------------------------------------------------------
### GENERAL CHARACTERISTICS OF THE STANDARD RUN - part 1 ###
# Anja Katzeberger, anja.katzenberger@pik-potsdam.de
#-------------------------------------------------------------

# This script creates the figures as displayed in the subchapter "General characteristics of the standard run (10-60°N)"
# The directories of the input files as well as directories for saving the figures must be adapted. 

#%%
### SETTING PATH TO ENVIRONMENT
#------------------------
import sys
sys.path.insert(0, r"C:\Users\anjaka\Nextcloud\PhD\01_Monsoon_Planet\Codes\Publication\venv_mp\lib\site-packages")

#%%
### INSTALLING PACKAGES
#------------------------

import netCDF4 as nc
import statistics as stat
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import argrelextrema


#%%
### DIRECTORIES
#------------------------

fn_50 = 'C:/Users/anjaka/Nextcloud/PhD/01_Monsoon_Planet/A_data//MonsoonPlanet_slab_50.nc'
#mas_50 ='/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/data/MonsoonPlanet_10N_mastrfu.nc'

savedir = 'C:/Users/anjaka/Nextcloud/PhD/01_Monsoon_Planet/Figures/00_Publication_figures'

#%%
### LODING DATA
#------------------------

month_list = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
clrs = ['#6EA6CD','#98CAE1','#C2E4EF','#ACD39E','#5AAE61','#FDB366','#F67E4B','#DD3D2D','#A50026','#762A83','#36489A','#4A7BB7']

ds_50 = nc.Dataset(fn_50)
#ds_mas_50 = nc.Dataset(mas_50)

lats = ds_50.variables['lat'][:]
lons = ds_50.variables['lon'][:]
pfull = ds_50.variables['pfull'][:]

evap_50 = ds_50['evap'][:]
precip_50 = ds_50['precip'][:]
tsurf_50 = ds_50['t_surf'][:]
ps_50 = ds_50['ps'][:]
olr_50 = ds_50['olr'][:]
vcomp_50 = ds_50['vcomp'][:]
ucomp_50 = ds_50['ucomp'][:]
#mastrfu_50 = np.squeeze(ds_mas_50['mastrfu'][:])
sw_50 = ds_50['swdn_toa'][:]

evap_z_50 = evap_50.mean(axis=2)*86400
precip_z_50 = precip_50.mean(axis = 2)*86400
olr_z_50 = olr_50.mean(axis=2)
vcomp_z_50 = vcomp_50.mean(axis=3)
ucomp_z_50 = ucomp_50.mean(axis=3)
ps_z_50 = ps_50.mean(axis=2)
tsurf_z_50 = tsurf_50.mean(axis=2)
sw_z_50 = sw_50.mean(axis=2)


#%%
#------------------------
### General characteristics of the 10N configuration
#------------------------

# Precipitation (contour)
plt.figure()
plt.title('Precipitation')
plt.contourf(range(1,13),lats,precip_z_50.transpose(),cmap="Blues",extend="max",levels = [0,2,4,6,8,10,12,14,16,18,20])
plt.axhline(0, linewidth = 0.7, color = 'black', linestyle = '--')
plt.axhline(10, linewidth = 0.7, color = 'black')
plt.axhline(60, linewidth = 0.7, color = 'black')
plt.xlabel('Months')
plt.xlim([1,12])
plt.ylabel('Latitude')
plt.ylim([-90,90])
clb = plt.colorbar()
clb.ax.set_title('mm/day',fontsize=8)
#plt.savefig(savedir + '/precip_10N_280ppm.pdf')

# Evaporation (contour)
plt.figure()
plt.title('Evaporation')
plt.contourf(range(1,13),lats,evap_z_50.transpose(),cmap="Blues",extend="max",levels = [0,1,2,3,4,5,6,7,8,9,10])
plt.axhline(0, linewidth = 0.7, color = 'black',linestyle = '--')
plt.axhline(10, linewidth = 0.7, color = 'black')
plt.axhline(60, linewidth = 0.7, color = 'black')
plt.xlabel('Months')
plt.xlim([1,12])
plt.ylabel('Latitude')
plt.ylim([-90,90])
clb = plt.colorbar()
clb.ax.set_title('mm/day',fontsize=8)
#plt.savefig(savedir + '/evap_10N_280ppm.pdf')

# Meridional P-E distribution
fig = plt.figure()
for m in range(0,12):
    plt.xlabel('latitude')
    plt.xlim([-91,91])
    plt.ylim([-15,30])
    plt.xticks(np.arange(-90, 91, step=30))
    plt.ylabel( 'precipitation - evaporation (mm/day)')
    plt.axvspan(10, 60, facecolor='grey', alpha=0.0022)
    plt.axvline(0, linewidth = 0.7, color = 'black', linestyle = '--')
    plt.axvline(10, linewidth = 0.7, color = 'black')
    plt.axvline(60, linewidth = 0.7, color = 'black')
    plt.axhline(0, linewidth = 0.7, color = 'black')
    plt.plot(lats,precip_z_50[m]-evap_z_50[m],label=month_list[m], color = clrs[m], linewidth = 1) #alpha = 1 
plt.legend(loc= 'center left', frameon = False, bbox_to_anchor = (1,0.5))
#plt.savefig(savedir + '/4_Sensitivity_slab/fig_A1_precip_minus_evap_10N_280ppm.pdf', bbox_to_anchor = 'tight' )

# Surface temperature (contour)
plt.figure()
plt.title('Surface Temperature')
plt.contourf(range(1,13),lats,tsurf_z_50.transpose()-273.15,cmap="RdBu_r",extend="both",levels = range(-35,35))
plt.axhline(0, linewidth = 0.7, color = 'black',linestyle = '--')
plt.axhline(10, linewidth = 0.7, color = 'black')
plt.axhline(60, linewidth = 0.7, color = 'black')
plt.xlabel('Months')
plt.xlim([1,12])
plt.ylabel('Latitude')
plt.ylim([-90,90])
clb = plt.colorbar()
clb.ax.set_title('°C',fontsize=8)
#plt.savefig(savedir + '/tsurf_10N_280ppm.pdf')

# Surface pressure (contour)
plt.figure()
plt.title('Surface Pressure')
plt.contourf(range(1,13),lats,ps_z_50.transpose()/100,cmap="Greens",extend="both",levels = range(960,1020,5))
plt.axhline(0, linewidth = 0.7, color = 'black',linestyle = '--')
plt.axhline(10, linewidth = 0.7, color = 'black')
plt.axhline(60, linewidth = 0.7, color = 'black')
plt.xlabel('Months')
plt.xlim([1,12])
plt.ylabel('Latitude')
plt.ylim([-90,90])
clb = plt.colorbar()
clb.ax.set_title('hPa',fontsize=8)
#plt.savefig(savedir + '/ps_10N_280ppm.pdf')

# Meridional winds at 5°N
plt.figure()
vcomp_ad = np.squeeze(vcomp_z_50[:,:,47]) # 47 refers to 5.1°N, check lats[47]
plt.title('Meridional wind')
ax = plt.contourf(range(1,13),pfull,vcomp_ad.transpose(),extend="both",cmap='RdBu',levels = [-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5])#drange(0, 100, '0.1')
ax = plt.gca()
ax.set_ylim(ax.get_ylim()[::-1])
clb = plt.colorbar()
clb.ax.set_title('m/s',fontsize=8)
plt.xlabel('Months')
plt.xlim([1,12])
plt.ylabel('Altitude [hPa]')
plt.ylim([1000,0])
#plt.savefig(savedir + '/vcomp_11S.pdf')


# Incoming shortwave radiation at the TOA for selected months
plt.figure()
plt.plot(lats,sw_z_50[12-1], label = "December")
plt.plot(lats,sw_z_50[6-1], label = "June")
plt.plot(lats,sw_z_50[3-1], label = "March/ September")
#plt.axvline(0,linewidth = 0.7,color="black",linestyle = '--')
plt.axhline(0,linewidth = 0.7,color="black")
plt.xlim([-90,90])
plt.ylim([0,600])
plt.xlabel("Latitudes")
plt.ylabel("Incoming SW radiation - TOA [$W/m^2$]")
plt.legend(frameon = False)
#plt.savefig(savedir + '/2_General_characteristics_10N_280ppm/sw_10N_280ppm.pdf')




#%% 
### MASS STREAM FUNCTION
#------------------------

# Mass stream function
# for i in range(0,12):
#     plt.figure()
#     #plt.contourf(lats,pfull,np.squeeze(mastrfu_50[i]),cmap='PuOr',levels = [-4e9,-3.5e9,-3e9,-2.5e9,-2e9,-1.5e9,-1e9,-0.5e9,0,0.5e9,1e9,1.5e9,2e9,2.5e9,3e9,3.5e9,4e9],extend = 'both')
#     plt.contour(lats,pfull,np.squeeze(mastrfu_50[i]),levels = [-4e9,-3.5e9,-3e9,-2.5e9,-2e9,-1.5e9,-1e9,-0.5e9,0,0.5e9,1e9,1.5e9,2e9,2.5e9,3e9,3.5e9,4e9],linewidths = 1,colors = "darkgrey",extend = 'both')
#    # plt.contour(lats,pfull,np.squeeze(mastrfu_50[i]),levels = [-4e9,-3.5e9,-3e9,-2.5e9,-2e9,-1.5e9,-1e9,-0.5e9,0,0.5e9,1e9,1.5e9,2e9,2.5e9,3e9,3.5e9,4e9],linewidths = 1,cmap = "PuOr",extend = 'both')
    
#     plt.colorbar(extend = 'both')
#     plt.axvline(0, linewidth = 0.7, color = 'black',linestyle = '--')
#     plt.axvline(10, linewidth = 0.7, color = 'black')
#     plt.axvline(60, linewidth = 0.7, color = 'black')
#     plt.axhline(500, linewidth = 0.7, color = 'black')
#     plt.ylabel('Height [hPa]')
#     plt.xlabel('Latitude')
#     plt.xlim([-91,91])
#     plt.ylim([1000,0])
#     plt.title(month_list[i])
#     plt.tight_layout()
#     #plt.savefig(savedir + '/Slab_depth/streamfkt_50_' + str(i) +'.pdf' )
#     plt.close()

