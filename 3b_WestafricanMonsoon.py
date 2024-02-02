#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 18:53:34 2023

@author: anjakatzenberger
"""

#-------------------------------------------------------------
###         REALITY CHECK: WESTAFRICAN MONSOON             ###
# Anja Katzeberger, anja.katzenberger@pik-potsdam.de
#-------------------------------------------------------------

# This script creates the figures as created for the comparison with the Westafrican monsoon and some more figures used in the context of the analysis
# The directories of the input files as well as directories for saving the figures must be adapted. 


import netCDF4 as nc
import statistics as stat
import matplotlib.pyplot as plt
import numpy as np
import iris.plot as iplt
import pylab as pl
from scipy.signal import argrelextrema
from matplotlib.legend_handler import HandlerTuple
import matplotlib

month_list = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
clrs = ['#6EA6CD','#98CAE1','#C2E4EF','#ACD39E','#5AAE61','#FDB366','#F67E4B','#DD3D2D','#A50026','#762A83','#36489A','#4A7BB7']


### Loading data
data_pr = '/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/A_data/gswp3-w5e5_obsclim_pr_global_daily_1981_2010_ymonmean.nc'
data_ps ='/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/A_data/gswp3-w5e5_obsclim_ps_global_daily_1981_2010_ymonmean.nc'

ds_pr = nc.Dataset(data_pr)
ds_ps = nc.Dataset(data_ps)

lats = ds_pr.variables['lat'][:]
lons = ds_pr.variables['lon'][:]

pr = ds_pr['pr'][:]
ps = ds_ps['ps'][:]

pr_z = pr.mean(axis = 2)*86400
ps_z = ps.mean(axis=2)/100



#%%

# REGION: West africa region from 5.25°E to 15.25°W

# Meridional rainfall distribution 
pr_a = pr[:,:,370:390] 
pr_z = pr_a.mean(axis = 2)*86400 # to mm/day

ps_a = ps[:,:,370:390]
ps_z = ps_a.mean(axis = 2)/100 # to hPa

# Precipitation
fig = plt.figure()
for m in range(5,9):
    plt.xlabel('Latitude')
    plt.ylim([0,12])
    plt.xticks(np.arange(0, 17, step=2.5))
    plt.ylabel( 'Precipitation (mm/day)')
    plt.axvline(0, linewidth = 0.7, color = 'black',linestyle = '--')
   # plt.axvspan(5,7.3,alpha = 0.02)
   # plt.axvspan(10,12.3,alpha = 0.02)
    exec("plt.plot(lats,pr_z[m],label=month_list[m], color = clrs[m], linewidth = 1)") #alpha = 1 
    plt.xlim([0,15])
plt.legend(loc= 'center left', bbox_to_anchor = (1,0.5), frameon = False)
plt.tight_layout()
plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/00_Publication_figures/Reality_check/westafrica_pr.pdf", bbox_inches='tight')

# Surface pressure
fig = plt.figure()
for m in range(5,9):
    plt.xlabel('Latitude')
    plt.xlim([0,16])
    plt.ylim([935,990])
    #plt.xticks(np.arange(-90, 90, step=15))
    plt.ylabel( 'Surface Pressure (hPa)')
    plt.axvline(0, linewidth = 0.7, color = 'black',linestyle = '--')
    exec("plt.plot(lats,ps_z[m],label=month_list[m], color = clrs[m], linewidth = 1)") #alpha = 1 
    plt.axvspan(5,6.6,alpha = 0.02)
    plt.axvspan(10,12.2,alpha = 0.02)
plt.legend(loc= 'center left', bbox_to_anchor = (1,0.5), frameon = False)
plt.tight_layout()
#plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/00_Publication_figures/Reality_check/westafrica_ps.pdf", bbox_inches='tight')

# Surface pressure - zoom
fig = plt.figure()
for m in range(5,9):
    plt.xlabel('Latitude')
    plt.xlim([4,14])
    plt.ylim([940,980])
    #plt.xticks(np.arange(0, 17, step=5))
    plt.ylabel( 'Surface Pressure (hPa)')
    plt.axvline(0, linewidth = 0.7, color = 'black',linestyle = '--')
    exec("plt.plot(lats,ps_z[m],label=month_list[m], color = clrs[m], linewidth = 1)") #alpha = 1 
plt.legend(loc= 'center left', bbox_to_anchor = (1,0.5), frameon = False)
plt.tight_layout()
#plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/00_Publication_figures/Reality_check/westafrica_ps_zoom.pdf", bbox_inches='tight')

# Surface pressure - zoom and norm
fig = plt.figure()
for m in range(5,9):
    plt.xlabel('Latitude')
    plt.xlim([4,13])
    plt.ylim([-5,35])
    #plt.xticks(np.arange(0, 17, step=5))
    plt.ylabel( 'Surface Pressure (hPa)')
    plt.axvline(0, linewidth = 0.7, color = 'black',linestyle = '--')
    exec("plt.plot(lats,ps_z[m,:]-ps_z[m,166],label=month_list[m], color = clrs[m], linewidth = 1)") #alpha = 1 
plt.legend(loc= 'center left', bbox_to_anchor = (1,0.5), frameon = False)
plt.tight_layout()
#plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/00_Publication_figures/Reality_check/westafrica_ps_zoom_norm.pdf", bbox_inches='tight')


# Twin plot: Precipitation and surface pressure
for m in range(4,9):
    fig,ax = plt.subplots()
    plt.title(month_list[m])
    plt.xlabel("Latitude")
    ax.plot(lats,ps_z[m],color = 'darkorange')
    plt.ylim([935,1000])
    ax2=ax.twinx()
    ax2.plot(lats,pr_z[m],color = 'blue')
    plt.ylim([0,13])
    plt.axvline(5,color='black',linewidth = 0.7)    
    plt.xlim([0,17])
    ax2.set_ylabel(" Precipitation [mm]",color="blue")
    ax.set_ylabel("Surface Pressure [hPa]",color="darkorange")
    #plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/00_Publication_figures/Reality_check/westafrica_ps_pr_" + str(m) + ".pdf", bbox_inches='tight')

