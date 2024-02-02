#!/usr/bin/env python3
# -*- coding: utf-8 -*-


#-------------------------------------------------------------
###                     BIMODALITY              
# Anja Katzeberger, anja.katzenberger@pik-potsdam.de
#-------------------------------------------------------------

# This script creates the figures as displayed in the subchapter "Bimodality" and some more figures used in the context of the analysis
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
fn_50 = '/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/A_data/MonsoonPlanet_slab_50.nc'

ds_50 = nc.Dataset(fn_50)


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
wvp_50 = ds_50['WVP'][:]

evap_z_50 = evap_50.mean(axis=2)*86400
precip_z_50 = precip_50.mean(axis = 2)*86400
olr_z_50 = olr_50.mean(axis=2)
vcomp_z_50 = vcomp_50.mean(axis=3)
ucomp_z_50 = ucomp_50.mean(axis=3)
ps_z_50 = ps_50.mean(axis=2)
tsurf_z_50 = tsurf_50.mean(axis=2)
wvp_z_50 = wvp_50.mean(axis=2)


#%%
### Rainfall distribution

# Meridional rainfall distribution 
fig = plt.figure()
for m in range(0,12):
    plt.xlabel('Latitude')
    plt.xlim([-91,91])
    plt.ylim([-5,35])
    plt.xticks(np.arange(-90, 91, step=30))
    plt.ylabel( 'Precipitation (mm/day)')
    plt.axvline(10, linewidth = 0.7, color = 'black')
    plt.axvline(60, linewidth = 0.7, color = 'black')
    plt.axhline(0, linewidth = 0.7, color = 'black')
    exec("plt.plot(lats,precip_z_50[m],label=month_list[m], color = clrs[m], linewidth = 1)") #alpha = 1 
    #plt.legend(loc= 'center left', bbox_to_anchor = (1,0.5), frameon = False)
plt.tight_layout()
plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/00_Publication_figures/2-Bimodality/precip_10N.pdf", bbox_inches='tight')


#%%
### Surface pressure 
    
fig = plt.figure()
for m in range(4,9):
    plt.xlabel('Latitude')
    plt.xticks(np.arange(-90, 91, step=30))
    plt.ylabel( 'Surface pressure [hPa]')
    plt.axvline(10, linewidth = 0.7, color = 'black')
    plt.axvline(60, linewidth = 0.7, color = 'black')
    plt.axvline(0, linewidth = 0.7, color = 'black', linestyle = '--')
    plt.plot(lats,np.array(ps_z_50[m])/100-np.array(ps_z_50[m][49])/100,label=month_list[m], color = clrs[m], linewidth = 1) #alpha = 1 
    plt.xlim([-15,15])
    plt.xticks(np.arange(-15, 15, step=5))
   # plt.ylim([1001.5,1006])
    plt.ylim([-1.1,0.6])
#plt.legend(loc= 'center left', bbox_to_anchor = (1,0.5), frameon = False)
plt.tight_layout()
plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/00_Publication_figures/2-Bimodality/ps_10N_280ppm_norm.pdf', bbox_to_anchor = 'tight' )
    

#%%
### WVP

# Meridional WVP distribution 
fig = plt.figure()
for m in range(4,9):
    plt.xlabel('Latitude')
    plt.xlim([-5,20])
    plt.ylim([27,77])
    plt.xticks(np.arange(0, 20, step=5))
    plt.ylabel( 'Transportable water [mm]')
    plt.axvline(10, linewidth = 0.7, color = 'black')
    plt.axvline(60, linewidth = 0.7, color = 'black')
    plt.axhline(0, linewidth = 0.7, color = 'black')
    exec("plt.plot(lats,wvp_z_50[m],label=month_list[m], color = clrs[m], linewidth = 1)") #alpha = 1 
    #plt.legend(loc= 'center left', bbox_to_anchor = (1,0.5), frameon = False)
    plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/00_Publication_figures/2-Bimodality/wvp.pdf", bbox_inches='tight')



#%% 
### Relative barrier height/ left slope

# Calculate left slope 

# MAY
ps = ps_z_50[4]
maxl = argrelextrema(np.array(ps), np.greater)
minl = argrelextrema(np.array(ps), np.less)
ls_5 = ps[maxl[0][1]] - ps[minl[0][1]]
plt.plot(lats,ps, label = "JAN")
plt.xlim([-15,30])
plt.ylim([100000,101000])
plt.axhline(ps[maxl[0][1]], color = "red")
plt.axhline(ps[minl[0][1]], color = "blue")
plt.show()
plt.close()

# JUNE
ps = ps_z_50[5]
maxl = argrelextrema(np.array(ps), np.greater)
minl = argrelextrema(np.array(ps), np.less)
ls_6 = ps[maxl[0][1]] - ps[minl[0][1]]
plt.plot(lats,ps, label = "JAN")
plt.xlim([-15,30])
plt.ylim([100000,101000])
plt.axhline(ps[maxl[0][1]], color = "red")
plt.axhline(ps[minl[0][1]], color = "blue")
plt.show()
plt.close()

# JUL
ps = ps_z_50[6]
maxl = argrelextrema(np.array(ps), np.greater)
minl = argrelextrema(np.array(ps), np.less)
ls_7 = ps[maxl[0][1]] - ps[minl[0][1]]
plt.plot(lats,ps, label = "JAN")
plt.xlim([-15,30])
plt.ylim([100000,101000])
plt.axhline(ps[maxl[0][1]], color = "red")
plt.axhline(ps[minl[0][1]], color = "blue")
plt.show()
plt.close()

# AUG
ps = ps_z_50[7]
maxl = argrelextrema(np.array(ps), np.greater)
minl = argrelextrema(np.array(ps), np.less)
ls_8 = ps[maxl[0][1]] - ps[minl[0][1]]
plt.plot(lats,ps, label = "JAN")
plt.xlim([-15,30])
plt.ylim([100000,101000])
plt.axhline(ps[maxl[0][1]], color = "red")
plt.axhline(ps[minl[0][1]], color = "blue")
plt.show()
plt.close()

# SEP
ps = ps_z_50[8]
maxl = argrelextrema(np.array(ps), np.greater)
minl = argrelextrema(np.array(ps), np.less)
ls_9 = ps[maxl[0][1]] - ps[minl[0][1]]
plt.plot(lats,ps, label = "JAN")
plt.xlim([-15,30])
plt.ylim([100000,101000])
plt.axhline(ps[maxl[0][1]], color = "red")
plt.axhline(ps[minl[0][1]], color = "blue")
plt.show()
plt.close()

# Precipitation in monsoon region from May-September
# extracted via processing.py
precip_10N = [0.11082077617174946, 0.8364636974874884, 1.8688858661334962, 3.6348289693705738, 6.374520971439779]

plt.figure()
ls = [ls_5,ls_6,ls_7,ls_8,0]
fig,ax = plt.subplots()
plt.plot([5,6,7,8,9],ls, color = "darkorange")
plt.xlabel("Months")
ax.tick_params(axis='y', colors='darkorange')
plt.ylabel('Barrier height [Pa]', color = "darkorange")
ax2=ax.twinx()
ax2.plot([5,6,7,8,9],precip_10N, color = "darkblue")
plt.ylabel('Monsoon rainfall [mm/day]', color = "darkblue")
ax2.tick_params(axis='y', colors='darkblue')
plt.xlim([5,9])
plt.xticks(np.arange(5, 10, step=1))
#plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/barrier_height_monsoon.pdf')

# Monsoon rainfall per barrier height
plt.figure()
plt.plot(ls,precip_10N, '.-')
plt.ylabel("Monsoon rainfal [mm/day]")
plt.xlabel("Barrier height [Pa]")
# #plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/monsoon_per_barrier.pdf')


#%%
# Evolution of bimodal structure during annual cycle


cmap = matplotlib.cm.get_cmap("seismic_r")
rgba = cmap(0.5) # get middle value of color map

norm = matplotlib.colors.Normalize(vmin=-5, vmax=5) # translates vcomp value to a value between 0 and 1 
print(norm(0)) # 0.5


of = 1 #factor to avoid overlapping coloured background

    
# Meridional ps distribution
for m in range(1,13):
    vcompp_10 = vcomp_z_50[m-1]
    vcomppp_10 = vcompp_10[23] 
    fig,ax  =  plt.subplots()
    plt.title(month_list[m-1])
    plt.xlabel('Latitude')
    plt.ylabel( 'Surface pressure [hPa]', color = "darkorange")
    plt.ylim([995,1015])
    plt.xticks(np.arange(-30, 31, step=15))
    plt.axvline(10, linewidth = 0.7, color = 'black')
    plt.axvline(60, linewidth = 0.7, color = 'black')
    plt.axvline(0, linewidth = 0.7, color = 'black',linestyle = '--')
    plt.plot(lats,ps_z_50[m-1]/100,label=month_list[m-1], color = 'darkorange', linewidth = 1)
    plt.xlim([-20,30])
    ax.tick_params(axis='y', colors='darkorange')
    
    ax2 = ax.twinx()
    plt.bar(lats,precip_z_50[m-1], color = '#3b57b4', width = 1, alpha = 0.5)
    plt.ylabel("Precipitation [mm]", color = '#3b57b4') 
    plt.ylim([0,30])
    ax2.tick_params(axis='y', colors='#3b57b4') 
    
    for i in range(1,len(lats)-1): # leave out first and last step
        plt.axvspan(lats[i] - of*(lats[i]-lats[i-1])/2, lats[i] + of*(lats[i+1]-lats[i])/2, color=cmap(norm(vcomppp_10[i])), alpha = 0.1)

    #plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/3_Sensitivity_stripe_position/bimodality_10_' + str(m) + '.pdf' )



