#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#-------------------------------------------------------------
### ANALYSIS OF SLAB OCEAN DEPTH ###
# Anja Katzeberger, anja.katzenberger@pik-potsdam.de
#-------------------------------------------------------------

# This script creates the figures as displayed in the subchapter "Sensitivty to different slab ocean depths" and some more figures used in the context of the analysis
# The directories of the input files as well as directories for saving the figures must be adapted. 


import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from matplotlib.legend_handler import HandlerTuple


month_list = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
clrs = ['#6EA6CD','#98CAE1','#C2E4EF','#ACD39E','#5AAE61','#FDB366','#F67E4B','#DD3D2D','#A50026','#762A83','#36489A','#4A7BB7']
#%%
### Loading data

fn_500 ='/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/A_data/MonsoonPlanet_slab_500.nc'
fn_200 = '/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/A_data/MonsoonPlanet_slab_200.nc'
fn_50 = '/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/A_data/MonsoonPlanet_slab_50.nc'

ds_500 = nc.Dataset(fn_500)
ds_200 = nc.Dataset(fn_200)
ds_50 = nc.Dataset(fn_50)


lats = ds_500.variables['lat'][:]
lons = ds_500.variables['lon'][:]
pfull = ds_500.variables['pfull'][:]

evap_500 = ds_500['evap'][:]
precip_500 = ds_500['precip'][:]
tsurf_500 = ds_500['t_surf'][:]
ps_500 = ds_500['ps'][:]
olr_500 = ds_500['olr'][:]
vcomp_500 = ds_500['vcomp'][:]
ucomp_500 = ds_500['ucomp'][:]

evap_z_500 = evap_500.mean(axis=2)*86400
precip_z_500 = precip_500.mean(axis = 2)*86400
olr_z_500 = olr_500.mean(axis=2)
vcomp_z_500 = vcomp_500.mean(axis=3)
ucomp_z_500 = ucomp_500.mean(axis=3)
ps_z_500 = ps_500.mean(axis=2)
tsurf_z_500 = tsurf_500.mean(axis=2)

evap_200 = ds_200['evap'][:]
precip_200 = ds_200['precip'][:]
tsurf_200 = ds_200['t_surf'][:]
ps_200 = ds_200['ps'][:]
olr_200 = ds_200['olr'][:]
vcomp_200 = ds_200['vcomp'][:]
ucomp_200 = ds_200['ucomp'][:]

evap_z_200 = evap_200.mean(axis=2)*86400
precip_z_200 = precip_200.mean(axis = 2)*86400
olr_z_200 = olr_200.mean(axis=2)
vcomp_z_200 = vcomp_200.mean(axis=3)
ucomp_z_200 = ucomp_200.mean(axis=3)
ps_z_200 = ps_200.mean(axis=2)
tsurf_z_200 = tsurf_200.mean(axis=2)


evap_50 = ds_50['evap'][:]
precip_50 = ds_50['precip'][:]
tsurf_50 = ds_50['t_surf'][:]
ps_50 = ds_50['ps'][:]
olr_50 = ds_50['olr'][:]
vcomp_50 = ds_50['vcomp'][:]
ucomp_50 = ds_50['ucomp'][:]
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
### GENERAL CHARACTERISTICS OF THE SLAB RUNS

# Wind field for all months - slab50
vcomp_ad = np.squeeze(vcomp_z_50[:,:,50]) # 50 = 10N / See lats[50]
plt.figure()
plt.title('Meridional wind - slab 50m')
plt.contourf(range(1,13),pfull,vcomp_ad.transpose(),extend="both",cmap='RdBu',levels = [-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5])
plt.xlabel('Months')
plt.xlim([1,12])
plt.ylabel('Altitude [hPa]')
plt.ylim([1000,0])
clb = plt.colorbar()
clb.ax.set_title('m/s',fontsize=8)
#plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/4_Sensitivity_slab/vcomp_10N_280ppm.pdf')

# Wind field for all months - slab500
plt.figure()
vcomp_ad = np.squeeze(vcomp_z_500[:,:,50])
plt.title('Meridional wind - slab 500m')
plt.contourf(range(1,13),pfull,vcomp_ad.transpose(),extend="both",cmap='RdBu',levels = [-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5])
plt.xlabel('Months')
plt.xlim([1,12])
plt.ylabel('Altitude [hPa]')
plt.ylim([1000,0])
clb = plt.colorbar()
clb.ax.set_title('m/s',fontsize=8)
#plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/4_Sensitivity_slab/vcomp_10N_280ppm_500.pdf')




### Surface temperature (contour)
plt.figure()
plt.title('Slab 50m')
plt.contourf(range(1,13),lats,tsurf_z_50.transpose()-273.15,cmap="OrRd",extend="both",levels = range(0,50))
plt.axhline(0, linewidth = 0.7, color = 'black',linestyle = '--')
plt.axhline(10, linewidth = 0.7, color = 'black')
plt.axhline(60, linewidth = 0.7, color = 'black')
plt.xlabel('Months')
plt.xlim([1,12])
plt.ylabel('Latitude')
plt.ylim([-35,35])
clb = plt.colorbar()
clb.ax.set_title('°C',fontsize=8)
plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/00_Publication_figures/4_Sensitivity_slab/tsurf_10N_280ppm_50.pdf')

plt.figure()
plt.title('Slab 200m')
plt.contourf(range(1,13),lats,tsurf_z_200.transpose()-273.15,cmap="OrRd",extend="both",levels = range(0,50))
plt.axhline(0, linewidth = 0.7, color = 'black',linestyle = '--')
plt.axhline(10, linewidth = 0.7, color = 'black')
plt.axhline(60, linewidth = 0.7, color = 'black')
plt.xlabel('Months')
plt.xlim([1,12])
plt.ylabel('Latitude')
plt.ylim([-35,35])
clb = plt.colorbar()
clb.ax.set_title('°C',fontsize=8)
plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/00_Publication_figures/4_Sensitivity_slab/tsurf_10N_280ppm_200.pdf')

plt.figure()
plt.title('Slab 500m')
plt.contourf(range(1,13),lats,tsurf_z_500.transpose()-273.15,cmap="OrRd",extend="both",levels = range(0,50))
plt.axhline(0, linewidth = 0.7, color = 'black',linestyle = '--')
plt.axhline(10, linewidth = 0.7, color = 'black')
plt.axhline(60, linewidth = 0.7, color = 'black')
plt.xlabel('Months')
plt.xlim([1,12])
plt.ylabel('Latitude')
plt.ylim([-35,35])
clb = plt.colorbar()
clb.ax.set_title('°C',fontsize=8)
plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/00_Publication_figures/4_Sensitivity_slab/tsurf_10N_280ppm_500.pdf')



### Precipitation (contour)
plt.figure()
plt.title("Slab 50m")
plt.contourf(range(1,13),lats,precip_z_50.transpose(),cmap="Blues",extend="max",levels = [0,2,4,6,8,10,12,14,16,18,20])
plt.axhline(0, linewidth = 0.7, color = 'black',linestyle = '--')
plt.axhline(10, linewidth = 0.7, color = 'black')
plt.axhline(60, linewidth = 0.7, color = 'black')
plt.xlabel('Months')
plt.xlim([1,12])
plt.ylabel('Latitude')
plt.ylim([-35,35])
clb = plt.colorbar()
clb.ax.set_title('mm/day',fontsize=8)
#plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/4_Sensitivity_slab/precip_10N_280ppm_50.pdf')

plt.figure()
plt.title("Slab 200m")
plt.contourf(range(1,13),lats,precip_z_200.transpose(),cmap="Blues",extend="max",levels = [0,2,4,6,8,10,12,14,16,18,20])
plt.axhline(0, linewidth = 0.7, color = 'black',linestyle = '--')
plt.axhline(10, linewidth = 0.7, color = 'black')
plt.axhline(60, linewidth = 0.7, color = 'black')
plt.xlabel('Months')
plt.xlim([1,12])
plt.ylabel('Latitude')
plt.ylim([-35,35])
clb = plt.colorbar()
clb.ax.set_title('mm/day',fontsize=8)
#plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/4_Sensitivity_slab/precip_10N_280ppm_200.pdf')

plt.figure()
plt.title("Slab 500m")
plt.contourf(range(1,13),lats,precip_z_500.transpose(),cmap="Blues",extend="max",levels = [0,2,4,6,8,10,12,14,16,18,20])
plt.axhline(0, linewidth = 0.7, color = 'black',linestyle = '--')
plt.axhline(10, linewidth = 0.7, color = 'black')
plt.axhline(60, linewidth = 0.7, color = 'black')
plt.xlabel('Months')
plt.xlim([1,12])
plt.ylabel('Latitude')
plt.ylim([-35,35])
clb = plt.colorbar()
clb.ax.set_title('mm/day',fontsize=8)
#plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/4_Sensitivity_slab/precip_10N_280ppm_500.pdf')

# Meridional rainfall distribution
fig = plt.figure()
for m in range(0,12):
    plt.xlabel('latitude')
    plt.xlim([-91,91])
    plt.ylim([-5,35])
    plt.xticks(np.arange(-90, 91, step=30))
    plt.ylabel( 'Precipitation [mm/day]')
    plt.axvline(10, linewidth = 0.7, color = 'black')
    plt.axvline(60, linewidth = 0.7, color = 'black')
    plt.axhline(0, linewidth = 0.7, color = 'black')
    plt.plot(lats,precip_z_50[m],label=month_list[m], color = clrs[m], linewidth = 1) #alpha = 1 
plt.legend(loc= 'center left', frameon = False, bbox_to_anchor = (1,0.5))
#plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/4_Sensitivity_slab/precip_evap_10N_280ppm.pdf' )
    

#%%
### Pressure barrier dynamics


### Examine area between 1-7°N - surface temperature and surface winds
fig, axes = plt.subplots(nrows=1, ncols=2)
plt.subplot(1, 2, 1)
plt.xlim([1,7])
plt.ylim([22,32])
plt.xlabel("Latitudes")
plt.ylabel("Surface temperature [°C]")
for m in range(3,11):
    plt.plot(lats,tsurf_z_50[m-1]-273.15, color = clrs[m-1], label = month_list[m-1])
    plt.plot(lats,tsurf_z_500[m-1]-273.15, color = clrs[m-1], linestyle = '--')

plt.subplot(1, 2, 2)
plt.xlim([1,7])
plt.xlabel("Latitudes")
plt.ylabel("Meridonal surface winds [m/s]")
plt.axhline(0, linewidth = 1, color = "black")
for m in range(8,11):
    vcompp_50 = vcomp_z_50[m-1]
    vcompp_500 = vcomp_z_500[m-1]
    plt.plot(lats,vcompp_50[23], color = clrs[m-1])
    plt.plot(lats,vcompp_500[23], color = clrs[m-1], linestyle = '--')

# Put a legend to the right of the current axis
plt.subplot(1, 2, 2)
plt.figlegend(loc='center left', bbox_to_anchor=(1, 0.5), frameon = False)
plt.tight_layout()
#plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/4_Sensitivity_slab/slab_tsurf_vcomp_mar_oct.pdf', bbox_inches='tight')



### Figures as individuals 
# fig = plt.figure(figsize = (6,4))
# plt.xlim([1,7])
# plt.ylim([22,32])
# plt.xlabel("Latitudes")
# plt.ylabel("Surface temperature [°C]")
# for m in range(3,11):
#     plt.plot(lats,tsurf_z_50[m-1]-273.15, color = clrs[m-1], label = month_list[m-1])
#     plt.plot(lats,tsurf_z_500[m-1]-273.15, color = clrs[m-1], linestyle = '--')
# plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon = False)
# plt.tight_layout()
# #plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/4_Sensitivity_slab/slab_tsurf_vcomp_mar_oct1.pdf', bbox_inches='tight')

# fig = plt.figure(figsize = (6,4))
# plt.xlim([1,7])
# plt.xlabel("Latitudes")
# plt.ylabel("Meridonal surface winds [m/s]")
# plt.axhline(0, linewidth = 1, color = "black")
# for m in range(8,11):
#     vcompp_50 = vcomp_z_50[m-1]
#     vcompp_500 = vcomp_z_500[m-1]
#     plt.plot(lats,vcompp_50[23], color = clrs[m-1], label = month_list[m-1])
#     plt.plot(lats,vcompp_500[23], color = clrs[m-1], linestyle = '--')

# plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon = False)
# plt.tight_layout()
# #plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/4_Sensitivity_slab/slab_tsurf_vcomp_mar_oct2.pdf', bbox_inches='tight')


### Show how wind direction in 0-10°N affect the rainfall distribution during monsoon 
for m in range(1,13):
    vcompp_50 = vcomp_z_50[m-1]
    vcompp_500 = vcomp_z_500[m-1]
    
    plt.figure()
    fig,ax = plt.subplots()
    plt.title(month_list[m-1])
    plt.axvline(10,linewidth = 0.7, color = "black")
    plt.axvline(0,linewidth = 0.7, color = "black", linestyle = '--')
    plt.xlim([-15,30])
    plt.xlabel(["Latitudes"])
    plt.ylabel("Precipitation [mm/day]", color = "darkblue")
    ax.tick_params(axis='y', colors='darkblue') 
    plt.ylim([0,25])
    ax.plot(lats,precip_z_50[m-1], color = "darkblue")
    ax.plot(lats,precip_z_500[m-1], color = "darkblue",linestyle = '--')
    ax.axhline(0,linewidth = 0.7, color = "black")
    
    ax2=ax.twinx()
    ax2.plot(lats,vcompp_50[23], color = "mediumvioletred")
    ax2.plot(lats,vcompp_500[23], color = "mediumvioletred",linestyle = '--')
    ax2.set_ylim([-7,7])
    ax2.axhline(0,linewidth = 0.7, color = "black")
    ax2.tick_params(axis='y', colors='mediumvioletred')
    plt.ylabel("Meridional surface wind [m/s]", color = 'mediumvioletred')


# tsurf and ps distribution for all months 
for m in range(5,11):
    plt.figure()
    fig,ax = plt.subplots()
    plt.title(month_list[m-1])
    plt.axvline(10,linewidth = 1, color = "black")
    plt.axvline(0,linewidth = 1, color = "black", linestyle = '--')
    plt.xlim([-15,30])
    plt.xlabel("Latitudes")
    plt.ylabel("Surface temperature [°C]", color = "darkorange")
    ax.tick_params(axis='y', colors='darkorange') 
    plt.ylim([0,45])
    ax.plot(lats,tsurf_z_50[m-1]-273.15, color = "darkorange")
    ax.plot(lats,tsurf_z_500[m-1]-273.15, color = "darkorange",linestyle = '--')
    plt.axhline(0,color = "black", linewidth = 0.7)
    ax2=ax.twinx()
    ax2.plot(lats,ps_z_50[m-1]/100, color = "darkgreen")
    ax2.plot(lats,ps_z_500[m-1]/100, color = "darkgreen",linestyle = '--')
    ax2.set_ylim([995,1015])
    ax2.axhline(0,linewidth = 0.7, color = "black")
    ax2.tick_params(axis='y', colors='darkgreen')
    plt.ylabel("Surface pressure [hPa]", color = 'darkgreen')
    #plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/4_Sensitivity_slab/slab50_tsurf_ps_' + str(m) + '.pdf', bbox_inches='tight')
    

# Create figure with meridiodonal precip distribution, surface pressure and winds

cmap = matplotlib.cm.get_cmap("seismic_r")
rgba = cmap(0.5) # get middle value of color map

norm = matplotlib.colors.Normalize(vmin=-5, vmax=5) # translates vcomp value to a value between 0 and 1 
print(norm(0)) # 0.5

of = 1 #factor to avoid overlapping coloured background

for m in range(9,10):
    vcompp_50 = vcomp_z_50[m-1]
    vcompp_500 = vcomp_z_500[m-1]
        
    plt.figure(figsize = (6,4))
    fig,ax = plt.subplots()
    plt.title(month_list[m-1])
    plt.axvline(10,linewidth = 1, color = "black")
    plt.axvline(0,linewidth = 1, color = "black", linestyle = '--')
    plt.xlim([-15,30])
    plt.xlabel("Latitudes")
    plt.ylabel("Precipitation [mm/day]", color = "darkblue")
    ax.tick_params(axis='y', colors='darkblue') 
    plt.ylim([-10,25])
    ax.plot(lats,precip_z_50[m-1], color = "darkblue")
    vcomppp_50 = vcompp_50[23]   
    vcomppp_500 = vcompp_500[23] 
    for i in range(1,len(lats)-1): # leave out first and last step
        plt.axvspan(lats[i] - of*(lats[i]-lats[i-1])/2, lats[i] + of*(lats[i+1]-lats[i])/2, color=cmap(norm(vcomppp_50[i])), alpha = 0.3)
    plt.axhline(0,color = "black", linewidth = 0.7)
    ax2=ax.twinx()
    ax2.plot(lats,ps_z_50[m-1]/100, color = "darkgreen")
    ax2.set_ylim([995,1015])
    ax2.axhline(0,linewidth = 0.7, color = "black")
    ax2.tick_params(axis='y', colors='darkgreen')
    plt.ylabel("Surface pressure [hPa]", color = 'darkgreen')
    #plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/4_Sensitivity_slab/slab50_tsurf_ps_wind_' + str(m) + '.pdf', bbox_inches='tight')
    
    
    plt.figure(figsize = (6,4))
    fig,ax = plt.subplots()
    plt.title(month_list[m-1])
    plt.axvline(10,linewidth = 1, color = "black")
    plt.axvline(0,linewidth = 1, color = "black", linestyle = '--')
    plt.xlim([-15,30])
    plt.xlabel("Latitudes")
    plt.ylabel("Precipitation [mm/day]", color = "darkblue")
    ax.tick_params(axis='y', colors='darkblue') 
    plt.ylim([-10,25])
    ax.plot(lats,precip_z_500[m-1], color = "darkblue")
    vcomppp_50 = vcompp_50[23]   
    vcomppp_500 = vcompp_500[23] 
    for i in range(1,len(lats)-1): # leave out first and last step
        plt.axvspan(lats[i] - of*(lats[i]-lats[i-1])/2, lats[i] + of*(lats[i+1]-lats[i])/2, color=cmap(norm(vcomppp_500[i])), alpha = 0.3)
    plt.axhline(0,color = "black", linewidth = 0.7)
    ax2=ax.twinx()
    ax2.plot(lats,ps_z_500[m-1]/100, color = "darkgreen")
    ax2.set_ylim([995,1015])
    ax2.axhline(0,linewidth = 0.7, color = "black")
    ax2.tick_params(axis='y', colors='darkgreen')
    plt.ylabel("Surface pressure [hPa]", color = 'darkgreen')
    #plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/4_Sensitivity_slab/slab500_tsurf_ps_wind_' + str(m) + '.pdf', bbox_inches='tight')




#%%### 
### HADLEY CELL and ITCZ

# A for slab50: calculate ITCZ as P-E max
ITCZ_50 = []
for m in range(0,12):
    ITCZ_50.append(lats[np.argmax(precip_z_50[m]-evap_z_50[m])])


# A for slab200: calculate ITCZ as P-E max
ITCZ_200 = []
for m in range(0,12):
    ITCZ_200.append(lats[np.argmax(precip_z_200[m]-evap_z_200[m])])


# A for slab500: calculate ITCZ as P-E max
ITCZ_500 = []
for m in range(0,12):
    ITCZ_500.append(lats[np.argmax(precip_z_500[m]-evap_z_500[m])])



# calculate Hadley cell edges via OLR

# slab50:  calculate upper and lower end of Hadley cell (most poleward latitude where OLR = 0 )
trans_nh_50 = [] # northern hemisphere
lats_nh_50 = []
lats_sh_50 = []
trans_sh_50 = [] # southern hemisphere
for m in range(0,12):
    olr = olr_z_50[m]
    
    trans = []
    for i in range(0,len(olr)-1):
        if (olr[i] < 250 and olr[i+1] > 250) or (olr[i] > 250 and olr[i+1] < 250):
            trans.append(i)
            
    if m +1 != 8:   
        print(m+1)
        print(trans)
        trans_nh_50.append(trans[3]) 
        trans_sh_50.append(trans[0])
        lats_sh_50.append(lats[trans[0]])
        lats_nh_50.append(lats[trans[3]])
    else:
        print(m+1)
        print(len(trans))
        trans_nh_50.append(trans[5]) 
        trans_sh_50.append(trans[0])
        lats_sh_50.append(lats[trans[0]])
        lats_nh_50.append(lats[trans[5]])
        
# for m in range(0,12):
#     plt.figure()
#     plt.title(month_list[m])
#     plt.plot(lats,olr_z_50[m])
#     plt.axvline(lats_sh_50[m], color = "orange")
#     plt.axvline(lats_nh_50[m], color = "orange")
       
    




# B slab200:  calculate upper and lower end of Hadley cell (most poleward latitude where OLR = 0 )
trans_nh_200 = [] # northern hemisphere
lats_nh_200 = []
lats_sh_200 = []
trans_sh_200 = [] # southern hemisphere
for m in range(0,12):
    olr = olr_z_200[m]
    
    trans = []
    for i in range(0,len(olr)-1):
        if (olr[i] < 250 and olr[i+1] > 250) or (olr[i] > 250 and olr[i+1] < 250):
            trans.append(i)
            
    if m +1 != 9:   
        print(m+1)
        print(trans)
        trans_nh_200.append(trans[3]) 
        trans_sh_200.append(trans[0])
        lats_sh_200.append(lats[trans[0]])
        lats_nh_200.append(lats[trans[3]])
    else:
        print(m+1)
        print(len(trans))
        trans_nh_200.append(trans[5]) 
        trans_sh_200.append(trans[0])
        lats_sh_200.append(lats[trans[0]])
        lats_nh_200.append(lats[trans[5]])

# for m in range(0,12):
#     plt.figure()
#     plt.title(month_list[m])
#     plt.plot(lats,olr_z_200[m])
#     plt.axvline(lats_sh_200[m], color = "orange")
#     plt.axvline(lats_nh_200[m], color = "orange")




# B slab500:  calculate upper and lower end of Hadley cell (most poleward latitude where OLR = 0 )
trans_nh_500 = [] # northern hemisphere
lats_nh_500 = []
lats_sh_500 = []
trans_sh_500 = [] # southern hemisphere
for m in range(0,12):
    olr = olr_z_500[m]
    print(len(olr))
    
    trans = []
    for i in range(0,len(olr)-1):
        if (olr[i] < 250 and olr[i+1] > 250) or (olr[i] > 250 and olr[i+1] < 250):
            trans.append(i)
            
    print(trans)     
    trans_nh_500.append(trans[3]) 
    trans_sh_500.append(trans[0])
    lats_sh_500.append(lats[trans[0]])
    lats_nh_500.append(lats[trans[3]])
        
# for m in range(0,12):
#     plt.figure()
#     plt.title(month_list[m])
#     plt.plot(lats,olr_z_500[m])
#     plt.axvline(lats_sh_500[m], color = "orange")
#     plt.axvline(lats_nh_500[m], color = "orange")
       
    
    
plt.figure()
plt.plot(range(1,13),lats_nh_50, color = "darkblue", linestyle = '--')
plt.plot(range(1,13),lats_sh_50, color = "darkblue", linestyle = '--')   
plt.plot(range(1,13),ITCZ_50, color = "darkblue", label = "50m")
plt.plot(range(1,13),ITCZ_200,color = "darkgreen", label = "200m")
plt.plot(range(1,13),lats_nh_200,color = "darkgreen", linestyle = '--')
plt.plot(range(1,13),lats_sh_200, color = "darkgreen", linestyle = '--')
plt.plot(range(1,13),ITCZ_500,color = "darkorange", label = "500m")
plt.plot(range(1,13),lats_nh_500,color = "darkorange", linestyle = '--')
plt.plot(range(1,13),lats_sh_500, color = "darkorange", linestyle = '--')
plt.xlabel("Months")
plt.xlim([1,12])
plt.ylim([-70,70])
plt.ylabel("Latitude")
leg = plt.legend(frameon = False, title = "Slab Depth",loc='center left', bbox_to_anchor=(1, 0.5))
leg._legend_box.align = "left"
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/4_Sensitivity_slab/slab_hadley_olr_pminuse.pdf', bbox_inches='tight')
   



#%%
### Monthly evolution of monsoon rain

# The grid cell weighted average was calculated via processing.py

# Monthly monsoon rainfall for each of the simulations
precip_m = [[1.3660164e-08,4.7876936e-10,3.1295400e-09,7.9958681e-08,1.1295950e-06,9.1309057e-06,2.1973536e-05,4.3431541e-05,7.4565265e-05,4.2340656e-05,1.5969730e-06,2.9456896e-08],[1.0326717e-08,9.7325825e-10,5.4131948e-09,9.3833414e-08,6.7169202e-07,1.8932678e-06,1.3223824e-05,2.5032767e-05,3.2737309e-05,1.8370407e-05,3.4031908e-07,7.6871132e-08],[2.2326665e-08,8.8336654e-09,1.0766793e-08,2.0907267e-07,1.0199841e-06,2.4448550e-06,1.2612470e-05,2.2486911e-05,2.4886969e-05,1.5546675e-05,2.3946666e-07,3.7047815e-08]]
colorss = ["darkblue", "darkgreen", "darkorange"]
slab = [50,200,500]

plt.figure()
for i in range(0,len(precip_m)):
    plt.plot(range(1,13),np.array(precip_m[i])*86400, label = '' + str(slab[i]) + 'm', color = colorss[i])
plt.ylabel("Monsoon Rainfall [mm/day]")
plt.xlabel('Months')
plt.ylim([0,7])
plt.axhline(0,color="black",linewidth = 0.7)
plt.legend(frameon = False, title = "Slab Depth", loc='center left', bbox_to_anchor=(1, 0.5))
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/4_Sensitivity_slab/pr_per_month_slab.pdf", bbox_inches='tight')
  
    
  
#%%
### In dependence of slab depth

# mr = monsoon region (10-30°N)
# JJAS = June-September

# The grid cell weighted average was calculated via processing.py

slab = [50,200,500]
precip_mr_annual = [1.6191263e-05,7.7047498e-06,6.6271145e-06]
precip_global_annual = [3.0008980e-05,2.9380648e-05,2.9112593e-05]
tsurf_global_annual = [290.76,291.99,292.09]
tsurf_mr_annual = [300.95865,303.82892,304.17676]

fig,ax = plt.subplots()
pr_mr, = plt.plot(slab,np.array(precip_mr_annual)*86400*365, color = "darkblue",marker = ".")
pr_g, = plt.plot(slab,np.array(precip_global_annual)*86400*365, color = "darkblue",linestyle = '--',marker = ".")
plt.xlabel("Slab ocean depth [m]")
plt.ylabel("Annual precipitation [mm]", color = "darkblue")
ax.tick_params(axis='y', colors='darkblue')
plt.legend(frameon = False,loc = "center left")
ax2=ax.twinx()
tsurf_mr, = ax2.plot(slab,np.array(tsurf_mr_annual)-273.15, color = 'darkorange', label = "Monsoon region",marker = ".")
tsurf_g, = ax2.plot(slab,np.array(tsurf_global_annual)-273.15, color = 'darkorange',linestyle = '--', label = "Global",marker = ".")
ax2.set_ylabel("Surface Temperature [°C]", color = "darkorange")
ax2.tick_params(axis='y', colors='darkorange')
plt.legend([(pr_g, tsurf_g), (pr_mr, tsurf_mr)], ['Global', 'Monsoon region'], numpoints=1,handler_map={tuple: HandlerTuple(ndivide=None)}, handlelength=3, frameon = False)
#plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/4_Sensitivity_slab/pr_per_slab.pdf", bbox_inches='tight')



    
    
    
    
    
    
    
    
    
