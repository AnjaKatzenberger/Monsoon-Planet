#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 08:51:57 2022

@author: anjakatzenberger
"""


#-------------------------------------------------------------
###          ANALYSIS OF STRIPE POSITION                ###
#    Anja Katzeberger, anja.katzenberger@pik-potsdam.de
#-------------------------------------------------------------

# This script creates the figures as displayed in the subchapter "Sensitivty to different stripe positions" and some more figures used in the context of the analysis
# The directories of the input files as well as directories for saving the figures must be adapted. 

import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import statistics as stat


month_list = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
clrs = ['#6EA6CD','#98CAE1','#C2E4EF','#ACD39E','#5AAE61','#FDB366','#F67E4B','#DD3D2D','#A50026','#762A83','#36489A','#4A7BB7']

#%%
### LOADING DATA

latlow = [0,2,4,6,8,10,12,14,16,89]
latlowi = [0,2,4,6,8,10,12,14,16,99]
latlow2 = [0,2,4,6,8,10,12,14,16,'Aquaplanet']
latlow3 = ['0°N','2°N','4°N','6°N','8°N','10°N','12°N','14°N','16°N','Aquaplanet']
year = [71,41,41,71,41,91,91,91,41,41]

dir = []
for i in range(0,len(latlow)):
    print(latlow[i])
    
    if i == len(latlow)-1:
        dir.append("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/A_data/Aquaplanet.nc")
    else:
        dir.append("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/A_data/MonsoonPlanet_stripe_" + str(latlow[i]) +"N.nc")
    
    
    exec("ds_" + str(latlow[i]) + " = nc.Dataset(dir[i])")
    exec("evap_" + str(latlow[i]) + " = ds_" + str(latlow[i]) + "['evap'][:]")
    exec("precip_" + str(latlow[i]) + " = ds_" + str(latlow[i]) + "['precip'][:]")
    exec("tsurf_" + str(latlow[i]) + " = ds_" + str(latlow[i]) + "['t_surf'][:]")
    exec("tref_" + str(latlow[i]) + " = ds_" + str(latlow[i]) + "['t_ref'][:]")
    exec("ps_" + str(latlow[i]) + " = ds_" + str(latlow[i]) + "['ps'][:]")
    exec("olr_" + str(latlow[i]) + " = ds_" + str(latlow[i]) + "['olr'][:]")
    exec("ucomp_" + str(latlow[i]) + " = ds_" + str(latlow[i]) + "['ucomp'][:]")
    exec("vcomp_" + str(latlow[i]) + " = ds_" + str(latlow[i]) + "['vcomp'][:]")
    exec("sw_" + str(latlow[i]) + " = ds_" + str(latlow[i]) + "['swdn_toa'][:]")
    exec("swsur_" + str(latlow[i]) + " = ds_" + str(latlow[i]) + "['swdn_sfc'][:]") # shortwave down top of atmosphere
    exec("cld_" + str(latlow[i]) + " = ds_" + str(latlow[i]) + "['cld_amt'][:]") # cloud amount
    exec("cldt_" + str(latlow[i]) + " = ds_" + str(latlow[i]) + "['tot_cld_amt'][:]") # total cloud amount
    exec("rh_" + str(latlow[i]) + " = ds_" + str(latlow[i]) + "['rh'][:]") # relative humidity
    exec("sphum_" + str(latlow[i]) + " = ds_" + str(latlow[i]) + "['sphum'][:]") # specific humidity
    exec("wvp_" + str(latlow[i]) + " = ds_" + str(latlow[i]) + "['WVP'][:]") # column integrated water vapour



    exec("evap_z_" + str(latlow[i]) + " = evap_" + str(latlow[i]) + ".mean(axis=2)*86400")
    exec("precip_z_" + str(latlow[i]) + " = precip_" + str(latlow[i]) + ".mean(axis = 2)*86400")
    exec("olr_z_" + str(latlow[i]) + " = olr_" + str(latlow[i]) + ".mean(axis=2)")
    exec("vcomp_z_" + str(latlow[i]) + " = vcomp_" + str(latlow[i]) + ".mean(axis=3)")
    exec("ucomp_z_" + str(latlow[i]) + " = ucomp_" + str(latlow[i]) + ".mean(axis=3)")
    exec("ps_z_" + str(latlow[i]) + " = ps_" + str(latlow[i]) + ".mean(axis=2)/100")
    exec("tsurf_z_" + str(latlow[i]) + " = tsurf_" + str(latlow[i]) + ".mean(axis=2)-273.15")
    exec("tref_z_" + str(latlow[i]) + " = tref_" + str(latlow[i]) + ".mean(axis=2)-273.15")
    exec("sw_z_" + str(latlow[i]) + " = sw_" + str(latlow[i]) + ".mean(axis=2)")
    exec("swsur_z_" + str(latlow[i]) + " = swsur_" + str(latlow[i]) + ".mean(axis=2)")
    exec("cld_z_" + str(latlow[i]) + " = cld_" + str(latlow[i]) + ".mean(axis=3)")
    exec("cldt_z_" + str(latlow[i]) + " = cldt_" + str(latlow[i]) + ".mean(axis=2)")
    exec("rh_z_" + str(latlow[i]) + " = rh_" + str(latlow[i]) + ".mean(axis=3)")
    exec("sphum_z_" + str(latlow[i]) + " = sphum_" + str(latlow[i]) + ".mean(axis=3)")
    exec("wvp_z_" + str(latlow[i]) + " = wvp_" + str(latlow[i]) + ".mean(axis=2)")
                 

lats = ds_10.variables['lat'][:]
lons = ds_10.variables['lon'][:]
pfull = ds_10.variables['pfull'][:]


#%% 
### GENERAL CHARACTERISTICS

# Meridional rainfall distribution
for i in range(0,len(latlow)):
    fig = plt.figure()
    for m in range(5,9):
        plt.title(latlow3[i])
        plt.xlabel('latitude')
        plt.xlim([-20,20])
        plt.ylim([0,15])
        plt.xticks(np.arange(-20, 20, step=5))
        plt.ylabel( 'Precipitation (mm/day)')
        plt.axvline(latlowi[i], linewidth = 0.7, color = 'black')
        plt.axvline(latlow[i]+50, linewidth = 0.7, color = 'black')
        plt.axhline(0, linewidth = 0.7, color = 'black')
        exec("plt.plot(lats,precip_z_" + str(latlow[i]) + "[m],label=month_list[m], color = clrs[m], linewidth = 1)") #alpha = 1 
    plt.legend(loc= 'center left', bbox_to_anchor = (1,0.5),frameon = False)
    #plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/00_Publication_figures/3_Sensitivity_stripe_position/precip_' + str(latlow[i]) + '_JJAS.pdf', bbox_inches='tight' )

# Meridional rainfall distribution
for i in range(0,len(latlow)):
    fig = plt.figure()
    for m in range(5,9):
        plt.title(latlow3[i])
        plt.xlabel('latitude')
        plt.xlim([-20,20])
        plt.ylim([1000,1010])
        plt.xticks(np.arange(-20, 20, step=5))
        plt.ylabel( 'Precipitation (mm/day)')
        plt.axvline(latlowi[i], linewidth = 0.7, color = 'black')
        plt.axvline(latlow[i]+50, linewidth = 0.7, color = 'black')
        plt.axhline(0, linewidth = 0.7, color = 'black')
        exec("plt.plot(lats,ps_z_" + str(latlow[i]) + "[m],label=month_list[m], color = clrs[m], linewidth = 1)") #alpha = 1 
    plt.legend(loc= 'center left', bbox_to_anchor = (1,0.5),frameon = False)
    #plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/00_Publication_figures/3_Sensitivity_stripe_position/precip_' + str(latlow[i]) + '_JJAS.pdf', bbox_inches='tight' )
    

# Meridional precipitation distribution (contour)
for i in range(0,len(latlow)):
    plt.figure()
    if i == len(latlow)-1:
        plt.title(str(latlow2[i]))
    else:
        plt.title(str(latlow2[i])+'°N')
    exec("plt.contourf(range(1,13),lats,precip_z_" + str(latlow[i]) + ".transpose(),cmap='Blues',extend='max',levels = [0,2,4,6,8,10,12,14,16,18,20])")
    plt.axhline(0, linewidth = 0.7, color = 'black',linestyle = '--')
    plt.axhline(latlow[i], linewidth = 0.7, color = 'black')
    plt.axhline(60, linewidth = 0.7, color = 'black')
    plt.xlabel('Months')
    plt.xlim([1,12])
    plt.ylabel('Latitude')
    plt.ylim([-90,90])
    clb = plt.colorbar()
    clb.ax.set_title('mm/day',fontsize=8)
    #plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/3_Sensitivity_stripe_position/precip_global_" + str(latlow[i]) + "N_280ppm.pdf", bbox_inches='tight')

# Meridional evap distribution (evaporation)
for i in range(0,len(latlow)):
    plt.figure()
    plt.title(str(latlow[i])+'°N')
    exec("plt.contourf(range(1,13),lats,evap_z_" + str(latlow[i]) + ".transpose(),cmap='Blues',extend='max',levels = [0,1,2,3,4,5,6,7,8,9,10])")
    plt.axhline(0, linewidth = 0.7, color = 'black',linestyle = '--')
    plt.axhline(latlow[i], linewidth = 0.7, color = 'black')
    plt.axhline(60, linewidth = 0.7, color = 'black')
    plt.xlabel('Months')
    plt.xlim([1,12])
    plt.ylabel('Latitude')
    plt.ylim([-30,30])
    clb = plt.colorbar()
    clb.ax.set_title('mm/day',fontsize=8)
    #plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/3_Sensitivity_stripe_position/evap_" + str(latlow[i]) + "N_280ppm.pdf", bbox_inches='tight')

# Meridional tsurf distribution (contour)
for i in range(0,len(latlow)):
    plt.figure()
    plt.title(str(latlow[i])+'°N')
    exec("plt.contourf(range(1,13),lats,tsurf_z_" + str(latlow[i]) + ".transpose(),cmap='OrRd',extend='both',levels = range(0,30))")
    plt.axhline(0, linewidth = 0.7, color = 'black',linestyle = '--')
    plt.axhline(latlow[i], linewidth = 0.7, color = 'black')
    plt.axhline(60, linewidth = 0.7, color = 'black')
    plt.xlabel('Months')
    plt.xlim([1,12])
    plt.ylabel('Latitude')
    plt.ylim([-35,35])
    clb = plt.colorbar()
    clb.ax.set_title('°C',fontsize=8)
    #plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/3_Sensitivity_stripe_position/tsurf_' + str(latlow[i]) + 'N_280ppm.pdf')

# Surface pressure
for i in range(0,len(latlow)):
    plt.figure()
    plt.title(str(latlow[i])+'°N')
    exec("plt.contourf(range(1,13),lats,ps_z_" + str(latlow[i]) + ".transpose(),cmap='Greens',extend='both',levels = range(960,1030,2))")
    plt.axhline(0, linewidth = 0.7, color = 'black',linestyle = '--')
    plt.axhline(latlow[i], linewidth = 0.7, color = 'black')
    plt.axhline(60, linewidth = 0.7, color = 'black')
    plt.xlabel('Months')
    plt.xlim([1,12])
    plt.ylabel('Latitude')
    plt.ylim([-35,35])
    clb = plt.colorbar()
    clb.ax.set_title('hPa',fontsize=8)
    #plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/3_Sensitivity_stripe_position/ps_' + str(latlow[i]) + 'N_280ppm.pdf')

#%%   
### Pressure barrier dynamics
    

# Meridional ps distribution
for i in range(0,len(latlowi)):
    fig = plt.figure(figsize=(6,4))
    for m in range(5,10):
        if i == len(latlow)-1:
            plt.title(str(latlow2[i]))
        else:
            plt.title(str(latlow2[i])+'°N')
        plt.xlabel('Latitude')
        plt.xlim([-90,90])
        plt.ylim([975,1025])
        plt.xticks(np.arange(-90, 91, step=30))
        plt.ylabel( 'Surface pressure [hPa]')
        plt.axvline(latlow3[i], linewidth = 0.7, color = 'black')
       # plt.axvline(latlow3[i] + 50 , linewidth = 0.7, color = 'black')
        plt.axvline(0, linewidth = 0.7, color = 'black',linestyle = '--')
        exec("plt.plot(lats,ps_z_" + str(latlow[i]) + "[m],label=month_list[m], color = clrs[m], linewidth = 1)") #alpha = 1 
    plt.legend(loc= 'center left',frameon = False,bbox_to_anchor = (1,0.5))
    #plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/3_Sensitivity_stripe_position/ps_' + str(latlow4[i]) + '.pdf', bbox_inches='tight' )
      


# Meridional tsurf distribution

fig = plt.figure(figsize=(6,4))
for m in range(0,12):
    if i == len(latlow)-1:
        plt.title(str(latlow2[i]))
    else:
        plt.title(str(latlow2[i])+'°N')
    plt.xlabel('Latitude')
    plt.xlim([-88,88])
    #plt.ylim([975,1025])
    plt.xticks(np.arange(-90, 90, step=30),np.arange(-90, 90, step=30))
    plt.ylabel( 'Surface temperature [°C]')
    plt.axvline(latlow3[i], linewidth = 0.7, color = 'black')
   # plt.axvline(latlow3[i] + 50 , linewidth = 0.7, color = 'black')
    plt.axvline(0, linewidth = 0.7, color = 'black',linestyle = '--')
    exec("plt.plot(lats,tsurf_z_89[m],label=month_list[m], color = clrs[m], linewidth = 1)") #alpha = 1 
    plt.legend(loc= 'center left',frameon = False,bbox_to_anchor = (1,0.5))
    #plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/00_Publication_figures/3_Sensitivity_stripe_position/tsurf_aquaplanet.pdf', bbox_inches='tight' )
      

     
#%% 
### Hadley cell and ITCZ

#Calculate ITCZ 
for i in range(0,len(latlow)):
    exec("ITCZ_" + str(latlow[i]) + " = []")
    for m in range(0,12):
        exec("ITCZ_" + str(latlow[i]) + ".append(lats[np.argmax(precip_z_" + str(latlow[i]) + "[m])])")#-evap_z_" + str(latlow[i]) + "[m]

#Show ITCZ
plt.figure()
plt.xlabel("Months")
plt.ylabel("Latitude")
plt.axhline(0,linewidth = 0.7, linestyle = '--', color = "black")
for i in range(0,len(latlow)):
    exec("plt.plot(range(1,13),ITCZ_" + str(latlow[i]) + ",color = clrs[i],label = str(latlow[i]) + '°N')")
plt.legend(loc= 'center left',frameon = False,bbox_to_anchor = (1,0.5), title = "South Coast")
#plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/3_Sensitivity_stripe_position/ITCZ_stripes_280ppm.pdf')


# calculate upper and lower end of Hadley cell
for i in range(0,len(latlow)):
    exec("trans_nh_" + str(latlow[i]) + " = []") # northern hemisphere
    exec("lats_nh_" + str(latlow[i]) + " = []")
    exec("trans_sh_" + str(latlow[i]) + " = []") # southern hemisphere
    exec("lats_sh_" + str(latlow[i]) + " = []")
    for m in range(0,12):
        exec("olr = olr_z_" + str(latlow[i]) + "[m]")
        trans = []
        for s in range(0,len(olr)-1):
            if (olr[s] < 250 and olr[s+1] > 250) or (olr[s] > 250 and olr[s+1] < 250):
                trans.append(s)
                
        print(m+1)
        print(trans)
        exec("trans_nh_" + str(latlow[i]) + ".append(trans[-1])")
        exec("trans_sh_" + str(latlow[i]) + ".append(trans[0])")
        exec("lats_sh_" + str(latlow[i]) + ".append(lats[trans[0]])")
        exec("lats_nh_" + str(latlow[i]) + ".append(lats[trans[-1]])")
    
   
# OLR (plot) and marked ends of Hadley cells
# i= 5
# print(latlow[i])
# for m in range(0,12):
#     plt.figure()
#     plt.title(month_list[m])
#     exec("plt.plot(lats,olr_z_" + str(latlow[i]) + "[m])")
#     exec("plt.axvline(lats_sh_" + str(latlow[i]) + "[m], color = 'orange')")
#     exec("plt.axvline(lats_nh_" + str(latlow[i]) + "[m], color = 'orange')")    
    
# Hadley cell
color = ['purple','red','orange','yellow','green','darkblue', "darkred"]
plt.figure(figsize = (6,4))
f = 0
for i in range(0,len(latlow)):
    print(str(latlow[i])) 
    exec("plt.plot(range(1,13),lats_nh_" + str(latlow[i]) + ", color = clrs[i], linestyle = '--')")
    exec("plt.plot(range(1,13),lats_sh_" + str(latlow[i]) + ", color = clrs[i], linestyle = '--')")
    exec("plt.plot(range(1,13),ITCZ_" + str(latlow[i]) + ", color = clrs[i], label = str(latlow[i]) + '°N')")
    f = f+1
plt.xlabel("Months")
plt.xlim([1,12])
plt.ylim([-70,70])
plt.ylabel("Latitude")
plt.legend(frameon = False, bbox_to_anchor = (1,0.5), loc = "center left", title = "South Coast")
#plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/3_Sensitivity_stripe_position/slab_hadley_olr_pminuse.pdf', bbox_inches = 'tight')



#%%
### analysing differences between Aquaplanet and 16N 

# compare ITCZ position
dif = np.array(ITCZ_89)-np.array(ITCZ_16)

# Compare precipitation
for i in range(0,12):
    plt.figure()
    plt.title(month_list[i])
    plt.plot(lats,precip_z_89[i],color = 'red')
    plt.plot(lats,precip_z_16[i], color = 'blue', linestyle = '--')
    plt.axvline(0,color='black',linewidth = 0.7,linestyle = '--')
    plt.axvline(16,color='black',linewidth = 0.7,linestyle = '-')

# # Compare incoming shortwave radiation at the surface from 16N and Aquaplanet
# for i in range(0,12):
#     plt.figure()
#     plt.title(month_list[i])
#     plt.plot(lats,swsur_z_89[i],color = 'red')
#     plt.plot(lats,swsur_z_16[i], color = 'blue', linestyle = '--')
#     plt.axvline(0,color='black',linewidth = 0.7,linestyle = '--')

# # Compare surface temperature
# for i in range(0,12):
#     plt.figure()
#     plt.title(month_list[i])
#     plt.plot(lats,tsurf_z_89[i],color = 'red')
#     plt.plot(lats,tsurf_z_16[i], color = 'blue', linestyle = '--')
#     plt.axvline(0,color='black',linewidth = 0.7,linestyle = '--')

# # Compare surface temperature and incoming shortwave radiation at the surface/toa
# for i in range(0,12):
#     fig,ax = plt.subplots()
#     plt.title(month_list[i])
#     ax.plot(lats,sw_z_89[i],color = 'red')
#     ax.plot(lats,sw_z_16[i], color = 'red', linestyle = '--')
#     plt.ylim([100,500])
#     ax2=ax.twinx()
#     ax2.plot(lats,tref_z_89[i],color = 'darkorange')
#     ax2.plot(lats,tref_z_16[i], color = 'darkorange', linestyle = '--')
#     plt.ylim([10,35])
#     plt.axvline(0,color='black',linewidth = 0.7,linestyle = '--')   
#     plt.axvline(16,color='black',linewidth = 0.7,linestyle = '-')   
#     plt.xlim([-60,60])
#     ax2.set_ylabel("Surface temperature [°C]",color="darkorange")
#     ax.set_ylabel("Incoming SW surface [$W/m^2$]",color="red")

# # Compare surface temperature and precipition
# for i in range(0,12):
#     fig,ax = plt.subplots()
#     plt.title(month_list[i])
#     ax.plot(lats,precip_z_89[i],color = 'darkblue')
#     ax.plot(lats,precip_z_16[i], color = 'darkblue', linestyle = '--')
#     plt.ylim([0,35])
#     ax2=ax.twinx()
#     ax2.plot(lats,tsurf_z_89[i],color = 'darkorange')
#     ax2.plot(lats,tsurf_z_16[i], color = 'darkorange', linestyle = '--')
#     plt.ylim([15,35])
#     plt.axvline(0,color='black',linewidth = 0.7,linestyle = '--')   
#     plt.axvline(16,color='black',linewidth = 0.7,linestyle = '-')   
#     plt.xlim([-60,60])
#     ax2.set_ylabel("Surface temperature [°C]",color="darkorange")
#     ax.set_ylabel("Precipitation [mm/day]",color="darkblue")
    
# # Compare surface temperature and clouds tot 
# for i in range(0,12):
#     fig,ax = plt.subplots()
#     plt.title(month_list[i])
#     ax.plot(lats,cldt_z_89[i],color = 'darkblue')
#     ax.plot(lats,cldt_z_16[i], color = 'darkblue', linestyle = '--')
#     plt.ylim([0,100])
#     ax2=ax.twinx()
#     ax2.plot(lats,tsurf_z_89[i],color = 'darkorange')
#     ax2.plot(lats,tsurf_z_16[i], color = 'darkorange', linestyle = '--')
#     plt.ylim([15,35])
#     plt.axvline(0,color='black',linewidth = 0.7,linestyle = '--')   
#     plt.axvline(16,color='black',linewidth = 0.7,linestyle = '-')   
#     plt.xlim([-60,60])
#     ax2.set_ylabel("Surface temperature [°C]",color="darkorange")
#     ax.set_ylabel("Clouds",color="darkblue")

# # Compare surface temperature and wvp
# for i in range(0,12):
#     fig,ax = plt.subplots()
#     plt.title(month_list[i])
#     ax.plot(lats,wvp_z_89[i],color = 'darkblue')
#     ax.plot(lats,wvp_z_16[i], color = 'darkblue', linestyle = '--')
#     plt.ylim([0,100])
#     ax2=ax.twinx()
#     ax2.plot(lats,tsurf_z_89[i],color = 'darkorange')
#     ax2.plot(lats,tsurf_z_16[i], color = 'darkorange', linestyle = '--')
#     plt.ylim([15,35])
#     plt.axvline(0,color='black',linewidth = 0.7,linestyle = '--')   
#     plt.axvline(16,color='black',linewidth = 0.7,linestyle = '-')   
#     plt.xlim([-60,60])
#     ax2.set_ylabel("Surface temperature [°C]",color="darkorange")
#     ax.set_ylabel("Cloumn integrated water vapor",color="darkblue")

# # Compare surface temperature and relative humidity
# for i in range(0,12):
#     rh_z_89_i = rh_z_89[i]
#     rh_z_16_i = rh_z_16[i]
#     fig,ax = plt.subplots()
#     plt.title(month_list[i])
#     ax.plot(lats,sum(rh_z_89_i),color = 'darkblue')
#     ax.plot(lats,sum(rh_z_16_i), color = 'darkblue', linestyle = '--')
#    # plt.ylim([0,6])
#     ax2=ax.twinx()
#     ax2.plot(lats,tsurf_z_89[i],color = 'darkorange')
#     ax2.plot(lats,tsurf_z_16[i], color = 'darkorange', linestyle = '--')
#     plt.ylim([15,35])
#     plt.axvline(0,color='black',linewidth = 0.7,linestyle = '--')   
#     plt.axvline(16,color='black',linewidth = 0.7,linestyle = '-')   
#     plt.xlim([-60,60])
#     ax2.set_ylabel("Surface temperature [°C]",color="darkorange")
#     ax.set_ylabel("Relative humidity",color="darkblue")

# # Compare surface temperature and specific humidity
# for i in range(0,12):
#     sphum_z_89_i = sphum_z_89[i]
#     sphum_z_16_i = sphum_z_16[i]
#     fig,ax = plt.subplots()
#     plt.title(month_list[i])
#     ax.plot(lats,sum(sphum_z_89_i),color = 'darkblue')
#     ax.plot(lats,sum(sphum_z_16_i), color = 'darkblue', linestyle = '--')
#    # plt.ylim([0,6])
#     ax2=ax.twinx()
#     ax2.plot(lats,tsurf_z_89[i],color = 'darkorange')
#     ax2.plot(lats,tsurf_z_16[i], color = 'darkorange', linestyle = '--')
#     plt.ylim([15,35])
#     plt.axvline(0,color='black',linewidth = 0.7,linestyle = '--')   
#     plt.axvline(16,color='black',linewidth = 0.7,linestyle = '-')   
#     plt.xlim([-60,60])
#     ax2.set_ylabel("Surface temperature [°C]",color="darkorange")
#     ax.set_ylabel("Specific humidity",color="darkblue")

# # Compare surface pressure
# for i in range(0,12):
#     plt.figure()
#     plt.title(month_list[i])
#     plt.plot(lats,ps_z_89[i],color = 'red')
#     plt.plot(lats,ps_z_16[i], color = 'blue', linestyle = '--')
#     plt.axvline(0,color='black',linewidth = 0.7,linestyle = '--')

# # Compare meridional surface wind
# for i in range(0,12):
#     vcomp_z_89_i = vcomp_z_89[i]
#     vcomp_z_16_i = vcomp_z_16[i]
#     plt.figure()
#     plt.title(month_list[i])
#     plt.plot(lats,vcomp_z_89_i[23],color = 'red', label = '89')
#     plt.plot(lats,vcomp_z_16_i[23], color = 'blue', linestyle = '--', label = '16')
#     plt.axvline(0,color='black',linewidth = 0.7,linestyle = '--')
#     plt.axhline(0,color='black',linewidth = 0.7,linestyle = '-')
#     plt.ylim([-7,7])
#     plt.legend() 


#%%
### Monthly evolution of monsoon rain


# unit: mm/day = l m^(-2) day^(-1) 
# precip of latlow to landstripe, averaged over area
precip_0N = [0.08290639634651598, 0.2928551839431748, 1.1204153357539326, 2.159426937578246, 3.5442412365227938, 5.560370162129402, 6.452375859953463, 6.790469354018569, 6.937280693091452, 6.104630185291171, 0.7195206679170951, 0.17604656604817137]
precip_2N = [0.025207684848282952, 0.16705664093024097, 0.7882387057179585, 1.9245955569203943, 3.1473052105866373, 4.9919727724045515, 6.770646991208196, 7.4714233400300145, 7.695877086371183, 6.286217411980033, 0.4326584137743339, 0.0699011237884406]
precip_4N = [0.003263481063869259, 0.044823488037733006, 0.5697431951563404, 1.8395106145089597, 3.24384714639331, 5.373744255951438, 7.1606842175557555, 7.865640307171303, 8.133224930862587, 5.906198803499925, 0.2915878439165468, 0.021458222813895183]
precip_5N = [0.004682395399413508, 0.010514111181691987, 0.15549072450085077, 1.0017870721640065, 2.3102674051187932, 3.168900974560529, 4.885573103092611, 7.537623192183673, 8.521117782220244, 7.120785280130804, 0.8208966959500685, 0.018743236159934895]
precip_7N = [0.0016093022338736773, 0.0030903205242793774, 0.017194450765600777, 0.2584024776297156, 1.312448483076878, 2.1432288573123515, 3.400016133673489, 6.162905995734036, 8.297486929222941, 6.630139681510627, 0.6897941173519939, 0.01043863844643056]
precip_10N = [0.0014880401806749433, 4.1953710905318076e-05, 0.00028091526544926637, 0.00752542941881984, 0.11082077617174946, 0.8364636974874884, 1.8688858661334962, 3.6348289693705738, 6.374520971439779, 3.5414205165579915, 0.1328200894931797, 0.003222451834972162]
precip_11N = [0.00287709430040195, 0.0014346420755373401, 0.0006706978012971376, 0.0025519173334487277, 0.015007786987553118, 0.10866363772947807, 0.5637162015773356, 1.8511920876335353, 3.2603461761027575, 2.468434546608478, 0.2574966427346226, 0.014505243598250672]
precip_12N = [0.0005096683011629466, 6.597385198148231e-05, 5.2057711741326784e-05, 0.0006985881441323727, 0.010114334304489603, 0.0871271036885446, 0.3701335284858942, 1.4194883464369923, 3.152483026497066, 1.8284525081980973, 0.09900330442178529, 0.010588401528366376]
precip_13N = [0.00998725515728438, 0.0027367765142116696, 0.0005947767974134877, 0.0005388167352293749, 0.008106537802632374, 0.014780020637772395, 0.057953159557655454, 0.14628168792114593, 0.8805394114460796, 1.038077988778241, 0.07815148510417202, 0.026908448307949584]
precip_14N = [0.004286899002181599, 0.00027999568743553027, 2.5992975771771398e-05, 0.00017441517030647447, 0.0038318868973874487, 0.010095175434798875, 0.044783708108298015, 0.10386463800386991, 0.7663866912480444, 0.7635568559635431, 0.036094178449275205, 0.01886997520159639]
precip_15N = [0.018198444786321488, 0.007567096645288984, 0.0015383234654109401, 0.0005503962086095271, 0.004419212950779183, 0.008891054471860116, 0.015206752414087532, 0.03524914754962083, 0.2383291292062495, 0.22019991738488898, 0.05042119937570533, 0.045050232984067407]


# plot monthly averages over year for different configurations
plt.figure()
plt.plot(range(1,13),precip_0N, color = clrs[0],label  = '0°N')
plt.plot(range(1,13),precip_2N, color = clrs[1],label = "2°N")
plt.plot(range(1,13),precip_4N, color = clrs[2],label = "4°N")
plt.plot(range(1,13),precip_5N, color = clrs[3],label = "6°N")
plt.plot(range(1,13),precip_7N, color = clrs[4], label = "8°N")
plt.plot(range(1,13),precip_10N, color = clrs[5], label = "10°N")
plt.plot(range(1,13),precip_12N, color = clrs[6], label = "12°N")
plt.plot(range(1,13),precip_14N, color = clrs[7], label = "14°N")
plt.plot(range(1,13),precip_15N, color = clrs[8], label = "16°N")
plt.ylabel("Monsoon rainfall [mm/day]")
plt.xlim(1,12,1)
plt.xlabel("Months")
leg = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),frameon=False, title = "South coast")
leg._legend_box.align = "left"
#plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/3_Sensitivity_stripe_position/precip_season.pdf', bbox_inches='tight')
    

#%%
### In dependence of Land stripe position

# mr = monsoon region (10-30°N)
# JJAS = June-September

latlow2 = [0,2,4,6,8,10,12,14,16]

# The grid cell weighted average was calculated via processing.py

# from latlow to 30
precip89_mr2_annual=[4.9486633e-05,3.4891047e-05,2.3742454e-05,1.7618226e-05,1.4975964e-05,1.4016481e-05,1.3739592e-05,1.3731747e-05,1.3881478e-05]
precip_mr2_annual=[3.8446255e-05,3.8601502e-05,3.6998073e-05,3.3309545e-05,2.6212037e-05,1.6191263e-05,6.8018435e-06,1.6864690e-06,5.9399946e-07]

plt.figure(figsize = (6,4))
plt.plot(latlow2,np.array(precip_mr2_annual)*86400*365, label = 'Monsoon Planet')
plt.plot(latlow2,np.array(precip89_mr2_annual)*86400*365, linestyle = '-', label = 'Aquaplanet', color = 'darkred')
plt.axvspan(1.5,10.5,linewidth = 0.7, alpha = 0.1)
plt.axvspan(0,1.5,linewidth = 0.7, alpha = 0.1, color = "darkred")
plt.axvspan(10.5,16,linewidth = 0.7, alpha = 0.1, color = "darkred")
plt.xlabel("Latitude of equatorward coast [°N]")
plt.ylabel("Annual precipitation [mm]")
plt.xlim([0,16])
plt.legend(frameon = False)#, bbox_to_anchor = (1,0.5), loc = "center left")
#plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/3_Sensitivity_stripe_position/pr_per_land.pdf", bbox_inches='tight')

#####
# global precipitation for each month and each simulation
p14N =  [2.7034755796194077, 2.6245842163916677, 2.623438515001908, 2.5710851477924734, 2.4044542748015374, 2.372928784461692, 2.2896138194482774, 2.3050548566970974, 2.4937449605204165, 2.8128003235906363, 2.992598747368902, 2.8948649181984365]
p89N =  [3.5950981196947396, 3.6334619857370853, 3.701452841050923, 3.7524357670918107, 3.8122152909636497, 3.7349393824115396, 3.6884811124764383, 3.6781933740712702, 3.6982505349442363, 3.686141304206103, 3.7011602078564465, 3.6461266223341227]
p15N =  [2.790954045485705, 2.7444222243502736, 2.7100643841549754, 2.6913138572126627, 2.5078120990656316, 2.4074959626886994, 2.326228015590459, 2.3086533648893237, 2.4252511910162866, 2.7471615350805223, 2.984119614120573, 2.919180190656334]
p12N =  [2.5841153401415795, 2.5399491016287357, 2.4788560287561268, 2.407949371263385, 2.3517419525887817, 2.2923432290554047, 2.3154332768172026, 2.4369771068450063, 2.8003906016238034, 3.0190265737473965, 2.9070499003864825, 2.7673733420670033]
p02N =  [2.097973966738209, 2.031351038021967, 2.0994358754251152, 2.3299572814721614, 2.550054219318554, 2.824580459855497, 3.0449835467152297, 3.1859610811807215, 3.196975844912231, 3.0845584929920733, 2.620977535843849, 2.0917185000143945]
p05N =  [2.3003901704214513, 2.2005008824635297, 2.1140901662874967, 2.233013656223193, 2.534883812768385, 2.6974399806931615, 2.8710047830827534, 3.148884989786893, 3.2774307997897267, 3.2010051305405796, 2.7788790757767856, 2.3016664723400027]
p07N =  [2.393161179497838, 2.342344372300431, 2.2601017146371305, 2.213314507389441, 2.383782929973677, 2.5127310713287443, 2.621096820803359, 2.989908470772207, 3.245345188770443, 3.2143964781425893, 2.8566937311552465, 2.4680350441485643]
p10N =  [2.4828051624353975, 2.4315157730598003, 2.36140230554156, 2.300508040934801, 2.292630518786609, 2.34755597775802, 2.4446353898383677, 2.65898056095466, 3.0983213684521616, 3.17082084948197, 2.887167187873274, 2.6107360026799142]
p00N =  [2.0162853004876524, 1.997180376201868, 2.1029298717621714, 2.3379048972856253, 2.5656036974396557, 2.861376490909606, 3.0304958461783826, 3.1080719316378236, 3.1185702653601766, 2.9919182416051626, 2.605484949890524, 2.0575997012201697]

precip_mm_list = [p00N,p02N,p05N,p07N,p10N,p12N,p14N,p15N]

latlow_list = [0,2,5,7,10,12,14,15]
precip_mm_list_global = []
for p in precip_mm_list:
    precip_mm = p[0]*31 + p[1]*28 + p[2]*31 + p[3]*30 + p[4]*31 + p[5]*30 + p[6]*31 + p[7]*31 + p[8]*30 + p[9]*31 + p[10]*30 + p[11]*31
    precip_mm_list_global.append(precip_mm)

p89N_year_global = p89N[0]*31 + p89N[1]*28 + p89N[2]*31 + p89N[3]*30 + p89N[4]*31 + p89N[5]*30 + p89N[6]*31 + p89N[7]*31 + p89N[8]*30 + p89N[9]*31 + p89N[10]*30 + p89N[11]*31


#NH
p02N =  [0.519361330952961, 0.49263305263593793, 0.7672556326724589, 1.3090770720737055, 1.8426840368192643, 2.8782599489204586, 3.872717753984034, 4.3867325177416205, 4.575996729545295, 4.208645888138562, 1.3343482743948698, 0.7252218289067969]
p10N =  [0.4007744006230496, 0.16782156162662432, 0.17905088898260146, 0.3119492837868165, 0.5695118947187439, 0.9913031972246245, 1.5624555351678282, 2.870563475880772, 4.292107955552638, 4.584664141293615, 3.8598758401349187, 1.595712307607755]
p89N =  [4.393182706553489, 3.5257773706689477, 2.749512030277401, 2.306204487103969, 2.1253816899843514, 2.2695577587001026, 2.7716521988622844, 3.5640183370560408, 4.577416519168764, 5.005709559191018, 5.215514358133078, 4.974767134990543]
p15N =  [0.602694605186116, 0.34786346659529954, 0.2609134142403491, 0.285721287218621, 0.3475955469184555, 0.3414949227590114, 0.3217103767383378, 0.5313064859365113, 1.0679653147235513, 2.2554239840246737, 3.2307191868312657, 1.98573813540861]
p12N =  [0.659074900613632, 0.231320862803841, 0.19155683257849887, 0.2951001064502634, 0.4629246541298926, 0.5507923677214421, 1.0418137768283486, 2.192161465063691, 3.4778568777255714, 4.17875989805907, 3.966071200557053, 2.60296079213731]
p05N =  [0.362095505988691, 0.18553529953351244, 0.29626376344822347, 0.8053012483287603, 1.4663766080047935, 1.842021447373554, 2.9964799876324832, 4.304583999328315, 4.762615612708032, 4.670108633581549, 2.855560602620244, 0.6458513613324612]
p07N =  [0.3478927770629525, 0.1635698052268708, 0.1924646712723188, 0.4301324486732483, 1.0385245608631521, 1.2985851819394156, 2.103411569260061, 3.835047595202923, 4.647492268122733, 4.714148514904082, 3.581879020202905, 0.9168531949399039]
p14N =  [0.8435965137323365, 0.2988762847962789, 0.22578597781830467, 0.27057916959165595, 0.3701431938679889, 0.3783676089369692, 0.5223917818511836, 1.2024384399410337, 2.416716265724972, 3.4967274754308164, 4.042437148746103, 3.023242251947522]
p00N =  [0.5927251564571634, 0.5733958852943033, 0.9342592908069491, 1.441855455050245, 2.1372570656239986, 3.251897217705846, 3.8437304086983204, 4.18510950403288, 4.355911107268184, 3.960237081628293, 1.2309968151384965, 0.7941815012600273]


precip_mm_list = [p00N,p02N,p05N,p07N,p10N,p12N,p14N,p15N]

precip_mm_list_NH = []
for p in precip_mm_list:
    precip_mm = p[0]*31 + p[1]*28 + p[2]*31 + p[3]*30 + p[4]*31 + p[5]*30 + p[6]*31 + p[7]*31 + p[8]*30 + p[9]*31 + p[10]*30 + p[11]*31
    precip_mm_list_NH.append(precip_mm)

p89N_year_NH = p89N[0]*31 + p89N[1]*28 + p89N[2]*31 + p89N[3]*30 + p89N[4]*31 + p89N[5]*30 + p89N[6]*31 + p89N[7]*31 + p89N[8]*30 + p89N[9]*31 + p89N[10]*30 + p89N[11]*31


#SH
p02N =  [3.676586563233286, 3.5700690234079957, 3.4316161181777716, 3.350837412290275, 3.257424558978528, 2.7709006564691663, 2.217249182285741, 1.9851893302984536, 1.8179546459577978, 1.9604709406848997, 3.9076067972928286, 3.4582152497023344]
p14N =  [4.563354409765452, 4.950292187277228, 5.02109118970111, 4.8715911456383765, 4.438765277154744, 4.367489763535559, 4.056835896335542, 3.4076711162924767, 2.570773655315861, 2.128873171750456, 1.9427605031523854, 2.7664878987707198]
p05N =  [4.238684638403356, 4.215466347523034, 3.9319167262874544, 3.660725906956941, 3.6033908603712916, 3.5528586711734533, 2.7455295785330236, 1.993185665924102, 1.7922461440321058, 1.7319017846602947, 2.7021972346119583, 3.9574817405082285]
p10N =  [4.564835806377232, 4.695210023783147, 4.543753643520176, 4.2890668963082135, 4.015749064274132, 3.7038089940324426, 3.3268154016695917, 2.4473976460285485, 1.904534624191001, 1.756977871991694, 1.914458378450945, 3.625759540591389]
p00N =  [3.439845365937799, 3.4209647099487484, 3.2716004527173936, 3.2339541823603213, 2.9939503292553127, 2.4708560784347355, 2.217261283658445, 2.0310342020820826, 1.8812295806128532, 2.023599558742717, 3.9799730060622096, 3.3210180583409965]
p89N =  [2.7970132185146213, 3.741146600805223, 4.6533936518244445, 5.198666732758284, 5.499049206264317, 5.2003210061229765, 4.605310026090592, 3.7923684110864997, 2.819084550719708, 2.3665725777391344, 2.1868057432584465, 2.3174859525170177]
p15N =  [4.9792132107540965, 5.140981217846274, 5.1592156291008, 5.096906132530421, 4.668028769083321, 4.473497159779072, 4.330745595507324, 4.086000204551965, 3.782537067309022, 3.238898771815002, 2.737519727088511, 3.852622245904058]
p07N =  [4.4384295819327235, 4.521118733100593, 4.327738797292113, 3.996496566105634, 3.729041141923517, 3.726877039298415, 3.1387820723466575, 2.1447693463414907, 1.8431982665788382, 1.7146442842204124, 2.131508756428957, 4.019216971937567]
p12N =  [4.509155661799014, 4.84857716364786, 4.7661551856435835, 4.52079875394702, 4.240559251047671, 4.0338942082598805, 3.5890527768060565, 2.6817927486263216, 2.12292448268272, 1.8592934065964073, 1.8480287573765963, 2.931786049157381]

precip_mm_list = [p00N,p02N,p05N,p07N,p10N,p12N,p14N,p15N]

precip_mm_list_SH = []
for p in precip_mm_list:
    precip_mm = p[0]*31 + p[1]*28 + p[2]*31 + p[3]*30 + p[4]*31 + p[5]*30 + p[6]*31 + p[7]*31 + p[8]*30 + p[9]*31 + p[10]*30 + p[11]*31
    precip_mm_list_SH.append(precip_mm)

p89N_year_SH = p89N[0]*31 + p89N[1]*28 + p89N[2]*31 + p89N[3]*30 + p89N[4]*31 + p89N[5]*30 + p89N[6]*31 + p89N[7]*31 + p89N[8]*30 + p89N[9]*31 + p89N[10]*30 + p89N[11]*31




# global SH NH
latlow_list = [0,2,6,8,10,12,14,16]

# plot precip per latlow
fig = plt.figure(figsize = (6,4))
ax = plt.subplot(111)
plt.plot(latlow_list, precip_mm_list_global, "-b", label = "Global - Monsoon Planet", color = '#1f77b4')
plt.plot(latlow_list, precip_mm_list_NH, "--b", label = "NH", color = '#1f77b4')
plt.plot(latlow_list, precip_mm_list_SH,  ":b", label = "SH", color = '#1f77b4')
plt.axhline(p89N_year_global, color = "darkred", linestyle = "-", label = " Global - Aquaplanet")
plt.ylabel("Annual precipitation [mm]")
plt.xlabel("Latitude of equatorward coast [°N]")
plt.xlim(0,16)

# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])

# Put a legend to the right of the current axis
ax.legend(frameon = False)#,loc='center left', bbox_to_anchor=(1, 0.5))

#plt.legend()
fig.tight_layout()
#plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/3_Sensitivity_stripe_position/pr_global.pdf', bbox_inches='tight')
       

