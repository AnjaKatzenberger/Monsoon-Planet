#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#-------------------------------------------------------------
### ANALYSIS OF CO2 SIMULATIONS ###
# Anja Katzeberger, anja.katzenberger@pik-potsdam.de
#-------------------------------------------------------------

# This script creates the figures as displayed in the subchapter "Sensitivty to changes in CO2 Concentration" and some more figures used in the context of the analysis
# The directories of the input files as well as directories for saving the figures must be adapted. 


import math
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import statistics as stat
from matplotlib.legend_handler import HandlerTuple


#%%
### Loading data 

ppm = [70,140,200,280,400,560,750,1120]

dir = []
for i in range(0,len(ppm)):
    print(i)
    dir.append("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/A_data/MonsoonPlanet_co2_" + str(ppm[i]) + "ppm.nc")
    print(dir[i])
    exec("ds_" + str(ppm[i]) + " = nc.Dataset(dir[i])")
    
    exec("evap_" + str(ppm[i]) + " = ds_" + str(ppm[i]) + "['evap'][:]")
    exec("precip_" + str(ppm[i]) + " = ds_" + str(ppm[i]) + "['precip'][:]")
    exec("tsurf_" + str(ppm[i]) + " = ds_" + str(ppm[i]) + "['t_surf'][:]")
    exec("tref_" + str(ppm[i]) + " = ds_" + str(ppm[i]) + "['t_ref'][:]")
    exec("ps_" + str(ppm[i]) + " = ds_" + str(ppm[i]) + "['ps'][:]")
    exec("olr_" + str(ppm[i]) + " = ds_" + str(ppm[i]) + "['olr'][:]")
    exec("ucomp_" + str(ppm[i]) + " = ds_" + str(ppm[i]) + "['ucomp'][:]")
    exec("vcomp_" + str(ppm[i]) + " = ds_" + str(ppm[i]) + "['vcomp'][:]")
    exec("sw_" + str(ppm[i]) + " = ds_" + str(ppm[i]) + "['swdn_toa'][:]")
    exec("swsur_" + str(ppm[i]) + " = ds_" + str(ppm[i]) + "['swdn_sfc'][:]") # shortwave down top of atmosphere
    exec("cld_" + str(ppm[i]) + " = ds_" + str(ppm[i]) + "['cld_amt'][:]") # cloud amount
    exec("cldt_" + str(ppm[i]) + " = ds_" + str(ppm[i]) + "['tot_cld_amt'][:]") # total cloud amount
    exec("rh_" + str(ppm[i]) + " = ds_" + str(ppm[i]) + "['rh'][:]") # relative humidity
    exec("sphum_" + str(ppm[i]) + " = ds_" + str(ppm[i]) + "['sphum'][:]") # specific humidity
    exec("wvp_" + str(ppm[i]) + " = ds_" + str(ppm[i]) + "['WVP'][:]") # column integrated water vapour



    exec("evap_z_" + str(ppm[i]) + " = evap_" + str(ppm[i]) + ".mean(axis=2)*86400")
    exec("precip_z_" + str(ppm[i]) + " = precip_" + str(ppm[i]) + ".mean(axis = 2)*86400")
    exec("olr_z_" + str(ppm[i]) + " = olr_" + str(ppm[i]) + ".mean(axis=2)")
    exec("vcomp_z_" + str(ppm[i]) + " = vcomp_" + str(ppm[i]) + ".mean(axis=3)")
    exec("ucomp_z_" + str(ppm[i]) + " = ucomp_" + str(ppm[i]) + ".mean(axis=3)")
    exec("ps_z_" + str(ppm[i]) + " = ps_" + str(ppm[i]) + ".mean(axis=2)/100")
    exec("tsurf_z_" + str(ppm[i]) + " = tsurf_" + str(ppm[i]) + ".mean(axis=2)-273.15")
    exec("tref_z_" + str(ppm[i]) + " = tref_" + str(ppm[i]) + ".mean(axis=2)-273.15")
    exec("sw_z_" + str(ppm[i]) + " = sw_" + str(ppm[i]) + ".mean(axis=2)")
    exec("swsur_z_" + str(ppm[i]) + " = swsur_" + str(ppm[i]) + ".mean(axis=2)")
    exec("cld_z_" + str(ppm[i]) + " = cld_" + str(ppm[i]) + ".mean(axis=3)")
    exec("cldt_z_" + str(ppm[i]) + " = cldt_" + str(ppm[i]) + ".mean(axis=2)")
    exec("rh_z_" + str(ppm[i]) + " = rh_" + str(ppm[i]) + ".mean(axis=3)")
    exec("sphum_z_" + str(ppm[i]) + " = sphum_" + str(ppm[i]) + ".mean(axis=3)")
    exec("wvp_z_" + str(ppm[i]) + " = wvp_" + str(ppm[i]) + ".mean(axis=2)")
                 

lats = ds_280.variables['lat'][:]
lons = ds_280.variables['lon'][:]
pfull = ds_280.variables['pfull'][:]


month_list = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
clrs = ['#6EA6CD','#98CAE1','#C2E4EF','#ACD39E','#5AAE61','#FDB366','#F67E4B','#DD3D2D','#A50026','#762A83','#36489A','#4A7BB7']


#%%
### General plots


# Meridional precipitation distribution (contour)
for i in range(0,len(ppm)):
    plt.figure()
    plt.title(str(ppm[i])+' $ppm$')
    exec("plt.contourf(range(1,13),lats,precip_z_" + str(ppm[i]) + ".transpose(),cmap='Blues',extend='max',levels = [0,2,4,6,8,10,12,14,16,18,20])")
    plt.axhline(0, linewidth = 0.7, color = 'black',linestyle = '--')
    plt.axhline(10, linewidth = 0.7, color = 'black')
    plt.axhline(60, linewidth = 0.7, color = 'black')
    plt.xlabel('Months')
    plt.xlim([1,12])
    plt.ylabel('Latitude')
    plt.ylim([-35,35])
    clb = plt.colorbar()
    clb.ax.set_title('mm/day',fontsize=8)
    #plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/6_Sensitivity_ppm/precip_10N_" + str(ppm[i]) + "ppm_contour.pdf", bbox_inches='tight')

# Meridional rainfall distribution 
for i in range(0,len(ppm)):
    fig = plt.figure()
    for m in range(0,12):
        plt.title(str(ppm[i])+' $ppm$')
        plt.xlabel('latitude')
        plt.xlim([-91,91])
        plt.ylim([-5,35])
        plt.xticks(np.arange(-90, 91, step=30))
        plt.ylabel( 'Precipitation (mm/day)')
        plt.axvline(10, linewidth = 0.7, color = 'black')
        plt.axvline(60, linewidth = 0.7, color = 'black')
        plt.axhline(0, linewidth = 0.7, color = 'black')
        exec("plt.plot(lats,precip_z_" + str(ppm[i]) + "[m],label=month_list[m], color = clrs[m], linewidth = 1)") #alpha = 1 
    plt.legend(loc= 'center left', bbox_to_anchor = (1,0.5), frameon = False)
    #plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/6_Sensitivity_ppm/precip_10N_" + str(ppm[i]) + "ppm.pdf", bbox_inches='tight')


# Meridional tsurf distribution (contour)
for i in range(0,len(ppm)):
    plt.figure()
    plt.title(str(ppm[i])+'ppm')
    exec("plt.contourf(range(1,13),lats,tsurf_z_" + str(ppm[i]) + ".transpose(),cmap='OrRd',extend='both',levels = range(0,50))")
    plt.axhline(0, linewidth = 0.7, color = 'black',linestyle = '--')
    plt.axhline(10, linewidth = 0.7, color = 'black')
    plt.axhline(60, linewidth = 0.7, color = 'black')
    plt.xlabel('Months')
    plt.xlim([1,12])
    plt.ylabel('Latitude')
    plt.ylim([-35,35])
    clb = plt.colorbar()
    clb.ax.set_title('°C',fontsize=8)
    plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/00_Publication_figures/6_Sensitivity_ppm/tsurf_10N_" + str(ppm[i]) + "ppm_contour.pdf", bbox_inches='tight')

# Meridional evaporation distribution
for i in range(0,len(ppm)):
    fig = plt.figure()
    for m in range(0,12):
        plt.title(str(ppm[i])+'ppm')
        plt.xlabel('latitude')
        plt.xlim([-91,91])
        plt.ylim([-5,35])
        plt.xticks(np.arange(-90, 91, step=30))
        plt.ylabel( 'Evaporation (mm/day)')
        plt.axvline(10, linewidth = 0.7, color = 'black')
        plt.axvline(60, linewidth = 0.7, color = 'black')
        plt.axhline(0, linewidth = 0.7, color = 'black')
        exec("plt.plot(lats,evap_z_" + str(ppm[i]) + "[m],label=month_list[m], color = clrs[m], linewidth = 1)") #alpha = 1 
    plt.legend(loc= 'center left', bbox_to_anchor = (1,0.5), frameon = False)
    #plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/6_Sensitivity_ppm/evap_10N_" + str(ppm[i]) + "ppm.pdf", bbox_inches='tight')
    
    
    
 # Incoming shortwave radiation at the surface for selected months
plt.figure()
plt.plot(lats,swsur_z_140[-1], label = "140 ppm",linestyle = '--', color = "blue")
plt.plot(lats,swsur_z_280[8-1], label = "280 ppm",linestyle = '--', color = "green")
plt.plot(lats,swsur_z_560[8-1], label = "560 ppm",linestyle = '--', color = "darkorange")
#plt.axvline(0,linewidth = 0.7,color="black",linestyle = '--')
plt.axhline(0,linewidth = 0.7,color="black")
#plt.ylim([100,300])
plt.xlim([-20,20])
plt.xlabel("Latitudes")
plt.ylabel("Incoming SW radiation - surface [$ppm$]")
plt.legend(frameon = False)
#plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/2_General_characteristics_10N_280ppm/sw_10N_280ppm.pdf')


# Clouds
plt.figure()
plt.plot(lats,cldt_z_140[-1], label = "140 ppm", color = "blue")
plt.plot(lats,cldt_z_280[8-1], label = "280 ppm", color = "green")
plt.plot(lats,cldt_z_560[8-1], label = "560 ppm", color = "darkorange")
#plt.axvline(0,linewidth = 0.7,color="black",linestyle = '--')
plt.axhline(0,linewidth = 0.7,color="black")
#plt.ylim([100,300])
plt.xlim([-20,20])
plt.xlabel("Latitudes")
plt.ylabel("Incoming SW radiation - surface [$ppm$]")
plt.legend(frameon = False)
#plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/2_General_characteristics_10N_280ppm/sw_10N_280ppm.pdf')

    
# Surface pressure and clouds
fig, ax = plt.subplots()
plt.plot(lats,ps_z_140[-1], label = "140 ppm", color = "blue")
plt.plot(lats,ps_z_280[8-1], label = "280 ppm", color = "green")
plt.plot(lats,ps_z_560[8-1], label = "560 ppm", color = "darkorange")
#plt.axvline(0,linewidth = 0.7,color="black",linestyle = '--')
plt.axhline(0,linewidth = 0.7,color="black")
plt.ylim([980,1040])
plt.xlim([-20,20])
plt.xlabel("Latitudes")
plt.ylabel("Incoming SW radiation - surface [$ppm$]")
plt.legend(frameon = False)

ax2 = ax.twinx()
ax2.plot(lats,cldt_z_140[8-1], label = "140 ppm", color = "blue")
ax2.plot(lats,cldt_z_280[8-1], label = "280 ppm", color = "green")
ax2.plot(lats,cldt_z_560[8-1], label = "560 ppm", color = "darkorange")
#plt.axvline(0,linewidth = 0.7,color="black",linestyle = '--')
plt.axhline(0,linewidth = 0.7,color="black")
#plt.ylim([100,300])
plt.xlim([-20,20])
plt.xlabel("Latitudes")
plt.ylabel("Incoming SW radiation - surface [$ppm$]")
plt.legend(frameon = False)
#plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/2_General_characteristics_10N_280ppm/sw_10N_280ppm.pdf')
   
    
    
# August surface winds (meridional)
plt.figure()
plt.title("August")
vcompp_140 = vcomp_z_140[8-1]
vcompp_280 = vcomp_z_280[8-1]
vcompp_560 = vcomp_z_560[8-1]
vcompp_1120 = vcomp_z_1120[8-1]
plt.plot(lats,vcompp_140[23], label = "140 ppm", color = clrs[1])
plt.plot(lats,vcompp_280[23], label = "280 ppm", color = clrs[3])
plt.plot(lats,vcompp_560[23], label = "560 ppm", color = clrs[5])
plt.plot(lats,vcompp_1120[23], label = "1120 ppm", color = clrs[7])
plt.legend(frameon = False,loc='center left', bbox_to_anchor=(1, 0.5), title = "Carbon Dioxid")
plt.xlim([-15,30])
plt.xlabel("Latitude")
plt.ylabel("August surface winds [m/s]")
plt.axvline(0, color = "black", linewidth = 0.7, linestyle = '--')
plt.axvline(10, color = "black", linewidth = 0.7)
plt.axhline(0,color = "black", linewidth = 0.7)
#plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/6_Sensitivity_ppm/vcomp_aug.pdf", bbox_inches='tight')

    
# September surface winds (meridional)
plt.figure()
plt.title("September")
vcompp_140 = vcomp_z_140[9-1]
vcompp_280 = vcomp_z_280[9-1]
vcompp_560 = vcomp_z_560[9-1]
vcompp_1120 = vcomp_z_1120[9-1]
plt.plot(lats,vcompp_140[23], label = "140 ppm", color = clrs[1])
plt.plot(lats,vcompp_280[23], label = "280 ppm", color = clrs[3])
plt.plot(lats,vcompp_560[23], label = "560 ppm", color = clrs[5])
plt.plot(lats,vcompp_1120[23], label = "1120 ppm", color = clrs[7])
plt.legend(frameon = False,loc='center left', bbox_to_anchor=(1, 0.5), title = "Carbon Dioxid")
plt.xlim([-15,30])
plt.xlabel("Latitude")
plt.ylabel("September surface winds [m/s]")
plt.axvline(0, color = "black", linewidth = 0.7, linestyle = '--')
plt.axvline(10, color = "black", linewidth = 0.7)
plt.axhline(0,color = "black", linewidth = 0.7)
#plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/6_Sensitivity_ppm/vcomp_sep.pdf", bbox_inches='tight')




#%%
### Pressure barrier dynamics

# August surface pressure
plt.figure()
plt.plot(lats,ps_z_140[8-1], label = "140 ppm", color = clrs[1])
plt.plot(lats,ps_z_280[8-1], label = "280 ppm", color = clrs[3])
plt.plot(lats,ps_z_560[8-1], label = "560 ppm", color = clrs[5])
plt.plot(lats,ps_z_1120[8-1], label = "1120 ppm", color = clrs[7])
plt.legend(frameon = False,loc='center left', bbox_to_anchor=(1, 0.5), title = "Carbon Dioxid")
plt.xlim([-15,30])
plt.ylim([995,1015])
plt.xlabel("Latitude")
plt.ylabel("August surface pressure [hPa]")
plt.axvline(0, color = "black", linewidth = 0.7, linestyle = '--')
plt.axvline(10, color = "black", linewidth = 0.7)
#plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/6_Sensitivity_ppm/ps_aug.pdf", bbox_inches='tight')


# August rainfall
plt.figure()
plt.plot(lats,precip_z_140[8-1], label = "140 ppm", color = clrs[1])
plt.plot(lats,precip_z_280[8-1], label = "280 ppm", color = clrs[3])
plt.plot(lats,precip_z_560[8-1], label = "560 ppm", color = clrs[5])
plt.plot(lats,precip_z_1120[8-1], label = "1120 ppm", color = clrs[7])
plt.legend(frameon = False,loc='center left', bbox_to_anchor=(1, 0.5), title = "Carbon Dioxid")
plt.xlim([-15,30])
plt.ylim([0,14])
plt.xlabel("Latitude")
plt.ylabel("August rainfall [mm/day]")
plt.axvline(0, color = "black", linewidth = 0.7, linestyle = '--')
plt.axvline(10, color = "black", linewidth = 0.7)
#plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/6_Sensitivity_ppm/precip_aug.pdf", bbox_inches='tight')

# September surface pressure
plt.figure()
plt.plot(lats,ps_z_140[9-1], label = "140 ppm", color = clrs[1])
plt.plot(lats,ps_z_280[9-1], label = "280 ppm", color = clrs[3])
plt.plot(lats,ps_z_560[9-1], label = "560 ppm", color = clrs[5])
plt.plot(lats,ps_z_1120[9-1], label = "1120 ppm", color = clrs[7])
plt.legend(frameon = False,loc='center left', bbox_to_anchor=(1, 0.5), title = "Carbon Dioxid")
plt.xlim([-7,12.5])
plt.ylim([1001,1007])
plt.xlabel("Latitude")
plt.ylabel("September surface pressure [hPa]")
plt.axvline(0, color = "black", linewidth = 0.7, linestyle = '--')
plt.axvline(10, color = "black", linewidth = 0.7)
#plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/6_Sensitivity_ppm/ps_sep.pdf", bbox_inches='tight')


# September rainfall
plt.figure()
plt.plot(lats,precip_z_140[9-1], label = "140 ppm", color = clrs[1])
plt.plot(lats,precip_z_280[9-1], label = "280 ppm", color = clrs[3])
plt.plot(lats,precip_z_560[9-1], label = "560 ppm", color = clrs[5])
plt.plot(lats,precip_z_1120[9-1], label = "1120 ppm", color = clrs[7])
plt.legend(frameon = False,loc='center left', bbox_to_anchor=(1, 0.5), title = "Carbon Dioxid")
plt.xlim([-15,30])
plt.ylim([0,16])
plt.xlabel("Latitude")
plt.ylabel("September rainfall [mm/day]")
plt.axvline(0, color = "black", linewidth = 0.7, linestyle = '--')
plt.axvline(10, color = "black", linewidth = 0.7)
#plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/6_Sensitivity_ppm/precip_sep.pdf", bbox_inches='tight')


#%%
### Hadley cell and ITCZ


#Calculate ITCZ 
for i in range(0,len(ppm)):
    exec("ITCZ_" + str(ppm[i]) + " = []")
    for m in range(0,12):
        exec("ITCZ_" + str(ppm[i]) + ".append(lats[np.argmax(precip_z_" + str(ppm[i]) + "[m])])")#-evap_z_" + str(ppm[i]) + "[m]

# Show ITCZ
plt.figure()
plt.xlabel("Months")
plt.ylabel("Latitude")
plt.axhline(0,linewidth = 0.7, linestyle = '--', color = "black")
for i in range(0,len(ppm)):
    exec("plt.plot(range(1,13),ITCZ_" + str(ppm[i]) + ",label = " + str(ppm[i]) + " )")
plt.legend(frameon = False)
#plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/5_Sensitivity_ppm/ITCZ_stripes_ppm.pdf')


# calculate upper and lower end of Hadley cell 
for i in range(0,len(ppm)):
    exec("trans_nh_" + str(ppm[i]) + " = []") # northern hemisphere
    exec("lats_nh_" + str(ppm[i]) + " = []")
    exec("trans_sh_" + str(ppm[i]) + " = []") # southern hemisphere
    exec("lats_sh_" + str(ppm[i]) + " = []")
    for m in range(0,12):
        exec("olr = olr_z_" + str(ppm[i]) + "[m]")
        trans = []
        for s in range(0,len(olr)-1):
            if (olr[s] < 250 and olr[s+1] > 250) or (olr[s] > 250 and olr[s+1] < 250):
                trans.append(s)
                
        exec("trans_nh_" + str(ppm[i]) + ".append(trans[-1])")
        exec("trans_sh_" + str(ppm[i]) + ".append(trans[0])")
        exec("lats_sh_" + str(ppm[i]) + ".append(lats[trans[0]])")
        exec("lats_nh_" + str(ppm[i]) + ".append(lats[trans[-1]])")
    
    
# OLR (plot) and marked ends of Hadley cells
# i= 5
# print(ppm[i])
# for m in range(0,12):
#     plt.figure()
#     plt.title(month_list[m])
#     exec("plt.plot(lats,olr_z_" + str(ppm[i]) + "[m])")
#     exec("plt.axvline(lats_sh_" + str(ppm[i]) + "[m], color = 'orange')")
#     exec("plt.axvline(lats_nh_" + str(ppm[i]) + "[m], color = 'orange')")
       
    
# Plot ITCZ and Hadley cell extension 
plt.figure()
for i in range(0,len(ppm)):
    exec("plt.plot(range(1,13),lats_nh_" + str(ppm[i]) + ", color = clrs[i], linestyle = '--')")
    exec("plt.plot(range(1,13),lats_sh_" + str(ppm[i]) + ", color = clrs[i], linestyle = '--')")
    exec("plt.plot(range(1,13),ITCZ_" + str(ppm[i]) + ", color = clrs[i], label = str(ppm[i]) + 'ppm')")
plt.xlabel("Months")
plt.xlim([1,12])
plt.ylim([-70,70])
plt.ylabel("Latitude")
plt.axhline(0,color= "black", linewidth = 0.7,linestyle = '--')
plt.axhline(10,color= "black", linewidth = 0.7)
leg = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),frameon=False, title= "Carbon Dioxide")
leg.get_title().set_ha("left") 
#plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/6_Sensitivity_ppm/hadley_ppm.pdf', bbox_inches='tight')

    

#%%
### Monthly evolution of monsoon rain

# The grid cell weighted average was calculated via processing.py

# Monthly monsoon rainfall for each of the simulations
precip_m=[[6.4199335e-09,5.0032434e-10,3.5851697e-09,9.7042090e-08,1.4742824e-06
,7.8892626e-06,1.8152374e-05,3.4551183e-05,5.7049674e-05,2.1263673e-05
,5.3469029e-07,4.5584827e-08]
,[1.0564064e-08,5.0026067e-10,3.0060043e-09,5.0247277e-08,9.7841803e-07
,7.4703080e-06,2.0072952e-05,3.8057409e-05,6.5999266e-05,3.0377874e-05
,8.0059073e-07,5.0394330e-08]
,[2.6895783e-09,6.8969008e-10,2.9988285e-09,6.0667602e-08,1.1285196e-06
,8.8589541e-06,2.0973081e-05,3.9816350e-05,6.9490539e-05,3.5824036e-05
,1.0256909e-06,4.0709573e-08]
,[1.3660164e-08,4.7876936e-10,3.1295400e-09,7.9958681e-08,1.1295950e-06
,9.1309057e-06,2.1973536e-05,4.3431541e-05,7.4565265e-05,4.2340656e-05
,1.5969730e-06,2.9456896e-08]
,[7.9104465e-09,4.5139190e-10,3.3680558e-09,6.1642766e-08,1.1622068e-06
,9.7975662e-06,2.4622628e-05,4.9278366e-05,7.8924393e-05,4.8464084e-05
,1.9752661e-06,5.2929419e-08]
,[1.1568752e-08,5.1290944e-10,3.9284505e-09,4.5581473e-08,9.7300699e-07
,1.0245581e-05,2.6694232e-05,5.3077132e-05,8.5328880e-05,6.0334900e-05
,3.0231095e-06,9.3089895e-08]
,[7.6370110e-09,5.4870758e-10,3.4822134e-09,5.4210055e-08,7.8976137e-07
,1.2162588e-05,3.1422667e-05,5.7509398e-05,9.1322516e-05,6.7451430e-05
,4.7402718e-06,1.5845217e-07]
,[2.9153679e-08,5.6962374e-10,4.5894066e-09,5.9455306e-08,1.0426485e-06
,1.4793726e-05,3.4876499e-05,6.4090818e-05,9.3765208e-05,7.6016106e-05
,7.4993827e-06,2.1835483e-07]]


plt.figure()
for i in range(0,len(ppm)):
    plt.plot(range(1,13),np.array(precip_m[i])*86400, label =  str(ppm[i]) + 'ppm', color = clrs[i])
plt.ylabel("Monsoon rainfall [mm/day]")
plt.xlabel('Months')
plt.ylim([0,9])
plt.axhline(0,color="black",linewidth = 0.7)
plt.legend(frameon = False, bbox_to_anchor = (1,0.5), loc = "center left", title = "Carbon Dioxide")
#plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/6_Sensitivity_ppm/pr_per_month_ppm.pdf", bbox_inches='tight')



#%%
### In dependence of CO2 Concentration

# mr = monsoon region (10-30°N)
# JJAS = June-September

# The grid cell weighted average was calculated via processing.py

# Monsoon rainfall, evaporation and transportable water
precip_mr_annual=[1.1755689e-05,1.3655961e-05,1.4768744e-05,1.6191263e-05,1.7862567e-05,1.9985961e-05,2.2135247e-05,2.4366376e-05]
evap_global_annual=[2.7610276e-05,2.8834087e-05,2.9465933e-05,3.0006197e-05,3.0601092e-05,3.1140251e-05,3.1679629e-05,3.2310218e-05]
wvp_global_annual=[17.396332,20.108015,21.8638,23.44319,25.454134,27.591684,29.629019,32.534977]

fig,ax = plt.subplots()
plt.ylim([200,1200])
pr, = ax.plot(ppm,np.array(precip_mr_annual[:13])*86400*365, color = "darkblue", label= "Monsoon rainfall (10-30°N)", marker = ".")
evap, = plt.plot(ppm,np.array(evap_global_annual[:13])*86400*365, color = "darkgreen",linestyle = '--',label= "Evaporation (global)", marker = ".")
plt.ylabel("Water [mm]")
plt.xlabel("$CO_2$ concentration [ppm]")
ax2=ax.twinx()
wvp, = plt.plot(ppm,np.array(wvp_global_annual[:13]), color = "darkorange",linestyle = '--',label= "Transportable water (global)", marker = ".")
plt.ylim([15,35])
ax2.set_ylabel("Transportable Water [mm]", color = "darkorange")
ax2.tick_params(axis='y', colors='darkorange')
plt.legend(frameon = False)
plt.legend([pr,evap,wvp], ['Monsoon rainfall (10-30°N)', 'Evaporation (global)','Transportable water (global)'], numpoints=1,handler_map={tuple: HandlerTuple(ndivide=None)}, handlelength=3, frameon = False, loc = "lower right")
#plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/6_Sensitivity_ppm/evap_contribution2.pdf', bbox_inches='tight')


# Overview figure
precip_global_annual=[2.7611337e-05,2.8835879e-05,2.9467377e-05,3.0008980e-05,3.0602685e-05,3.1142197e-05,3.1681437e-05,3.2312710e-05]
precip_mr_annual=[1.1755689e-05,1.3655961e-05,1.4768744e-05,1.6191263e-05,1.7862567e-05,1.9985961e-05,2.2135247e-05,2.4366376e-05]
tsurf_global_annual=[286.33853,288.58524,289.7542,290.76453,291.89136,293.01608,293.93185,295.249]
tsurf_mr_annual=[294.97186,297.93857,299.5395,300.95865,302.4874,304.06342,305.2409,307.09006]

fig,ax = plt.subplots()
pr_mr, = plt.plot(ppm,np.array(precip_mr_annual)*86400*365, color = "darkblue", marker = ".")
pr_g, = plt.plot(ppm,np.array(precip_global_annual)*86400*365, color = "darkblue",linestyle = '--', marker = ".")
plt.xlabel("Carbon dioxide [ppm]")
plt.ylabel("Annual precipitation [mm]", color = "darkblue")
ax.tick_params(axis='y', colors='darkblue')
plt.legend(frameon = False)
ax2=ax.twinx()
tsurf_mr, = ax2.plot(ppm,np.array(tsurf_mr_annual)-273.15, color = 'darkorange', label = "Monsoon region", marker = ".")
tsurf_g, = ax2.plot(ppm,np.array(tsurf_global_annual)-273.15, color = 'darkorange',linestyle = '--', label = "Global", marker = ".")
ax2.set_ylabel("Surface Temperature [°C]", color = "darkorange")
ax2.tick_params(axis='y', colors='darkorange')
plt.legend([(pr_g, tsurf_g), (pr_mr, tsurf_mr)], ['Global', 'Monsoon region'], numpoints=1,handler_map={tuple: HandlerTuple(ndivide=None)}, handlelength=3, frameon = False, loc = "lower right")
#plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/6_Sensitivity_ppm/pr_per_ppm.pdf", bbox_inches='tight')















