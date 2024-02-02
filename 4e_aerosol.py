#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#-------------------------------------------------------------
### ANALYSIS OF AEROSOL SIMULATIONS ###
# Anja Katzeberger, anja.katzenberger@pik-potsdam.de
#-------------------------------------------------------------

# This script creates the figures as displayed in the subchapter "Sensitivty to aerosol changes" and some more figures used in the context of the analysis
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

aerosol = [0, 10**(-7), 10**(-6), 2*10**(-6), 4*10**(-6), 6*10**(-6), 8*10**(-6), 10**(-5), 2*10**(-5), 10**(-4)]
aerosol_name = [0,7, 6, '6_2', '6_4', '6_6', '6_8', 5, '5_2',4] # 6_2 refers to 2*10**(-6) etc.
aerosol_name2 = ['0','$10^{-7}$','$10^{-6}$','$2x10^{-6}$','$4x10^{-6}$','$6x10^{-6}$','$8x10^{-6}$','$10^{-5}$','$2x10^{-5}$','$10^{-4}$']

dir = []
for i in range(0,len(aerosol_name)):
    print(i)
    dir.append("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/A_data/MonsoonPlanet_aerosol_" + str(aerosol_name[i]) + ".nc")
    print(dir[i])
    exec("ds_" + str(aerosol_name[i]) + " = nc.Dataset(dir[i])")
    
    exec("evap_" + str(aerosol_name[i]) + " = ds_" + str(aerosol_name[i]) + "['evap'][:]")
    exec("precip_" + str(aerosol_name[i]) + " = ds_" + str(aerosol_name[i]) + "['precip'][:]")
    exec("tsurf_" + str(aerosol_name[i]) + " = ds_" + str(aerosol_name[i]) + "['t_surf'][:]")
    exec("tref_" + str(aerosol_name[i]) + " = ds_" + str(aerosol_name[i]) + "['t_ref'][:]")
    exec("ps_" + str(aerosol_name[i]) + " = ds_" + str(aerosol_name[i]) + "['ps'][:]")
    exec("olr_" + str(aerosol_name[i]) + " = ds_" + str(aerosol_name[i]) + "['olr'][:]")
    exec("ucomp_" + str(aerosol_name[i]) + " = ds_" + str(aerosol_name[i]) + "['ucomp'][:]")
    exec("vcomp_" + str(aerosol_name[i]) + " = ds_" + str(aerosol_name[i]) + "['vcomp'][:]")
    exec("sw_" + str(aerosol_name[i]) + " = ds_" + str(aerosol_name[i]) + "['swdn_toa'][:]")
    exec("swsur_" + str(aerosol_name[i]) + " = ds_" + str(aerosol_name[i]) + "['swdn_sfc'][:]") # shortwave down top of atmosphere
    exec("cld_" + str(aerosol_name[i]) + " = ds_" + str(aerosol_name[i]) + "['cld_amt'][:]") # cloud amount
    exec("cldt_" + str(aerosol_name[i]) + " = ds_" + str(aerosol_name[i]) + "['tot_cld_amt'][:]") # total cloud amount
    exec("rh_" + str(aerosol_name[i]) + " = ds_" + str(aerosol_name[i]) + "['rh'][:]") # relative humidity
    exec("sphum_" + str(aerosol_name[i]) + " = ds_" + str(aerosol_name[i]) + "['sphum'][:]") # specific humidity
    exec("wvp_" + str(aerosol_name[i]) + " = ds_" + str(aerosol_name[i]) + "['WVP'][:]") # column integrated water vapour



    exec("evap_z_" + str(aerosol_name[i]) + " = evap_" + str(aerosol_name[i]) + ".mean(axis=2)*86400")
    exec("precip_z_" + str(aerosol_name[i]) + " = precip_" + str(aerosol_name[i]) + ".mean(axis = 2)*86400")
    exec("olr_z_" + str(aerosol_name[i]) + " = olr_" + str(aerosol_name[i]) + ".mean(axis=2)")
    exec("vcomp_z_" + str(aerosol_name[i]) + " = vcomp_" + str(aerosol_name[i]) + ".mean(axis=3)")
    exec("ucomp_z_" + str(aerosol_name[i]) + " = ucomp_" + str(aerosol_name[i]) + ".mean(axis=3)")
    exec("ps_z_" + str(aerosol_name[i]) + " = ps_" + str(aerosol_name[i]) + ".mean(axis=2)/100")
    exec("tsurf_z_" + str(aerosol_name[i]) + " = tsurf_" + str(aerosol_name[i]) + ".mean(axis=2)-273.15")
    exec("tref_z_" + str(aerosol_name[i]) + " = tref_" + str(aerosol_name[i]) + ".mean(axis=2)-273.15")
    exec("sw_z_" + str(aerosol_name[i]) + " = sw_" + str(aerosol_name[i]) + ".mean(axis=2)")
    exec("swsur_z_" + str(aerosol_name[i]) + " = swsur_" + str(aerosol_name[i]) + ".mean(axis=2)")
    exec("cld_z_" + str(aerosol_name[i]) + " = cld_" + str(aerosol_name[i]) + ".mean(axis=3)")
    exec("cldt_z_" + str(aerosol_name[i]) + " = cldt_" + str(aerosol_name[i]) + ".mean(axis=2)")
    exec("rh_z_" + str(aerosol_name[i]) + " = rh_" + str(aerosol_name[i]) + ".mean(axis=3)")
    exec("sphum_z_" + str(aerosol_name[i]) + " = sphum_" + str(aerosol_name[i]) + ".mean(axis=3)")
    exec("wvp_z_" + str(aerosol_name[i]) + " = wvp_" + str(aerosol_name[i]) + ".mean(axis=2)")
                 

lats = ds_0.variables['lat'][:]
lons = ds_0.variables['lon'][:]
pfull = ds_0.variables['pfull'][:]


month_list = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
clrs = ['#6EA6CD','#98CAE1','#C2E4EF','#ACD39E','#5AAE61','#FDB366','#F67E4B','#DD3D2D','#A50026','#762A83','#36489A','#4A7BB7']



#%%
### General plots

# Meridional precipitation distribution (contourf)
for i in range(0,len(aerosol)):
    print(aerosol_name[i])
    plt.figure()
    plt.title(str(aerosol_name2[i])+' $kg/m^2$')
    exec("plt.contourf(range(1,13),lats,precip_z_" + str(aerosol_name[i]) + ".transpose(),cmap='Blues',extend='max',levels = [0,2,4,6,8,10,12,14,16,18,20])")
    plt.axhline(0, linewidth = 0.7, color = 'black',linestyle = '--')
    plt.axhline(10, linewidth = 0.7, color = 'black')
    plt.axhline(60, linewidth = 0.7, color = 'black')
    plt.xlabel('Months')
    plt.xlim([1,12])
    plt.ylabel('Latitude')
    plt.ylim([-35,35])
    clb = plt.colorbar()
    clb.ax.set_title('mm/day',fontsize=8)
    #plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/7_Sensitivity_aerosol/precip_10N_280ppm_aerosol" + str(aerosol_name[i]) + "_contour.pdf", bbox_inches='tight')


# Meridional tsurf distribution (contourf)
for i in range(0,len(aerosol)):
    plt.figure()
    plt.title(str(aerosol_name2[i])+' $kg/m^2$')
    exec("plt.contourf(range(1,13),lats,tsurf_z_" + str(aerosol_name[i]) + ".transpose(),cmap='OrRd',extend='both',levels = range(0,50))")
    plt.axhline(0, linewidth = 0.7, color = 'black',linestyle = '--')
    plt.axhline(10, linewidth = 0.7, color = 'black')
    plt.axhline(60, linewidth = 0.7, color = 'black')
    plt.xlabel('Months')
    plt.xlim([1,12])
    plt.ylabel('Latitude')
    plt.ylim([-35,35])
    clb = plt.colorbar()
    clb.ax.set_title('°C',fontsize=8)
    plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/00_Publication_figures/7_Sensitivity_aerosol/tsurf_10N_280ppm_aerosol" + str(aerosol_name[i]) + "contour.pdf", bbox_inches='tight')



# Meridional vcomp distribution (plot)
for i in range(0,len(aerosol)):
    fig = plt.figure()
    for m in range(0,12):
        exec("vcompp = vcomp_z_" + str(aerosol_name[i]) + "[m]")
        plt.title(str(aerosol_name2[i])+' $kg/m^2$')
        plt.xlabel('latitude')
        plt.xlim([-91,91])
        plt.ylim([-6,6])
        plt.xticks(np.arange(-90, 91, step=30))
        plt.ylabel( 'Meridional surfce winds (m/s)')
        plt.axvline(10, linewidth = 0.7, color = 'black')
        plt.axvline(60, linewidth = 0.7, color = 'black')
        plt.axhline(0, linewidth = 0.7, color = 'black')
        exec("plt.plot(lats,vcompp[23],label=month_list[m], color = clrs[m], linewidth = 1)") #alpha = 1 
    plt.legend(loc= 'center left',frameon = False, bbox_to_anchor = (1,0.5))
    #plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/7_Sensitivity_aerosol/vcomp_10N_280ppm_aerosol" + str(aerosol_name[i]) + "_contour.pdf", bbox_inches='tight')




#%%
### Pressure barrier dynamics

# September surface pressure zoom
plt.figure()
plt.plot(lats,ps_z_4[9-1], label = str(aerosol_name2[9])+' $kg/m^2$', color = clrs[1])
plt.plot(lats,ps_z_5[9-1], label = str(aerosol_name2[7])+' $kg/m^2$', color = clrs[3])
plt.plot(lats,ps_z_6[9-1], label = str(aerosol_name2[2])+' $kg/m^2$', color = clrs[5])
plt.plot(lats,ps_z_7[9-1], label = str(aerosol_name2[1])+' $kg/m^2$', color = clrs[7])
plt.legend(frameon = False,loc='center left', bbox_to_anchor=(1, 0.5), title = "Aerosols")
plt.xlim([-7,12.5])
plt.ylim([996,1006])
plt.xlabel("Latitude")
plt.ylabel("September surface pressure [hPa]")
plt.axvline(0, color = "black", linewidth = 0.7, linestyle = '--')
plt.axvline(10, color = "black", linewidth = 0.7)
#plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/7_Sensitivity_aerosol/ps_sep_zoom.pdf", bbox_inches='tight')

# September precipitation
plt.figure()
plt.plot(lats,precip_z_4[9-1], label = "$ 10^{-4} kg/m^2$", color = clrs[1])
plt.plot(lats,precip_z_5[9-1], label = "$ 10^{-5} kg/m^2$", color = clrs[3])
plt.plot(lats,precip_z_6[9-1], label = "$ 10^{-6} kg/m^2$", color = clrs[5])
plt.plot(lats,precip_z_7[9-1], label = "$10^{-7} kg/m^2$", color = clrs[7])
plt.legend(frameon = False,loc='center left', bbox_to_anchor=(1, 0.5), title = "Aerosol concentration")
plt.xlim([-15,30])
plt.ylim([0,16])
plt.xlabel("Latitude")
plt.ylabel("September precipitation [mm/day]")
plt.axvline(0, color = "black", linewidth = 0.7, linestyle = '--')
plt.axvline(10, color = "black", linewidth = 0.7)
#plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/7_Sensitivity_aerosol/precip_sep_zoom.pdf", bbox_inches='tight')



#%%
### Hadley cell and ITCZ

#Calculate ITCZ 
for i in range(0,len(aerosol)):
    exec("ITCZ_" + str(aerosol_name[i]) + " = []")
    for m in range(0,12):
        exec("ITCZ_" + str(aerosol_name[i]) + ".append(lats[np.argmax(precip_z_" + str(aerosol_name[i]) + "[m])])")#-evap_z_" + str(sol[i]) + "[m]

# Show ITCZ
plt.figure()
plt.xlabel("Months")
plt.ylabel("Latitude")
plt.axhline(0,linewidth = 0.7, linestyle = '--', color = "black")
for i in range(0,len(aerosol)):
    exec("plt.plot(range(1,13),ITCZ_" + str(aerosol_name[i]) + ",color = clrs[i], label = " + str(aerosol_name[i]) + " )")
plt.legend(frameon = False, bbox_to_anchor = (1,0.5), loc = "center left", title = "Aerosols")
#plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/6_Sensitivity_aerosol/ITCZ_stripes_aerosol.pdf', bbox_inches = 'tight' )


# calculate upper and lower end of Hadley cell 
# 10^-4 aerosol run does not fulfill definition of Hadley cell edges (does not reach 250 W/m^2)
for i in range(0,len(aerosol)-1):
    print(aerosol_name[i])
    exec("trans_nh_" + str(aerosol_name[i]) + " = []") # northern hemisphere
    exec("lats_nh_" + str(aerosol_name[i]) + " = []")
    exec("trans_sh_" + str(aerosol_name[i]) + " = []") # southern hemisphere
    exec("lats_sh_" + str(aerosol_name[i]) + " = []")
    for m in range(0,12):
        exec("olr = olr_z_" + str(aerosol_name[i]) + "[m]")
        trans = []
        for s in range(0,len(olr)-1):
            if (olr[s] < 250 and olr[s+1] > 250) or (olr[s] > 250 and olr[s+1] < 250):
                trans.append(s)
        exec("trans_nh_" + str(aerosol_name[i]) + ".append(trans[-1])")
        exec("trans_sh_" + str(aerosol_name[i]) + ".append(trans[0])")
        exec("lats_sh_" + str(aerosol_name[i]) + ".append(lats[trans[0]])")
        exec("lats_nh_" + str(aerosol_name[i]) + ".append(lats[trans[-1]])")
    

# OLR (plot) and marked ends of Hadley cells
# i= 5
# print(aerosol[i])
# for m in range(0,12):
#     plt.figure()
#     plt.title(month_list[m])
#     exec("plt.plot(lats,olr_z_" + str(aerosol_name[i]) + "[m])")
#     exec("plt.axvline(lats_sh_" + str(aerosol_name[i]) + "[m], color = 'orange')")
#     exec("plt.axvline(lats_nh_" + str(aerosol_name[i]) + "[m], color = 'orange')")
       
    
# Plot ITCZ and Hadley cell extension 
plt.figure()
for i in range(0,len(aerosol)-1):
    exec("plt.plot(range(1,13),lats_nh_" + str(aerosol_name[i]) + ", color = clrs[i], linestyle = '--')")
    exec("plt.plot(range(1,13),lats_sh_" + str(aerosol_name[i]) + ", color = clrs[i], linestyle = '--')")
    exec("plt.plot(range(1,13),ITCZ_" + str(aerosol_name[i]) + ", color = clrs[i], label = str(aerosol_name2[i]) + '$ kg/m^2$')")
#plt.plot(range(1,13),ITCZ_89, color = clrs[i], label = 'Aquaplanet')
plt.xlabel("Months")
plt.xlim([1,12])
plt.ylim([-70,70])
plt.ylabel("Latitude")
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),frameon=False, title = "Sulfate Aerosols")
#plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/7_Sensitivity_aerosol/hadley_aerosol.pdf', bbox_inches='tight')



#%%
### Monthly evolution of monsoon rain

# The grid cell weighted average was calculated via processing.py

# Monthly monsoon rainfall for each of the simulations
precip_m=[[1.36601637e-08,4.78769357e-10,3.12953996e-09,7.99586815e-08
,1.12959503e-06,9.13090571e-06,2.19735357e-05,4.34315407e-05
,7.45652651e-05,4.23406564e-05,1.59697299e-06,2.94568956e-08]
,[3.41349038e-09,4.88320606e-10,3.06426373e-09,5.57806814e-08
,1.09780831e-06,8.55950384e-06,2.10945836e-05,4.09420245e-05
,7.27461666e-05,4.06755025e-05,1.21922574e-06,2.14375877e-08]
,[8.25431190e-09,4.60900762e-10,2.75245449e-09,4.28608651e-08
,8.32389162e-07,7.27495535e-06,1.90321371e-05,3.70038506e-05
,6.32174560e-05,3.25809815e-05,1.02211709e-06,3.85913310e-08]
,[4.84194596e-09,5.80272108e-10,2.36333686e-09,2.53796664e-08
,5.81516076e-07,5.56431769e-06,1.59834708e-05,2.98907453e-05
,5.20882531e-05,2.52850696e-05,6.53404868e-07,1.12016858e-08]
,[3.40526651e-09,3.78173937e-10,1.79287574e-09,1.59668883e-08
,2.62306116e-07,3.75578065e-06,1.27234625e-05,2.29630987e-05
,3.79100493e-05,1.40556558e-05,2.20580787e-07,2.89916251e-08]
,[7.29319438e-10,2.10505458e-09,1.51413326e-09,1.05827231e-08
,1.28207176e-07,2.09403765e-06,9.66691823e-06,1.75429486e-05
,2.77222025e-05,8.60444652e-06,9.13821623e-08,2.86204109e-08]
,[3.77332210e-09,2.70849371e-10,1.28898403e-09,7.63268559e-09
,7.62502737e-08,1.13926399e-06,6.90089428e-06,1.34173952e-05
,2.16975986e-05,5.47503623e-06,2.57785917e-08,2.66311293e-08]
,[1.53767177e-09,2.37161818e-10,1.02778286e-09,6.54236088e-09
,4.74106443e-08,7.97876453e-07,4.78161610e-06,9.72806993e-06
,1.78119553e-05,3.64891730e-06,1.82732123e-08,2.99932132e-08]
,[2.99527603e-09,8.11816170e-11,1.20898208e-10,4.62899746e-10
,1.24500179e-08,1.15194496e-07,9.07293952e-07,4.05217361e-06
,7.20640219e-06,2.50155722e-07,2.16809593e-08,9.21277632e-09]
,[7.98288102e-10,2.05461759e-11,1.88726899e-11,4.96812487e-12
,1.66143228e-12,1.65673794e-11,4.94541630e-11,1.13923718e-10
,1.82689974e-10,3.37139865e-11,6.72382462e-12,3.92896566e-12]]

plt.figure()
for i in range(0,len(aerosol)):
    plt.plot(range(1,13),np.array(precip_m[i])*86400, label =  str(aerosol[i]) + '$kg/m^2$', color = clrs[i])
plt.ylabel("Monsoon rainfall [mm/day]")
plt.xlabel('Months')
plt.ylim([0,7])
plt.axhline(0,color="black",linewidth = 0.7)
plt.legend(frameon = False, bbox_to_anchor = (1,0.5),loc = "center left",title = "Sulfate Aerosols")
#plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/7_Sensitivity_aerosol/pr_per_month_aerosol.pdf", bbox_inches='tight')
  


#%%
### In dependence of aerosols

# mr = monsoon region (10-30°N)
# JJAS = June-September

# The grid cell weighted average was calculated via processing.py
precip_mr_annual = [1.6191263e-05,1.5534917e-05,1.3421401e-05,1.0840929e-05,7.6617889e-06,5.4911411e-06,4.0643176e-06,3.0727881e-06,1.0481854e-06,1.0427821e-10]
precip_global_annual = [3.0008980e-05,2.9977744e-05,2.9738982e-05,2.9374294e-05,2.8781433e-05,2.8154403e-05,2.7674276e-05,2.7267221e-05,2.5944059e-05,2.1938815e-05]
tsurf_land_annual = [293.54468,293.56046,292.9768,292.24304,291.0919,289.76193,288.48123,287.1858,281.3394]
tsurf_10S_10N_annual = [300.86884,300.84454,300.7271,300.39047,300.01437,299.46335,298.81555,298.2871,295.77118]
tsurf_30S_10N_annual = [298.81635,298.83368,298.77325,298.552,298.3516,298.0355,297.6334,297.31152,295.70026]
tsurf_90S_10N_annual = [291.79892,291.7631,291.69363,291.57126,291.37946,291.14612,290.89545,290.68445,289.7397]
tsurf_mr_jun_oct=[309.73355,309.60574,309.52783,309.1793,308.7015,308.0566,307.085,306.14166,300.4314,269.95694]
tsurf_10S_10N_jun_oct=[301.70868,301.6804,301.56512,301.22806,300.85092,300.32092,299.69696,299.19037,296.7388,283.99985]
tsurf_global_annual=[290.76453,290.75244,290.4779,290.124,289.56592,288.9163,288.28098,287.64954,284.87186,270.27216]
tsurf_mr_annual=[300.95865,300.83334,300.54657,299.90344,299.00665,297.88428,296.5743,295.4022,289.74847,264.77393]



# Decrease in temperature with increasing sulfate aerosol concentration
plt.figure()
plt.plot(aerosol,np.array(tsurf_mr_jun_oct)-273.15, label = "Monsoon region (10-30°N)", color = "darkorange")
plt.plot(aerosol,np.array(tsurf_10S_10N_jun_oct)-273.15, label = "10°S-10°N", color = "darkblue")
plt.legend(frameon = False)
plt.ylabel("Surface temperatures [°C]")
plt.xlabel("Sulfate Aerosol [$kg/m^2$]")
#plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/7_Sensitivity_aerosol/tsurf_gradient.pdf", bbox_inches='tight')

# Monsoon rainfall and global precipitation in dependence of sulfate aerosols
plt.figure()
plt.plot(aerosol,np.array(precip_mr_annual)*86400*365, label = "monsoon region")
plt.plot(aerosol,np.array(precip_global_annual)*86400*365, label = "global")
plt.ylabel("Annual monsoon rainfall [mm]")
plt.xlabel("Sulfate Aerosol [$kg/m^2$]")
plt.legend(frameon = False)

# Overview figure
fig,ax = plt.subplots()
pr_mr, = plt.plot(aerosol,np.array(precip_mr_annual)*86400*365, color = "darkblue", marker = ".")
pr_g, = plt.plot(aerosol,np.array(precip_global_annual)*86400*365, color = "darkblue",linestyle = '--', marker = ".")
plt.xlabel("Aerosol concentration [$kg/m^2$]")
plt.ylabel("Annual precipitation [mm]", color = "darkblue")
ax.tick_params(axis='y', colors='darkblue')
plt.legend(frameon = False)
ax2=ax.twinx()
tsurf_mr, = ax2.plot(aerosol,np.array(tsurf_mr_annual)-273.15, color = 'darkorange', label = "Monsoon region", marker = ".")
tsurf_g, = ax2.plot(aerosol,np.array(tsurf_global_annual)-273.15, color = 'darkorange',linestyle = '--', label = "Global", marker = ".")
ax2.set_ylabel("Surface Temperature [°C]", color = "darkorange")
ax2.tick_params(axis='y', colors='darkorange')
plt.legend([(pr_g, tsurf_g), (pr_mr, tsurf_mr)], ['Global', 'Monsoon region'], numpoints=1,handler_map={tuple: HandlerTuple(ndivide=None)}, handlelength=3, frameon = False)
#plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/7_Sensitivity_aerosol/pr_per_aerosol.pdf", bbox_inches='tight')






