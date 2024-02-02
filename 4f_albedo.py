#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#-------------------------------------------------------------
### ANALYSIS OF ALBEDO SIMULATIONS ###
# Anja Katzeberger, anja.katzenberger@pik-potsdam.de
#-------------------------------------------------------------

# This script creates the figures as displayed in the subchapter "Sensitivty to albedo changes" and some more figures used in the context of the analysis
# The directories of the input files as well as directories for saving the figures must be adapted. 

import math
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import statistics as stat


#%%
### Loading data 

albedo = [10,12,14,16,18,20,22,24,26,28,30,32,34]

dir = []
for i in range(0,len(albedo)):
    print(i)
    dir.append("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/A_data/MonsoonPlanet_albedo_" + str(albedo[i]) + ".nc")
    print(dir[i])
    exec("ds_" + str(albedo[i]) + " = nc.Dataset(dir[i])")
    
    exec("evap_" + str(albedo[i]) + " = ds_" + str(albedo[i]) + "['evap'][:]")
    exec("precip_" + str(albedo[i]) + " = ds_" + str(albedo[i]) + "['precip'][:]")
    exec("tsurf_" + str(albedo[i]) + " = ds_" + str(albedo[i]) + "['t_surf'][:]")
    exec("tref_" + str(albedo[i]) + " = ds_" + str(albedo[i]) + "['t_ref'][:]")
    exec("ps_" + str(albedo[i]) + " = ds_" + str(albedo[i]) + "['ps'][:]")
    exec("olr_" + str(albedo[i]) + " = ds_" + str(albedo[i]) + "['olr'][:]")
    exec("ucomp_" + str(albedo[i]) + " = ds_" + str(albedo[i]) + "['ucomp'][:]")
    exec("vcomp_" + str(albedo[i]) + " = ds_" + str(albedo[i]) + "['vcomp'][:]")
    exec("sw_" + str(albedo[i]) + " = ds_" + str(albedo[i]) + "['swdn_toa'][:]")
    exec("swsur_" + str(albedo[i]) + " = ds_" + str(albedo[i]) + "['swdn_sfc'][:]") # shortwave down top of atmosphere
    exec("cld_" + str(albedo[i]) + " = ds_" + str(albedo[i]) + "['cld_amt'][:]") # cloud amount
    exec("cldt_" + str(albedo[i]) + " = ds_" + str(albedo[i]) + "['tot_cld_amt'][:]") # total cloud amount
    exec("rh_" + str(albedo[i]) + " = ds_" + str(albedo[i]) + "['rh'][:]") # relative humidity
    exec("sphum_" + str(albedo[i]) + " = ds_" + str(albedo[i]) + "['sphum'][:]") # specific humidity
    exec("wvp_" + str(albedo[i]) + " = ds_" + str(albedo[i]) + "['WVP'][:]") # column integrated water vapour


    exec("evap_z_" + str(albedo[i]) + " = evap_" + str(albedo[i]) + ".mean(axis=2)*86400")
    exec("precip_z_" + str(albedo[i]) + " = precip_" + str(albedo[i]) + ".mean(axis = 2)*86400")
    exec("olr_z_" + str(albedo[i]) + " = olr_" + str(albedo[i]) + ".mean(axis=2)")
    exec("vcomp_z_" + str(albedo[i]) + " = vcomp_" + str(albedo[i]) + ".mean(axis=3)")
    exec("ucomp_z_" + str(albedo[i]) + " = ucomp_" + str(albedo[i]) + ".mean(axis=3)")
    exec("ps_z_" + str(albedo[i]) + " = ps_" + str(albedo[i]) + ".mean(axis=2)/100")
    exec("tsurf_z_" + str(albedo[i]) + " = tsurf_" + str(albedo[i]) + ".mean(axis=2)-273.15")
    exec("tref_z_" + str(albedo[i]) + " = tref_" + str(albedo[i]) + ".mean(axis=2)-273.15")
    exec("sw_z_" + str(albedo[i]) + " = sw_" + str(albedo[i]) + ".mean(axis=2)")
    exec("swsur_z_" + str(albedo[i]) + " = swsur_" + str(albedo[i]) + ".mean(axis=2)")
    exec("cld_z_" + str(albedo[i]) + " = cld_" + str(albedo[i]) + ".mean(axis=3)")
    exec("cldt_z_" + str(albedo[i]) + " = cldt_" + str(albedo[i]) + ".mean(axis=2)")
    exec("rh_z_" + str(albedo[i]) + " = rh_" + str(albedo[i]) + ".mean(axis=3)")
    exec("sphum_z_" + str(albedo[i]) + " = sphum_" + str(albedo[i]) + ".mean(axis=3)")
    exec("wvp_z_" + str(albedo[i]) + " = wvp_" + str(albedo[i]) + ".mean(axis=2)")
                 

lats = ds_10.variables['lat'][:]
lons = ds_10.variables['lon'][:]
pfull = ds_10.variables['pfull'][:]


month_list = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
clrs = ['#6EA6CD','#98CAE1','#C2E4EF','#ACD39E','#5AAE61','#FDB366','#F67E4B','#DD3D2D','#A50026','#762A83','#36489A','#4A7BB7']




#%%
### General plots 

# Meridional rainfall distribution
for i in range(0,len(albedo)):
    fig = plt.figure()
    for m in range(0,12):
        plt.title(str(albedo[i])+' %')
        plt.xlabel('latitude')
        plt.xlim([-91,91])
        plt.ylim([-5,35])
        plt.xticks(np.arange(-90, 91, step=30))
        plt.ylabel( 'Precipitation (mm/day)')
        plt.axvline(10, linewidth = 0.7, color = 'black')
        plt.axvline(60, linewidth = 0.7, color = 'black')
        plt.axhline(0, linewidth = 0.7, color = 'black')
        exec("plt.plot(lats,precip_z_" + str(albedo[i]) + "[m],label=month_list[m], color = clrs[m], linewidth = 1)") #alpha = 1 
    plt.legend(loc= 'center left', bbox_to_anchor = (1,0.5), frameon = False)
    #plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/8_Sensitivity_albedo/precip_' + str(albedo[i]) + 'albedo.pdf' )
  
    
# Meridional precipitation distribution - Contour
for i in range(0,len(albedo)):
    plt.figure()
    plt.title(str(albedo[i])+' %')
    exec("plt.contourf(range(1,13),lats,precip_z_" + str(albedo[i]) + ".transpose(),cmap='Blues',extend='max',levels = [0,2,4,6,8,10,12,14,16,18,20])")
    plt.axhline(0, linewidth = 0.7, color = 'black',linestyle = '--')
    plt.axhline(10, linewidth = 0.7, color = 'black')
    plt.axhline(60, linewidth = 0.7, color = 'black')
    plt.xlabel('Months')
    plt.xlim([1,12])
    plt.ylabel('Latitude')
    plt.ylim([-35,35])
    clb = plt.colorbar()
    clb.ax.set_title('mm/day',fontsize=8)
    #plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/8_Sensitivity_albedo/precip_10N_280ppm_" + str(albedo[i]) + "albedo_contour.pdf", bbox_inches='tight')


# Meridional tsurf distribution - Contour
for i in range(0,len(albedo)):
    plt.figure()
    plt.title(str(albedo[i])+' %')
    exec("plt.contourf(range(1,13),lats,tsurf_z_" + str(albedo[i]) + ".transpose(),cmap='OrRd',extend='both',levels = range(0,50))")
    plt.axhline(0, linewidth = 0.7, color = 'black',linestyle = '--')
    plt.axhline(10, linewidth = 0.7, color = 'black')
    plt.axhline(60, linewidth = 0.7, color = 'black')
    plt.xlabel('Months')
    plt.xlim([1,12])
    plt.ylabel('Latitude')
    plt.ylim([-35,35])
    clb = plt.colorbar()
    clb.ax.set_title('°C',fontsize=8)
    plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/00_Publication_figures/8_Sensitivity_albedo/tsurf_10N_280ppm_" + str(albedo[i]) + "albedo_contour.pdf", bbox_inches='tight')


# Surface temperature (plot) + Surface pressure (plot)
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
    ax.plot(lats,tsurf_z_30[m-1], color = "darkorange")
    ax.plot(lats,tsurf_z_10[m-1], color = "darkorange",linestyle = '--')
    plt.axhline(0,color = "black", linewidth = 0.7)
    ax2=ax.twinx()
    ax2.plot(lats,ps_z_30[m-1], color = "darkgreen")
    ax2.plot(lats,ps_z_10[m-1], color = "darkgreen",linestyle = '--')
    ax2.set_ylim([995,1015])
    ax2.axhline(0,linewidth = 0.7, color = "black")
    ax2.tick_params(axis='y', colors='darkgreen')
    plt.ylabel("Surface pressure [hPa]", color = 'darkgreen')
    #plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/8_Sensitivity_albedo/tsurf_ps_' + str(m) + '.pdf', bbox_inches='tight')



#%%
### Pressure barrier dynamics


# September surface pressure zoom
plt.figure()
plt.plot(lats,ps_z_10[9-1], label = "10%", color = clrs[1])
plt.plot(lats,ps_z_20[9-1], label = "20%", color = clrs[3])
plt.plot(lats,ps_z_30[9-1], label = "30%", color = clrs[5])
plt.legend(frameon = False,loc='center left', bbox_to_anchor=(1, 0.5), title = "Albedo")
plt.xlim([-7,12.5])
plt.ylim([1000,1006])
plt.xlabel("Latitude")
plt.ylabel("September surface pressure [hPa]")
plt.axvline(0, color = "black", linewidth = 0.7, linestyle = '--')
plt.axvline(10, color = "black", linewidth = 0.7)
#plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/8_Sensitivity_albedo/ps_sep_zoom.pdf", bbox_inches='tight')

# September precipitation
plt.figure()
plt.plot(lats,precip_z_10[9-1], label = " 10%", color = clrs[1])
plt.plot(lats,precip_z_20[9-1], label = " 20%", color = clrs[3])
plt.plot(lats,precip_z_30[9-1], label = " 30%", color = clrs[5])
plt.legend(frameon = False,loc='center left', bbox_to_anchor=(1, 0.5), title = "Albedo")
plt.xlim([-15,30])
plt.ylim([0,16])
plt.xlabel("Latitude")
plt.ylabel("September precipitation [mm/day]")
plt.axvline(0, color = "black", linewidth = 0.7, linestyle = '--')
plt.axvline(10, color = "black", linewidth = 0.7)
#plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/8_Sensitivity_albedo/precip_sep_zoom.pdf", bbox_inches='tight')


#%% 
### Hadley cell and ITCZ

#Calculate ITCZ 
for i in range(0,len(albedo)):
    exec("ITCZ_" + str(albedo[i]) + " = []")
    for m in range(0,12):
        exec("ITCZ_" + str(albedo[i]) + ".append(lats[np.argmax(precip_z_" + str(albedo[i]) + "[m])])")

# Show ITCZ
plt.figure()
plt.xlabel("Months")
plt.ylabel("Latitude")
plt.axhline(0,linewidth = 0.7, linestyle = '--', color = "black")
for i in range(0,len(albedo)):
    exec("plt.plot(range(1,13),ITCZ_" + str(albedo[i]) + ",label = " + str(albedo[i]) + " )")
plt.legend(frameon = False, bbox_to_anchor = (1,0.5), loc = "center left", title = 'Albedo')
#plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/8_Sensitivity_albedo/ITCZ_stripes_albedo.pdf')


# calculate upper and lower end of Hadley cell
for i in range(0,len(albedo)-5): 
    print(albedo[i])
    exec("trans_nh_" + str(albedo[i]) + " = []") # northern hemisphere
    exec("lats_nh_" + str(albedo[i]) + " = []")
    exec("trans_sh_" + str(albedo[i]) + " = []") # southern hemisphere
    exec("lats_sh_" + str(albedo[i]) + " = []")
    for m in range(0,12):
        print(m)
        exec("olr = olr_z_" + str(albedo[i]) + "[m]")
        trans = []
        for s in range(0,len(olr)-1):
            if (olr[s] < 250 and olr[s+1] > 250) or (olr[s] > 250 and olr[s+1] < 250):
                trans.append(s)
                
        exec("trans_nh_" + str(albedo[i]) + ".append(trans[-1])")
        exec("trans_sh_" + str(albedo[i]) + ".append(trans[0])")
        exec("lats_sh_" + str(albedo[i]) + ".append(lats[trans[0]])")
        exec("lats_nh_" + str(albedo[i]) + ".append(lats[trans[-1]])")
    
# OLR (plot) and marked ends of Hadley cells
# i= 5
# print(albedo[i])
# for m in range(0,12):
#     plt.figure()
#     plt.title(month_list[m])
#     exec("plt.plot(lats,olr_z_" + str(albedo[i]) + "[m])")
#     exec("plt.axvline(lats_sh_" + str(albedo[i]) + "[m], color = 'orange')")
#     exec("plt.axvline(lats_nh_" + str(albedo[i]) + "[m], color = 'orange')")
       
    

plt.figure()
for i in range(0,len(albedo)-5):
    exec("plt.plot(range(1,13),lats_nh_" + str(albedo[i]) + ", color = clrs[i], linestyle = '--')")
    exec("plt.plot(range(1,13),lats_sh_" + str(albedo[i]) + ", color = clrs[i], linestyle = '--')")
    exec("plt.plot(range(1,13),ITCZ_" + str(albedo[i]) + ", color = clrs[i], label = str(albedo[i]) + '%')")
#plt.plot(range(1,13),ITCZ_89, color = clrs[i], label = 'Aquaplanet')
plt.xlabel("Months")
plt.xlim([1,12])
plt.ylim([-70,70])
plt.ylabel("Latitude")
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),frameon=False, title = "Land Albedo")
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/8_Sensitivity_albedo/hadley_albedo.pdf', bbox_inches='tight')



#%%
### Monthly evolution of monsoon rain

# The grid cell weighted average was calculated via processing.py

# Monthly monsoon rainfall for each of the simulations
precip_m=[[9.43049550e-09,2.82568835e-09,5.33765920e-08,1.33819913e-06
,8.33280683e-06,2.43013455e-05,5.24299503e-05,8.60066502e-05
,1.05121188e-04,9.21075261e-05,1.03351149e-05,1.45996182e-07]
,[1.77583388e-08,1.73732073e-09,2.27490293e-08,7.27433928e-07
,5.75952618e-06,1.95070261e-05,4.41161428e-05,7.78838948e-05
,1.01091071e-04,8.20655405e-05,6.66853475e-06,1.95100228e-07]
,[7.01873226e-09,5.87375371e-10,1.08105942e-08,3.25580118e-07
,3.47746823e-06,1.63134428e-05,3.60423619e-05,6.81187885e-05
,9.43657506e-05,7.02711914e-05,4.59343846e-06,7.42866249e-08]
,[5.72319747e-09,1.45338119e-09,5.58726532e-09,1.51231689e-07
,2.29254420e-06,1.34029542e-05,2.88074189e-05,5.62712339e-05
,8.77404600e-05,5.71578057e-05,2.65771882e-06,4.50655371e-08]
,[2.03942943e-08,2.98443470e-09,3.30643668e-09,6.53145804e-08
,1.37682741e-06,9.53620656e-06,2.26912871e-05,4.37604540e-05
,7.44520657e-05,4.27498271e-05,1.53501435e-06,7.45865734e-08]
,[7.45375583e-09,4.28527547e-10,2.21252638e-09,3.48012819e-08
,7.72043734e-07,6.92036065e-06,1.78625232e-05,3.46775414e-05
,5.95542078e-05,2.92829354e-05,7.71819884e-07,5.76674708e-08]
,[1.15086154e-08,3.11396270e-10,1.69880154e-09,2.10588826e-08
,3.63079607e-07,4.70175019e-06,1.44472842e-05,2.64608771e-05
,4.53786997e-05,1.98478319e-05,3.52886417e-07,3.89030497e-08]
,[2.48817345e-09,7.78090425e-10,1.21141575e-09,9.47583167e-09
,1.87187041e-07,2.71576891e-06,1.17471627e-05,2.12450159e-05
,3.40186089e-05,1.15436760e-05,1.11913529e-07,3.07490104e-08]
,[4.31384928e-09,2.15465645e-10,9.06726916e-10,6.08124529e-09
,9.17842442e-08,1.72928856e-06,9.26178291e-06,1.83648408e-05
,2.88427254e-05,7.01796671e-06,3.71665294e-08,4.54375524e-08]
,[3.32435746e-10,1.80156307e-10,5.99976901e-10,3.98782518e-09
,7.84444865e-08,1.31988884e-06,6.53882080e-06,1.38778069e-05
,2.16518874e-05,4.87229181e-06,3.31650760e-08,3.01933198e-08]
,[1.59421587e-09,1.27784422e-10,3.89734550e-10,2.58243271e-09
,3.89633286e-08,7.97259190e-07,5.36536072e-06,1.27250587e-05
,1.96610526e-05,3.40649080e-06,3.57799053e-08,3.79108833e-08]
,[3.79264176e-09,6.53918933e-11,2.70216904e-10,1.74901749e-09
,1.92925764e-08,4.81754398e-07,3.80580991e-06,1.18337712e-05
,1.70893363e-05,1.59733600e-06,6.67159270e-08,3.54965479e-08]
,[2.15381907e-10,3.77113687e-11,1.16050426e-10,9.30045985e-10
,9.45610878e-09,3.73939855e-07,2.80889049e-06,1.03198336e-05
,1.50235555e-05,9.56674739e-07,2.58955080e-08,4.61734011e-08]]

plt.figure()
for i in range(0,len(albedo)-1):
    plt.plot(range(1,13),np.array(precip_m[i])*86400, label =  str(albedo[i]) + '%', color = clrs[i])
plt.ylabel("Monsoon rainfall [mm/day]")
plt.xlabel('Months')
plt.ylim([0,10])
plt.axhline(0,color="black",linewidth = 0.7)
plt.legend(frameon = False)
#plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/8_Sensitivity_albedo/pr_per_month_albedo.pdf", bbox_inches='tight')
  
#%%
### In dependence of albedo

# mr = monsoon region (10-30°N)
# JJAS = June-September

# The grid cell weighted average was calculated via processing.py
albedo = [10,12,14,16,18,20,22,24,26,28,30,32,34]
precip_mr_annual = [3.1682033e-05,2.8171376e-05,2.4466726e-05,2.0711599e-05,1.6355689e-05,1.2495333e-05,9.3021572e-06,6.8011695e-06,5.4502093e-06,4.0339664e-06,3.5060475e-06,2.9112825e-06,2.4638098e-06]
precip_global_annual = [3.1538941e-05,3.1187523e-05,3.0821695e-05,3.0504614e-05,3.0017558e-05,2.9603627e-05,2.9066858e-05,2.8469110e-05,2.8000233e-05,2.7445763e-05,2.7075515e-05,2.6687521e-05,2.6332264e-05]
tsurf_global_annual = [292.59238,292.14255,291.71094,291.2827,290.7962,290.34195,289.7692,289.21555,288.4838,287.68384,286.8092,285.9841,285.14462]
tsurf_mr_annual = [303.65204,302.989,302.4095,301.82495,301.07227,300.21036,299.1172,298.02628,296.49387,294.87775,293.10013,291.3877,289.62616]
tsurf_mr_JJAS = [312.16772,312.12958,311.98196,311.95822,311.6611,311.22034,310.56607,309.71533,308.49274,307.2664,305.60715,303.93744,302.19662]
tsurf_20S_10N_annual = [300.76074,300.5874,300.45844,300.37393,300.16962,299.96222,299.60568,299.25034,298.7516,298.15,297.55856,296.99545,296.391]
tsurf_20S_10N_JJAS = [301.09338,300.91476,300.74734,300.60074,300.3604,300.07532,299.6776,299.28943,298.7676,298.1501,297.53815,296.96994,296.35336]

# tsurf in the monsoon region and a region south of the landstripe (entire year)
plt.figure()
plt.plot(albedo,np.array(tsurf_mr_annual)-273-15, label = "Monsoon region")
plt.plot(albedo,np.array(tsurf_20S_10N_annual)-273.15, label = "South")
plt.legend(frameon = False)

# tsurf in the monsoon region and a region south of the landstripe (JJAS)
plt.figure()
plt.plot(albedo,np.array(tsurf_mr_JJAS)-273-15, label = "Monsoon region")
plt.plot(albedo,np.array(tsurf_20S_10N_JJAS)-273.15, label = "South")
plt.legend(frameon = False)

# Overview figure
fig,ax = plt.subplots()
plt.plot(albedo,np.array(precip_mr_annual)*86400*365, color = "darkblue", label = "Monsoon region", marker = '.')
plt.plot(albedo,np.array(precip_global_annual)*86400*365, color = "darkblue",linestyle = '--', label = 'Global', marker = '.')
plt.xlabel("Albedo [%]")
plt.ylabel("Annual precipitation [mm]", color = "darkblue")
ax.tick_params(axis='y', colors='darkblue')
plt.legend(frameon = False)
ax2=ax.twinx()
ax2.plot(albedo,np.array(tsurf_mr_annual)-273.15, color = 'darkorange',  marker = '.')
ax2.plot(albedo,np.array(tsurf_global_annual)-273.15, color = 'darkorange',linestyle = '--',  marker = '.')
ax2.set_ylabel("Surface Temperature [°C]", color = "darkorange")
ax2.tick_params(axis='y', colors='darkorange')
plt.legend(frameon = False)
#plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/8_Sensitivity_albedo/pr_per_albedo.pdf", bbox_inches='tight')



#%%
### Summary figure

slab = [50,200,500]
precip_mr_annual_slab = [1.6191263e-05,7.7047498e-06,6.6271145e-06]
precip_global_annual_slab = [3.0008980e-05,2.9380648e-05,2.9112593e-05]
tsurf_global_annual_slab = [290.76,291.99,292.09]
tsurf_mr_annual_slab = [300.95865,303.82892,304.17676]

sol =  [1400, 1340, 1320, 1300, 1280, 1260, 1240, 1220, 1200, 1180,1160, 1140, 1120, 1100, 1080, 1060, 1040, 1020, 1000]
sol2 = [1400, 1340, 1320, 1300, 1280, 1260, 1240, 1220, 1200, 1180,1160, 1140, 1120]
precip_global_annual_sol=[3.2474160e-05,2.8455828e-05,2.7057653e-05,2.5746720e-05,2.4620618e-05,2.3629911e-05,2.2477425e-05,2.1166337e-05,2.0031399e-05,1.8814941e-05,1.7509223e-05,1.6284472e-05,1.4775650e-05,1.2362219e-05,1.0786938e-05,5.6803237e-06,5.3893764e-06,5.0107005e-06,4.7828939e-06]
precip_mr_annual_sol=[2.2514931e-05,1.6191263e-05,1.3706567e-05,1.2151729e-05,1.1044361e-05,1.0543251e-05,9.9233202e-06,9.4784600e-06,9.0642925e-06,8.6720893e-06,8.3828591e-06,8.0624013e-06,8.0787577e-06,8.1644612e-06,7.8298262e-06,1.0846883e-05,3.4762204e-06,3.2997664e-06,3.0319554e-06,3.0103740e-06]
tsurf_global_annual_sol=[293.96896,290.76453,288.5692,286.58844,284.75378,282.80124,280.88815,278.747,276.35776,273.90945,271.16425,268.2663,265.47302,261.8714,256.12714,252.24355,244.4436,243.57248,241.95451,240.92987]
tsurf_mr_annual_sol=[305.3422,297.8474,295.15826,292.83908,290.46027,288.18854,285.85425,283.531,281.27203,279.0063,276.51257,274.3762,271.6137,267.39008,264.77637,252.5583,251.75221,250.63507,250.01729]

ppm = [70,140,200,280,400,560,750,1120]
precip_global_annual_ppm=[2.7611337e-05,2.8835879e-05,2.9467377e-05,3.0008980e-05,3.0602685e-05,3.1142197e-05,3.1681437e-05,3.2312710e-05]
precip_mr_annual_ppm=[1.1755689e-05,1.3655961e-05,1.4768744e-05,1.6191263e-05,1.7862567e-05,1.9985961e-05,2.2135247e-05,2.4366376e-05]
tsurf_global_annual_ppm=[286.33853,288.58524,289.7542,290.76453,291.89136,293.01608,293.93185,295.249]
tsurf_mr_annual_ppm=[294.97186,297.93857,299.5395,300.95865,302.4874,304.06342,305.2409,307.09006]

aerosol = [0, 10**(-7), 10**(-6), 2*10**(-6), 4*10**(-6), 6*10**(-6), 8*10**(-6), 10**(-5), 2*10**(-5), 10**(-4)]
aerosol_name = [0,7, 6, '6_2', '6_4', '6_6', '6_8', 5, '5_2',4]
aerosol_name2 = ['0','$10^{-7}$','$10^{-6}$','$2x10^{-6}$','$4x10^{-6}$','$6x10^{-6}$','$8x10^{-6}$','$10^{-5}$','$2x10^{-5}$','$10^{-4}$']
precip_global_annual_aerosol=[3.0008980e-05,2.9977744e-05,2.9738982e-05,2.9374294e-05,2.8781433e-05,2.8154403e-05,2.7674276e-05,2.7267221e-05,2.5944059e-05,2.1938815e-05]
precip_mr_annual_aerosol=[1.6191263e-05,1.5534917e-05,1.3421401e-05,1.0840929e-05,7.6617889e-06,5.4911411e-06,4.0643176e-06,3.0727881e-06,1.0481854e-06,1.0427821e-10]
tsurf_global_annual_aerosol=[290.76453,290.75244,290.4779,290.124,289.56592,288.9163,288.28098,287.64954,284.87186,270.27216]
tsurf_mr_annual_aerosol=[300.95865,300.83334,300.54657,299.90344,299.00665,297.88428,296.5743,295.4022,289.74847,264.77393]

albedo = [10,12,14,16,18,20,22,24,26,28,30,32,34]
precip_mr_annual_albedo = [3.1682033e-05,2.8171376e-05,2.4466726e-05,2.0711599e-05,1.6355689e-05,1.2495333e-05,9.3021572e-06,6.8011695e-06,5.4502093e-06,4.0339664e-06,3.5060475e-06,2.9112825e-06,2.4638098e-06]
precip_global_annual_albedo = [3.1538941e-05,3.1187523e-05,3.0821695e-05,3.0504614e-05,3.0017558e-05,2.9603627e-05,2.9066858e-05,2.8469110e-05,2.8000233e-05,2.7445763e-05,2.7075515e-05,2.6687521e-05,2.6332264e-05]
tsurf_global_annual_albedo = [292.59238,292.14255,291.71094,291.2827,290.7962,290.34195,289.7692,289.21555,288.4838,287.68384,286.8092,285.9841,285.14462]
tsurf_mr_annual_albedo = [303.65204,302.989,302.4095,301.82495,301.07227,300.21036,299.1172,298.02628,296.49387,294.87775,293.10013,291.3877,289.62616]

plt.figure()
#plt.plot(np.array(tsurf_global_annual_slab)-273.15,np.array(precip_mr_annual_slab)*86400*365, label = "Slab" )
plt.plot(np.array(tsurf_global_annual_sol[:13])-273.15,np.array(precip_mr_annual_sol[:13])*86400*365, label = "Solar Constant", marker = ".")
plt.plot(np.array(tsurf_global_annual_ppm)-273.15,np.array(precip_mr_annual_ppm)*86400*365, label = "Carbon Dioxide", marker = ".")
plt.plot(np.array(tsurf_global_annual_aerosol)-273.15,np.array(precip_mr_annual_aerosol)*86400*365, label = "Sulfate Aerosols", marker = ".")
plt.plot(np.array(tsurf_global_annual_albedo)-273.15,np.array(precip_mr_annual_albedo)*86400*365, label = "Land Albedo", marker = ".")
plt.legend(frameon = False)
plt.ylabel("Monsoon rainfall [mm]")
plt.xlabel("Global mean surface temperature [°C]")
plt.xlim([0,25])
#plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/pr_per_GMT.pdf", bbox_inches='tight')


