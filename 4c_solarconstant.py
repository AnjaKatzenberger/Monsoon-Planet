#!/usr/bin/env python3
# -*- coding: utf-8 -*-


#-------------------------------------------------------------
### ANALYSIS OF SOLAR CONSTANT ###
#-------------------------------------------------------------
# Anja Katzeberger, anja.katzenberger@pik-potsdam.de

# This script creates the figures as displayed in the subchapter "Sensitivty to different solar constant values" and some more figures used in the context of the analysis
# The directories of the input files as well as directories for saving the figures must be adapted. 


import math
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import statistics as stat
from matplotlib.legend_handler import HandlerTuple
from scipy.signal import argrelextrema


sol =  [1400, 1340, 1320, 1300, 1280, 1260, 1240, 1220, 1200, 1180,1160, 1140, 1120, 1100, 1080, 1060, 1040, 1020, 1000] # the simulations below 1100 W/m^2 are unintresting for our monsoon rainfall studies since the planet is increasingly dominated by snow fall rather than rainfall 
sol2 = [1400, 1340, 1320, 1300, 1280, 1260, 1240, 1220, 1200, 1180,1160, 1140, 1120]

dir = []
for i in range(0,len(sol2)):
    print(i)
    dir.append("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/A_data/MonsoonPlanet_sol_" + str(sol[i]) + ".nc")
    print(dir[i])
    exec("ds_" + str(sol[i]) + " = nc.Dataset(dir[i])")
    
    exec("evap_" + str(sol[i]) + " = ds_" + str(sol[i]) + "['evap'][:]")
    exec("precip_" + str(sol[i]) + " = ds_" + str(sol[i]) + "['precip'][:]")
    exec("tsurf_" + str(sol[i]) + " = ds_" + str(sol[i]) + "['t_surf'][:]")
    exec("tref_" + str(sol[i]) + " = ds_" + str(sol[i]) + "['t_ref'][:]")
    exec("ps_" + str(sol[i]) + " = ds_" + str(sol[i]) + "['ps'][:]")
    exec("olr_" + str(sol[i]) + " = ds_" + str(sol[i]) + "['olr'][:]")
    exec("ucomp_" + str(sol[i]) + " = ds_" + str(sol[i]) + "['ucomp'][:]")
    exec("vcomp_" + str(sol[i]) + " = ds_" + str(sol[i]) + "['vcomp'][:]")
    exec("sw_" + str(sol[i]) + " = ds_" + str(sol[i]) + "['swdn_toa'][:]")
    exec("swsur_" + str(sol[i]) + " = ds_" + str(sol[i]) + "['swdn_sfc'][:]") # shortwave down top of atmosphere
    exec("cld_" + str(sol[i]) + " = ds_" + str(sol[i]) + "['cld_amt'][:]") # cloud amount
    exec("cldt_" + str(sol[i]) + " = ds_" + str(sol[i]) + "['tot_cld_amt'][:]") # total cloud amount
    exec("rh_" + str(sol[i]) + " = ds_" + str(sol[i]) + "['rh'][:]") # relative humidity
    exec("sphum_" + str(sol[i]) + " = ds_" + str(sol[i]) + "['sphum'][:]") # specific humidity
    exec("wvp_" + str(sol[i]) + " = ds_" + str(sol[i]) + "['WVP'][:]") # column integrated water vapour



    exec("evap_z_" + str(sol[i]) + " = evap_" + str(sol[i]) + ".mean(axis=2)*86400")
    exec("precip_z_" + str(sol[i]) + " = precip_" + str(sol[i]) + ".mean(axis = 2)*86400")
    exec("olr_z_" + str(sol[i]) + " = olr_" + str(sol[i]) + ".mean(axis=2)")
    exec("vcomp_z_" + str(sol[i]) + " = vcomp_" + str(sol[i]) + ".mean(axis=3)")
    exec("ucomp_z_" + str(sol[i]) + " = ucomp_" + str(sol[i]) + ".mean(axis=3)")
    exec("ps_z_" + str(sol[i]) + " = ps_" + str(sol[i]) + ".mean(axis=2)/100")
    exec("tsurf_z_" + str(sol[i]) + " = tsurf_" + str(sol[i]) + ".mean(axis=2)-273.15")
    exec("tref_z_" + str(sol[i]) + " = tref_" + str(sol[i]) + ".mean(axis=2)-273.15")
    exec("sw_z_" + str(sol[i]) + " = sw_" + str(sol[i]) + ".mean(axis=2)")
    exec("swsur_z_" + str(sol[i]) + " = swsur_" + str(sol[i]) + ".mean(axis=2)")
    exec("cld_z_" + str(sol[i]) + " = cld_" + str(sol[i]) + ".mean(axis=3)")
    exec("cldt_z_" + str(sol[i]) + " = cldt_" + str(sol[i]) + ".mean(axis=2)")
    exec("rh_z_" + str(sol[i]) + " = rh_" + str(sol[i]) + ".mean(axis=3)")
    exec("sphum_z_" + str(sol[i]) + " = sphum_" + str(sol[i]) + ".mean(axis=3)")
    exec("wvp_z_" + str(sol[i]) + " = wvp_" + str(sol[i]) + ".mean(axis=2)")
                 

lats = ds_1300.variables['lat'][:]
lons = ds_1300.variables['lon'][:]
pfull = ds_1300.variables['pfull'][:]


month_list = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
clrs = ['#6EA6CD','#98CAE1','#C2E4EF','#ACD39E','#5AAE61','#FDB366','#F67E4B','#DD3D2D','#A50026','#762A83','#36489A','#4A7BB7']


#%% 
### GENERAL CHARACTERISTICS

# Meridional rainfall distribution
for i in range(0,len(sol2)):
    fig = plt.figure()
    for m in range(0,12):
        plt.xlabel('latitude')
        plt.xlim([-91,91])
        plt.ylim([-5,35])
        plt.xticks(np.arange(-90, 91, step=30))
        plt.ylabel( 'Precipitation (mm/day)')
        plt.axvline(10, linewidth = 0.7, color = 'black')
        plt.axvline(60, linewidth = 0.7, color = 'black')
        plt.axhline(0, linewidth = 0.7, color = 'black')
        exec("plt.plot(lats,precip_z_" + str(sol[i]) + "[m],label=month_list[m], color = clrs[m], linewidth = 1)") #alpha = 1 
    plt.legend(loc= 'center left', frameon = False, bbox_to_anchor = (1,0.5))
    #plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/5_Sensitivity_sol/precip_' + str(sol[i]) + 'rad.pdf' )
  
    
    
# Meridional precipitation distribution (contour)
for i in range(0,len(sol2)):
    print(sol[i])
    plt.figure()
    plt.title(str(sol[i])+' $W/m^2$')
    exec("plt.contourf(range(1,13),lats,precip_z_" + str(sol[i]) + ".transpose(),cmap='Blues',extend='max',levels = [0,2,4,6,8,10,12,14,16,18,20])")
    plt.axhline(0, linewidth = 0.7, color = 'black',linestyle = '--')
    plt.axhline(10, linewidth = 0.7, color = 'black')
    plt.axhline(60, linewidth = 0.7, color = 'black')
    plt.xlabel('Months')
    plt.xlim([1,12])
    plt.ylabel('Latitude')
    plt.ylim([-35,35])
    clb = plt.colorbar()
    clb.ax.set_title('mm/day',fontsize=8)
    #plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/5_Sensitivity_sol/precip_10N_280ppm" + str(sol[i]) + "rad_contour.pdf", bbox_inches='tight')

# Meridional evap distribution
for i in range(0,len(sol2)):
    fig = plt.figure()
    for m in range(0,12):
        plt.xlabel('latitude')
        plt.xlim([-91,91])
        plt.ylim([-5,15])
        plt.xticks(np.arange(-90, 91, step=30))
        plt.ylabel( 'Evaporation (mm/day)')
        plt.axvline(10, linewidth = 0.7, color = 'black')
        plt.axvline(60, linewidth = 0.7, color = 'black')
        plt.axhline(0, linewidth = 0.7, color = 'black')
        exec("plt.plot(lats,evap_z_" + str(sol[i]) + "[m],label=month_list[m], color = clrs[m], linewidth = 1)") #alpha = 1 
    plt.legend(loc= 'center left', frameon = False, bbox_to_anchor = (1,0.5))
    #plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/5_Sensitivity_sol/evap_' + str(sol[i]) + 'rad.pdf' )
  
# Meridional tsurf distribution for all months
for i in range(0,len(sol2)):
    plt.figure()
    plt.title(str(sol[i])+'$W/m^2$')
    exec("plt.contourf(range(1,13),lats,tsurf_z_" + str(sol[i]) + ".transpose(),cmap='OrRd',extend='both',levels = range(0,50))")
    plt.axhline(0, linewidth = 0.7, color = 'black',linestyle = '--')
    plt.axhline(10, linewidth = 0.7, color = 'black')
    plt.axhline(60, linewidth = 0.7, color = 'black')
    plt.xlabel('Months')
    plt.xlim([1,12])
    plt.ylabel('Latitude')
    plt.ylim([-35,35])
    clb = plt.colorbar()
    clb.ax.set_title('°C',fontsize=8)
    plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/00_Publication_figures/5_Sensitivity_sol/tsurf_10N_280ppm" + str(sol[i]) + "rad_contour.pdf", bbox_inches='tight')



#%%
### Pressure barrier dynamics

# August surface pressure zoom
plt.figure()
plt.plot(lats,ps_z_1200[8-1], label = "1200", color = clrs[1])
plt.plot(lats,ps_z_1300[8-1], label = "1300", color = clrs[3])
plt.plot(lats,ps_z_1400[8-1], label = "1400", color = clrs[5])
plt.legend(frameon = False,loc='center left', bbox_to_anchor=(1, 0.5), title = "Solar Constant")
plt.xlim([-7,12.5])
plt.ylim([1000,1007])
plt.xlabel("Latitude")
plt.ylabel("August surface pressure [hPa]")
plt.axvline(0, color = "black", linewidth = 0.7, linestyle = '--')
plt.axvline(10, color = "black", linewidth = 0.7, linestyle = '-')
plt.axvline(10, color = "black", linewidth = 0.7)
#plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/5_Sensitivity_sol/ps_aug_zoom.pdf", bbox_inches='tight')

# September surface pressure zoom
plt.figure()
plt.plot(lats,ps_z_1200[9-1], label = "1200", color = clrs[1])
plt.plot(lats,ps_z_1300[9-1], label = "1300", color = clrs[3])
plt.plot(lats,ps_z_1400[9-1], label = "1400", color = clrs[5])
plt.legend(frameon = False,loc='center left', bbox_to_anchor=(1, 0.5), title = "Solar Constant")
plt.xlim([-7,12.5])
plt.ylim([998,1007])
plt.xlabel("Latitude")
plt.ylabel("September surface pressure [hPa]")
plt.axvline(0, color = "black", linewidth = 0.7, linestyle = '--')
plt.axvline(10, color = "black", linewidth = 0.7, linestyle = '-')
plt.axvline(10, color = "black", linewidth = 0.7)
#plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/5_Sensitivity_sol/ps_sep_zoom.pdf", bbox_inches='tight')

# August precipitation
plt.figure()
plt.plot(lats,precip_z_1200[8-1], label = "1200", color = clrs[1])
plt.plot(lats,precip_z_1300[8-1], label = "1300", color = clrs[3])
plt.plot(lats,precip_z_1400[8-1], label = "1400", color = clrs[5])
plt.legend(frameon = False,loc='center left', bbox_to_anchor=(1, 0.5), title = "Solar Constant")
plt.xlim([-15,30])
plt.ylim([0,15])
plt.xlabel("Latitude")
plt.ylabel("August precipitation [mm/day]")
plt.axvline(0, color = "black", linewidth = 0.7, linestyle = '--')
plt.axvline(10, color = "black", linewidth = 0.7)
#plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/5_Sensitivity_sol/precip_aug_zoom.pdf", bbox_inches='tight')

# September precipitation
plt.figure()
plt.plot(lats,precip_z_1200[9-1], label = "1200 $W/m^2$", color = clrs[1])
plt.plot(lats,precip_z_1300[9-1], label = "1300 $W/m^2$", color = clrs[3])
plt.plot(lats,precip_z_1400[9-1], label = "1400 $W/m^2$", color = clrs[5])
plt.legend(frameon = False,loc='center left', bbox_to_anchor=(1, 0.5), title = "Solar Constant")
plt.xlim([-15,30])
plt.ylim([0,15])
plt.xlabel("Latitude")
plt.ylabel("September precipitation [mm/day]")
plt.axvline(0, color = "black", linewidth = 0.7, linestyle = '--')
plt.axvline(10, color = "black", linewidth = 0.7)
#plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/5_Sensitivity_sol/precip_sep_zoom.pdf", bbox_inches='tight')



### Calculating barrier height in August for different solar constant simulations for discussion

# ps = ps_z_1200[7]
# maxl = argrelextrema(np.array(ps), np.greater)
# minl = argrelextrema(np.array(ps), np.less)
# ls_1200 = ps[maxl[0][1]] - ps[minl[0][0]]
# plt.plot(lats,ps, label = "AUG")
# plt.xlim([-7.5,12.5])
# plt.ylim([1000,1007])
# plt.axhline(ps[maxl[0][1]], color = "red")
# plt.axhline(ps[minl[0][0]], color = "blue")
# plt.show()
# plt.close()

# ps = ps_z_1300[7]
# maxl = argrelextrema(np.array(ps), np.greater)
# minl = argrelextrema(np.array(ps), np.less)
# ls_1300 = ps[maxl[0][1]] - ps[minl[0][1]]
# plt.plot(lats,ps, label = "AUG")
# plt.xlim([-7.5,12.5])
# plt.ylim([1000,1007])
# plt.axhline(ps[maxl[0][1]], color = "red")
# plt.axhline(ps[minl[0][1]], color = "blue")
# plt.show()
# plt.close()

# ps = ps_z_1400[7]
# maxl = argrelextrema(np.array(ps), np.greater)
# minl = argrelextrema(np.array(ps), np.less)
# ls_1400 = ps[maxl[0][1]] - ps[minl[0][1]]
# plt.plot(lats,ps, label = "AUG")
# plt.xlim([-7.5,12.5])
# plt.ylim([1002,1007])
# plt.axhline(ps[maxl[0][1]], color = "red")
# plt.axhline(ps[minl[0][1]], color = "blue")
# plt.show()
# plt.close()



#%%

# Calculate ITCZ 
for i in range(0,len(sol2)):
    exec("ITCZ_" + str(sol[i]) + " = []")
    for m in range(0,12):
        exec("ITCZ_" + str(sol[i]) + ".append(lats[np.argmax(precip_z_" + str(sol[i]) + "[m])])")#-evap_z_" + str(sol[i]) + "[m]

# # Show ITCZ
# plt.figure()
# plt.xlabel("Months")
# plt.ylabel("Latitude")
# plt.axhline(0,linewidth = 0.7, linestyle = '--', color = "black")
# for i in range(0,len(sol2)):
#     exec("plt.plot(range(1,13),ITCZ_" + str(sol[i]) + ",label = " + str(sol[i]) + " )")
# plt.legend(frameon = False,loc='center left', bbox_to_anchor=(1, 0.5), title = "Solar constant")
# #plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/5_Sensitivity_sol/ITCZ_stripes_sol.pdf')


# calculate upper and lower end of Hadley cell 
for i in range(0,len(sol2)-5):
    print(sol2[i])
    exec("trans_nh_" + str(sol[i]) + " = []") # northern hemisphere
    exec("lats_nh_" + str(sol[i]) + " = []")
    exec("trans_sh_" + str(sol[i]) + " = []") # southern hemisphere
    exec("lats_sh_" + str(sol[i]) + " = []")
    for m in range(0,12):
        print(m)
        exec("olr = olr_z_" + str(sol[i]) + "[m]")
        trans = []
        for s in range(0,len(olr)-1):
            if (olr[s] < 250 and olr[s+1] > 250) or (olr[s] > 250 and olr[s+1] < 250):
                trans.append(s)
                
       # print(m+1)
       # print(trans)
        exec("trans_nh_" + str(sol[i]) + ".append(trans[-1])")
        exec("trans_sh_" + str(sol[i]) + ".append(trans[0])")
        exec("lats_sh_" + str(sol[i]) + ".append(lats[trans[0]])")
        exec("lats_nh_" + str(sol[i]) + ".append(lats[trans[-1]])")
    
# OLR (plot) and marked ends of Hadley cells
i= 4
print(sol[i])
for m in range(0,12):
    plt.figure()
    plt.title(month_list[m])
    exec("plt.plot(lats,olr_z_" + str(sol[i]) + "[m])")
    exec("plt.axvline(lats_sh_" + str(sol[i]) + "[m], color = 'orange')")
    exec("plt.axvline(lats_nh_" + str(sol[i]) + "[m], color = 'orange')")


plt.figure()
for i in range(0,len(sol2)-8):
    exec("plt.plot(range(1,13),lats_nh_" + str(sol[i]) + ", color = clrs[i+3], linestyle = '--')")
    exec("plt.plot(range(1,13),lats_sh_" + str(sol[i]) + ", color = clrs[i+3], linestyle = '--')")
    exec("plt.plot(range(1,13),ITCZ_" + str(sol[i]) + ", color = clrs[i+3], label = str(sol2[i]) + '$ W/m^2$')")
#plt.plot(range(1,13),ITCZ_89, color = clrs[i], label = 'Aquaplanet')
plt.xlabel("Months")
plt.xlim([1,12])
plt.ylim([-70,70])
plt.ylabel("Latitude")
leg = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),frameon=False, title = "Solar Constant")
leg.get_title().set_ha("left") 
#plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/5_Sensitivity_sol/hadley_solar.pdf', bbox_inches='tight')

#%%
### Monthly evolution of monsoon rain

# The grid cell weighted average was calculated via processing.py

# Monthly monsoon rainfall for each of the simulations
precip_m = [[1.54424544e-08,6.68710032e-10,4.76042805e-09,4.71379344e-08
,6.92601702e-07,1.16960418e-05,3.00673364e-05,5.95437341e-05
,9.40945320e-05,6.86983767e-05,5.19373043e-06,1.24799200e-07]
,[2.21093810e-09,3.80141280e-10,2.27029839e-09,4.98724262e-08
,9.90895842e-07,7.65865025e-06,2.04978660e-05,3.88597182e-05
,6.54957403e-05,3.02141980e-05,6.57623843e-07,4.93802688e-08]
,[1.74514145e-08,2.77233903e-10,1.65648539e-09,6.74257947e-08
,1.38727319e-06,8.22772563e-06,1.93463729e-05,3.56640012e-05
,5.80518936e-05,2.23622992e-05,6.19810180e-07,7.45557571e-08]
,[1.26741941e-08,3.57799901e-10,1.22410249e-09,1.14527367e-07
,1.95650159e-06,7.12362544e-06,1.79896451e-05,3.42631611e-05
,5.22275441e-05,1.83542652e-05,4.46495477e-07,4.23116866e-08]
,[6.38846576e-09,2.54141402e-10,1.05048481e-09,1.38737050e-07
,2.76129094e-06,8.03487728e-06,1.81667692e-05,3.36827070e-05
,4.81531497e-05,1.49571924e-05,5.71644080e-07,4.49525110e-08]
,[3.35327037e-08,5.90368254e-10,1.00190123e-09,1.84246574e-07
,3.53963583e-06,8.74451416e-06,1.83937791e-05,3.20453400e-05
,4.29624524e-05,1.26003652e-05,5.14735689e-07,5.96508016e-08]
,[3.01206420e-08,8.98830010e-10,1.18086818e-09,2.06409851e-07
,3.64533844e-06,8.73069348e-06,1.82993463e-05,3.08774252e-05
,4.00196986e-05,1.12293574e-05,6.40581106e-07,6.04722104e-08]
,[5.32788498e-08,2.92385094e-09,1.55660451e-09,2.37929484e-07
,4.25108692e-06,9.13307395e-06,1.74285524e-05,3.04822315e-05
,3.65756241e-05,9.91059460e-06,5.78813911e-07,1.15848096e-07]
,[6.17147720e-08,4.85015139e-09,3.47737461e-09,3.00653369e-07
,5.05663093e-06,9.91127217e-06,1.72841264e-05,2.92749755e-05
,3.23695422e-05,8.97563677e-06,7.20936782e-07,1.01257363e-07]
,[8.89991867e-08,1.23281332e-08,9.11535736e-09,5.37620622e-07
,5.69721351e-06,1.02339736e-05,1.68086644e-05,2.76570972e-05
,3.00945685e-05,8.55689450e-06,8.04053229e-07,9.37744531e-08]
,[1.15362532e-07,1.16696421e-08,1.19351498e-08,5.05548769e-07
,5.54304597e-06,1.06020671e-05,1.61199459e-05,2.60399302e-05
,2.79067299e-05,8.52908033e-06,1.22837207e-06,1.35132495e-07]
,[9.55892361e-08,1.88846965e-08,1.38016336e-08,3.84866524e-07
,4.58299155e-06,1.06783073e-05,1.65770198e-05,2.61509103e-05
,2.73835140e-05,9.25752283e-06,1.62529716e-06,1.76391012e-07]
,[1.87349698e-07,8.79774618e-08,2.59179025e-08,3.83608125e-07
,3.19235414e-06,9.31664908e-06,1.67167054e-05,2.60372362e-05
,2.81773137e-05,1.06410580e-05,2.67533665e-06,5.32026093e-07]]

clrs = ['#4A7BB7',
 '#36489A',
 '#762A83',
 '#A50026',
 '#DD3D2D',
 '#F67E4B',
 '#FDB366',
 '#5AAE61',
 '#ACD39E',
 '#C2E4EF',
 '#98CAE1',
 '#6EA6CD', '#4A7BB7','#36489A', '#762A83', 'darkblue']

plt.figure()
for i in [0,1,2,3,4,6,8,10,12]:
    plt.plot(range(1,13),np.array(precip_m[i])*86400, label =  str(sol[i]) + '$ W/m^2$', color = clrs[i+3])
plt.ylabel("Monsoon rainfall [mm/day]")
plt.xlabel('Months')
plt.ylim([0,9])
plt.axhline(0,color="black",linewidth = 0.7)
leg = plt.legend(frameon = False,loc='center left', bbox_to_anchor=(1, 0.5),title = "Solar Constant")
leg.get_title().set_ha("left")
#plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/5_Sensitivity_sol/pr_per_month_sol.pdf", bbox_inches='tight')
  
#%%
### In dependence of solar constant values

# mr = monsoon region (10-30°N)
# JJAS = June-September

# The grid cell weighted average was calculated via processing.py

# evaporation from 90S to 10N
evap_o = [4.8830414e-05,4.3229949e-05,4.1095336e-05,3.9139624e-05,3.7338690e-05,3.5782174e-05,3.3944962e-05,3.1829466e-05,3.0072097e-05,2.8204558e-05,2.6141242e-05,2.4105100e-05,2.1520103e-05]

# evap over land
evap_l = [9.2683558e-06,7.5171074e-06,7.1939753e-06,6.7973579e-06,6.6748225e-06,6.5224081e-06,6.3740140e-06,6.2419022e-06,6.0101397e-06,5.7165821e-06,5.4624143e-06,5.3971830e-06,5.4276043e-06]

# evap north of land 
evap_n = [8.4609028e-06,6.6610792e-06,6.1874371e-06,5.8326395e-06,5.4432066e-06,5.0903900e-06,4.7737221e-06,4.4355538e-06,4.0996297e-06,3.8431203e-06,3.7676341e-06,3.6901733e-06,3.6889612e-06]

# evap global
evap_g = [3.2471329e-05,2.8454389e-05,2.7056820e-05,2.5746442e-05,2.4620242e-05,2.3629655e-05,2.2477787e-05,2.1166637e-05,2.0031446e-05,1.8815070e-05,1.7508926e-05,1.6284404e-05,1.4775448e-05]

plt.figure()
plt.plot(sol2,np.array(evap_o)*86400, label = "90S - 10N")
plt.plot(sol2,np.array(evap_l)*86400, label = "10N - 60N")
plt.plot(sol2,np.array(evap_n)*86400, label = "60N - 90N")
plt.plot(sol2,np.array(evap_g)*86400, label = "90S - 90N")
plt.legend(frameon = False)



# Monsoon rainfall, evaporation and transportable water
evap_global_annual = [3.2471329e-05,2.8454389e-05,2.7056820e-05,2.5746442e-05,2.4620242e-05,2.3629655e-05,2.2477787e-05,2.1166637e-05,2.0031446e-05,1.8815070e-05,1.7508926e-05,1.6284404e-05,1.4775448e-05,1.2362077e-05,1.0786991e-05,5.6802605e-06,5.3893677e-06,5.0106214e-06,4.7827366e-06]
wvp_global_annual = [29.839777,20.022635,17.421024,15.349231,13.506909,11.89889,10.367379,8.956019,7.7106495,6.554062,5.500548,4.6998725,3.835738,2.7155876,2.030468,0.96507204,0.9055032,0.8130128,0.76363003]
precip_mr_annual = [2.2514931e-05,1.3706567e-05,1.2151729e-05,1.1044361e-05,1.0543251e-05,9.9233202e-06,9.4784600e-06,9.0642925e-06,8.6720893e-06,8.3828591e-06,8.0624013e-06,8.0787577e-06,8.1644612e-06,7.8298262e-06,1.0846883e-05,3.4762204e-06,3.2997664e-06,3.0319554e-06,3.0103740e-06]

fig,ax = plt.subplots()
plt.ylim([200,1200])
pr, = ax.plot(sol2,np.array(precip_mr_annual[:13])*86400*365, color = "darkblue", label= "Monsoon rainfall (10-30°N)",marker = ".")
evap, = ax.plot(sol2,np.array(evap_global_annual[:13])*86400*365, color = "darkgreen",linestyle = '--',label= "Evaporation (global)",marker = ".")
plt.ylabel("Water [mm]")
plt.xlabel("Solar Constant [$W/m^2$]")
ax2=ax.twinx()
wvp, = ax2.plot(sol2,np.array(wvp_global_annual[:13]), color = "darkorange",linestyle = '--',label= "Water vapour content (global)",marker = ".")
plt.ylim([0,30])
plt.ylabel("Water [mm]")
plt.xlabel("Solar constant [$W/m^2$]")
ax2.set_ylabel("Transportable water [mm]", color = "darkorange")
ax2.tick_params(axis='y', colors='darkorange')
plt.legend(frameon = False)
plt.legend([pr,evap,wvp], ['Monsoon rainfall (10-30°N)', 'Evaporation (global)','Transportable water (global)'], numpoints=1,handler_map={tuple: HandlerTuple(ndivide=None)}, handlelength=3, frameon = False)
#plt.savefig('/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/5_Sensitivity_sol/evap_contribution.pdf', bbox_inches='tight')



precip_global_annual=[3.2474160e-05,2.8455828e-05,2.7057653e-05,2.5746720e-05,2.4620618e-05,2.3629911e-05,2.2477425e-05,2.1166337e-05,2.0031399e-05,1.8814941e-05,1.7509223e-05,1.6284472e-05,1.4775650e-05,1.2362219e-05,1.0786938e-05,5.6803237e-06,5.3893764e-06,5.0107005e-06,4.7828939e-06]
precip_mr_annual=[2.2514931e-05,1.3706567e-05,1.2151729e-05,1.1044361e-05,1.0543251e-05,9.9233202e-06,9.4784600e-06,9.0642925e-06,8.6720893e-06,8.3828591e-06,8.0624013e-06,8.0787577e-06,8.1644612e-06,7.8298262e-06,1.0846883e-05,3.4762204e-06,3.2997664e-06,3.0319554e-06,3.0103740e-06]
tsurf_global_annual=[293.96896,288.5692,286.58844,284.75378,282.80124,280.88815,278.747,276.35776,273.90945,271.16425,268.2663,265.47302,261.8714,256.12714,252.24355,244.4436,243.57248,241.95451,240.92987]
tsurf_mr_annual=[305.3422,297.8474,295.15826,292.83908,290.46027,288.18854,285.85425,283.531,281.27203,279.0063,276.51257,274.3762,271.6137,267.39008,264.77637,252.5583,251.75221,250.63507,250.01729]

fig,ax = plt.subplots()
pr_mr, = plt.plot(sol2,np.array(precip_mr_annual[:13])*86400*365, color = "darkblue", marker = ".")
pr_g, = plt.plot(sol2,np.array(precip_global_annual[:13])*86400*365, color = "darkblue",linestyle = '--', marker = ".")
plt.xlabel("Solar Constant [$W/m^2$]")
plt.ylabel("Annual precipitation [mm]", color = "darkblue")
ax.tick_params(axis='y', colors='darkblue')
plt.legend(frameon = False)
ax2=ax.twinx()
tsurf_mr,  = ax2.plot(sol2,np.array(tsurf_mr_annual[:13])-273.15, color = 'darkorange', label = "Monsoon region", marker = ".")
tsurf_g, = ax2.plot(sol2,np.array(tsurf_global_annual[:13])-273.15, color = 'darkorange',linestyle = '--', label = "Global", marker = ".")
ax2.set_ylabel("Surface Temperature [°C]", color = "darkorange")
ax2.tick_params(axis='y', colors='darkorange')
plt.legend([(pr_g, tsurf_g), (pr_mr, tsurf_mr)], ['Global', 'Monsoon region'], numpoints=1,handler_map={tuple: HandlerTuple(ndivide=None)}, handlelength=3, frameon = False)

#plt.savefig("/Users/anjakatzenberger/Dokumente/PhD/01_Monsoon_Planet/Figures/Publication_figures/5_Sensitivity_sol/pr_per_sol.pdf", bbox_inches='tight')






