#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#-------------------------------------------------------------
###  PROCESSING DATA WITH CDOs ###
# Anja Katzeberger, anja.katzenberger@pik-potsdam.de
#-------------------------------------------------------------

import subprocess
import netCDF4 as nc
import numpy as np

import numpy.ma as ma
import os.path

# For the different parts of the sensitivity analysis, this code can be adapted in order to 
# A. process the row data by applying CDOs (Climate Data Operators)
# B. Extracting data for a specific area and a specific variable eg- in order to reproduce the data extracts given in the other codes
# Note in this context, that the code needs to be adapted, e.g. regarding the region or  regarding the CDOs, e.g. including timmean

 
#%% 
### 1 LAND POSITION 

latlow = [0,2,4,5,7,10,12,14,15,89]
latlow2 = [0,2,4,6,8,10,12,14,16,89]
width = [50,50,50,50,50,50,50,50,50,1]
year = [71,41,41,71,41,91,91,91,41,41]

latlow = [0,2,4,5,7,10,12,14,15]
latlow2 = [0,2,4,6,8,10,12,14,16]
width = [50,50,50,50,50,50,50,50,50]
year = [71,41,41,71,41,91,91,91,41]


# for i in range(0,len(latlow)):
#     filename = "MonsunPlanet_" + str(latlow[i]) + "N_" + str(width[i]) + "_280ppm/history/00" + str(year[i]) + "0101.atmos_month.nc"
#     outname = "MonsunPlanet_" + str(latlow[i]) + "N_" + str(width[i]) + "_280ppm/history/00" + str(year[i]) + "0101.atmos_month_selyear20_ymonmean_selname_" + str(latlow2[i]) +"N_280ppm.nc"
#     cdo = "cdo -ymonmean -selyear," + str(year[i]-20) + "/" + str(year[i]) + " -selname,evap,precip,olr,ucomp,vcomp,t_surf,t_ref,WVP,ps,swdn_toa,swdn_sfc,tot_cld_amt,cld_amt,rh,sphum " + filename + " " + outname
#     print(cdo)
#     subprocess.check_call(cdo,shell=True)
#     cp = "cp " + outname + " data"
#     subprocess.check_call(cp,shell=True)
   
# for i in range(0,len(latlow)):
#       filename = "data/00" + str(year[i]) + "0101.atmos_month_selyear20_ymonmean_selname_" + str(latlow2[i]) + "N_280ppm.nc"
#       outname = "data/Mon_" + str(latlow2[i]) + "N_280ppm_processed.nc"
#       cdo = "cdo -fldmean -sellonlatbox,0,360,10,30 -selname,precip " + filename + " " + outname 
#       print(cdo)
#       subprocess.check_call(cdo,shell=True)
    
# dir = []
# precip = []
# for i in range(0,len(latlow)):
#     print(latlow[i])
#     exec("dir.append('/p/projects/climber3/anjaka/POEM/work/slab/data/00' + str(year[i]) + '0101.atmos_month_selyear20_ymonmean_selname_' + str(latlow2[i]) + 'N_280ppm_processed.nc')")
#     exec("ds_" + str(latlow[i]) + " = nc.Dataset(dir[i])")
#     exec("precip_ = ds_" + str(latlow[i]) + "['precip'][:]")
#     precip.append(ma.getdata(precip_.data))
#     print(precip_)
# print(np.squeeze(precip))
    



#%% 
# 2 SLAB DEPTHS

slab = [50,200,500]
width = [50,50,50]
latlow = [10,10,10]
year = [91,91,91]


# for i in range(0,len(slab)):
#       filename = "data/00" + str(year[i]) + "0101.atmos_month_selyear20_ymonmean_selname_" + str(latlow[i]) + "N_280ppm_slab" + str(slab[i]) + ".nc"
#       outname = "data/00" + str(year[i]) + "0101.atmos_month_selyear20_ymonmean_selname_" + str(latlow[i]) + "N_280ppm_slab" + str(slab[i]) + "_processed.nc"
#       cdo = "cdo -fldmean -sellonlatbox,0,360,10,30 -selname,precip " + filename + " " + outname 
#       print(cdo)
#       subprocess.check_call(cdo,shell=True)
    
# dir = []
# precip = []
# for i in range(0,len(slab)):
#     print(slab[i])
#     exec("dir.append('/p/projects/climber3/anjaka/POEM/work/slab/data/00' + str(year[i]) + '0101.atmos_month_selyear20_ymonmean_selname_' + str(latlow[i]) + 'N_280ppm_slab' + str(slab[i]) + '_processed.nc')")
    
#     exec("ds_" + str(slab[i]) + " = nc.Dataset(dir[i])")
#     exec("precip_ = ds_" + str(slab[i]) + "['precip'][:]")
#     precip.append(ma.getdata(precip_.data))

#     print(precip_)

# print(np.squeeze(precip))
    




#%% 
# 3 SOLAR CONSTANT
 
sol = [1400, 1340, 1320, 1300, 1280, 1260, 1240, 1220, 1200, 1180,1160, 1140, 1120]


# for i in range(0,len(sol)):
#     filename = "MonsunPlanet_10N_50_280ppm_radiation" + str(sol[i]) + "/history/00910101.atmos_month.nc"
#     outname = "MonsunPlanet_10N_50_280ppm_radiation" + str(sol[i]) + "/history/00910101.atmos_month_selyear20_ymonmean_selname_10N_280ppm_" + str(sol[i]) + "rad.nc"
#     cdo = "cdo -ymonmean -selyear," + str(91-20) + "/" + str(91) + " -selname,evap,precip,olr,ucomp,vcomp,t_surf,t_ref,WVP,ps,swdn_toa,swdn_sfc,tot_cld_amt,cld_amt,rh,sphum " + filename + " " + outname
#     print(cdo)
#     subprocess.check_call(cdo,shell=True)
    
#     cp = "cp " + outname + " data"
#     subprocess.check_call(cp,shell=True)
    

# for i in range(0,len(sol)):
#       filename = "data/00910101.atmos_month_selyear20_ymonmean_selname_10N_280ppm_" + str(sol[i]) + "rad.nc"
#       outname_help = "data/00910101.atmos_month_selyear20_ymonmean_selname_10N_280ppm_" + str(sol[i]) + "rad_0_10N.nc"
#       outname = "data/00910101.atmos_month_selyear20_ymonmean_selname_10N_280ppm_" + str(sol[i]) + "rad_wvp_sum.nc"
#       cdo = "cdo -fldmean -sellonlatbox,0,360,10,30 -selname,t_surf " + filename + " " + outname # -sellonlatbox,0,360,-90,90
#       print(cdo)
#       subprocess.check_call(cdo,shell=True)
    
# dir = []
# precip = []
# for i in range(0,len(sol)):
#     print(sol[i])
#     exec("dir.append('/p/projects/climber3/anjaka/POEM/work/slab/data/00910101.atmos_month_selyear20_ymonmean_selname_10N_280ppm_" + str(sol[i]) + "rad_wvp_sum.nc')")
    
#     exec("ds_" + str(sol[i]) + " = nc.Dataset(dir[i])")
#     exec("precip_ = ds_" + str(sol[i]) + "['t_surf'][:]")
#     precip.append(ma.getdata(precip_.data))

#     print(precip_)

# print(np.squeeze(precip))
    


#%%
# 4 CARBON DIOXIDE

# ppm = [70,140,200,280,400,560,750,1120]


# for i in range(0,len(ppm)):
#     filename = "MonsunPlanet_10N_50_" + str(ppm[i]) + "ppm/history/00610101.atmos_month.nc"
#     outname = "MonsunPlanet_10N_50_" + str(ppm[i]) + "ppm/history/00610101.atmos_month_selyear20_ymonmean_selname_10N_" + str(ppm[i]) + "ppm.nc"
#     cdo = "cdo -ymonmean -selyear," + str(61-20) + "/" + str(61) + " -selname,evap,precip,olr,ucomp,vcomp,t_surf,t_ref,WVP,ps,swdn_toa,swdn_sfc,tot_cld_amt,cld_amt,rh,sphum " + filename + " " + outname
#     print(cdo)
#     subprocess.check_call(cdo,shell=True)
#     cp = "cp " + outname + " data"
#     subprocess.check_call(cp,shell=True) 

# for i in range(0,len(ppm)):
#       filename = "data/00610101.atmos_month_selyear20_ymonmean_selname_10N_" + str(ppm[i]) + "ppm.nc"
#       outname = "data/00610101.atmos_month_selyear20_ymonmean_selname_10N_" + str(ppm[i]) + "ppm_processed.nc"
#       cdo = "cdo -fldmean -sellonlatbox,0,360,10,30 -selname,t_surf " + filename + " " + outname 
#       subprocess.check_call(cdo,shell=True)
      
# dir = []
# precip = []
# for i in range(0,len(ppm)):
#     print(ppm[i])
#     exec("dir.append('data/00610101.atmos_month_selyear20_ymonmean_selname_10N_' + str(ppm[i]) + 'ppm_processed.nc')")
#     exec("ds_" + str(ppm[i]) + " = nc.Dataset(dir[i])")
#     exec("precip_ = ds_" + str(ppm[i]) + "['t_surf'][:]")
#     precip.append(ma.getdata(precip_.data))
#     print(precip_)
# print(np.squeeze(precip))
    
    

#%% 
# 5 SULFATE AEROSOLS

aerosol = [0, 10**(-7), 10**(-6), 2*10**(-6), 4*10**(-6), 6*10**(-6), 8*10**(-6), 10**(-5), 2*10**(-5),10**(-4)]
aerosol_name = [0,7, 6, '6_2', '6_4', '6_6', '6_8', 5, '5_2',4]
year = [101,111,111,111,111,111,111,111,111]


aerosol = [ 10**(-7), 10**(-6), 2*10**(-6), 4*10**(-6), 6*10**(-6), 8*10**(-6), 10**(-5), 2*10**(-5),10**(-4)]
aerosol_name = [7, 6, '6_2', '6_4', '6_6', '6_8', 5, '5_2',4]
year = [101,111,111,111,111,111,111,111,111]



# for i in range(0,len(aerosol_name)):
#     filename = "MonsunPlanet_10N_50_280ppm_aerosol" + str(aerosol_name[i]) + "/history/0" + str(year[i]) + "0101.atmos_month.nc"
#     outname = "MonsunPlanet_10N_50_280ppm_aerosol" + str(aerosol_name[i]) + "/history/0" + str(year[i]) + "0101.atmos_month_selyear20_ymonmean_selname_10N_280ppm_aerosol" + str(aerosol_name[i]) + ".nc"
#     cdo = "cdo -ymonmean -selyear," + str(year[i]-20) + "/" + str(year[i]) + " -selname,evap,precip,olr,ucomp,vcomp,t_surf,t_ref,WVP,ps,swdn_toa,swdn_sfc,tot_cld_amt,cld_amt,rh,sphum " + filename + " " + outname
#     print(cdo)
#     subprocess.check_call(cdo,shell=True)
#     cp = "cp " + outname + " data"
#     subprocess.check_call(cp,shell=True) 


# for i in range(0,len(aerosol)):
#       filename = "data/01110101.atmos_month_selyear20_ymonmean_selname_10N_280ppm_aerosol" + str(aerosol_name[i]) + ".nc"
#       outname = "data/01110101.atmos_month_selyear20_ymonmean_selname_10N_280ppm_aerosol" + str(aerosol_name[i]) + "_processed.nc"
#       cdo = "cdo -fldmean -sellonlatbox,0,360,10,30 -selname,precip " + filename + " " + outname
#       print(cdo)
#       subprocess.check_call(cdo,shell=True)
      
# dir = []
# precip = []
# for i in range(0,len(aerosol)):
#     print(aerosol[i])
#     exec("dir.append('data/01110101.atmos_month_selyear20_ymonmean_selname_10N_280ppm_aerosol' + str(aerosol_name[i]) + '_processed.nc')")
#     exec("ds_" + str(aerosol_name[i]) + " = nc.Dataset(dir[i])")
#     exec("precip_ = ds_" + str(aerosol_name[i]) + "['t_surf'][:]")
#     precip.append(ma.getdata(precip_.data))
#     print(precip_)
# print(np.squeeze(precip))



#%%

# 6 ALBEDO 
albedo = [10,12,14,16,18,20,22,24,26,28,30,32,34]


# for i in range(0,len(albedo)):
#     filename = "MonsunPlanet_10N_50_280ppm_albedo" + str(albedo[i]) + "/history/00910101.atmos_month.nc"
#     outname = "MonsunPlanet_10N_50_280ppm_albedo" + str(albedo[i]) + "/history/00910101.atmos_month_selyear20_ymonmean_selname_10N_280ppm_albedo" + str(albedo[i]) + ".nc"
#     cdo = "cdo -ymonmean -selyear," + str(91-20) + "/" + str(91) + " -selname,evap,precip,olr,ucomp,vcomp,t_surf,t_ref,WVP,ps,swdn_toa,swdn_sfc,tot_cld_amt,cld_amt,rh,sphum " + filename + " " + outname
#     print(cdo)
#     subprocess.check_call(cdo,shell=True)
#     cp = "cp " + outname + " data"
#     subprocess.check_call(cp,shell=True) 

# for i in range(0,len(albedo)):
#       filename = "data/00910101.atmos_month_selyear20_ymonmean_selname_10N_280ppm_albedo" + str(albedo[i]) + ".nc"
#       outname = "data/00910101.atmos_month_selyear20_ymonmean_selname_10N_280ppm_albedo" + str(albedo[i]) + "_processed.nc"
#       cdo = "cdo  -fldmean -sellonlatbox,0,360,-10,10 -selname,t_surf " + filename + " " + outname # -sellonlatbox,0,360,-90,90
#       print(cdo)
#       subprocess.check_call(cdo,shell=True)
      
# dir = []
# precip = []
# for i in range(0,len(albedo)):
#     print(albedo[i])
#     exec("dir.append('data/00910101.atmos_month_selyear20_ymonmean_selname_10N_280ppm_albedo' + str(albedo[i]) + '_processed.nc')")
#     exec("ds_" + str(albedo[i]) + " = nc.Dataset(dir[i])")
#     exec("precip_ = ds_" + str(albedo[i]) + "['t_surf'][:]")
#     precip.append(ma.getdata(precip_.data))
#     print(precip_)
# print(np.squeeze(precip))

