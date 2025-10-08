######################################################################
### LIBRARIES ########################################################
######################################################################
import sys
import numpy as np
import xarray as xr
sys.path.append('/home/cls/projects/ctb-itan/cls/pylib/')
from cls_tools import make_DA

######################################################################
### PROPERTIES #######################################################
######################################################################
bot_min,bot_max,hgt_min,hgt_max = 160., 1670., 538., 1612. # int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]) # 160., 1670., 538., 1612., # 
years,months=np.arange(2014,2020),np.arange(1,13)
nyrs,nmth, = len(years),len(months)

######################################################################
### PROCESS DATA #####################################################
######################################################################
### bin bounds
ds_bins = xr.open_dataset('./stratiform_dataset_histogram_bins_20142019.nc')
### particle size data
ds_dmv = xr.open_dataset('/home/cls/projects/ctb-itan/cls/DATASETS/nsaarsclkazr1kolliasC1.c0_hourly/cloud_dataset_kazr_2000-2020.nc')[['dmv_gcpex']]
### water path data
ds_wp = xr.open_dataset('/home/cls/projects/ctb-itan/cls/CLIMATOLOGY/datasets/DOE_ARM_NSA_atmosphere_dataset.nc')[['iwps_shup','lwps_shup']]
### stratiform data
ds_strat = xr.open_dataset('/home/cls/projects/ctb-itan/cls/CLIMATOLOGY/datasets/DOE_ARM_NSA_stratiform_dataset_only_layer_05p.nc')[['strat_mask','cbot_armb','ctop_armb']]
### surface temperature datasets
ds_sfc = xr.open_dataset('/home/cls/projects/ctb-itan/cls/CLIMATOLOGY/datasets/DOE_ARM_NSA_surfaceobs_2000-2020.nc')[['temp_srfc']]
### big dataset
ds = xr.merge([ds_bins,ds_dmv,ds_wp,ds_strat,ds_sfc]);ds_bins.close();ds_dmv.close();ds_wp.close();ds_strat.close();ds_sfc.close()
ds = ds.rename({'dmv_gcpex':'dmv','iwps_shup':'iwp','lwps_shup':'lwp','cbot_armb':'cbot','ctop_armb':'ctop',})

### COMPUTE DERIVED VARIABLES ########################################
ds['twp']  = ds['iwp']+ds['lwp']
ds['slf']  = ds['lwp']/(ds['lwp']+ds['iwp'])
ds['chgt'] = ds['ctop']-ds['cbot']
ds['slf_bounds_regular'] = np.arange(0,1.1,0.2)

### 2014--2019, CTOP/CHEIGHT #########################################
ds = ds.isel(time=np.where((ds.time.dt.year>=2014)&(ds.time.dt.year<=2019)&\
                                    (ds.cbot.values>=bot_min)&(ds.cbot.values<=bot_max)&\
                                        (ds.chgt.values>=hgt_min)&(ds.chgt.values<=hgt_max))[0]).dropna(dim='time')

### MONTHLY CLOUD NUMBERS ############################################
cloud_number,strat_number = np.ones((nyrs,nmth))*np.nan,np.ones((nyrs,nmth))*np.nan
for yidx,year in enumerate(years):
    for midx,month in enumerate(months):
        cloud_number[yidx,midx] = len(ds.isel(time=np.where((ds.time.dt.year==year)&(ds.time.dt.month==month)&(ds.twp.values>0))[0]).dropna(dim='time')['time'].values)
        strat_number[yidx,midx] = len(ds.isel(time=np.where((ds.time.dt.year==year)&(ds.time.dt.month==month)&(ds.strat_mask.values==1)&(ds.lwp.values>0)&(ds.iwp.values>0))[0]).dropna(dim='time')['time'].values)

### MIXED-PHASE STRATIFORM CLOUDS ####################################
Good_strat = (ds.strat_mask.values==1)&(ds.lwp.values>0)&(ds.iwp.values>0)

######################################################################
### MAKE AND SAVE CLOUD HISTOGRAM ####################################
######################################################################
nslf,ntwp,ndmv=len(ds.slf_bounds_regular.values)-1,len(ds.twp_bounds.values)-1,len(ds.dmv_bounds.values)-1
cloud_histogram = np.zeros((nyrs, nmth, nslf, ndmv, ntwp)) # yaxis: Dmv, xaxis: TWP
# loop through years
for yidx,year in enumerate(years):
    Good_yrs = (ds.time.dt.year==year)
    # loop through months
    for midx,month in enumerate(months):
        Good_mth = (ds.time.dt.month==month)
        # loop through SLF bins
        for sidx in range(nslf):
            if sidx == 5: Good_slf = (ds.slf.values>ds.slf_bounds_regular.values[sidx]) & (ds.slf.values<ds.slf_bounds_regular.values[sidx+1])
            else: Good_slf = (ds.slf.values>ds.slf_bounds_regular.values[sidx]) & (ds.slf.values<=ds.slf_bounds_regular.values[sidx+1])
            # loop through DMV bins
            for didx in range(ndmv):
                if didx == 0: Good_dmv = (ds.dmv.values>=ds.dmv_bounds.values[didx]) & (ds.dmv.values<=ds.dmv_bounds.values[didx+1])
                else: Good_dmv = (ds.dmv.values>ds.dmv_bounds.values[didx]) & (ds.dmv.values<=ds.dmv_bounds.values[didx+1])
                # loop through TWP
                for tidx in range(ntwp):
                    if tidx == 0: Good_twp = (ds.twp.values>=ds.twp_bounds.values[tidx]) & (ds.twp.values<=ds.twp_bounds.values[tidx+1])
                    else: Good_twp = (ds.twp.values>ds.twp_bounds.values[tidx]) & (ds.twp.values<=ds.twp_bounds.values[tidx+1])
                    # monthly cloud histogram count
                    cloud_histogram[yidx,midx,sidx,didx,tidx] = np.nansum(ds.strat_mask.values[Good_strat & Good_yrs & Good_mth & Good_slf & Good_dmv & Good_twp])#/cloud_number[yidx,midx]

### XARRAY DATASET ###################################################
ds_new = xr.Dataset({
    'cloud_number':make_DA(cloud_number,'cloud number','',['year','month'],{'year':years,'month':months},'cloud number'),
    'strat_number':make_DA(strat_number,'strat number','',['year','month'],{'year':years,'month':months},'strat number'),
    'slf_bounds':make_DA(ds.slf_bounds_regular.values,'supercooled liquid fraction','',['nbins'],{'nbins':np.arange(nslf+1)},'supercooled liquid fraction bins'),
    'dmv_bounds':ds.dmv_bounds,
    'twp_bounds':ds.twp_bounds,
    'cloud_count':make_DA(cloud_histogram,'cloud count','',['year','month','slf','dmv','twp'],{'year':years,'month':months,'slf':np.arange(nslf),'dmv':np.arange(ndmv),'twp':np.arange(ntwp)},'cloud count for each bin'),
    # 'annual_mean_cloud_count':make_DA(cloud_histogram_annual_mean,'annual mean cloud count','',['slf','dmv','twp'],{'slf':np.arange(nslf),'dmv':np.arange(ndmv)},'cloud count for each bin'),
    # 'annual_sum_cloud_count':make_DA(cloud_histogram_annual_sum,'annual sum cloud count','',['slf','dmv','twp'],{'slf':np.arange(nslf),'dmv':np.arange(ndmv)},'cloud count for each bin'),
    })
ds_new.to_netcdf('./stratiform_cloud_histogram_20142019_BOT'+str(bot_min).zfill(4)+'_TOP'+str(bot_min+hgt_min).zfill(4)+'.nc')
ds_new.close()

######################################################################
### END ##############################################################
######################################################################
"""
python make_cloud_histograms_annual.py 160 1670 538 1612;python make_cloud_histograms_annual.py 1670 3180 1285 3120;python make_cloud_histograms_annual.py 3180 4690 353 593
"""
