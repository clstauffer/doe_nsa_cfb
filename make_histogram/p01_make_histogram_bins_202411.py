######################################################################
### LIBRARIES ########################################################
######################################################################
import sys
import numpy as np
import pandas as pd
import xarray as xr
from scipy.stats import norm
sys.path.append('/home/cls/projects/ctb-itan/cls/pylib/')
from cls_tools import make_DA

def fit_gaussian(data):
    # mean, std = norm.fit(data[~np.isnan(data)])
    mean = np.nanmean(data)
    std = np.nanstd(data)
    # return np.random.normal(mean, std, 5)
    return std,mean,np.array([np.round(mean-3*std,4),
                     np.round(mean-2*std,4),
                     np.round(mean-1*std,4),
                     np.round(mean,4),
                     np.round(mean+1*std,4),
                     np.round(mean+2*std,4),
                     np.round(mean+3*std,4)])

def gaussian_ypdf(x, mu, sigma):
    return 1 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-(x - mu)**2 / (2 * sigma**2))

def generate_xval(mean, std, num_points=100):
    lower_bound = mean - 3 * std
    upper_bound = mean + 3 * std
    x_values = np.linspace(lower_bound, upper_bound, num_points)
    return x_values

######################################################################
### PROPERTIES #######################################################
######################################################################
start_year,end_year = 2000,2020
years = np.arange(start_year,end_year+1)
months = np.arange(1,13)

######################################################################
### PROCESS DATA #####################################################
######################################################################
### OPEN CLOUD DATASET ###############################################
ds = xr.open_dataset('/home/cls/projects/ctb-itan/cls/DATASETS/DOE_ARM_NSA_cloud_dataset_20002020_mixed_phase_stratiform_clouds.nc')
time  = pd.to_datetime(ds.time.values)
dmv = ds.dmv_avg_gcpex.values
iwp = ds.iwp_avg5km.values
lwp = ds.lwp_avg5km.values
iwc = ds.iwc_avg5km.values
lwc = ds.lwc_avg5km.values
ice = ds.rice_avg5km.values
liq = ds.rliq_avg5km.values
cbot_armb = ds.cbot_armb.values
ctop_armb = ds.ctop_armb.values
chgt_armb = ctop_armb-cbot_armb # ds.cdepth_armb.values
msk = ds.strat_mask.values
ds.close()

### COMPUTE DERIVED VARIABLES ########################################
twp = iwp+lwp
slf = lwp/(lwp+iwp)
eff = (liq*slf) + (ice*(1-slf))

### LIMIT TO YEARS 2014--2019 ########################################
dmv = dmv[np.where((time.year>=2014) & (time.year<=2019))[0]]
ice = ice[np.where((time.year>=2014) & (time.year<=2019))[0]]
liq = liq[np.where((time.year>=2014) & (time.year<=2019))[0]]
msk = msk[np.where((time.year>=2014) & (time.year<=2019))[0]]
twp = twp[np.where((time.year>=2014) & (time.year<=2019))[0]]
slf = slf[np.where((time.year>=2014) & (time.year<=2019))[0]]
eff = eff[np.where((time.year>=2014) & (time.year<=2019))[0]]
bot = cbot_armb[np.where((time.year>=2014) & (time.year<=2019))[0]]
hgt = chgt_armb[np.where((time.year>=2014) & (time.year<=2019))[0]]
top = ctop_armb[np.where((time.year>=2014) & (time.year<=2019))[0]]

time = time[np.where((time.year>=2014) & (time.year<=2019))[0]]

### LIMIT TO STRATIFORM LAYERS #######################################
Good_data = (msk==1) & (twp>0)
dmv = dmv[Good_data]
ice = ice[Good_data]
liq = liq[Good_data]
twp = twp[Good_data]
slf = slf[Good_data]
eff = eff[Good_data]
bot = bot[Good_data]
top = top[Good_data]
hgt = hgt[Good_data]

######################################################################
### VARIABLE DISTRIBUTIONS ###########################################
######################################################################
# 20th percentile of TWP/SLF
twp_bounds_array = np.array([np.nanpercentile(twp,iperc) for iperc in np.arange(0,101,20)])
slf_bounds_array = np.array([np.nanpercentile(slf,iperc) for iperc in np.arange(0,101,20)])

# 20th percentiles of DMV
dmv_bounds_array = np.array([np.nanpercentile(dmv,iperc) for iperc in np.arange(0,101,20)])
dmv_bounds_array[0] = 0

# 5th percentiles of DMV
dmv_p05_array = np.array([np.nanpercentile(dmv,iperc) for iperc in np.arange(0,101,5)])
dmv_p05_array[0] = 0

# 30th percentiles of CBOT
bot_bounds,hgt_bounds = np.ones((2,3))*np.nan,np.ones((2,3))*np.nan
bot_p30_array = np.histogram(bot,bins=3)[1]#np.sort(np.histogram(bot,bins=3)[1])
hgt_p30_array = np.histogram(hgt,bins=3)[1]#np.sort(np.histogram(hgt,bins=3)[1])
for i_cb in range(len(bot_p30_array)-1): 
    Good_cb = (bot>=bot_p30_array[i_cb]) & (bot<bot_p30_array[i_cb+1])
    hist_ch, bins_ch = np.histogram(hgt[Good_cb], bins = np.logspace(np.log10(np.min(hgt[Good_cb])), np.log10(np.max(hgt[Good_cb])),5))
    Good_ch = (hgt>=bins_ch[np.argmax(hist_ch)]) & (hgt<=bins_ch[np.argmax(hist_ch)+1])
    bot_bounds[:,i_cb] = np.array([int(bot_p30_array[i_cb]),int(bot_p30_array[i_cb+1])]) # np.append(bot_bounds,int(bot_p30_array[i_cb]))
    hgt_bounds[:,i_cb] = np.array([int(bins_ch[np.argmax(hist_ch)]),int(bins_ch[np.argmax(hist_ch)+1])]) # np.append(hgt_bounds,int(bins_ch[np.argmax(hist_ch)]))

rei_dmv_0_min = np.nanmin(ice[(ice>0) & (dmv<=dmv_p05_array[2])])
rei_dmv_1_min = np.nanmin(ice[(ice>0) & (dmv>dmv_p05_array[2])  & (dmv<=dmv_p05_array[6]) ])
rei_dmv_2_min = np.nanmin(ice[(ice>0) & (dmv>dmv_p05_array[6])  & (dmv<=dmv_p05_array[10])])
rei_dmv_3_min = np.nanmin(ice[(ice>0) & (dmv>dmv_p05_array[10]) & (dmv<=dmv_p05_array[14])])
rei_dmv_4_min = np.nanmin(ice[(ice>0) & (dmv>dmv_p05_array[14]) & (dmv<=dmv_p05_array[18])])
rei_dmv_5_min = np.nanmin(ice[(ice>0) & (dmv>dmv_p05_array[18])])

rel_dmv_0_min = np.nanmin(liq[(liq>0) & (dmv<=dmv_p05_array[2])])
rel_dmv_1_min = np.nanmin(liq[(liq>0) & (dmv>dmv_p05_array[2])  & (dmv<=dmv_p05_array[6]) ])
rel_dmv_2_min = np.nanmin(liq[(liq>0) & (dmv>dmv_p05_array[6])  & (dmv<=dmv_p05_array[10])])
rel_dmv_3_min = np.nanmin(liq[(liq>0) & (dmv>dmv_p05_array[10]) & (dmv<=dmv_p05_array[14])])
rel_dmv_4_min = np.nanmin(liq[(liq>0) & (dmv>dmv_p05_array[14]) & (dmv<=dmv_p05_array[18])])
rel_dmv_5_min = np.nanmin(liq[(liq>0) & (dmv>dmv_p05_array[18])])

rei_dmv_0_max = np.nanmax(ice[dmv<=dmv_p05_array[2]])
rei_dmv_1_max = np.nanmax(ice[(dmv>dmv_p05_array[2])  & (dmv<=dmv_p05_array[6]) ])
rei_dmv_2_max = np.nanmax(ice[(dmv>dmv_p05_array[6])  & (dmv<=dmv_p05_array[10])])
rei_dmv_3_max = np.nanmax(ice[(dmv>dmv_p05_array[10]) & (dmv<=dmv_p05_array[14])])
rei_dmv_4_max = np.nanmax(ice[(dmv>dmv_p05_array[14]) & (dmv<=dmv_p05_array[18])])
rei_dmv_5_max = np.nanmax(ice[(dmv>dmv_p05_array[18])])

rel_dmv_0_max = np.nanmax(liq[dmv<=dmv_p05_array[2]])
rel_dmv_1_max = np.nanmax(liq[(dmv>dmv_p05_array[2])  & (dmv<=dmv_p05_array[6]) ])
rel_dmv_2_max = np.nanmax(liq[(dmv>dmv_p05_array[6])  & (dmv<=dmv_p05_array[10])])
rel_dmv_3_max = np.nanmax(liq[(dmv>dmv_p05_array[10]) & (dmv<=dmv_p05_array[14])])
rel_dmv_4_max = np.nanmax(liq[(dmv>dmv_p05_array[14]) & (dmv<=dmv_p05_array[18])])
rel_dmv_5_max = np.nanmax(liq[(dmv>dmv_p05_array[18])])

# effective radii "bins" for each Dmv bin
eff_dmv_0_mean,eff_dmv_0_std,eff_dmv_0 = fit_gaussian(eff[dmv<=dmv_p05_array[2]])
eff_dmv_1_mean,eff_dmv_1_std,eff_dmv_1 = fit_gaussian(eff[(dmv>dmv_p05_array[2])  & (dmv<=dmv_p05_array[6]) ])
eff_dmv_2_mean,eff_dmv_2_std,eff_dmv_2 = fit_gaussian(eff[(dmv>dmv_p05_array[6])  & (dmv<=dmv_p05_array[10])])
eff_dmv_3_mean,eff_dmv_3_std,eff_dmv_3 = fit_gaussian(eff[(dmv>dmv_p05_array[10]) & (dmv<=dmv_p05_array[14])])
eff_dmv_4_mean,eff_dmv_4_std,eff_dmv_4 = fit_gaussian(eff[(dmv>dmv_p05_array[14]) & (dmv<=dmv_p05_array[18])])
eff_dmv_5_mean,eff_dmv_5_std,eff_dmv_5 = fit_gaussian(eff[(dmv>dmv_p05_array[18])])

# restrict to realistic values
eff_dmv_0[eff_dmv_0<0] = np.nan
eff_dmv_1[eff_dmv_1<0] = np.nan
eff_dmv_2[eff_dmv_2<0] = np.nan
eff_dmv_3[eff_dmv_3<0] = np.nan
eff_dmv_4[eff_dmv_4<0] = np.nan
eff_dmv_5[eff_dmv_5<0] = np.nan

# final arrays for output
rei_min    = np.array([rei_dmv_0_min,rei_dmv_1_min,rei_dmv_2_min,rei_dmv_3_min,rei_dmv_4_min,rei_dmv_5_min])
rel_min    = np.array([rel_dmv_0_min,rel_dmv_1_min,rel_dmv_2_min,rel_dmv_3_min,rel_dmv_4_min,rel_dmv_5_min])
rei_max    = np.array([rei_dmv_0_max,rei_dmv_1_max,rei_dmv_2_max,rei_dmv_3_max,rei_dmv_4_max,rei_dmv_5_max])
rel_max    = np.array([rel_dmv_0_max,rel_dmv_1_max,rel_dmv_2_max,rel_dmv_3_max,rel_dmv_4_max,rel_dmv_5_max])
eff_bounds = np.array([eff_dmv_0    ,eff_dmv_1    ,eff_dmv_2    ,eff_dmv_3    ,eff_dmv_4,eff_dmv_5])

######################################################################
### SAVE DATASET #####################################################
######################################################################
dmv_bounds_array = make_DA(dmv_bounds_array,'particle size',              'cm',    ['nbins'],{'nbins':np.arange(len(dmv_bounds_array))},'particle size bins')
ice_min = make_DA(rei_min,'minimum ice effective radii',        'um',    ['nbins'],{'nbins':np.arange(len(rei_min))},'minimum ice effective radii bins')
ice_max = make_DA(rei_max,'maximum ice effective radii',        'um',    ['nbins'],{'nbins':np.arange(len(rei_max))},'maximum ice effective radii bins')
liq_min = make_DA(rel_min,'minimum liq effective radii',        'um',    ['nbins'],{'nbins':np.arange(len(rel_min))},'minimum liq effective radii bins')
liq_max = make_DA(rel_max,'maximum liq effective radii',        'um',    ['nbins'],{'nbins':np.arange(len(rel_max))},'maximum liq effective radii bins')
twp_bounds_array = make_DA(twp_bounds_array,'total water path',           'g m-2', ['nbins'],{'nbins':np.arange(len(twp_bounds_array))},'total water path bins')
slf_bounds_array = make_DA(slf_bounds_array,'supercooled liquid fraction','',      ['nbins'],{'nbins':np.arange(len(slf_bounds_array))},'supercooled liquid fraction bins')
eff_dmv_bounds_array = make_DA(eff_bounds,'effective radii',    'um',    ['nbins','nbins2'],{'nbins':np.arange(np.shape(eff_bounds)[0]),'nbins2':np.arange(np.shape(eff_bounds)[1])},'effective radii bins')
bot_bounds = make_DA(bot_bounds,'cloud bottom altitude bins','m',['bnds','calt'],{'bnds':[0,1],'calt':[0,1,2]},'cloud bottom altitude bins')
hgt_bounds = make_DA(hgt_bounds,'cloud height altitude bins','m',['bnds','calt'],{'bnds':[0,1],'calt':[0,1,2]},'cloud height altitude bins')
ds = xr.Dataset({'dmv_bounds':dmv_bounds_array,
                 'ice_min':ice_min,
                 'ice_max':ice_max,
                 'liq_min':liq_min,
                 'liq_max':liq_max,
                 'twp_bounds':twp_bounds_array,
                 'slf_bounds':slf_bounds_array,
                 'eff_bounds':eff_dmv_bounds_array,
                 'bot_bounds':bot_bounds,
                 'hgt_bounds':hgt_bounds,
})
ds.to_netcdf('./stratiform_dataset_histogram_bins_20002020_mixed_phase_stratiform_clouds_v2.nc')
ds.close()

######################################################################
### END ##############################################################
######################################################################
"""
dmv = [dmv_p05_array[0],dmv_p05_array[2],dmv>dmv_p05_array[6],dmv>dmv_p05_array[10],dmv>dmv_p05_array[14],dmv_p05_array[18],dmv_p05_array[-1]]
"""