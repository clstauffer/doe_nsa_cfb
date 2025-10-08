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
years  = np.arange(2014,2019+1)
months = np.arange(1,13)
limbot = int(sys.argv[1])

cbot_bounds = np.array([160.0,1450.,2740.,4030.])

if limbot >= 0: 
    print(limbot)
    filename = 'BOT_'+str(int(cbot_bounds[limbot])).zfill(4)+'_'+str(int(cbot_bounds[limbot+1])).zfill(4)
else: 
    print(limbot)
    filename = 'BOT_all_altitudes'

######################################################################
### PROCESS DATA #####################################################
######################################################################
# open dataset
ds = xr.open_dataset('/home/cls/projects/ctb-itan/cls/DATASETS/DOE_ARM_NSA_cloud_dataset_20002020.nc')[['dmv_avg_gcpex','iwp_avglayer','lwp_avglayer','cbot_armb','ctop_armb','strat_mask','precip_mask']]
# limit to 2014-2019, drop NaN
ds = ds.isel(time=np.where((ds.time.dt.year>=2014)&(ds.time.dt.year<=2019))[0]).dropna('time')
if limbot==0: 
    print(limbot)
    ds = ds.isel(time=np.where((ds.cbot_armb.values>=cbot_bounds[limbot])&(ds.cbot_armb.values<=cbot_bounds[limbot+1]))[0]).dropna('time')
if limbot> 0: 
    print(limbot)
    ds = ds.isel(time=np.where((ds.cbot_armb.values> cbot_bounds[limbot])&(ds.cbot_armb.values<=cbot_bounds[limbot+1]))[0]).dropna('time')

# hours with mixed-phase stratiform clouds
Good_strat      = ( ds.strat_mask.values  == 1 ) # stratiform layer
Good_precip     = ( ds.precip_mask.values == 0 ) # no precipitation
Good_lwp        = ( ds.lwp_avglayer.values > 0 ) # liquid water is present
Good_iwp        = ( ds.iwp_avglayer.values > 0 ) # ice water is present
ds['mps_cloud'] = ( Good_strat & Good_precip & Good_lwp & Good_iwp ) * 1

lwp = ds.lwp_avglayer.values
iwp = ds.iwp_avglayer.values
dmv = ds.dmv_avg_gcpex.values
mps = ds.mps_cloud.values
bot = ds.cbot_armb.values

twp = iwp + lwp
slf = lwp / twp

Good_cloud = ( twp > 0 ) # liquid or ice water is present

######################################################################
### HISTOGRAM BOUNDS AND BINS DATA ###################################
######################################################################
dmv_bounds = np.array([0,0.04115, 0.06845, 0.0958 , 0.12315,0.145]) # YI RELATIONSHIP V2
# dmv_bounds = np.array([0.000,0.026,0.031,0.042,0.062,0.145])
twp_bounds = np.array([0.000,27.00,58.00,105.0,180.0,2250.])
slf_bounds = np.array([0.000,0.300,0.750,0.950,0.999,1.000])

slf_bins = slf_bounds[:-1]+(slf_bounds[1:]-slf_bounds[:-1])/2
twp_bins = twp_bounds[:-1]+(twp_bounds[1:]-twp_bounds[:-1])/2
# dmv_bins = dmv_bounds[:-1]+(dmv_bounds[1:]-dmv_bounds[:-1])/2
dmv_bins = np.array([0.0275, 0.0548, 0.0821, 0.1095, 0.1368]) # YI RELATIONSHIP

######################################################################
### MAKE HISTOGRAM ###################################################
######################################################################
nyrs,nmth,nslf,ntwp,ndmv = len(years),len(months),len(slf_bins),len(twp_bins),len(dmv_bins)
cloud_histogram = np.ones((nyrs, nmth, nslf, ntwp, ndmv))*np.nan

for yidx,year in enumerate(years): # loop through years
    Good_yrs = (ds.time.dt.year==year)        
    for midx,month in enumerate(months): # loop through months
        if month in np.unique(ds.isel(time=np.where(ds.time.dt.year==year)[0]).time.dt.month):
            Good_mth = (ds.time.dt.month==month)        
            for sidx in range(nslf): # loop through SLF bins (only mixed-phase)
                Good_slf = (slf>slf_bounds[sidx]) & (slf<slf_bounds[sidx+1])
                for tidx in range(ntwp): # loop through TWP bins
                    if tidx == 0: Good_twp = (twp>=twp_bounds[tidx]) & (twp<=twp_bounds[tidx+1])
                    else:         Good_twp = (twp> twp_bounds[tidx]) & (twp<=twp_bounds[tidx+1])
                    for didx in range(ndmv): # loop through DMV bins
                        if didx == 0: Good_dmv = (dmv>=dmv_bounds[didx]) & (dmv<=dmv_bounds[didx+1])
                        else:         Good_dmv = (dmv> dmv_bounds[didx]) & (dmv<=dmv_bounds[didx+1])
                        
                        cloud_histogram[yidx,midx,sidx,tidx,didx] = 100*np.nansum(
                            mps[Good_yrs & Good_mth & Good_slf & Good_twp & Good_dmv]
                            )/np.nansum((Good_yrs & Good_mth & Good_cloud)*1) # count
        else:
            print(yidx,midx,year,month)
cloud_histogram_annual = np.nanmean(cloud_histogram,axis=(0,1))
cloud_histogram_annual_anomaly = cloud_histogram-cloud_histogram_annual[np.newaxis,np.newaxis,...]
cloud_histogram_annual_anomaly_climo = np.nanmean(cloud_histogram_annual_anomaly,axis=0)


######################################################################
### SAVE HISTOGRAM ###################################################
######################################################################
ds = xr.Dataset({
    'year':make_DA(years,'Year','',['year'],{'year':years},'Year'),
    'month':make_DA(months,'Month','',['month'],{'month':months},'Month'),
    'slf':make_DA(slf_bins,'SLF','',['slf'],{'slf':slf_bins},'SLF'),
    'twp':make_DA(twp_bins,'TWP','kg m$^{-3}$',['twp'],{'twp':twp_bins},'TWP'),
    'dmv':make_DA(dmv_bins,'D$_{mv}$','cm',['dmv'],{'dmv':dmv_bins},'D$_{mv}$'),
    'slf_bounds':make_DA(slf_bounds,'SLF','',['slf_bounds'],{'slf_bounds':slf_bounds},'SLF'),
    'twp_bounds':make_DA(twp_bounds,'TWP','',['twp_bounds'],{'twp_bounds':twp_bounds},'TWP'),
    'dmv_bounds':make_DA(dmv_bounds,'D$_{mv}$','',['dmv_bounds'],{'dmv_bounds':dmv_bounds},'D$_{mv}$'),
    'cloud_amount':make_DA(cloud_histogram,
        'Mixed-Phase Stratiform Cloud Amount',
        '%',
        ['year','month','slf','twp','dmv'],
        {'year':years,'month':months,'slf':slf_bins,'twp':twp_bins,'dmv':dmv_bins},
        'Mixed-Phase Stratiform Cloud Amount',
    ),
    'cloud_amount_average':make_DA(cloud_histogram_annual,
        'Average Mixed-Phase Stratiform Cloud Amount',
        '%',
        ['slf','twp','dmv'],
        {'slf':slf_bins,'twp':twp_bins,'dmv':dmv_bins},
        'Average Mixed-Phase Stratiform Cloud Amount',
    ),
    'cloud_amount_anomaly_monthly_climatology':make_DA(cloud_histogram_annual_anomaly_climo,
        'Monthly Climatology\nAnomalous Mixed-Phase Stratiform Cloud Amount',
        '%',
        ['month','slf','twp','dmv'],
        {'month':months,'slf':slf_bins,'twp':twp_bins,'dmv':dmv_bins},
        'Monthly Climatology\nAnomalous Mixed-Phase Stratiform Cloud Amount',
    ),
    'cloud_amount_anomaly':make_DA(cloud_histogram_annual_anomaly,
        'Anomalous Mixed-Phase Stratiform Cloud Amount',
        '%',
        ['year','month','slf','twp','dmv'],
        {'year':years,'month':months,'slf':slf_bins,'twp':twp_bins,'dmv':dmv_bins},
        'Anomalous Mixed-Phase Stratiform Cloud Amount',
    ),
})
ds.to_netcdf('../data/stratiform_cloud_histogram_20142019_'+filename+'_v4_APR2025.nc')
ds.close()

######################################################################
### END ##############################################################
######################################################################
"""
python p02_make_cloud_histograms_202411.py 0;python p02_make_cloud_histograms_202411.py 1;python p02_make_cloud_histograms_202411.py 2;python p02_make_cloud_histograms_202411.py -1
"""
