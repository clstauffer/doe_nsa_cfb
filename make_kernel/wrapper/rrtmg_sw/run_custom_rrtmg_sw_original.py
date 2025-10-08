######################################################################
### LIBRARIES ########################################################
######################################################################
import sys
import numpy as np
import xarray as xr
# import functions
sys.path.insert(1,'/home/cls/projects/ctb-itan/cls/KERNELS/RRTMG_wrapper/rrtmg_sw/wrapper/') # replace with full path to ./wrapper
import run_rrtmg_cloud_allsky
import solar_zenithal_angle

def make_DA(data,longname,units,dims,coords,description):
    da = xr.DataArray(
        data=data,
        dims=dims,
        coords=coords,
        attrs=dict(
            long_name=longname,
            description=description,
            units=units,
        ),)
    return da

sw_bins  = np.array([2925.,  3625.,  4325.,  4900.,  5650.,  6925.,  7875., 10450., 14425., 19325., 25825., 33500., 44000.,  1710.])
def main_run_rrtmg(month,cloud_bot,cloud_hgt,twp,slf,rei,rel,outfile):
    ######################################################################
    ### PROCESS PROFILES -- All profiles are bottom-up ###################
    ######################################################################
    # pmid_tmp... pressure profile ..................... hPa
    # tmid_tmp... temperature profile .................. K
    # h20_raw.... mass mixing ratio profile ............ g/g
    # o3_raw..... ozone mixing ratio profile ........... g/g
    # psfc....... surface pressure ..................... hPa
    # tsfc....... surface temperature .................. K
    # cc_raw..... cloud profile ........................ array of 1s and 0s, 1 for cloudy level
    # ci_raw..... cloud ice water content profile ...... g/g
    # cl_raw..... cloud liquid water content profile ... g/g
    # rei_raw.... ice effective radii .................. um, must be between 2.5-60
    # rel_raw.... liquid effective radii ............... um, must be between 13-130

    ds_sfc = xr.open_dataset('/home/cls/projects/ctb-itan/cls/KERNELS/input_files/temp_pres_surface/temp_pres_surface_'+str(month)+'.nc')
    psfc = ds_sfc['pres'].values       # hPa
    tsfc = ds_sfc['temp'].values + 273 # K
    ds_sfc.close()

    ds_prf = xr.open_dataset('/home/cls/projects/ctb-itan/cls/KERNELS/input_files/temp_rh_profiles/temp_rh_profile_'+str(month)+'.nc')
    pmid_tmp = ds_prf['lev'].values[::-1]        # hPa
    tmid_tmp = ds_prf['temp'].values[::-1] + 273 # K
    h2o_raw  = ds_prf['mmr'].values[::-1]        # g/g
    ds_prf.close()

    ds_o3p = xr.open_dataset('/home/cls/projects/ctb-itan/cls/KERNELS/input_files/o3_profiles/o3_profile_'+str(month)+'.nc')
    o3_raw = ds_o3p['o3'].values[::-1] # g/g
    ds_o3p.close()

    ds = xr.open_dataset('/home/cls/projects/ctb-itan/cls/KERNELS/input_files/z_pres_profiles/z_pres_profile_'+str(month)+'.nc')
    altitude = ds.alt.values
    ds.close()

    nlev = len(altitude)
    pmid_tmp = pmid_tmp[:nlev]
    tmid_tmp = tmid_tmp[:nlev]
    h2o_raw  = h2o_raw[:nlev]
    o3_raw   = o3_raw[:nlev]

    ##################################################################
    ### SLAB CLOUD ###################################################
    ##################################################################

    cloud_top = cloud_bot+cloud_hgt
    cloud_bot_idx = int(np.argmin(abs(altitude-cloud_bot)))
    cloud_top_idx = int(np.argmin(abs(altitude-cloud_top)))

    cc_raw = np.zeros(nlev)
    cc_raw[cloud_bot_idx:cloud_top_idx] = 1 # 1/0 slab cloud fraction
    ci_raw   = cc_raw*(twp*(1-slf))         # g/m2 ice water path
    cl_raw   = cc_raw*(twp*slf)             # g/m2 liquid water path
    rei_raw  = cc_raw*rei                   # um ice effective radii
    rel_raw  = cc_raw*rel                   # um liquid effective radii

    ######################################################################
    ### RUN SHORTWAVE RRTMG ##############################################
    ######################################################################
    sza_range = solar_zenithal_angle.range_theta(month,71)
    if len(sza_range)>0:
        for sidx2,sza in enumerate(sza_range): 
            # RUN RRTMG ######################################
            fluxbroad_all, fluxband_all, nlev_all, fluxbroad_clr,fluxband_clr,nlev_clr = run_rrtmg_cloud_allsky.run_rrtmg_specific(pmid_tmp,tmid_tmp,h2o_raw,o3_raw,psfc,tsfc,cc_raw,ci_raw,cl_raw,rei_raw,rel_raw,nlev,sza,month)
            # print(fluxband_clr.fdn_band,np.nansum(fluxband_clr.fdn_band,axis=0),fluxbroad_clr.fdn)
            if sidx2==0:
                fluxband_all_fdn_band = fluxband_all.fdn_band[np.newaxis,...]
                fluxband_all_fup_band = fluxband_all.fup_band[np.newaxis,...]
                fluxband_clr_fdn_band = fluxband_clr.fdn_band[np.newaxis,...]
                fluxband_clr_fup_band = fluxband_clr.fup_band[np.newaxis,...]
            else:
                fluxband_all_fdn_band = np.append(fluxband_all_fdn_band,fluxband_all.fdn_band[np.newaxis],axis=0)
                fluxband_all_fup_band = np.append(fluxband_all_fup_band,fluxband_all.fup_band[np.newaxis],axis=0)
                fluxband_clr_fdn_band = np.append(fluxband_clr_fdn_band,fluxband_clr.fdn_band[np.newaxis],axis=0)
                fluxband_clr_fup_band = np.append(fluxband_clr_fup_band,fluxband_clr.fup_band[np.newaxis],axis=0)
        fluxband_all_fdn_band = np.nanmean(fluxband_all_fdn_band,axis=0)
        fluxband_all_fup_band = np.nanmean(fluxband_all_fup_band,axis=0)
        fluxband_clr_fdn_band = np.nanmean(fluxband_clr_fdn_band,axis=0)
        fluxband_clr_fup_band = np.nanmean(fluxband_clr_fup_band,axis=0)
    else:
        fluxband_all_fdn_band = np.ones((len(sw_bins),len(pmid_tmp)))*np.nan
        fluxband_all_fup_band = np.ones((len(sw_bins),len(pmid_tmp)))*np.nan
        fluxband_clr_fdn_band = np.ones((len(sw_bins),len(pmid_tmp)))*np.nan
        fluxband_clr_fup_band = np.ones((len(sw_bins),len(pmid_tmp)))*np.nan

    ######################################################################
    ### SAVE DATA TO NETCDF ##############################################
    ######################################################################
    fluxband_all_fdn_band = make_DA(fluxband_all_fdn_band,'downwelling shortwave flux','W m-2',['sw','pres'],{'sw':sw_bins,'pres':pmid_tmp[::-1]},'all sky downwelling shortwave broadband flux')
    fluxband_all_fup_band = make_DA(fluxband_all_fup_band,  'upwelling shortwave flux','W m-2',['sw','pres'],{'sw':sw_bins,'pres':pmid_tmp[::-1]},'all sky upwelling shortwave broadband flux')
    fluxband_clr_fdn_band = make_DA(fluxband_clr_fdn_band,'downwelling shortwave flux','W m-2',['sw','pres'],{'sw':sw_bins,'pres':pmid_tmp[::-1]},'clear sky downwelling shortwave broadband flux')
    fluxband_clr_fup_band = make_DA(fluxband_clr_fup_band,  'upwelling shortwave flux','W m-2',['sw','pres'],{'sw':sw_bins,'pres':pmid_tmp[::-1]},'clear sky upwelling shortwave broadband flux')
    ds = xr.Dataset({
                    'sw_flux_down_all':fluxband_all_fdn_band,
                    'sw_flux_down_clr':fluxband_clr_fdn_band,
                    'sw_flux_up_all':fluxband_all_fup_band,
                    'sw_flux_up_clr':fluxband_clr_fup_band,
                    })
    ds.attrs['description'] = "RRTM shortwave flux output" # add to description profiles tested
    ds.to_netcdf(outfile) # add to file name the type of profiles used
    ds.close()

######################################################################
### USER PARAMETERS ##################################################
######################################################################
month     = int(str(sys.argv[1]))
cloud_bot = int(str(sys.argv[2]))
cloud_hgt = int(str(sys.argv[3]))

######################################################################
### HISTOGRAM BIN DATA ###############################################
######################################################################
dmv_bounds = np.array([0.0275, 0.0548, 0.0821, 0.1095, 0.1368]) # YI RELATIONSHIP V2
eff_bounds = np.array([0.64,13.53,26.42,39.31,52.20]) # YI RELATIONSHIP V2
# dmv_bounds = np.array([0.,0.026, 0.031, 0.042, 0.062, 0.145])
twp_bounds = np.array([0.,  27.,   58.,  105.,  180., 2250.])
slf_bounds = np.array([0., 0.30,  0.75,  0.95,  0.999, 1.00])
"""
eff_bounds = np.array([
    [ np.nan,  np.nan,  np.nan,  1.0,  5.6,  9.8, 14. ],
    [ np.nan,  np.nan,  np.nan,  2.4,  8.0, 13.6, 20. ],
    [ np.nan,  np.nan,  np.nan,  3.4,  9.2, 15.1, 21. ],
    [ np.nan,  np.nan,  np.nan,  6.0, 13.7, 21.3, 30. ],
    [ np.nan,  np.nan,     2.0, 13.9, 25.3, 36.7, 50. ]
    ])
"""
# eff_bounds = np.array([
#     [    np.nan,     np.nan,     np.nan,  1.7139,  6.2081, 10.7024, 15.1966],
#     [    np.nan,     np.nan,     np.nan,  1.7691,  6.1515, 10.5339, 14.9162],
#     [    np.nan,     np.nan,     np.nan,  3.0717,  9.0238, 14.9759, 20.928 ],
#     [    np.nan,     np.nan,     np.nan,  4.5141, 11.2884, 18.0626, 24.8369],
#     [    np.nan,     np.nan,     np.nan,  8.0161, 16.3709, 24.7256, 33.0803],
#     [    np.nan,     np.nan,     7.2131, 18.875 , 30.5369, 42.1987, 53.8606]])

######################################################################
### THEORETICAL ICE/LIQ EFFECTIVE RADII PROPERTIES ###################
######################################################################
rel_bins = np.arange(2.5,60,.1)
rei_bins = np.arange(13,130,.1)
unique_combinations = np.array(np.meshgrid(rel_bins, rei_bins)).T.reshape(-1,2)

######################################################################
### LOOP THROUGH BINS ################################################
######################################################################
for sidx, slf in enumerate(slf_bounds): # SLF BINS
# for sidx, slf in zip([9],[0.999]): # SLF BINS
    for tidx, twp in enumerate(twp_bounds): # TWP BINS
        for didx,dmv in enumerate(dmv_bounds): # DMV BINS
            eff_bins = slf*unique_combinations[:,0] + unique_combinations[:,1]*(1-slf)
            for ridx,eff in enumerate(eff_bounds[didx,:]): # EFF BINS
                if (~np.isnan(eff)) & (len(eff_bins[((eff-2.5) <= eff_bins) & ((eff+2.5)>=eff_bins)]) > 1):
                    # THEORETICAL EFFECTIVE RADII FOR RRTMG INPUT ####
                    where_min = np.argwhere(np.abs(eff_bins - eff)==np.nanmin(np.abs(eff_bins - eff)))
                    rei = unique_combinations[where_min[0][0]][1]
                    rel = unique_combinations[where_min[0][0]][0]
                    # RUN RRTMG ######################################
                    outfile   = '../data/rrtmg_sw/rrtmg_sw_output_MONTH'+str(month).zfill(2)+'_CBOT'+str(cloud_bot).zfill(4)+'_CHGT'+str(cloud_hgt).zfill(4)+'_SLF'+str(sidx).zfill(2)+'_TWP'+str(tidx).zfill(2)+'_DMV'+str(didx).zfill(2)+'_REF'+str(ridx).zfill(2)+'_v3_MAR2025.nc'
                    main_run_rrtmg(month,cloud_bot,cloud_hgt,twp,slf,rei,rel,outfile)
                    print(outfile)
                else:
                    ds_prf = xr.open_dataset('/home/cls/projects/ctb-itan/cls/KERNELS/input_files/temp_rh_profiles/temp_rh_profile_'+str(month)+'.nc')
                    pmid_tmp = ds_prf['lev'].values[::-1]        # hPa
                    ds_prf.close()
                    ds = xr.open_dataset('/home/cls/projects/ctb-itan/cls/KERNELS/input_files/z_pres_profiles/z_pres_profile_'+str(month)+'.nc')
                    altitude = ds.alt.values
                    ds.close()
                    nlev = len(altitude)
                    pmid_tmp = pmid_tmp[:nlev]
                    outfile = '../data/rrtmg_sw/rrtmg_sw_output_MONTH'+str(month).zfill(2)+'_CBOT'+str(cloud_bot).zfill(4)+'_CHGT'+str(cloud_hgt).zfill(4)+'_SLF'+str(sidx).zfill(2)+'_TWP'+str(tidx).zfill(2)+'_DMV'+str(didx).zfill(2)+'_REF'+str(ridx).zfill(2)+'_v3_MAR2025.nc'
                    fluxband_all_fdn_band = make_DA(np.ones((len(sw_bins),nlev))*np.nan,'downwelling shortwave flux','W m-2',['sw','pres'],{'sw':sw_bins,'pres':pmid_tmp[::-1]},'all sky downwelling shortwave broadband flux')
                    fluxband_all_fup_band = make_DA(np.ones((len(sw_bins),nlev))*np.nan,  'upwelling shortwave flux','W m-2',['sw','pres'],{'sw':sw_bins,'pres':pmid_tmp[::-1]},'all sky upwelling shortwave broadband flux')
                    fluxband_clr_fdn_band = make_DA(np.ones((len(sw_bins),nlev))*np.nan,'downwelling shortwave flux','W m-2',['sw','pres'],{'sw':sw_bins,'pres':pmid_tmp[::-1]},'clear sky downwelling shortwave broadband flux')
                    fluxband_clr_fup_band = make_DA(np.ones((len(sw_bins),nlev))*np.nan,  'upwelling shortwave flux','W m-2',['sw','pres'],{'sw':sw_bins,'pres':pmid_tmp[::-1]},'clear sky upwelling shortwave broadband flux')
                    ds = xr.Dataset({
                                    'sw_flux_down_all':fluxband_all_fdn_band,
                                    'sw_flux_down_clr':fluxband_clr_fdn_band,
                                    'sw_flux_up_all':fluxband_all_fup_band,
                                    'sw_flux_up_clr':fluxband_clr_fup_band,
                                    })
                    ds.attrs['description'] = "RRTM shortwave flux output" # add to description profiles tested
                    ds.to_netcdf(outfile) # add to file name the type of profiles used
                    ds.close()
                    print(outfile)

######################################################################
### END ##############################################################
######################################################################
# hist_bins  = xr.open_dataset('/home/cls/projects/ctb-itan/cls/HISTOGRAMS/scripts/stratiform_dataset_histogram_bins_20002020.nc')
# dmv_bounds = hist_bins.dmv_bounds.values
# twp_bounds = hist_bins.twp_bounds.values
# slf_bounds = hist_bins.slf_bounds.values
# eff_bounds = hist_bins.eff_bounds.values
# hist_bins.close()
