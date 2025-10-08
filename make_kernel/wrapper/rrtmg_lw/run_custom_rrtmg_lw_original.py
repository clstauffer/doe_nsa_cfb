######################################################################
### LIBRARIES ########################################################
######################################################################
import sys
import numpy as np
import xarray as xr
# import functions
sys.path.insert(1,'/home/cls/projects/ctb-itan/cls/KERNELS/RRTMG_wrapper/rrtmg_lw/wrapper/') # replace with full path to ./wrapper
import run_rrtmg_cloud_allsky

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

lw_bins    = np.array([350, 500, 630, 700, 820, 980, 1080, 1180, 1390, 1480, 1800, 2080, 2250, 2380, 2600, 3250])
def main_run_rrtmg(month,cloud_bot,cloud_hgt,twp,slf,rei,rel,outfile):
    ##################################################################
    ### PROCESS PROFILES -- All profiles are bottom-up ###############
    ##################################################################
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
    pmid_tmp = ds_prf['lev'].values[::-1]        # hPa TOA --> SFC
    tmid_tmp = ds_prf['temp'].values[::-1] + 273 # K   TOA --> SFC
    h2o_raw  = ds_prf['mmr'].values[::-1]        # g/g TOA --> SFC
    ds_prf.close()

    ds_o3p = xr.open_dataset('/home/cls/projects/ctb-itan/cls/KERNELS/input_files/o3_profiles/o3_profile_'+str(month)+'.nc')
    o3_raw = ds_o3p['o3'].values[::-1] # g/g
    ds_o3p.close()

    ds = xr.open_dataset('/home/cls/projects/ctb-itan/cls/KERNELS/input_files/z_pres_profiles/z_pres_profile_'+str(month)+'.nc')
    altitude = ds.alt.values # m FSC --> TOA
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

    ### test 16 March 2025: convert TWP to TWC using g/m3/density
    # lay=[0]
    # for idx,alt in enumerate(ds['alt'].values):
    #     lay = np.append(lay,lay[idx]+((alt-lay[idx])*2))
    # twc = np.zeros(nlev)
    # twc[cloud_bot_idx:cloud_top_idx] = twp/len(twc[cloud_bot_idx:cloud_top_idx])
    # twc = twc*len(twc[cloud_bot_idx:cloud_top_idx])/(lay[1:]-lay[:-1])/(100*1000*(pmid_tmp/(287*tmid_tmp)))
    # ci_raw   = cc_raw*(twc*(1-slf))         # g/g ice water content
    # cl_raw   = cc_raw*(twc*slf)             # g/g liquid water content

    ci_raw   = cc_raw*(twp*(1-slf))         # g/m2 ice water path
    cl_raw   = cc_raw*(twp*slf)             # g/m2 liquid water path
    ### end test

    rei_raw  = cc_raw*rei                   # um ice effective radii
    rel_raw  = cc_raw*rel                   # um liquid effective radii

    ##################################################################
    ### RUN LONGWAVE RRTMG ###########################################
    ##################################################################
    fluxbroad_all, fluxband_all, nlev_all, fluxbroad_clr,fluxband_clr,nlev_clr = run_rrtmg_cloud_allsky.run_rrtmg_specific(pmid_tmp,tmid_tmp,h2o_raw,o3_raw,psfc,tsfc,cc_raw,ci_raw,cl_raw,rei_raw,rel_raw,nlev)
    fluxband_all_fdn_band = fluxband_all.fdn_band
    fluxband_all_fup_band = fluxband_all.fup_band
    fluxband_clr_fdn_band = fluxband_clr.fdn_band
    fluxband_clr_fup_band = fluxband_clr.fup_band
    print(fluxband_all_fup_band)

    ##################################################################
    ### SAVE DATA TO NETCDF ##########################################
    ##################################################################
    fluxband_all_fdn_band = make_DA(fluxband_all_fdn_band,'downwelling longwave flux','W m-2',['lw','pres'],{'lw':lw_bins,'pres':pmid_tmp[::-1]},'all sky downwelling longwave broadband flux')
    fluxband_all_fup_band = make_DA(fluxband_all_fup_band,  'upwelling longwave flux','W m-2',['lw','pres'],{'lw':lw_bins,'pres':pmid_tmp[::-1]},'all sky upwelling longwave broadband flux')
    fluxband_clr_fdn_band = make_DA(fluxband_clr_fdn_band,'downwelling longwave flux','W m-2',['lw','pres'],{'lw':lw_bins,'pres':pmid_tmp[::-1]},'clear sky downwelling longwave broadband flux')
    fluxband_clr_fup_band = make_DA(fluxband_clr_fup_band,  'upwelling longwave flux','W m-2',['lw','pres'],{'lw':lw_bins,'pres':pmid_tmp[::-1]},'clear sky upwelling longwave broadband flux')
    ds = xr.Dataset({
                    'lw_flux_down_all':fluxband_all_fdn_band,
                    'lw_flux_down_clr':fluxband_clr_fdn_band,
                    'lw_flux_up_all':fluxband_all_fup_band,
                    'lw_flux_up_clr':fluxband_clr_fup_band,
                    })
    ds.attrs['description'] = "RRTM longwave flux output" # add to description profiles tested
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
dmv_bounds = np.array([0.,0.026, 0.031, 0.042, 0.062, 0.145])
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
eff_bounds = np.array([
    [    np.nan,     np.nan,     np.nan,  1.7139,  6.2081, 10.7024, 15.1966],
    [    np.nan,     np.nan,     np.nan,  1.7691,  6.1515, 10.5339, 14.9162],
    [    np.nan,     np.nan,     np.nan,  3.0717,  9.0238, 14.9759, 20.928 ],
    [    np.nan,     np.nan,     np.nan,  4.5141, 11.2884, 18.0626, 24.8369],
    [    np.nan,     np.nan,     np.nan,  8.0161, 16.3709, 24.7256, 33.0803],
    [    np.nan,     np.nan,     7.2131, 18.875 , 30.5369, 42.1987, 53.8606]])

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
                    outfile   = '../data/rrtmg_lw/rrtmg_lw_output_MONTH'+str(month).zfill(2)+'_CBOT'+str(cloud_bot).zfill(4)+'_CHGT'+str(cloud_hgt).zfill(4)+'_SLF'+str(sidx).zfill(2)+'_TWP'+str(tidx).zfill(2)+'_DMV'+str(didx).zfill(2)+'_REF'+str(ridx).zfill(2)+'_v3_MAR2025.nc'
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
                    outfile = '../data/rrtmg_lw/rrtmg_lw_output_MONTH'+str(month).zfill(2)+'_CBOT'+str(cloud_bot).zfill(4)+'_CHGT'+str(cloud_hgt).zfill(4)+'_SLF'+str(sidx).zfill(2)+'_TWP'+str(tidx).zfill(2)+'_DMV'+str(didx).zfill(2)+'_REF'+str(ridx).zfill(2)+'_v3_MAR2025.nc'
                    fluxband_all_fdn_band = make_DA(np.ones((len(lw_bins),nlev))*np.nan,'downwelling longwave flux','W m-2',['lw','pres'],{'lw':lw_bins,'pres':pmid_tmp[::-1]},'all sky downwelling longwave broadband flux')
                    fluxband_all_fup_band = make_DA(np.ones((len(lw_bins),nlev))*np.nan,  'upwelling longwave flux','W m-2',['lw','pres'],{'lw':lw_bins,'pres':pmid_tmp[::-1]},'all sky upwelling longwave broadband flux')
                    fluxband_clr_fdn_band = make_DA(np.ones((len(lw_bins),nlev))*np.nan,'downwelling longwave flux','W m-2',['lw','pres'],{'lw':lw_bins,'pres':pmid_tmp[::-1]},'clear sky downwelling longwave broadband flux')
                    fluxband_clr_fup_band = make_DA(np.ones((len(lw_bins),nlev))*np.nan,  'upwelling longwave flux','W m-2',['lw','pres'],{'lw':lw_bins,'pres':pmid_tmp[::-1]},'clear sky upwelling longwave broadband flux')
                    ds = xr.Dataset({
                                    'lw_flux_down_all':fluxband_all_fdn_band,
                                    'lw_flux_down_clr':fluxband_clr_fdn_band,
                                    'lw_flux_up_all':fluxband_all_fup_band,
                                    'lw_flux_up_clr':fluxband_clr_fup_band,
                                    })
                    ds.attrs['description'] = "RRTM longwave flux output" # add to description profiles tested
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
