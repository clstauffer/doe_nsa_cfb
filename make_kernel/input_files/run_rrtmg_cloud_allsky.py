# import libraries
import sys
import time
import xarray
import numpy as np
# import functions
sys.path.insert(1,'/home/cls/projects/ctb-itan/cls/kernels/rrtmg_lw/wrapper/')
import single_driver_band as driver_clr
import single_driver_band_allsky as driver

t0 = time.time()

### Helpers
class Dummy:
    pass

def run_rrtmg_specific(data_dir, month,cloud_base, cloud_top, rei_base, rel_base, total_water_path,slf):

    params = Dummy()
    gas    = Dummy()
    atmprofile    = Dummy()

    outname        = str(month).zfill(2)+'_by_monthlymean_clearsky.nc'
    outname_flx    = 'flux_'+ outname
    outname_tau    = 'tau_' + outname
    outname_taucld = 'tau_cloud_' + outname
    outname_frc    = 'frac_' + outname

    # load data
    PTsfcroot   = xarray.open_dataset(data_dir+'temp_pres_surface/'+'temp_pres_surface_%s.nc'%(month))
    h2oTsroot   = xarray.open_dataset(data_dir+'temp_rh_profiles/'+'temp_rh_profile_%s.nc'%(month))
    o3root      = xarray.open_dataset(data_dir+ 'o3_profiles/'+'o3_profile_%s.nc'%(month))

    level       = h2oTsroot['lev'][::-1].data
    pmid_tmp    = h2oTsroot['lev'][::-1].data # make coordinate bottom-up
    tmid_tmp    = h2oTsroot['temp'][::-1].data + 273 # use .data to retrieve data but not mask. Can also use .compressed()
    h2o_raw     = h2oTsroot['mmr'][::-1].data
    o3_raw      = o3root['o3'][::-1].data

    psfc        = PTsfcroot['pres'].data
    tsfc        = PTsfcroot['temp'].data + 273

    not_nan = (np.isnan(tmid_tmp)==False) & (np.isnan(pmid_tmp)==False)& (np.isnan(level)==False)  &  (np.isnan(h2o_raw)==False) & (np.isnan(o3_raw)==False)

    nlev = len(level[not_nan])
    if nlev==33: nlev=32

    emis = 1 # QC A Changer

    # set well-mixed GHG
    gas.n2o = 0.323 * 1e-6
    gas.co  = 0.15  * 1e-6
    gas.ch4 = 1.797 * 1e-6
    gas.o2  = 0.209
    gas.co2 = 380*1* 1e-6

    atmprofile.h2o_raw  = h2o_raw[not_nan]
    atmprofile.o3_raw   = o3_raw[not_nan]
    atmprofile.tsfc     = np.array([tsfc])
    atmprofile.psfc     = np.array([psfc])
    atmprofile.pmid_tmp = pmid_tmp[not_nan]
    atmprofile.tmid_tmp = tmid_tmp[not_nan]

    atmprofile.semiss = np.array([emis]) #Check what is it ?!

    fluxbroad_clr,fluxband_clr,nlev_clr = driver_clr.rrtmg_lw_driver(atmprofile,gas) # to check

    ci_base = total_water_path*(1-slf)
    cl_base = total_water_path*slf

    #----- Cloudy situations
    slab    = cloud_base # QC
    nfrac   = 1 # To check
    cc_raw  = np.full( nlev, 0. )
    ci_raw  = np.full( nlev, 0. )
    cl_raw  = np.full( nlev, 0. )
    rei_raw = np.full( nlev, 0. )
    rel_raw = np.full( nlev, 0. )

    cc_raw[slab:cloud_top]  = nfrac # *0.1  # cloud fraction Y-T why 27
    ci_raw[slab:cloud_top]  = nfrac*ci_base # cloud ice content at slab cloud layer
    cl_raw[slab:cloud_top]  = nfrac*cl_base
    rei_raw[slab:cloud_top] = rei_base
    rel_raw[slab:cloud_top] = rel_base

    atmprofile.cfrc = cc_raw
    atmprofile.ciwc = ci_raw 
    atmprofile.clwc = cl_raw

    atmprofile.crei = rei_raw
    atmprofile.crel = rel_raw

    fluxbroad,fluxband,nlev = driver.rrtmg_lw_driver_allsky( atmprofile,gas )

    return fluxbroad, fluxband, nlev, fluxbroad_clr,fluxband_clr,nlev_clr

# end