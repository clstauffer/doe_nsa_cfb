# import libraries
import sys
import time
import xarray
import numpy as np
# import functions
sys.path.insert(1,'/home/cls/projects/ctb-itan/cls/kernels/rrtmg_sw/wrapper/') # replace with full path to ./wrapper
import single_driver_band as driver_clr
import single_driver_band_allsky as driver

t0 = time.time()

### Helpers
class Dummy:
    pass

def run_rrtmg_specific(pmid_tmp,tmid_tmp,h2o_raw,o3_raw,psfc,tsfc,cc_raw,ci_raw,cl_raw,rei_raw,rel_raw,nlev,sza,month):

    """
    All profiles are bottom-up
    pmid_tmp... pressure profile ..................... hPa
    tmid_tmp... temperature profile .................. K
    h20_raw.... mass mixing ratio profile ............ g/g
    o3_raw..... ozone mixing ratio profile ........... g/g
    psfc....... surface pressure ..................... hPa
    tsfc....... surface temperature .................. K
    cc_raw..... cloud profile ........................ array of 1s and 0s, 1 for cloudy pixel
    ci_raw..... cloud ice water content profile ...... g/g
    cl_raw..... cloud liquid water content profile ... g/g
    rei_raw.... ice effective radii .................. um, must be between 2.5-60
    rel_raw.... liquid effective radii ............... um, must be between 13-130
    sza........
    month...... number 1-12
    """

    # MARC 2025 -- removing non_nan stuff

    params = Dummy()
    gas    = Dummy()
    atmprofile    = Dummy()

    # not_nan = (np.isnan(tmid_tmp)==False) & (np.isnan(pmid_tmp)==False)& (np.isnan(pmid_tmp)==False)  &  (np.isnan(h2o_raw)==False) & (np.isnan(o3_raw)==False)

    # nlev = len(pmid_tmp[not_nan])

    emis = 1 # QC A Changer

    # set well-mixed GHG
    gas.n2o = 0.323 * 1e-6
    gas.co  = 0.15  * 1e-6
    gas.ch4 = 1.797 * 1e-6
    gas.o2  = 0.209
    gas.co2 = 380*1* 1e-6

    atmprofile.h2o_raw  = h2o_raw # h2o_raw[not_nan]
    atmprofile.o3_raw   = o3_raw # o3_raw[not_nan]
    atmprofile.tsfc     = np.array([tsfc])
    atmprofile.psfc     = np.array([psfc])
    atmprofile.pmid_tmp = pmid_tmp # pmid_tmp[not_nan]
    atmprofile.tmid_tmp = tmid_tmp # tmid_tmp[not_nan]

    atmprofile.semiss = np.array([emis])
    atmprofile.sza    = np.array([sza])
    atmprofile.month  = np.array([month])

    fluxbroad_clr,fluxband_clr,nlev_clr = driver_clr.rrtmg_lw_driver(atmprofile,gas) # to check

    atmprofile.cfrc = cc_raw
    atmprofile.ciwc = ci_raw 
    atmprofile.clwc = cl_raw

    atmprofile.crei = rei_raw
    atmprofile.crel = rel_raw

    atmprofile.sza = np.array([sza])
    atmprofile.month = np.array([month])

    fluxbroad,fluxband,nlev = driver.rrtmg_lw_driver_allsky( atmprofile,gas )

    return fluxbroad, fluxband, nlev, fluxbroad_clr,fluxband_clr,nlev_clr

# end