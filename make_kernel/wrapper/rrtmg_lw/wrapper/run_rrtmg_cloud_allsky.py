# import libraries
import sys
import time
import numpy as np
# import functions
sys.path.insert(1,'/home/cls/projects/ctb-itan/cls/KERNELS/RRTMG_wrapper/rrtmg_lw/wrapper/') # replace with full path to ./wrapper
import single_driver_band as driver_clr
import single_driver_band_allsky as driver

t0 = time.time()

### Helpers
class Dummy:
    pass

def run_rrtmg_specific(pmid_tmp,tmid_tmp,h2o_raw,o3_raw,psfc,tsfc,cc_raw,ci_raw,cl_raw,rei_raw,rel_raw,nlev):

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
    """

    params     = Dummy()
    gas        = Dummy()
    atmprofile = Dummy()

    emis = 1

    # set well-mixed GHG
    gas.n2o = 0.323 * 1e-6
    gas.co  = 0.15  * 1e-6
    gas.ch4 = 1.797 * 1e-6
    gas.o2  = 0.209
    gas.co2 = 380*1 * 1e-6

    atmprofile.h2o_raw  = h2o_raw
    atmprofile.o3_raw   = o3_raw
    atmprofile.tsfc     = np.array([tsfc])
    atmprofile.psfc     = np.array([psfc])
    atmprofile.pmid_tmp = pmid_tmp
    atmprofile.tmid_tmp = tmid_tmp

    atmprofile.semiss = np.array([emis])

    fluxbroad_clr,fluxband_clr,nlev_clr = driver_clr.rrtmg_lw_driver(atmprofile,gas)

    atmprofile.cfrc = cc_raw
    atmprofile.ciwc = ci_raw 
    atmprofile.clwc = cl_raw

    atmprofile.crei = rei_raw
    atmprofile.crel = rel_raw

    fluxbroad,fluxband,nlev = driver.rrtmg_lw_driver_allsky( atmprofile,gas )

    return fluxbroad, fluxband, nlev, fluxbroad_clr,fluxband_clr,nlev_clr

# end