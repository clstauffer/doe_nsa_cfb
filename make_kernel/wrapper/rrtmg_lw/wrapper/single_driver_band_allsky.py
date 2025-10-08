"""
This is python version of main_lw.m for rrtmg

%%% This is for *clear-sky* calculation with skin temperature for land %%% 
Temporarily set main_lw for single column with all functions done in a single script
Then put functions in corresponding functions
"""
# import libraries
import numpy as np
# import functions
import rrtmg_cld_func_band_N as rrtmg

def rrtmg_lw_driver_allsky( atmprofile,gas ):
    nlev    = np.size(atmprofile.pmid_tmp)
    if nlev==33:
        nlev=32

    atm,cflag = rrtmg.read_profile(nlev,atmprofile )

    # update nlev after cut
    nlev    = np.size(atm.pint)
    co2     = gas.co2 * np.ones(nlev) # * 1e-6
    n2o     = gas.n2o * np.ones(nlev) # * 1e-6
    co      = gas.co  * np.ones(nlev) # * 1e-6
    ch4     = gas.ch4 * np.ones(nlev) # * 1e-6
    o2      = gas.o2  * np.ones(nlev)

    atm.co2 = co2
    atm.n2o = n2o
    atm.co  = co
    atm.ch4 = ch4
    atm.o2  = o2
    
    atm.semiss = atmprofile.semiss

    wbroadl = rrtmg.broad(atm,gas,nlev)

    rrtmg.rrtmg_tape5_writer_htr_lw_band(atm,gas,cflag,wbroadl,nlev)
    rrtmg.rrtmg_lw_htr(atm,nlev)

    fluxbroad,fluxband,nlev = rrtmg.rrtmg_lw_output_read_band(nlev)
    return fluxbroad,fluxband,nlev

# end