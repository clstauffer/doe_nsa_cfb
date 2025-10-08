"""
This is the python version of related functions for rrtmg all-sky
"""
# import libraries
import os
import numpy as np
from collections import namedtuple
from scipy.interpolate import interp1d

class Dummy:
    pass

def read_profile( nlev,atmprofile ):
    atm   = Dummy()
    cflag = Dummy()
    # cut unrealistic layers
    n_to_cut = np.size(np.argwhere(atmprofile.pmid_tmp > atmprofile.psfc))

    atm.pmid   = atmprofile.pmid_tmp[n_to_cut:nlev]
    after_size = np.size(atm.pmid)
    atm.tmid   = atmprofile.tmid_tmp[n_to_cut:nlev]

    # add surface layer
    # cut the outermost layer 0.5hpa since I could't do extrap.
    # cut the layer in tape writer
    atm.pint = np.concatenate(((atmprofile.psfc)[0], (atm.pmid[0:after_size-1]+atm.pmid[1:after_size])/2.))#, [0.5]))  # set model top at 0.5hPa
    t_itp    = interp1d(np.log(atm.pmid),atm.tmid)
    atm.tint = np.concatenate(((atmprofile.tsfc)[0],t_itp(np.log(atm.pint[1:]))))
    atm.h2o  = atmprofile.h2o_raw*28.960/18.016
    atm.o3   = atmprofile.o3_raw*28.960/47.998
    atm.h2o  = atm.h2o[n_to_cut:]
    atm.o3   = atm.o3[n_to_cut:]
    atm.cfrc = atmprofile.cfrc[n_to_cut:]
    atm.ciwc = atmprofile.ciwc[n_to_cut:]
    atm.clwc = atmprofile.clwc[n_to_cut:]

    atm.crei = atmprofile.crei[n_to_cut:]
    atm.crel = atmprofile.crel[n_to_cut:]


    cflag.imca = 0
    cflag.icld = 0

    if any(atmprofile.cfrc)>0:
        cflag.icld = 2
        cld_writer( atm.pint,np.size(atm.pint),atm.cfrc,atm.ciwc,atm.clwc,  atm.crei, atm.crel )

    return atm,cflag

def broad(atm,gas,nlev):
    """
    Calculate column density molecules/cm**2
    wbroadl 1 x nlev-1
    from bottom to top
    """

    dp = atm.pint[0:nlev-1] - atm.pint[1:nlev]
    Na = 6.02e23
    gravity = 9.8
    wbroadl = np.zeros(nlev-1)

    for ilayer in range(0,nlev-1):
        amm = (1 - atm.h2o[ilayer]) * 28.966 + atm.h2o[ilayer] * 18.016 # The molecular weight of moist air g/mol
        dry_air = dp[ilayer] * 1e3 * Na /(100 * gravity * amm * (1+atm.h2o[ilayer]))
        summol = gas.co2 + atm.o3[ilayer] + gas.n2o + gas.co + gas.ch4 + gas.o2
        wbroadl[ilayer] = dry_air * (1-summol)

    return wbroadl

def cld_writer( pint,nlev,cfrc,ciwc,clwc, effradice, effradliq):
    inflag = 2 # 0 direct specification of cloud optical depth 2 calculation separate ice and liquid cloud optical depths
    iceflag, liqflag = 2, 1

    dp      = pint[0:nlev-1] - pint[1:nlev]
    indyc   = np.argwhere((ciwc+clwc) > 1e-10)
    ncld    = len(indyc)
    cldfrc  = cfrc[indyc]
    fracice = ciwc[indyc]/(ciwc[indyc]+clwc[indyc])
    crice   = effradice[indyc]
    crliq   = effradliq[indyc]
    # effradice,effradliq = 20,5
    # g = 9.8
    # Han Huang's validation
    # cwp = (1/g) * 1e2 * 1e3 * (ciwc[indyc]+clwc[indyc])*dp[indyc]
    # cwp = (ciwc+clwc)/(float(ncld))*np.ones(ncld)
    cwp =  (ciwc[indyc]+clwc[indyc])/(float(ncld))

    # open IN_CLD_RRTM file for RRTMG
    cfile = 'IN_CLD_RRTM'
    cldID = open(cfile,'w')

    # record C1.1
    cldID.write('%5i%5i%5i\n' % (inflag, iceflag, liqflag))

    # record C1.2
    for i in range(ncld):
        cldID.write('%5i%10.5f%10.5f%10.5f%10.5f%10.5f\n' % (indyc[i]+1,cldfrc[i],cwp[i],fracice[i],crice[i],crliq[i]))

    cldID.write('%s' % '%')
    cldID.close()

def rrtmg_tape5_writer_htr_lw_band(atm,gas,cflag,wbroadl,nlev):
    # write INPUT_RRTMG file for RRTMG
    filename = 'INPUT_RRTM'
    fileID   = open(filename,'w')


    # long wave set record 1.1~basic setting
    rdef1    = namedtuple('record1','iaer,iatm,ixsect,numangs,iout,idrv,imca,icld \
        tbound,iemis,ireflect,semiss')
    record1  = rdef1(0,0,0,0,99,0,cflag.imca,cflag.icld,-99.9,1,0,atm.semiss)
    rdef2    = namedtuple('record2','iform,nlayrs,nmol,pave,tave,pz,tz')
    record2  = rdef2(0,nlev-1,7,atm.pmid,atm.tmid,atm.pint,atm.tint)

    # record 1.1
    fileID.write('%s\n' % '$ ATMOSPHERE PROFILE')

    # record 1.2
    fileID.write('%20d%30d%20d%15d%5d%2d%2d%1d\n' % (record1.iaer, record1.iatm, \
        record1.ixsect, record1.numangs, record1.iout,record1.idrv, \
        record1.imca, record1.icld))

    # record 1.4
    fileID.write('%10.3f %1d %2d %5.3f\n' % (record1.tbound, record1.iemis, \
        record1.ireflect, record1.semiss))

    # record 2.1
    fileID.write('%2d%3d%5d\n' % (record2.iform, record2.nlayrs, record2.nmol))

    # record 2.1.1
    # record 2.1.2

    # first layer
    ilayer = 0
    fileID.write('%10.4f%10.4f%31.3f%7.2f%15.3f%7.2f\n' % (record2.pave[ilayer], record2.tave[ilayer], \
        record2.pz[ilayer], record2.tz[ilayer], \
        record2.pz[ilayer+1], record2.tz[ilayer+1]))
    fileID.write('%10.3e%10.3e%10.3e%10.3e%10.3e%10.3e%10.3e%10.3e\n' % (atm.h2o[ilayer], \
        gas.co2, atm.o3[ilayer], gas.n2o,gas.co, \
        gas.ch4, gas.o2, wbroadl[ilayer]))

    # higher layer
    for ilayer in range(1,nlev-1):
        fileID.write('%10.4f%10.4f%53.3f%7.2f\n' % (record2.pave[ilayer], record2.tave[ilayer], \
            record2.pz[ilayer+1], record2.tz[ilayer+1]))

        fileID.write('%10.3e%10.3e%10.3e%10.3e%10.3e%10.3e%10.3e%10.3e\n' % (atm.h2o[ilayer], \
            gas.co2, atm.o3[ilayer], gas.n2o,gas.co, \
            gas.ch4, gas.o2, wbroadl[ilayer]))

    fileID.write('%s' % '%%%%%')
    fileID.close()

def rrtmg_lw_htr(atm,nlev):
    os.system('! rm -rf TAPE5,TAPE6,OUTPUT_RRTM')
    os.system('/home/cls/projects/ctb-itan/cls/RRTMG/RRTMG_LW/run_examples_std_atm/rrtmg_lw') # replace with full path to compiled RRTMG_LW

def rrtmg_lw_output_read_band(nlev):
    fluxbroad = Dummy()
    fluxband  = Dummy()
    fluxbroad.fup = np.full(nlev,np.nan)
    fluxbroad.fdn = np.full(nlev,np.nan)
    fluxbroad.fnt = np.full(nlev,np.nan)
    fluxbroad.htr = np.full(nlev,np.nan)

    fi  = open('OUTPUT_RRTM','r')
    lines = fi.readlines()[3:]

    for j in range(0,nlev):
        tmp1 = lines[j]
        tmp2 = tmp1.split()
        for idx in range(len(tmp2)):
            if tmp2[idx][0]=='*':
                tmp2[idx]=np.nan
        fluxbroad.fup[j] = tmp2[2]
        fluxbroad.fdn[j] = tmp2[3]
        fluxbroad.fnt[j] = tmp2[4]
        fluxbroad.htr[j] = tmp2[5]

    fluxband.fup_band = np.full([16,nlev],np.nan)
    fluxband.fdn_band = np.full([16,nlev],np.nan)
    fluxband.fnt_band = np.full([16,nlev],np.nan)
    fluxband.htr_band = np.full([16,nlev],np.nan)

    for k in range(1,17):
        fi = open('OUTPUT_RRTM','r')
        lines = fi.readlines()[(nlev+4)*k+3:(nlev+4)*k+3+nlev]
        for j in range(0,nlev):
            tmp1 = lines[j]
            tmp2 = tmp1.split()
            for idx in range(len(tmp2)):
                if tmp2[idx][0]=='*':
                    tmp2[idx]=np.nan
            fluxband.fup_band[k-1,j] = tmp2[2]
            fluxband.fdn_band[k-1,j] = tmp2[3]
            fluxband.fnt_band[k-1,j] = tmp2[4]
            fluxband.htr_band[k-1,j] = tmp2[5]

    fi.close()
    return fluxbroad,fluxband,nlev

# end