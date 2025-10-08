# IMPORT LIBRARIES
import sys
import numpy as np
import xarray as xr
sys.path.insert(1,'/home/cls/projects/ctb-itan/cls/pylib/')
from cls_tools import make_DA

######################################################################
### HISTOGRAM BIN DATA ###############################################
######################################################################

bot_bounds = np.array([['0160', '1450', '2740'],['1450', '2740', '4030']])
hgt_bounds = np.array([['0538', '1499', '0509'],['1612', '3479', '0741']])
months     = np.arange(1,13)
slf        = np.array([ 0.15 , 0.525 , 0.85  ,  0.9745,   0.9995])
twp        = np.array([13.5  ,42.5   ,81.5   ,142.5   ,1215.    ])
dmv        = np.array([ 0.013, 0.0285, 0.0365,  0.052 ,   0.1035])
slf_bounds = np.array([0., 0.3  , 0.75 ,  0.95 ,  0.999,  1.    ])
twp_bounds = np.array([0.,27.   ,58.   ,105.   ,180.   ,2250.   ])
dmv_bounds = np.array([0., 0.026, 0.031,  0.042,  0.062,   0.145])

variables = ['lw_kernel_sfc','lw_kernel_toa','sw_kernel_sfc','sw_kernel_toa']

coords,dims = ['month','slf','twp','dmv'], {'month':months,'slf':slf,'twp':twp,'dmv':dmv}

for bidx in range(3):
    filename = 'BOT_'+bot_bounds[0,bidx]+'_'+bot_bounds[1,bidx]+'.nc'

    ds1 = xr.open_dataset('/home/cls/projects/ctb-itan/cls/KERNELS/RRTMG_wrapper/data/kernel_histograms/cloud_radiative_kernel_CBOT'+bot_bounds[0,bidx]+'_CHGT'+hgt_bounds[0,bidx]+'_v4_MAR2025.nc')
    ds2 = xr.open_dataset('/home/cls/projects/ctb-itan/cls/KERNELS/RRTMG_wrapper/data/kernel_histograms/cloud_radiative_kernel_CBOT'+bot_bounds[1,bidx]+'_CHGT'+hgt_bounds[1,bidx]+'_v4_MAR2025.nc')

    ds_new = {}
    for variable in variables:
        ds_new[variable] = make_DA(np.nanmean(np.stack((ds1[variable].values,ds2[variable].values),axis=0),axis=0),
                                    ds1[variable].long_name,
                                    ds1[variable].units,
                                    coords,
                                    dims,
                                    ds1[variable].description
                                    )
    ds = xr.Dataset(ds_new)
    ds.to_netcdf('/home/cls/projects/ctb-itan/cls/HISTOGRAMS/data/cloud_radiative_kernel_'+filename)

    ds.close()
    ds1.close()
    ds2.close()

ds_all = []
for bidx1 in range(3):
    for bidx2 in range(2):
        ds_all.append(xr.open_dataset('/home/cls/projects/ctb-itan/cls/KERNELS/RRTMG_wrapper/data/kernel_histograms/cloud_radiative_kernel_CBOT'+bot_bounds[bidx2,bidx1]+'_CHGT'+hgt_bounds[bidx2,bidx1]+'_v4_MAR2025.nc'))

ds_new = {}
for variable in variables:
    ds_new[variable] = make_DA(np.nanmean(np.stack((
                                            ds_all[0][variable].values,
                                            ds_all[1][variable].values,
                                            ds_all[2][variable].values,
                                            ds_all[3][variable].values,
                                            ds_all[4][variable].values,
                                            ds_all[5][variable].values),axis=0),axis=0),
                                ds_all[0][variable].long_name,
                                ds_all[0][variable].units,
                                coords,
                                dims,
                                ds_all[0][variable].description
                                )

ds_new = xr.Dataset(ds_new)
ds_new.to_netcdf('/home/cls/projects/ctb-itan/cls/HISTOGRAMS/data/cloud_radiative_kernel_BOT_all_altitudes.nc')
ds_new.close()

# END