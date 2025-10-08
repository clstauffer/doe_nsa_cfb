#####################################
#   for information the level for clouds in hPa:
#[1014.3088  987.5     962.5     937.5     912.5     887.5     862.5    837.5     812.5     787.5     762.5     725.      675.      625.    575.      525.      475.      425.      375.      325.      275.    237.5     212.5     187.5     162.5     137.5     112.5      85.    60.       40.       25.       15.    ]
#####################################################

# import libraries
import os
import sys
import time
import argparse
import traceback
import numpy as np
import netCDF4 as nc
# import functions
sys.path.insert(1,'/home/cls/projects/ctb-itan/cls/kernels/rrtmg_lw/wrapper/')
import run_rrtmg_cloud_allsky


t0 = time.time()

try: ncfile.close()  # just to be safe, make sure dataset is not already open.
except: pass
global month, cloud_base, cloud_top

data_dir = '/home/cls/projects/ctb-itan/cls/kernels/input_files/'
# args = ['armb',117,268,2533,442,359,1040]
# args = ['ceil',133,379,2624,561,511,462]
args = ['extr',117,268,2533,561,511,1040]

# Cloud_Base_bounds   = np.array([int(str(sys.argv[1])),int(str(sys.argv[2])),int(str(sys.argv[3])),int(str(sys.argv[4]))])
# Cloud_Height_bounds = np.array([int(str(sys.argv[5])),int(str(sys.argv[6])),int(str(sys.argv[7])),int(str(sys.argv[8]))])

Cloud_Base   = np.array([args[1],args[2],args[3]]) # (Cloud_Base_bounds[:-1]+Cloud_Base_bounds[1:])/2.
Cloud_Height = np.array([args[4],args[5],args[6]]) # (Cloud_Height_bounds[:-1]+Cloud_Height_bounds[1:])/2.
Cloud_Top    = Cloud_Base + Cloud_Height

# Cloud_Base_bounds      = np.array([35., 232., 527., 4721.])#[180, 254, 356, 1262]
# Cloud_Base      = (Cloud_Base_bounds[:-1]+Cloud_Base_bounds[1:])/2.
# Cloud_Height    = np.array([561,513,463])
# Cloud_Top = Cloud_Base + Cloud_Height

dname = args[0]
cld_percentiles = ['033','066','100']

lw_bins  = np.array([350, 500, 630, 700, 820, 980, 1080, 1180, 1390, 1480, 1800, 2080, 2250, 2380, 2600, 3250])
slf_bins = np.arange(0,1.1,0.2)
re_bins  = [ 1., 11., 23., 32., 39., 55.]
wp_bins  = [0.004, 0.7, 2., 15., 119., 1027.]

for month in np.arange(1,13): # loop over months

    for cp,cloud_base,cloud_top in zip(cld_percentiles,Cloud_Base,Cloud_Top): # loop over cloud layers
        start = time.process_time()

        nc_profile = nc.Dataset(data_dir+'z_pres_profiles/z_pres_profile_%s.nc'%(int(month)))
        altitude   = nc_profile['alt'][:]

        cloud_base_arg = int(np.argmin(abs(altitude-cloud_base)))
        cloud_top_arg  = int(np.argmin(abs(altitude-cloud_top)))#+1

        print(dname,month,cp,altitude[cloud_base_arg],altitude[cloud_top_arg])
