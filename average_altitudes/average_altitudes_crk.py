import xarray as xr

ds = []
files = ['cloud_radiative_kernel_CBOT0160_CHGT0538_v2.nc','cloud_radiative_kernel_CBOT1450_CHGT1499_v2.nc','cloud_radiative_kernel_CBOT1450_CHGT1612_v2.nc','cloud_radiative_kernel_CBOT2740_CHGT0509_v2.nc','cloud_radiative_kernel_CBOT2740_CHGT3479_v2.nc','cloud_radiative_kernel_CBOT4030_CHGT0741_v2.nc']
for file in files:
    ds.append(xr.open_dataset(file))
ds = xr.concat(ds, dim='ensemble')  # new dimension: 'ensemble'
ds = ds.mean(dim='ensemble',skipna=True)
ds.to_netcdf('cloud_radiative_kernel_all_altitudes_v2.nc')
ds.close()
