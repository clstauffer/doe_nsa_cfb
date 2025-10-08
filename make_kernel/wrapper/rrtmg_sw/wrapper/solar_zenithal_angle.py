# import libraries
from math import * 
import numpy as np 

def subsolar(utc):
    ye, mo, da, ho, mi, se = utc
    ta = pi * 2
    ut = ho + mi / 60 + se / 3600
    t = 367 * ye - 7 * (ye + (mo + 9) // 12) // 4
    dn = t + 275 * mo // 9 + da - 730531.5 + ut / 24
    sl = dn * 0.01720279239 + 4.894967873
    sa = dn * 0.01720197034 + 6.240040768
    t = sl + 0.03342305518 * sin(sa)
    ec = t + 0.0003490658504 * sin(2 * sa)
    ob = 0.4090877234 - 0.000000006981317008 * dn
    st = 4.894961213 + 6.300388099 * dn
    ra = atan2(cos(ob) * sin(ec), cos(ec))
    de = asin(sin(ob) * sin(ec))
    la = degrees(de)
    lo = degrees(ra - st) % 360
    lo = lo - 360 if lo > 180 else lo
    return [round(la, 6), round(lo, 6)]

def theta_min_max(lat_subsolar, lat_point):
    return [abs(lat_point - lat_subsolar),180-(abs(lat_point - lat_subsolar))]

def range_theta(month_ref, lat_point): 
    date_ref = (2008,month_ref, 15,0,0,0)
    lat_s, lon_s = subsolar(date_ref)
    theta_min, theta_max = theta_min_max(lat_s, lat_point)
    if theta_min<0: 
        theta_min=0
    if theta_max>90:
        theta_max = 90 
    return np.arange(theta_min, theta_max)

# end