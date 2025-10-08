import numpy as np 
import matplotlib.pyplot as plt








def Plt_profile(fluxbroad_clr, fluxband_clr, fluxbroad, fluxband, lev, name):
    plt.clf()
    plt.figure(1)

    lev = lev[::-1]

#    print(lev)
    plt.title('longwave clear sky')
    plt.plot(np.sum(fluxband_clr.fup_band, axis=0),lev, '-o', color = 'r', label= 'Longwave up')
    plt.plot(np.sum(fluxband_clr.fdn_band, axis=0),lev, '-o', color = 'Orange', label= 'Longwave down')
    plt.xlabel(r'flux (W m$^{-2}$)')
    plt.legend(loc='best')
    plt.gca().invert_yaxis()
    plt.savefig('./Graph/vert_updn_%s.png'%(name))

    plt.clf()
    plt.figure(2)
    plt.title('longwave')
    plt.plot(np.sum(fluxband_clr.fdn_band, axis=0), lev, '-o', color = 'r', label= 'clear sky')
    plt.plot(np.sum(fluxband.fdn_band, axis=0), lev,'-o', color = 'k', label= 'cloudy sky')
    plt.xlabel(r'flux (W m$^{-2}$)')
    plt.legend(loc='best')
    plt.gca().invert_yaxis()
    plt.savefig('./Graph/vert_clr_cld_%s.png'%(name))

    plt.clf()
    plt.figure(3)
    plt.title('longwave downward at the surface')
    LW_bands = [350, 500, 630, 700, 820, 980, 1080, 1180, 1390, 1480, 1800, 2080, 2250, 2380, 2600, 3250]
#    print(fluxband_clr.fdn_band[:,:])
    plt.plot(LW_bands, fluxband_clr.fdn_band[:,-1], '-o', color = 'r', label= 'clear sky')
    plt.plot(LW_bands, fluxband.fdn_band[:,-1], '-o', color = 'k', label= 'cloudy sky')
    plt.ylabel(r'flux (W m$^{-2}$)')
    plt.xlabel(r'Longwave (cm$^{-1}$)')
    plt.legend(loc='best')
    plt.savefig('./Graph/longwave_clr_cld_%s.png'%(name))


    





















