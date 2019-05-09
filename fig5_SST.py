from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker
from matplotlib.mlab import bivariate_normal
import math

FILEDIR = '/n/home12/hongwei/TC/NCEP_data/'

geos_nc = Dataset(FILEDIR+'sst.mnmean.nc','r',format='NETCDF4_CLASSIC')

SST = geos_nc.variables['sst']

iyear = np.linspace(1988, 2017, num=30)

# small: 1995, 1996, 1999, 2004, 2007, 2010, 2011
# big: 1990, 1992, 2006, 2014, 2015

SST_small = SST[0:7,:,:]
SST_big = SST[0:5,:,:]
print(SST_small.shape)

SST_small[0] = np.mean(SST[162:166,:,:], axis=0)
SST_small[1] = np.mean(SST[174:178,:,:], axis=0)
SST_small[2] = np.mean(SST[210:214,:,:], axis=0)
SST_small[3] = np.mean(SST[270:274,:,:], axis=0)
SST_small[4] = np.mean(SST[306:310,:,:], axis=0)
SST_small[5] = np.mean(SST[342:346,:,:], axis=0)
SST_small[6] = np.mean(SST[354:358,:,:], axis=0)

small_mean = np.mean(SST_small, axis=0)

SST_big[0] = np.mean(SST[102:106,:,:], axis=0)
SST_big[1] = np.mean(SST[126:130,:,:], axis=0)
SST_big[2] = np.mean(SST[294:394,:,:], axis=0)
SST_big[3] = np.mean(SST[390:394,:,:], axis=0)
SST_big[4] = np.mean(SST[402:406,:,:], axis=0)

big_mean = np.mean(SST_big, axis=0)

a_SST = big_mean - small_mean

del geos_nc

# t test:
# x - small
# y - big
small_num	= 7
big_num		= 5
small_std	= np.std(SST_small, axis=0)
big_std		= np.std(SST_big, axis=0)

t_alpha = 2.228
aSST_alpha = np.sqrt(np.power(big_std,2)/big_num + np.power(small_std,2)/small_num ) * t_alpha

# n = 10;  t_alpha: 1.812(0.1), 2.228(0.05)
t = np.divide((big_mean - small_mean),np.sqrt(np.power(big_std,2)/big_num + np.power(small_std,2)/small_num ))
print(t.shape)

#------------------------------------------------
# plot  -----------------------------------------
#------------------------------------------------

plt.figure(figsize=(7,9))

m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=0,urcrnrlon=360)
m.drawcoastlines()
m.drawparallels(np.arange(-90,91.,30.))
m.drawmeridians(np.arange(0.0,361.,60.))
m.drawmapboundary(fill_color='white')

parallels = np.arange(90,90,30.)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=8)
meridians = np.arange(0.0,361.,60.)
m.drawmeridians(meridians,labels=[1,0,0,0],fontsize=8)
#m.drawcountries()

ny = SST.shape[1]; nx = SST.shape[2]
print(ny)
lons, lats = m.makegrid(nx, ny) # get lat/lons of ny by nx evenly space grid.
x, y = m(lons, lats) # compute map proj coordinates.

# uneven bounds changes the colormapping:
bounds = np.array([-1.0,-0.7,-0.4,-0.1,-0.01,0,0.01,0.1,0.4,0.7,1.0])
#bounds = np.array([20.0,22.5,25.0,27.5,30.0,32.5])
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)

cs = m.pcolormesh(x, y, a_SST[::-1,:], norm=norm, cmap='bwr')
cs.cmap.set_under('w')
cs.set_clim(-10)

levs = [-2.228,2.228]
plt.contour(x, y, t[::-1,:], levels=levs, linewidths=2, colors='y')

# add colorbar.
fmt = matplotlib.ticker.ScalarFormatter(useMathText=True)
fmt.set_powerlimits((0, 0))

cbar = m.colorbar(cs,location='bottom',pad="9%",format=fmt)
cbar.set_label('degC')

plt.title('SST', fontsize=10)

plt.savefig('SST.png')
plt.clf()
plt.cla()


