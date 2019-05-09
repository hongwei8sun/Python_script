from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker
from matplotlib.mlab import bivariate_normal
import math

FILEDIR = '/n/home12/hongwei/TC/NCEP_data/'

geos_nc = Dataset(FILEDIR+'mslp.mon.mean.nc','r',format='NETCDF4_CLASSIC')

SLP = geos_nc.variables['mslp']

iyear = np.linspace(1988, 2017, num=30)

# small: 1995, 1996, 1999, 2004, 2007, 2010, 2011
# big: 1990, 1992, 2006, 2014, 2015

SLP_small = SLP[0:7,:,:]
SLP_big = SLP[0:5,:,:]
print(SLP_small.shape)

SLP_small[0] = np.mean(SLP[198:202,:,:], axis=0)
SLP_small[1] = np.mean(SLP[210:214,:,:], axis=0)
SLP_small[2] = np.mean(SLP[246:250,:,:], axis=0)
SLP_small[3] = np.mean(SLP[306:310,:,:], axis=0)
SLP_small[4] = np.mean(SLP[342:346,:,:], axis=0)
SLP_small[5] = np.mean(SLP[378:382,:,:], axis=0)
SLP_small[6] = np.mean(SLP[390:394,:,:], axis=0)

small_mean = np.mean(SLP_small, axis=0)

SLP_big[0] = np.mean(SLP[138:142,:,:], axis=0)
SLP_big[1] = np.mean(SLP[162:166,:,:], axis=0)
SLP_big[2] = np.mean(SLP[330:334,:,:], axis=0)
SLP_big[3] = np.mean(SLP[426:430,:,:], axis=0)
SLP_big[4] = np.mean(SLP[438:442,:,:], axis=0)

big_mean = np.mean(SLP_big, axis=0)

a_SLP = big_mean - small_mean

del geos_nc

# t test:
# x - small
# y - big
small_num	= 7
big_num		= 5
small_std	= np.std(SLP_small, axis=0)
big_std		= np.std(SLP_big, axis=0)

t_alpha = 2.228
aSLP_alpha = np.sqrt(np.power(big_std,2)/big_num + np.power(small_std,2)/small_num ) * t_alpha

# n = 10;  t_alpha: 1.812(0.1), 2.228(0.05)
t = np.divide( a_SLP , np.sqrt(np.power(big_std,2)/big_num + np.power(small_std,2)/small_num ) )
print(t.shape)

#------------------------------------------------
# plot  -----------------------------------------
#------------------------------------------------

plt.figure(figsize=(7,9))

m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=0,urcrnrlon=360)
m.drawcoastlines()
m.drawparallels(np.arange(-90.,91.,30.))
m.drawmeridians(np.arange(-0.01,361.,60.))
m.drawmapboundary(fill_color='white')

parallels = np.arange(-90.,90,30.)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=8)
meridians = np.arange(-0.01,361.,60.)
m.drawmeridians(meridians,labels=[1,0,0,0],fontsize=8)
#m.drawcountries()

ny = SLP.shape[1]; nx = SLP.shape[2]
print(ny)
lons, lats = m.makegrid(nx, ny) # get lat/lons of ny by nx evenly space grid.
x, y = m(lons, lats) # compute map proj coordinates.

# uneven bounds changes the colormapping:
# bounds = np.array([100000,100750,101000,102000])
bounds = np.array([-700,-400,-100,-10,0,10,100,400,700])
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)

# plot the variable:
cs = m.pcolormesh(x, y, a_SLP[::-1,:], norm=norm, cmap='bwr')  # Blues_r') #bwr
cs.cmap.set_under('w')
# cs.set_clim(-10000000)

levs = [-2.228,2.228]
plt.contour(x, y, t[::-1,:], levels=levs, linewidths=2, colors='y')

# add colorbar.
fmt = matplotlib.ticker.ScalarFormatter(useMathText=True)
fmt.set_powerlimits((0, 0))

cbar = m.colorbar(cs,location='bottom',pad="9%",format=fmt)
cbar.set_label('degC')

plt.title('SLP', fontsize=10)

plt.savefig('SLP.png')
plt.clf()
plt.cla()


