import numpy as np
import sys
import pygrib
import os
import time as timey
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, interp
from dateutils import hrstodate, daterange, dayofyear, \
     splitdate, datetohrs, dateshift, dateto_hrs_since_day1CE
from netCDF4 import Dataset
import cPickle

# --- Read input data from command line.  Get beginning and ending lead time 
#     in hours and the name of the file with the dates we want

cleadb = sys.argv[1]
cleade = sys.argv[2]
cmean_or_control = sys.argv[3] # C or c for control, m or M for mean
cyyyymmddhh_end = sys.argv[4]
iyyyymmddhh_end = int(cyyyymmddhh_end)
rthresh = [0.1,1.0,2.5,5.0,10.0,15.0,20.0,25.0,30.0,35.0,40.0] # list of event thresholds
nthresh = len(rthresh)

# ---- read in ccpa conus mask

infile = '/mac_home/jwhitaker/python/NLDASmask_CONUS.nc'
cfile = Dataset(infile,'r')
iccpa_mask = cfile.variables['mask'][:,:]
iccpa_mask_t = np.transpose(iccpa_mask)
cfile.close()

# ---- open netcdf file with coincident precipitation analysis data and reforecast data
#      and read in list of dates that are in the file as well as lat/lon arrays for 
#      analysis and forecast grid.  Note that python does its arrays backwards from 
#      fortran, and clasically the data is stored as (y,x) not (x,y).

infilename = '/Rf2_tests/ccpa/netcdf/refcstv2_precip_ccpav3_'+cleadb+'_to_'+cleade+'.nc'
print infilename
rfile = Dataset(infilename, 'r')
yyyymmddhh_init  = rfile.variables["yyyymmddhh_init"][:]
yyyymmddhh_fcstb = rfile.variables["yyyymmddhh_fcstb"][:]
yyyymmddhh_fcste = rfile.variables["yyyymmddhh_fcste"][:]

print 'yyyymmddhh_init[0:4] = ',yyyymmddhh_init[0:4]
print 'yyyymmddhh_fcstb[0:4] = ',yyyymmddhh_fcstb[0:4]
print 'yyyymmddhh_fcste[0:4] = ',yyyymmddhh_fcste[0:4]

print np.shape(yyyymmddhh_fcste)
lendates = len(yyyymmddhh_fcste)
print 'lendates = ',lendates
lons_anal = rfile.variables["lons_anal"][:,:]
lats_anal = rfile.variables["lats_anal"][:,:]
lons_fcst = rfile.variables["lons_fcst"][:,:]
lats_fcst = rfile.variables["lats_fcst"][:,:]
rlons_1d = lons_fcst[0,:]
rlats_1d = lats_fcst[:,0]
lonsf, latsf = np.meshgrid(rlons_1d,rlats_1d)

nya, nxa = np.shape(lons_anal)
nyf, nxf = np.shape(lons_fcst)
# uncomment to do ETS calcs on the analysis grid
#ones  = np.ones((nya,nxa),dtype=np.float32)
#zeros = np.zeros((nya,nxa),dtype=np.float32)
ones  = np.ones((nyf,nxf),dtype=np.float32)
zeros = np.zeros((nyf,nxf),dtype=np.float32)


print 'iyyyymmddhh_end = ',iyyyymmddhh_end
print 'np.min,max (yyyymmddhh_fcste)= ', np.min(yyyymmddhh_fcste),np.max(yyyymmddhh_fcste)
itemindex = np.where(yyyymmddhh_fcste == iyyyymmddhh_end)
ids = int(itemindex[0])
apcp_anal_today = rfile.variables['apcp_anal'][ids,:,:]
apcp_ensfcst_today = rfile.variables['apcp_fcst_ens'][ids,:,:,:]

if cmean_or_control == 'c' or cmean_or_control == 'C':
    apcp_today = apcp_ensfcst_today[0,:,:] # control member
else:
    apcp_today = np.mean(apcp_ensfcst_today, axis=0) # ensemble mean
  
# ----- make plots of precipitation analyses and forecasts

colorst = ['White','BurlyWood','LightGreen','Green','Gold',\
           'Orange','Red','Plum','LightSkyBlue','LightGray','Gray','Brown']
colorstblack=['White','Black','Black','Black','Black','Black','Black','Black','Black','Black','Black']
colorstwhite=['Black','White','White','White','White','White','White','White','White','White','White']

fig1 = plt.figure(figsize=(6.5,9))

# ---- plot 1: analysis data

clevs = [0.0, 0.1, 0.5, 1, 2.5, 5, 10, 25, 50, 100, 150, 200]
ax = fig1.add_axes([0.05,.6,0.9,0.35])
ctitle = '(a) 24-h accumulated precipitation analysis\nending '+cyyyymmddhh_end 
ax.set_title(ctitle,fontsize=13)
#m = Basemap(projection='mill',llcrnrlon=-95.,llcrnrlat=25.,\
#        urcrnrlon=-75,urcrnrlat=39.,resolution='l')  # resolution next step = i
m = Basemap(projection='mill',llcrnrlon=-125.,llcrnrlat=25.,\
        urcrnrlon=-65,urcrnrlat=50.,resolution='l')  # resolution next step = i
x,y = m(lons_anal,lats_anal)
CS1 = m.contour(x,y,apcp_anal_today,clevs,colors=colorstblack,cmap=None,linewidths=0.3)
CS2 = m.contourf(x,y,apcp_anal_today,clevs,colors=colorst,cmap=None,extend='neither')
m.drawcoastlines(linewidth=.5)
m.drawstates(linewidth=.5)
m.drawcountries(linewidth=.5)

cax = fig1.add_axes([0.2,0.56,0.6,0.02])
cbar = fig1.colorbar(CS2,extend='neither', \
   orientation='horizontal',cax=cax,drawedges=True,ticks=clevs,format='%g')
cax.set_xlabel('Analyzed precipitation amount (mm)',fontsize=10)

# ---- plot 2: forecast data

clevs = [0.0, 0.1, 0.5, 1, 2.5, 5, 10, 25, 50, 100, 150, 200]
ax = fig1.add_axes([0.05,.1,0.9,0.35])
ctitle = '(b) '+cleadb+'-h to '+cleade+'-h accumulated precipitation\n forecast ending '+cyyyymmddhh_end
ax.set_title(ctitle,fontsize=13)
#m = Basemap(projection='mill',llcrnrlon=-95.,llcrnrlat=25.,\
#        urcrnrlon=-75,urcrnrlat=39.,resolution='l')  # resolution next step = i
m = Basemap(projection='mill',llcrnrlon=-125.,llcrnrlat=25.,\
        urcrnrlon=-65,urcrnrlat=50.,resolution='l')  # resolution next step = i
x,y = m(lonsf,latsf)
CS1 = m.contour(x,y,apcp_today,clevs,colors=colorstblack,cmap=None,linewidths=0.3)
CS2 = m.contourf(x,y,apcp_today,clevs,colors=colorst,cmap=None,extend='neither')
m.drawcoastlines(linewidth=.5)
m.drawstates(linewidth=.5)
m.drawcountries(linewidth=.5)

cax = fig1.add_axes([0.2,0.06,0.6,0.02])
cbar = fig1.colorbar(CS2,extend='neither', \
   orientation='horizontal',cax=cax,drawedges=True,ticks=clevs,format='%g')
cax.set_xlabel('Forecast precipitation amount (mm)',fontsize=10)

# ---- set plot title

plot_title = 'refcstv2_precip_'+cleadb+'_to_'+cleade+'_endtime='+cyyyymmddhh_end+'.pdf'
fig1.savefig(plot_title)
print 'saving plot to file = ',plot_title
#plt.show()


