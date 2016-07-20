"""this script operates on the output of the analog method ppn_analog_ccpa_supplocns.x 
   and provides plots in this case for probabilities of exceeding the input threshold
   entered in millimeters.  User must also input:
   
   (a) the initial date/time of the forecast of interest.
   (b) the beginning forecast lead time range (e.g., 024).
   (c) the end forecast lead time (e.g., 036)
   (d) the threshold (POP,1mm,2p5mm,5mm,10mm,25mm,50mm,q50,q67,q80,q90,q95)
   (e) the number of supplemental locations (01 to 20)
   
   Running this program is predicated on you having previously run the 
   fortran routine ppn_analog_ccpa_supplocns.x for this particular lead time, 
   and for this particular number of supplemental locations.
   
   Coded by: Tom Hamill, tom.hamill@noaa.gov, Mar 2016
   
"""
# --- import the following python modules

import numpy as np
import numpy.ma as ma
import sys, os
from dateutils import hrstodate, daterange, dayofyear, splitdate, \
    datetohrs, dateshift, dateto_hrs_since_day1CE  # homegrown code
import matplotlib
matplotlib.use('Agg')

from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
from matplotlib import rcParams
from netCDF4 import Dataset
from mpl_toolkits.basemap import interp, Basemap, shiftgrid, addcyclic

# --- import f2py fortran routines 

from read_analog_probabilities import read_analog_probabilities
from get_quantiles_fromfile import get_quantiles_fromfile

# ---- internal functions

def average_analysis_to_fcstgrid(nya, nxa, latsa, lonsa, nyf, nxf, \
    latsf, lonsf, thresh_mediana, thresh_uppertercilea, \
    thresh_upperquintilea, thresh_upperdecilea, thresh_95a):

    thresh_medianf = np.zeros((nyf,nxf), dtype=np.float32)
    thresh_uppertercilef = np.zeros((nyf,nxf), dtype=np.float32)
    thresh_upperquintilef = np.zeros((nyf,nxf), dtype=np.float32)
    thresh_upperdecilef = np.zeros((nyf,nxf), dtype=np.float32)
    thresh_95f = np.zeros((nyf,nxf), dtype=np.float32)

    nearest_lon = np.zeros(nxa,dtype=np.float32)
    nearest_lat = np.zeros(nya,dtype=np.float32)
    
    numer_medianf = np.zeros((nyf,nxf), dtype=np.float32)
    denom_medianf = np.zeros((nyf,nxf), dtype=np.float32)
    numer_uppertercilef = np.zeros((nyf,nxf), dtype=np.float32)
    denom_uppertercilef = np.zeros((nyf,nxf), dtype=np.float32)
    numer_upperquintilef = np.zeros((nyf,nxf), dtype=np.float32)
    denom_upperquintilef = np.zeros((nyf,nxf), dtype=np.float32)
    numer_upperdecilef = np.zeros((nyf,nxf), dtype=np.float32)
    denom_upperdecilef = np.zeros((nyf,nxf), dtype=np.float32)
    numer_95f = np.zeros((nyf,nxf), dtype=np.float32)
    denom_95f = np.zeros((nyf,nxf), dtype=np.float32)
    
    # ---- determine indices of nearest forecast lat/lon grid point for every analysis

    for ixa in range(nxa):
        dmin = 99999.
        ilon = -99
        for ixf in range(nxf):
            dist = np.abs(lonsf[0,ixf] - lonsa[0,ixa])
            if dist < dmin:
                dmin = dist
                ilon = ixf
        nearest_lon[ixa] = ilon
            
    for jya in range(nya):
        dmin = 99999.
        ilat = -99
        for jyf in range(nyf):
            dist = np.abs(latsf[jyf,0] - latsa[jya,0])
            if dist < dmin:
                dmin = dist
                ilat = jyf
        nearest_lat[jya] = ilat           

    # ---- tally information needed to calculate average

    for jya in range (nya):
        jyf = nearest_lat[jya]
        for ixa in range(nxa):
            ixf = nearest_lon[ixa]
            if thresh_mediana[jya,ixa] >= 0.:
                numer_medianf[jyf,ixf] = numer_medianf[jyf,ixf] + thresh_mediana[jya,ixa]
                denom_medianf[jyf,ixf] = denom_medianf[jyf,ixf] + 1.0
            if thresh_uppertercilea[jya,ixa] >= 0.:
                numer_uppertercilef[jyf,ixf] = numer_uppertercilef[jyf,ixf] + thresh_uppertercilea[jya,ixa]
                denom_uppertercilef[jyf,ixf] = denom_uppertercilef[jyf,ixf] + 1.0
            if thresh_upperquintilea[jya,ixa] >= 0.:
                numer_upperquintilef[jyf,ixf] = numer_upperquintilef[jyf,ixf] + thresh_upperquintilea[jya,ixa]
                denom_upperquintilef[jyf,ixf] = denom_upperquintilef[jyf,ixf] + 1.0
            if thresh_upperdecilea[jya,ixa] >= 0.:
                numer_upperdecilef[jyf,ixf] = numer_upperdecilef[jyf,ixf] + thresh_upperdecilea[jya,ixa]
                denom_upperdecilef[jyf,ixf] = denom_upperdecilef[jyf,ixf] + 1.0
            if thresh_95a[jya,ixa] >= 0.:
                numer_95f[jyf,ixf] = numer_95f[jyf,ixf] + thresh_95a[jya,ixa]
                denom_95f[jyf,ixf] = denom_95f[jyf,ixf] + 1.0

    # ---- calculate average for forecast grid box if data available there.

    for ixf in range(nxf):
        for jyf in range(nyf):
            if denom_medianf[jyf,ixf] > 0.:
                thresh_medianf[jyf,ixf] = numer_medianf[jyf,ixf] / denom_medianf[jyf,ixf]
            else:
                thresh_medianf[jyf,ixf] = -99.99
            if denom_uppertercilef[jyf,ixf] > 0.:
                thresh_uppertercilef[jyf,ixf] = numer_uppertercilef[jyf,ixf] / denom_uppertercilef[jyf,ixf]
            else:
                thresh_uppertercilef[jyf,ixf] = -99.99
            if denom_upperquintilef[jyf,ixf] > 0.:
                thresh_upperquintilef[jyf,ixf] = numer_upperquintilef[jyf,ixf] / denom_upperquintilef[jyf,ixf]
            else:
                thresh_upperquintilef[jyf,ixf] = -99.99
            if denom_upperdecilef[jyf,ixf] > 0.:
                thresh_upperdecilef[jyf,ixf] = numer_upperdecilef[jyf,ixf] / denom_upperdecilef[jyf,ixf]
            else:
                thresh_upperdecilef[jyf,ixf] = -99.99
            if denom_95f[jyf,ixf] > 0.:
                thresh_95f[jyf,ixf] = numer_95f[jyf,ixf] / denom_95f[jyf,ixf]
            else:
                thresh_95f[jyf,ixf] = -99.99
  
    return thresh_medianf, thresh_uppertercilef, thresh_upperquintilef, thresh_upperdecilef, thresh_95f


# --- set some display properties

rcParams['legend.fontsize']='small'
rcParams['xtick.labelsize']='x-small' # for colorbar tick labels
rcParams['axes.labelsize']='small'  # for colorbar wording
cmonths = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

# --- from command line get the initial date and forecast lead time in hours

cyyyymmddhh_init = sys.argv[1] # YYYYMMDDHH format of initial forecast date
iyyyymmddhh_init = int(cyyyymmddhh_init)
cleadb = sys.argv[2]  # enter beginning lead time in hours
cleade = sys.argv[3]  # enter endlead time in hours
cthresh = sys.argv[4] # threshold values must take one from below:
    # 'POP', '1mm', '2p5mm', '5mm', '10mm', '25mm', '50mm', 'q50', 'q67', 'q80', 'q90', 'q95'
csupp = sys.argv[5] # number of supplemental locations, 01 to 20

ilead = int(cleade)
cyyyymmddhh_fcst = dateshift(cyyyymmddhh_init,ilead)
iday = ilead/24
clead_dayb = str(int(cleadb)/24)
clead_daye = str(iday)
cmonth = cyyyymmddhh_init[4:6]
imonth = int(cmonth)
cmonthname = cmonths[imonth-1]

# --- set up inputs and constants

input_data_directory = '/Projects/Reforecast2/netcdf' # where input CCPA and reforecast data stored
analog_data_directory = '/Rf2_tests/ccpa'  # where analog forecast output was stored
nxa = 464  # ccpa grid dimensions
nya = 224
nxf = 128  # reforecast grid surrounding CONUS
nyf = 61
nmembers = 11

# ---- read in the raw ensemble from netcdf file and calculate raw event probabilities

infile = input_data_directory+'/refcstv2_precip_ccpav3_'+cleadb+'_to_'+cleade+'.nc'
print infile
nc = Dataset(infile,'r')
yyyymmddhh_init = nc.variables["yyyymmddhh_init"][:]
yyyymmddhh_fcst = nc.variables["yyyymmddhh_fcste"][:]
conusmask = nc.variables["conusmask"][:]
latsf = nc.variables['lats_fcst'][:]
lonsf = nc.variables['lons_fcst'][:]
latsa = nc.variables['lats_anal'][:]
lonsa = nc.variables['lons_anal'][:]
itemindex = np.where(yyyymmddhh_init == iyyyymmddhh_init)
id = int(itemindex[0])

print '******* id = ',id
apcp_anal = nc.variables["apcp_anal"][id,:,:]
apcp_fcst_ens = nc.variables["apcp_fcst_ens"][id,:,:,:]
apcp_fcst_ens = np.squeeze(apcp_fcst_ens)
apcp_fcst_mean = np.mean(apcp_fcst_ens,axis=0)
nc.close()
print 'np.shape(apcp_fcst_ens) = ',np.shape(apcp_fcst_ens) 

# ---- in some circumstances we may be plottting the probability of exceeding 
#      a quantile of the observed distribution.  To prepare for this contingency
#      read in the threshold precipitation amounts associated with these quantiles

infile_quantiles = analog_data_directory + '/quantiles_'+cmonths[imonth-1]+'_'+\
	 cleadb+'_to_'+cleade+'.dat'
print infile_quantiles
thresh_median_t = np.zeros((nxa,nya),dtype=np.float32)
thresh_uppertercile_t = np.zeros((nxa,nya),dtype=np.float32)
thresh_upperquintile_t = np.zeros((nxa,nya),dtype=np.float32)
thresh_upperdecile_t = np.zeros((nxa,nya),dtype=np.float32)
thresh_95_t = np.zeros((nxa,nya),dtype=np.float32)
print get_quantiles_fromfile.__doc__
print infile_quantiles, nxa, nya
print np.shape(thresh_upperquintile_t)
thresh_median_t, thresh_uppertercile_t, \
    thresh_upperquintile_t,thresh_upperdecile_t, thresh_95_t=\
    get_quantiles_fromfile(infile_quantiles, nxa, nya)    
    
thresh_mediana = np.transpose(thresh_median_t)
thresh_uppertercilea= np.transpose(thresh_uppertercile_t)
thresh_upperquintilea = np.transpose(thresh_upperquintile_t)
thresh_upperdecilea = np.transpose(thresh_upperdecile_t)
thresh_95a = np.transpose(thresh_95_t)

thresh_medianf, thresh_uppertercilef, thresh_upperquintilef, \
    thresh_upperdecilef, thresh_95f = average_analysis_to_fcstgrid(nya, \
    nxa, latsa, lonsa, nyf, nxf, latsf, lonsf, thresh_mediana, \
    thresh_uppertercilea, thresh_upperquintilea, thresh_upperdecilea, thresh_95a)

rthreshf = np.zeros((nmembers,nyf, nxf), dtype=np.float32)
rthresha = np.zeros((nya, nxa), dtype=np.float32)
# 'POP', '1mm', '2p5mm', '5mm', '10mm', '25mm', '50mm', 'q50', 'q67', 'q80', 'q90', 'q95'
    
if cthresh == 'POP':
    for i in range(nmembers):
        rthreshf[i,:,:] = 0.4
    rthresha[:,:] = 0.4
elif cthresh == '1mm':
    for i in range(nmembers):
        rthreshf[i,:,:] = 1.0
    rthresha[:,:] = 1.0
elif cthresh == '2p5mm':
    for i in range(nmembers):
        rthreshf[i,:,:] = 2.5
    rthresha[:,:] = 2.5
elif cthresh == '5mm':
    for i in range(nmembers):
        rthreshf[i,:,:] = 5.0
    rthresha[:,:] = 5.0
elif cthresh == '10mm':
    for i in range(nmembers):
        rthreshf[i,:,:] = 10.0
    rthresha[:,:] = 10.0
elif cthresh == '25mm':
    for i in range(nmembers):
        rthreshf[i,:,:] = 25.0
    rthresha[:,:] = 25.0
elif cthresh == '50mm':
    for i in range(nmembers):
        rthreshf[i,:,:] = 50.0
    rthresha[:,:] = 50.0
elif cthresh == 'q50':
    for i in range(nmembers):
        rthreshf[i,:,:] = thresh_medianf[:,:]
    rthresha[:,:] = thresh_mediana[:,:]
elif cthresh == 'q67':
    for i in range(nmembers):
        rthreshf[i,:,:] = thresh_uppertercilef[:,:]
    rthresha[:,:] = thresh_uppertercilea[:,:]
elif cthresh == 'q80':
    for i in range(nmembers):
        rthreshf[i,:,:] = thresh_upperquintilef[:,:]
    rthresha[:,:] = thresh_upperquintilea[:,:]
elif cthresh == 'q90':
    for i in range(nmembers):
        rthreshf[i,:,:] = thresh_upperdecilef[:,:]
    rthresha[:,:] = thresh_upperdecilea[:,:]
elif cthresh == 'q95':
    for i in range(nmembers):
        rthreshf[i,:,:] = thresh_95f[:,:]
    rthresha[:,:] = thresh_95a[:,:]
else:
    print 'inappropriate threshold entered.  Exiting'
    sys.exit()


prob_raw = np.zeros((nyf, nxf),dtype=np.float32)
oneorzero = np.zeros((nmembers,nyf,nxf),dtype=np.float32)
ones = np.ones((nmembers,nyf,nxf),dtype=np.float32)
zeros = np.zeros((nmembers,nyf,nxf),dtype=np.float32)

oneorzero = np.where(apcp_fcst_ens >= rthreshf, ones, zeros)
sumonezero = np.sum(oneorzero,axis=0)
prob_raw = sumonezero[:,:]/np.real(nmembers)
missing = -99.99*np.ones((nyf,nxf),dtype=np.float32)
prob_raw = np.where(rthreshf[0,:,:] >=0, prob_raw, missing)

# ---- read in post-processed info, previously created by ppn_analog_ccpa_realtime_supplocns.x
#      as well as vector of date information.

iccpa_mask_t = np.zeros((nxa,nya),dtype=np.int16)
prob_t  = np.zeros((nxa,nya),dtype=np.float32)

infile_dates = analog_data_directory + '/ppn_pwat_ccpa_supp' + csupp + '_analog_dates_' + \
    cmonthname + '_' + cleadb + '_to_' + cleade + '.dat'
infile_data =  analog_data_directory + '/ppn_pwat_ccpa_supp' + csupp + '_analog_data_' + \
    cmonthname + '_' + cleadb + '_to_' + cleade + '.dat'

date_anal_in = int(cyyyymmddhh_fcst)
probfcst_t, iccpa_mask_t = read_analog_probabilities(infile_dates, infile_data, \
    date_anal_in, nxa, nya, nxf, nyf, cthresh)

probfcst = np.transpose(probfcst_t)
print 'max,min probfcst',np.max(probfcst)
iccpa_mask = np.transpose(iccpa_mask_t)

# ---- define color tables for plotting

colors_ahps = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8','#A6ECA6','#42F742','Yellow','Gold',\
  'Orange','#FCD5D9','#F6A3AE','#FA5257','#E1C4FF','#AD8ADB','#A449FF']

colorstblack=['White','Black','Black','Black','Black', 'Black','Black','Black',\
  'Black','Black','Black','Black','Black','Black','Black','Black']

fig1 = plt.figure(figsize=(9.,6.5))
clead_dayb = str(int(cleadb)/24)
clead_daye = str(iday)

# ---- plot 1: raw forecast ensemble mean

clevs = [0.0, 0.1, 0.5, 1, 2.5, 5, 10, 15, 20, 25, 50, 75, 100, 125, 150]
ax = fig1.add_axes([0.02,.1,0.96,0.85])
ctitle = '(a) '+cleadb+'-'+cleade+'-h raw mean forecast, initialized '+\
    cyyyymmddhh_init+', valid '+cyyyymmddhh_fcst
ax.set_title(ctitle,fontsize=15)
m = Basemap(projection='mill',llcrnrlon=lonsa[0,0],llcrnrlat=latsa[0,0],\
        urcrnrlon=lonsa[-1,-1],urcrnrlat=latsa[-1,-1],resolution='l')
x,y = m(lonsf,latsf)
#ppnfcst_today_m = ma.array(apcp_fcst_mean_realtime ,mask=1-np.int(iccpa_mask))
CS2 = m.contourf(x,y,apcp_fcst_mean ,clevs,colors=colors_ahps,cmap=None,extend='neither')
CS1 = m.contour(x,y,apcp_fcst_mean,clevs,colors=colorstblack,cmap=None,linewidths=0.25)
m.drawcoastlines(linewidth=.5)
m.drawstates(linewidth=.5)
m.drawcountries(linewidth=.5)

cax = fig1.add_axes([0.05,0.09,0.9,0.04])
print 'clevs = ',clevs
cbar = fig1.colorbar(CS2,extend='neither', \
   orientation='horizontal',cax=cax,drawedges=True,ticks=clevs,format='%g')
cax.set_xlabel('Accumulated precipitation (mm)',fontsize=12)

# ---- set plot title and save output to a pdf file

plot_title = 'refcst_rawmean_flead'+cleadb+'_to_'+cleade+'_IC'+cyyyymmddhh_init+'.pdf'
fig1.savefig(plot_title)
print 'saving plot to file = ',plot_title

# ---- plot 2: post-processed probability

fig1 = plt.figure(figsize=(9.,6.5))
clevs = [0.01, 0.03, 0.05,0.1, 0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99,1.0]

ax = fig1.add_axes([0.02,.1,0.96,0.85])
ctitle = cleadb+'-'+cleade+'-h S-G smoothed analog P(precip > '+cthresh+\
    '),\ninitialized '+cyyyymmddhh_init+', valid '+cyyyymmddhh_fcst
ax.set_title(ctitle,fontsize=14)
m = Basemap(projection='mill',llcrnrlon=lonsa[0,0],llcrnrlat=latsa[0,0],\
        urcrnrlon=lonsa[-1,-1],urcrnrlat=latsa[-1,-1],resolution='l')
x,y = m(lonsa,latsa)
probfcst_m = ma.array(probfcst, mask=1-iccpa_mask.astype(int))
CS2 = m.contourf(x,y,probfcst_m,clevs,colors=colors_ahps,cmap=None,extend='neither')
CS1 = m.contour(x, y,probfcst_m,clevs,colors=colorstblack,cmap=None,linewidths=0.25)
m.drawcoastlines(linewidth=.5)
m.drawstates(linewidth=.5)
m.drawcountries(linewidth=.5)

cax = fig1.add_axes([0.05,0.09,0.9,0.04])
cbar = fig1.colorbar(CS2,extend='neither', \
   orientation='horizontal',cax=cax,drawedges=True,ticks=clevs,format='%g')
cax.set_xlabel('Probability',fontsize=12)

# ---- set plot title and save output to a pdf file

plot_title = 'refcst_'+cthresh+'_supp'+csupp+'_smoothed_flead'+cleadb+\
    '_to_'+cleade+'_IC='+cyyyymmddhh_init+'.pdf'
fig1.savefig(plot_title)
print 'saving plot to file = ',plot_title

# ---- plot 3: raw ensemble probability

fig1 = plt.figure(figsize=(9.,6.5))
clevs = [0.01, 0.03, 0.05,0.1, 0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99,1.0]

ax = fig1.add_axes([0.02,.1,0.96,0.85])
ctitle = '(c) '+cleadb+'-'+cleade+'-h raw P(precip > '+cthresh+\
    '),\ninitialized '+cyyyymmddhh_init+', valid '+cyyyymmddhh_fcst
ax.set_title(ctitle,fontsize=14)
m = Basemap(projection='mill',llcrnrlon=lonsa[0,0],llcrnrlat=latsa[0,0],\
        urcrnrlon=lonsa[-1,-1],urcrnrlat=latsa[-1,-1],resolution='l')
x,y = m(lonsf,latsf)
print np.shape(x), np.shape(prob_raw)
CS2 = m.contourf(x,y,prob_raw,clevs,colors=colors_ahps,cmap=None,extend='neither')
CS1 = m.contour(x,y,prob_raw,clevs,colors=colorstblack,cmap=None,linewidths=0.25)
m.drawcoastlines(linewidth=.5)
m.drawstates(linewidth=.5)
m.drawcountries(linewidth=.5)

cax = fig1.add_axes([0.05,0.09,0.9,0.04])
cbar = fig1.colorbar(CS2,extend='neither', \
   orientation='horizontal',cax=cax,drawedges=True,ticks=clevs,format='%g')
cax.set_xlabel('Probability',fontsize=12)

# ---- set plot title and save output to a pdf file.

plot_title = 'refcst_raw_p'+cthresh+'_flead'+cleadb+'_to_'+cleade+'_IC='+cyyyymmddhh_init+'.pdf'
fig1.savefig(plot_title)
print 'saving plot to file = ',plot_title




