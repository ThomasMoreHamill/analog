""" this script is designed to read in grib data of 1/8-degree CCPA precipitation analyses
    as well as precipitation and precipitable water reforecast data for a particular forecast
    lead time, and save them to a netcdf file.  Coded by Tom Hamill, Mar 2016, 
    tom.hamill@noaa.gov """
    
import numpy as np
import sys
import pygrib
import os
import time as timey
from netCDF4 import Dataset
from dateutils import hrstodate, daterange, dayofyear, \
     splitdate, datetohrs, dateshift, dateto_hrs_since_day1CE
from mpl_toolkits.basemap import interp
import cPickle

def read_gribdata(idxfilename, gribfilename, validityDate, validityTime):
    """ reads in data from a particular grib file with a particular date and time"""
    grb=[]
    istat = -1
    for iter in range(10):
        fexist_grib = False
        fexist_grib = os.path.exists(gribfilename)
        fexist_idx  = False
        fexist_idx  = os.path.exists(idxfilename)
        if fexist_grib and fexist_idx:
            try:
                fcstfile = pygrib.index(idxfilename)
                grb = fcstfile.select(validityDate=validityDate, \
                    validityTime=validityTime)[0]
                istat = 0
                fcstfile.close()
                return istat,grb
            except (IOError, ValueError, RuntimeError):
                print 'Error reading ', gribfilename
                istat = -1
    return istat,grb

def form_filename(week1, cyear, cyearmo, date_IC, g_or_latlon, mbr_mean_or_sprd, variable):
    """compose the correct grib file name and index file name to read for the 
    forecast data of interest"""
    istat = -1
    if week1 == True:
        filename_d = '/Projects/Reforecast2/'+cyear+'/'+cyearmo+'/'+ date_IC+ \
            '/'+mbr_mean_or_sprd+'/'+g_or_latlon+'/'+variable+date_IC+'_'+mbr_mean_or_sprd+'.grib2'
        filename = '/Projects/Reforecast2/'+cyear+'/'+cyearmo+'/'+ date_IC+ \
            '/'+mbr_mean_or_sprd+'/'+g_or_latlon+'/'+variable+date_IC+'_'+mbr_mean_or_sprd+'.grib2.pyidx'
    else:
        filename_d = '/Projects/Reforecast2/'+cyear+'/'+cyearmo+'/'+ \
            date_IC+ '/'+mbr_mean_or_sprd+'/'+g_or_latlon+'/'+variable+date_IC+'_'+mbr_mean_or_sprd+'_t190.grib2'
        filename = '/Projects/Reforecast2/'+cyear+'/'+cyearmo+'/'+ \
            date_IC+ '/'+mbr_mean_or_sprd+'/'+g_or_latlon+'/'+variable+date_IC+'_'+mbr_mean_or_sprd+'_t190.grib2.pyidx'
    return istat, filename, filename_d

# =====================================================================================

cleadb = sys.argv[1]  # beginning lead time in hours
cleade = sys.argv[2]  # end lead time in hours

ileadb = int(cleadb)  # converted to integers
ileade = int(cleade)
nmembers = 11 # reforecast data has 11 members

# ---- initialize

date_begin = '2002010200'  # beginning date
date_end   = '2015123100'  # ending date
date_list  = daterange(date_begin, date_end, 24)
ndates = len(date_list)
cmembers = ['c00','p01','p02','p03','p04','p05','p06','p07','p08','p09','p10']

# ---- read in conus mask

infile = '/Users/thamill/refcst2/git_archive/conusmask_ccpa.nc'
nc = Dataset(infile)
conusmask_in = nc.variables['conusmask'][:,:]
lons_anal = nc.variables['lons'][:,:]
lats_anal = nc.variables['lats'][:,:]
latsa1d = lats_anal[:,0]
lonsa1d = lons_anal[0,:]
nlatsa, nlonsa = lons_anal.shape
#print 'nlonsa, nlatsa = ', nlonsa, nlatsa
coslatsa = np.cos((np.pi/180.)*lats_anal)
zeros = np.zeros((nlatsa,nlonsa),dtype=np.float)
ones = np.ones((nlatsa,nlonsa),dtype=np.float)
nc.close()

# ---- read in sample forecast lat/lon information for week +1 Gaussian grid

infile = '/Users/thamill/refcst2/apcp_sfc_2010010100_c00.grib2'
flatlon = pygrib.open(infile)
idateno00 = 20100101
fcst = flatlon.select(shortName='tp', level=0, dataDate=idateno00)[0]
latsf_week1, lonsf_week1 = fcst.latlons()
nlatsf_week1, nlonsf_week1 = lonsf_week1.shape
flatlon.close()
nib_week1 = 245 # latitude to chop out
nie_week1 = 306
njb_week1 = 500 # longitude
nje_week1 = 628
ni_week1 = nie_week1 - nib_week1 
nj_week1 = nje_week1 - njb_week1 
lats1f_week1 = latsf_week1[:,0]
if lats1f_week1[0] > 0:
    flipit = True
    lats1f_week1 = lats1f_week1[::-1]  # note flipping of data to make ascending
    latsf_week1 = np.flipud(latsf_week1)
else:
    flipit = False
lonsf_week1 = lonsf_week1 - 360.
lons1f_week1 = lonsf_week1[0,:] 

# ---- read in sample forecast lat/lon information for week +2 Gaussian grid

infile = '/Users/thamill/refcst2/apcp_sfc_2010010100_c00_t190.grib2'
flatlon = pygrib.open(infile)
idateno00 = 20100101
fcst = flatlon.select(shortName='tp',level=0,dataDate=idateno00)[0]
latsf_week2, lonsf_week2 = fcst.latlons()
nlatsf_week2, nlonsf_week2 = lonsf_week2.shape
flatlon.close()
lats1f_week2 = latsf_week2[:,0]
if lats1f_week2[0] > 0:
    lats1f_week2 = lats1f_week2[::-1]  # note flipping of data to make ascending
    latsf_week2 = np.flipud(latsf_week2)
lonsf_week2 = lonsf_week2 - 360.
lons1f_week2 = lonsf_week2[0,:] 

# ---- set the Gaussian lat / lon sub-array for week +1

lons_fsub_week1 = np.zeros((ni_week1,nj_week1), dtype=np.float)
lats_fsub_week1 = np.zeros((ni_week1,nj_week1), dtype=np.float)
lons_fsub_week1[:,:] = lonsf_week1[nib_week1:nie_week1,njb_week1:nje_week1]
lats_fsub_week1[:,:] = latsf_week1[nib_week1:nie_week1,njb_week1:nje_week1]

# ---- read in sample forecast lat/lon information for latlon grid

infile = '/Users/thamill/refcst2/apcp_sfc_2010010100_c00_latlon.grib2'
flatlon = pygrib.open(infile)
idateno00 = 20100101
fcst = flatlon.select(shortName='tp', level=0, dataDate=20100101)[0]
latsf_1deg, lonsf_1deg = fcst.latlons()
nlatsf_1deg, nlonsf_1deg = lonsf_1deg.shape
latsf_1deg_1d = np.zeros((nlatsf_1deg),dtype=np.float)
lonsf_1deg_1d = np.zeros((nlonsf_1deg),dtype=np.float)
latsf_1deg_1d[:] = latsf_1deg[:,0]
lonsf_1deg_1d[:] = lonsf_1deg[0,:]
latsf00 = latsf_1deg_1d[0]

# --- there was a grib routine change so that data can no longer be accessed
#     properly if lat order is descending as opposed to ascending.  Check and
#     correct if so.

if latsf00 > 0:
    latsf_1deg_1d = latsf_1deg_1d[::-1]  # note flipping of data

flatlon.close()
njb_1deg = 235 # lon
nje_1deg = 293
nib_1deg = 115 # lat
nie_1deg = 143


# ---- set up netCDF file particulars

outfilename = '/Users/thamill/refcst2/test/netcdf/refcstv2_precip_ccpav3_'+\
    cleadb+'_to_'+cleade+'.nc'

print 'writing netCDF data to ',outfilename
rootgrp = Dataset(outfilename,'w',format='NETCDF4_CLASSIC')

xa = rootgrp.createDimension('xa',nlonsa)
xva = rootgrp.createVariable('xa','f4',('xa',))
xva.long_name = "eastward grid point number, precip analysis" 
xva.units = "n/a" 

ya = rootgrp.createDimension('ya',nlatsa)
yva = rootgrp.createVariable('ya','f4',('ya',))
yva.long_name = "northward grid point number, precip analysis" 
yva.units = "n/a" 

xf = rootgrp.createDimension('xf',nj_week1)
xvf = rootgrp.createVariable('xf','f4',('xf',))
xvf.long_name = "eastward grid point number, precip forecast" 
xvf.units = "n/a" 

yf = rootgrp.createDimension('yf',ni_week1)
yvf = rootgrp.createVariable('yf','f4',('yf',))
yvf.long_name = "northward grid point number, precip forecast" 
yvf.units = "n/a" 

time = rootgrp.createDimension('time',None)
timev = rootgrp.createVariable('time','f4',('time',))
timev.units = "hours since 1-1-1 00:00:0.0" 

ens = rootgrp.createDimension('ens',11)
ensv = rootgrp.createVariable('ensv','i4',('ens',))
ensv.long_name = "Ensemble member number (control, perts 1-10)" 
ensv.units = " " 

lonsa = rootgrp.createVariable('lons_anal','f4',('ya','xa',))
lonsa.long_name = "longitude" 
lonsa.units = "degrees_east" 

latsa = rootgrp.createVariable('lats_anal','f4',('ya','xa',))
latsa.long_name = "latitude" 
latsa.units = "degrees_north" 

lonsf = rootgrp.createVariable('lons_fcst','f4',('yf','xf',))
lonsf.long_name = "longitude" 
lonsf.units = "degrees_east" 

latsf = rootgrp.createVariable('lats_fcst','f4',('yf','xf',))
latsf.long_name = "latitude" 
latsf.units = "degrees_north" 

conusmaska = rootgrp.createVariable('conusmask','i2',('ya','xa',))
conusmaska.long_name = "mask (1=CONUS, 0=outside of CONUS)"
conusmaska.units = "none"

yyyymmddhh_init2 = rootgrp.createVariable('yyyymmddhh_init','i4',('time',))
yyyymmddhh_init2.longname = "Initial condition date/time in yyyymmddhh format"

yyyymmddhh_fcst2b = rootgrp.createVariable('yyyymmddhh_fcstb','i4',('time',))
yyyymmddhh_fcst2b.longname = "Forecast valid date/time in yyyymmddhh format"

yyyymmddhh_fcst2e = rootgrp.createVariable('yyyymmddhh_fcste','i4',('time',))
yyyymmddhh_fcst2e.longname = "Forecast valid date/time in yyyymmddhh format"

apcp_anal = rootgrp.createVariable('apcp_anal','f4',('time','ya','xa',), \
   zlib=True,least_significant_digit=2)
apcp_anal.units = "mm" 
apcp_anal.long_name = "Analyzed Accumulated Precipitation from CCPA" 
apcp_anal.valid_range = [0.,1000.]
apcp_anal.missing_value = -99.99

apcp_fcst_ens = rootgrp.createVariable('apcp_fcst_ens','f4',('time','ens','yf','xf',),
   zlib=True,least_significant_digit=2)  
apcp_fcst_ens.units = "mm" 
apcp_fcst_ens.long_name = "Ensemble member forecast accumulated precipitation" 
apcp_fcst_ens.valid_range = [0.,1000.]
apcp_fcst_ens.missing_value = -99.99

pwat_fcst_mean = rootgrp.createVariable('pwat_fcst_mean','f4',('time','yf','xf',),
   zlib=True,least_significant_digit=2)  
pwat_fcst_mean.units = "mm" 
pwat_fcst_mean.long_name = "Ensemble precipitable water forecast mean, time averaged"
pwat_fcst_mean.valid_range = [0.,1000.]
pwat_fcst_mean.missing_value = -99.99

conusmaska[:] = conusmask_in

xva[:] = np.arange(nlonsa)
yva[:] = np.arange(nlatsa)

xvf[:] = np.arange(nj_week1)
yvf[:] = np.arange(ni_week1)

lonsa[:] = lons_anal
latsa[:] = lats_anal

lonsf[:] = lons_fsub_week1
latsf[:] = lats_fsub_week1 

ensv[:] = range(11)

rootgrp.stream = "s4" # ????
rootgrp.title = \
   "Reforecast V2 accumulated forecast and CCPA 1/8 degree precip analysis over CONUS"
rootgrp.Conventions = "CF-1.0"  # ????
rootgrp.history = "Created March 2016 by Tom Hamill" 
rootgrp.institution = \
   "CCPA from NCEP/EMC, Reforecast from ERSL/PSD using NCEP/EMC GEFS circa 2012"
rootgrp.platform = "Model" 
rootgrp.references = "http://www.esrl.noaa.gov/psd/forecasts/reforecast2/" 

ktr = 0
apcpf = np.zeros((nmembers,ni_week1,nj_week1),dtype=np.float)
pwatf = np.zeros((ni_week1,nj_week1),dtype=np.float)
missing_values_anal = -99.99 * np.ones((nlatsa,nlonsa),dtype=np.float)

# ---- loop thru all dates ...
ktr = 0

for idate, date_IC in zip(range(ndates),date_list):

    print 'processing record number ',idate,'  initial date = ',date_IC

    # ---- read in CCPA 6-hourly precipitation analysis data

    precipver = False
    date_ppend   = dateshift(date_IC, ileade+6)
    date_FCend   = dateshift(date_IC, ileade)
    date_FCbegin = dateshift(date_IC, ileadb)
    date_ppbegin = dateshift(date_IC, ileadb+6)
    date_lista   = daterange(date_ppbegin,date_ppend,6)
    tp_anal = np.zeros((nlatsa, nlonsa),dtype=np.float)
    nogood = False
    for datea in date_lista:
        cyear = datea[0:4]
        chour = datea[8:10]
        
        cyyyymmddhh = datea[0:10]
        cyyyymmdd = datea[0:8]
        #infile = '/Rf2_tests/ccpa_v1/0.125d/ccpa.'+cyyyymmdd+'/'+chour+ \
        #    '/ccpa.t'+chour+'z.06h.0p125.conus'
        infile = '/Users/thamill/ccpa_v1/ccpa.'+cyyyymmdd+'/'+chour+ \
            '/ccpa.t'+chour+'z.06h.0p125.conus'
        fexist_apcp = os.path.exists(infile)
        if fexist_apcp and nogood == False:
            ffcst = pygrib.open(infile)
            grb = ffcst.select(shortName='tp')[0]
            #tp_anal = tp_anal + np.flipud(grb.values)
            tp_anal = tp_anal + grb.values
            ffcst.close()
        else: 
            nogood = True
            tp_anal[:,:] = -99.99
            print 'missing precip analysis data for ', idate, ileade, date_IC, date_FCend
    
    tp_anal = np.where(tp_anal < 80000, tp_anal, missing_values_anal)
    
    # --- now read in the associated precipitation ensemble forecasts on the global grid.

    apcpf[:,:,:] = 0.0
    cyear = date_IC[0:4]
    chour = date_IC[8:10]
    cyyyymmddhh = date_IC
    cyearmo = date_IC[0:6]
    for imem, cmem in zip(range(11), cmembers): 

        # ---- load the forecast data

        date_ppend   = dateshift(date_IC, ileade+6)
        date_ppbegin = dateshift(date_IC, ileadb+6)
        date_lista   = daterange(date_ppbegin,date_ppend,6)
        leadrange = range(ileadb+6,ileade+6,6)

        for ilead, datea in zip(leadrange, date_lista):
            validityDate = int(datea)/100
            validityTime = (int(datea) - validityDate*100)*100
            g_or_latlon = 'gaussian'
            mbr_mean_or_sprd = cmem
            if ilead <= 192:
                week1 = True
            else:
                week1 = False
            istat_filename, filename_apcp, filename_apcp_d = \
                form_filename(week1, cyear, cyearmo, date_IC, g_or_latlon, mbr_mean_or_sprd,'apcp_sfc_')
            istat_apcp, grb_apcp = read_gribdata(filename_apcp, filename_apcp_d, validityDate, validityTime)
            if istat_apcp == 0:
                grbdata = grb_apcp.values
                if flipit == True: grbdata = np.flipud(grbdata)
                if week1 == False:
                    # --- interpolate to week +1 grid
                    grbdata_week1grid = interp(grbdata, lons1f_week2, lats1f_week2, \
                        lons_fsub_week1, lats_fsub_week1, checkbounds=False, masked=False, order=1)
                    apcpf[imem,:,:] = apcpf[imem,:,:] + grbdata_week1grid[:,:]
                else:
                    apcpf[imem,:,:] = apcpf[imem,:,:] + grbdata[nib_week1:nie_week1,njb_week1:nje_week1]
            else:
                apcpf[imem,:,:] = -99.99
                fexist_apcp = False
                
    # --- read in precipitable water forecast mean on latlon grid
    
    date_ppend   = dateshift(date_IC, ileade+6)
    date_ppbegin = dateshift(date_IC, ileadb)
    date_lista   = daterange(date_ppbegin,date_ppend,6)
    pwatf[:,:] = 0.0
    nda = len(date_lista)
    for datea in date_lista:
        g_or_latlon = 'latlon'
        mbr_mean_or_sprd = 'mean'
        istat_filename, filename_pwat, filename_pwat_d = \
            form_filename(week1, cyear, cyearmo, date_IC, g_or_latlon, mbr_mean_or_sprd,'pwat_eatm_')
        istat_pwat, grb_pwat = read_gribdata(filename_pwat, filename_pwat_d, validityDate, validityTime)
        if istat_pwat == 0:
            grbdata = grb_pwat.values
            if latsf00 > 0:
                grbdata = np.flipud(grbdata)
            pwat_meanfcst_gauss = interp(grbdata, lonsf_1deg_1d, latsf_1deg_1d, \
                lons_fsub_week1, lats_fsub_week1, checkbounds=False, masked=False, order=1)
            pwatf[:,:] = pwatf[:,:]+ pwat_meanfcst_gauss[:,:]/np.float(nda)
        else:
            pwatf[:,:] = -99.99
            fexist_pwat = False

    ds = dateshift( date_IC, ileade)
    timev[ktr]             = datetohrs(ds)  # in our file, since 0 AD
    yyyymmddhh_init2[ktr]  = date_IC
    yyyymmddhh_fcst2b[ktr] = date_FCbegin
    yyyymmddhh_fcst2e[ktr] = date_FCend
    apcp_anal[ktr]         = tp_anal[:,:]
    apcp_fcst_ens[ktr]     = apcpf[:,:,:]
    pwat_fcst_mean[ktr]    = pwatf[:,:]
    ktr = ktr + 1 

print 'Data was written to ',outfilename
rootgrp.close()
print 'Normal termination.'

