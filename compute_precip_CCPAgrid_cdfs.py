""" this python routine will take precipitation analysis and forecast data (analysis data on
    the CCPA 1/8-degree grid; forecast data on its native grid) and generate CDF information.
    This data is used by the program that calculates the supplemental data locations,
    compute_precip_analog_locations_ccpa5.f90
"""

import numpy as np
from netCDF4 import Dataset
import sys
import pygrib
from dateutils import hrstodate, daterange, dayofyear, \
     splitdate, datetohrs, dateshift, dateto_hrs_since_day1CE
from mpl_toolkits.basemap import interp
import os
#import scipy.stats as stats
import time as timey
from get_ceedeef_precip import get_ceedeef_precip
from get_quantiles_linear import get_quantiles_linear


input_data_directory  = '/Users/thamill/refcst2/test/netcdf'
#input_data_directory = '/Projects/Reforecast2/netcdf'

output_data_directory =  '/Projects/Reforecast2/netcdf'
cleadb = sys.argv[1]  # enter begin lead time in hours
cleade = sys.argv[2]  # enter end lead time in hours
ileadb = int(cleadb)
ileade = int(cleade)
nmembers = 11
npct = 99 + 8 # quantiles, from 0.01 to 0.99 by 0.01, also (.0001, .005, .001, .0005), 
              # and 1 - (.0001, .005, .001, .0005)

# --- initialize quantiles of the CDF for which we will calculate the associated
#     precipitation amount where this fraction has values below or equal to.

pctiles = [.0001, .0005, .001, .005, \
            .01, .02, .03, .04, .05, .06, .07, .08, .09, .10, \
            .11, .12, .13, .14, .15, .16, .17, .18, .19, .20, \
            .21, .22, .23, .24, .25, .26, .27, .28, .29, .30, \
            .31, .32, .33, .34, .35, .36, .37, .38, .39, .40, \
            .41, .42, .43, .44, .45, .46, .47, .48, .49, .50, \
            .51, .52, .53, .54, .55, .56, .57, .58, .59, .60, \
            .61, .62, .63, .64, .65, .66, .67, .68, .69, .70, \
            .71, .72, .73, .74, .75, .76, .77, .78, .79, .80, \
            .81, .82, .83, .84, .85, .86, .87, .88, .89, .90, \
            .91, .92, .93, .94, .95, .96, .97, .98, .99, \
            .995, .999, .9995, .9999]
            
# ---- define the precipitation amount thresholds that we will calculate CDF at.

xc = range(90)
thresh = [.001,.003,.005,.01,.03, .05,.07,.1,.2,.3,  .4,.5,.6,.7,.8,  \
    .9,1.0,1.2, 1.4, 1.6,    1.8, 2.0, 2.25, 2.5, 2.75,   3.0, 3.5, 4.0, 4.5, 5.0, \
    6.0,7.0,8.0,9.0,10.0,  11.0,12.0,13.0,14.0,15.0,  16.0,17.0,18.0,19.0,19.5,  \
    20.0,22.5,25.,27.5,30.0,   32.5,35.0,37.5,40.0,42.5,  45.0,50.0,55.0,60.0,65.0,  \
    70.0,75.0,80.0,85.0,90.0,   95.0,100.0,105.0,110.0,120.0,  130.0,140.0,150.0,160.0,170.0,  \
    180.0,190.0,200.0,220.0,240.0,   260.0,280.0,300.0,325.0,350.0,   400.0,500.0,600.0,700.0,1000.]
nthresh = len(xc)
print 'len xc, thresh = ',len(xc), len(thresh)

# ---- determine analysis grid lat/lon information from a sample file

filename = input_data_directory + '/refcstv2_precip_ccpav3_024_to_048.nc'
print filename
nc = Dataset(filename)
lats_anal = nc.variables['lats_anal'][:]
lons_anal = nc.variables['lons_anal'][:]
lats_fcst = nc.variables['lats_fcst'][:]
lons_fcst = nc.variables['lons_fcst'][:]
conusmask_ccpa = nc.variables['conusmask'][:]
xavals = nc.variables['xa'][:]
yavals = nc.variables['ya'][:]
xfvals = nc.variables['xf'][:]
yfvals = nc.variables['yf'][:]
nia = len(xavals)
nja = len(yavals)
nif = len(xfvals)
njf = len(yfvals)
nc.close()

# ---- define the corner lat/lons.

llcrnrlat = lats_anal[0,0]
llcrnrlon = lons_anal[0,0]
urcrnrlat = lats_anal[-1,-1]
urcrnrlon = lons_anal[-1,-1]

# ---- determine nearest forecast grid point for each analysis lat/lon

ixmin = np.zeros((nia),dtype=np.int)
jymin = np.zeros((nja),dtype=np.int)
for ixa in range(nia):
    dmin = 999999999.
    for ixf in range(nif):
        d = np.abs(lons_anal[0,ixa]-lons_fcst[0,ixf])
        if d < dmin:
            ixmin[ixa] = ixf
            dmin = d
for jya in range(nja):
    dmin = 999999999.
    for jyf in range(njf):
        d = np.abs(lats_anal[jya,0]-lats_fcst[jyf,0])
        if d < dmin:
            jymin[jya] = jyf
            dmin = d

#print 'samples nearest lon index = ',ixmin[:]
#print 'samples nearest lat index = ',jymin[:]


# --- loop thru and process each month, reading that data plus the surrounding two months.

cmonths = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
icenter = [14,45,73,104,134,165,195,226,257,287,318,348] # 15th of each month, in Julian days
for imonth, cmonth in zip(range(12), cmonths): 
    print 'processing month = ',cmonth

    # ---- initialize CDFs and counter

    CDFf = np.zeros((nthresh,njf,nif), dtype=np.float64) # forecast
    CDFa = np.zeros((nthresh,nja,nia), dtype=np.float64) # analyzed
    CDFworkf = np.zeros((nthresh,njf,nif), dtype=np.float64)
    CDFworka = np.zeros((nthresh,nja,nia), dtype=np.float64)
    icdf = 0

    # ---- open and initialize the netCDF file we'll be writing to.

    outfilename = output_data_directory+\
        '/refcstv2_apcp_CCPAgrid_CDF_fhour'+cleadb+'_to_'+cleade+'_'+cmonth+'.nc'

    print outfilename
    rootgrp = Dataset(outfilename,'w',format='NETCDF4_CLASSIC')
    xf = rootgrp.createDimension('xf',nif)
    xfv = rootgrp.createVariable('xf','f4',('xf',))
    xfv.long_name = "forecast grid eastward distance from southwest corner of domain in grid points" 
    xfv.units = "grid index (dimensionless)" 
    yf = rootgrp.createDimension('yf',njf)
    yfv = rootgrp.createVariable('yf','f4',('yf',))
    yfv.long_name = "forecast grid northward distance from southwest corner of domain in grid points" 
    yfv.units = "grid index (dimensionless)"
    xa = rootgrp.createDimension('xa',nia)
    xav = rootgrp.createVariable('xa','f4',('xa',))
    xav.long_name = "analysis grid eastward distance from southwest corner of domain in grid points" 
    xav.units = "grid index (dimensionless)" 
    ya = rootgrp.createDimension('ya',nja)
    yav = rootgrp.createVariable('ya','f4',('ya',))
    yav.long_name = "analysis grid northward distance from southwest corner of domain in grid points" 
    yav.units = "grid index (dimensionless)"
    pct = rootgrp.createDimension('pct',npct)
    pctv = rootgrp.createVariable('pct','f4',('pct'))
    pctv.long_name = "quantiles of the distribution"
    pctv.units = "fraction"
    thrnum = rootgrp.createDimension('thrnum',nthresh)
    thrnumv = rootgrp.createVariable('thrnum','i4',('thrnum',))
    thrnumv.long_name = "Threshold iterator (0:nthresh)" 
    thrnumv.units = " "  
    thrval = rootgrp.createDimension('thrval',nthresh)
    thrvalv = rootgrp.createVariable('thrval','f4',('thrval',))
    thrvalv.long_name = "Precip thresholds (mm) that precip_CDFs evaluated at" 
    thrvalv.units = "K"  
    time = rootgrp.createDimension('time',None)
    timev = rootgrp.createVariable('time','f4',('time',))
    timev.units = "hours since 1-1-1 00:00:0.0" 
    lonsa = rootgrp.createVariable('lonsa','f4',('ya','xa',))
    lonsa.long_name = "longitude" 
    lonsa.units = "degrees_east" 
    latsa = rootgrp.createVariable('latsa','f4',('ya','xa',))
    latsa.long_name = "latitude" 
    latsa.units = "degrees_north" 
    lonsf = rootgrp.createVariable('lonsf','f4',('yf','xf',))
    lonsf.long_name = "longitude" 
    lonsf.units = "degrees_east" 
    latsf = rootgrp.createVariable('latsf','f4',('yf','xf',))
    latsf.long_name = "latitude" 
    latsf.units = "degrees_north" 

    conusmask = rootgrp.createVariable('conusmask','i2',('ya','xa',))
    conusmask.long_name = "mask for grid points inside CONUS (1=yes,0=no)"
    conusmask.units=""

    pforecast_CDF = rootgrp.createVariable('pforecast_CDF','f4',('thrnum','yf','xf',),
        zlib=True,least_significant_digit=3)  
    pforecast_CDF.units = "" 
    pforecast_CDF.long_name = "Cumulative distribution function of forecast precip, defined at thrval" 
    pforecast_CDF.valid_range = [0.0,1.0]
    pforecast_CDF.missing_value = -99.99

    panal_CDF = rootgrp.createVariable('panal_CDF','f4',('thrnum','ya','xa',),
        zlib=True,least_significant_digit=3)  
    panal_CDF.units = "" 
    panal_CDF.long_name = "Cumulative distribution function of forecast precip, defined at thrval" 
    panal_CDF.valid_range = [0.0,1.0]
    panal_CDF.missing_value = -99.99

    pforecast_quantiles = rootgrp.createVariable('pforecast_quantiles','f4',('pct','yf','xf',),
        zlib=True,least_significant_digit=2)  
    pforecast_quantiles.units = "mm" 
    pforecast_quantiles.long_name = "Precip amount associated with quantiles of forecast distribution" 
    pforecast_quantiles.valid_range = [0.0,10000.]
    pforecast_quantiles.missing_value = -99.99

    panal_quantiles = rootgrp.createVariable('panal_quantiles','f4',('pct','ya','xa',),
        zlib=True,least_significant_digit=2)  
    panal_quantiles.units = "mm" 
    panal_quantiles.long_name = "Precip amount associated with quantiles of analyzed distribution" 
    panal_quantiles.valid_range = [0.0,10000.]
    panal_quantiles.missing_value = -99.99

    pctv[:] = pctiles[:]

    xav[:]   = np.arange(nia)
    yav[:]   = np.arange(nja)
    xfv[:]   = np.arange(nif)
    yfv[:]   = np.arange(njf)
    lonsa[:]  = lons_anal
    latsa[:]  = lats_anal
    lonsf[:]  = lons_fcst
    latsf[:]  = lats_fcst

    conusmask[:] = conusmask_ccpa

    thrnumv[:] = xc[:]
    thrvalv[:] = thresh[:]
 
    rootgrp.latcorners = [lats_anal[0,0], lats_anal[0,-1], lats_anal[-1,0], lats_anal[-1,-1]]
    rootgrp.loncorners = [lons_anal[0,0], lons_anal[0,-1], lons_anal[-1,0], lons_anal[-1,-1]]

    rootgrp.stream = "s4" # ????
    rootgrp.title = "Reforecast V2 accum. ensemble-mean precip forecast and analyzed CDF + rank correlation"
    rootgrp.Conventions = "CF-1.0"  # ????
    rootgrp.history = "Revised Mar 2016 by Hamill" 
    rootgrp.institution = \
        "Reforecast from ERSL/PSD using NCEP/EMC GEFS, circa 2012"
    rootgrp.platform = "Model" 
    rootgrp.references = "http://www.esrl.noaa.gov/psd/forecasts/reforecast2/" 
    
    # ---- open ensemble data file for each year, read in data, and augment cdf 
    #      information for that year if the sample is within the month of interest
    #      or the neighboring month

    rankcorr_fa = -99.99*np.ones((nja,nia),dtype=np.float)    
    nyears = len(range(2002,2016))
    print 'nsamps = ',92*nyears
    precipa = np.zeros((nja,nia),dtype=np.float)
    precipf3d = np.zeros((92*nyears,njf,nif),dtype=np.float)
    precipa3d = np.zeros((92*nyears,nja,nia),dtype=np.float)
    ipktr = 0

    infilename = input_data_directory+'/refcstv2_precip_ccpav3_'+\
       cleadb+'_to_'+cleade+'.nc'
    cfile1 = Dataset(infilename,"r")
    yyyymmddhh_init = cfile1.variables['yyyymmddhh_init'][:]
    ensv = cfile1.variables['ensv'][:]
    nmembers = len(ensv)
    precipf = np.zeros((njf,nif),dtype=np.float)
    precipfens = np.zeros((nmembers,njf,nif),dtype=np.float)

    print 'looping thru days'
    for iyyyymmddhh, iday in zip(yyyymmddhh_init,range(len(yyyymmddhh_init))):
        if iday % 90 == 0: print iyyyymmddhh, imonth
        cyyyymmddhh = str(iyyyymmddhh)
        cmonth2 = cyyyymmddhh[4:6]
        imonth2 = int(cmonth2) - 1
        if imonth2 == 11 and imonth == 0: imonth2 = -1
        if imonth2 == 0 and imonth == 11: imonth2 = 12
        if imonth2 == imonth-1 or imonth2 == imonth or imonth2 == imonth+1:
            precipfens[:,:,:] = cfile1.variables['apcp_fcst_ens'][iday,:,:,:]
            precipf = np.mean(precipfens,axis=0)
            precipf3d[ipktr,:,:] = precipf[:,:]
            CDFworkf, istatf = get_ceedeef_precip(nthresh,njf,nif,thresh,precipf)
            precipa[:,:] = cfile1.variables['apcp_anal'][iday,:,:]
            precipa3d[ipktr,:,:] = precipa[:,:]
            ipktr = ipktr + 1
            CDFworka, istata = get_ceedeef_precip(nthresh,nja,nia,thresh,precipa)
            if istatf == 1 and istata == 1:
                CDFf = CDFf + CDFworkf
                CDFa = CDFa + CDFworka
                icdf = icdf + 1
    cfile1.close()

    # --- write out the CDF record

    print 'writing out CDF record',timey.localtime()
    flicdf = np.float(icdf)
    for it in range(nthresh):
        CDFf[it,:,:] = CDFf[it,:,:] / flicdf
        CDFa[it,:,:] = CDFa[it,:,:] / flicdf
    pforecast_CDF[:] = CDFf
    panal_CDF[:] = CDFa
    print 'sample pforecast_CDF = ',CDFf[:,njf/2,nif/2]
    print 'sample panal_CDF     = ',CDFa[:,nja/2,nia/2]

    # --- get quantiles

    print 'getting precip amount associated with quantiles of distribution for forecast'
    precip_qret = np.zeros((npct,njf,nif),dtype=np.float)
    print 'nthresh, npct, njf, nif = ', nthresh, npct, njf, nif 
    print 'pctiles = ', pctiles
    print 'thresh = ', thresh
    print 'np.shape(CDFf) = ', np.shape(CDFf)
    precip_qret,istat = get_quantiles_linear(nthresh,npct,njf,nif,pctiles,thresh,CDFf)
    pforecast_quantiles[:] = precip_qret
    print 'sample pforecast_quantiles = ',precip_qret[:,njf/2,nif/2]

    print 'getting precip amount associated with quantiles of distribution for analyzed'
    precip_qret = np.zeros((npct,nja,nia),dtype=np.float)
    precip_qret,istat = get_quantiles_linear(nthresh,npct,nja,nia,pctiles,thresh,CDFa)
    panal_quantiles[:] = precip_qret
    print 'sample panal_quantiles     = ',precip_qret[:,nja/2,nia/2]

    rootgrp.close()
