""" to ensure that the supplemental data locations look realistic, plot them
    for a small number of locations around the CONUS, underplotting q95 of
    the precipitation climatology
"""

from netCDF4 import Dataset
import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap, addcyclic
import scipy.ndimage
import numpy as np
from numpy import ma
import math
import os, sys
from read_ppn_analog_locns_ccpa5 import read_ppn_analog_locns_ccpa5
from matplotlib import rcParams
rcParams['xtick.labelsize']='x-small'
rcParams['ytick.labelsize']='x-small'
rcParams['contour.negative_linestyle']='solid'
rcParams['path.simplify'] = False

def add_analog_locns(ycoord,xcoord,ylocs,xlocs,lons,lats,map,symbolshape,color):
    istat=0
    nsize = len(ylocs)
    #print '**** xcoord, ycoord = ',xcoord-1,ycoord-1
    # --- first plot the original point in a larger font
    x,y = map(lons[ycoord,xcoord],lats[ycoord,xcoord])
    map.scatter(x,y,s=35,color=color,marker=symbolshape,zorder=12)
    for i in range(nsize):
        #print 'xlocs, ylocs, i = ',xlocs[i]-1,ylocs[i]-1,i
        x,y = map(lons[ylocs[i]-1,xlocs[i]-1],lats[ylocs[i]-1,xlocs[i]-1])
        map.scatter(x,y,s=13,color=color,marker=symbolshape,zorder=12)
    return istat

def find_nearest_latlon(lonin, latin, rlon, rlat):
    istat = 0
    dmin = 999999.
    nx, ny = np.shape(rlon)
    #print 'ny, nx = ',ny,nx
    for i in range(nx):
        dist = np.abs(lonin - rlon[i,0])
        if dist < dmin:
            dmin = dist
            ilon = i
    dmin = 999999.
    for j in range(ny):
        dist = np.abs(latin - rlat[0,j])
        if dist < dmin:
            dmin = dist
            ilat = j
    return istat, ilon, ilat

# ---- hardcode the lead times for which supplemental locations calculated, but 
#      enter month on the command line.

cleadb = '024' #sys.argv[1]
cleade = '048' #sys.argv[2]
cmonth = sys.argv[1] # Jan, Feb, etc

# ---- read in analog locations

nxa = 464
nya = 224
nanalogs = 20

xlocationa = np.zeros((nxa,nya,nanalogs),dtype=np.int32)
ylocationa = np.zeros((nxa,nya,nanalogs),dtype=np.int32)
conusmaska = np.zeros((nxa,nya),dtype=np.int16)
rlona = np.zeros((nxa,nya),dtype=np.float)
rlata = np.zeros((nxa,nya),dtype=np.float)
xlocationa, ylocationa, conusmaska, rlona, rlata = \
    read_ppn_analog_locns_ccpa5(cmonth, cleadb, cleade, nxa, nya, nanalogs)
#print 'conusmaska(nxa/2,nya/2) = ',conusmaska[nxa/2, nya/2]
#print 'xlocationa[nxa/2,nya/2,:] = ',xlocationa[nxa/2,nya/2,:]
#print 'ylocationa[nxa/2,nya/2,:] = ',ylocationa[nxa/2,nya/2,:]
#print 'max, min rlona =',np.max(rlona), np.min(rlona)
#print 'max, min rlata =',np.max(rlata), np.min(rlata)

# ---- get nearest gridpoint to several major cities

istat, iber_a, jber_a = find_nearest_latlon(-122.0, 37.95, rlona, rlata)
istat, ipor_a, jpor_a = find_nearest_latlon(-122.6, 45.54, rlona, rlata)
istat, iboz_a, jboz_a = find_nearest_latlon(-111.05, 45.67, rlona, rlata)
istat, ibis_a, jbis_a = find_nearest_latlon(-101.7, 46.8, rlona, rlata)
istat, ibou_a, jbou_a = find_nearest_latlon(-105.23, 40.02, rlona, rlata)
istat, ilub_a, jlub_a = find_nearest_latlon(-101.88, 33.59, rlona, rlata)
istat, iphx_a, jphx_a = find_nearest_latlon(-112.1, 33.6, rlona, rlata)
istat, ilit_a, jlit_a = find_nearest_latlon(-92.23, 34.72, rlona, rlata)
istat, iatl_a, jatl_a = find_nearest_latlon(-84.42, 33.76, rlona, rlata)
istat, ichi_a, jchi_a = find_nearest_latlon(-87.73, 41.83, rlona, rlata)
istat, icin_a, jcin_a = find_nearest_latlon(-84.54, 39.13, rlona, rlata)
istat, inyc_a, jnyc_a = find_nearest_latlon(-73.97, 40.7, rlona, rlata)
istat, ioma_a, joma_a = find_nearest_latlon(-96.04, 41.29, rlona, rlata)

print 'iber_a, jber_a = ',iber_a,jber_a

# ---- read in the climatology of forecast and analyzed quantiles

infile = '/Projects/Reforecast2/netcdf/refcstv2_apcp_CCPAgrid_CDF_fhour'+\
  cleadb+'_to_'+cleade+'_'+cmonth+'.nc'

print infile
nc = Dataset(infile)
panal_quantiles     = nc.variables['panal_quantiles'][:]
rankcorr            = nc.variables['rankcorr'][:]
lonsa               = nc.variables['lonsa'][:]
latsa               = nc.variables['latsa'][:]
nc.close()
 

if cmonth == 'Jan' : imonth = 0
if cmonth == 'Feb' : imonth = 1
if cmonth == 'Mar' : imonth = 2
if cmonth == 'Apr' : imonth = 3
if cmonth == 'May' : imonth = 4
if cmonth == 'Jun' : imonth = 5
if cmonth == 'Jul' : imonth = 6
if cmonth == 'Aug' : imonth = 7
if cmonth == 'Sep' : imonth = 8
if cmonth == 'Oct' : imonth = 9
if cmonth == 'Nov' : imonth = 10
if cmonth == 'Dec' : imonth = 11

#  ---- prep arrays for display

xlocationa_t = np.transpose(xlocationa)
ylocationa_t = np.transpose(ylocationa)
rlatsa_t     = np.transpose(rlata)
rlonsa_t     = np.transpose(rlona)
conusmaska_t = np.transpose(conusmaska)

#print 'np.shape(panal_quantiles) = ',np.shape(panal_quantiles)
#print 'ylocationa_t[:,nya/2,nxa/2,:] = ',ylocationa_t[:,nya/2,nxa/2] 
#print 'xlocationa_t[:,nya/2,nxa/2,:] = ',xlocationa_t[:,nya/2,nxa/2] 

paquant = panal_quantiles[98,:,:]
paquant50 = panal_quantiles[53,:,:]
paquant_m = ma.array(paquant, mask = 1-conusmaska_t)

colorst = ['White','BurlyWood','LightGreen','DarkSeaGreen','Gold',\
           'Orange','Red','MediumOrchid','DeepSkyBlue','LightGray']
colorstblack=['White','Black','Black','Black','Black', 'Black','Black','Black','Black','Black']

colorst = ['White','#ECFFFF','#D9F7FF','#C4E8FF','#E8FBE8','#C7F4C7','#92F592','Yellow','Gold',\
  'Orange','#FFB2B2','#EC5B71','Red','Magenta','DarkOrchid','White']
colorstblack=['White','Black','Black','Black','Black', 'Black','Black','Black',\
  'Black','Black','Black','Black','Black','Black','Black','Black']

# --- Plot analyzed 95th percentile and analog locations on analysis grid

fig1 = plt.figure(figsize=(9.,5.7))

ax = fig1.add_axes([0.03,.12,0.94,.83])
ctitle = 'Supplemental locations and 95th percentile of analyses, '+\
   cleadb+' to '+cleade+'-h forecast, '+cmonth
ax.set_title(ctitle,fontsize=15)
m = Basemap(projection='mill',llcrnrlon=rlonsa_t[0,0],llcrnrlat=rlatsa_t[0,0],\
        urcrnrlon=rlonsa_t[-1,-1],urcrnrlat=50.,resolution='l')  # resolution next step = i
x,y = m(lonsa,latsa)
#pfquant_m = ma.array(pfquant,mask=1-conusmask_t)
clevs = [0.,1,2,4,7,10,15,20,25,30,40,50,60]
#CS1 = m.contour(x,y,paquant_m,clevs,colors=colorstblack,cmap=None,linewidths=0.5)
CS2 = m.contourf(x,y,paquant_m,clevs,colors=colorst,cmap=None,extend='max')

#print 'jber_a, iber_a= ',jber_a, iber_a
#print 'xlocationa_t[:,jber_a,iber_a]=',xlocationa_t[:,jber_a,iber_a]
#print 'ylocationa_t[:,jber_a,iber_a]=',ylocationa_t[:,jber_a,iber_a]

#print 'xlocationa_t[:,jbou_a,ibou_a]=',xlocationa_t[:,jbou_a,ibou_a]
#print 'ylocationa_t[:,jbou_a,ibou_a]=',ylocationa_t[:,jbou_a,ibou_a]

istat = add_analog_locns(ycoord=jber_a,xcoord=iber_a,\
    ylocs=ylocationa_t[:,jber_a,iber_a],xlocs=xlocationa_t[:,jber_a,iber_a],\
    lons=rlonsa_t,lats=rlatsa_t,map=m,symbolshape='o',color='Black') # Berkeley
istat = add_analog_locns(ycoord=jpor_a,xcoord=ipor_a,\
    ylocs=ylocationa_t[:,jpor_a,ipor_a],xlocs=xlocationa_t[:,jpor_a,ipor_a],\
    lons=rlonsa_t,lats=rlatsa_t,map=m,symbolshape='v',color='Black') # Portland
istat = add_analog_locns(ycoord=jboz_a,xcoord=iboz_a,\
    ylocs=ylocationa_t[:,jboz_a,iboz_a],xlocs=xlocationa_t[:,jboz_a,iboz_a],\
    lons=rlonsa_t,lats=rlatsa_t,map=m,symbolshape='+',color='Black') # Bozeman
istat = add_analog_locns(ycoord=jbou_a,xcoord=ibou_a,\
    ylocs=ylocationa_t[:,jbou_a,ibou_a],xlocs=xlocationa_t[:,jbou_a,ibou_a],\
    lons=rlonsa_t,lats=rlatsa_t,map=m,symbolshape='p',color='Black') # Boulder
istat = add_analog_locns(ycoord=jphx_a,xcoord=iphx_a,\
    ylocs=ylocationa_t[:,jphx_a,iphx_a],xlocs=xlocationa_t[:,jphx_a,iphx_a],\
    lons=rlonsa_t,lats=rlatsa_t,map=m,symbolshape='s',color='Black') # Phoenix
istat = add_analog_locns(ycoord=jlub_a,xcoord=ilub_a,\
    ylocs=ylocationa_t[:,jlub_a,ilub_a],xlocs=xlocationa_t[:,jlub_a,ilub_a],\
    lons=rlonsa_t,lats=rlatsa_t,map=m,symbolshape='D',color='Black') # Lubbock
istat = add_analog_locns(ycoord=joma_a,xcoord=ioma_a,\
    ylocs=ylocationa_t[:,joma_a,ioma_a],xlocs=xlocationa_t[:,joma_a,ioma_a],\
    lons=rlonsa_t,lats=rlatsa_t,map=m,symbolshape='*',color='Black') # Omaha
istat = add_analog_locns(ycoord=jbis_a,xcoord=ibis_a,\
    ylocs=ylocationa_t[:,jbis_a,ibis_a],xlocs=xlocationa_t[:,jbis_a,ibis_a],\
    lons=rlonsa_t,lats=rlatsa_t,map=m,symbolshape='^',color='Black') # Bismarck
istat = add_analog_locns(ycoord=jlit_a,xcoord=ilit_a,\
    ylocs=ylocationa_t[:,jlit_a,ilit_a],xlocs=xlocationa_t[:,jlit_a,ilit_a],\
    lons=rlonsa_t,lats=rlatsa_t,map=m,symbolshape='|',color='Black') # Little Rock
istat = add_analog_locns(ycoord=jchi_a,xcoord=ichi_a,\
    ylocs=ylocationa_t[:,jchi_a,ichi_a],xlocs=xlocationa_t[:,jchi_a,ichi_a],\
    lons=rlonsa_t,lats=rlatsa_t,map=m,symbolshape='x',color='Black') # Chicago
istat = add_analog_locns(ycoord=jcin_a,xcoord=icin_a,\
    ylocs=ylocationa_t[:,jcin_a,icin_a],xlocs=xlocationa_t[:,jcin_a,icin_a],\
    lons=rlonsa_t,lats=rlatsa_t,map=m,symbolshape='>',color='Black') # Cincinnati
istat = add_analog_locns(ycoord=jatl_a,xcoord=iatl_a,\
    ylocs=ylocationa_t[:,jatl_a,iatl_a],xlocs=xlocationa_t[:,jatl_a,iatl_a],\
    lons=rlonsa_t,lats=rlatsa_t,map=m,symbolshape='+',color='Black') # Atlanta
istat = add_analog_locns(ycoord=jnyc_a,xcoord=inyc_a,\
    ylocs=ylocationa_t[:,jnyc_a,inyc_a],xlocs=xlocationa_t[:,jnyc_a,inyc_a],\
    lons=rlonsa_t,lats=rlatsa_t,map=m,symbolshape='o',color='Black') # NYC

m.drawcoastlines(linewidth=.5)
m.drawstates(linewidth=.5)
m.drawcountries(linewidth=.5)
cax = fig1.add_axes([0.1,0.08,0.8,0.03])
print 'clevs = ',clevs
cbar = fig1.colorbar(CS2,extend='neither', \
   orientation='horizontal',cax=cax,drawedges=True,ticks=clevs,format='%g')
cax.set_xlabel('Forecast precipitation amount, 95th percentile (mm)')
plot_title = 'supp_locns_ppnfcst95_ccpagrid_'+cmonth+'.pdf'
fig1.savefig(plot_title)
print 'saving plot to file = ',plot_title
print 'Plot done' 

sys.exit()

