""" this python script assumes that you have already generated analog
forecasts and raw ensemble forecasts with the appropriate fortran programs
and now wish to plot skill scores as a function of lead time and month of 
the year.  Assuming that data has been generated, it will make plots for
the events > 1 mm, > 10 mm, > 25 mm, and > q95 (95th percentile of 
climatological distribution).

You may care about different lead times, or have different methodologies
to compare, or care about different event thresholds.  In such cases, 
you should be able to adapt this software accordingly.

Lead times for which you should have calculated forecasts for this 
script to work are:  012-024, 036-048, 060-072, 084-096, 108-120, 132-144, 
156-168, and 180-192.

Coded by Tom Hamill, tom.hamill@noaa.gov, latest version April 2016
"""
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np
from numpy import ma
import os, sys
from read_and_process_reliability_bssinfo_rawens \
    import read_and_process_reliability_bssinfo_rawens
from read_and_process_reliability_bssinfo_supp \
    import read_and_process_reliability_bssinfo_supp
    
rcParams['legend.fontsize']='x-small'
rcParams['legend.fancybox']=True
rcParams['xtick.labelsize']='small'
rcParams['ytick.labelsize']='small'

nmonths = 12
nclasses = 21
nxa = 464 # CCPA 1/8 degree grid over CONUS
nya = 224 # 

output_data_directory = '/Users/thamill/refcst2/test/git' # where your post-processing algorithms 
  # stashed their data

# ---- read in brier skill score information and reliability/frequency of usage
#      for a variety of lead times and each of the 12 months

bssyearly = 0.0
bssmonthly = np.zeros((12),dtype=np.float)
bssmap     = np.zeros((nxa,nya),dtype=np.float)
rlons_anal = np.zeros((nxa,nya),dtype=np.float)
rlats_anal = np.zeros((nxa,nya),dtype=np.float)
iccpa_mask = np.zeros((nxa,nya),dtype=np.int16)
frequse    = np.zeros((nclasses),dtype=np.float)
relia      = np.zeros((nclasses),dtype=np.float)

# ---- now read in the brier skill scores for original rank analog procedure

cleadb_arr = ['012','036','060','084','108','132','156','180']
cleade_arr = ['024','048','072','096','120','144','168','192']
cthresh_arr = ['1mm','10mm','25mm','q95']
cmonths = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
#cthresh_arr = ['POP']
csupp = '10'
csuppm1 = '09'

for cthresh in cthresh_arr:
    bssmonthly_raw = np.zeros((12,8),dtype=np.float)
    bssmonthly_analogsupp = np.zeros((12,8),dtype=np.float)
    bssmonthly_analogsupp = np.zeros((12,8),dtype=np.float)
    for cleadb, cleade, ileadb in zip(cleadb_arr,cleade_arr,range(8)):

        bssyearly, bssmonthly, bssmap, relia, frequse, rlons_anal, rlats_anal, iccpa_mask = \
            read_and_process_reliability_bssinfo_supp(cleadb, cleade, cthresh, \
            csupp, output_data_directory, nxa, nya)
        bssmonthly_analogsupp[:,ileadb] = bssmonthly[:]
        
        bssyearly, bssmonthly, bssmap, relia, frequse, rlons_anal, rlats_anal, iccpa_mask = \
            read_and_process_reliability_bssinfo_rawens(cleadb, cleade, cthresh, \
            output_data_directory, nxa, nya)
        bssmonthly_raw[:,ileadb] = bssmonthly[:]

    fig = plt.figure(figsize=(9.,5.))
    plt.suptitle('Brier skill scores and differences, > '+cthresh,fontsize=16)

    # --- lefthand top plot, analog with supplemental

    a2 = fig.add_axes([.07,.11,.24,.74])
    a2.set_title('(a) Analog,\n'+csuppm1+' supplemental locations',fontsize=12)
    a2.plot(range(12),bssmonthly_analogsupp[:,0],'o-',color='black',label='Day 0.5-1.0',markersize=5.0)
    a2.plot(range(12),bssmonthly_analogsupp[:,1],'^-',color='LimeGreen',label='Day 1.5-2.0',markersize=5.0)
    a2.plot(range(12),bssmonthly_analogsupp[:,2],'s-',color='red',label='Day 2.5-3.0',markersize=5.0)
    a2.plot(range(12),bssmonthly_analogsupp[:,3],'D-',color='DodgerBlue',label='Day 3.5-4.0',markersize=5.0)
    a2.plot(range(12),bssmonthly_analogsupp[:,4],'v-',color='Orchid',label='Day 4.5-5.0',markersize=5.0)
    a2.plot(range(12),bssmonthly_analogsupp[:,5],'o-',color='orange',label='Day 5.5-6.0',markersize=5.0)
    a2.plot(range(12),bssmonthly_analogsupp[:,6],'^-',color='silver',label='Day 6.5-7.0',markersize=5.0)
    a2.plot(range(12),bssmonthly_analogsupp[:,7],'s-',color='khaki',label='Day 7.5-8.0',markersize=5.0)
    a2.plot([-0.5,11.5],[0,0],'k--',linewidth=2)
    
    a2.set_ylabel('Brier Skill Score',fontsize=11)
    a2.set_xlabel('Month',fontsize=11)
    a2.set_ylim(-0.25,0.55)
    a2.set_xlim(-0.5,11.5)
    a2.set_xticks(range(12))
    a2.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
    a2.grid (True)

    # --- middle plot, raw ensemble

    a2 = fig.add_axes([.4,.11,.24,.74])
    a2.set_title('(b) Raw ensemble',fontsize=12)
    a2.plot(range(12),bssmonthly_raw[:,0],'o-',color='black',label='Day 0.5-1.0',markersize=5.0)
    a2.plot(range(12),bssmonthly_raw[:,1],'^-',color='LimeGreen',label='Day 1.5-2.0',markersize=5.0)
    a2.plot(range(12),bssmonthly_raw[:,2],'s-',color='red',label='Day 2.5-3.0',markersize=5.0)
    a2.plot(range(12),bssmonthly_raw[:,3],'D-',color='DodgerBlue',label='Day 3.5-4.0',markersize=5.0)
    a2.plot(range(12),bssmonthly_raw[:,4],'v-',color='Orchid',label='Day 4.5-5.0',markersize=5.0)
    a2.plot(range(12),bssmonthly_raw[:,5],'o-',color='orange',label='Day 5.5-6.0',markersize=5.0)
    a2.plot(range(12),bssmonthly_raw[:,6],'^-',color='silver',label='Day 6.5-7.0',markersize=5.0)
    a2.plot(range(12),bssmonthly_raw[:,7],'s-',color='khaki',label='Day 7.5-8.0',markersize=5.0)
    a2.plot([-0.5,11.5],[0,0],'k--',linewidth=2)

    a2.set_ylabel('Brier Skill Score',fontsize=11)
    a2.set_xlabel('Month',fontsize=11)
    a2.set_ylim(-0.25,0.55)
    a2.set_xlim(-0.5,11.5)
    a2.set_xticks(range(12))
    a2.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
    a2.grid (True)
    a2.legend(loc=9,markerscale=1.0,ncol=1)

    # --- right plot, Skill difference, analog - raw

    a2 = fig.add_axes([.72,.11,.24,.74])
    a2.set_title('(c) Skill difference,\nanalog - raw',fontsize=12)
    a2.plot(range(12),bssmonthly_analogsupp[:,0]-bssmonthly_raw[:,0],\
        'o-',color='black',label='Day 0.5-1.0',markersize=5.0)
    a2.plot(range(12),bssmonthly_analogsupp[:,1]-bssmonthly_raw[:,1],\
        '^-',color='LimeGreen',label='Day 1.5-2.0',markersize=5.0)
    a2.plot(range(12),bssmonthly_analogsupp[:,2]-bssmonthly_raw[:,2],\
        's-',color='red',label='Day 2.5-3.0',markersize=5.0)
    a2.plot(range(12),bssmonthly_analogsupp[:,3]-bssmonthly_raw[:,3],\
        'D-',color='DodgerBlue',label='Day 3.5-4.0',markersize=5.0)
    a2.plot(range(12),bssmonthly_analogsupp[:,4]-bssmonthly_raw[:,4],\
        'v-',color='Orchid',label='Day 4.5-5.0',markersize=5.0)
    a2.plot(range(12),bssmonthly_analogsupp[:,5]-bssmonthly_raw[:,5],\
        'o-',color='orange',label='Day 5.5-6.0',markersize=5.0)
    a2.plot(range(12),bssmonthly_analogsupp[:,6]-bssmonthly_raw[:,6],\
        '^-',color='silver',label='Day 6.5-7.0',markersize=5.0)
    a2.plot(range(12),bssmonthly_analogsupp[:,7]-bssmonthly_raw[:,7],\
        's-',color='khaki',label='Day 7.5-8.0',markersize=5.0)
    a2.plot([-0.5,11.5],[0,0],'k--',linewidth=2)
    
    a2.set_ylabel('Brier Skill Score difference',fontsize=11)
    a2.set_xlabel('Month',fontsize=11)
    a2.set_ylim(-0.05,0.5)
    a2.set_xlim(-0.5,11.5)
    a2.set_xticks(range(12))
    a2.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
    a2.grid (True)

    plot_title = 'BSS_monthly_ccpa_nsupp'+csupp+'_'+cthresh+'.pdf'
    print 'saving plot to file = ',plot_title
    plt.savefig(plot_title)
    print 'Plot done'

