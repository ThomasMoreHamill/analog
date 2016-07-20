""" this routine generates reliability diagrams for the 
    analog forecasts, along with the reliability diagrams
    for the raw ensemble.
    
    coded by Tom Hamill, tom.hamill@noaa.gov, Apr 2016
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
from fortran_routine4 import read_conusmask_ccpa
rcParams['legend.fontsize']='medium'
rcParams['legend.fancybox']=True
rcParams['xtick.labelsize']='medium'

cleadb = sys.argv[1] # beginning of period's lead time in hours, e.g., '096'
                     # (don't include '' marks)
cleade = sys.argv[2] # ending of periods lead time in hours, e.g., '108'
cthresh = sys.argv[3]
# acceptable cthresh values are 'POP','1mm','2p5mm','5mm','10mm','25mm','50mm'
#   'q50','q67','q80','q90','q95'
csupp = '20'

output_data_directory = '/Rf2_tests/ccpa'
nmonths = 12
nclasses = 21
nxa = 464 # CCPA 1/8 degree grid over CONUS
nya = 224 # 

# ---- read in brier skill score information and reliability/frequency of usage
#      for analog post-processed

bssyearly = 0.0
bssmonthly = np.zeros((12),dtype=np.float)
bssmap     = np.zeros((nxa,nya),dtype=np.float)
rlons_anal = np.zeros((nxa,nya),dtype=np.float)
rlats_anal = np.zeros((nxa,nya),dtype=np.float)
iccpa_mask = np.zeros((nxa,nya),dtype=np.int16)
frequse    = np.zeros((nclasses),dtype=np.float)
relia      = np.zeros((nclasses),dtype=np.float)

bssyearly, bssmonthly, bssmap, relia, frequse, rlons_anal, rlats_anal, iccpa_mask = \
    read_and_process_reliability_bssinfo(cleadb, cleade, cthresh, \
    output_data_directory, nxa, nya)

relia_m = ma.array(relia)
relia_m = ma.masked_where(relia_m < 0.0, relia_m)
bssmap_t = np.transpose(bssmap)

# ---- read in Brier skill score and reliability / frequency for supplemental analog

bssyearly_supp = 0.0
bssmonthly_supp = np.zeros((12),dtype=np.float)
bssmap_supp     = np.zeros((nxa,nya),dtype=np.float)
frequse_supp    = np.zeros((nclasses),dtype=np.float)
relia_supp      = np.zeros((nclasses),dtype=np.float)

bssyearly_supp, bssmonthly_supp, bssmap_supp, relia_supp, frequse_supp, \
   rlons_anal, rlats_anal, iccpa_mask = read_and_process_reliability_bssinfo_supp \
   (cleadb, cleade, cthresh, csupp, output_data_directory, nxa, nya)

relia_supp_m = ma.array(relia_supp)
relia_supp_m = ma.masked_where(relia_supp_m < 0.0, relia_supp_m)
bssmap_supp_t = np.transpose(bssmap_supp)

# ---- read in brier skill score information and reliability/frequency of usage
#      for raw ensemble

bssyearly_raw = 0.0
bssmonthly_raw = np.zeros((12),dtype=np.float)
bssmap_raw     = np.zeros((nxa,nya),dtype=np.float)
frequse_raw    = np.zeros((nclasses),dtype=np.float)
relia_raw      = np.zeros((nclasses),dtype=np.float)

bssyearly_raw, bssmonthly_raw, bssmap_raw, relia_raw, frequse_raw, rlons_anal, rlats_anal, iccpa_mask = \
    read_and_process_reliability_bssinfo_rawens(cleadb, cleade, cthresh, nxa, nya)

# --- only 11 members with raw data, so 12 relia diagram ticks; map them over from the original 21

relia_raw_12 = np.zeros((12),dtype=np.float)
frequse_raw_12 = np.zeros((12),dtype=np.float)
relia_raw_12[0] = relia_raw[0]
relia_raw_12[1] = relia_raw[2]
relia_raw_12[2] = relia_raw[4]
relia_raw_12[3] = relia_raw[5]
relia_raw_12[4] = relia_raw[7]
relia_raw_12[5] = relia_raw[9]
relia_raw_12[6] = relia_raw[11]
relia_raw_12[7] = relia_raw[13]
relia_raw_12[8] = relia_raw[15]
relia_raw_12[9] = relia_raw[16]
relia_raw_12[10] = relia_raw[18]
relia_raw_12[11] = relia_raw[20]

frequse_raw_12[0] = frequse_raw[0]
frequse_raw_12[1] = frequse_raw[2]
frequse_raw_12[2] = frequse_raw[4]
frequse_raw_12[3] = frequse_raw[5]
frequse_raw_12[4] = frequse_raw[7]
frequse_raw_12[5] = frequse_raw[9]
frequse_raw_12[6] = frequse_raw[11]
frequse_raw_12[7] = frequse_raw[13]
frequse_raw_12[8] = frequse_raw[15]
frequse_raw_12[9] = frequse_raw[16]
frequse_raw_12[10] = frequse_raw[18]
frequse_raw_12[11] = frequse_raw[20]

print 'relia_raw_12 = ',relia_raw_12
#sys.exit()

relia_raw_m = ma.array(relia_raw_12)
relia_raw_m = ma.masked_where(relia_raw_m < 0.0, relia_raw_m)
bssmap_raw_t = np.transpose(bssmap_raw)


# -- make reliability diagram (left: analog, right: raw)

fig = plt.figure(figsize=(9.,5.))
plt.suptitle('Reliability for '+cleadb+'-'+cleade+'-h forecasts, > '+cthresh,fontsize=23)

##### PANEL 1 #####

rcParams['xtick.labelsize']='small'
rcParams['ytick.labelsize']='small'
a3 = fig.add_axes([.09,.1,.38,.75])
a3.set_title('(a) Analog, '+csupp+' supplemental locations',fontsize=13)

# --- add basic reliability diagram and 5/95 confidence intervals

probs = np.arange(nclasses) * 100./np.real(nclasses-1)
a3.plot(probs,100.*relia_supp_m,'o-',color='r')
a3.plot([0,100],[0,100],'--',color='k')
a3.set_ylabel('Observed Relative Frequency (%)',fontsize=11)
a3.set_xlabel('Forecast Probability (%)',fontsize=11)
a3.set_ylim(-1,101)
a3.set_xlim(-1,101)

# -- BSS inserted here

strbss = 'BSS = %0.2f' %(bssyearly_supp)
a3.text(20,6,strbss,fontsize=13)
rcParams['xtick.labelsize']='x-small'
rcParams['ytick.labelsize']='x-small'

# --- Frequency of usage inset diagram

a4 = fig.add_axes([.13,.59,.15,.19])
a4.bar(probs,frequse_supp,width=5,bottom=0.00001,\
    log=True,color='red',edgecolor='black',align='center')
a4.set_xlim(-3,103)
#a2.set_ylim(.00001,1.)
a4.set_ylim(0.0,1.)
a4.set_title('Frequency of usage',fontsize=8)
a4.set_xlabel('Probability',fontsize=7)
a4.set_ylabel('Frequency',fontsize=7)
a4.hlines([.0001,.001,.01,.1],0,100,linestyles='dashed',colors='black')


##### PANEL 2 #####

rcParams['xtick.labelsize']='small'
rcParams['ytick.labelsize']='small'
a5 = fig.add_axes([.72,.1,.27,.75])
a5.set_title('(b) Raw 11-member ensemble',fontsize=13)

# --- add basic reliability diagram and 5/95 confidence intervals

nclasses12 = 12
probs = np.arange(nclasses12) * 100./np.real(nclasses12-1)
a5.plot(probs,100.*relia_raw_m,'o-',color='r')
a5.plot([0,100],[0,100],'--',color='k')
a5.set_ylabel('Observed Relative Frequency (%)',fontsize=11)
a5.set_xlabel('Forecast Probability (%)',fontsize=11)
a5.set_ylim(-1,101)
a5.set_xlim(-1,101)

# -- BSS inserted here

strbss = 'BSS = %0.2f' %(bssyearly_raw)
a5.text(20,6,strbss,fontsize=13)
rcParams['xtick.labelsize']='x-small'
rcParams['ytick.labelsize']='x-small'

# --- Frequency of usage inset diagram

a6 = fig.add_axes([.63,.59,.15,.19])
a6.bar(probs,frequse_raw_12,width=5*21./12.,bottom=0.00001,\
    log=True,color='red',edgecolor='black',align='center')
a6.set_xlim(-3,103)
#a2.set_ylim(.00001,1.)
a6.set_ylim(0.0,1.)
a6.set_title('Frequency of usage',fontsize=8)
a6.set_xlabel('Probability',fontsize=7)
a6.set_ylabel('Frequency',fontsize=7)
a6.hlines([.0001,.001,.01,.1],0,100,linestyles='dashed',colors='black')

plot_title = 'relia_rankan_supp_rawens_ccpa_nsupp'+csupp+'_'+\
    cleadb+'_to_'+cleade+'h_'+cthresh+'.pdf'

print 'saving plot to file = ',plot_title
plt.savefig(plot_title)
print 'Plot done'

