import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab
from scipy.stats import norm

import numpy as np 
import pandas as pd 

import sqlalchemy
from sqlalchemy import *
from sqlalchemy.orm import relationship

from dblogin import *
import tkp.db
from tkp.db.model import Runningcatalog
from tkp.db.model import Extractedsource
from tkp.db.model import Assocxtrsource
import dbtools 
import myplotting as myplt

global db 

# Inputs to edit
database = '' 
dataset_ids = [6]
window_size = 10
cutoff_val = 3
minDatapoints = 10
sigma =3. # the number of sigma the deviation has to be above the flux uncertainties
sigmaThresh = 2

##############################

def plothist(x, threshold, filename):
# Create a histogram of the data t
    plt.hist(x,bins=50,histtype='stepfilled')
    plt.axvline(x=threshold, linewidth=2, color='k',linestyle='--')
    plt.xlabel('Deviation')
    plt.ylabel('Number of sources')
    plt.savefig(filename)
    plt.close()
    return

def checkCand(row):
    if (row.flux - row.MAvg) > row.fluxErr:
        val = 1
    else:
        val = 0
    return val

def GetAllLightcurves(session,dataset_id):
    # Returns all the light curves for the unique sources in a given dataset
    x = session.query(Runningcatalog).filter(Runningcatalog.dataset_id == dataset_id)
    dx = pd.read_sql_query(x.statement,db.connection)
    dx = dx.rename(index=str,columns={'id' : 'runcat'})
    y = session.query(Extractedsource,Assocxtrsource).select_from(join(Extractedsource,Assocxtrsource)).filter(Assocxtrsource.runcat_id.in_(dx.runcat)).all()
    data = [[y[i].Extractedsource.id, y[i].Assocxtrsource.runcat_id, y[i].Extractedsource.image.taustart_ts, y[i].Extractedsource.f_int, y[i].Extractedsource.f_int_err] for i in range(len(y))]
    data=pd.DataFrame(data=data, columns =['srcID','runcat','time','flux','fluxErr'])
    data=data.sort_values(['runcat','time'])
    return(data)

def access(engine,host,port,user,password,database):
    """ Access the database using sqlalchemy"""
    # make db global in order to be used in GetPandaExtracted
    global db
    db = tkp.db.Database(engine=engine, host=host, port=port,
                     user=user, password=password, database=database)
    db.connect()
    session = db.Session()
    print ('connected!')
    return session

def SigmaFit(data):
    median = np.median(data)
    std_median = np.sqrt(np.mean([(i-median)**2. for i in data]))
    tmp_data = [a for a in data if a < 3.*std_median+median and a > median - 3.*std_median]
    param1 = norm.fit(tmp_data)
    param2 = norm.fit(data)
    return param1, param2

##############################

# Connect to the database and run the queries
session = access(engine,host,port,user,password,database)
for dataset_id in dataset_ids:
    lightcurves = GetAllLightcurves(session,dataset_id)  # get all the lightcurves
    lightcurves.to_csv('ds'+str(dataset_id)+'_lightcurves.csv', index=False)

data = pd.read_csv('ds'+str(dataset_ids[0])+'_lightcurves.csv')
dataset_ids.pop(0)

for dataset_id in dataset_ids:
    data=pd.concat([data,pd.read_csv('ds'+str(dataset_id)+'_lightcurves.csv')], ignore_index=True)
    

finalData = pd.DataFrame(columns = ["srcID","runcat","time","flux","fluxErr","dpts","MAvg","deviation","candidate"])
runcats = data.runcat.unique()

for runcat in runcats:
    lc = data.loc[data['runcat'] == runcat]
    lc = lc.reset_index()
    lc["dpts"] = lc.index+1
    lc['MAvg'] = lc['flux'].rolling(window=window_size, min_periods=1).mean()
    std = lc.flux.std()
    lc['deviation'] = (lc['flux'] - lc['MAvg']) / std
    lc['candidate'] = np.where(np.abs(lc['flux'] - lc['MAvg']) > (sigma * lc['fluxErr']), 1, 0)
    lc = lc.fillna(0)
    lc = lc.drop('index',axis=1)
    finalData = finalData._append(lc, ignore_index=True)

# remove all rows with a zero deviation
finalData = finalData.loc[finalData['deviation'] != 0]

# Find transients that have deviation values that lie above threshold value
all_deviations = finalData.deviation
params_med = np.median(all_deviations), np.sqrt(np.mean([(i-np.median(all_deviations))**2. for i in all_deviations]))
threshold = cutoff_val * params_med[1] + params_med[0]
# Save plot of histogram and threshold
plothist(all_deviations, threshold, 'LOFAR_deviation_hist.png')
print(cutoff_val, '*', params_med[1], '+',  + params_med[0], '=', threshold)
print(params_med[0])

candidates = finalData.loc[(finalData['deviation'] > threshold) & (finalData['dpts'] >= minDatapoints) & (finalData['candidate'] == 1)]

if len(candidates) == 0:
    print('No candidate variable sources found.')
else:
    print(candidates)
    candidates.to_csv('candidates.csv', index=False)


VarParams = dbtools.GetVarParams(session,dataset_id)
plotdata = [[VarParams[i].Runningcatalog.id, VarParams[i].Varmetric.eta_int, VarParams[i].Varmetric.v_int,
                 VarParams[i].Varmetric.lightcurve_max, VarParams[i].Varmetric.lightcurve_median,
                 (VarParams[i].Varmetric.band.freq_central/1e6), VarParams[i].Runningcatalog.datapoints,
                 VarParams[i].Varmetric.newsource, VarParams[i].Runningcatalog.wm_ra, VarParams[i].Runningcatalog.wm_decl] for i in range(len(VarParams))]
plotdata = pd.DataFrame(data=plotdata,columns=['runcat','eta','V','maxFlx','avgFlx','freq','dpts','newSrc','ra','dec'])
plotdata = plotdata.fillna('N')

# Create new columns in dataframe with log10 values
plotdata['logEta'] = plotdata.apply(lambda row: np.log10(row.eta), axis=1)
plotdata['logV'] = plotdata.apply(lambda row: np.log10(row.V), axis=1)


plotdata_candidates = pd.merge(candidates,plotdata,how='inner', on='runcat')
runcats=candidates.runcat
plotdata_rest = plotdata[~plotdata.runcat.isin(runcats)]
plotdata_rest = plotdata_rest.loc[(plotdata_rest['eta'] > 0) & (plotdata_rest['V'] > 0) & (plotdata_rest['dpts']>1) & (plotdata_rest['newSrc']=='N')]

# finding the old sigma thresholds for plotting later
paramx, paramx2 = SigmaFit(plotdata_rest['logEta'])
paramy, paramy2 = SigmaFit(plotdata_rest['logV'])
sigcutx = paramx[1]*sigmaThresh+paramx[0]
sigcuty = paramy[1]*sigmaThresh+paramy[0]

EtaVsVout, axveta = myplt.OutInPlot(np.array([plotdata_candidates.logEta, plotdata_candidates.logV]).T,np.array([plotdata_rest.logEta,plotdata_rest.logV]).T,'OutInEtavsV')
axveta.axvline(x=sigcutx, linewidth=2, color='k', linestyle='--')
axveta.axhline(y=sigcuty, linewidth=2, color='k', linestyle='--')

axveta.set_xlabel(r'$log_{10}(\eta_{\nu})$',fontsize=30)
axveta.set_ylabel(r'$log_{10}(V_{\nu})$',fontsize=30)
axveta.tick_params(axis='both', which='major', labelsize=25)
plt.savefig('EtavsV_flares')


exit()
