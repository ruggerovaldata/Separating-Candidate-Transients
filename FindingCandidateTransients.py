import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab

import numpy as np 
import pandas as pd 
import myplotting as myplt
import scipy.optimize as spopt
import scipy.stats as spstat
from scipy.stats import norm

import sqlalchemy
from sqlalchemy import *
from sqlalchemy.orm import relationship
from dblogin import *
import dbtools 
import tkp.db
import os
import os.path



def LinearFit(data,a,b):
    return a * data + b

def Params_distance(y_o,y_p):
    return (y_o-y_p)

def Probability(data, mean, sigma):
    prob=[]
    for x in data: 
        dist = np.dot((x-mean).transpose(),np.linalg.inv(sigma))
        dist = np.dot(dist,(x-mean))
        prob.append(1- spstat.chi2.cdf(dist,2))
    prob = np.array(prob)        

    return prob*100

def SigmaFit(data):
    median = np.median(data)
    std_median = np.sqrt(np.mean([(i-median)**2. for i in data]))
    tmp_data = [a for a in data if a < 3.*std_median+median and a > median - 3.*std_median]
    param1 = norm.fit(tmp_data)
    param2 = norm.fit(data)
    return param1, param2


p = 0.95 #Inserting here the percentage with which the source should be classified as inlier

database = 'antoniar' 
dataset_ids = [6]
sigmaThresh = 2

global db

# Connect to the database and run the queries
session = dbtools.access(engine,host,port,user,password,database)
dataset_ids_TMP=[]
for dataset_id in dataset_ids:
    if not os.path.isfile('ds'+str(dataset_id)+'.csv'):
        VarParams = dbtools.GetVarParams(session,dataset_id)   # Get the running catalogue and varmetric catalogues and combine
        plotdata = [[VarParams[i].Runningcatalog.id, VarParams[i].Varmetric.eta_int, VarParams[i].Varmetric.v_int,
                            VarParams[i].Varmetric.lightcurve_max, VarParams[i].Varmetric.lightcurve_median,
                            (VarParams[i].Varmetric.band.freq_central/1e6), VarParams[i].Runningcatalog.datapoints,
                            VarParams[i].Varmetric.newsource, VarParams[i].Runningcatalog.wm_ra, VarParams[i].Runningcatalog.wm_decl] for i in range(len(VarParams))]
        plotdata = pd.DataFrame(data=plotdata,columns=['runcat','eta','V','maxFlx','avgFlx','freq','dpts','newSrc','ra','dec'])
        plotdata = plotdata.fillna('N')
        plotdata = plotdata.loc[(plotdata['eta'] > 0) & (plotdata['V'] > 0) & (plotdata['dpts']>1) & (plotdata['newSrc']=='N')]
        if len(plotdata) != 0:
            # Save the data for plotting
            plotdata.to_csv('ds'+str(dataset_id)+'.csv', index=False)
        else:
            print('Dataset '+str(dataset_id)+' contains no sources. Remove from list.')
            dataset_ids_TMP.append(dataset_id)

for id in dataset_ids_TMP:
    dataset_ids.remove(id)

#Loading the First file and creating the array that will be used

data = pd.read_csv('ds'+str(dataset_ids[0])+'.csv') #Loading the first file
dataset_ids.pop(0)

for name in dataset_ids:
    data=pd.concat([data,pd.read_csv('ds'+str(name)+'.csv')], ignore_index=True)

data = data.loc[ (data['V']>0.) & (data['eta']>0.)]

# Create new columns in dataframe with log10 values
data['logEta'] = data.apply(lambda row: np.log10(row.eta), axis=1)
data['logV'] = data.apply(lambda row: np.log10(row.V), axis=1)
data['logFlux'] = data.apply(lambda row: np.log10(row.maxFlx), axis=1)

freq = data.freq.unique()  #Keeping track of the frequencies used eliminating duplicates

print('Number of sources analysed:', len(data))

# finding the old sigma thresholds for plotting later
paramx, paramx2 = SigmaFit(data['logEta'])
paramy, paramy2 = SigmaFit(data['logV'])
sigcutx = paramx[1]*sigmaThresh+paramx[0]
sigcuty = paramy[1]*sigmaThresh+paramy[0]

tmp_data=data.loc[ (data['logEta'] > sigcutx) & (data['logV'] > sigcuty)]

print('2sigma outliers are:')
print(tmp_data)

fig,ax1,ax2=myplt.EtaVscatter(data,freq,'EtavsVUn')

#Finding the line that best represents the two parameters of the data
best_params_eta, ml_cfcovar_linear = spopt.curve_fit(LinearFit,data.logFlux.values, data.logEta.values)
best_params_V, ml_cfcovar_linear = spopt.curve_fit(LinearFit, data.logFlux.values, data.logV.values)

m_eta, q_eta = best_params_eta[0],best_params_eta[1]
m_V, q_V = best_params_V[0],best_params_V[1]

print('eta linear fit parameters: m = '+str(np.around(m_eta,3))+' q = '+str(np.around(q_eta,3)))
print('V linear fit parameters: : m = '+str(np.around(m_V,3))+' q = '+str(np.around(q_V,3)))

y_eta = LinearFit(data.logFlux.values,m_eta,q_eta)
y_V = LinearFit(data.logFlux.values,m_V,q_V)

ax1.plot(data.logFlux.values,y_eta,label='Best fit',color='gray', ls='--')
ax2.plot(data.logFlux.values,y_V,label='Best fit',color='gray',ls='--')
#ax1.set_title(r'Linear Fit Stable and unstable sources',fontsize=25)
ax1.legend(fontsize=25,markerscale=1.5)
ax2.legend(fontsize=25,markerscale=1.5)
ax1.tick_params(axis='both', which='major', labelsize=25)
ax2.tick_params(axis='both', which='major', labelsize=25)
plt.show()
plt.savefig('ParametersLinearFit.png')
plt.close()


#Calculating the paramater distance from the line that has been found earlier

data['distsEta'] = data.apply(lambda row: Params_distance(row.logEta,LinearFit(row.logFlux,m_eta,q_eta)), axis=1)
data['distsV'] = data.apply(lambda row: Params_distance(row.logV,LinearFit(row.logFlux,m_V,q_V)), axis=1)


#Creating a matrix necessary to calculate the Gaussian distribution of the data.
print('\n WHOLE DATASET: \n')

ndim = 2
data_graph = np.vstack([data.distsEta.values,data.distsV.values])
data_graph = data_graph.T
mean_deta = np.mean(data.distsEta.values)
mean_dV = np.mean(data.distsV.values)
mu = [mean_deta,mean_dV]

cov_matrix = np.cov(data.distsEta.values,data.distsV.values)
likelihood = spstat.multivariate_normal.pdf(data_graph,[mean_deta,mean_dV],cov_matrix)

outliers_prob = Probability(data_graph,mu,cov_matrix) #Calculating the probability for every parameter of being associated to an "inlier" source

data['probability'] = 100.-outliers_prob
data = data.sort_values('probability',ascending=False)
data.to_csv('Wholedatasetoutput.csv', index=False)

figure, axes = myplt.MyCorner(data.distsEta,data.distsV,data.probability/100.,'CornerPlot') #Printing the corner plot both with the likelihood and without the likelihood

chi2 = spstat.chi2.ppf([p],2)[0]

# finding the outliers with probabilities >99% and with positive distances abovve the trend line
inliers = data.loc[ (data['probability'] <= p*100.) | ( (data['distsEta'] < 0) | (data['distsV'] < 0))]
outliers = data.loc[ (data['probability'] > p*100.) & (data['distsEta'] > 0) & (data['distsV'] > 0)]


    

# Plotting
fig,(ax1,ax2) = plt.subplots(2,1,figsize=(14,14))

ax1.scatter(outliers.logFlux,outliers.logEta,color='red',label='Outliers')
ax1.scatter(inliers.logFlux,inliers.logEta,color='blue',label='Inliers')
ax2.scatter(outliers.logFlux,outliers.logV,color='red',label='Outliers')
ax2.scatter(inliers.logFlux,inliers.logV,color='blue',label = 'Inliers')


ax1.set_ylabel(r'$log_{10}(\eta_{\nu}$)',fontsize=30)
ax2.legend(fontsize=25,markerscale=1.5)
ax1.legend(fontsize=25,markerscale=1.5)
ax2.set_ylabel(r'$log_{10}(V_{\nu}$)',fontsize=30)
ax2.set_xlabel(r'$log_{10}(Flux) (Jy)$',fontsize=30)
ax1.tick_params(labelsize=20)
ax2.tick_params(labelsize=20)
plt.savefig('EtavsVscatterinout')



figure,ax = myplt.OutInPlot(np.array([outliers.distsEta, outliers.distsV]).T,np.array([inliers.distsEta,inliers.distsV]).T,'OutIn_Unstable')
EtaVsVout, axveta = myplt.OutInPlot(np.array([outliers.logEta, outliers.logV]).T,np.array([inliers.logEta,inliers.logV]).T,'OutInEtavsV')
axveta.axvline(x=sigcutx, linewidth=2, color='k', linestyle='--')
axveta.axhline(y=sigcuty, linewidth=2, color='k', linestyle='--')

axveta.set_xlabel(r'$log_{10}(\eta_{\nu})$',fontsize=30)
axveta.set_ylabel(r'$log_{10}(V_{\nu})$',fontsize=30)
plt.legend()
axveta.tick_params(axis='both', which='major', labelsize=25)
plt.savefig('OutInEtavsV')

# Outputting variable candidates
print('Number of outliers : ', len(outliers))

if len(outliers) == 0:
    print('No candidate variable sources found.')
else:
    outliers = outliers.sort_values('probability',ascending=False)
    print(outliers)
    outliers.to_csv('Outliers.csv', index=False)

    
print('\n DATASET ABOVE THE LINE: \n')

ndim = 2

dataBest = data.loc[(data['distsEta'] > 0) & (data['distsV'] > 0)]

data_graph = np.vstack([dataBest.distsEta,dataBest.distsV])
data_graph = data_graph.T

mean_deta = np.mean(dataBest.distsEta)
mean_dV = np.mean(dataBest.distsV)
mu = [mean_deta,mean_dV]

cov_matrix = np.cov(dataBest.distsEta,dataBest.distsV)
likelihood = spstat.multivariate_normal.pdf(data_graph,[mean_deta,mean_dV],cov_matrix)

outliers_prob = Probability(data_graph,mu,cov_matrix) #Calculating the probability for every parameter of being associated to an "inlier" source

dataBest['probability'] = 100.-outliers_prob
dataBest = dataBest.sort_values('probability',ascending=False)
dataBest.to_csv('WholedatasetoutputBest.csv', index=False)

figure, axes = myplt.MyCorner(dataBest.distsEta,dataBest.distsV,dataBest.probability/100.,'CornerPlotBest') #Printing the corner plot both with the likelihood and without the likelihood

chi2 = spstat.chi2.ppf([p],2)[0]

# finding the outliers with probabilities >99% and with positive distances abovve the trend line
inliersBest = dataBest.loc[ (dataBest['probability'] <= p*100.) ]
outliersBest = dataBest.loc[ (dataBest['probability'] > p*100.) ]

# Plotting
fig,(ax1,ax2) = plt.subplots(2,1,figsize=(14,14))

ax1.scatter(outliersBest.logFlux,outliersBest.logEta,color='red',label='Outliers')
ax1.scatter(inliersBest.logFlux,inliersBest.logEta,color='blue',label='Inliers')
ax2.scatter(outliersBest.logFlux,outliersBest.logV,color='red',label='Outliers')
ax2.scatter(inliersBest.logFlux,inliersBest.logV,color='blue',label = 'Inliers')

ax1.set_ylabel(r'$log_{10}(\eta_{\nu}$)',fontsize=30)
ax2.legend(fontsize=25,markerscale=1.5)
ax1.legend(fontsize=25,markerscale=1.5)
ax2.set_ylabel(r'$log_{10}(V_{\nu}$)',fontsize=30)
ax2.set_xlabel(r'$log_{10}(Flux) (Jy)$',fontsize=30)
ax1.tick_params(labelsize=20)
ax2.tick_params(labelsize=20)
plt.savefig('EtavsVscatterinoutBest')

figure,ax = myplt.OutInPlot(np.array([outliersBest.distsEta, outliersBest.distsV]).T,np.array([inliersBest.distsEta,inliersBest.distsV]).T,'OutIn_UnstableBest')
EtaVsVout, axveta = myplt.OutInPlot(np.array([outliersBest.logEta, outliersBest.logV]).T,np.array([inliersBest.logEta,inliersBest.logV]).T,'OutInEtavsVBest')
axveta.axvline(x=sigcutx, linewidth=2, color='k', linestyle='--')
axveta.axhline(y=sigcuty, linewidth=2, color='k', linestyle='--')
axveta.set_xlabel(r'$log_{10}(\eta_{\nu})$',fontsize=30)
axveta.set_ylabel(r'$log_{10}(V_{\nu})$',fontsize=30)
axveta.tick_params(axis='both', which='major', labelsize=25)
plt.savefig('OutInEtavsVBest')

# Outputting variable candidates
print('Number of outliers : ', len(outliersBest))

if len(outliersBest) == 0:
    print('No candidate variable sources found.')
else:
    outliersBest = outliersBest.sort_values('probability',ascending=False)
    print(outliersBest)
    outliersBest.to_csv('OutliersBest.csv', index=False)
