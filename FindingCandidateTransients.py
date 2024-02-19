import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab

import numpy as np 
import pandas as pd 
import myplotting as myplt
import functions as func
import scipy.optimize as spopt
import scipy.stats as spstat

# Inputs to edit
names = ['GRBdata'] #Insert here the name of the .csv files that need to be analysed
p = 0.99 #Inserting here the percentage with which the source should be classified as inlier

##############################

# Reading the data files and creating the arrays to be used
data = pd.read_csv(names[0]+'.csv')
names.pop(0)
for name in names:
    tmp_data = pd.read_csv(name+'.csv')
    data = data.append(tmp_data, ignore_index=True)

# Filter out all sources with zero or negative variability parameters
data = data.loc[ (data['V']>0.) & (data['eta']>0.)]

# Create new columns in dataframe with log10 values
data['logEta'] = data.apply(lambda row: np.log10(row.eta), axis=1)
data['logV'] = data.apply(lambda row: np.log10(row.V), axis=1)
data['logFlux'] = data.apply(lambda row: np.log10(row.maxFlx), axis=1)

freq = data.freq.unique()  #Keeping track of the frequencies used eliminating duplicates

print('Number of sources analysed:', len(data))

fig,ax1,ax2=myplt.EtaVscatter(data,freq,'EtavsVUn')

#Finding the line that best represents the two parameters of the data
best_params_eta, ml_cfcovar_linear = spopt.curve_fit(func.LinearFit, data.logFlux, data.logEta)
best_params_V, ml_cfcovar_linear = spopt.curve_fit(func.LinearFit, data.logFlux, data.logV)

m_eta, q_eta = best_params_eta[0],best_params_eta[1]
m_V, q_V = best_params_V[0],best_params_V[1]

y_eta = func.LinearFit(data.logFlux,m_eta,q_eta)
y_V = func.LinearFit(data.logFlux,m_V,q_V)

ax1.plot(data.logFlux,y_eta,label='Best fit',color='gray', ls='--')
ax2.plot(data.logFlux,y_V,label='Best fit',color='gray',ls='--')
ax1.set_title(r'Linear Fit Stable and unstable sources',fontsize=20)
ax1.legend(fontsize=15,markerscale=1.5)
ax2.legend(fontsize=15,markerscale=1.5)
plt.show()
plt.savefig('ParametersLinearFit.png')
plt.close()


#Calculating the paramater distance from the line that has been found earlier

data['distsEta'] = data.apply(lambda row: func.Params_distance(row.logEta,func.LinearFit(row.logFlux,m_eta,q_eta)), axis=1)
data['distsV'] = data.apply(lambda row: func.Params_distance(row.logV,func.LinearFit(row.logFlux,m_V,q_V)), axis=1)


#Creating a matrix necessary to calculate the Gaussian distribution of the data.

ndim = 2
data_graph = np.vstack([data.distsEta,data.distsV])
data_graph = data_graph.T

mean_deta = np.mean(data.distsEta)
mean_dV = np.mean(data.distsV)
mu = [mean_deta,mean_dV]

cov_matrix = np.cov(data.distsEta,data.distsV)
likelihood = spstat.multivariate_normal.pdf(data_graph,[mean_deta,mean_dV],cov_matrix)

outliers_prob = func.Probability(data_graph,mu,cov_matrix) #Calculating the probability for every parameter of being associated to an "inlier" source

data['probability'] = 100.-outliers_prob
temp = data.sort_values('probability')
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
ax2.legend(fontsize=15,markerscale=1.5)
ax1.legend(fontsize=15,markerscale=1.5)
ax2.set_ylabel(r'$log_{10}(V_{\nu}$)',fontsize=30)
ax2.set_xlabel(r'$log_{10}(Flux) (Jy)$',fontsize=30)
ax1.tick_params(labelsize=15)
ax2.tick_params(labelsize=15)
plt.savefig('EtavsVscatterinout')

figure,ax = myplt.OutInPlot(np.array([outliers.distsEta, outliers.distsV]).T,np.array([inliers.distsEta,inliers.distsV]).T,'OutIn_Unstable')
EtaVsVout, axveta = myplt.OutInPlot(np.array([outliers.logEta, outliers.logV]).T,np.array([inliers.logEta,inliers.logV]).T,'OutInEtavsV')
axveta.set_xlabel(r'$log_{10}(\eta_{\nu})$',fontsize=30)
axveta.set_ylabel(r'$log_{10}(V_{\nu})$',fontsize=30)
plt.savefig('OutInEtavsV')

# Outputting variable candidates
print('Number of outliers : ', len(outliers))

if len(outliers) == 0:
    print('No candidate variable sources found.')
else:
    print(outliers)
    outliers.to_csv('Outliers.csv', index=False)

