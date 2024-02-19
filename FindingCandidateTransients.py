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


names = ['GRBdata'] #Insert here the name of the .csv files that need to be analysed
p = 0.99 #Inserting here the percentage with which the source should be classified as inlier

# Reading the data files and creating the arrays to be used
data = pd.read_csv(names[0]+'.csv')
names.pop(0)
for name in names:
    tmp_data = pd.read_csv(name+'.csv')
    data = data.append(tmp_data, ignore_index=True)

data = data.loc[ (data['V']>0.) & (data['eta']>0.)]

idx = np.array((np.log10(data['eta']),np.log10(data['V']),data['runcat'],np.log10(data['maxFlx']))) #Connecting the data to their running catalog indexes
eta = np.array((data['eta'],data['freq'])) #Collecting eta parameter of the data connected to the frequency
V = np.array((data['V'],data['freq'])) #Same as for eta
ra = np.array(data['ra'])#Getting the positions of the sources
dec = np.array(data['dec'])
max_flux = np.array((data['maxFlx'],data['freq'])) #Same as V and eta

freq = data.freq.unique()  #Keeping track of the frequencies used eliminating duplicates


#Loading the First file and creating the array that will be used

#data = pd.read_csv(names[0]+'.csv') #Loading the first file
#names.pop(0)
#data = data.loc[ (data['V']>0.) & (data['eta']>0.)]

#freq=[]
#freq.append(data['freq'][2])


#idx = np.array((np.log10(data['eta']),np.log10(data['V']),data['runcat'],np.log10(data['maxFlx']))) #Connecting the data to their running catalog indexes
#eta = np.array((data['eta'],data['freq'])) #Collecting eta parameter of the data connected to the frequency
#V = np.array((data['V'],data['freq'])) #Same as for eta

#ra = np.array(data['ra'])#Getting the positions of the sources
#dec = np.array(data['dec'])

#max_flux = np.array((data['maxFlx'],data['freq'])) #Same as V and eta
#i=1

#print('Sources in the file number ', i, ' : ', eta.shape[-1])
# Loading the rest of the files
#i=2
#for name in names: 
#    temp_eta, temp_V, temp_flux, temp_idx = func.Data_Load(name+'.csv')
#    print('Sources in the file number ', i, ' : ', temp_eta.shape[-1])
#    freq.append(temp_eta[1][0])
#    eta = np.concatenate((eta,temp_eta),axis=1)
#    V = np.concatenate((V,temp_V),axis=1)
#    max_flux = np.concatenate((max_flux,temp_flux),axis=1)
#    idx= np.concatenate((idx,temp_idx),axis=1)
#    i+=1

#freq=np.array(freq)
#freq=np.unique(freq) #Keeping track of the frequencies used eliminating duplicates


eta_log = np.array((np.log10(eta[0]),eta[1]))
V_log = np.array((np.log10(V[0]),V[1]))
flux_log = np.array((np.log10(max_flux[0]),max_flux[1]))

print('Number of sources analysed:', eta_log.shape[-1])

fig,ax1,ax2=myplt.EtaVscatter(eta_log,V_log,flux_log,freq,'EtavsVUn')

#Finding the line that best represents the two parameters of the data

best_params_eta, ml_cfcovar_linear = spopt.curve_fit(func.LinearFit, flux_log[0], eta_log[0])
best_params_V, ml_cfcovar_linear = spopt.curve_fit(func.LinearFit, flux_log[0], V_log[0])

m_eta, q_eta = best_params_eta[0],best_params_eta[1]
m_V, q_V = best_params_V[0],best_params_V[1]

y_eta = func.LinearFit(flux_log[0],m_eta,q_eta)
y_V = func.LinearFit(flux_log[0],m_V,q_V)

ax1.plot(flux_log[0],y_eta,label='Best fit',color='gray', ls='--')
ax2.plot(flux_log[0],y_V,label='Best fit',color='gray',ls='--')
ax1.set_title(r'Linear Fit Stable and unstable sources',fontsize=20)
ax1.legend(fontsize=15,markerscale=1.5)
ax2.legend(fontsize=15,markerscale=1.5)
plt.show()
plt.savefig('ParametersLinearFit.png')
plt.close()

#Calculating the paramater distance from the line that has been found earlier

dists_eta = func.Params_distance(eta_log[0],y_eta)
dists_V = func.Params_distance(V_log[0],y_V)

#Creating a matrix necessary to calculate the Gaussian distribution of the data.

ndim = 2
data_graph = np.vstack([dists_eta,dists_V])
data_graph = data_graph.T

mean_deta = np.mean(dists_eta)
mean_dV = np.mean(dists_V)
mu = [mean_deta,mean_dV]

cov_matrix = np.cov(dists_eta,dists_V)
likelihood = spstat.multivariate_normal.pdf(data_graph,[mean_deta,mean_dV],cov_matrix)

outliers_prob = func.Probability(data_graph,mu,cov_matrix) #Calculating the probability for every parameter of being associated to an "inlier" source

func.GetOutput(idx,ra,dec,eta,V,dists_eta,dists_V, outliers_prob,'Wholedatasetoutput') #Obtaining the output file

figure, axes = myplt.MyCorner(dists_eta,dists_V,likelihood,'CornerPlot') #Printing the corner plot both with the likelihood and without the likelihood


chi2 = spstat.chi2.ppf([p],2)

inliers = []
outliers = []

for x in data_graph: 
    if func.IntervalCalc(x,mu,cov_matrix) <= chi2:
        inliers.append(x)
    else: 
        outliers.append(x)

inliers = np.array(inliers)
outliers = np.array(outliers)


print('\n WHOLE DATASET: \n')
outliers_whole = func.FindOutliersIdx(outliers,eta_log, V_log, idx, y_eta, y_V)
outliers_eta=[]
outliers_V=[]
inliers_eta=[]
inliers_V=[]
flux_outliers=[]
flux_inliers=[]

#The next few lines are made to reconnect the distances to the original parameters value so that it is possible to obtain 
#the orgiinal graph with the prediction

for i,val in enumerate(idx[-2]):
    if val in outliers_whole:
        outliers_eta.append(idx[0][i])
        outliers_V.append(idx[1][i])
        flux_outliers.append(idx[-1][i])
    else:
        inliers_eta.append(idx[0][i])
        inliers_V.append(idx[1][i])
        flux_inliers.append(idx[-1][i])

outliers_EtavsV = np.array((outliers_eta, outliers_V))

inliers_EtavsV= np.array((inliers_eta,inliers_V))

flux_outliers=np.array(flux_outliers)
flux_inliers=np.array(flux_inliers)

outliers_EtavsV = outliers_EtavsV.T
inliers_EtavsV =inliers_EtavsV.T

fig,(ax1,ax2) = plt.subplots(2,1,figsize=(14,14))
ax1.scatter(flux_outliers,outliers_EtavsV[:,0],color='red',label='Outliers')
ax1.scatter(flux_inliers,inliers_EtavsV[:,0],color='blue',label='Inliers')
ax2.scatter(flux_outliers,outliers_EtavsV[:,1],color='red',label='Outliers')
ax2.scatter(flux_inliers,inliers_EtavsV[:,1],color='blue',label = 'Inliers')
ax1.set_ylabel(r'$log_{10}(\eta_{\nu}$)',fontsize=30)
ax2.legend(fontsize=15,markerscale=1.5)
ax1.legend(fontsize=15,markerscale=1.5)
ax2.set_ylabel(r'$log_{10}(V_{\nu}$)',fontsize=30)
ax2.set_xlabel(r'$log_{10}(Flux) (Jy)$',fontsize=30)
ax1.tick_params(labelsize=15)
ax2.tick_params(labelsize=15)
plt.savefig('EtavsVscatterinout')

figure,ax = myplt.OutInPlot(outliers,inliers,'OutIn_Unstable')
EtaVsVout, axveta = myplt.OutInPlot(outliers_EtavsV,inliers_EtavsV,'OutInEtavsV')
axveta.set_xlabel(r'$log_{10}(\eta_{\nu})$',fontsize=30)
axveta.set_ylabel(r'$log_{10}(V_{\nu})$',fontsize=30)
plt.savefig('OutInEtavsV')

#Finding points over the line and repeating the analysis 

eta_best,flux_best,V_best, idx_best, ra_best, dec_best = func.BothOverLine(m_eta,q_eta,m_V,q_V,eta_log,V_log,flux_log,idx, ra, dec)

fig,ax1,ax2=myplt.EtaVscatter(eta_best,V_best,flux_best,freq,'EtavsVUnBest')
ax1.set_title(r'Linear Fit Stable and unstable sources',fontsize=20)
ax1.plot(flux_log[0],y_eta,label='Best fit',color='gray', ls='--')
ax2.plot(flux_log[0],y_V,label='Best fit',color='gray',ls='--')
ax1.legend(fontsize=15,markerscale=1.5)
ax2.legend(fontsize=15,markerscale=1.5)
plt.savefig('EtavsVUnBest')

inliers_best = []
outliers_best = []

for x in data_graph: 
    if func.IntervalCalc(x,mu,cov_matrix) <= chi2:
        inliers_best.append(x)
    else: 
        outliers_best.append(x)

inliers_best = np.array(inliers_best)
outliers_best = np.array(outliers_best)

print('Number of outliers : ', len(outliers_best))

if outliers_best.shape[-1] == 0:
    print('No candidate variable sources found.')
    exit
else:
    outliers_idx=func.FindOutliersIdx(outliers_best,eta_log, V_log, idx, y_eta, y_V)
    print(outliers_idx)

    outliers_eta=[]
    outliers_V=[]
    inliers_eta=[]
    inliers_V=[]
    flux_outliers=[]
    flux_inliers=[]
    
    for i,val in enumerate(idx_best[-2]):
        if val in outliers_idx:
            outliers_eta.append(idx_best[0][i])
            outliers_V.append(idx_best[1][i])
            flux_outliers.append(idx_best[-1][i])
        else:
            inliers_eta.append(idx_best[0][i])
            inliers_V.append(idx_best[1][i])
            flux_inliers.append(idx_best[-1][i])

    outliers_EtavsV = np.array((outliers_eta, outliers_V))
    inliers_EtavsV= np.array((inliers_eta,inliers_V))

    flux_outliers=np.array(flux_outliers)
    flux_inlliers=np.array(flux_inliers)

    outliers_EtavsV = outliers_EtavsV.T
    inliers_EtavsV =inliers_EtavsV.T

    figure,ax = myplt.OutInPlot(outliers_best,inliers_best,'OutIn_Unstable_Best',dists_eta,dists_V)

    EtaVsVout, axveta = myplt.OutInPlot(outliers_EtavsV,inliers_EtavsV,'OutInEtavsVbest')
    axveta.set_xlabel(r'$log_{10}(\eta_{\nu})$',fontsize=30)
    axveta.set_ylabel(r'$log_{10}(V_{\nu})$',fontsize=30)
    plt.savefig('OutInEtavsVbest')

    fig,(ax1,ax2) = plt.subplots(2,1,figsize=(14,14))
    ax1.scatter(flux_outliers,outliers_EtavsV[:,0],color='red',label='Outliers')
    ax1.scatter(flux_inliers,inliers_EtavsV[:,0],color='blue',label='Inliers')
    ax2.scatter(flux_outliers,outliers_EtavsV[:,1],color='red',label='Outliers')
    ax2.scatter(flux_inliers,inliers_EtavsV[:,1],color='blue',label='Inliers')
    ax1.set_ylabel(r'$log_{10}(\eta_{\nu}$)',fontsize=30)
    ax2.legend(fontsize=15,markerscale=1.5)
    ax1.legend(fontsize=15,markerscale=1.5)
    ax2.set_ylabel(r'$log_{10}(V_{\nu}$)',fontsize=30)
    ax2.set_xlabel(r'$log_{10}(Flux) (Jy)$',fontsize=30)
    ax1.tick_params(labelsize=15)
    ax2.tick_params(labelsize=15)
    plt.savefig('EtavsVscatterinoutBest')



