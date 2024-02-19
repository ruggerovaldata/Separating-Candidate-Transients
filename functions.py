from scipy.stats import norm
from scipy.stats import chi2
import numpy as np
import pandas as pd 
from tabulate import tabulate

def GaussianFit(data, mu,sigma):
    range_data=np.linspace(min(data),max(data),len(data))
    fit=norm.pdf(range_data,loc=mu,scale=sigma)
    return fit 

def LinearFit(data,a,b):
    return a*data+b

def EuclidDistance(x1,y1,x2,y2):
    return np.sqrt((x1-x2)**2+(y1-y2)**2)

def LineDistance(x1,y1,m,q):
    n = np.abs(y1-(m*x1+q))
    d = np.sqrt(1+m**2)
    return n/d

def Params_distance(y_o,y_p):
    return (y_o-y_p)

def IntervalCalc(x,mu,sigma):
    diff = x-mu
    first = np.dot(diff.T,np.linalg.inv(sigma))
    return first.dot(diff)

def Data_Load(name):
    data = pd.read_csv(name)
    data = data.loc[(data['V']>0.) & (data['eta']>0.)]
    eta = np.array((data['eta'],data['freq']))
    V = np.array((data['V'],data['freq']))
    max_flux = np.array((data['maxFlx'],data['freq']))
    idx = np.array((np.log10(data['eta']),np.log10(data['V']),data['runcat'],np.log10(data['maxFlx'])))
    return eta,V,max_flux, idx

def OverLine(m,q,eta,flux):
    eta_best=[]
    flux_best_eta=[]
    for i,val in enumerate(eta[0]):
        if val - (m*flux[0,i]+q)>=0:
            flux_best_eta.append(flux[:,i])
            eta_best.append(eta[:,i])

    eta_best = np.array(eta_best)
    flux_best_eta = np.array(flux_best_eta)

    return eta_best.T, flux_best_eta.T

def BothOverLine(m_eta,q_eta,m_V,q_V,eta,V,flux,idx,ra,dec):
    eta_best=[]
    flux_best_eta=[]
    V_best=[]
    idx_best=[]
    ra_best = []
    dec_best = []
    for i,val in enumerate(eta[0]):
        if (val - (m_eta*flux[0,i]+q_eta)>=0) and (V[0,i]-(m_V*flux[0,i]+q_V)>=0):
            #print(idx[:,i],val,m_eta*flux[0,i]+q_eta,V[0,i],m_V*flux[0,i]+q_V)
            flux_best_eta.append(flux[:,i])
            eta_best.append(eta[:,i])
            V_best.append(V[:,i])
            idx_best.append(idx[:,i])
            ra_best.append(ra[i])
            dec_best.append(dec[i])
    eta_best = np.array(eta_best)
    flux_best_eta = np.array(flux_best_eta)
    V_best = np.array(V_best)
    idx_best = np.array(idx_best)
    ra_best = np.array(ra_best)
    dec_best = np.array(dec_best)

    return eta_best.T, flux_best_eta.T, V_best.T, idx_best.T, ra_best.T, dec_best.T

def Accuracy(pred_out,true_out):
    aset = set([tuple(x) for x in pred_out])
    bset = set([tuple(x) for x in true_out])
    common = np.array([x for x in aset & bset])
    acc = len(common)*100/len(pred_out)
    return common, acc

def FindOutliersIdx(outliers,eta_best,V_best,idx,y_eta_best,y_V_best):

    outliers_found_eta=[]
    outliers_found_V = []
    for val in outliers:
        for i,eta in enumerate(eta_best[0]):
            if val[0]==Params_distance(eta,y_eta_best[i]) and val[1] == Params_distance(V_best[0][i],y_V_best[i]):
                outliers_found_eta.append(eta)
                outliers_found_V.append(V_best[0][i])

    outliers_found=np.array((outliers_found_eta,outliers_found_V))

    outliers_idx=[]
    temp=idx.T
    for val in temp:
        if val[0] in outliers_found[0] and val[1] in outliers_found[1]:
            outliers_idx.append(int(val[2]))
    return outliers_idx

def Probability(data, mean, sigma):
    prob=[]
    for x in data: 
        dist = np.dot((x-mean).transpose(),np.linalg.inv(sigma))
        dist = np.dot(dist,(x-mean))
        prob.append(1- chi2.cdf(dist,2))
    prob = np.array(prob)        

    return prob*100

def GetOutput(idx,ra, dec, eta,V,dists_eta,dists_V,percentage, title):
    #output = np.vstack([idx[2],ra, dec, eta[0],V[0],dists_eta,dists_V,percentage])
    #output = output.T

    #sort = np.argsort(output[:,-1])
    #output = output[sort]

    output = [idx,ra, dec, eta,V,dists_eta,dists_V,percentage]

    file = open(title+'.txt','w')
    file.write(tabulate(output,headers=['Runcat','RA','DEC','Eta','V','Dists eta','Dists V','Percentage'],floatfmt=('.0f','.4f','.4f','.3f','.3f','.3f','.3f','.10f')))
    file.close()
