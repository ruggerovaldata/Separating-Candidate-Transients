import matplotlib.pyplot as plt
import numpy as np
import pylab
pylab.rcParams['legend.loc'] = 'best'


#def EtaVscatter(Eta,V,Flux,freq,name):
def EtaVscatter(data,freq,name):

    fig,(ax1,ax2) = plt.subplots(2,1,figsize=(14,14))

    eta_graph=[]
    flux_graph=[]
    V_graph=[]
    
    for i,val in enumerate(freq):
        tmp = data.loc[data['freq'] == val]
        eta_graph.append(tmp.logEta)
        flux_graph.append(tmp.logFlux)
        V_graph.append(tmp.logV)

    eta_graph = np.array(eta_graph,dtype=list)
    flux_graph = np.array(flux_graph,dtype=list)
    V_graph = np.array(V_graph,dtype=list)


    fig,(ax1,ax2) = plt.subplots(2,1,figsize=(14,14))

    for i,val in enumerate(freq):
        val = round(val,3)
        ax1.scatter(flux_graph[i],eta_graph[i],s=20, label=str(val))
        ax2.scatter(flux_graph[i],V_graph[i],s=20, label=str(val))
    ax1.legend(fontsize=15,markerscale=1.5)
    ax1.set_ylabel(r'$log_{10}(\eta_{\nu}$)',fontsize=30)
    ax2.legend(fontsize=15,markerscale=1.5)
    ax2.set_ylabel(r'$log_{10}(V_{\nu}$)',fontsize=30)
    ax2.set_xlabel(r'$log_{10}(Flux) (Jy)$',fontsize=30)
    ax1.tick_params(labelsize=20)
    ax2.tick_params(labelsize=20)
    plt.savefig( name +'.png')
    return fig, ax1,ax2
        
def HistPlot(dist1,dist2,name):
    fig,(ax1,ax2)= plt.subplots(1,2,figsize=(14,14))
    ax1.hist(dist1,bins=20, density= True,alpha=0.5)
    ax2.hist(dist2,bins=20, density= True,alpha=0.5)
    ax1.set_title('Distance eta distribution')
    ax2.set_title('Distance V distribution')
    plt.savefig(name+'.png')
    return fig, ax1,ax2

def CornerPlot(dist1,dist2,name):
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left+width+0.02
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]
    fig = plt.figure(1,figsize=(14,14))
    axS = fig.add_subplot(223,position = rect_scatter)
    axS.set_xlim(-6,2)
    axS.set_ylim(-6,2)
    plt.xlabel(r'$\eta_{\nu}$')
    plt.ylabel(r'$V_{\nu}$')
    axHistx=fig.add_subplot(221, position=rect_histx)
    axHisty=fig.add_subplot(224, position=rect_histy)
    sns.kdeplot(dist1,dist2,fill=True, color='k',ax = axS,alpha=1)
    _,bins1,_= axHistx.hist(dist1,bins=20,density=True,alpha=0.5)
    axHistx.set_xlim(-6,2)
    _,bins2,_=axHisty.hist(dist2,bins=20,density=True,alpha=0.5,orientation='horizontal')
    axHisty.set_ylim(-6,2)
    plt.savefig(name+'.png')

    return fig, axS,  axHistx,axHisty, bins1, bins2

def EtaVscatterover(Eta,V,Flux_eta, Flux_V,freq,name):
    eta_graph=[]
    flux_graph_eta=[]

    for val in freq: 
        temp_flux=[]
        temp_eta=[]
        for i,f in enumerate(Flux_eta[1,:]):
            if val == f:
                temp_eta.append(Eta[0,i])
                temp_flux.append(Flux_eta[0,i])
        eta_graph.append(temp_eta)
        flux_graph_eta.append(temp_flux)

    eta_graph = np.array(eta_graph,dtype=list)
    flux_graph_eta = np.array(flux_graph_eta,dtype=list)

    V_graph = []
    flux_graph_V = []
    for val in freq: 
        temp_flux=[]
        temp_V =[]
        for i,f in enumerate(Flux_V[1,:]):
            if val == f:
                temp_flux.append(Flux_V[0,i])
                temp_V.append(V[0,i])
        flux_graph_V.append(temp_flux)
        V_graph.append(temp_V)
    flux_graph_V = np.array(flux_graph_V,dtype=list)
    V_graph = np.array(V_graph,dtype=list)
    
    fig,(ax1,ax2) = plt.subplots(2,1,figsize=(14,14))

    for i,val in enumerate(freq):
        val = round(val,3)
        ax1.scatter(flux_graph_eta[i],eta_graph[i],s=20, label=str(int(val)))
        ax2.scatter(flux_graph_V[i],V_graph[i],s=20, label=str(int(val)))
    
    ax1.legend(fontsize='small')
    ax1.set_title('Stable sources')
    ax1.set_ylabel(r'$\eta_{\nu}$')
    ax2.legend(fontsize='small')
    ax2.set_ylabel(r'$V_{\nu}$')
    ax2.set_xlabel('Max Flux')
    ax1.tick_params(axis='both', which='major', labelsize=20)
    ax2.tick_params(axis='both', which='major', labelsize=20)

    plt.savefig( name +'.png')
    return fig, ax1,ax2

def MyCorner(dists_eta,dists_V,likelihood,name):
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    spacing = 0.005


    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom + height + spacing, width, 0.2]
    rect_histy = [left + width + spacing, bottom, 0.2, height]  

    fig = plt.figure(figsize=(14,14))
    ax = fig.add_axes(rect_scatter)
    ax_histx = fig.add_axes(rect_histx, sharex=ax)
    ax_histy = fig.add_axes(rect_histy, sharey=ax)

    fig.suptitle('Distances Corner Plot',fontsize=30)

    ax_histx.hist(dists_eta,bins=20,density=True)
    ax_histx.tick_params(labelsize=15)
    ax_histx.set_ylabel(r'Counts',fontsize=30)
    ax_histy.hist(dists_V,bins=20,orientation='horizontal',color='orange',density=True)
    ax.scatter(dists_eta,dists_V,s=50,color='grey')
    ax.tick_params(labelsize=20)
    ax.set_ylabel(r'Dists $V_{\nu}$',fontsize=30)
    ax.set_xlabel(r'Dists $\eta_{\nu}$',fontsize=30)
    ax_histy.set_xlabel(r'Dists $V_{\nu}$',fontsize=30)
    ax_histy.tick_params(labelsize=20)
    ax_histx.tick_params(axis='both', which='major', labelsize=20)
    

    plt.savefig(name)

    ax.scatter(dists_eta,dists_V,c=likelihood,s=50)
    ax.tick_params(axis='both', which='major', labelsize=20)
    plt.savefig(name+'_Likelihood.png')
    return fig,[ax,ax_histx,ax_histy]

def OutInPlot(outliers,inliers,name,dists_eta_un=[],dists_V_un=[]):
    figure = plt.figure(figsize=(14,14))
    ax = plt.axes()
    ax.set_title(r'Predicted Outliers', fontsize=30)
    ax.scatter(inliers[:,0],inliers[:,1],color='blue',label='Inliers')
    ax.scatter(outliers[:,0],outliers[:,1],color='red',label='Outliers')
    ax.set_xlabel(r'Dists $\eta_{\nu}$',fontsize=30)
    ax.set_ylabel(r'Dists $V_{\nu}$',fontsize=30)
    #ax.scatter(dists_eta_un,dists_V_un,facecolors='None',label='TO',edgecolors='black')
    ax.legend(fontsize=20,markerscale=1.5)
    ax.tick_params(axis='both', which='major', labelsize=20)
    plt.savefig(name)
    return figure, ax
