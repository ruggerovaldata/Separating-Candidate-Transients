import scipy as sp 
import numpy as np 
import pandas as pd 
import sqlalchemy
from sqlalchemy import *
from sqlalchemy.orm import relationship
import logging
logging.basicConfig(level=logging.INFO)
query_loglevel = logging.WARNING
import sys 
sys.path.append
from astropy import units as u
from astropy.coordinates import SkyCoord

from dblogin import *
from databaseTools import dbtools 
import tkp.db

from tkp.db.model import Skyregion
from tkp.db.model import Assocxtrsource

FilterRadius = 0.051 #Radius inside of which the sources will be extracted to the csv file
maxSigma = 5.56
BMaj = 0.002185 # in degrees
beamwidths=5.
SrcAssocRadius = 0
#from plotting import plot_varib_params as pltvp

outputname = 'freq3' #Name for the output .csv file
database = 'freq3rvaldata' 
websiteURL = 'http://banana.transientskp.org/r3/vlo_'+database+'/runningcatalog/'
dataset_id = 1

websiteURL = 'http://banana.transientskp.org/r3/vlo_'+database+'/runningcatalog/'
latexURL = '\url{http://banana.transientskp.org/r3/vlo_'+database+'/runningcatalog/'
latexHREF = '\href{http://banana.transientskp.org/r3/vlo_'+database+'/runningcatalog/'
query_loglevel = logging.WARNING  # Set to INFO to see queries, otherwise WARNING

global db 

# Connect to the database and run the queries
session = dbtools.access(engine,host,port,user,password,database)
NewSrcs = dbtools.GetNewSrcs(session,dataset_id)       # Get the new source table
VarParams = dbtools.GetVarParams(session,dataset_id)   # Get the running catalogue and varmetric catalogues and combine
Img1Srcs = dbtools.GetImg1Srcs(session,dataset_id)     # Get all the sources identified in the first image of the dataset


# Get co-ordinates, fluxes and other information for the new sources
NewSrcData=[[NewSrcs[i].Runningcatalog.id,NewSrcs[i].Runningcatalog.wm_ra,NewSrcs[i].Runningcatalog.wm_decl,NewSrcs[i].Newsource.trigger_xtrsrc.id,NewSrcs[i].Newsource.trigger_xtrsrc.image.taustart_ts,NewSrcs[i].Newsource.trigger_xtrsrc.image.band.id,NewSrcs[i].Newsource.trigger_xtrsrc.det_sigma] for i in range(len(NewSrcs))]
NewSrcDataFrame = pd.DataFrame(data=NewSrcData, columns=['RuncatID','ra','decl','xtrsrc','TimeDetect','Band','detSigma'])
NewSrcDataFrame = NewSrcDataFrame.sort_values(by=['RuncatID'])

print "Number of new sources: "+str(len(NewSrcDataFrame))


# Filter out all new sources that were below the detection threshold for transient sources
NewSrcDataFrame = NewSrcDataFrame.loc[NewSrcDataFrame['detSigma'] >= maxSigma]

# Make a list of all the positions of the new sources 
NewSrcs_coord = SkyCoord(ra=(NewSrcDataFrame.ra*u.degree).values,dec=(NewSrcDataFrame.decl*u.degree).values)

VarData = [[VarParams[i].Runningcatalog.id, VarParams[i].Varmetric.eta_int, VarParams[i].Varmetric.v_int, VarParams[i].Varmetric.lightcurve_max, VarParams[i].Varmetric.lightcurve_median, (VarParams[i].Varmetric.band.freq_central/1e6), VarParams[i].Runningcatalog.datapoints, VarParams[i].Varmetric.newsource] for i in range(len(VarParams))]
VarData = pd.DataFrame(data=VarData,columns=['RuncatID','eta','V','maxFlx','avgFlx','freq','dpts','newSrc'])
VarData = VarData.sort_values(by=['RuncatID'])
VarData = VarData.fillna('N')
NewSrcDataFrame = pd.merge(NewSrcDataFrame, VarData, on="RuncatID")

print "Number of new sources: "+str(len(NewSrcs_coord))


# Get co-ordinates of all sources in the first image
Img1VarParamsData=[[VarParams[i].Runningcatalog.id,VarParams[i].Runningcatalog.wm_ra,VarParams[i].Runningcatalog.wm_decl] for i in range(len(VarParams)) if VarParams[i].Runningcatalog.id in Img1Srcs]
Img1VarParamsDataFrame = pd.DataFrame(data=Img1VarParamsData, columns=['RuncatID','ra','decl'])
Img1VarParamsDataFrame = Img1VarParamsDataFrame.sort_values(by=['RuncatID'])
Img1VarParams_coord = SkyCoord(ra=(Img1VarParamsDataFrame.ra*u.degree).values,dec=(Img1VarParamsDataFrame.decl*u.degree).values)

# From CatalogueMatching.ipynb
# idx, d2d, d3d = aart_coord.match_to_catalog_sky(green_coord)
#idx is an array of indices for the green Dataframe, such that green.iloc[idx[i]] is the nearest source to aart.iloc[i].  len(idx) = len(aart)
#d2d is the pointwise 2D angular distances. len(d2d) = len(aart)
#d3d is the 3D distance, usefull if distance to each source is known, otherwise they're assumed to be on a unit sphere. To add physical distances to the catalogued soueces see http://docs.astropy.org/en/stable/coordinates/matchsep.html. len(d3d) = len(aart)


###### Filter 1 - remove any poorly associated sources
###### i.e. those that were detected in the first image but TraP incorrectly labeled as new sources
idx, d2d, d3d = NewSrcs_coord.match_to_catalog_sky(Img1VarParams_coord)
UnassocNewSrcs = NewSrcDataFrame[d2d.deg > SrcAssocRadius] # Unassociated sources - those not found in Img1
UnassocNewSrcs_coord = NewSrcs_coord[d2d.deg > SrcAssocRadius] # Unassociated sources - those not found in Img1

print "Number of new sources after Filter 1: "+str(len(UnassocNewSrcs_coord))

uniqueIDlist=list(UnassocNewSrcs.RuncatID)


###### Filter 2 - Reject sources too close to the source extraction radius.
# Not general enough at the moment...

skyreg = session.query(Skyregion).select_from(Skyregion).filter(Skyregion.dataset_id == dataset_id).one()
centre = SkyCoord(ra=(skyreg.centre_ra*u.degree),dec=(skyreg.centre_decl*u.degree))
xtrRadius = skyreg.xtr_radius

uniqueIDlistTMP=[]
for a in range(len(UnassocNewSrcs)):
    sep = UnassocNewSrcs_coord[a].separation(centre)
    if sep.degree < xtrRadius - FilterRadius:
        uniqueIDlistTMP.append(UnassocNewSrcs.iloc[a].RuncatID)

uniqueIDlist=uniqueIDlistTMP

print "Number of new sources after Filter 2: "+str(len(uniqueIDlist))

VarParams = dbtools.GetVarParams(session,dataset_id)  


plotdata = [[VarParams[i].Runningcatalog.id, VarParams[i].Varmetric.eta_int, VarParams[i].Varmetric.v_int, VarParams[i].Varmetric.lightcurve_max, VarParams[i].Varmetric.lightcurve_median, (VarParams[i].Varmetric.band.freq_central/1e6), VarParams[i].Runningcatalog.datapoints, VarParams[i].Varmetric.newsource] for i in range(len(VarParams))]
plotdata = pd.DataFrame(data=plotdata,columns=['runcat','eta','V','maxFlx','avgFlx','freq','dpts','newSrc'])
plotdata = plotdata.fillna('N')

plotdata.to_csv('~/Transients/DbMine/'+outputname+'.csv',index=False)

