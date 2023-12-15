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

global db 

# Connect to the database and run the queries
session = dbtools.access(engine,host,port,user,password,database)
VarParams = dbtools.GetVarParams(session,dataset_id)   # Get the running catalogue and varmetric catalogues and combine

plotdata = [[VarParams[i].Runningcatalog.id, VarParams[i].Varmetric.eta_int, VarParams[i].Varmetric.v_int, VarParams[i].Varmetric.lightcurve_max, 
            VarParams[i].Varmetric.lightcurve_median, (VarParams[i].Varmetric.band.freq_central/1e6), VarParams[i].Runningcatalog.datapoints, VarParams[i].Varmetric.newsource, 
            VarParams[i].Runningcatalog.wm_ra, VarParams[i].Runningcatalog.wm_decl] for i in range(len(VarParams))]  #if VarParams[i].Runningcatalog.id not in matchSrcs]
plotdata = pd.DataFrame(data=plotdata,columns=['runcat','eta','V','maxFlx','avgFlx','freq','dpts','newSrc','ra','dec'])
plotdata = plotdata.fillna('N')

# Save the data for plotting
plotdata.to_csv(outputname+'.csv', index=False)

