# plotLightcurve.py
#
# a code to plot a simple lightcurve for a source

# Import all the dependencies and generic setup
import os
import scipy as sp
import numpy as np
import pandas as pd
from sqlalchemy import *
from sqlalchemy.orm import relationship
import tkp.db
#import logging
#logging.basicConfig(level=logging.INFO)
#query_loglevel = logging.WARNING  # Set to INFO to see queries, otherwise WARNING
import sys
sys.path.append('../')
from dblogin import * # This file contains all the variables required to connect to the database
from databaseTools import dbtools
from tools import tools
from plotting import plot_varib_params as pltvp
import matplotlib
import matplotlib.pyplot as plt
import pylab
pylab.rcParams['legend.loc'] = 'best'
from matplotlib.ticker import NullFormatter
from matplotlib.font_manager import FontProperties
from astropy import units as u
from astropy.coordinates import SkyCoord
from datetime import datetime

from tkp.db.model import Runningcatalog
from tkp.db.model import Extractedsource
from tkp.db.model import Assocxtrsource
from tkp.db.model import Image

# The input database, dataset and thresholds
if len(sys.argv) == 1:
    sourceID=915
    database = 'antoniar'
elif len(sys.argv) ==3:
    database=sys.argv[1]
    sourceID=sys.argv[2]
else:
    print('Define variables in script or on command line')

outname = '/home/rvaldata/Lightcurves_MA/new/'+ str(sourceID)+'_lightcurve.png' # name of the plotdata datafile


def parse_txt(txt):
    try: 
        return datetime.strptime(txt, "%Y-%m-%d %H:%M:%S.%f")
    except ValueError:
        return datetime.strptime(txt, "%Y-%m-%d %H:%M:%S")

if os.path.exists(str(sourceID)+'_lightcurve.csv'):
    lightcurve = pd.read_csv(str(sourceID)+'_lightcurve.csv')
    lightcurve['Time'] = lightcurve.apply(lambda row: parse_txt(row.Time), axis=1)
    
else:
    # Connect to the database and run the queries
    session = dbtools.access(engine,host,port,user,password,database)

    runcatData = session.query(Runningcatalog).filter(Runningcatalog.id == sourceID).all()
    positions = [[runcatData[i].id,runcatData[i].wm_ra,runcatData[i].avg_ra_err,runcatData[i].wm_decl,runcatData[i].avg_decl_err] for i in range(len(runcatData))]
    print(positions)

    flxVals = session.query(Assocxtrsource,Extractedsource).select_from(join(Assocxtrsource,Extractedsource)).filter(Assocxtrsource.runcat_id == sourceID).all()
    lightcurve = pd.DataFrame(data=[[flxVals[x].Extractedsource.image.id, flxVals[x].Extractedsource.f_int, flxVals[x].Extractedsource.f_int_err, flxVals[x].Extractedsource.extract_type] for x in range(len(flxVals))], columns = ['Image','Flux','FluxErr','Type'])
    
    images= session.query(Image).all()
    images = pd.DataFrame(data=[[images[x].id,images[x].taustart_ts,np.around(images[x].band.freq_central/1e9, decimals=3)] for x in range(len(images))], columns=['Image','Time','Freq'])

    TMP= pd.merge(lightcurve, images, on="Image")
    lcurve = TMP.to_numpy()
    lightcurve = pd.DataFrame(lcurve, columns=['Image','Flux','FluxErr','Type','Time','Freq'])
    lightcurve.to_csv(str(sourceID)+'_lightcurve.csv')


# setting up the plot
nullfmt   = NullFormatter()         # no labels
fontP = FontProperties()
fontP.set_size('xx-large')

plt.figure(figsize=(10,5))
for index, row in lightcurve.iterrows():
    if row.Type==0:
        plt.errorbar(row.Time,row.Flux*1e3, yerr=row.FluxErr*1e3, fmt='o', markersize=7, linestyle='-',color='b')
    if row.Type==1:
        plt.errorbar(row.Time,row.Flux*1e3, yerr=row.FluxErr*1e3, fmt='v', markersize=7, linestyle='-',color='r')
    if row.Type==2:
        plt.errorbar(row.Time,row.Flux*1e3, yerr=row.FluxErr*1e3, fmt='o', markersize=7, linestyle='-',color='k')

plt.ylabel('Flux Density (mJy)', fontsize=14)
plt.xlabel('Time', fontsize=14)
plt.savefig(outname)
