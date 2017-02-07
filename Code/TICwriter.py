# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 13:03:56 2016

@author: Max Karlsson
"""
import matplotlib
#Forces matplotlib to not show plot:
matplotlib.use('Agg')
import matplotlib.pyplot as pl    
import os

def TICwriter(TIC, dataFile, saveDirectory):
    """ Saves a TIC file as a .png in the specified saveDirectory under a 
        output subdirectory.
    """
    #Create savename from data file name:
    savefile =  dataFile.split('/')[-1].split('.')[0] + '_TIC.png'
    #Create ouput directory:
    saveDirectory = os.path.join(saveDirectory, 'output/')
    os.makedirs(os.path.dirname(saveDirectory), exist_ok=True)
    #Plot figure:
    Plot = pl.figure()
    TICplot = Plot.add_subplot(111)
    TICplot.plot([d[0] for d in TIC], [d[1] for d in TIC])
    
    #Save and close plot:
    pl.savefig(saveDirectory + savefile)
    pl.close(Plot)
    
