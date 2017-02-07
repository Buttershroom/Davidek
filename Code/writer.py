# -*- coding: utf-8 -*-
"""
Created on Sun Dec 18 19:09:40 2016

@author: Max Karlsson
"""


from pickle import dump, load, HIGHEST_PROTOCOL
from os import path, makedirs
import time

def pickleExperiment(variable, saveName, saveDirectory):
    """ Creates a file with a given save name in a given save directory with
        given contents. Can be opened by using the unpickleExperiment function.
    """
    with open(path.join(saveDirectory, saveName+'.pkl'), 'wb') as saveFile:
        dump(variable, saveFile, HIGHEST_PROTOCOL)

def unpickleExperiment(loadPath):
    """ Opens a .pkl file from a given path.
    """
    with open(loadPath, 'rb') as loadfile:
        content = load(loadfile)
        return content
        
def saveSettingsSummary(peptideSettings, settings, saveName, saveDirectory):
    with open(path.join(saveDirectory, saveName+'_summary.txt'), 'w') as saveFile:
        saveFile.write('Peptides included:\n')
        for peptide in peptideSettings:
            saveFile.write(peptide+'\n')
        saveFile.write('\nSettings:\n')
        for setting, value in zip(settings, settings.values()):
            saveFile.write(setting + '\t' + str(value) + '\n')
            
    
def pickler(indexedData, peptideSettings, settings, saveDirectory, dataFilePath):
    """ Pickles experiment data from one file.
    """
    fileName = dataFilePath.split('/')[-1].split('.')[0]
    #If no save directory was chosen, save in file location:
    if saveDirectory == '':
        saveDirectory = '/'.join(dataFilePath.split('/')[:-1])
    #Create time specific file name:
    saveName = '-'.join((time.strftime('%Y-%m-%d_%H-%M'),fileName))
    #Create output folder if non-existent:
    makedirs(path.dirname(saveDirectory + '/output/'), exist_ok=True)
    #Elongate save directory to include output folder:
    saveDirectory = saveDirectory + '/output/'
    #Save data:
    pickleExperiment((indexedData, peptideSettings, settings), saveName, saveDirectory)
    saveSettingsSummary(peptideSettings, settings, saveName, saveDirectory)
    return path.join(saveDirectory, saveName+'.pkl')

