# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 12:55:20 2016

@author: Max Karlsson
"""

from setting_handler import settingsGUI
from mass_calculator import massCalculator
from mzML_reader import mzMLdataExtractor
from TICwriter import TICwriter
from compensator import compensateData, create_compensationLibrary
from writer import pickler
from analyser import experimentAnalyser
from Printicek import daviPrint
from progress_bar import statusWindow



def main():
    def fileLoop(index=0):
        dataFilePath = dataFilePaths[index]
        #Extract data according to isolation list from data file:
        indexedData = mzMLdataExtractor(dataFilePath, settings, peptideSettings)
        #Create TIC file in save directory/output:
        TICwriter(indexedData['TIC'], dataFilePath, saveDirectory)
        #Compensate MS1 & MS2 peaks:
        if settings['Co-isolation'] and settings['MS2 compensation'] or settings['MS1 compensation']:
            indexedData = compensateData(compensationLibrary, indexedData, settings)
        #Save data:
        saveFilePath = pickler(indexedData, peptideSettings, settings, saveDirectory, dataFilePath)
        experimentFiles.append(saveFilePath)
        statWin.updating(index+1, dataFilePath[index])
        if index+1<len(dataFilePaths):
            index += 1
            statWin.after(1, fileLoop, index)
    
    experimentFiles = []
    #load Settings:
    settings, isolationListPath, dataFilePaths, saveDirectory, experimentFiles, experimentName = settingsGUI()
    #Cover up for if user clicked 'Browse experiment files' without choosing 
    #any:
    if experimentFiles == '':
        experimentFiles = []
    #Only run peak calling if no experiment files were selected:
    if len(experimentFiles) == 0:
        #Calculate peptide fragment masses:
        peptideSettings = massCalculator(isolationListPath, settings)
        #Create library of compensation constants for peptide in isolation list:
        if settings['Co-isolation'] and settings['MS2 compensation']:
            compensationLibrary = create_compensationLibrary(peptideSettings, settings)
        #Parse data files:
        daviPrint('Data indexing and peak calling:', line=True)
        
        #Create status window:
        statWin = statusWindow(len(dataFilePaths))
        statWin.after(1, fileLoop, 0)
        statWin.mainloop()
            
        
    #Create result files:
    experimentAnalyser(experimentFiles, saveDirectory, experimentName, settings)
    

    
    
    
if __name__ == '__main__':
    main()