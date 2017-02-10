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
from progress_bar import bar
import time
from statistics import median


def main():
    #--------------------------------------------------------------------------
    def fileLoop(index, fileTime):
        """
        """
        showTime = ''
        if not fileTime == 0:
            times.append(time.time()-fileTime)
            secTime = median(times)*(len(dataFilePaths)-index)
            showTime = ' '.join([t for t in time.strftime('%Hh %Mm %Ss', time.gmtime(secTime)).split(' ') if '00' not in t])
    
        fileTime = time.time()
        dataFilePath = dataFilePaths[index]
        window.currentText('Current file: ' + dataFilePath.split('/')[-1] + '\n(' + str(index+1) + '/' + str(len(dataFilePaths)) + ')\n' + showTime)
        
        #File Analysis
        #----------------------------------------------------------------------
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
        #----------------------------------------------------------------------
        window.add()
        window.update_idletasks()
        index+=1
        if not index == len(dataFilePaths):
            window.after(1, fileLoop, index, fileTime)
    #--------------------------------------------------------------------------
            
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
        
        #Create status window
        #----------------------------------------------------------------------
        fileTime = 0
        times = []
        index = 0
        #Create a progress bar window:
        window = bar(len(dataFilePaths))
        window.wm_title("Analysis status")
        window.after(10, fileLoop, index, fileTime)
        window.mainloop()
        #----------------------------------------------------------------------
        
    #Create result files:
    experimentAnalyser(experimentFiles, saveDirectory, experimentName, settings)
    

    
    
    
if __name__ == '__main__':
    main()