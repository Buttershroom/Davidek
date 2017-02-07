# -*- coding: utf-8 -*-
"""
Created on Sun Dec 11 12:22:34 2016

@author: Max Karlsson
"""

from pymzml import run
from Printicek import daviPrint

def getSpectrumInfo(spectrum):
    """ Inputs a spectrum object and outputs spectrum retention time and mass
        spec level, to classify spectrum type.
    """
    MSlevel = spectrum['ms level']
    
    if MSlevel == 1:
        MSlevel = 'MS1'
    elif MSlevel == 2:
        MSlevel = 'MS2'
    else:
        daviPrint('ERROR: Unexpected MS type in mzMLdataExtractor. MS type: ' + str(MSlevel), line=True, pause=True)
        
    
    return MSlevel, spectrum['scan start time']

def isTIC(spectrum):
    """ Checks if spectrum is total ion current spectrum.
    """
    if spectrum['id']=='TIC':
        return True
    return False
    
def callPeaksFromSpectrum(spectrum, settings, peptideSettings):
    """ Calls peaks from a spectrum object
    """
    indexedData = {}
    #Get spectrum retention time and MS level:
    MSlevel, RT = getSpectrumInfo(spectrum)
    if MSlevel == 'MS2':
        selectedIonMz = spectrum['selected ion m/z']
    #Parse through peptide objects:
    for peptide in peptideSettings:
        peptideObject = peptideSettings[peptide]
        #If the run is scheduled peptides outside of the retention time span 
        #are skipped:
        if settings['Scheduled run'] and not (peptideObject.scanPeriod[0] < RT < peptideObject.scanPeriod[1]):
            continue
        
        if MSlevel == 'MS1':
            referencePeakRanges = peptideObject.precursorRanges
        if MSlevel == 'MS2':
            referencePeakRanges = peptideObject.fragmentRanges
            precursorRanges = peptideObject.precursorRanges
            #Check if selected ion is this peptide. OBS: Relies on that the 
            #selected ion value has the mz of one of the precursor isotope
            #cluster peaks:
            if not any([limits[0] < selectedIonMz < limits[1] for limits in precursorRanges.values()]):
                continue
            
        for peak in spectrum.peaks:
            for referencePeak in referencePeakRanges:
                limits = referencePeakRanges[referencePeak]
                if limits[0] < peak[0] < limits[1]:
                    #b-type ions are heavy or light depending on if the 
                    #precursor ion is heavy or light, because b-ions loose
                    #their heavy isotope marker in fragmentation:
                    if referencePeak[0]=='b':
                        lightCount = 0
                        heavyCount = 0
                        #Check how many light vs heavy precursors match:
                        for precursor in precursorRanges:
                            precursorLimits = precursorRanges[precursor]
                            if precursorLimits[0] < selectedIonMz < precursorLimits[1]:
                                if precursor[1] == 'l':
                                    lightCount+=1
                                elif precursor[0] == 'h':
                                    heavyCount+=1
                        #Change peak name to state whether light or heavy:
                        if lightCount > heavyCount:
                            referencePeak = referencePeak.split(' ')[0] + ' l' + referencePeak.split(' ')[1]
                        elif lightCount < heavyCount:
                            referencePeak = referencePeak.split(' ')[0] + ' h' + referencePeak.split(' ')[1]
                        else:
                            #If it is not possible to state whether it is a
                            #heavy or light fragment, skip it:
                            continue
                    
                    #Create sub-dictionaries if non-existant:
                    if not peptide in indexedData:
                        indexedData[peptide] = {}
                    if not MSlevel in indexedData[peptide]:
                        indexedData[peptide][MSlevel] = {}
                    if not RT in indexedData[peptide][MSlevel]:
                        indexedData[peptide][MSlevel][RT] = {}
                    
                    #If the peak is not already called, add peak to indexed
                    #data. Possible implementation could be to call the 
                    #data peak which is closest to the reference peak.
                    if not referencePeak in indexedData[peptide][MSlevel][RT]:
                        indexedData[peptide][MSlevel][RT][referencePeak] = peak[1]
                        
    return indexedData
                                
def dataMerger(indexedData, newData):
    """ Merges indexedData with new Data
    """
    for peptide in newData:
        for MSlevel in newData[peptide]:
            indexedData[peptide][MSlevel].update(newData[peptide][MSlevel])
    return indexedData
            
def mzMLdataExtractor(dataFilePath, settings, peptideSettings):
    """ Extracts data from a .mzml file into a dictionary indexedData:
        -[peptide]
            -[MS level]
                -[retention time]
                    -[peak identifyer]
                        intensities
    """
    #Load file: (cannot use with statement with this module)
    dataFile = run.Reader(dataFilePath)
    numberSpectra = dataFile.getSpectrumCount()
    #Declare variables for storage of extracted data:
    indexedData = {'TIC': []}
    for peptide in peptideSettings: 
        indexedData[peptide] = {'MS1':{}, 'MS2':{}}
    daviPrint(dataFilePath.split('/')[-1])
    #Parse through spectra:
    for i, spectrum in enumerate(dataFile):
        print("\rProgress {:2.1%}".format(i/numberSpectra), end="\r")
        if settings['Reduce noise']:
            spectrum.removeNoise()
        #If the spectrum is the total ion current, save it separately and 
        #continue:
        if isTIC(spectrum):
            indexedData['TIC'] = spectrum.peaks
            continue
        #Get data from spectrum:
        spectrumData = callPeaksFromSpectrum(spectrum, settings, peptideSettings)
        #Merge spectrumData with rest of data:
        indexedData = dataMerger(indexedData, spectrumData)
    
    return indexedData

    


