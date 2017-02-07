# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 17:26:21 2016

@author: Max Karlsson
"""

from pickle import load
from sheetTools import make_resultSheet, make_detailSheet, make_summarySheet, outputFile
import matplotlib
#Turns off figure showing to save memory. Make line comment to show figures:
matplotlib.use('Agg') 
import statsmodels.api as sm
from os import makedirs
from Printicek import daviPrint
from plotTools import create_figures

    
def unpickleExperiment(loadPath):
    """ Opens a .pkl file from a given path.
    """
    with open(loadPath, 'rb') as loadfile:
        content = load(loadfile)
        return content

def coIsolatedPeakPairing(lightChannel, heavyChannel):
    """ Creates a dictionary with retention times as keys and tuples of light
        and heavy intensities (light intensity, heavy intensity) for each 
        retention time where there are both a heavy and light peak  from paired
        channels detected.
    """
    ratioChannel = {}
    
    for RT, heavyIntensity in heavyChannel:
        #Search for retention time in light channel:
        for lightRT, lightIntensity in lightChannel:
            #Break loop if found:
            if lightRT == RT:
                break
            else:
                lightIntensity = None
        
        #If retention time found in light data, add to outdata:
        if not lightIntensity==None:
            ratioChannel[RT] = (lightIntensity, heavyIntensity)
    return ratioChannel


def inRange(lightRT, RT, retentionTimeTolerance):
    """ Checks if a light channel retention time is within limits 
    """
    if RT < lightRT < lightRT + retentionTimeTolerance:
        return True
    else:
        return False
        
def asynchronalPeakPairing(lightChannel, heavyChannel):
    """ Creates a dictionary with retention times as keys and tuples of light
        and heavy intensities (light intensity, heavy intensity) where there 
        are both a heavy and light peak in range from paired channels.
    """
    
    ratioChannel = {}
    #Return empty channel if input channels are empty:
    if len(lightChannel)<2 or len(heavyChannel)<2:
        return ratioChannel
        
    RTshift = min([heavyChannel[i+1][0]-heavyChannel[i][0] for i in range(len(heavyChannel)-1)])
    for RT, heavyIntensity in heavyChannel:
        #Search for retention time in light channel:
        for lightRT, lightIntensity in lightChannel:
            #Only take light retention times that are later than the heavy:
            if lightRT<RT:
                continue
            #If the peak cannot be paired with a point with a retention time of
            #at most two transitions further, skip the peak.
            if lightRT>RT+RTshift*8:
                lightIntensity = None
                break
            #Break loop if found:
            if inRange(lightRT, RT, 0.2): #Consider having different retention time tolerance than 0.2
                break
            else:
                lightIntensity = None
        
        #If retention time found in light data, add to outdata:
        if not lightIntensity==None:
            ratioChannel[RT] = (lightIntensity, heavyIntensity)
    return ratioChannel
    
def pair_heavyAndLight(sortedData, coIsolated, settings):
    """ Pairs heavy and light data points in order to enable linear regression
        for quantification. 
    """
    pairedData = {}
    for peptide, peptideData in zip(sortedData, sortedData.values()):
        pairedData[peptide] = {}
        for channel in peptideData:
            #Continue if light channel. Only pair heavy channels:
            if 'l' in channel:
                continue
            refChannel = channel.replace('h', 'l')
            if refChannel in peptideData:
                #Name of channel containing ratio:              
                if len(channel)>2 and settings['MS2 quantification']:
                    ratioChannel = channel.split(' ')[0]
                    charge = channel.count('+')
                    ratioChannel += ' ' + channel.strip('+')[-1] + '+'*charge
                elif settings['MS1 quantification']:
                    ratioChannel = channel.replace('h', 'Prec ')
                else:
                    #If the channel is an MSlevel that is not included in 
                    #quantigication, continue:
                    continue
                
                newRatioChannel = {}
                
                #Get new chanel with paired intensities:
                if coIsolated:
                    newRatioChannel = coIsolatedPeakPairing(peptideData[refChannel], peptideData[channel])
                    
                else:
                    newRatioChannel =  asynchronalPeakPairing(peptideData[refChannel], peptideData[channel])
                    
                #Append data if not empty:
                if not len(newRatioChannel)==0:
                    pairedData[peptide][ratioChannel] = newRatioChannel
                

        #Delete empty elements:
        if len(pairedData[peptide])==0:
            del pairedData[peptide]
        
                
    return pairedData
        
def create_fragmentChannels(indexedData):
    """ Creates a "channel" for each fragment in data by making a 
        sorted list of tuples containing retention time and intensity. 
    """
    sortedData = {}
    for peptide, peptideData in zip(indexedData, indexedData.values()):
        #Skip the TIC:
        if peptide=='TIC':
            continue
        fragmentChannel = {}
        #Transform data to channel format:
        for MSlevel in peptideData:
            for RT in peptideData[MSlevel]:
                for fragment, intensity in zip(peptideData[MSlevel][RT], peptideData[MSlevel][RT].values()):
                    if fragment not in fragmentChannel:
                        fragmentChannel[fragment] = []
                    fragmentChannel[fragment].append((RT, intensity))
                    
        #Sort data from each channel according to retention time:
        for fragment in fragmentChannel:
            fragmentChannel[fragment].sort(key=lambda x: x[0])
        
        sortedData[peptide] = fragmentChannel
    return sortedData

def get_saveNames(experimentFiles):
    """ Gets savenames from the experiment file names to name result workbook
        sheets.
    """
    saveNames = {}
    for filePath in experimentFiles:
        #Removes path and file ending:
        saveName = filePath.split('/')[-1].split('.')[0]
        saveName = saveName[17:]
        saveNames[filePath] = saveName
    
    #Trims away starting substring common to all files:
    refSaveName = str(saveName)
    for i in reversed(range(1,len(refSaveName)+1)):
        try:
            if all([refSaveName[:i]==saveName[:i] for saveName in saveNames.values()]):
                for filePath, saveName in zip(saveNames, saveNames.values()):
                    saveName = saveName[i:]
                    saveNames[filePath] = saveName
                break
        except:
            pass
    
    #Trims the name down to 23 characters, the maximum sheet name + '_details' 
    #(total 31) in Excel:
    for filePath, saveName in zip(saveNames, saveNames.values()):
        if len(saveName)>23:
            saveName = saveName[-23:]
            saveNames[filePath] = saveName
    return saveNames

def calculate_linearRegressionParameters(pairedData):
    """ Calculates linear regression parameters for all channels in paired Data
    """
    linearRegressionParameters = {}
            
    for peptide in pairedData:
        linearRegressionParameters[peptide] = {}
        
        #Add data to fragment by fragment linear regression:
        for fragment, channel in zip(pairedData[peptide], pairedData[peptide].values()):
            lightInts, heavyInts = zip(*list(channel.values()))
            #At least 3 data points are needed for quantification:
            if len(lightInts)>2:
                #Make linear regression:
                fittedValues, m, k, r2 = linearRegression(heavyInts, lightInts)
                #Save parameters:
                linearRegressionParameters[peptide][fragment] = {'k':k, 'm':m, 'r2':r2, 'fitted values': fittedValues}
                    
            #Add data to total linear regression:
            allLightInts, allHeavyInts = [],[]
            for fragment, channel in zip(pairedData[peptide], pairedData[peptide].values()):
                lightInts, heavyInts = zip(*list(channel.values()))
                allLightInts += lightInts
                allHeavyInts += heavyInts
            
            #At least 3 data points are needed for quantification:
            if len(allLightInts)>2:
                #Make linear regression:
                fittedValues, m, k, r2 = linearRegression(allHeavyInts, allLightInts)
                #Save parameters:
                linearRegressionParameters[peptide]['Total'] = {'k':k, 'm':m, 'r2':r2, 'fitted values': fittedValues}
            
        #Remove empty elements:
        if len(linearRegressionParameters[peptide])==0:
            del linearRegressionParameters[peptide]
    
    return linearRegressionParameters
                

def linearRegression(heavyInts, lightInts):
    """ Makes a robust linear regression of paired light and heavy intensities
        and returns values fitted to the regression curve that can be plotted 
        to the heavy intensities to plot the regression curve. It also outputs
        the m, k and r2 values of the regression curve.
    """
    Y = lightInts
    X = sm.add_constant(heavyInts)
    #Create the regression model:
    robustLinearModel = sm.RLM(Y, X)
    #Calculate r2 value:
    r2 = sm.WLS(robustLinearModel.endog, robustLinearModel.exog, 
                weights=robustLinearModel.fit().weights).fit().rsquared
    #Get fitted values and parameters:
    results = robustLinearModel.fit()
    fittedValues = results.fittedvalues
    m, k = results.params[:2]
    
    return fittedValues, m, k, r2

    
def experimentAnalyser(experimentFiles, saveDirectory, experimentName, settings):
    """ Goes through experimetal data with heavy and light called peaks, pairs
        data points to enable linear regression for quantification of ratios
        between heavy and light peptides. Saves an xlsx file with results and
        figures of the chromatograms and linear regressions.
    """
    #Get names of output file sheets from experiment file names:
    sheetNames = get_saveNames(experimentFiles)
    #Summary variable:
    allLinRegParams = {}
    #Create image folder:
    makedirs((saveDirectory + '/output/images'), exist_ok=True)
    #Create output file:
    with outputFile(experimentName + '_results', saveDirectory) as workbook:
        daviPrint('Saving result file...',line=True)
        #Parse experiment files:
        for filePath in experimentFiles:
            daviPrint('Saving sheet for ' + filePath.split('/')[-1], line=True)
            sortedData = {}
            #Get data from pickled experiment file:
            indexedData, peptideSettings, settingsDummy = unpickleExperiment(filePath) #settingsDummy is an unused variable but necessary for unpacking.
            coIsolated = settings['Co-isolation (analysis)']
            #Sort data and put in new variable:
            sortedData = create_fragmentChannels(indexedData)
            #Delete indexed data to save memory:
            del indexedData
            #Pair heavy and light data:
            pairedData = pair_heavyAndLight(sortedData, coIsolated, settings)
            #Make quantification:
            linearRegressionParameters = calculate_linearRegressionParameters(pairedData)
            
            if settings['Save figures']:
                #Create figures:
                daviPrint('Saving figures:')
                create_figures(pairedData, sortedData, linearRegressionParameters, sheetNames[filePath], saveDirectory + '/output/images', settings)
            
            #Add data to summary variable:
            allLinRegParams[sheetNames[filePath]] = linearRegressionParameters
            #Make sheet with results in outfile:
            make_resultSheet(workbook, sheetNames[filePath], linearRegressionParameters)
            #Make sheet with detailed results in outfile:
            make_detailSheet(workbook, sheetNames[filePath], linearRegressionParameters) 
        #Make sheet with experiment summary:
        make_summarySheet(workbook, allLinRegParams, sheetNames) 


if __name__ == '__main__':
    import doctest
    experimentFiles = ('C:/Users/bjorn.forsstrom/Documents/Max/Kod/AppBio project/testData/output/2016-12-29_18-33-wellness_prm(pilot)_vis(2)_pla(1)_160522_A2_117055.pkl', 'C:/Users/bjorn.forsstrom/Documents/Max/Kod/AppBio project/testData/output/2016-12-19_16-28-Max_semi_heavy_160826_1p6_vial_6.pkl')
    saveDirectory = 'C:/Users/bjorn.forsstrom/Documents/Max/Kod/AppBio project/testData'
    experimentName = 'TESTNAME_X'
    settings = {'Save figures': 0, 'Peaks in isotope clusters': 3, 'Fragment peptide calling': 0, 'Reduce noise': 1, 'Co-isolation (analysis)': 0, 'Max fragment charge': 2, 'MS1 quantification': 0, 'MS2 quantification': 1, 'Co-isolation': 1, 'IAA': 1, 'Include b-fragments': 0, 'Isolation window offset': 0.55, 'MS1 tolerance': 10.0, 'Heavy labels': 'K:0C2N, R:0C2N', 'Isolation window size': 1.2, 'Scheduled run': 1, 'MS2 tolerance': 7.0, 'MS2 compensation': 1}
    doctest.testmod()
