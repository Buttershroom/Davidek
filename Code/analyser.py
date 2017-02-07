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
from analyser_tools import coIsolatedPeakPairing, asynchronalPeakPairing

    
def unpickleExperiment(loadPath):
    """ Opens a .pkl file from a given path.
    """
    with open(loadPath, 'rb') as loadfile:
        content = load(loadfile)
        return content

def pair_heavyAndLight(sortedData, coIsolated, settings):
    """ Pairs heavy and light data points in order to enable linear regression
        for quantification. 
    
    #Coisolated peak pairing:
    >>> settings = {'MS1 quantification': 0, 'MS2 quantification': 1}
    >>> sortedData = {'PEPTIDE': {'y5 l0++': [(1, 2039.6504953772994), (2, 4032.683751226058), (3, 6255.661783188084), (4, 8347.654545887279), (5, 10054.777187936974), (6, 12753.377561975742), (7, 15092.835379643253), (8, 16652.87979203035), (9, 18867.037445433903), (10, 20916.81055594802), (11, 23504.82612261323), (12, 24972.72925124599), (13, 26472.038514797292), (14, 28390.543370616673), (15, 30618.0717340172), (16, 33551.51213570254), (17, 34937.784961234545), (18, 39487.16106527799), (19, 41246.70408864039)], 'y6 h1+': [(1, 3070.1689165838934), (2, 6094.287004261431), (3, 9211.112248596273), (4, 12458.82061284305), (5, 15254.654903792823), (6, 18746.3742951612), (7, 21748.681539527115), (8, 24593.22009446358), (9, 28326.061598204593), (10, 30318.59822593895), (11, 35139.9220556041), (12, 38064.56620613969), (13, 39924.235532685954), (14, 42229.62225453523), (15, 45409.58193912281), (16, 48342.084556080285), (17, 54182.45680023217), (18, 54979.0150075922), (19, 58984.11976070764)], 'y5 h0++': [(1, 2103.6257080054893), (2, 4271.006552628103), (3, 6057.294057732438), (4, 8200.195969676499), (5, 10538.835798123931), (6, 12568.909234004188), (7, 14182.205017528417), (8, 16846.730166385205), (9, 18308.568712833963), (10, 21729.835659048294), (11, 22591.90586934931), (12, 25833.00555839655), (13, 26142.22285991264), (14, 30761.863755720227), (15, 31388.09973191338), (16, 32315.478061519214), (17, 34616.669175314004), (18, 39261.41576636129), (19, 38364.20519613705)], 'y6 l1+': [(1, 3079.9754872403346), (2, 6139.815394766756), (3, 9391.44233848008), (4, 12185.89393762708), (5, 15560.601327188167), (6, 18779.84036533217), (7, 21618.405529013562), (8, 25361.201518735656), (9, 28256.250147745304), (10, 30575.983639655773), (11, 34704.867078652904), (12, 37911.769596891056), (13, 39989.67256772389), (14, 42199.69793882989), (15, 47152.1860310885), (16, 48969.69786360805), (17, 52537.4470058678), (18, 55729.34257574856), (19, 57540.53285786968)]}}
    >>> pair_heavyAndLight(sortedData, 1, settings)=={'PEPTIDE': {'y5 0++': [(2039.6504953772994, 2103.6257080054893), (4032.683751226058, 4271.006552628103), (6255.661783188084, 6057.294057732438), (8347.654545887279, 8200.195969676499), (10054.777187936974, 10538.835798123931), (12753.377561975742, 12568.909234004188), (15092.835379643253, 14182.205017528417), (16652.87979203035, 16846.730166385205), (18867.037445433903, 18308.568712833963), (20916.81055594802, 21729.835659048294), (23504.82612261323, 22591.90586934931), (24972.72925124599, 25833.00555839655), (26472.038514797292, 26142.22285991264), (28390.543370616673, 30761.863755720227), (30618.0717340172, 31388.09973191338), (33551.51213570254, 32315.478061519214), (34937.784961234545, 34616.669175314004), (39487.16106527799, 39261.41576636129), (41246.70408864039, 38364.20519613705)], 'y6 1+':[(3079.9754872403346, 3070.1689165838934), (6139.815394766756, 6094.287004261431), (9391.44233848008, 9211.112248596273), (12185.89393762708, 12458.82061284305), (15560.601327188167, 15254.654903792823), (18779.84036533217, 18746.3742951612), (21618.405529013562, 21748.681539527115), (25361.201518735656, 24593.22009446358), (28256.250147745304, 28326.061598204593), (30575.983639655773, 30318.59822593895), (34704.867078652904, 35139.9220556041), (37911.769596891056, 38064.56620613969), (39989.67256772389, 39924.235532685954), (42199.69793882989, 42229.62225453523), (47152.1860310885, 45409.58193912281), (48969.69786360805, 48342.084556080285), (52537.4470058678, 54182.45680023217), (55729.34257574856, 54979.0150075922), (57540.53285786968, 58984.11976070764)]}}
    True
    
    #Not coisolated peak pairing:
    >>> settings = {'MS1 quantification': 0, 'MS2 quantification': 1}
    >>> sortedData = {'PEPTIDE': {'y5 l0++': [(1, 2039.6504953772994), (2, 4032.683751226058), (3, 6255.661783188084), (4, 8347.654545887279), (5, 10054.777187936974), (6, 12753.377561975742), (7, 15092.835379643253), (8, 16652.87979203035), (9, 18867.037445433903), (10, 20916.81055594802), (11, 23504.82612261323), (12, 24972.72925124599), (13, 26472.038514797292), (14, 28390.543370616673), (15, 30618.0717340172), (16, 33551.51213570254), (17, 34937.784961234545), (18, 39487.16106527799), (19, 41246.70408864039)], 'y6 h1+': [(1.02, 3070.1689165838934), (2.02, 6094.287004261431), (3.02, 9211.112248596273), (4.02, 12458.82061284305), (5.02, 15254.654903792823), (6.02, 18746.3742951612), (7.02, 21748.681539527115), (8.020000000000001, 24593.22009446358), (9.020000000000001, 28326.061598204593), (10.020000000000001, 30318.59822593895), (11.020000000000001, 35139.9220556041), (12.020000000000001, 38064.56620613969), (13.020000000000001, 39924.235532685954), (14.020000000000001, 42229.62225453523), (15.020000000000001, 45409.58193912281), (16.02, 48342.084556080285), (17.02, 54182.45680023217), (18.02, 54979.0150075922), (19.02, 58984.11976070764)], 'y5 h0++': [(1.02, 2103.6257080054893), (2.02, 4271.006552628103), (3.02, 6057.294057732438), (4.02, 8200.195969676499), (5.02, 10538.835798123931), (6.02, 12568.909234004188), (7.02, 14182.205017528417), (8.020000000000001, 16846.730166385205), (9.020000000000001, 18308.568712833963), (10.020000000000001, 21729.835659048294), (11.020000000000001, 22591.90586934931), (12.020000000000001, 25833.00555839655), (13.020000000000001, 26142.22285991264), (14.020000000000001, 30761.863755720227), (15.020000000000001, 31388.09973191338), (16.02, 32315.478061519214), (17.02, 34616.669175314004), (18.02, 39261.41576636129), (19.02, 38364.20519613705)], 'y6 l1+': [(1, 3079.9754872403346), (2, 6139.815394766756), (3, 9391.44233848008), (4, 12185.89393762708), (5, 15560.601327188167), (6, 18779.84036533217), (7, 21618.405529013562), (8, 25361.201518735656), (9, 28256.250147745304), (10, 30575.983639655773), (11, 34704.867078652904), (12, 37911.769596891056), (13, 39989.67256772389), (14, 42199.69793882989), (15, 47152.1860310885), (16, 48969.69786360805), (17, 52537.4470058678), (18, 55729.34257574856), (19, 57540.53285786968)]}}
    >>> pair_heavyAndLight(sortedData, 0, settings)=={'PEPTIDE': {'y5 0++': [(2039.6504953772994, 2103.6257080054893), (4032.683751226058, 4271.006552628103), (6255.661783188084, 6057.294057732438), (8347.654545887279, 8200.195969676499), (10054.777187936974, 10538.835798123931), (12753.377561975742, 12568.909234004188), (15092.835379643253, 14182.205017528417), (16652.87979203035, 16846.730166385205), (18867.037445433903, 18308.568712833963), (20916.81055594802, 21729.835659048294), (23504.82612261323, 22591.90586934931), (24972.72925124599, 25833.00555839655), (26472.038514797292, 26142.22285991264), (28390.543370616673, 30761.863755720227), (30618.0717340172, 31388.09973191338), (33551.51213570254, 32315.478061519214), (34937.784961234545, 34616.669175314004), (39487.16106527799, 39261.41576636129), (41246.70408864039, 38364.20519613705)], 'y6 1+':[(3079.9754872403346, 3070.1689165838934), (6139.815394766756, 6094.287004261431), (9391.44233848008, 9211.112248596273), (12185.89393762708, 12458.82061284305), (15560.601327188167, 15254.654903792823), (18779.84036533217, 18746.3742951612), (21618.405529013562, 21748.681539527115), (25361.201518735656, 24593.22009446358), (28256.250147745304, 28326.061598204593), (30575.983639655773, 30318.59822593895), (34704.867078652904, 35139.9220556041), (37911.769596891056, 38064.56620613969), (39989.67256772389, 39924.235532685954), (42199.69793882989, 42229.62225453523), (47152.1860310885, 45409.58193912281), (48969.69786360805, 48342.084556080285), (52537.4470058678, 54182.45680023217), (55729.34257574856, 54979.0150075922), (57540.53285786968, 58984.11976070764)]}}
    True
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
        
    >>> indexedData = {'peptide': {'MS1':{},'MS2':{1:{"y6 h1+":3070.1689165838934},2:{"y6 h1+":6094.287004261431},3:{"y6 h1+":9211.112248596273},4:{"y6 h1+":12458.82061284305},5:{"y6 h1+":15254.654903792823},6:{"y6 h1+":18746.3742951612},7:{"y6 h1+":21748.681539527115},8:{"y6 h1+":24593.22009446358},9:{"y6 h1+":28326.061598204593},10:{"y6 h1+":30318.59822593895},11:{"y6 h1+":35139.9220556041},12:{"y6 h1+":38064.56620613969},13:{"y6 h1+":39924.235532685954},14:{"y6 h1+":42229.62225453523},15:{"y6 h1+":45409.58193912281},16:{"y6 h1+":48342.084556080285},17:{"y6 h1+":54182.45680023217},18:{"y6 h1+":54979.0150075922},19:{"y6 h1+":58984.11976070764}}}}
    >>> create_fragmentChannels(indexedData)
    {'peptide': {'y6 h1+': [(1, 3070.1689165838934), (2, 6094.287004261431), (3, 9211.112248596273), (4, 12458.82061284305), (5, 15254.654903792823), (6, 18746.3742951612), (7, 21748.681539527115), (8, 24593.22009446358), (9, 28326.061598204593), (10, 30318.59822593895), (11, 35139.9220556041), (12, 38064.56620613969), (13, 39924.235532685954), (14, 42229.62225453523), (15, 45409.58193912281), (16, 48342.084556080285), (17, 54182.45680023217), (18, 54979.0150075922), (19, 58984.11976070764)]}}
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
            lightInts, heavyInts = zip(*channel)
            #At least 3 data points are needed for quantification:
            if len(lightInts)>2:
                #Make linear regression:
                fittedValues, m, k, r2 = linearRegression(heavyInts, lightInts)
                #Save parameters:
                linearRegressionParameters[peptide][fragment] = {'k':k, 'm':m, 'r2':r2, 'fitted values': fittedValues}
                    
            #Add data to total linear regression:
            allLightInts, allHeavyInts = [],[]
            for fragment, channel in zip(pairedData[peptide], pairedData[peptide].values()):
                lightInts, heavyInts = zip(*channel)
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
    
    experimentFiles = ('C:/Users/max.karlsson/Documents/Max/Kod/AppBio project/testData/output/2016-12-29_18-33-wellness_prm(pilot)_vis(2)_pla(1)_160522_A2_117055.pkl', 'C:/Users/bjorn.forsstrom/Documents/Max/Kod/AppBio project/testData/output/2016-12-19_16-28-Max_semi_heavy_160826_1p6_vial_6.pkl')
    saveDirectory = 'C:/Users/max.karlsson/Documents/Max/Kod/AppBio project/testData'
    experimentName = 'TESTNAME_X'
    settings = {'Save figures': 0, 'Peaks in isotope clusters': 3, 'Fragment peptide calling': 0, 'Reduce noise': 1, 'Co-isolation (analysis)': 0, 'Max fragment charge': 2, 'MS1 quantification': 0, 'MS2 quantification': 1, 'Co-isolation': 1, 'IAA': 1, 'Include b-fragments': 0, 'Isolation window offset': 0.55, 'MS1 tolerance': 10.0, 'Heavy labels': 'K:0C2N, R:0C2N', 'Isolation window size': 1.2, 'Scheduled run': 1, 'MS2 tolerance': 7.0, 'MS2 compensation': 1}
    #experimentAnalyser(experimentFiles, saveDirectory, experimentName, settings)
    import doctest
    doctest.testmod()
