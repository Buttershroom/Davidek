# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 13:44:04 2016

@author: Max Karlsson
"""

def findMaxMS2Intensity(spectrum, peak):
    """ Finds the peak with the most intensity in the light cluster in the 
        spectrum.
    >>> spectrum = {'y3 l0++': 2000,'y3 l1++': 1000,'y3 l2++': 1000,'y3 h0++': 20000,'y4 h0++': 10000,'y2 l0++': 5000,}
    >>> findMaxMS2Intensity(spectrum, 'y3 h0++')
    """
    maxPeak = None
    maxIntensity = 0
    #Get the identifyer of the cluster (e.g. 'y5'):
    clusterID = peak.split(' ')[0]
    charge = peak.count('+')
    #Parse through spectrum to find maximum intensity:
    for refPeak, intensity in zip(spectrum, spectrum.values()):
        if not refPeak.startswith(clusterID):
            continue
        if 'h' in refPeak:
            continue
        if not charge == refPeak.count('+'):
            continue
        if peak == refPeak:
            continue
        if intensity>maxIntensity:
            maxPeak = refPeak
            maxIntensity = intensity
            
    if not maxPeak==None and 'h' in maxPeak:
        print('HEAVY')
    return maxPeak, maxIntensity

def calculateTheoreticalMS2Intensity(peak, refPeak, refIntensity, precursorDistribution):
    """ Calculated the theoretical intensity from the light peak with the 
        highest intensity.
    
    >>> calculateTheoreticalMS2Intensity('y8 l2++', 'y8 l0++', 1000, [0.6, 0.3, 0.1])
    166.66666666666669
    """
    print(peak)
    print(refPeak)
    print(precursorDistribution)
    #Get peak numbers to find right abundance:
    peakNumber = int(peak.strip('+')[-1])
    refPeakNumber = int(refPeak.strip('+')[-1])
    peakAbundance = precursorDistribution[peakNumber]
    refPeakAbundance = precursorDistribution[refPeakNumber]
    theoreticalIntensity = peakAbundance*(refIntensity/refPeakAbundance)
    return theoreticalIntensity
    
def compensatePeakOverlap(overlappingPeaks, lightMS1distribution, spectrum, peak):
    """
    >>> spectrum = {'y3 l0++': 2000,'y3 l1++': 1000,'y3 l2++': 1000,'y3 h0++': 20000,'y4 h0++': 10000,'y2 l0++': 5000}
    >>> overlappingPeaks = {'y3 l2++': 'y3 h0++', 'y3 h0++': 'y3 l2++'}
    >>> lightMS1distribution = [0.6, 0.3, 0.1]
    >>> compensatePeakOverlap(overlappingPeaks, lightMS1distribution, spectrum, 'y3 h0++')
    """
    overlappingPeak = overlappingPeaks[peak]
    #Check if there is a overlapping peak:
    if not overlappingPeak==None:
        heavyLightIdentifyer = peak.split(' ')[1][0]
        if heavyLightIdentifyer == 'l':
            #Rename peak to the overlapping peak. Always prioritize
            #heavy peaks as these should constitute a larger 
            #intensity:
            peak.replace('l','h')
            heavyLightIdentifyer = 'h'
            
        if heavyLightIdentifyer == 'h':
            #Find largest light peak:
            maxPeak, maxIntensity = findMaxMS2Intensity(spectrum, peak)
            #If no light peaks were found continue:
            if maxPeak == None:
                return None
            precursorDistribution = lightMS1distribution
            #Calculate theoretical intensity of overlapping peak:
            theoreticalIntensity = calculateTheoreticalMS2Intensity(
                    overlappingPeak, maxPeak, maxIntensity, 
                    precursorDistribution)
            #Subtract theoretical overlapping intensity from peak
            #intensity:
            intensity = spectrum[peak]-theoreticalIntensity
    else:
        intensity = spectrum[peak]
    
    return intensity
        
def compensateMS2data(MS2data, peptide, compObj, compensatedData, compensationConstants):
    """
    
    """
    for RT in MS2data:
            spectrum = MS2data[RT]
            for peak in spectrum:
                #Ignore if b-fragment:
                if peak[0] == 'b':
                    continue
                
                intensity = compensatePeakOverlap(compObj.overlappingPeaks, compObj.lightMS1distribution, spectrum, peak)
                if intensity == None:
                    continue
                
                #Multiply peak intensity with compensation constant:
                intensity *= compensationConstants[peak]
                #Check if peak intensity is less than zero. If true, ignore:
                if intensity<=0:
                    continue
                
                #Add compensated peak to new data variable:
                try:
                    compensatedData[peptide]['MS2'][RT][peak] = intensity
                except:
                    compensatedData[peptide]['MS2'][RT] = {}
                    compensatedData[peptide]['MS2'][RT][peak] = intensity
    
    return compensatedData

'''
def findMaxMS1Intensity(spectrum, peak):
    """ Finds the peak with the most intensity in the light cluster in the 
        spectrum.
    """
    maxPeak = None
    maxIntensity = 0
    for refPeak, intensity in zip(spectrum, spectrum.values()):
        heavyLightIdentifyer = refPeak[0]
        if heavyLightIdentifyer=='h':
            continue
        if intensity>maxIntensity:
            maxPeak = peak
            maxIntensity = intensity
    return maxPeak, maxIntensity

def calculateTheoreticalMS1Intensity(peak, refPeak, refIntensity, precursorDistribution):
    """ Calculated the theoretical intensity from the light peak with the 
        highest intensity.
    """
    #Get peak numbers to find right abundance:
    peakNumber = int(peak.strip('+')[-1])
    refPeakNumber = int(refPeak.strip('+')[-1])
    print(peakNumber, precursorDistribution, peak, refPeak)
    peakAbundance = precursorDistribution[peakNumber]
    refPeakAbundance = precursorDistribution[refPeakNumber]
    theoreticalIntensity = peakAbundance*(refIntensity/refPeakAbundance)
    return theoreticalIntensity

def compensateMS1data(MS1data, peptide, compObj, compensatedData):
    for RT in MS1data:
            spectrum = MS1data[RT]
            print(spectrum)
            for peak in spectrum:
                overlappingPeak = compObj.overlappingPeaks[peak]
                #Check if there is a overlapping peak:
                if not overlappingPeak==None:
                    heavyLightIdentifyer = peak[0]
                    if heavyLightIdentifyer == 'l':
                        #Rename peak to the overlapping peak. Always prioritize
                        #heavy peaks as these should constitute a larger 
                        #intensity:
                        peak.replace('l','h')
                        heavyLightIdentifyer = 'h'
                        
                    if heavyLightIdentifyer == 'h':
                        #Find largest light peak:
                        maxPeak, maxIntensity = findMaxMS1Intensity(spectrum, peak)
                        #If no light peaks were found continue:
                        if maxPeak == None:
                            continue
                        #Calculate theoretical intensity of overlapping peak:
                        theoreticalIntensity = calculateTheoreticalMS1Intensity(
                                overlappingPeak, maxPeak, maxIntensity, 
                                compObj.MS1distribution)
                        #Subtract theoretical overlapping intensity from peak
                        #intensity:
                        intensity = spectrum[peak]-theoreticalIntensity
                else:
                    intensity = spectrum[peak]
                
                #Check if peak intensity is less than zero. If true, ignore:
                if intensity<=0:
                    continue
                
                #Add compensated peak to new data variable:
                try:
                    compensatedData[peptide]['MS1'][RT][peak] = intensity
                except:
                    compensatedData[peptide]['MS1'][RT] = {}
                    compensatedData[peptide]['MS1'][RT][peak] = intensity
    
    return compensatedData

'''
def compensateData(compensationLibrary, indexedData, settings):
    compensatedData = {}
    #Create compensation library:    
    for peptide in indexedData:
        if peptide == 'TIC':
            continue
        #Add MS1 data to compensated data:
        MS1data = indexedData[peptide]['MS1']
        MS2data = indexedData[peptide]['MS2']
        
        compensatedData[peptide] = {'MS1':{}, 'MS2':{}}
        #Get compensation object:
        compObj = compensationLibrary[peptide]
        compensationConstants = compObj.compensationConstants
        
        if settings['MS2 compensation']:
            #Compensate MS2data
            compensatedData = compensateMS2data(MS2data, peptide, compObj, compensatedData, compensationConstants)
        
        '''
        if settings['MS1 compensation']:
            #Compensate MS1data
            compensatedData = compensateMS1data(MS1data, peptide, compObj, compensatedData)
        '''
                    
    return compensatedData
                    
            

if __name__ == '__main__':
    
    import doctest
    doctest.testmod()
