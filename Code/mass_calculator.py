# -*- coding: utf-8 -*-
"""
Created on Sat Dec 10 15:28:30 2016

@author: Max Karlsson
"""
import re
import mass_calculator_tools as tools
"""
tools.calculate_fragmentPeakCluster(sequence, charge, IAA, maxClusterPeaks, labels, ionType, includedPeaks)
tools.calculate_heavyShift(sequence, labels)
tools.calculate_precursorPeakCluster(sequence, charge, maxClusterPeaks, labels, IAA=False)
tools.calculate_mass(sequence, charge, ionType, IAA)
"""

class peptide:
    """ Peptide object containing peptide specific variables such as sequence,
        mass, scan period, fragment masses.
    """
    def __init__(self, sequence, charge, scanStart, scanEnd, isoWinSize, isoWinOffset, labels, coIsolation, maxPeaksInCluster, IAA):
        self.sequence = sequence
        self.charge = charge
        self.scanPeriod = (scanStart, scanEnd)
        self.fragments = {}
        self.includedPeaks = self.__calculate_included_peaks(isoWinSize, isoWinOffset, labels, coIsolation, maxPeaksInCluster)
        self.precursorMass = tools.calculate_precursorPeakCluster(sequence, charge, maxPeaksInCluster, labels, IAA=IAA)
    
    def __calculate_included_peaks(self, isoWinSize, isoWinOffset, labels, coIsolation, maxPeaksInCluster):
        """ Function to calculate which peaks in a cluster are isolated by the 
            isolation window size and offset.
        
        #Co-isolated:
        >>> a = peptide('PEPTIDECK', 1, 0, 10, 3, 1.3, 'K:0C2N, R:0C2N', 1, 3, 1)
        >>> a.includedPeaks
        ['l0', 'h0', 'l1', 'l2']
        
        #Not co-isolated:
        >>> a = peptide('PEPTIDECK', 1, 0, 10, 3, 1.3, 'K:0C2N, R:0C2N', 0, 3, 1)
        >>> a.includedPeaks
        ['l0', 'l1', 'l2', 'h0', 'h1', 'h2']
        """
        includedPeaks = []
        limits = (isoWinOffset-isoWinSize/2, isoWinOffset+isoWinSize/2)
        heavyShift = tools.calculate_heavyShift(self.sequence, labels)
        
        mz = 0
        #Check if peaks are within isolation window
        peakNumber = 0
        while mz < limits[1]:
            #light peak
            if limits[0] < mz < limits[1]:
                includedPeaks.append('l'+str(peakNumber))
            if coIsolation:
                #heavy peak
                if limits[0] < mz + heavyShift/self.charge < limits[1]:
                    includedPeaks.append('h'+str(peakNumber))
            mz += tools.Cmass/self.charge
            peakNumber+=1
            if peakNumber == maxPeaksInCluster:
                break
        
        if not coIsolation:
            #heavy peaks are the same as light if not co-isolated:
            heavyIncludedPeaks = []
            for peak in includedPeaks:
                heavyIncludedPeaks.append('h'+peak[-1])
            includedPeaks += heavyIncludedPeaks
        return includedPeaks
            
    def add_yFragments(self, maxCharge, maxClusterPeaks, labels, IAA = False):
        """ Adds y-fragments to the fragment dictionary when called.
        
        >>> def inrange(number, reference): return reference-1E-5 < number < reference+1E-5 
        >>> a = peptide('PEPTIDECK', 1, 0, 10, 3, 1.3, 'K:0C2N, R:0C2N', 1, 3, 1)
        >>> a.add_yFragments(2, 3, 'K:0C2N, R:0C2N', IAA=True)
        >>> inrange(a.fragments['y8 l0+'], 991.4400844433401)
        True
        >>> inrange(a.fragments['y8 l0++'], 496.22368045505505)
        True
        >>> inrange(a.fragments['y9 l0+'], 1088.49284829219)
        True
        >>> 'y1 l0+' in a.fragments
        False
        >>> 'y2 l0+' in a.fragments
        True
        
        """
        for charge in range(1, maxCharge+1):
            for i in range(len(self.sequence)-1):
                sequence = self.sequence[i:]
                #Calculate cluster:
                isotopeCluster = tools.calculate_fragmentPeakCluster(sequence, charge, IAA, maxClusterPeaks, labels, 'y', self.includedPeaks)
                #Update fragment dictionary
                self.fragments.update(isotopeCluster)
                
    def add_bFragments(self, maxCharge, maxClusterPeaks, labels, IAA = False):
        """ Adds b-fragments to the fragment dictionary when called.
        
        >>> a = peptide('PEPTIDECK', 2, 0, 10, 3, 1.3, 'K:0C2N, R:0C2N', 1, 3, 1)
        >>> a.add_bFragments(2, 3, 'K:0C2N, R:0C2N', IAA=True)
        >>> def inrange(number, reference): return reference-1E-5 < number < reference+1E-5 
        >>> inrange(a.fragments['b8 0+'], 942.38732059449)
        True
        >>> inrange(a.fragments['b8 0++'], 471.69729853063)
        True
        >>> inrange(a.fragments['b9 0+'], 1070.48228360849)
        True
        >>> 'b1 0+' in a.fragments
        False
        >>> 'b2 0+' in a.fragments
        True
        """
        for charge in range(1, maxCharge+1):
            for i in range(0, len(self.sequence)-1):
                sequence = self.sequence[:len(self.sequence)-i]
                #Calculate cluster:
                isotopeCluster = tools.calculate_fragmentPeakCluster(sequence, charge, IAA, maxClusterPeaks, labels, 'b', self.includedPeaks)
                #Update fragment dictionary
                self.fragments.update(isotopeCluster)
    
    def add_fragmentRanges(self, tolerance):
        """ Calculates tuples of min and max m/z whithin which a peak will be 
            called for the specific peak type.
        
        >>> def inrange(number, reference): return reference-1E-5 < number < reference+1E-5 
        >>> a = peptide('PEPTIDECK', 1, 0, 10, 3, 1.3, 'K:0C2N, R:0C2N', 1, 3, 1)
        >>> a.add_yFragments(2, 3, 'K:0C2N, R:0C2N', IAA=True)
        >>> a.add_fragmentRanges(5)
        >>> range = a.fragmentRanges['y7 l0+']
        >>> inrange(range[0], 862.3931793679134)
        True
        >>> inrange(range[1], 862.4018033428268)
        True
        >>> range = a.fragmentRanges['y7 l0++']
        >>> inrange(range[0], 431.7002253991505)
        True
        >>> inrange(range[1], 431.7045424229896)
        True
        """
        self.fragmentRanges = {}
        for peak in self.fragments:
            mz = self.fragments[peak]
            halfSpan = mz*tolerance/1000000
            self.fragmentRanges[peak] = (mz-halfSpan, mz+halfSpan)
    
    def add_precursorRanges(self, tolerance):
        """ Calculates tuples of min and max m/z whithin which a peak will be 
            called for the specific peak type.
        """
        self.precursorRanges = {}
        for peak in self.precursorMass:
            mz = self.precursorMass[peak]
            halfSpan = mz*tolerance/1000000
            self.precursorRanges[peak] = (mz-halfSpan, mz+halfSpan)
            
    def __add_MS1peaks(self):
        """Placeholder for addition of all MS1 peaks.
        """
        pass
  
def fileParser(filePath):
    """ Creates a generator that outputs lines in forms of lists, splitted by
        commas and semicolons.
    """
    #Open file:
    with open(filePath, 'r') as file:
        for line in file:
            #Make list of csv row elements:
            line = re.split('\n|,|;', line.strip())
            
            yield line

def peptideSettingsImporter(isolationListPath, settings):
    """ Imports settings from and isolation list and returns a dictionary with
        scan period and charge for each peptide in isolation list.
    """
    peptideSettings = {}
    first = True
    for line in fileParser(isolationListPath):
        if first:
            #Create headers
            headers = line
            first = False
            continue
        
        if 'heavy' in line[headers.index('Comment')]:
            #Heavy peptide weights are automatically calculated. 
            #Addition of heavy peptides is therefore skipped.
            continue

        #Extract settings:
        sequence = line[headers.index('Comment')].split(' ')[0]
        scanStart = float(line[headers.index('Start [min]')])
        scanEnd = float(line[headers.index('End [min]')])
        charge = int(line[headers.index('CS [z]')])
        
        #Save settings as peptide objects:
        peptideSettings[sequence+'+'*charge] = peptide(sequence, charge, scanStart, scanEnd, settings['Isolation window size'], settings['Isolation window offset'], settings['Heavy labels'], settings['Co-isolation'], settings['Peaks in isotope clusters'], settings['IAA'])
        
    return peptideSettings

def massCalculator(isolationListPath, settings):
    """ Calculates mzs for precursors and fragments and returns a list of 
        peptide objects.
    >>> from setting_handler import templateLoader
    >>> isolationListPath = "../DocTests/DocTest IsoList.csv"
    >>> settings = templateLoader("../DocTests/standard settings.pkl").get_settings()
    
    >>> peptides = ['YGVSDYHK++', 'NLINNAK++', 'DLQAQVVESAK++', 'ISEATDGLSDFLK++', 'SQTPAEDTVK++', 'SIELAEAK++', 'YGVSDYYK++']
    >>> peptideSettings = massCalculator(isolationListPath, settings)
    >>> all([peptide in peptideSettings for peptide in peptides])
    True
    >>> peptideSettings['DLQAQVVESAK++'].precursorMass=={'h1': 595.8156720775, 'h0': 595.3139946586, 'l1': 594.8186371841, 'l2': 595.3203146029999, 'h2': 596.3173494964, 'l0': 594.3169597652}
    True
    """
    #Import sequence, scan period and charge from isolation list to dictionary
    #of peptide objects:
    peptideSettings = peptideSettingsImporter(isolationListPath, settings)
    
    #Do mass calculation for every peptide object in dictionary:
    for peptideObject in peptideSettings.values():
        #Add precursor mzs:
        '''
        peptideObject.calculate_precursorMass(settings['Peaks in isotope clusters'], 
                                     settings['Heavy labels'])
        '''
        #Add y-fragment mzs:
        peptideObject.add_yFragments(settings['Max fragment charge'],
                            settings['Peaks in isotope clusters'], 
                            settings['Heavy labels'],
                            IAA=settings['IAA'])
                            
        if settings['Include b-fragments']:
            #Add b-fragment mzs:
            peptideObject.add_bFragments(settings['Max fragment charge'],
                                settings['Peaks in isotope clusters'], 
                                settings['Heavy labels'],
                                IAA=settings['IAA']) 
        
        #Add ranges for peakCalling:
        peptideObject.add_precursorRanges(settings['MS1 tolerance'])
        peptideObject.add_fragmentRanges(settings['MS2 tolerance'])
            
    return peptideSettings


if __name__ == "__main__":
    import doctest
    doctest.testmod()