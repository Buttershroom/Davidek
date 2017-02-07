# -*- coding: utf-8 -*-
"""
Created on Sat Dec 10 15:28:30 2016

@author: Max Karlsson
"""
import re
from pyteomics import mass

class peptide:
    """ Peptide object containing peptide specific variables such as sequence,
        mass, scan period, fragment masses.
    """
    def __init__(self, sequence, charge, scanStart, scanEnd, isoWinSize, isoWinOffset, labels, coIsolation, maxPeaksInCluster):
        self.sequence = sequence
        self.charge = charge
        self.scanPeriod = (scanStart, scanEnd)
        self.fragments = {}
        self.includedPeaks = self.__calculate_included_peaks(isoWinSize, isoWinOffset, labels, coIsolation, maxPeaksInCluster)
    
    def __calculate_included_peaks(self, isoWinSize, isoWinOffset, labels, coIsolation, maxPeaksInCluster):
        """ Function to calculate which peaks in a cluster are isolated by the 
            isolation window size and offset.
        
        #Co-isolated:
        >>> a = peptide('PEPTIDECK', 1, 0, 10, 3, 1.3, 'K:0C2N, R:0C2N', 1, 3)
        >>> a.includedPeaks
        ['l0', 'h0', 'l1', 'l2']
        
        #Not co-isolated:
        >>> a = peptide('PEPTIDECK', 1, 0, 10, 3, 1.3, 'K:0C2N, R:0C2N', 0, 3)
        >>> a.includedPeaks
        ['l0', 'l1', 'l2', 'h0', 'h1', 'h2']
        """
        Cmass = 13.0033548378 - 12
        includedPeaks = []
        limits = (isoWinOffset-isoWinSize/2, isoWinOffset+isoWinSize/2)
        heavyShift = self.__calculate_heavyShift(self.sequence, labels)
        
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
            mz += Cmass/self.charge
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
    
    def __calculate_heavyShift(self, sequence, labels):
        """ Calculate the mass shift of heavy peptides.
        """
        Cmass = 13.0033548378 - 12
        Nmass = 15.0001088982 - 14.0030740048
        #Count number of K and R in sequence:
        Knum = sequence.count('K')
        Rnum = sequence.count('R')
        #Calculate mass shift:
        Klabel = int(labels.split(',')[0].split(':')[1][0])*Cmass + int(labels.split(',')[0].split(':')[1][2])*Nmass
        Rlabel = int(labels.split(' ')[1].split(':')[1][0])*Cmass + int(labels.split(' ')[1].split(':')[1][2])*Nmass
        return Klabel*Knum + Rlabel*Rnum
    
    def __monoisotopicMass(self, sequence, IAA, charge, ionType):
        monoisotopicMass = mass.fast_mass(sequence, charge=charge, ion_type=ionType)
        if IAA: monoisotopicMass += sequence.count('C')*57.02146/charge
        return monoisotopicMass
    
    def calculate_precursorMass(self, maxClusterPeaks, labels, IAA=False):
        """ Calculate the masses of the precursor's isotope cluster.
        
        #No IAA/CAA:
        >>> a = peptide('PEPTIDECK', 1, 0, 10, 3, 1.3, 'K:0C2N, R:0C2N', 1, 3)
        >>> a.calculate_precursorMass(3, 'K:0C2N, R:0C2N', IAA=False)
        >>> a.precursorMass['l0']
        1031.47137115047
        >>> a.precursorMass['l1']
        1032.4747259882702
        >>> a.precursorMass['h0']
        1033.4654409372702
        
        #With IAA/CAA:
        >>> a = peptide('PEPTIDECK', 1, 0, 10, 3, 1.3, 'K:0C2N, R:0C2N', 1, 3)
        >>> a.calculate_precursorMass(3, 'K:0C2N, R:0C2N', IAA=True)
        >>> a.precursorMass['l0']
        1088.49283115047
        
        #With IAA/CAA and charge 2:
        >>> a = peptide('PEPTIDECK', 2, 0, 10, 3, 1.3, 'K:0C2N, R:0C2N', 1, 3)
        >>> a.calculate_precursorMass(3, 'K:0C2N, R:0C2N', IAA=True)
        >>> a.precursorMass['l0']
        544.75005380862
        
        """
        Cmass = 13.0033548378 - 12
        self.precursorMass = {}
        #Calculate the monoisotopic mass over charge:
        monoisotopicMass = self.__monoisotopicMass(self.sequence, IAA, self.charge, 'y')
        #Calculate peak mass over charges for every peaknumber:
        for peakNumber in range(maxClusterPeaks):
            #Calculate light peak m/z:
            lightPeakMass = monoisotopicMass + peakNumber*Cmass/self.charge
            self.precursorMass['l'+str(peakNumber)] = lightPeakMass
            #Calculate heavy peak m/z
            heavyShift = self.__calculate_heavyShift(self.sequence, labels)
            heavyPeakMass = lightPeakMass + heavyShift/self.charge
            self.precursorMass['h'+str(peakNumber)] = heavyPeakMass       
    
    
    def __calculate_fragment_masses(self, sequence, charge, IAA, maxClusterPeaks, labels, ionType, fragmentNumber):
        """ Calculate masses of a specified fragment's isotope cluster.
        """
        Cmass = 13.0033548378 - 12
        monoisotopicMass = self.__monoisotopicMass(sequence, IAA, charge, ionType)
        fragmentMass = {}
        lightID = ' l'
        if ionType=='b':
            lightID=' '
        for peakNumber in range(maxClusterPeaks):
            #Calculate light peak m/z:
            
            if 'l'+str(peakNumber) in self.includedPeaks:
                lightPeakMass = monoisotopicMass + peakNumber*Cmass/charge
                fragmentMass[ionType + str(fragmentNumber) + lightID +str(peakNumber) + '+'*charge] = lightPeakMass
            #Calculate heavy peak m/z. Not for 'b'-ions
            if not ionType=='b':
                if 'h'+str(peakNumber) in self.includedPeaks:
                    heavyShift = self.__calculate_heavyShift(sequence, labels)
                    heavyPeakMass = lightPeakMass + heavyShift/charge
                    fragmentMass[ionType + str(fragmentNumber) + ' h' +str(peakNumber) + '+'*charge] = heavyPeakMass
            
        return fragmentMass
            
    
            
    def add_yFragments(self, maxCharge, maxClusterPeaks, labels, IAA = False):
        """ Adds y-fragments to the fragment dictionary when called.
        
        >>> a = peptide('PEPTIDECK', 1, 0, 10, 3, 1.3, 'K:0C2N, R:0C2N', 1, 3)
        >>> a.add_yFragments(2, 3, 'K:0C2N, R:0C2N', IAA=True)
        >>> a.fragments['y8 l0+']
        991.44007115047
        >>> a.fragments['y8 l0++']
        496.22367380862
        >>> a.fragments['y9 l0+']
        1088.49283115047
        >>> 'y1 l0+' in a.fragments
        False
        >>> 'y2 l0+' in a.fragments
        True
        
        """
        for charge in range(1, maxCharge+1):
            for i in range(len(self.sequence)-1):
                fragmentNumber = len(self.sequence)-i
                sequence = self.sequence[i:]
                #Calculate cluster:
                isotopeCluster = self.__calculate_fragment_masses(sequence, charge, IAA, maxClusterPeaks, labels, 'y', fragmentNumber)
                #Update fragment dictionary
                self.fragments.update(isotopeCluster)
                
    def add_bFragments(self, maxCharge, maxClusterPeaks, labels, IAA = False):
        """ Adds b-fragments to the fragment dictionary when called.
        
        >>> a = peptide('PEPTIDECK', 2, 0, 10, 3, 1.3, 'K:0C2N, R:0C2N', 1, 3)
        >>> a.add_bFragments(2, 3, 'K:0C2N, R:0C2N', IAA=True)
        
        >>> a.fragments['b8 0+']
        942.38730646677
        >>> a.fragments['b8 0++']
        471.69729146677
        >>> a.fragments['b9 0+']
        1070.4822664667702
        >>> 'b1 0+' in a.fragments
        False
        >>> 'b2 0+' in a.fragments
        True
        """
        for charge in range(1, maxCharge+1):
            for i in range(0, len(self.sequence)-1):
                fragmentNumber = len(self.sequence)-i
                sequence = self.sequence[:fragmentNumber]
                #Calculate cluster:
                isotopeCluster = self.__calculate_fragment_masses(sequence, charge, IAA, maxClusterPeaks, labels, 'b', fragmentNumber)
                #Update fragment dictionary
                self.fragments.update(isotopeCluster)
    
    def add_fragmentRanges(self, tolerance):
        """ Calculates tuples of min and max m/z whithin which a peak will be 
            called for the specific peak type.
        
        >>> a = peptide('PEPTIDECK', 1, 0, 10, 3, 1.3, 'K:0C2N, R:0C2N', 1, 3)
        >>> a.add_yFragments(2, 3, 'K:0C2N, R:0C2N', IAA=True)
        >>> a.add_fragmentRanges(5)
        >>> a.fragmentRanges['y7 l0+']
        (862.3931691630642, 862.4017931378758)
        >>> a.fragmentRanges['y7 l0++']
        (431.70022029672594, 431.704537320514)
        
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
        peptideSettings[sequence+'+'*charge] = peptide(sequence, charge, scanStart, scanEnd, settings['Isolation window size'], settings['Isolation window offset'], settings['Heavy labels'], settings['Co-isolation'], settings['Peaks in isotope clusters'])
            
    return peptideSettings

def massCalculator(isolationListPath, settings):
    """ Calculates mzs for precursors and fragments and returns a list of 
        peptide objects.
    """
    #Import sequence, scan period and charge from isolation list to dictionary
    #of peptide objects:
    peptideSettings = peptideSettingsImporter(isolationListPath, settings)
    
    #Do mass calculation for every peptide object in dictionary:
    for peptideObject in peptideSettings.values():
        #Add precursor mzs:
        peptideObject.calculate_precursorMass(settings['Peaks in isotope clusters'], 
                                     settings['Heavy labels'])
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