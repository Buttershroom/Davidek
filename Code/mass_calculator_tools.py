# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 11:41:59 2017

@author: max.karlsson

"""
from pyteomics import mass

#Carbon isotopic label weight:
Cmass = 13.0033548378 - 12
Nmass = 15.0001088982 - 14.0030740048
#IAA added weight:
IAAmass = 57.02146
    
    
def calculate_mass(sequence, charge, ionType, IAA):
    """ Calculates and returns the monoisotopic mass of an ion.
    
    >>> calculate_mass('PEPTIDECK', 1, 'y', 1)
    1088.4928482921898
    >>> calculate_mass('PEPTIDECK', 1, 'b', 1)
    1070.48228360849
    >>> calculate_mass('PEPTIDECK', 1, 'y', 0)
    1031.4713882921899
    >>> calculate_mass('PEPTIDECK', 2, 'y', 1)
    544.7500623794799
    """
    monoisotopicMass = mass.calculate_mass(sequence, charge=charge, ion_type=ionType)
    if IAA: 
        monoisotopicMass += IAAmass*sequence.count('C')/charge
    
    return monoisotopicMass
    
def calculate_heavyShift(sequence, labels):
    """ Calculates and returns the mass shift of a heavy peptide given the sequence.
    
    >>> calculate_heavyShift("PEPTIDECK", 'K:0C2N, R:0C2N')
    1.9940697868000008
    """
    #Count number of K and R in sequence:
    Knum = sequence.count('K')
    Rnum = sequence.count('R')
    #Calculate mass shift:
    Kmass = int(labels.split(',')[0].split(':')[1][0])*Cmass + int(labels.split(',')[0].split(':')[1][2])*Nmass
    Rmass= int(labels.split(' ')[1].split(':')[1][0])*Cmass + int(labels.split(' ')[1].split(':')[1][2])*Nmass
    return Kmass*Knum + Rmass*Rnum
    
def calculate_fragmentPeakCluster(sequence, charge, IAA, maxClusterPeaks, labels, ionType, includedPeaks):
    """ Calculates and returns a dictionary of the peak cluster of a certain fragment.
    
    >>> calculate_fragmentPeakCluster('PEPTIDECK', 1, 1, 3, 'K:0C2N, R:0C2N', 'y', ['l0','l1','l2','h0'])=={'y9 l1+': 1089.49620312999, 'y9 l0+': 1088.4928482921898, 'y9 l2+': 1090.4995579677898, 'y9 h0+': 1090.48691807899}
    True
    
    """
    monoisotopicMass = calculate_mass(sequence, charge, ionType, IAA)
    heavyShift = calculate_heavyShift(sequence, labels)
    fragmentNumber = len(sequence)
    
    fragmentCluster = {}
    
    lightID = ' l'
    #b-fragments are the same mass whether from heavy or light origins:
    if ionType=='b':
        lightID=' '
    
    #Parse through peaks in isotopic order:
    for peakNumber in range(maxClusterPeaks):
        #Calculate light peak m/z:
        
        if 'l'+str(peakNumber) in includedPeaks:
            lightPeakMass = monoisotopicMass + peakNumber*Cmass/charge
            fragmentCluster[ionType + str(fragmentNumber) + lightID +str(peakNumber) + '+'*charge] = lightPeakMass
        #Calculate heavy peak m/z. Not for 'b'-ions
        if not ionType=='b':
            if 'h'+str(peakNumber) in includedPeaks:
                heavyPeakMass = lightPeakMass + heavyShift/charge
                fragmentCluster[ionType + str(fragmentNumber) + ' h' +str(peakNumber) + '+'*charge] = heavyPeakMass
    
    return fragmentCluster
 
 

def calculate_precursorPeakCluster(sequence, charge, maxClusterPeaks, labels, IAA=False):
    """ Calculate and returns the m/z of the precursor's isotope cluster.
    
    #No IAA/CAA:
    >>> calculate_precursorPeakCluster("PEPTIDECK", 1, 3, 'K:0C2N, R:0C2N', IAA=False)=={'l0': 1031.4713882921899, 'l1': 1032.47474312999, 'h1': 1034.46881291679, 'h2': 1035.47216775459, 'l2': 1033.4780979677898, 'h0': 1033.46545807899}
    True
    
    #With IAA/CAA:
    >>> calculate_precursorPeakCluster("PEPTIDECK", 1, 3, 'K:0C2N, R:0C2N', IAA=True)=={'l0': 1088.4928482921898, 'l1': 1089.49620312999, 'h1': 1091.49027291679, 'h2': 1092.4936277545899, 'l2': 1090.4995579677898, 'h0': 1090.48691807899}
    True
    
    #With IAA/CAA and charge 2:
    >>> calculate_precursorPeakCluster("PEPTIDECK", 2, 3, 'K:0C2N, R:0C2N', IAA=True)=={'l0': 544.7500623794799, 'l1': 545.25173979838, 'h1': 546.24877469178, 'h2': 546.7504521106799, 'l2': 545.7534172172799, 'h0': 545.74709727288}
    True
    
    """
    precursorCluster = {}
    #Calculate the monoisotopic mass over charge:
    monoisotopicMass = calculate_mass(sequence, charge, 'y', IAA)
    #Calculate peak mass over charges for every peaknumber:
    for peakNumber in range(maxClusterPeaks):
        #Calculate light peak m/z:
        lightPeakMass = monoisotopicMass + peakNumber*Cmass/charge
        precursorCluster['l'+str(peakNumber)] = lightPeakMass
        #Calculate heavy peak m/z
        heavyShift = calculate_heavyShift(sequence, labels)
        heavyPeakMass = lightPeakMass + heavyShift/charge
        precursorCluster['h'+str(peakNumber)] = heavyPeakMass       
    
    return precursorCluster
    

if __name__ == '__main__':
    import doctest
    doctest.testmod()