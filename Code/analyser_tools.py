# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 11:35:34 2017

@author: max.karlsson
"""

import collections


def coIsolatedPeakPairing(lightChannel, heavyChannel):
    """ Creates a dictionary with retention times as keys and tuples of light
        and heavy intensities (light intensity, heavy intensity) for each 
        retention time where there are both a heavy and light peak  from paired
        channels detected.
    
    >>> heavyChannel = [(26.784027, 76426.3046875), (26.838704, 71134.5546875), (26.896054, 39654.8125), (26.947822, 18284.69140625), (27.004714, 8455.0673828125), (27.059382, 8262.1708984375), (27.239498, 5855.65234375)]
    >>> lightChannel = [(26.784027, 76426.3046875), (26.838704, 71134.5546875), (26.896054, 39654.8125), (26.947822, 18284.69140625), (27.004714, 8455.0673828125), (27.059382, 8262.1708984375), (27.239498, 5855.65234375)]
    >>> coIsolatedPeakPairing(lightChannel, heavyChannel)
    [(76426.3046875, 76426.3046875), (71134.5546875, 71134.5546875), (39654.8125, 39654.8125), (18284.69140625, 18284.69140625), (8455.0673828125, 8455.0673828125), (8262.1708984375, 8262.1708984375), (5855.65234375, 5855.65234375)]
    """
    ratioChannel = []
    
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
            ratioChannel.append((lightIntensity, heavyIntensity))
    return ratioChannel
    
def asynchronalPeakPairing(lightChannel, heavyChannel):
    """ Creates a dictionary with retention times as keys and tuples of light
        and heavy intensities (light intensity, heavy intensity) where there 
        are both a heavy and light peak in range from paired channels.
        
    >>> heavyChannel = [(26.784027, 76426.3046875), (26.838704, 71134.5546875), (26.896054, 39654.8125), (26.947822, 18284.69140625), (27.004714, 8455.0673828125), (27.059382, 8262.1708984375), (27.239498, 5855.65234375)]
    >>> lightChannel = [(26.677597, 57905.5859375), (26.72487, 266683.28125), (26.781769, 735656.5), (26.836444, 803626.0), (26.893793, 410744.34375), (26.945562, 124779.6015625), (27.002453, 25164.228515625), (27.057121, 10775.19140625), (27.119456, 12546.4140625), (27.175771, 10831.6982421875), (27.237237, 9189.2158203125)]
    >>> asynchronalPeakPairing(lightChannel, heavyChannel)
    [(735656.5, 76426.3046875), (803626.0, 71134.5546875), (410744.34375, 39654.8125), (124779.6015625, 18284.69140625), (25164.228515625, 8455.0673828125), (10775.19140625, 8262.1708984375), (9189.2158203125, 5855.65234375)]

    >>> heavyChannel = [(26.32832, 6616.873046875), (26.357566, 12807.7353515625), (26.382969, 12354.1728515625), (26.406697, 10298.384765625), (26.431453, 11595.12109375), (26.784027, 28867.45703125), (26.838704, 18894.296875), (26.896054, 14050.2841796875)]
    >>> lightChannel = [(26.357566, 4081.56591796875), (26.677597, 29509.734375), (26.72487, 123921.0546875), (26.781769, 421344.34375), (26.836444, 419462.34375), (26.893793, 208354.171875), (26.945562, 57058.8046875), (27.057121, 5221.6513671875), (27.175771, 4774.64453125)]
    >>> asynchronalPeakPairing(lightChannel, heavyChannel)
    [(4081.56591796875, 12807.7353515625), (421344.34375, 28867.45703125), (419462.34375, 18894.296875), (208354.171875, 14050.2841796875)]

    
    """
    ratioChannel = []
    RTs = []
    
    #Return empty channel if input channels are empty:
    if len(lightChannel)<2 or len(heavyChannel)<2:
        return ratioChannel
    
    lightRTs = [x[0] for x in lightChannel]
    heavyRTs = [x[0] for x in heavyChannel]
    
    i = 0
    for heavyRT in heavyRTs:
        mindist = 1
        lightIntensity = None
        for lightRT in lightRTs[i:]:
            #Calculate distance in retention time:
            dist = abs(lightRT-heavyRT)
            #If the retention time is not within 1.8 seconds, continue:
            if dist > 0.03:
                continue
            #If the distance is smaller than the smallest recorded distance, 
            #choose this data point as the one to pair:
            if dist < mindist:
                i = lightRTs.index(lightRT)
                mindist = dist
                lightIntensity = lightChannel[i][1]
                chosenLightRT = lightRT
                
            else:
                break
        
        if not lightIntensity==None:
            heavyIntensity = heavyChannel[heavyRTs.index(heavyRT)][1]
            ratioChannel.append(((chosenLightRT, lightIntensity), (heavyRT, heavyIntensity)))
            RTs.append(chosenLightRT)
    
    #Remove duplicates that are furthest regarding RT:
    duplicateRTs = [item for item, count in collections.Counter(RTs).items() if count > 1]
    for RT in duplicateRTs:
        
        duplicateIndices = [i for i, pairedDataPoint in enumerate(ratioChannel) if pairedDataPoint[0][0]==RT]
        duplicateData = {}
        for index in duplicateIndices:
            duplicateData[index] = ratioChannel[index]
        
        mindist = 1
        for index, pairedDataPoint in duplicateData.items():
            dist = abs(pairedDataPoint[0][0]-pairedDataPoint[1][0])
            if dist<mindist:
                minIndex = index
                mindist = dist
        indicesToRemove = [index for index in duplicateIndices if not index == minIndex]
        
        ratioChannel = [value for i, value in enumerate(ratioChannel) if not i in indicesToRemove]
    
    #Remove RTs from ratioChannel, so it only includes paired intensities:
    ratioChannel = [(x[0][1],x[1][1]) for x in ratioChannel]
                
    return ratioChannel
    
if __name__=="__main__":
    import doctest
    doctest.testmod()