# -*- coding: utf-8 -*-
"""
Created on Fri Dec 30 15:30:56 2016

@author: Max Karlsson
"""
import time
from matplotlib.font_manager import FontProperties
import matplotlib.pyplot as pl
from Printicek import daviPrint


class savedPlot:
    """ A plot object that automatically saves when closed. Should be used with 
        the with statement.
    """
    def __init__(self, saveName, saveDirectory, figsize=(10,10), title=''):
        self.__figsize=figsize
        self.__saveFile = saveDirectory + '/' + saveName + '.png'
        self.__title = title
        self.subplots = {}
        
    def __enter__(self):
        self.__figure = pl.figure(figsize=self.__figsize)
        return self
        
    def __exit__(self, type, value, traceback):
        #Save and exit:
        pl.title(self.__title)
        self.__figure.savefig(self.__saveFile)
        pl.close(self.__figure)
    
    def __colors(self, colorIndex):
        """ A function that returns colors from an index. If index exceeds list
            indices, it loops.
        """
        colors = ['cornflowerblue', 'royalblue', 'darkcyan', 
        'orangered', 'darkturquoise', 'crimson', 'midnightblue', 
        'mediumvioletred', 'yellowgreen', 'slategray', 'darkorchid', 'gold', 
        'dimgray', 'purple', 'steelblue', 'brown', 'orchid', 'skyblue', 
        'black', 'turquoise', 'sandybrown', 'darksalmon', 'coral', 
        'darkorange', 'silver', 'darkmagenta', 'darkviolet', 'fuchsia', 
        'burlywood', 'darkslateblue', 'orange', 'darkolivegreen', 
        'forestgreen', 'rosybrown', 'mediumaquamarine', 'darkred', 
        'springgreen', 'chocolate', 'mediumpurple', 'teal', 'blue', 'sienna', 
        'deepskyblue', 'hotpink', 'dodgerblue', 'seagreen', 'powderblue', 
        'slateblue', 'gray', 'violet', 'mediumorchid', 'blueviolet', 'cyan', 
        'indigo', 'cadetblue', 'mediumblue', 'salmon', 'saddlebrown', 
        'darkslategray', 'darkgreen', 'maroon', 'darkblue', 'tomato', 
        'mediumturquoise', 'navy', 'green', 'mediumslateblue', 'peru']
        while True:
            if colorIndex>=len(colors):
                colorIndex-=len(colors)
            else:
                break
        return colors[colorIndex]
        
    def add_subplot(self, plotNumber, tall=1, wide=1, xlabel='', ylabel=''):
        """ Adds a subplot to the figure. 
        """
        self.subplots[plotNumber] = self.__figure.add_subplot(str(tall) + str(wide) + str(plotNumber))
        self.subplots[plotNumber].set_xlabel(xlabel)
        self.subplots[plotNumber].set_ylabel(ylabel)
    
    def add_dataToSubplot(self, x, y, plotNumber, label, colorIndex, marker='-'):
        self.subplots[plotNumber].plot(x, y, marker, label=label, color=self.__colors(colorIndex))

    def add_legend(self, plotNumber, location='upper right'):
        fontP = FontProperties()
        fontP.set_size('xx-small')
        handles, labels = self.subplots[plotNumber].get_legend_handles_labels()
        #If there are legend labels, order them:
        if not len(labels)==0:
            labels, handles = zip(*sorted(zip(labels, handles)))
        #Sort the label and handles in the legend:
        #labels, handles = zip(*sorted(zip(labels, handles)))
        self.__figure.legend(handles, labels, loc=location, prop = fontP)
            

def calculate_colorIndex(fragment, settings):
    """ Calculates the colorIndex so that each fragment gets a unique color 
        based on its peak number and fragment number.
    """
    fragment = fragment.strip('+')
    if len(fragment)==2 or "Prec" in fragment:
        colorIndex = int(fragment[-1])
    else:
        colorIndex = settings['Peaks in isotope clusters']
        fragmentNumber = int(fragment.split(' ')[0][1:])
        colorIndex += settings['Peaks in isotope clusters']*fragmentNumber
        colorIndex += int(fragment[-1])
    return colorIndex

def create_figures(pairedData, sortedData, linearRegressionParameters, saveName, saveDirectory, settings):
    """ Creates a figure containing the MS chromatogram, a linear regression
        plot fragment by fragment, and a linear regression of all fragments'
        data points.
    """
    saveName = '-'.join((time.strftime('%Y-%m-%d_%H-%M'),saveName))
            
    for peptide in pairedData:
        #Make figurename:
        figureName = '-'.join((saveName, peptide))
        #Create figure:
        with savedPlot(figureName, saveDirectory, figsize=(20,10), title=peptide) as figure:
            daviPrint(peptide)
            #Create chromatogram:
            figure.add_subplot(0, tall=2, wide=1, xlabel='Retention times [min]', ylabel='Intensity')
            #Create fragment by fragment linear regression plot:
            figure.add_subplot(1, tall=2, wide=2, xlabel='Heavy intensity', ylabel='Light intensity')
            #Create total linear regression plot:
            figure.add_subplot(2, tall=2, wide=2, xlabel='Heavy intensity', ylabel='Light intensity')
            #Add data to chromatogram:
            for fragment, channel in zip(sortedData[peptide], sortedData[peptide].values()):
                #Do not include MS1 data in this plot:
                if len(fragment)==2:
                    continue
                #Choose marker depending on heavy or light fragment:
                if 'h' in fragment:
                    marker='-d'
                else:
                    marker='-o'
                #Calculate color index:
                colorIndex = calculate_colorIndex(fragment, settings)
                RTs, intensities = zip(*channel)
                #Add points to plot:
                figure.add_dataToSubplot(RTs, intensities, 0, fragment, colorIndex, marker=marker)
            
            
            #Add data to fragment by fragment linear regression:
            for fragment, channel in zip(pairedData[peptide], pairedData[peptide].values()):
                #Calculate color index:
                colorIndex = calculate_colorIndex(fragment, settings)
                lightInts, heavyInts = zip(*channel)
                #At least 3 data points are needed for quantification:
                if len(lightInts)>2:
                    #Add points to plot:
                    figure.add_dataToSubplot(heavyInts, lightInts, 1, fragment, colorIndex, marker='.')
                    #Plot linear regression:
                    fittedValues = linearRegressionParameters[peptide][fragment]['fitted values']
                    figure.add_dataToSubplot(heavyInts, fittedValues, 1, '', colorIndex, marker='-')
                    
            #Add data to total linear regression:
            allLightInts, allHeavyInts = [],[]
            for fragment, channel in zip(pairedData[peptide], pairedData[peptide].values()):
                lightInts, heavyInts = zip(*channel)
                allLightInts += lightInts
                allHeavyInts += heavyInts
            
            #At least 3 data points are needed for quantification:
            if len(allLightInts)>2:
                #Add points to plot:
                figure.add_dataToSubplot(allHeavyInts, allLightInts, 2, '', 0, marker='.')
                #Plot linear regression:
                fittedValues = linearRegressionParameters[peptide]['Total']['fitted values']
                
                figure.add_dataToSubplot(allHeavyInts, fittedValues, 2, '', 0, marker='-')
            
            #Add legend to figure:
            figure.add_legend(0)
            figure.add_legend(1, location='upper left')