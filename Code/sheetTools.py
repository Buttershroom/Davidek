# -*- coding: utf-8 -*-
"""
Created on Fri Dec 30 15:24:02 2016

@author: Max Karlsson
"""
import xlsxwriter

class outputFile:
    """ An object containing the outdata savefile.
    """
    def __init__(self, saveName, saveDirectory):
        self.__saveName = saveName
        self.__saveDir = saveDirectory
        self.worksheets = {}
    
    def __enter__(self):
        #Creates the file:
        self.wb = xlsxwriter.Workbook(self.__saveDir + '/output/' + self.__saveName + '.xlsx')
        return self
        
    def __exit__(self, type, value, traceback):
        #Closes the file:
        self.wb.close()
        
    def add_worksheet(self, worksheetName):
        worksheet = self.wb.add_worksheet(name=worksheetName)
        self.worksheets[worksheetName] = worksheet
        
    def write_column(self, iterator, worksheetName, row, column):
        for element in iterator:
            worksheet = self.worksheets[worksheetName]
            try:
                worksheet.write(row, column, element)
            except:
                worksheet.write(row, column, 'NA')
            row += 1
            
    def write_row(self, iterator, worksheetName, row, column):
        for element in iterator:
            worksheet = self.worksheets[worksheetName]
            worksheet.write(row, column, element)
            column += 1
            
def sortFragments(fragments):
    """ Sorts list of fragments to have precursor peaks first, and fragment 
        peaks last in order of fragment number. The total regression is not
        returned.
    """
    #Add precursors:
    sortedFragments = sorted([fragment for fragment in fragments if 'Prec' in fragment])
    #Add y-fragments:
    sortedFragments += sorted([fragment for fragment in fragments if 'y' in fragment])
    #Add b-fragments:
    sortedFragments += sorted([fragment for fragment in fragments if 'b' in fragment])
    
    return sortedFragments
    
def make_resultSheet(workbook, saveName, linearRegressionParameters):
    """ Saves a sheet in the result file with a quantification summary of the
        data file.
    """
    LRP = linearRegressionParameters
    peptides = LRP.keys()
    #Extract output columns:
    col1, col2, col3, col4 = (peptides, 
                            [LRP[peptide]['Total']['k'] for peptide in peptides],
                            [LRP[peptide]['Total']['m'] for peptide in peptides],
                            [LRP[peptide]['Total']['r2'] for peptide in peptides])
    #Create sheet:
    workbook.add_worksheet(saveName)
    #Create  titles:
    workbook.write_row(('Peptide', 'k', 'm', 'r2'), saveName, 0, 0)
    #Write peptides:
    workbook.write_column(col1, saveName, 1, 0)
    #Write ratios:
    workbook.write_column(col2, saveName, 1, 1)
    #Write intersects:
    workbook.write_column(col3, saveName, 1, 2)
    #Write r2 values:
    workbook.write_column(col4, saveName, 1, 3)
    
def make_detailSheet(workbook, saveName, linearRegressionParameters):
    """ Saves a sheet in the result file with quantification data of all the
        fragments quantified in the data file.
    """
    #Create sheet:
    saveName = saveName+'_details'
    workbook.add_worksheet(saveName)
    LRP = linearRegressionParameters
    
    table = []
    #Create output table:
    for peptide in LRP:
        sortedFragments = sortFragments(LRP[peptide].keys())
        for fragment in sortedFragments:
            k, m, r2 = (LRP[peptide][fragment]['k'], 
                        LRP[peptide][fragment]['m'], 
                        LRP[peptide][fragment]['r2'])
            table.append((peptide, fragment, k, m, r2))
    #If no values are found, don't write anything:
    if table==[]:
        return
    #Extract output columns:
    col1, col2, col3, col4, col5 = zip(*table)
    #Create  titles:
    workbook.write_row(('Peptide', 'Peak', 'k', 'm', 'r2'), saveName, 0, 0)
    
    #Write peptides:
    workbook.write_column(col1, saveName, 1, 0)
    #Write peaks:
    workbook.write_column(col2, saveName, 1, 1)
    #Write ratios:
    workbook.write_column(col3, saveName, 1, 2)
    #Write intersects:
    workbook.write_column(col4, saveName, 1, 3)
    #Write r2 values:
    workbook.write_column(col5, saveName, 1, 4)

def make_summarySheet(workbook, allLinRegParams, sheetNames):
    """ Saves a sheet in the result file with a quantification summary of all
        data files included in the experiment.
    """
    #Create sheet:
    saveName = 'Summary'
    workbook.add_worksheet(saveName)
    #Get full filenames:
    fileNames = dict( (v,k.split('/')[-1]) for k, v in sheetNames.items() )
    
    table = []
    #Create output table:
    for fileID in allLinRegParams:
        fileName = fileNames[fileID]
        for peptide in allLinRegParams[fileID]:
            k, m, r2 = (allLinRegParams[fileID][peptide]['Total']['k'], 
                        allLinRegParams[fileID][peptide]['Total']['m'], 
                        allLinRegParams[fileID][peptide]['Total']['r2'])
            table.append((fileName, peptide, k, m, r2))
    #If table is empty, skip
    if table==[]:
        return
    #Extract output columns:
    col1, col2, col3, col4, col5 = zip(*table)
    
    #Create  titles:
    workbook.write_row(('Filename', 'Peptide', 'k', 'm', 'r2'), saveName, 0, 0)
    #Write file IDs
    workbook.write_column(col1, saveName, 1, 0)
    #Write peptides:
    workbook.write_column(col2, saveName, 1, 1)
    #Write ratios:
    workbook.write_column(col3, saveName, 1, 2)
    #Write intersects:
    workbook.write_column(col4, saveName, 1, 3)
    #Write r2 values:
    workbook.write_column(col5, saveName, 1, 4)
    

