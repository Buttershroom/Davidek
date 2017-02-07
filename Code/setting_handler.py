# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 13:04:46 2016

@author: Max Karlsson
"""
from Printicek import daviPrint
from pickle import dump, load, HIGHEST_PROTOCOL
import os
import tkinter as tk
from tkinter.filedialog import askopenfilename, askopenfilenames, asksaveasfilename, askdirectory
import sys
import textwrap

def templateSaver(settings, savename):
    """ Saves a settings variable to a pickled file in ./Settings templates.
    """
    #Create save directory if it does not exist:
    if not os.path.exists(os.path.join(os.getcwd(),'Setting templates\\')):
        os.mkdir(os.path.join(os.getcwd(),'Setting templates\\'))
    
    with open(savename +'.pkl', 'wb') as savefile:
        dump(settings, savefile, HIGHEST_PROTOCOL)
    daviPrint('Settings file "' + savename + '.pkl" was saved.')
    
def templateLoader(loadname):
    """ Loads a settings variable from a pickled file in ./Settings templates.
    """
    with open(loadname, 'rb') as loadfile:
        settings = load(loadfile)
    
        return settings

class GUI:
    """ A GUI object that can have buttons, fields and labels added.
    """
    def __init__(self, settings):
        self.settings = settings
        self.top = tk.Tk()
        self.top.title('Davidek')
        self.entries = {}
        self.buttons = {}
        self.isolationListFilename = ''
        self.dataFilenames = ''
        self.saveDir = ''
        self.expFilenames = []
        
    def add_separator(self, row, column):
        tk.Label(self.top,text='_____________________________________').grid(row=row, column=column, sticky=tk.NW)
        
    def add_label(self, text, row, column, bold=False):
        """ Adds a label to the GUI when called.
        """
        if bold:
            tk.Label(self.top, text=text, justify=tk.LEFT, font=('bold')).grid(row=row, column=column, sticky=tk.W)
        else:
            tk.Label(self.top, text=text, justify=tk.LEFT).grid(row=row, column=column, sticky=tk.W)
    
    #Inputs:
    def add_checkbutton(self, setting, text, state, row, column):
        """ Adds a checkbutton to the GUI when called.
        """
        #Create variable to store user input:
        var = tk.IntVar(self.top, value=int(state)) 
        #Create and place button:
        self.button = tk.Checkbutton(self.top, text=text, variable=var, onvalue=True, offvalue=False)
        self.button.grid(row=row, column=column, sticky=tk.W)
        #Save button variable in dictionary:
        self.buttons[setting] = var
    
    def add_entry(self, setting, state, row, column):
        """ Adds an entry field to the GUI when called.
        """
        #Create variable to store user input:
        var = tk.StringVar(self.top, value=str(state))
        #Create and place entry field:
        entry = tk.Entry(self.top, textvariable = var)
        entry.grid(row=row, column=column, sticky=tk.W)
        #Save entry field object in dictionary
        self.entries[setting] = var
    
    #Buttons:
    def add_experimentNameField(self, row, column):
        """ Adds an entry field to the GUI when called.
        """
        #Create variable to store user input:
        self.experimentName = tk.StringVar(self.top, value='')
        #Create and place entry field:
        self.nameEntry = tk.Entry(self.top, textvariable = self.experimentName)
        self.add_label('Experiment name:', row, column, bold=True)
        self.nameEntry.grid(row=row+1, column=column)
        
    def add_destroy_button(self, text, row, column):
        """ Adds a button to the GUI when called.
        """
        def destroy():
            #Close window and program:
            self.top.destroy()
            sys.exit()
        button = tk.Button(self.top, text=text, command=destroy)
        button.grid(row=row, column=column)
        
    def add_OK_button(self, text, row, column):
        """ Adds a button to the GUI when called.
        """
        def destroy():
            #Collect settings and close window:
            self.settingsCollector()
            self.top.destroy()
        button = tk.Button(self.top, text=text, command=destroy)
        button.grid(row=row, column=column)
        
    def add_browse_isolationList_button(self, row, column):
        """ Adds a browse button to the GUI when called.
        """
        self.isolationList_labelText = tk.StringVar('')
        tk.Label(self.top, textvariable=self.isolationList_labelText, justify=tk.LEFT).grid(row=row+1, column=column)
        button = tk.Button(self.top, text='Browse isolation list', command=self.browseIsolationList)
        button.grid(row=row, column=column)
        
    def add_browse_dataFiles_button(self, row, column):
        """ Adds a browse button to the GUI when called.
        """
        self.dataFile_labelText = tk.StringVar('')
        tk.Label(self.top, textvariable=self.dataFile_labelText, justify=tk.LEFT).grid(row=row+1, column=column)
        button = tk.Button(self.top, text='Browse data files', command=self.browseDataFiles)
        button.grid(row=row, column=column)
    
    def add_browse_saveDir_button(self, row, column):
        """ Adds a browse button to the GUI when called.
        """
        self.saveDir_labelText = tk.StringVar('')
        tk.Label(self.top, textvariable=self.saveDir_labelText, justify=tk.LEFT).grid(row=row+1, column=column)
        button = tk.Button(self.top, text='Browse save directory', command=self.browseSaveDir)
        button.grid(row=row, column=column)
    
    def add_browse_experimentFiles_button(self, row, column):
        """ Adds a browse button to the GUI when called.
        """
        self.exp_labelText = tk.StringVar('')
        self.exp_labelText2 = tk.StringVar('')
        tk.Label(self.top, textvariable=self.exp_labelText, justify=tk.LEFT).grid(row=row+1, column=column)
        tk.Label(self.top, textvariable=self.exp_labelText2, justify=tk.LEFT, fg='red').grid(row=row+2, column=column)
        button = tk.Button(self.top, text='Browse experiment files', command=self.browseExperimentFiles)
        button.grid(row=row, column=column)
    
    def browseExperimentFiles(self):
        """ Function to browse for .mzML data files.
        """
        self.expFilenames = askopenfilenames(filetypes=[('Experiment file', '*.pkl')])
        self.exp_labelText.set(str(len(self.expFilenames)) + ' experiment files selected.' )
        self.exp_labelText2.set('No peak calling will be performed')
        
    def browseIsolationList(self):
        """ Function to browse for .csv isolation lists.
        """
        self.isolationListFilename = askopenfilename(filetypes=[('csv file','*.csv'), ('All files','*.*')])
        self.isolationList_labelText.set('Isolation list: ' + self.isolationListFilename.split('/')[-1])
        
    def browseDataFiles(self):
        """ Function to browse for .mzML data files.
        """
        self.dataFilenames = askopenfilenames(filetypes=[('mzML file', '*.mzML')])
        self.dataFile_labelText.set(str(len(self.dataFilenames)) + ' data files selected.' )
        
    def browseSaveDir(self):
        """ Function to browse for .mzML data files.
        """
        self.saveDir = askdirectory()
        self.saveDir_labelText.set('Save directory: .../' + '/'.join(str(self.saveDir).split('/')[-2:]))
    
    #Help window
    def add_help_button(self, row, column):
        """ Adds a browse button to the GUI when called
        """
        button = tk.Button(self.top, text='What. Help.', command=self.displayHelp)
        button.grid(row=row, column=column)
        
    def displayHelp(self):
        """ Function to display help window.
        """

        helpWindow = tk.Toplevel()
        for i, setting in enumerate(self.settings.settings.values()):
            
            name = setting.name
            description = '\n'.join(textwrap.wrap(setting.description, 50))
            
            tk.Label(helpWindow, text=name, justify=tk.LEFT).grid(row=i, column=0, sticky=tk.NW)
            tk.Label(helpWindow, text=description, justify=tk.LEFT).grid(row=i, column=1, sticky=tk.NW)
        
    #Saving and loading settings:
    def add_browse_settingFiles_button(self, row, column):
        """ Adds a browse button to the GUI when called.
        """
        button = tk.Button(self.top, text='Browse setting files', command=self.browseSettingFile)
        button.grid(row=row, column=column)
        
    def browseSettingFile(self):
        """ Function to browse for .pkl settings files.
        """
        #Open file:
        self.settingsFilename = askopenfilename(filetypes=[('settings file','*.pkl')])
        self.settings.settings = templateLoader(self.settingsFilename)
        
        #Update settings data:
        for key in self.entries:
            var = self.entries[key]
            #Supposed to update fields:
            var.set(self.settings.settings[key].value)
        for key in self.buttons:
            var = self.buttons[key]
            #Supposed to update button states:
            var.set(self.settings.settings[key].value)
        
    
    def add_save_settings_button(self, row, column):
        """ Adds a browse button to the GUI when called
        """
        button = tk.Button(self.top, text='Save settings as...', command=self.saveSettingFile)
        button.grid(row=row, column=column)
        
    def saveSettingFile(self):
        """ Function to save settings file. Asks user for save file name.
        """
        self.saveSettingsFilename = asksaveasfilename()
        self.settingsCollector()
        templateSaver(self.settings.settings, self.saveSettingsFilename)
        
    def settingsCollector(self):
        """ Collects all values for settings when settings have been inputted and
            user have pressed OK-button
        """
        def converter(value, varType):
            if varType == float:
                value = float(value)
            elif varType == int:
                value = int(value)
            return value
            
        for setting in self.entries:
            settingType = self.settings.settings[setting].type
            value = self.entries[setting].get()
            self.settings.settings[setting].value = converter(value, settingType)
        for setting in self.buttons:
            value = self.buttons[setting].get()
            self.settings.settings[setting].value = value
        self.experimentName = self.nameEntry.get()
    
class setting:
    """ A setting class containing a specific setting's name, unit, value, 
        type, and description.
    """
    def __init__(self, name, unit, value, settingType, description, classification):
        self.name = name
        self.unit = unit
        self.value = value
        self.type = settingType
        self.description = description
        self.classification = classification

class settingsClass:
    """ A class which contains a dictionary of setting objects and a function 
        to add new settings.
    """
    def __init__(self, name):
        self.settings = {}
    
    def add_setting(self, name, unit, value, settingType, description, classification):
        self.settings[name] = setting(name, unit, value, settingType, description, classification)
        
    def get_settings(self):
        settingsDict = {}
        for setting in self.settings:
            settingsDict[setting] = self.settings[setting].value
        return settingsDict

def addSettingsToGUI(settings, classification, window, i):
    for setting in sorted(list(settings.settings.keys())):
        if classification in settings.settings[setting].classification:
            i += 1
            #Obtain set value:
            setSetting = settings.settings[setting]
            #Add labels for each setting:
            window.add_label(setting, i, 0)
            
            #Add entry or checkbutton depending on setting type:
            if settings.settings[setting].type == bool:
                window.add_checkbutton(setting, '', setSetting.value, i, 1)
            else:
                window.add_entry(setting, setSetting.value, i, 1)
                
            #Add units for all settings:
            window.add_label(setSetting.unit, i, 2)
    return i
            
def designGUI(window, settings):
    """ A function that sets all the buttons and entries in place in the GUI
    """
    window.add_label('Peak calling options:', 0, 0, bold=True)
    i = addSettingsToGUI(settings, 'callingType', window, 0)
    i+=1
    window.add_separator(i,0)
    i+=1
    window.add_label('Analysis options:', i, 0, bold=True)
    i = addSettingsToGUI(settings, 'analysisType', window, i)
    
    #Add buttons:
    window.add_OK_button('OK', i+1, 1)
    window.add_destroy_button('Cancel', i+1, 2)
    window.add_browse_isolationList_button(0, 3)
    window.add_browse_dataFiles_button(2, 3)
    window.add_browse_settingFiles_button(4, 3)
    window.add_save_settings_button(5, 3)
    window.add_browse_saveDir_button(7,3)
    window.add_browse_experimentFiles_button(9,3)
    window.add_experimentNameField(12,3)
    window.add_help_button(15, 3)


def standardSettings():
    """ Returns standard settings dictionary
    """
    try:
        #Tries to access previously saved standard settings:
        standardSettings = templateLoader(os.path.join(os.getcwd(),'Setting templates\\standard settings.pkl'))
    except:
        #Creates a standard settings file if it is not found:
        standardSettingsList = [('Reduce noise', '', False, bool, 'Whether or not noise reduction should be done in each spectrum. This may increase analysis speed and may reduce FDR while also might increase false negatives.', 'callingType'),
                            ('MS1 quantification', '', False, bool, 'Whether or not to quantify on MS1 level. OBS: True value severely increases indexing time.', 'analysisType'),
                            ('MS2 quantification', '', True, bool, 'Whether or not quantification on MS2 level should be done.', 'analysisType'),
                            ('Scheduled run', '', True, bool, 'Whether or not an isolationList should be imported to limit data importation to schedule.', 'callingType'),
                            ('Fragment peptide calling', '', False, bool, 'Whether or not chromatogram peak calling based on peptide fragment presence should be done or not.', 'callingType'),
                            ('Heavy labels', '', 'K:0C2N, R:0C2N', str, 'The mass label composition. Is written in the form K:*C*N, R:*C*N. Currently only supports K and R labels.', 'callingType'),
                            ('MS1 tolerance', 'ppm', 10.0, float, 'The tolerance for calling precursor peaks at MS1 level. This defines the largest accepted deviation from the predicted MZ value for a peak to still be called.', 'callingType'),
                            ('MS2 tolerance', 'ppm', 7.0, float, 'The tolerance used in peak calling MS2 peaks. This defines the largest accepted deviation from the predicted MZ value for a peak to still be called.', 'callingType'),
                            ('Co-isolation', '', True, bool, 'A boolean depicting the nature of signal collecting. If heavy and light peaks are collected in concert. If heavy and light peaks are collected in alternating transitions this should be False.', 'callingType'),
                            ('Co-isolation (analysis)', '', True, bool, 'A boolean depicting the nature of signal collecting. If heavy and light peaks are collected in concert. If heavy and light peaks are collected in alternating transitions this should be False.', 'analysisType'),
                            ('IAA', '', True, bool, 'Boolean of whether iodoacetamide has been used or not. Peptide masses will be increased by 57 Da per Cysteine while True. Other modifications are not yet available.', 'callingType'),
                            ('Include b-fragments', '', False, bool, 'Boolean of wether to use b-fragments for quantification. This should only be activated when having asynchronous heavy and light collection.', 'callingType'),
                            ('Peaks in isotope clusters', '', 3, int, 'The maximum number of peaks in the isotope cluster that are to be searched for in data.', 'callingType'),
                            ('Isolation window size', 'Da', 1.2, float, 'The size of the isolation window.', 'callingType'),
                            ('Isolation window offset', 'Da', 0.55, float, 'The offset of the isolation window center compared to the monoisotopic peak.', 'callingType'),
                            ('Max fragment charge', '', 2, int, 'The maximum charge of MS2 ions to search for.', 'callingType'),
                            ('MS1 compensation', '', False, bool, 'Whether or not compensation should be done in MS1.', 'callingType'),
                            ('MS2 compensation', '', False, bool, 'Whether or not compensation should be done in MS2. Will only have effect if Co-isolation is activated.', 'callingType'),
                            ('Save figures', '', False, bool, 'Whether or not to save figures in analysis.', 'analysisType')
                            ]
        standardSettings = settingsClass('standard settings')
        for s in standardSettingsList:
            standardSettings.add_setting(s[0], s[1], s[2], s[3], s[4], s[5])
            
        daviPrint('Standard settings file could not be found.')
        templateSaver(standardSettings, os.path.join(os.getcwd(),'Setting templates\\standard settings'))
                            
    return standardSettings

def checkSettings(settings, isolationList, dataFiles, saveDirectory, experimentFiles):
    """ Function to check if inputted settings are correct.
    """
    if len(experimentFiles)>0:
        return True
    if isolationList == '':
        daviPrint('No isolation list was chosen.', line=True)
        return False
    if saveDirectory == '':
        daviPrint('No save directory was chosen.', line=True)
        return False
    if len(dataFiles) == 0:
        daviPrint('No data files were selected.', line = True)
        return False
    else: return True
    
def settingsGUI():
    # Load standard settings:
    settings = standardSettings()
    while True:
        #Create GUI window:
        window = GUI(settings)
        #Place entry fields, buttons and text:
        designGUI(window, settings)
        #Display interactive GUI:
        window.top.mainloop()
        #Get settings:
        settings = window.settings
        settingsValues = window.settings.get_settings()
        isolationList = window.isolationListFilename
        saveDirectory = window.saveDir
        dataFiles = window.dataFilenames
        experimentFiles = window.expFilenames
        experimentName = window.experimentName
        #Delete GUI for clean up:
        del window
        
        #Check if settings are correctly inputted:
        if checkSettings(settingsValues, isolationList, dataFiles, saveDirectory, experimentFiles):
            break
        
    return settingsValues, isolationList, dataFiles, saveDirectory, experimentFiles, experimentName



if __name__=="__main__":
    settingsGUI()
