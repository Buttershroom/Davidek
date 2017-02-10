# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 11:57:00 2017

@author: max.karlsson
"""
import tkinter as tk
from tkinter import ttk
import time
from statistics import median

class bar(tk.Tk):
    def __init__(self, maxval):
        tk.Tk.__init__(self)
        #Create progress bar:
        self.progressBar =  ttk.Progressbar(self, orient="horizontal", length=200, mode="determinate")
        self.progressBar.grid(column = 0, row = 3)
        
        #Set progress bar values:
        self.progressBar["maximum"] = maxval
        self.progressBar["value"] = 0
        
        #Create text in window:
        self.label()
    
    def label(self):
        self.text = tk.StringVar('')
        label = tk.Label(self, textvariable=self.text)
        label.grid(row=0, column=0)
        
    def add(self):
        self.progressBar['value'] += 1
        if self.progressBar['value'] == self.progressBar["maximum"]:
            self.destroy()
            
    def currentText(self, text):
        self.text.set(text)
        
def fileLoop(index, fileTime):
    """ For development and testing purposes:
    """
    time.sleep(5)
    showTime = ''
    #Do stuff:
    if not fileTime == 0:
        times.append(time.time()-fileTime)
        secTime = median(times)*(len(files)-index)
        showTime = ' '.join([t for t in time.strftime('%Hh %Mm %Ss', time.gmtime(secTime)).split(' ') if '00' not in t])

    fileTime = time.time()
    file = files[index]
    window.currentText('Current file: ' + file.split('/')[-1] + '\n(' + str(index+1) + '/' + str(len(files)) + ')\n' + showTime)
    #####
    window.add()
    index+=1
    if not index == len(files):
        window.after(1, fileLoop, index, fileTime)

if __name__=='__main__':
    fileTime = 0
    times = []
    index = 0
    files = 10*('C:/Users/max.karlsson/Documents/Max/Wellness project/wellness_prm(pilot)_vis(2)_pla(1)_160522_A1_116157.mzML', 'C:/Users/max.karlsson/Documents/Max/Wellness project/wellness_prm(pilot)_vis(2)_pla(1)_160522_A2_117055.mzML')
    window = bar(20)
    window.wm_title("Analysis status")
    window.after(10, fileLoop, index, fileTime)
    window.mainloop()

   