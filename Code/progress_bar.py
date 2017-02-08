# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 11:57:00 2017

@author: max.karlsson
"""
import tkinter as tk
from tkinter import ttk
import time
from statistics import median

class statusWindow(tk.Tk):
    def __init__(self, fileNumber):
        tk.Tk.__init__(self)
        self.progress =  ttk.Progressbar(self, orient="horizontal", length=200, mode="determinate")
        self.progress.grid(column = 0, row = 3)
        self.val = 0
        self.maxval = fileNumber
        self.progress["maximum"] = fileNumber
        self.fileLabel()
        self.timeLabel()
        self.times = []
        
    def updating(self, val, file):
        self.val = val
        self.progress["value"] = self.val
        self.labelText.set('Current file: ' + file.split('/')[-1] + '\n(' + str(val) + '/' + str(self.maxval) + ')')
        self.timeUpdate()
        if self.val == self.maxval:
            self.destroy()
    
    def timeUpdate(self):
        try:    
            self.times.append(time.time()-self.time)
            secTime = median(self.times)*(self.maxval-self.val)
            showTime = ' '.join([t for t in time.strftime('%Hh %Mm %Ss', time.gmtime(secTime)).split(' ') if '00' not in t])
            self.timeLeft.set(showTime)
        except:
            pass
        
        self.time = time.time()
        
    def fileLabel(self):
        self.labelText = tk.StringVar('')
        label = tk.Label(self, textvariable=self.labelText)
        label.grid(row=0, column=0)
        
    def timeLabel(self):
        self.timeLeft = tk.StringVar('')
        label = tk.Label(self, textvariable=self.timeLeft)
        label.grid(row=2, column=0)
        

        
def __fileLoop(index=0):
    """ This module is for development only.
    """
    time.sleep(5)
    window.updating(index+1, files[index])
    if index+1<len(files):
        index += 1
        window.after(1, __fileLoop, index)

    
if __name__=='__main__':
    files = 10*('C:/Users/max.karlsson/Documents/Max/Wellness project/wellness_prm(pilot)_vis(2)_pla(1)_160522_A1_116157.mzML', 'C:/Users/max.karlsson/Documents/Max/Wellness project/wellness_prm(pilot)_vis(2)_pla(1)_160522_A2_117055.mzML')
    window = statusWindow(len(files))
    
    window.after(1, __fileLoop, 0)
    window.mainloop()
    
    
    