# -*- coding: utf-8 -*-
"""
Created on Sun Dec 11 12:05:48 2016

@author: Max Karlsson
"""

def daviPrint(text, pause=False, line=False):
    """ Basic printer function with functionality to pause script and add a 
        separator line.
    """
    if line: 
        print('----------------------------------------------------\n')
    
    print(text)
    
    if pause: 
        input('Press ENTER to continue.')