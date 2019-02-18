# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 15:25:50 2019

@author: HYEJEONG
"""

class defineModel(object):
    
    def __init__(self, secStr, aminoAcid, trans = None, emiss = None): 
    # initializing 
    ''' 
        input
                secStr : 상태(states), helix/coil/sheet
                aminoAcid : 관측치(observation), types of amino acid (A~Z except B,J,O,U,X,Z)   
                trans : 전이확률
                
    '''
        self._secStr = set(secStr)
        self._aminoAcid = set(aminoAcid)
        self._trans = 
        self._emiss = 
        
    def _forward(self, )