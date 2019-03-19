# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 15:41:52 2019

@author: HYEJEONG
"""
from mathematics import *
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def frequency(protein, secondstr):
        statecount = {}
        symbolcount = {}         
        statesymbolcount = {}
      
        startcount = {}
        transcount = {} 
        for structure, chain in zip(secondstr, protein):
            pre_state = None   
            count1d(startcount, structure[0])    
            for state, symbol in zip(structure, chain): #single protein
                count1d(statecount, state)
                count1d(symbolcount, symbol)
                count2d(statesymbolcount, state, symbol) # for transmission probability           
                if pre_state != None: 
                    count2d(transcount, pre_state, state) # for emittance probability            
                pre_state = state
        
        statelist = list(statecount.keys())
        symbollist = list(symbolcount.keys()) 
        
        print(startcount)
        print(transcount)
        print(statesymbolcount)
        

        freq_start = pd.DataFrame(startcount, index=[0])
        freq_trans = pd.DataFrame(transcount)
        freq_emit = pd.DataFrame(statesymbolcount)
        print(freq_start)
        print(freq_trans)
        print(freq_emit)
        
        prob_start = pd.DataFrame(norm1d(startcount, statelist), index=[0])
        prob_trans = pd.DataFrame(norm2d(transcount, statelist, statelist))
        prob_emit = pd.DataFrame(norm2d(statesymbolcount, statelist, symbollist))
        print(prob_start)
        print(prob_trans)
        print(prob_emit)
            
        '''
        flights = sns.load_dataset("flights")
        flights = flights.pivot("month", "year", "passengers")
        plt.figure(figsize=(10, 10))
        ax = sns.heatmap(flights, annot=True, fmt="d")
        '''