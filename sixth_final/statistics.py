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
'''
def probplot(prob_start, prob_trans, prob_emit): 
    

def heatplot(data1, date2):
     
    flights = sns.load_dataset("flights")
    flights = flights.pivot("month", "year", "passengers")
    plt.figure(figsize=(10, 10))
    ax = sns.heatmap(flights, annot=True, fmt="d")
'''
def singleprobplot(orig_protein, orig_str, protein2, secondstr2, protein3, secondstr3):
    first = accuracy(orig_protein, orig_str, protein2, secondstr2)
    second = accuracy(orig_protein, orig_str, protein3, secondstr3)
    # [[correct in single, total in single, percent], [], ... , []]
    
  
    print('First(test1) Q3 by protein length\n', pd.DataFrame(first))
    print(first)
    
    print('Second(test2) Q3 by protein length\n', pd.DataFrame(second))
    print(second)
    
    w = 0.3
    plot.bar(first[:][0]-w, first[:][2], )
    plot.bar(second[:][0]+w, second[:][2], )
    plt.show()
    
    return first, second
    
def accuracy(protein1, secondstr1, protein2, secondstr2):
#return 
    singlelist = []
    correctcount = 0
    totalcount = 0
    crosscount = {}
    
    for chain1, structure1, structure2 in zip(protein1, secondstr1, secondstr2):
        single_correctcount = 0
        for symbol1, state1, state2 in zip(chain1, structure1, structure2):
            totalcount += 1
            count2d(crosscount, state1, state2)
            count1d(statecount, state1)
            count1d(symbolcount, symbol1)
            if state1 == state2 :
                correctcount += 1
                single_correctcount += 1
  
        single_percent = 100 * single_correctcount / len(structure1)
        singlelist.append([single_correctcount, len(structure1), single_percent]) #protein 마다 q3 를 구함 
    
    statelist = list(statecount.keys())
    symbollist = list(symbolcount.keys()) 

    q3 = 100 * correctcount / totalcount    
    print('Correct, Total, Q3\n', correctcount, totalcount, q3)
    
    prob_cross = norm2d(crosscount, statelist, statelist)

    print('Frequency table of states for Q3a, Q3b, Q3c. original >> prediction\n', pd.DataFrame(crosscount)) 
    print('Probabilities for Q3a, Q3b, Q3c. original >> prediction\n', pd.DataFrame(norm2d(crosscount, statelist, statelist))) 

    
    return q3, prob_cross, singlelist
                

def summary(protein, secondstr):
        statecount = {}
        symbolcount = {}         
        statesymbolcount = {}
      
        startcount = {}
        transcount = {} 
        for chain, structure in zip(protein, secondstr):
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

        print('Frequency table of start states\n', pd.DataFrame(startcount, index=[0]))
        print('Frequency table of transition, pre>>post\n', pd.DataFrame(transcount))
        print('Frequency table of emission, states>>symbols\n', pd.DataFrame(statesymbolcount))
        
        prob_start = norm1d(startcount, statelist)
        prob_trans = norm2d(transcount, statelist, statelist)
        prob_emit = norm2d(statesymbolcount, statelist, symbollist)
        print('Start probabilities\n', pd.DataFrame(prob_start, index=[0]))
        print('Transition probabilities, pre>>post\n', pd.DataFrame(prob_trans))
        print('Emission probabilities, states>>symbols\n', pd.DataFrame(prob_emit))

        return prob_start, prob_trans, prob_emit, statelist, symbollist
            
