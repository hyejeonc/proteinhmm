# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 15:41:52 2019

@author: HYEJEONG
"""
from mathematics import *
import pandas as pd
    
def accuracy(origset, proteinset):
    """
    This method is for analysing accuracy depends on states by comparing two sets.
    Input : [ [Protein1, Protein2, Protein3 ... ], [Structure1, Structure2, Structure3 ... ] ]
    """
    proteins1 = origset[0]
    structures1 = origset[1]
    structures2 = proteinset[1]
    
    statecount = {}
    symbolcount = {}     
    crosscount = {}
    
    single_crosscount = {}    
    q3_single = []
    
    correctcount = 0
    totalcount = 0

    for protein1, structure1, structure2 in zip(proteins1, structures1, structures2):
        single_correctcount = 0
        single_crosscount = {}
        for symbol1, state1, state2 in zip(protein1, structure1, structure2):
            totalcount += 1
            count2d(crosscount, state1, state2)
            count2d(single_crosscount, state1, state2)
            count1d(statecount, state1)
            count1d(symbolcount, symbol1)
        '''
        q3h = 100 * single_crosscount['h']['h']/ sum(crosscount['h'].values())
        q3e = 100 * single_crosscount['e']['e']/ sum(crosscount['e'].values())
        q3c = 100 * single_crosscount['_']['_']/ sum(crosscount['_'].values())
        '''
        single_percent = 100 * single_correctcount / len(structure1)
        q3_single.append([len(structure1), single_correctcount, single_percent])
      #  q3_single_dic = single_crosscount] 
    
    statelist = statecount.keys()
    symbollist = symbolcount.keys()

    q3 = 100 * correctcount / totalcount
    '''
    q3h_tot = 100 * crosscount['h']['h']/ sum(crosscount['h'].values)
    q3e_tot = 100 * crosscount['e']['e']/ sum(crosscount['e'].values)
    q3c_tot = 100 * crosscount['_']['_']/ sum(crosscount['_'].values)
    '''
    q3_tot = [totalcount, correctcount, q3]
    q3_prob = norm2d(crosscount, statelist, statelist)
 
    return q3_tot, q3_prob, q3_single
                

def summary(proteinset):
    """
    This method is for analysing a amino acid sequence - secondary structure sequence data set.
    Input : [ [Protein1, Protein2, Protein3 ... ], [Structure1, Structure2, Structure3 ... ] ]

    """
    proteins = proteinset[0]
    structures = proteinset[1]
    
    statecount = {}
    symbolcount = {}         
    statesymbolcount = {}
  
    startcount = {}
    transcount = {} 
    for protein, structure in zip(proteins, structures):
        pre_state = None   
        count1d(startcount, structure[0])    
        for state, symbol in zip(structure, protein): #single protein
            count1d(statecount, state)
            count1d(symbolcount, symbol)
            count2d(statesymbolcount, state, symbol) # for transmission probability           
            if pre_state != None: 
                count2d(transcount, pre_state, state) # for emittance probability            
            pre_state = state
    
    statelist = list(statecount.keys())
    symbollist = list(symbolcount.keys()) 
    prob_start = norm1d(startcount, statelist)
    prob_trans = norm2d(transcount, statelist, statelist)
    prob_emit = norm2d(statesymbolcount, statelist, symbollist)
    '''
    print('this is emit count for simple', statesymbolcount)
    print('Frequency table of start states\n', pd.DataFrame(startcount, index=[0]))
    print('Frequency table of transition, pre>>post\n', pd.DataFrame(transcount))
    print('Frequency table of emission, states>>symbols\n', pd.DataFrame(statesymbolcount))
    
    
    print('Start probabilities\n', pd.DataFrame(prob_start, index=[0]))
    print('Transition probabilities, pre>>post\n', pd.DataFrame(prob_trans))
    print('Emission probabilities, states>>symbols\n', pd.DataFrame(prob_emit))
    '''
    return [[startcount, transcount, statesymbolcount], [prob_start, prob_trans, prob_emit], [statelist, symbollist]]
            
