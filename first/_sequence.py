# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 16:49:57 2019

@author: HYEJEONG
"""
from collections import Counter

#aminoacid = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
#structure = ['h', 'e', '_']


symbollist = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'}
statelist = {'h', 'e', '_'}



def readSeq(s):
    """ Determine if argument is a string representing a numeric value. """ 
    for kind in (str, str, str): 
     try: 
      kind(s) 
     except (TypeError, ValueError): 
      pass 
     else: 
      return True 
    else: 
     return False 

#def lineSeq()
allstring = [] 
with open('dataset/protein-secondary-structure.train') as f: 
    for line in (line.strip() for line in f):  ## line 말고 통째로 split 
     fields = line.split() 
     if fields: # non-blank line? 
         if readSeq(fields[0]):
             allstring += fields

protein_single = []
secondstr_single = []
protein = []
secondstr = []
i = 0
for j in range(len(allstring)):  ## <> end 
    if allstring[j] in symbollist:
#        print("########$$$$$$$$$$$$$$$$$$$$")            
        if allstring[j+2] in symbollist:
            protein_single.extend(allstring[j])
            secondstr_single.extend(allstring[j+1])  
        else:
            protein_single.extend(allstring[j])
            secondstr_single.extend(allstring[j+1])            
            protein.append(protein_single)
            secondstr.append(secondstr_single)
            protein_single = []
            secondstr_single = []         
    else: pass

def _count1d(count, item):
    if item not in count:
        count[item] = 0
    count[item] += 1

def _count2d(count, item1, item2):
    if item1 not in count:
        count[item1] = {}
    _count1d(count[item1], item2)

def _norm1d(prob, item_set):  
    result = {} 
    prob_sum = 0.0
    
    for item in item_set:
        prob_sum += prob.get(item, 0)
    
    if prob_sum == 0:
        result[item] = 0
    else:    
        for item in item_set:
            result[item] = prob.get(item, 0) / prob_sum            
    return result
                       
def _norm2d(prob, item_set1, item_set2):
    result = {}
    
    if prob is None:
        for item in item_set1:
            result[item] = _norm1d(None, item_set2)
    for item in item_set1:
        result[item] = _norm1d( prob.get(item, 0), item_set2)
    
    return result    
'''
sum = 0.0
k = 0
l = 0

statecount = {}
symbolcount = {}
statesymbolcount = {}
prob_emit = {}
      
for structure, chain in zip(secondstr, protein):
    for state, symbol in zip(structure, chain):
#        _count1d(symbolcount, symbol)        
#        _count1d(statecount, state)        
        _count2d(statesymbolcount, state, symbol)
prob_emit = _norm2d(statesymbolcount, statelist, symbollist)
'''


statecount = {}
symbolcount = {}
statesymbolcount = {}

sum_emitcount = {}
prob_emit = {}

sum_transcount = {}
prob_trans = {}

for structure, chain in zip(secondstr, protein):
    transcount = {}
    pre_statecount = {}
    
    for state, symbol in zip(structure, chain):
        statecount = {}
        _count1d(statecount, state) 
        _count2d(statesymbolcount, state, symbol) # for transmission probability           
        if pre_statecount != {}: 
            _count2d(transcount, pre_state, state) # for emittance probability            
        pre_state = state
        _count1d(pre_statecount, pre_state)

    if sum_transcount == {}:
        sum_transcount = transcount 
    else:
        for key1 in list(sum_transcount.keys()):
            for key2 in list(transcount.keys()):
                if key1 == key2:                    
                    b = Counter(sum_transcount[key1]) + Counter(transcount[key2])
                    sum_transcount[key1] = dict(Counter(b))

prob_trans = _norm2d(sum_transcount, statelist, statelist)
prob_emit = _norm2d(statesymbolcount, statelist, symbollist)
print(prob_trans)

startcount = {}
prob_start = {}
        

for structure, chain in zip(secondstr, protein):
    _count1d(startcount, structure[0])
prob_start = _norm1d(startcount, statelist)
print(prob_start)
