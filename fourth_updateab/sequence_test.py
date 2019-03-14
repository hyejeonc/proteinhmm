# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 16:49:57 2019

@author: HYEJEONG
"""
from collections import Counter

#aminoacid = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
#structure = ['h', 'e', '_']


#symbollist = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'}
#statelist = {'h', 'e', '_'}

def _readfile(string):
        """ Determine if argument is a string representing a numeric value. """ 
        for kind in (str, str, str): 
            try: 
                kind(string) 
            except (TypeError, ValueError): 
                pass 
            else: 
                return True 
        else: 
            return False 

def count1d(count, item):
    if item not in count:
        count[item] = 0
    count[item] += 1

def count2d(count, item1, item2):
    if item1 not in count:
        count[item1] = {}
    count1d(count[item1], item2)

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

        
class Seq(object):
    
    def __init__(self, path, statelist, symbollist):
        self._path = path
        self._statelist = set(statelist)
        self._symbollist = set(symbollist)
        
    def path(self):
        return self._path

    def states(self):
   #     if statelist not in self._statelist
        return set(self._statelist)

    def symbols(self):
        return set(self._symbollist)    
                
    
   
    def lineseq(self): #lineseq(p, )
        allstring = [] 
        with open(self._path) as f: 
            for line in (line.strip() for line in f):  ## line 말고 통째로 split 
             fields = line.split() 
             if fields: # non-blank line? 
                 
                 if _readfile(fields[0]):
                     allstring += fields
    
        protein_single = []
        secondstr_single = []
        protein = []
        secondstr = []      
        for j in range(len(allstring)):  ## <> end 
            if allstring[j] in self._symbollist:
        #        print("########$$$$$$$$$$$$$$$$$$$$")            
                if allstring[j+2] in self._symbollist:
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
        
        return [protein, secondstr]
    
 
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
    #        count1d(symbolcount, symbol)        
    #        count1d(statecount, state)        
            count2d(statesymbolcount, state, symbol)
    prob_emit = _norm2d(statesymbolcount, statelist, self._symbollist)
    '''
    
    def prob(self, protein, secondstr):

        symbolcount = {}
        statesymbolcount = {}
        
        #sum_emitcount = {}
        prob_emit = {}
        
        sum_transcount = {}
        prob_trans = {}
        
        for structure, chain in zip(secondstr, protein):
            transcount = {}
            pre_statecount = {}
            
            startcount = {}
            prob_start = {}        
            count1d(startcount, structure[0])    
            
            for state, symbol in zip(structure, chain): #single protein
                statecount = {}
                count1d(statecount, state) 
                count2d(statesymbolcount, state, symbol) # for transmission probability           
                if pre_statecount != {}: 
                    count2d(transcount, pre_state, state) # for emittance probability            
                pre_state = state
                count1d(pre_statecount, pre_state)
        
            if sum_transcount == {}:
                sum_transcount = transcount 
            else:
                for key1 in list(sum_transcount.keys()):
                    for key2 in list(transcount.keys()):
                        if key1 == key2:                    
                            b = Counter(sum_transcount[key1]) + Counter(transcount[key2])
                            sum_transcount[key1] = dict(Counter(b))
        
        prob_trans = _norm2d(sum_transcount, self._statelist, self._statelist)
        prob_emit = _norm2d(statesymbolcount, self._statelist, self._symbollist)
        #print(prob_trans)        

        prob_start = _norm1d(startcount, self._statelist)
        #print(prob_start)
        return [prob_start, prob_trans, prob_emit]
