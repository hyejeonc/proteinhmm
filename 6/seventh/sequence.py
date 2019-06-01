# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 16:49:57 2019

@author: HYEJEONG
"""
from collections import Counter
import mathematics
#http://www.compbio.dundee.ac.uk/jpred4/about.shtml
#aminoacid = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
#structure = ['h', 'e', '_']


#symbollist = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'}
#statelist = {'h', 'e', '_'}

def readfile(string):
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


def lineseq(path): #lineseq(p, )
    allstring = [] 
    with open(path) as f: 
        for line in (line.strip() for line in f):  ## line 말고 통째로 split 
         fields = line.split() 
         if fields: # non-blank line? 
             if readfile(fields[0]):
                 allstring += fields
    return allstring 


def proteinseq(allstring, dbtype = 1):
    if dbtype == 0:
        None #FASTA file reading method will be here       
    else:
        protein = []
        secondstr = []   
        i = None
        for j in range(len(allstring)):    
            if allstring[j] == '<>' or allstring[j] == '>':
                protein_single = []
                secondstr_single = []
                i = j+1
                while i < len(allstring):
                    if allstring[i] == 'end' or allstring[i] == '<>' or allstring[i] == '<end>'  :
                        protein.append(protein_single)
                        secondstr.append(secondstr_single)
                        protein_single = []
                        secondstr_single = []
                        break
                    else: #allstring[i] != '<' or allstring[i] != '>' or allstring[i] != '<>':
                        protein_single.extend(allstring[i])
                        secondstr_single.extend(allstring[i+1])
                        i += 2

    return protein, secondstr         
                        
                    
                    
'''                    
                for i in range
                protein_single.extend(allstring[j+1])
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
prob_emit = norm2d(statesymbolcount, statelist, self._symbollist)
'''

def prob(protein, secondstr):

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
    
    prob_trans = norm2d(sum_transcount, self._statelist, self._statelist)
    prob_emit = norm2d(statesymbolcount, self._statelist, self._symbollist)
    #print(prob_trans)        

    prob_start = norm1d(startcount, self._statelist)
    #print(prob_start)
    return [prob_start, prob_trans, prob_emit]
