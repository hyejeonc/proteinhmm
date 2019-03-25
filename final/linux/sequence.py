#!/usr/bin/env python[3.6]
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 16:49:57 2019

@author: HYEJEONG
"""
from collections import Counter
from mathematics import *

#aminoacid = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
#structure = ['h', 'e', '_']


#symbollist = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'}
#statelist = {'h', 'e', '_'}

def readfile(string):
        """
        This is a method to determine if argument is a string representing a numeric value. 
        """ 
        for kind in (str, str, str): 
            try: 
                kind(string) 
            except (TypeError, ValueError): 
                pass 
            else: 
                return True 
        else: 
            return False 

def lineseq(path): 
    """
    This method is for seperating string to word for getting symbols and states.
    """
    allstring = [] 
    with open(path) as f: 
        for line in (line.strip() for line in f):  ## line 말고 통째로 split 
         fields = line.split() 
         if fields: # non-blank line? 
             if readfile(fields[0]):
                 allstring += fields
    return allstring 

def proteinseq(path, dbtype = 1):
    """
    This is a method for getting a protein set. 
    output = [ [ [ protein1 ], [ protein2 ], [ protein3 ], [ protein4 ], ... ],
               [ [ structure1 ], [ structure2 ], [ structure3 ], [ structure4 ], ... ] ]
    """   
    if dbtype == 0:
        None #FASTA file reading method will be here... To be continued...      
    else:
        allstring = lineseq(path)
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

def getproteinset(path):
    proteinset = proteinseq(path) 
    return proteinset                        
                    
