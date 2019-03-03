# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 16:49:57 2019

@author: HYEJEONG
"""



aminoAcid = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
structure = ['h', 'e', '_']

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

def readPath(p):
    allStr = [] 
    with open('p') as f: 
        for line in (line.strip() for line in f): 
         fields = line.split() 
         if fields: # non-blank line? 
             if readSeq(fields[0]):
                 allStr += fields
        
    proteinSingle = []
    secStrSingle = []
    protein = []
    secStr = []

    for j in range(len(allStr)):
        if allStr[j] in aminoAcid:
    #        print("########$$$$$$$$$$$$$$$$$$$$")            
            if allStr[j+2] in aminoAcid:
                proteinSingle.extend(allStr[j])
                secStrSingle.extend(allStr[j+1])  
            else:
                proteinSingle.extend(allStr[j])
                secStrSingle.extend(allStr[j+1])            
                protein.append(proteinSingle)
                secStr.append(secStrSingle)
                proteinSingle = []
                secStrSingle = []         
        else: pass

def readEmitProb(proteinList, secStrList):
    
    prob = {}
    sum = 0.0
    
    for state in structure:
        prob[state] = 0
        for i in range(len(proteinList[i])):
            for j in range(len(proteinList[i][j])):        
                if protein[i][j] == symbol :
                    prob[state][symbol] += 1                
#                else : 
        print(prob)                
        sum += prob[symbol]   
    
    prob[symbol] = prob[symbol] 
    
    
def readEmitP:
        




p = dataset/protein-secondary-structure.train


'''

trans_prob = {
        'h' : { 'h' : 0.3, 'e' : 0.2, '_' : 0.5 },
        'e' : { 'h' : 0.1, 'e' : 0.3, '_' : 0.6 },
        '_' : { 'h' : 0.2, 'e' : 0.3, '_' : 0.5 }
}

emit_prob = {
    'h': { 'A' : 0.09, 'C' : 0.01, 'D' : 0.05, 'E' : 0.05, 'F' : 0.01, 'G' : 0.09, 'H' : 0.05, 'I' : 0.05, 'K' : 0.05, 'L' : 0.05, 'M' : 0.05, 'N' : 0.05, 'P' : 0.05, 'Q' : 0.05, 'R' : 0.05, 'S' : 0.05, 'T' : 0.05, 'V' : 0.05, 'W' : 0.05, 'Y' : 0.05 },
    'e': { 'A' : 0.05, 'C' : 0.06, 'D' : 0.09, 'E' : 0.01, 'F' : 0.05, 'G' : 0.05, 'H' : 0.05, 'I' : 0.05, 'K' : 0.05, 'L' : 0.05, 'M' : 0.05, 'N' : 0.05, 'P' : 0.05, 'Q' : 0.05, 'R' : 0.05, 'S' : 0.05, 'T' : 0.05, 'V' : 0.05, 'W' : 0.05, 'Y' : 0.05 },
    '_': { 'A' : 0.05, 'C' : 0.05, 'D' : 0.06, 'E' : 0.01, 'F' : 0.09, 'G' : 0.05, 'H' : 0.05, 'I' : 0.05, 'K' : 0.05, 'L' : 0.05, 'M' : 0.05, 'N' : 0.05, 'P' : 0.05, 'Q' : 0.05, 'R' : 0.05, 'S' : 0.05, 'T' : 0.05, 'V' : 0.05, 'W' : 0.05, 'Y' : 0.05 }
}
'''
 