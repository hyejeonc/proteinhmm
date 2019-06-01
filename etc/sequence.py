# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 16:49:57 2019

@author: HYEJEONG
"""


aminoAcid = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

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


allStr = [] 
with open('../dataset/protein-secondary-structure.train') as f: 
    for line in (line.strip() for line in f): 
     fields = line.split() 
     if fields: # non-blank line? 
         if readSeq(fields[0]):
             allStr += fields

proteinSingle = []
secstrSingle = []
protein = []
secstr = []
i = 0
for j in range(len(allStr)):
    if allStr[j] in aminoAcid:
#        print("########$$$$$$$$$$$$$$$$$$$$")            
        if allStr[j+2] in aminoAcid:
            proteinSingle.extend(allStr[j])
            secstrSingle.extend(allStr[j+1])  
        else:
            proteinSingle.extend(allStr[j])
            secstrSingle.extend(allStr[j+1])            
            protein.append(proteinSingle)
            secstr.append(secstrSingle)
            proteinSingle = []
            secstrSingle = []         
    else: pass



#should I define amino sequence that is two in-a
"""
for j in range(len(allStr)):
    if allStr[j] in aminoAcid:
        print("########$$$$$$$$$$$$$$$$$$$$")
        proteinSingle.extend(allStr[j])
        secstrSingle.extend(allStr[j+1])              
        if allStr[j+2] in aminoAcid:

        else:
            proteinSingle.extend(allStr[j+2])
            secstrSingle.extend(allStr[j+3])            
            protein.append(proteinSingle)
            secstr.append(secstrSingle)
            proteinSingle = []
            secstrSingle = []         
    else: pass
"""                 
'''            
        allStr[i+1] in aminoacid:
            protein[i][j] = allStr[i]
        elif 
        
    protein[i] = protein[j]
             
            ind = data.find(b'f i n a l   r e s u l t')
            data.seek(ind)

        data.seek(data.find(headstr.encode()))             
 '''
             
           
         
         
         