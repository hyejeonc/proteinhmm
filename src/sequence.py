# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 16:49:57 2019

@author: HYEJEONG
"""
"""
f = open("protein-secondary-structure.train", 'r')
all = f.read()

while True:
    line = f.readlines()
    if not line: break
    print(line)
f.close()

lines = all.split('\n')
"""

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
with open('protein-secondary-structure.train') as f: 
    for line in (line.strip() for line in f): 
     fields = line.split() 
     if fields: # non-blank line? 
         if readSeq(fields[0]):
             allStr += fields
             
           
         
         
         