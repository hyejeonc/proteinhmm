# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 16:49:57 2019

@author: HYEJEONG
"""
from collections import Counter

aminoacid = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
structure = ['h', 'e', '_']


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
    if allstring[j] in aminoacid:
#        print("########$$$$$$$$$$$$$$$$$$$$")            
        if allstring[j+2] in aminoacid:
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

def emit_prob(self, state, symbol):
    """
    Return the emission probability for `symbol` associated with the `state`.

    If either the `state` or the `symbol` are not contained in this model,
    0 is returned.
    """
    if state not in self._states or symbol not in self._symbols:
        return 0
    return self._emit_prob[state][symbol]
  '''
    

sum = 0.0
k = 0
l = 0
statesymbolcount = {}

for structure, chain in zip(secondstr, protein):
    for state, symbol in zip(structure, chain):
#        _count1d(symbolcount, symbol)        
#        _count1d(statecount, state)        
        _count2d(statesymbolcount, state, symbol)
prob_emit = _norm2d(statesymbolcount, statelist, symbollist)


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
        if pre_statecount != {}: _count2d(transcount, pre_state, state)        #important            
        pre_state = state
        _count1d(pre_statecount, pre_state)

    if sum_transcount == {}:
        print('**all sum and,. this is sum_transcount : ', sum_transcount)
        sum_transcount = transcount 
    else:
        print('what the hell is the sum_transcount? : ' , sum_transcount)
        print('what the hell is the transcount? : ' , transcount)
        for key1 in list(sum_transcount.keys()):
            for key2 in list(transcount.keys()):
                if key1 == key2:                    
                    b = Counter(sum_transcount[key1]) + Counter(transcount[key2])
                    sum_transcount[key1] = dict(Counter(b))

prob_trans = _norm2d(sum_transcount, statelist, statelist)
print(prob_trans)

'''           
if sum_emitcount == {}:
    sum_emitcount = statesymbolcount 
else:
    a = Counter(sum_emitcount)
#    print(a)
    
    b = Counter(statesymbolcount)
#    print(b)
    sum_emitcount = dict(Counter(a + b))
#    print(sum_emitcount)

            
            
            
        for key, value in sum_transcount.items():
            for k, v in transcount.items():
                if key == k:
                    a = value.keys
                    b = v.keys
                    print(a)
                    print(b)
                #    a.keys = a.keys + b.keys
                #    = value + v
        
        z = sum_transcount.copy()
        z.update(transcount)

        
#        z = {**sum_transcount, **transcount}    
        print('is this sum of sum+trans? : ', z)
        for key, value in transcount.item():
            print(key)
            print(value)
            for k, v in value.item():
                print(k)
                print(v)


        for k, v in sum_transcount{}, sum_transcount{}{}:
            sum_transcount[k][v] += transcount[k][v] 
        print('!!!!the hell is the sum_transcount? : ' ,sum_transcount)    

        a = Counter({item[0] : item[2] for item in sum_transcount})
        b = Counter({item[0] : item[2] for item in transcount})
        
        sum_transcount = [[key, value] for key, value in (a + b).items()]
        
        print(' sum_transcount? : ' ,sum_transcount)
       
        aa = Counter(a)   
        print(a)
        print(aa)
        b = Counter(transcount)
        bb = Counter(b)
        print(a)
        print(b) # Counter({'_': {'_': 62, 'e': 8, 'h': 2}, 'e': {'e': 35, '_': 8}, 'h': {'h': 11, '_': 2}})
        c = aa + bb
        print(c)
        
        
        
        from collections import Counter
list1= [['some',2],['other',1],['thing',5]]

list2= [['some',1],['thing',5]]
c1 = Counter({item[0]: item[1] for item in list1})
c2 = Counter({item[0]: item[1] for item in list2})
result = [[key, value] for key, value in (c1-c2).items()]

result
[['other', 1], ['some', 1]]
'''

'''       
#    print(b)
        sum_transcount = dict(Counter(a + b))
#    print(sum_emitcount)

#print(_norm2d(sum_emitcount, statelist, symbollist))
prob_trans = _norm2d(sum_transcount, statelist, statelist)


#        print('this is statecount', statecount)
#        print('this is pre statecount', pre_statecount)
'''
'''    
for structure, chain in zip(secondstr, protein):
    for state, symbol in zip(structure, chain):
        _count1d(symbolcount, symbol)        
        _count1d(statecount, state)        
        _count2d(statesymbolcount, state, symbol)
        
if sum_emitcount == {}:
    sum_emitcount = statesymbolcount 
else:
    a = Counter(sum_emitcount)
#    print(a)
    
    b = Counter(statesymbolcount)
#    print(b)
    sum_emitcount = dict(Counter(a + b))
#    print(sum_emitcount)

#print(_norm2d(sum_emitcount, statelist, symbollist))
prob_emit = _norm2d(sum_emitcount, statelist, symbollist)
'''



