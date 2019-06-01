# -*- coding: utf-8 -*-
"""
Created on Thu Feb 21 14:44:36 2019

@author: HYEJEONG
"""
#import hmmexample as hmmex
import hiddemarkovmodel_solveprob_butlog as hmm
from sequence import *

#states = b # hidden states, secondary structure
#symbols = a # observable, amino acid sequences

statelist = ['h', 'e', '_']
symbollist = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

start_prob = prob_start
#print('prob start : ', prob_start)
trans_prob = prob_trans
#print(prob_trans)
emit_prob = prob_emit
#print(prob_emit)
model1 = hmmex.Model(statelist, symbollist, start_prob, trans_prob, emit_prob)
model2 = hmm.Hmm(start_prob, trans_prob, emit_prob, statelist, symbollist)

#singleresult = model2.decode(sequence[0])

result = []
for sequence, structure in zip(protein, secondstr):
#sequence = protein[0] 
#structure = secondstr[0] #['G', 'C', 'A', 'G', 'C', 'A', 'G', 'C', 'A', 'G', 'C', 'A']
    #print(sequence)
    #print(structure)
   # print('This is sequence in a interest : ', sequence)   
    model2.check(sequence)
    #print('This is decoded states : \n', model2.decode(sequence))
    result.append(model2.decode(sequence))
    # print('This is real states : \n', structure)
'''    
    
    print(model1.evaluate(sequence))
    print(model1.decode(sequence))
    
    model2 = hmm.Hmm(start_prob, trans_prob, emit_prob, statelist, symbollist)
    print('This is model \n', model2)
    
    print('This is evaluate prob \n', model2.probability(sequence))

    
    decoded = model2.decode(sequence)
    print('This is prob :' )
    
    print('This is real states : \n', structure)
    print(len(structure))
   
    
'''