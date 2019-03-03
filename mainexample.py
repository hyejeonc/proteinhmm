# -*- coding: utf-8 -*-
"""
Created on Thu Feb 21 14:44:36 2019

@author: HYEJEONG
"""
import hmmexample as hmmex

from sequence import *

a = protein[2]
b = secondstr[2]
#states = b # hidden states, secondary structure
#symbols = a # observable, amino acid sequences

states = ('h', 'e', '_')
symbols = ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')

start_prob = {
    'h' : 0.33,
    'e' : 0.33,
    '_' : 0.34
}
#aminoAcid = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
'''
tranMatrix = [[0.06, 0.05, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.4], 
              [0.4, 0.6, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
              [0.5, 0.4, 0.6, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
              [0.5, 0.6, 0.4, 0.6, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
              [0.5, 0.6, 0.5, 0.4, 0.6, 0.5, 0.5, 0.5, 0.5, 0.5],
              [0.5, 0.6, 0.5, 0.5, 0.4, 0.6, 0.5, 0.5, 0.5, 0.5],
              [0.5, 0.6, 0.5, 0.5, 0.5, 0.4, 0.6, 0.5, 0.5, 0.5],
              [0.5, 0.6, 0.5, 0.5, 0.5, 0.5, 0.4, 0.6, 0.5, 0.5],
              [0.5, 0.6, 0.5, 0.5, 0.5, 0.5, 0.5, 0.4, 0.6, 0.5],
              [0.5, 0.6, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.4, 0.6]]
'''

trans_prob = prob_trans
'''
{
        'h' : { 'h' : 0.3, 'e' : 0.2, '_' : 0.5 },
        'e' : { 'h' : 0.1, 'e' : 0.3, '_' : 0.6 },
        '_' : { 'h' : 0.2, 'e' : 0.3, '_' : 0.5 }
}
'''
emit_prob = prob_emit
'''
{
    'h': { 'A' : 0.09, 'C' : 0.01, 'D' : 0.05, 'E' : 0.05, 'F' : 0.01, 'G' : 0.09, 'H' : 0.05, 'I' : 0.05, 'K' : 0.05, 'L' : 0.05, 'M' : 0.05, 'N' : 0.05, 'P' : 0.05, 'Q' : 0.05, 'R' : 0.05, 'S' : 0.05, 'T' : 0.05, 'V' : 0.05, 'W' : 0.05, 'Y' : 0.05 },
    'e': { 'A' : 0.05, 'C' : 0.06, 'D' : 0.09, 'E' : 0.01, 'F' : 0.05, 'G' : 0.05, 'H' : 0.05, 'I' : 0.05, 'K' : 0.05, 'L' : 0.05, 'M' : 0.05, 'N' : 0.05, 'P' : 0.05, 'Q' : 0.05, 'R' : 0.05, 'S' : 0.05, 'T' : 0.05, 'V' : 0.05, 'W' : 0.05, 'Y' : 0.05 },
    '_': { 'A' : 0.05, 'C' : 0.05, 'D' : 0.06, 'E' : 0.01, 'F' : 0.09, 'G' : 0.05, 'H' : 0.05, 'I' : 0.05, 'K' : 0.05, 'L' : 0.05, 'M' : 0.05, 'N' : 0.05, 'P' : 0.05, 'Q' : 0.05, 'R' : 0.05, 'S' : 0.05, 'T' : 0.05, 'V' : 0.05, 'W' : 0.05, 'Y' : 0.05 }
}
'''
sequence = a #['G', 'C', 'A', 'G', 'C', 'A', 'G', 'C', 'A', 'G', 'C', 'A']
print('This is sequence in a interest : ', sequence)
model = hmmex.Model(states, symbols, start_prob, trans_prob, emit_prob)
print('This is model \n', model)
print('This is evaluate prob \n', model.evaluate(sequence))



print('This is decoded states : \n', model.decode(sequence))
print(b)