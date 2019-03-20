#! /usr/bin/env python[3.6]
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 10:00:43 2019

@author: HYEJEONG
"""

#exe 만들기 리눅스 윈도우 http://codingdojang.com/scode/272?orderby=time
# pyinstaller https://hyrama.com/?p=579
import sequence as seq

import statistics as stat

import hiddenmarkovmodel as hmm
import os.path 
import sys

path_train = None
rawlines = None

rawlines = seq.lineseq('.\dataset\protein-secondary-structure.train')
testlines = seq.lineseq('.\dataset\protein-secondary-structure.test')
#proteinseq(rawlines)

trainset = seq.proteinseq(rawlines) 
testset = seq.proteinseq(testlines)

initialset = stat.summary(trainset[0], trainset[1])
#print statistical information about data sets, ex. frequency table, probabilities ... 
print(initialset)
old = hmm.initialize(initialset)
new = hmm.train(initialset, trainset[0], trainset[1], 0.001, 100)
print(new.get())

decode = []
likelihood = []
decode_em = []
likelihood_em = []

for i in range(len(testset[0])):
    for chain1, structure1 in zip(testset[0][i], testset[1][i]):
        viterbi = old.viterbi(chain1)
        decode.append(viterbi[:][1])
        likelihood.append(viterbi[:][2])
    
        viterbi_em = new.viterbi(chain1)
        decode_em.append(viterbi_em[:][1])
        likelihood_em.append(viterbi_em[:][2])

oldset = [testset[0], decode]
newset = [testset[0], decode_em] 

stat.summary(oldset[0], oldset[1])
stat.summary(newset[0], newset[1])

stat.singleprobplot(testset[0], testset[1], oldset[0], oldset[1], newset[0], newset[1])






'''
while True:
    try:
        path_train = input('Insert input file for training \n(for example: .\dataset\protein-secondary-structure.train) : ')
        rawlines = seq.lineseq(path_train)
    except FileNotFoundError:
        print('No available input file. Select options with integer number below (1, 2, 3)', \
                        '1. Try another file', \
                        '2. Use example input file (.\dataset\protein-secondary-structure.train)', \
                        '3. Exit', \
                        sep='\n')
        trytype = input()
        if trytype == 2:
            path_train = '.\dataset\protein-secondary-structure.train'
            break
        elif trytype == 3:
            sys.exit(0)
    if rawlines != None:
        break
        
 #
'''

#stat.frequency()
'''
path_train = Path(path)

if path_train == None:
    path_train = '.\dataset\protein-secondary-structure.test'
    
a = seq.lineseq(path_train) #
print(a)
#simple prob & training prob 출력 


path_test = input('Insert file directory for testing \n \
             For example : .\dataset\protein-secondary-structure.test \n \
             : ')
'''
