#! /usr/bin/env python[3.6]
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 10:00:43 2019

@author: HYEJEONG
"""

#exe 만들기 리눅스 윈도우 http://codingdojang.com/scode/272?orderby=time
# pyinstaller https://hyrama.com/?p=579
import statistics as stat
import sequence as seq
import os.path 
import sys

path_train = None
rawlines = None

rawlines = seq.lineseq('.\dataset\protein-secondary-structure.train')
a = seq.proteinseq(rawlines)
stat.frequency(seq.proteinseq(rawlines)[0], seq.proteinseq(rawlines)[1])

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
