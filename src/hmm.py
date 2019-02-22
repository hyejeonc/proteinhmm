# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 15:25:50 2019

@author: HYEJEONG
"""

#홑밑줄(single underscore) : protected 보통 내부적으로 사용하는 변수일 때 사용합니다.
#곁밑줄(double underscore) : private 클래스 외부에서 접근할 수 없도록 내부 변수로 만듭니다.

'''
tranMatrix = [[0.6, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.4], 
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
transMatrix = [[0.2, 0.3, 0.5], 
               [0.3, 0.5, 0.2],
               [0.5, 0.2, 0.3]]

emitMatrix = [[0.4, 0.6, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
              [0.5, 0.4, 0.6, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
              [0.5, 0.5, 0.4, 0.6, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]]

aminoAcid = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
aminoAcidNum = []
            


def _normprob(prob, set):
    
    if :
else:
    sum = 0 
    for s in set:
        sum += prob.get()
    
else: 


class defineModel(object):
    
    def __init__(self, secStr, aminoAcid, trans = None, emiss = None): 
    # initializing 
    ''' 
        input
                secStr : 상태(states), helix/coil/sheet
                aminoAcid : 관측치(observation), types of amino acid (A~Z except B,J,O,U,X,Z)   
                trans : 전이확률
                
    '''
        self._secStr = set(secStr)
        self._aminoAcid = set(aminoAcid)
        self._trans = 
        self._emiss = 
        
    def _forward(self, )