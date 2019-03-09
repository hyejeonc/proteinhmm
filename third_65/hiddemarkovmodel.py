# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 17:25:48 2019

@author: HYEJEONG

Problem with float --> change to log 
"""
#from math import log
 
class Hmm(object):
    
    def __init__(self, prob_start, prob_trans, prob_emit, statelist, symbollist):

        self._prob_start = prob_start
        self._prob_trans = prob_trans
        self._prob_emit = prob_emit
        self._statelist = statelist
        self._symbollist = symbollist
           
        
    def check(self, sequence):
        if sequence not in self._symbollist:
            print('Not available sequence')
            
            
    def prob_start(self, state):
        if state not in self._statelist:
            print('Not available states, not in state list')
            return 0
        #print(state)
        #print(self._prob_start)    
        #print(self._prob_start[state])
        return self._prob_start[state]
 

    def prob_trans(self, pre_state, post_state):
        if pre_state not in self._statelist:
            print('Not available pre state, not in state list')
            return 0
        if post_state not in self._statelist:
            print('Not available post state, not in state list')
            return 0
            
        return self._prob_trans[pre_state][post_state]    
    
    
    def prob_emit(self, state, symbol):
        if state not in self._statelist:
            print('Not available states, not in state list')
            return 0
        if symbol not in self._symbollist:
            print('Not available symbols, not in symbol list')
            return 0
            
        return self._prob_emit[state][symbol]
       
                         
    def forward(self, sequence):
        
        a = [{}] # forward probability, alpha (list for [(t+1) states]; dict for probabilities of states {'h', 'e', '_'})            
                 # alpha is saved for next alpha, dynamics
        c = [{}] # scaling factor for very small number   
        
        for t in range(0, len(sequence)+1, +1): # state number : 0, 1, ... , (length)               
            if t == 0: 
                for state in self._statelist:                                    
                    a[0][state] = 0.0
            elif t == 1:
                for state in self._statelist:                 
                    a.append({})
                    a[1][state] = self.prob_start(state) * self.prob_emit(state, sequence[0])  #t = 0, first state. sequence must be list (or tuple)  함수냐 변수냐 그것이 문제                               
                
#            elif t == len(sequence)+1:
#                a[t][state] = 1.0
                
            else:
                for state in self._statelist:
                    a.append({})
                    sum_prob = 0.0
                    for pre_state in self._statelist:                     
                        sum_prob += a[t-1][pre_state] * self.prob_trans(pre_state, state)
                        #pre_state = state

                    #print('this is sum_prob ', sum_prob)
                    a[t][state] = sum_prob * self.prob_emit(state, sequence[t-1])

            if t != 0:
                print('this is sum a : ', sum(list(a[t].values())))
                
                c[t] = 1 / sum(list(a[t].values()))
                for state in self._statelist:
                    a[t][state] = c[t] * a[t][state]
            
            print('this is a[t][state] ', a[t][state])
        return a 
    
    
    ''' 
    output must be a list as below 
    
    a = [ {'h':... , 'e':... , '_':... },    a_{0} = first state * prob_emit = prob_start * prob_emit
                        ...    
          {'h':... , 'e':... , '_':... }, ]  a_{length} = last state * prob_emit
    '''    
    
    def probability(self, sequence):
        a = self._forward(sequence)
        sum_prob = 0.0
        
        for state in a[len(sequence)]:
            sum_prob += a[len(sequence)][state]
        
        return sum_prob
    
    
    def decode(self, sequence): #viterbi 
    
        v = [{}]
        dec_state = []
        dec_prob = []

        for t in range(0, len(sequence)+1): # state number : 0, 1, ... , (length)               
            
            v.append({})
            
            if t == 0: 
                for state in self._statelist:                                    
                    v[0][state] = 0.0
                    print('this is v[0] when t=0', v[0])
            elif t == 1: 
                print('this t : ', t)
                for state in self._statelist:                 

                    v[1][state] = self.prob_start(state) * self.prob_emit(state, sequence[0]) 

            else:                
                for state in self._statelist:       
                   # print('this is vmax after t = 2 :', vmax[1])
                    #print('this is prob : ', self.prob_trans)
                    #print('this t : ', t)
                    #print('this is length :', len(sequence))
                    v[t][state] = vmax[0] * self.prob_trans(vmax[1], state) * self.prob_emit(state, sequence[t-1])
                    #print('this t : ', t)
                    #print('this is length :', len(sequence))# pre_state = state
            vmax = list(max(zip(v[t].values(), v[t].keys())))
            
            if vmax[0] != 0:
                c = 1 / sum(list(v[t].values()))
                for state in self._statelist:
                   v[t][state] = c * v[t][state]
               
                vmax[0] = c * vmax[0]#decode.append(vmax)
                dec_prob.append(vmax[0])
                print('this is vmax after t>1 : ', vmax)            
                dec_state.append(vmax[1])
            
        return v[t], dec_state, vmax[0]
                   # 


    def backward(self, sequence):        
        b = [] # backward probability, beta (list for [(t+1) states]; dict for probabilities of states {'h', 'e', '_'})            
                  # beta is saved for next beta, dynamics programming
        c = [] # scaling factor for very small number   
        
        for t in range(len(sequence), 0, -1): # state number : 0, 1, ... , (length)               
            print(t)
            b.insert(0, {})
                    
            for state in self._statelist: 
                if t == len(sequence):
                    b[0][state] = 1.0        
                else:    
                    
                    sum_prob = 0.0
                    for post_state in self._statelist:                     
                        sum_prob += b[1][post_state] * self.prob_trans(state, post_state)
                  
                    b[0][state] = sum_prob * self.prob_emit(state, sequence[t])

            c = 1 / sum(list(b[0].values()))
            for state in self._statelist:
                b[0][state] = c * b[0][state]
            

        return b
''' 
    output must be a list as below 
    
    a = [ {'h':... , 'e':... , '_':... },    a_{0} = first state * prob_emit = _prob_start * prob_emit
                        ...    
          {'h':... , 'e':... , '_':... }, ]  a_{length} = last state * prob_emit
''' 