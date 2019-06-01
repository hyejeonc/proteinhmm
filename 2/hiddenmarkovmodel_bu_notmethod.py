# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 17:25:48 2019

@author: HYEJEONG
"""

#from math import log
 
class Hmm(object):
    
    def __init__(self, prob_start, prob_trans, prob_emit, statelist, symbollist):
        #
        # for 
        #
        #self._state = state.values()
        #self._symbol = symbol.values()
        self._prob_start = prob_start
        self._prob_trans = prob_trans
        self._prob_emit = prob_emit
        self._statelist = statelist
        self._symbollist = symbollist
        
    def check(self, sequence):
        if sequence not in self._statelist:
            print('Not available sequence')
                            
    def _forward(self, sequence):
        a = [{}] # forward probability, alpha (list for [(t+1) states]; dict for probabilities of states {'h', 'e', '_'})            
                 # alpha is saved for next alpha, dynamics 
        for t in range(0, len(sequence), +1): # state number : 0, 1, ... , (length), (length+1)                  
            for state in self._statelist: 
                if t == 0:                    
                    a[0][state] = 0.0  
                elif t == 1:
                    a[1][state] = self.prob_start[state] * self.prob_emit[state][sequence[0]]  #t = 0, first state. sequence must be list (or tuple)  함수냐 변수냐 그것이 문제                             
                else:
                    sum_prob = 0.0
                    for pre_state in self._statelist:       
                        sum_prob += a[t - 1][pre_state] * self.prob_trans[pre_state][state]
                    a[t][state] = sum_prob * self.prob_emit[state][sequence[t]]

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
    
    def decode(self, sequence):
        a = [{}]
        dec_state = []
        
        for t in range(0, len(sequence), 1): # state number : 0, 1, ... , (length), (length+1)                  
            
            for state in self._statelist: 
                if t == 0:                    
                    a[t][state] = 0.0  
                elif t == 1:
                    a[t][state] = self.prob_start[state] * self.prob_emit[state][sequence[0]]  #t = 0, first state. sequence must be list (or tuple)  함수냐 변수냐 그것이 문제                             
                else:
                    sum_prob = 0.0
                    for pre_state in self._statelist:       
                        sum_prob += a[t - 1][pre_state] * self.prob_trans[pre_state][state]
                    a[t][state] = sum_prob * self.prob_emit[state][sequence[t]]
            
            #여기까지 오면 각 t 상태에서 h, e, _ 중에 최대 값을 알 수 있다.                 
            #그럼 여기서 다음 state 에서 있는 우도를 각각 구한 다음에, 제일 큰걸 현 step 의 state로 정한다  
            if t == 0:
                dec_state[0] = None
            elif t == len(sequence):
                max_state = max(zip(a[t].values(), a[t].keys()))
                dec_state[t] = max_state[1]
            else:
                for post_state in self._statelist:
                    post_prob = a[t][state] * self.prob_trans[state][post_state]
                    max_state = max(zip(post_prob.values(), post_prob.keys()))
                    dec_state[t] = max_state[1]
            '''                
            * self.prob_emit[post_state][sequence[t+1]]
            #다음 스텝에 가는 것 중에 
            
            제일 큰
            maxstate = max(zip(a[t].values(), a[t].keys()))
            for state in self._statelist:            
                dec_temp = self.prob_trans[pre_maxstate][state] * self.prob_emit[state][sequence[t]]
                max(dec_temp) maxstate = max(zip(a[t].values(), a[t].keys())) #셋 state 중에서 가장 큰 것은? h/e/_
            pre_maxstate = maxstate
            '''            
        return dec_state

        

    '''
    def _backward(self, sequence):
        b = [{}] # backward probability, beta (list for [(t+1) states]; dict for probabilities of states {'h', 'e', '_'})            
        
        for t in range(len(sequence) + 1, 1, -1): # state number : 0, 1, ... , (length), (length+1)                  
            for state in self._statelist: 
                if t == 0:                    
                    a[0][state] = 0.0  
                elif t == 1:
                    a[1][state] = self.prob_start[state] * self.prob_emit[state][sequence[0]]  #t = 0, first state. sequence must be list (or tuple)  함수냐 변수냐 그것이 문제                             
                else:
                    sum_prob = 0.0
                    for pre_state in self._statelist:       
                        sum_prob += a[t - 1][pre_state] * self.prob_trans[pre_state][state]
                    a[t][state] = sum_prob * self.prob_emit[state][sequence[t]]

        return a 
    '''    
    ''' 
    output must be a list as below 
    
    a = [ {'h':... , 'e':... , '_':... },    a_{0} = first state * prob_emit = prob_start * prob_emit
                        ...    
          {'h':... , 'e':... , '_':... }, ]  a_{length} = last state * prob_emit
    '''  