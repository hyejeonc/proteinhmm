# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 17:25:48 2019

@author: HYEJEONG

Problem with float --> change to log 
"""
#from math import log

#def _gamma(sequence, statelist, symbollist):
#    a = 
    
def train(proteins)
    statelist = ['h', 'e', '_']
    symbollist = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    
    pre_prob = 0
    for protein in proteins:
        prob += log(pre_prob)
    

 
#def em_prob_emit(sequence, statelist, symbollist):    
    
class Hmm(object):
    
    def __init__(self, prob_start, prob_trans, prob_emit, statelist, symbollist):

        self._prob_start = prob_start
        self._prob_trans = prob_trans
        self._prob_emit = prob_emit
        self._statelist = statelist
        self._symbollist = symbollist
        
    def initialstart(self):   
        for state in self._statelist:
            for post_state in self._statelist:
                self._prob_trans[state][post_state] = 0.0
        
            for symbol in self._symbollist:
                self._prob_emit[state][symbol] = 0.0
                
            self._prob_start[state] = 0.0
  
    def emupdate(self, sequence):
        print('this is prob_trans ', self._prob_trans)
        print('this is prob_emit ', self._prob_emit) 
        a = self.forward(sequence)
        b = self.backward(sequence)
        print('this is a ', a)
        print('this is b ', b)
        g = []
        sum_prob = 0.0
        for t in range(0, len(sequence)): # 0, 1, ... l, 
            g.append({})
            for state in self._statelist:
                g[t][state] = a[t][state] * b[t][state]
                
            sum_prob += g[t][state]            

            for state in self._statelist:
                if sum_prob > 1e-300:
                    g[t][state] /= sum_prob
              #  else:
               #     g[t][state] = 0
        #여기까지 n 번째 state 서  gamma 가 몇인지를 구함 
        #print('this is g', g)
        
        x = []
        for t in range(0, len(sequence)-1): #until t-1 state
            x.append({})
            sum_prob = 0.0
            for pre_state in self._statelist:            
                x[t][pre_state] = {}
                for post_state in self._statelist:
                    
                    '''
                    print('this is t : ', t)
                    print(a[t][pre_state])
                    print(self._prob_trans)
                    print(self._prob_trans[pre_state][post_state])
                    print(self._prob_emit[post_state][sequence[t+1]])
                    print(b[t+1][post_state])
                    '''
                    x[t][pre_state][post_state] = a[t][pre_state]  \
                        * self._prob_trans[pre_state][post_state] \
                        * self._prob_emit[post_state][sequence[t+1]] \
                        * b[t+1][post_state] 
                    sum_prob += x[t][pre_state][post_state]    
            #print('this is sum_prob', sum_prob)         
            if sum_prob > 1e-300:
                for pre_state in self._statelist:
                    for post_state in self._statelist:
                        x[t][pre_state][post_state] /= sum_prob 
        #여기까지 해서 xi 를 구함. 
        #print()
        #print('this is x', x)
        # print('this is g', g)
        
        '''
       #for Transition probabilities 
        for state in self._statelist:
            sum_g = 0.0
            for t in range(0, len(sequence)-1):         
                sum_g += g[t][state] 
            print('this is sum g t state in trans update', sum_g)# sigma gamma
            
            if sum_g > 1e-300:
                #sum_sum_x = 0.0
                for post_state in self._statelist:
                    sum_x = 0.0
                    for t in range(len(sequence)-1):
                        sum_x += x[t][state][post_state]   #분자
                    
                    print('this is sum x in trans update', sum_x)
                    print('this is sum_g in trans update', sum_g)
                    
                #sum_sum_x += sum_x 
                    
                    self._prob_trans[state][post_state] = sum_x / sum_g #왜 gamma 로 나눌까!      
            else:  
                for post_state in self._statelist:
                    self._prob_trans[state][post_state] = 0
        '''                
        #for state in self._statelist:
        #    self._prob_start[state] = g[0][state] #이 부분에 따라 start 필요하나 안필요하냐가 달라진다 

        
        for state in self._statelist:
               
            for post_xstate in self._statelist:   
                sum_sum_prob = 0.0 
                sum_x = 0.0
                for t in range(0, len(sequence)-1): #분모 
                    sum_prob= 0.0
                    for post_state in self._statelist:                
                        sum_prob += x[t][state][post_state]     
                    sum_sum_prob += sum_prob 

                    sum_x += x[t][state][post_xstate] #분자
                    
                if sum_sum_prob > 1e-300:
                    self._prob_trans[state][post_xstate] = sum_x / sum_sum_prob
                else: 
                    self._prob_trans[state][post_xstate] = sum_x 
                    
                
                #print('this is prob trans in update :', self._prob_trans)
            
            #for Emission probabilities
            #for state in self._statelist:
            sum_g = 0.0
            for t in range(0, len(sequence)):         
                sum_g += g[t][state]
            sum_g += g[len(sequence)-1][state]

            sum_g_emit = {}
            for symbol in self._symbollist:
                sum_g_emit[symbol] = 0.0
            #print(sum_g_emit)
                
            for t in range(0, len(sequence)):
                sum_g_emit[sequence[t]] += g[t][state] 
               #
                #print('this is g state', g[t][state])# 이 부분이 이해가 안된다. 골라서 더하는건가.  각 symbol 같은게 나올 때 마다 해당 state 의 gamma 를 더한다. 
            #print(sum_g_emit) 
            if sum_g > 1e-300:
                for symbol in self._symbollist:
                    self._prob_emit[state][symbol] = sum_g_emit[symbol] / sum_g       
            else:
                for symbol in self._symbollist:               
                    self._prob_emit[state][symbol] = 0.0
                
        return self._prob_emit, self.prob_trans
        
                    
        
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
        
        a = [] # forward probability, alpha (list for [(t+1) states]; dict for probabilities of states {'h', 'e', '_'})            
                 # alpha is saved for next alpha, dynamics
        c = [{}] # scaling factor for very small number   
        
        for t in range(0, len(sequence), +1): # state number : 0, 1, ... , (length)               
            a.append({})
            if t == 0: 
                for state in self._statelist:                 
                    
                    a[t][state] = self.prob_start(state) * self.prob_emit(state, sequence[0])  #t = 0, first state. sequence must be list (or tuple)  함수냐 변수냐 그것이 문제                               
                
#            elif t == len(sequence)+1:
#                a[t][state] = 1.0
                
            else:
                for state in self._statelist:
                    
                    sum_prob = 0.0
                    for pre_state in self._statelist:                     
                        sum_prob += a[t-1][pre_state] * self.prob_trans(pre_state, state)
                        #pre_state = state

                    #print('this is sum_prob ', sum_prob)
                    a[t][state] = sum_prob * self.prob_emit(state, sequence[t-1])
           #print('sum _prob  in forward', sum(list(a[t] .values())))
          #  print('a ', a[t])

            
            if sum(list(a[t] .values())) > 1e-300:
               # print('this is sum a : ', sum(list(a[t].values())))

                c = 1 / sum(list(a[t].values()))
                for state in self._statelist:
                    a[t][state] = c * a[t][state]
            #print('a after scaling', a[t])
            # print('this is a sum in forward', sum(list(a[t].values())))        
           # print('this is a in forward', a[t])
                    
            
           # print('this is t in forward ', t)
           # print('this is a[t][state] ', a[t][state])
        return a 

    def backward(self, sequence):        
        b = [] # backward probability, beta (list for [(t+1) states]; dict for probabilities of states {'h', 'e', '_'})            
                  # beta is saved for next beta, dynamics programming
        c = [] # scaling factor for very small number   
        
        for t in range(len(sequence)-1, -1, -1): # state number : 0, 1, ... , (length)               
            #print('this is t in back', t)
            b.insert(0, {})
                    
            for state in self._statelist: 
                if t == (len(sequence)-1):
                    b[0][state] = 1.0        
                else:      
                    sum_prob = 0.0
                    for post_state in self._statelist:                     
                        sum_prob += b[1][post_state] * self.prob_trans(state, post_state)
                  
                    b[0][state] = sum_prob * self.prob_emit(state, sequence[t])

            if sum(list(b[0].values())) != 0:
                c = 1 / sum(list(b[0].values()))
                for state in self._statelist:
                    b[0][state] = c * b[0][state]
            else:
                continue        
            
        return b    
    
    ''' 
    output must be a list as below 
    
    a = [ {'h':... , 'e':... , '_':... },    a_{0} = first state * prob_emit = prob_start * prob_emit
                        ...    
          {'h':... , 'e':... , '_':... }, ]  a_{length} = last state * prob_emit
    '''    
    '''    
        def probability(self, sequence): #a 전방확률의 나중에 구한 것 합 
            a = self._forward(sequence)
            sum_prob = 0.0
            
            for state in a[len(sequence)]:
                sum_prob += a[len(sequence)][state]
            
            return sum_prob
    
    '''     
    
    def decode(self, sequence): #viterbi 
    
        v = []
        dec_state = []
        dec_prob = []

        for t in range(0, len(sequence)): # state number : 0, 1, ... , (length)               
            
            v.append({})
            #print('this is start prob in decode', self._prob_start)
            if t == 0: 
                for state in self._statelist:                 
                    v[t][state] = self.prob_start(state) * self.prob_emit(state, sequence[0]) 
            else:                
                for state in self._statelist:       
                    v[t][state] = vmax[0] * self.prob_trans(vmax[1], state) * self.prob_emit(state, sequence[t-1])

            vmax = list(max(zip(v[t].values(), v[t].keys())))            
            if vmax[0] != 0:
                c = 1 / sum(list(v[t].values()))
                for state in self._statelist:
                   v[t][state] = c * v[t][state]
               
                vmax[0] = c * vmax[0]#decode.append(vmax)
                dec_prob.append(vmax[0])           
                dec_state.append(vmax[1])
           
        return v[t], dec_state, vmax[0]

                   


''' 
    output must be a list as below 
    
    a = [ {'h':... , 'e':... , '_':... },    a_{0} = first state * prob_emit = _prob_start * prob_emit
                        ...    
          {'h':... , 'e':... , '_':... }, ]  a_{length} = last state * prob_emit
''' 