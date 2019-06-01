# -*- coding: utf-8 -*-
#ppt : https://people.cs.umass.edu/~mccallum/courses/inlp2004a/lect10-hmm2.pdf
#https://sambaiga.github.io/2017/05/03/hmm-intro.html
#http://www.biostat.jhsph.edu/bstcourse/bio638/notes/HMMs_BaumWelch.pdf

#filtering smoothing, ..?https://danieltakeshi.github.io/2015-07-25-hidden-markov-models-and-particle-filtering/
"""
Created on Mon Mar  4 17:25:48 2019

@author: HYEJEONG

Problem with float --> change to log 
"""
import sequence as seq
import statistics as stat

from math import log












def initialize(initialset):
    prob_start = initialset[0]
    prob_trans = initialset[1]
    prob_emit = initialset[2]
    statelist = initialset[3]
    symbollist = initialset[4]
    return Hmm(prob_start, prob_trans, prob_emit, statelist, symbollist)

#MLE(Maximum Likelihood Estimation)    
def train(initialset, proteins, secondstrs, d=0.001, maxiter = 100): #, prob_start=None, prob_trans=None, prob_emit=None):
    #statelist = ['h', 'e', '_']
    #symbollist = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    update = initialize(initialset)

    pre_prob = 0.0
    for protein in proteins:
        prob = update.viterbi(protein)
        pre_prob += log(prob[2]) #expectation
    pre_prob /= len(proteins)
    print(pre_prob/len(proteins))
    
    for i in range(maxiter): 
        post_prob = 0.0
        print(i)
        for protein in proteins:
            #E-M update transition probabilites and emission probabilities 
            update.emupdate(protein)
            prob = update.viterbi(protein) #viterbi probability of 
            post_prob += log(prob[2])
        post_prob /= len(proteins)
        print(pre_prob, post_prob)
        #if (post_prob > pre_prob) and (abs(post_prob - pre_prob) < d):
        #    break
        #if post_prob <= pre_prob:   
        pre_prob = post_prob 
    #a = (update.prob_start, update.prob_trans, update.prob_emit, update.statelist, update.symbollist)    
    #Hmm(update.prob_start, update.prob_trans, update.prob_emit, update.statelist, update.symbollist)    
        
    '''
    for i in range(maxiter): 
        post_prob = 0.0
        print(i)
        for protein in proteins:
            update.emupdate(protein)#E-M update transition probabilites and emission probabilities 
            prob = update.viterbi(protein) #viterbi probability of 
            post_prob += log(prob[2]) #prob = [v[t], dec_state, vmax[0]]
            
#https://www.cs.cmu.edu/~epxing/Class/10701-08s/recitation/em-hmm.pdf      ->   
    
        print((pre_prob - post_prob)/len(proteins))
        if abs(pre_prob - post_prob)/len(proteins) < d:
            print('this is post', post_prob)
            print('this is pre', pre_prob)
            break
        post_prob = pre_prob
        
    '''    
    return update
    #return Hmm(update.prob_start, update.prob_trans, update.prob_emit, update.statelist, update.symbollist)
 
#def em_prob_emit(sequence, statelist, symbollist):    
    
class Hmm(object):
    
    def __init__(self, prob_start, prob_trans, prob_emit, statelist, symbollist):
        self._prob_start = prob_start
        self._prob_trans = prob_trans
        self._prob_emit = prob_emit
        self._statelist = list(statelist)
        self._symbollist = list(symbollist)

    def get(self):
        return self._prob_start, self._prob_trans, self._prob_emit, self._statelist, self._symbollist
        
    def zero(self):   
        self._prob_trans = []
        self._prob_emit = []
        self._prob_start = []
        
        for state in self._statelist:
            self._prob_trans.append({})
            self._prob_emit.append({})
            self._prob_start.append({})
            for post_state in self._statelist:
                self._prob_trans[state][post_state] = 0.0
                
            for symbol in self._symbollist:
                self._prob_emit[state][symbol] = 0.0
                
            self._prob_start[state] = 0.0
  
    def check(self, sequence):
        if sequence not in self._symbollist:
            print('Not available sequence')
 
    def prob_start(self):
      #  if state not in self._statelist:
     #       print('Not available states, not in state list')
        if self._prob_start == None:    
            return 0
        #print(state)
        #print(self._prob_start)    
        #print(self._prob_start[state])
        return self._prob_start
 
    def prob_trans(self):
        if self._prob_trans == None:    
            return 0   
        return self._prob_trans    
    
    def prob_emit(self):
        if self._prob_emit == None:    
            return 0      
        return self._prob_emit
    
    def statelist(self):
        if self._statelist == None:    
            return None     
        return self._statelist
    
    def symbollist(self):
        if self._symbollist == None:    
            return None     
        return self._symbollist    
    
    
            
    def emupdate(self, sequence):  #update, Baum-Welch
        a = self.forward(sequence)
        b = self.backward(sequence)
        
        # Expectation value (E-step)
        g = [] # get gamma 
        g_smooth = []
        sum_prob = 0.0
        for t in range(0, len(sequence)): # 0, 1, ... l, 
            g.append({})
            for state in self._statelist:
                g[t][state] = a[t][state] * b[t][state]
                sum_prob += g[t][state]    

            g_smooth.append(sum_prob)
            if sum_prob > 1e-300:
                for state in self._statelist:
                    g[t][state] /= sum_prob

        x = [] # get xi
        x_smooth = []
        for t in range(0, len(sequence)-1): #until t-1 state
            x.append({})
            sum_prob = 0.0
            for pre_state in self._statelist:            
                x[t][pre_state] = {}
                for post_state in self._statelist:
                    x[t][pre_state][post_state] = a[t][pre_state]  \
                        * self._prob_trans[pre_state][post_state] \
                        * self._prob_emit[post_state][sequence[t+1]] \
                        * b[t+1][post_state] 
                    sum_prob += x[t][pre_state][post_state]    
                x_smooth.append(sum_prob)
            
            if sum_prob > 1e-300:
                for pre_state in self._statelist:
                    for post_state in self._statelist:
                        x[t][pre_state][post_state] /= sum_prob 


        #Maximization step (M-step)
        #Laplace smoothing 
        p = 0.001
        for state in self._statelist:
            #for start probabilities
            #p_smooth = 1 / (len(self._statelist) + )
            self._prob_start[state] = (g[0][state] + p) / (1 + p * len(self._statelist))   
            #     self._prob_start[state] = (g[0][state] + 1) / (g_smooth[0] + len(self._statelist))
            
            #for transition probabilities   
            for post_xstate in self._statelist:   
                sum_sum_prob = 0.0 
                sum_x = 0.0
                for t in range(0, len(sequence)-1): #분모 
                    sum_prob= 0.0
                    for post_state in self._statelist:                
                        sum_prob += x[t][state][post_state]     
                    sum_sum_prob += sum_prob 

                    sum_x += x[t][state][post_xstate] #분자
                    
               # if sum_sum_prob > 1e-300:
                    self._prob_trans[state][post_xstate] = (sum_x + p) / (sum_sum_prob + p*len(self._statelist))
               # else: 
               #     self._prob_trans[state][post_xstate] = 1e-300
                    
                
                #print('this is prob trans in update :', self._prob_trans)
            
            #for emission probabilities
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
            
            for symbol in self._symbollist:
                self._prob_emit[state][symbol] = (sum_g_emit[symbol] + p) / (sum_g + p*len(self._statelist))  

                
            #print(self._prob_emit, self.prob_trans)
        return self._prob_emit, self._prob_trans
                
    def forward(self, sequence):
        ''' 
        output must be a list as below 
        
        a = [ {'h':... , 'e':... , '_':... },    a_{0} = first state * prob_emit = prob_start * prob_emit
                            ...    
              {'h':... , 'e':... , '_':... }, ]  a_{length} = last state * prob_emit
        '''    
         
        a = [] # forward probability, alpha (list for [(t+1) states]; dict for probabilities of states {'h', 'e', '_'})            
                 # alpha is saved for next alpha, dynamics
        c = [{}] # scaling factor for very small number   
        
        for t in range(0, len(sequence), +1): # state number : 0, 1, ... , (length)               
            a.append({})
            if t == 0: 
                for state in self._statelist:                 
                    
                    a[t][state] = self._prob_start[state] * self._prob_emit[state][sequence[0]]  #t = 0, first state. sequence must be list (or tuple)  함수냐 변수냐 그것이 문제                               
                
#            elif t == len(sequence)+1:
#                a[t][state] = 1.0
                
            else:
                for state in self._statelist:
                    
                    sum_prob = 0.0
                    for pre_state in self._statelist:                     
                        sum_prob += a[t-1][pre_state] * self._prob_trans[pre_state][state]
                        #pre_state = state

                    #print('this is sum_prob ', sum_prob)
                    a[t][state] = sum_prob * self._prob_emit[state][sequence[t-1]]
           #print('sum _prob  in forward', sum(list(a[t] .values())))
          #  print('a ', a[t])

            
            if sum(list(a[t] .values())) > 1e-300:
               # print('this is sum a : ', sum(list(a[t].values())))

                c = 1 / sum(list(a[t].values()))
                for state in self._statelist:
                    a[t][state] *= c 
   
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
                        sum_prob += b[1][post_state] * self._prob_trans[state][post_state]
                  
                    b[0][state] = sum_prob * self._prob_emit[state][sequence[t]]

            if sum(list(b[0].values())) > 1e-300:
                c = 1 / sum(list(b[0].values()))
                for state in self._statelist:
                    b[0][state] *= c 
       
            
        return b    
    
 
    
    def viterbi(self, sequence): #viterbi 

        v = []
        dec_state = []
        dec_prob = []

        for t in range(0, len(sequence)): # state number : 0, 1, ... , (length)               
            
            v.append({})
            #print('this is start prob in viterbi', self._prob_start)
            if t == 0: 
                for state in self._statelist:                 
                    v[t][state] = self._prob_start[state] * self._prob_emit[state][sequence[0]] 
            else:                
                for state in self._statelist:       
                    v[t][state] = vmax[0] * self._prob_trans[vmax[1]][state] * self._prob_emit[state][sequence[t-1]]

            vmax = list(max(zip(v[t].values(), v[t].keys())))            
            
            if vmax[0] > 1e-300 :
                c = 1 / sum(list(v[t].values()))
                for state in self._statelist:
                   v[t][state] *= c 
                vmax[0] *= c #viterbi.append(vmax)
             
            dec_prob.append(vmax[0])           
            dec_state.append(vmax[1])        

        return [v[t], dec_state, vmax[0]]
    

        ''' 
            output must be a list as below 
            
            a = [ {'h':... , 'e':... , '_':... },    a_{0} = first state * prob_emit = _prob_start * prob_emit
                                ...    
                  {'h':... , 'e':... , '_':... }, ]  a_{length} = last state * prob_emit
        ''' 