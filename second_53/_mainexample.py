# -*- coding: utf-8 -*-
"""
Created on Thu Feb 21 14:44:36 2019

@author: HYEJEONG
"""
#import hmmexample as hmmex
import hiddemarkovmodel_solveprob_butlog as hmm
import sequence_test as seq

#states = b # hidden states, secondary structure
#symbols = a # observable, amino acid sequences

statelist = ['h', 'e', '_']
symbollist = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

pathtrain = '../dataset/protein-secondary-structure.train'
pathtest = '../dataset/protein-secondary-structure.train'

seqtrain = seq.Seq(pathtrain, statelist, symbollist)
seqtest = seq.Seq(pathtest, statelist, symbollist)

[protein, secondstr] = seqtrain.lineseq()
[proteintest, secondstrtest] = seqtest.lineseq()

[prob_start, prob_trans, prob_emit] = seqtrain.prob(protein, secondstr)


modeltrain = hmm.Hmm(prob_start, prob_trans, prob_emit, statelist, symbollist)
modeltest = hmm.Hmm(prob_start, prob_trans, prob_emit, statelist, symbollist)

print('This is a start probabilities : ', prob_start)
print('This is a transition probabilities : ', prob_trans)
print('This is a emission probabilities : ', prob_emit)
# decimal.Decimal
result = []
for sequence, structure in zip(proteintest, secondstrtest):
    modeltest.check(sequence)
    result.append(modeltest.decode(sequence)) 

correctcount = 0
wrongcount = 0
for i in range(len(secondstrtest)):
    for j in range(len(secondstrtest[i])):  
       # print(result[i])
       # print(result[i][1])
        if result[i][1][j] == secondstrtest[i][j]:
            correctcount += 1
        else:
            wrongcount += 1
            
correctpercent = 100 * correctcount / (correctcount + wrongcount)
print(correctpercent)        #52.8%    
'''
with open('textresult.txt', 'w', encoding='utf8') as f:
   f.write(') 
'''    

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