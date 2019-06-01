# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 10:23:37 2019

@author: HYEJEONG
"""

def count1d(count, item):
    if item not in count:
        count[item] = 0
    count[item] += 1

def count2d(count, item1, item2):
    if item1 not in count:
        count[item1] = {}
    count1d(count[item1], item2)

def norm1d(prob, item_set):  
    result = {} 
    prob_sum = 0.0
    
    for item in item_set:
        prob_sum += prob.get(item, 0)
    
    if prob_sum == 0:
        result[item] = 0
    else:    
        for item in item_set:
            result[item] = prob.get(item, 0) / prob_sum            
    return result
                       
def norm2d(prob, item_set1, item_set2):
    result = {}
    
    if prob is None:
        for item in item_set1:
            result[item] = norm1d(None, item_set2)
    for item in item_set1:
        result[item] = norm1d( prob.get(item, 0), item_set2)
    
    return result