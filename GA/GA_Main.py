# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 11:41:27 2021

@author: Ding
"""

from GA_Function import *

ga = GA()

inipop = ga.initpop_binary(200, 100)

for i in range(10):
    
    [x, y] = ga.binary_to_decimal_2(inipop, 0, 10, 0, 10)
    
    scord = (x+y)
    
    [maxpop , maxscord] = ga.best_ans(inipop, scord)
    
    secpop = ga.selection(inipop, scord)
    
    cropop = ga.crossover(secpop, 1.0)
    
    mutpop = ga.mutation(cropop, 1.0)

    inipop = mutpop
    
[x, y] = ga.binary_to_decimal_2(maxpop, 0, 10, 0, 10)
print('\n',x[0,0], y[0,0])