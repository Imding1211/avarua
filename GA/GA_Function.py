# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 11:30:45 2021

@author: Ding
"""

import numpy as np
import random

def Min_Max(oldscord):
        
        max_scord = max(oldscord)
        min_scord = min(oldscord)
        ave_scord = sum(oldscord)
        [x, y] = np.shape(oldscord)
        ave_scord = sum(oldscord)/float(x)
        newscord = np.zeros(x)
        for i in range(x):
            newscord[i] = (oldscord[i] - min_scord) / (max_scord - min_scord) + ave_scord
        newscord = newscord.reshape(-1,1)
        return newscord
    
class GA:
    def __init__(self):
        return
        
    def initpop_binary(self, popsize, chromlength):
            
        ans = [np.random.randint(0, 2, chromlength).tolist() for _ in range(popsize)]
        iniput = np.array(ans)
            
        return iniput
        
    def binary_to_decimal_2(self, pop_binary, LB_1, UB_1, LB_2, UB_2):
    
        [x,y] = np.shape(pop_binary)
        temp_1 = np.zeros((x,1))
        output_1 = np.zeros((x,1))
        temp_2 = np.zeros((x,1))
        output_2 = np.zeros((x,1))
        for i in range(x):
            for j in range(0,int(y/2),1):
                temp_1[i] = temp_1[i] + (2**j)*pop_binary[i,j]
                output_1[i] = (temp_1[i]*(UB_1-LB_1))/(2**(y/2)-1)+LB_1
            n = 0
            for j in range(int(y/2),y,1):
                temp_2[i] = temp_2[i] + (2**n)*pop_binary[i,j]
                output_2[i] = (temp_2[i]*(UB_2-LB_2))/(2**(y/2)-1)+LB_2
                n = n + 1
                
        return output_1 , output_2
    
    def best_ans(self, bespop, besscord):
        
        [x,y] = np.shape(bespop)
        temppop = bespop[0,:]
        tempscord = besscord[0]
        for i in range(1,x,1):
            if (tempscord < besscord[i]):
                temppop = bespop[i,:]
                tempscord = besscord[i]
        minpop = temppop.reshape(1,-1)
        minscord = tempscord
    
        return minpop , minscord
    
    def selection(self, selpop, selscord):
        
        newscord = Min_Max(selscord)
        pickrate = newscord/sum(newscord)
        pickrate = np.cumsum(pickrate)
        [x,y] = np.shape(newscord)
        newselpop = np.zeros(np.shape(selpop))
        ms = np.sort(np.random.rand(x,1), axis=0)
        fitin = 0
        newin = 0
        while ((newin) < x):
            if (ms[newin] < pickrate[fitin]):
                newselpop[newin,:] = selpop[fitin,:];
                newin = newin+1
            else:
                fitin = fitin+1
                 
        return newselpop
    
    def crossover(self, cropop, crosoverrate):
    
        [x,y] = np.shape(cropop)
        newcropop = np.zeros((x,y))
        for i in range(0,int(x-1),2):
            if (crosoverrate > random.random()):
                cropoint = round(random.random()*y)
                if (cropoint == 0):
                    cropoint = 1;
                newcropop[i,:] = np.hstack(([cropop[i,0:cropoint],cropop[i+1,cropoint:y]]))
                newcropop[i+1,:] = np.hstack(([cropop[i+1,0:cropoint],cropop[i,cropoint:y]]))
            else:
                newcropop[i,:] = cropop[i,:]
                newcropop[i+1,:] = cropop[i+1,:]
        
        return newcropop
    
    def mutation(self, mutpop, mutationrate):
        
        [x,y] = np.shape(mutpop)
        newmutpop = mutpop+1-1
        for i in range(x):
            if (mutationrate > random.random()):
                mutationpoint = round(random.random()*y)          
                if (mutationpoint == y):
                    mutationpoint = y-1
                if (newmutpop[i,mutationpoint] == 0):
                    newmutpop[i,mutationpoint] = newmutpop[i,mutationpoint] + 1
                else:
                    newmutpop[i,mutationpoint] = newmutpop[i,mutationpoint] - 1

        return newmutpop