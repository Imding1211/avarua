#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 14:13:15 2021

@author: ding
"""

from PDE_env import AI
from RL_brain import QLearningTable
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

class agent:
    
    def __init__(self, ITMAX, isave):

        self.env = AI(ITMAX, isave)
        self.RL_t = QLearningTable(actions=list(range(self.env.n_actions_t)))
        self.RL_time = QLearningTable(actions=list(range(self.env.n_actions_time)))
    
        return 
        
#==============================================================================
        
    def update(self, epoch):
        
        T_lis=[]
        t_lis=[]
        time_lis=[]
        
        for episode in range(epoch):
            
            #print("Epoch "+str(episode+1))
            
            observation = 0
            
            a_t = self.RL_t.choose_action(str(observation))
            a_time = self.RL_time.choose_action(str(observation))
        
            observation_, reward, thrust, t, time = self.env.step(a_t, a_time)
            
            #print(thrust, t, time)
                    
            T_lis.append(thrust)
            t_lis.append(t)
            time_lis.append(time)
                    
            self.RL_t.learn(str(observation), a_t, reward, str(observation_))
            self.RL_time.learn(str(observation), a_time, reward, str(observation_))
    
            observation = observation_
        
        #print(self.env.max_thrust, self.env.max_t, self.env.max_time)
        
        return T_lis, t_lis, time_lis, self.env.max_thrust, self.env.max_t, self.env.max_time
