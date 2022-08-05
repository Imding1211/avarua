#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 14:13:15 2021

@author: ding
"""

from PDE_env import AI
from RL_DQbrain import DQN
import numpy as np

class agent:
    
    def __init__(self, ITMAX, isave):

        self.env = AI(ITMAX, isave)
        self.RL_t = DQN(self.env.n_actions_t)
        self.RL_time = DQN(self.env.n_actions_time)
    
        return 
        
#==============================================================================
        
    def update(self, epoch):
        
        MEMORY_CAPACITY = 2000

        T_lis=[]
        t_lis=[]
        time_lis=[]
        
        for episode in range(epoch):
            
            print("Epoch "+str(episode+1))
            
            s = np.array([0])
            
            a_t = self.RL_t.choose_action(s)
            a_time = self.RL_time.choose_action(s)
        
            s_, reward, thrust, t, time = self.env.step(a_t, a_time)
            
            print(thrust, t, time)
                    
            T_lis.append(thrust)
            t_lis.append(t)
            time_lis.append(time)
                    
            self.RL_t.store_transition(s, a_t, reward, s_)
            self.RL_time.store_transition(s, a_time, reward, s_)

            if self.RL_t.memory_counter > MEMORY_CAPACITY:
                self.RL_t.learn()
                self.RL_time.learn()
        
        print(self.env.max_thrust, self.env.max_t, self.env.max_time)
        
        return T_lis, t_lis, time_lis, self.env.max_thrust, self.env.max_t, self.env.max_time
