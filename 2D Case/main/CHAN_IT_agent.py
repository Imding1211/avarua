#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 14:13:15 2021

@author: ding
"""

from CHAN_IT_env import AI
import numpy as np
from RL_brain import QLearningTable
from RL_DQbrain import DQN

#==============================================================================
        
class agent_Q:
    
    def __init__(self, ITMAX, ngrd, imethod, alpha_ini, alpha_end, alpha_range, beta_ini, beta_end, beta_range):

        self.ITMAX    = ITMAX
        self.env      = AI(ngrd, imethod, alpha_ini, alpha_end, alpha_range, beta_ini, beta_end, beta_range)
        self.RL_alpha = QLearningTable(actions=list(range(self.env.n_actions_alpha)))
        self.RL_beta  = QLearningTable(actions=list(range(self.env.n_actions_beta)))
    
        return 
        
#==============================================================================
        
    def update(self, epoch):
        
        best_reward = 0
                
        for episode in range(epoch):
            
            alpha_lis  = []
            beta_lis   = []
            reward_lis = []
            
            Q = self.env.inicon()
            
            print("Epoch "+str(episode+1))
            
            for IT in range(self.ITMAX):
                
                print(f"\rIT = [{IT+1}/{self.ITMAX}]",end='')
            
                a_alpha = self.RL_alpha.choose_action(str(IT))
                a_beta = self.RL_beta.choose_action(str(IT))
            
                reward, alpha, beta = self.env.step(Q, a_alpha, a_beta)
                        
                alpha_lis.append(alpha)
                beta_lis.append(round(beta,3))
                reward_lis.append(reward)
                        
                self.RL_alpha.learn(str(IT), a_alpha, reward, str(IT+1))
                self.RL_beta.learn(str(IT), a_beta, reward, str(IT+1))
                
            sum_reward = sum(reward_lis)
            
            print("\n",round(sum_reward,2))
            
            if sum_reward > best_reward:
                best_reward     = sum_reward
                best_alpha_lis  = alpha_lis
                best_beta_lis   = beta_lis
                best_reward_lis = reward_lis
        
        print("Complete Best Result = ",round(best_reward,2))
        
    
        return best_alpha_lis, best_beta_lis, best_reward_lis

#==============================================================================

class agent_DQ:

    def __init__(self, ITMAX, ngrd, imethod, alpha_ini, alpha_end, alpha_range, beta_ini, beta_end, beta_range):

        self.ITMAX    = ITMAX
        self.env      = AI(ngrd, imethod, alpha_ini, alpha_end, alpha_range, beta_ini, beta_end, beta_range)
        self.RL_alpha = DQN(self.env.n_actions_alpha)
        self.RL_beta  = DQN(self.env.n_actions_beta)
    
        return 
        
#==============================================================================
        
    def update(self, epoch):
        
        MEMORY_CAPACITY = 2000
        
        best_reward = 0
                
        for episode in range(epoch):
            
            alpha_lis  = []
            beta_lis   = []
            reward_lis = []
            
            Q = self.env.inicon()
            
            print("Epoch "+str(episode+1))
            
            for IT in range(1,self.ITMAX+1):
                
                print(f"\rIT = [{IT}/{self.ITMAX}]",end='')
                
                s = np.array([IT])
            
                a_alpha = self.RL_alpha.choose_action(s)
                a_beta = self.RL_beta.choose_action(s)
            
                reward, alpha, beta = self.env.step(Q, a_alpha, a_beta)
                
                s_ = np.array([IT+1])
                        
                alpha_lis.append(alpha)
                beta_lis.append(round(beta,3))
                reward_lis.append(reward)
                        
                self.RL_alpha.store_transition(s, a_alpha, reward, s_)
                self.RL_beta.store_transition(s, a_beta, reward, s_)
                
                if self.RL_alpha.memory_counter > MEMORY_CAPACITY:
                    self.RL_alpha.learn()
                    self.RL_beta.learn()
                
            sum_reward = sum(reward_lis)
            
            print("\n",round(sum_reward,2))
            
            if sum_reward > best_reward:
                best_reward     = sum_reward
                best_alpha_lis  = alpha_lis
                best_beta_lis   = beta_lis
                best_reward_lis = reward_lis
        
        print("Complete Best Result = ",round(best_reward,2))
        
    
        return best_alpha_lis, best_beta_lis, best_reward_lis
       
#==============================================================================