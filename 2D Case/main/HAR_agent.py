#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 14:13:15 2021

@author: ding
"""

from HAR_env import AI
from RL_brain import QLearningTable
from RL_DQbrain import DQN
import numpy as np

#==============================================================================

class agent_Q:
    
    def __init__(self, ISAVE, ITMAX, IMETHOD, CFL, alpha_ini, alpha_end, alpha_range, beta_ini, beta_end, beta_range):

        self.env = AI(ISAVE, ITMAX, IMETHOD, CFL, alpha_ini, alpha_end, alpha_range, beta_ini, beta_end, beta_range)
        self.RL_alpha = QLearningTable(actions=list(range(self.env.n_actions_alpha)))
        self.RL_beta = QLearningTable(actions=list(range(self.env.n_actions_beta)))
    
        return 
        
#==============================================================================
        
    def update(self, epoch):
        
        alpha_lis=[]
        beta_lis=[]
        error_lis=[]
                
        for episode in range(epoch):
            
            print("Epoch "+str(episode+1))
            
            observation = 0
            
            a_alpha = self.RL_alpha.choose_action(str(observation))
            a_beta = self.RL_beta.choose_action(str(observation))
        
            observation_, reward, error, alpha, beta = self.env.step(a_alpha, a_beta)
            
            print(round(alpha,3), round(beta,3), round(error,6))
                    
            alpha_lis.append(alpha)
            beta_lis.append(beta)
            error_lis.append(error)
                    
            self.RL_alpha.learn(str(observation), a_alpha, reward, str(observation_))
            self.RL_beta.learn(str(observation), a_beta, reward, str(observation_))
    
            observation = observation_
        
        print("Best result")
        print(round(self.env.best_alpha,3), round(self.env.best_beta,3), round(self.env.max_score,6))
        
        return alpha_lis, beta_lis, error_lis, self.env.best_alpha, self.env.best_beta, self.env.max_score
        
#==============================================================================

class agent_DQ:
    
    def __init__(self, ISAVE, ITMAX, IMETHOD, CFL, alpha_ini, alpha_end, alpha_range, beta_ini, beta_end, beta_range):

        self.env = AI(ISAVE, ITMAX, IMETHOD, CFL, alpha_ini, alpha_end, alpha_range, beta_ini, beta_end, beta_range)
        self.RL_alpha = DQN(self.env.n_actions_alpha)
        self.RL_beta  = DQN(self.env.n_actions_beta)
    
        return 
        
#==============================================================================
        
    def update(self, epoch):

        MEMORY_CAPACITY = 2000
        
        alpha_lis=[]
        beta_lis=[]
        error_lis=[]
                
        for episode in range(epoch):
            
            print("Epoch "+str(episode+1))
            
            s = np.array([0])
            
            a_alpha = self.RL_alpha.choose_action(s)
            a_beta = self.RL_beta.choose_action(s)
        
            s_, reward, error, alpha, beta = self.env.step(a_alpha, a_beta)
            
            print(round(alpha,3), round(beta,3), round(error,6))
                    
            alpha_lis.append(alpha)
            beta_lis.append(beta)
            error_lis.append(error)
                    
            self.RL_alpha.store_transition(s, a_alpha, reward, s_)
            self.RL_beta.store_transition(s, a_beta, reward, s_)
        
        print("Best result")
        print(round(self.env.best_alpha,3), round(self.env.best_beta,3), round(self.env.max_score,6))
        
        return alpha_lis, beta_lis, error_lis, self.env.best_alpha, self.env.best_beta, self.env.max_score
        
#==============================================================================