#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 14:01:18 2021

@author: ding
"""

from solver_CHAN_IT import solver,inic

class AI:
    
    def __init__(self, ngrd, imethod, alpha_ini, alpha_end, alpha_range, beta_ini, beta_end, beta_range):
        
        self.n_actions_alpha = int(abs(alpha_end-alpha_ini)/alpha_range)+1
        self.n_actions_beta  = int(abs(beta_end-beta_ini)/beta_range)+1
        
        self.ngrd            = ngrd
        self.imethod         = imethod
        
        self.alpha_ini       = alpha_ini
        self.beta_ini        = beta_ini
        
        self.alpha_range     = alpha_range
        self.beta_range      = beta_range

#==============================================================================
          
    def inicon(self):
        
        Q = inic(self.ngrd)
        
        return Q
        
#==============================================================================

    def step(self, Q, act_alpha, act_beta):

        alpha = self.alpha_ini + act_alpha*self.alpha_range
        beta = self.beta_ini + act_beta*self.beta_range
        
        Q = solver(Q, alpha, beta, self.imethod)
        
        rew = 1
        
        for QQ in Q[0,1:401,42]:
            if QQ > 0.00121 and QQ < 0.0018:
                rew = rew + 1
        
        reward = round((1/rew*100),2)

        return reward, alpha, beta
