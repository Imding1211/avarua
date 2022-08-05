#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 14:01:18 2021

@author: ding
"""

from solver_CHAN_ITMAX import solver

class AI:
    
    def __init__(self, ISAVE, ITMAX, IMETHOD, ngrd, alpha_ini, alpha_end, alpha_range, beta_ini, beta_end, beta_range):
        
        self.n_actions_alpha = int(abs(alpha_end-alpha_ini)/alpha_range)+1
        self.n_actions_beta  = int(abs(beta_end-beta_ini)/beta_range)+1
        self.alpha_ini       = alpha_ini
        self.beta_ini        = beta_ini
        self.alpha_range     = alpha_range
        self.beta_range      = beta_range
        self.best_alpha      = 0
        self.best_beta       = 0
        self.min_error       = 100
        self.alpha           = 100
        self.beta            = 1.3
        self.ISAVE           = ISAVE
        self.ITMAX           = ITMAX
        self.IMETHOD         = IMETHOD
        self.ngrd            = ngrd
    
    def step(self, act_alpha, act_beta):

        self.alpha = self.alpha_ini + act_alpha*self.alpha_range
        self.beta = self.beta_ini + act_beta*self.beta_range
        
        Q = solver(self.alpha, self.beta, self.ISAVE, self.ITMAX, self.IMETHOD, self.ngrd)
        
        error = sum(Q[0,int(130+(self.ngrd-1)*130):int(138+(self.ngrd-1)*130),int(42*self.ngrd)])
        
        if error < self.min_error:
            reward = 1
            self.min_error = error
            self.best_alpha = self.alpha
            self.best_beta = self.beta
        else:
            reward = -1
        
        state_ = 1
        
        return state_, reward, error, self.alpha, self.beta
