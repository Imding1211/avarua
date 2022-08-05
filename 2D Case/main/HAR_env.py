#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 14:01:18 2021

@author: ding
"""

from solver_HAR import solver

class AI:
    
    def __init__(self, ISAVE, ITMAX, IMETHOD, CFL, alpha_ini, alpha_end, alpha_range, beta_ini, beta_end, beta_range):
        
        self.n_actions_alpha = int(abs(alpha_end-alpha_ini)/alpha_range)+1
        self.n_actions_beta  = int(abs(beta_end-beta_ini)/beta_range)+1
        self.alpha_ini       = alpha_ini
        self.beta_ini        = beta_ini
        self.alpha_range     = alpha_range
        self.beta_range      = beta_range
        self.best_alpha      = 0
        self.best_beta       = 0
        self.max_score       = 0
        self.alpha           = 0
        self.beta            = 0
        self.ISAVE           = ISAVE
        self.ITMAX           = ITMAX
        self.IMETHOD         = IMETHOD
        self.CFL             = CFL
        
        self.Roe             = solver(round(0.9,1), round(0.125,3), round(0.05,2), 2, self.ISAVE, self.ITMAX)[0,1:44,1]
        self.Ausm            = solver(round(0.9,1), round(0.125,3), round(0.1,1), 3, self.ISAVE, self.ITMAX)[0,11:36,19]

    def step(self, act_alpha, act_beta):

        self.alpha = self.alpha_ini + act_alpha*self.alpha_range
        self.beta = self.beta_ini + act_beta*self.beta_range

        Q = solver(round(self.alpha,3), round(self.beta,3), round(self.CFL,3), self.IMETHOD, self.ISAVE, self.ITMAX)
        
        error = abs(sum(Q[0,1:44,1]-self.Roe))+abs(sum(Q[0,11:36,19]-self.Ausm))

        scord = 1/(error+1)*100

        if scord > self.max_score:
            reward = 1
            self.max_score = scord
            self.best_alpha = self.alpha
            self.best_beta = self.beta
        else:
            reward = -1
        
        state_ = 1
        
        return state_, reward, scord, self.alpha, self.beta
