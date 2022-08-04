#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 14:01:18 2021

@author: ding
"""

import numpy as np
import math
from solver import solver, ini

class AI:
    
    def __init__(self, icase, mesh, imethod):
        self.mesh64 = np.loadtxt('12800.csv',delimiter=",",dtype=np.float32,skiprows=1)[:,1]
        self.n_actions_alpha = 31
        self.n_actions_beta = 101
        self.icase = icase
        self.imethod = imethod
        self.mesh = mesh
        
    def reset(self):
        self.mesh8, self.T, _ = ini(self.icase, self.mesh)
        self.T -= 1
        self.done = False
        self.alpha = 100
        self.beta = 1.3
        self.nt = 0
        state_ = self.nt
        return state_, self.nt, self.T
        
    def step(self, act_alpha, act_beta):
        
        self.alpha = 200 + act_alpha*10
        self.beta = 1.2 + act_beta*0.001
        
        self.mesh8 = solver(self.icase, self.imethod, self.mesh, self.alpha, self.beta, self.mesh8)
        
        x1 = self.mesh64[0:12800:16]
        x2 = self.mesh8[0:self.mesh*100,0]
        reward = -np.sqrt(np.mean((x1-x2)**2))
        
        if self.nt == self.T:
            state_ = 'terminal'
            self.done = True
        else:
            self.nt += 1
            state_ = self.nt
            self.done = False
        
        return state_, reward, self.done, self.nt, x2, self.alpha, self.beta
