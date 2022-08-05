#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 14:01:18 2021

@author: ding
"""

from solver import solver

class AI:
    
    def __init__(self, ITMAX, isave):
        
        self.n_actions_t    = 11
        self.n_actions_time = 31
        self.ITMAX          = ITMAX
        self.isave          = isave
        self.max_thrust     = 0
        self.max_t          = 0
        self.max_time       = 0
    
    def step(self, act_t, act_time):

        self.ini_t = 200 + act_t*10
        self.ini_time = 0.5 + act_time*0.01
        
        self.thrust = solver(self.ITMAX, self.isave, self.ini_t, self.ini_time)
        
        if self.thrust > self.max_thrust:
            reward = 1
            self.max_thrust = self.thrust
            self.max_t = self.ini_t
            self.max_time = self.ini_time
        else:
            reward = -1
        
        state_ = 1
        
        return state_, reward, self.thrust, self.ini_t, self.ini_time
