#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 13:27:26 2021

@author: ding
"""

from solver import solver
import pandas as pd

# =============================================================================

def creat_csv(value):

    input_value = pd.DataFrame(columns=['value'])
    
    input_value.to_csv('Thrust.csv')
    
    return

# =============================================================================

ITMAX = 20
isave = 10

ini_t = 300.1
ini_time = 0.7

max_thrust = solver(ITMAX,isave,ini_t,ini_time)

print(max_thrust)