#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 14:00:38 2021

@author: ding
"""

from PDE_agent import agent
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

#==============================================================================

epoch = 10
ITMAX = 177000
isace = 177000

#==============================================================================   

P = agent(ITMAX, isace)

thrust_P, t_P, time_P, max_thrust_P, max_t_P, max_time_P = P.update(epoch)

data = {
        "Temperature" : t_P,
        "Cycle Time"  : time_P,
        "Thrust"      : thrust_P,
        "MAX Temperature" : max_t_P,
        "MAX Cycle Time"  : max_time_P,
        "MAX Thrust"      : max_thrust_P
        }

output = pd.DataFrame(data).astype(float)
output.to_csv("Error.csv")
