#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 14:00:38 2021

@author: ding
"""

from CHAN_ITMAX_agent import agent_DQ,agent_Q
from solver_CHAN_ITMAX import solver
from load import loadsetting
import pandas as pd
import os
import shutil
import glob

"""
IMETHOD = 1 ORI
IMETHOD = 2 JCP
IMETHOD = 3 ATM
IMETHOD = 4 MUSCL
IMETHOD = 5 THINCEM

ngrd = 1 400x 80
ngrd = 2 800x160

Brain   = 1 Q
Brain   = 2 DQ
"""

#================參數設定======================================================

path = 'CHAN_setting.txt'
par = loadsetting(path,0)

Brain   = int(par[0])
epoch   = int(par[1])
IMETHOD = int(par[2])
ISAVE   = int(par[3])
ITMAX   = int(par[4])
ngrd    = int(par[5])

alpha_ini   = float(par[6])
alpha_end   = float(par[7])
alpha_range = float(par[8])

beta_ini    = float(par[9])
beta_end    = float(par[10])
beta_range  = float(par[11])

#==============================================================================

if Brain == 1:
    P = agent_Q(1075*ngrd, 1075*ngrd, IMETHOD, ngrd, alpha_ini, alpha_end, alpha_range, beta_ini, beta_end, beta_range)
elif Brain == 2:
    P = agent_DQ(1075*ngrd, 1075*ngrd, IMETHOD, ngrd, alpha_ini, alpha_end, alpha_range, beta_ini, beta_end, beta_range)

alpha_lis, beta_lis, error_lis, best_alpha, best_beta, min_error = P.update(epoch)

data = {
        "Alpha" : alpha_lis,
        "Beta"  : beta_lis,
        "Error"      : error_lis,
        "Best Alpha" : best_alpha,
        "Best Beta"  : best_beta,
        "Min Error"      : min_error
        }

output = pd.DataFrame(data).astype(float)
output.to_csv("Error.csv")

#==============================================================================

_ = solver(best_alpha, best_beta, ISAVE, ITMAX, IMETHOD, ngrd)

#==============================================================================

path = "Result"

if not os.path.isdir(path):
    os.mkdir(path)

files = glob.glob("*.plt")

for file in files:
    shutil.move(str(os.path.join(os.getcwd(),file)),str(os.path.join(os.getcwd(),path,file)))

filecsv = glob.glob("*.csv")[0]
shutil.move(str(os.path.join(os.getcwd(),filecsv)), str(os.path.join(os.getcwd(),path,filecsv)))

#==============================================================================