#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 15:46:32 2021

@author: ding
"""
import pandas as pd
from solver import solver, ini

u64, step64, lengh_x= ini(1,128)

for i in range(1,step64+1):

    print (f'\rStep [{i}/{step64}]',end='')
    u64 = solver(1,1,128,100,1.2,u64)

print(" ")
print(" Done")

u1 = u64[0:12800,0]
u2 = pd.DataFrame(u1).astype(float)
u2.to_csv("12800.csv")
