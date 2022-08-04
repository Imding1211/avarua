#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 13:27:26 2021

@author: ding
"""

from solver import solver, ini
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import numpy as np
import math

# =============================================================================
# Parameters
# =============================================================================
"""
1 MUSCL
2 Thincem
3 Thinc
4 ATM
5 ATM2
6 ATM3
7 ATM4
8 ATM5
9 ATM6
10 ATM7
11 ATM8
12 ATM(2021)
13 ATM(2021)2
"""  
icase               = 1
imethod             = [1,2,3,4,5,6,7,8,9,10,11,12,13]
mesh                = [4]
alpha               = [500]
beta                = [1.25]
#beta = np.linspace(1.2, 1.6, 5)

# =============================================================================
# Function
# =============================================================================

def sett(lengh, mesh): 
    
    dx = lengh/(100*mesh)
    indexmax = 100*mesh+1
    axis_x = np.arange(0, lengh, dx)
    
    return indexmax, axis_x

# =============================================================================

def plot1(n, u64, u8, imax64, imax8, x64, x8, a, b, m): 
    
    #Error = math.sqrt(sum((u64[0:12800:16]-u8[0:800,0])**2)/800)
    max64 = np.argmax(u64[:])
    max8 = np.argmax(u8[:,0])
    Error = 0
    for nn in [-1,0,1]:
        Error = Error + (u64[max64+nn]-u8[max8+nn,0])**2
    Error = math.sqrt(Error/3)
    #plt.title(case_name[icase-1]+"\nRMSE="+str(round(Error,3)))
    plt.xlabel('x')
    plt.ylabel('Density')
    plt.plot(x64,u64[0:imax64-1],'-',color='black',label="Ref.")
    if n == 1 or n == 2:
        plt.plot(x8,u8[0:imax8-1,0],'o',color='red',markersize=3,markerfacecolor="white",label="MUSCL-"+str(m*100))
    elif n == 3 or n == 4:
        plt.plot(x8,u8[0:imax8-1,0],'o',color='red',markersize=3,markerfacecolor="white",label="Thincem Beta="+str(b))
    else:
        plt.plot(x8,u8[0:imax8-1,0],'o',color='red',markersize=3,markerfacecolor="white",label="Alpha="+str(a)+" Beta="+str(b)+"-"+str(m*100))
    #plt.plot(x8,u81[0:imax8-1,0],'-.',color='red',label="Thincem Beta="+str(b))
    #plt.plot(x8,u82[0:imax8-1,0],'-',color='gray',label="Thinc Beta="+str(b))
    plt.axis([0,1,0,7])
    plt.legend(loc=2)
    plt.savefig(str(n)+"-Alpha="+str(a)+"-Beta="+str(b)+'-MESH'+str(m)+'.png',dpi=300)
    plt.clf()
    
# =============================================================================

def plot2(n, u64, u81, u82, imax64, imax8, x64, x8, a, b, m): 
    
    plt.title(case_name[icase-1])
    plt.xlabel('x')
    plt.ylabel('Density')
    plt.plot(x64,u64[0:imax64-1,0],'-',color='black',label="Ref.")
    #plt.plot(x8,u8[0:imax8-1,0],'-.',color='red',label="Alpha="+str(a)+" Beta="+str(b))
    #plt.plot(x8,u81[0:imax8-1,0],'-.',color='red',label="Thincem Beta="+str(b))
    #plt.plot(x8,u82[0:imax8-1,0],'-',color='gray',label="Thinc Beta="+str(b))
    plt.axis([0,1,0,7])
    plt.legend(loc=2)
    plt.savefig(str(n)+"-Alpha="+str(a)+"-Beta="+str(b)+'-MESH'+str(m)+'.png',dpi=300)
    plt.clf()

# =============================================================================
# Setting
# =============================================================================

case_name = ['Two interacting blast waves problem','Detonation Wave Case']
method_name = ['MUSCL_PRIMITIVE','MUSCL_CONSERVATIVE','THINCEM_PRIMITIVE','THINCEM_CONSERVATIVE','HMT_PRIMITIVE','HMT_CONSERVATIVE','HMT(2021)_PRIMITIVE','HMT(2021)_CONSERVATIVE']
imethod = 1+(np.array(imethod)-1)*2

# =============================================================================
# Ref
# =============================================================================

u64, step64, lengh_x= ini(icase,128)

imax64, x64 = sett(lengh_x,128)

u64 = np.loadtxt('12800.csv',delimiter=",",dtype=np.float32,skiprows=1)[:,1]

for method in imethod:
    
# =============================================================================
# 58
# =============================================================================

    for m in mesh:
        imax5, x5 = sett(1, m)
        for a in alpha:
            for b in beta:
                u5, step5, _= ini(icase,m)
                for i in range(1,step5+1):
                    print (f'\rMesh={m}, Alpha={a}, Beta={b}, Step [{i}/{step5}]',end='')
                    u5 = solver(icase,method,m,a,b,u5)
                plot1(method, u64, u5, imax64, imax5, x64, x5, a, b, m)
    print(" ")
    print(" Done")

# =============================================================================
# 68
# =============================================================================

    for m in mesh:
        imax6, x6 = sett(1, m)
        for a in alpha:
            for b in beta:
                u6, step6, _= ini(icase,m)
                for i in range(1,step6+1):
                    print (f'\rMesh={m}, Alpha={a}, Beta={b}, Step [{i}/{step6}]',end='')
                    u6 = solver(icase,method+1,m,a,b,u6)
                plot1(method+1, u64, u6, imax64, imax6, x64, x6, a, b, m)
    print(" ")
    print(" Done")
    
