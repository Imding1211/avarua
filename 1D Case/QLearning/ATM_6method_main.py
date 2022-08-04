#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 14:00:38 2021

@author: ding
"""

from ATM_6method_agent import agent
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

#==============================================================================
"""
1 MUSCL
2 Thincem
3 Thinc
4 ATM1
5 ATM2
6 ATM3
"""        
epoch  = 5000
method = 4
period = 1
Mesh   = 8

#==============================================================================

def sett(lengh, mesh):
    dx = lengh/(100*mesh)
    indexmax = 100*mesh+1
    axis_x = np.arange(0, lengh, dx)
    return indexmax, axis_x

#==============================================================================

def plot(m, r, u1, u2, n, nname, a1, a2, b1, b2):

    msize=3
    if n != 1:    
        for j in range(n):
            plt.title("Two interacting blast waves problem\n"+nname)
            plt.xlabel('x')
            plt.ylabel('Density')
            plt.plot(x128,r[:],'-',color="black",label='Ref.')
            if m == 1 or m == 2:
                plt.plot(x,u1[j,:],'o',color="red",markerfacecolor="white",markersize=msize,label='Pri.')
                plt.plot(x,u2[j,:],'v',color="gray",markerfacecolor="white",markersize=msize,label='Con.')
            elif m == 3 or m == 4:
                plt.plot(x,u1[j,:],'o',color="red",markerfacecolor="white",markersize=msize,label='Pri.-Beta='+str(round(b1[j],3)))
                plt.plot(x,u2[j,:],'v',color="gray",markerfacecolor="white",markersize=msize,label='Con.-Beta='+str(round(b2[j],3)))
            elif m == 5 or m == 6:
                plt.plot(x,u1[j,:],'o',color="red",markerfacecolor="white",markersize=msize,label='Pri.-Beta='+str(round(b1[j],3)))
                plt.plot(x,u2[j,:],'v',color="gray",markerfacecolor="white",markersize=msize,label='Con.-Beta='+str(round(b2[j],3)))
            elif m == 11 or m == 12:
                plt.plot(x,u1[j,:],'o',color="red",markerfacecolor="white",markersize=msize,label='Pri.-Beta='+str(round(b1[j],3)))
                plt.plot(x,u2[j,:],'v',color="gray",markerfacecolor="white",markersize=msize,label='Con.-Beta='+str(round(b2[j],3)))
            else:
                plt.plot(x,u1[j,:],'o',color="red",markerfacecolor="white",markersize=msize,label='Pri.-Alpha='+str(round(a1[j],1))+' Beta='+str(round(b1[j],3)))
                plt.plot(x,u2[j,:],'v',color="gray",markerfacecolor="white",markersize=msize,label='Con.-Alpha='+str(round(a2[j],1))+' Beta='+str(round(b2[j],3)))
            plt.axis([0,1.0,0,30])
            plt.legend(loc=2)
            plt.savefig(file+'//'+str(j+1)+'-compare.png', dpi=300)
            plt.clf()
    
    plt.title("Two interacting blast waves problem\n"+name[int((method-1)/2)])
    plt.xlabel('x')
    plt.ylabel('Density')
    plt.plot(x128,r[:],'-',color="black",label='Ref.')
    if m == 1 or m == 2:
        plt.plot(x,u1[-1,:],'o',color="red",markerfacecolor="white",markersize=msize,label='Pri.')
        plt.plot(x,u2[-1,:],'v',color="gray",markerfacecolor="white",markersize=msize,label='Con.')
    elif m == 3 or m == 4:
        plt.plot(x,u1[-1,:],'o',color="red",markerfacecolor="white",markersize=msize,label='Pri.-Beta='+str(round(b1[-1],3)))
        plt.plot(x,u2[-1,:],'v',color="gray",markerfacecolor="white",markersize=msize,label='Con.-Beta='+str(round(b2[-1],3)))
    elif m == 5 or m == 6:
        plt.plot(x,u1[-1,:],'o',color="red",markerfacecolor="white",markersize=msize,label='Pri.-Beta='+str(round(b1[-1],3)))
        plt.plot(x,u2[-1,:],'v',color="gray",markerfacecolor="white",markersize=msize,label='Con.-Beta='+str(round(b2[-1],3)))
    elif m == 11 or m == 12:
        plt.plot(x,u1[-1,:],'o',color="red",markerfacecolor="white",markersize=msize,label='Pri.-Beta='+str(round(b1[-1],3)))
        plt.plot(x,u2[-1,:],'v',color="gray",markerfacecolor="white",markersize=msize,label='Con.-Beta='+str(round(b2[-1],3)))
    else:
        plt.plot(x,u1[-1,:],'o',color="red",markerfacecolor="white",markersize=msize,label='Pri.-Alpha='+str(round(a1[-1],1))+' Beta='+str(round(b1[-1],3)))
        plt.plot(x,u2[-1,:],'v',color="gray",markerfacecolor="white",markersize=msize,label='Con.-Alpha='+str(round(a2[-1],1))+' Beta='+str(round(b2[-1],3)))
    plt.legend(loc=2)
    plt.savefig(file+'//'+'compare_paper_size.png',dpi=300)
    plt.clf()
    return 

#==============================================================================

name = ["MUSCL",
        "Thincem",
        "Thinc",
        "ATM1",
        "ATM2",
        "ATM3",]
[imax, x] = sett(1.0, Mesh)
[_, x128] = sett(1.0, 128)
ref = np.loadtxt('12800.csv',delimiter=",",dtype=np.float32,skiprows=1)[:,1]
method = 1+(method-1)*2
data = np.zeros([2,epoch])

#==============================================================================   

aP = agent(1, method , Mesh, period)
aC = agent(1, method+1, Mesh, period)

dataP, uP, alphaP, betaP = aP.update(epoch)
dataC, uC, alphaC, betaC = aC.update(epoch)

file = name[int((method-1)/2)]+"-compare"+str(period)
if not os.path.isdir(file):
    os.mkdir(file)
    
plot(method, ref, uP, uC, period, name[int((method-1)/2)], alphaP, alphaC, betaP, betaC)

data[0,:]=dataP
data[1,:]=dataC

output = pd.DataFrame(data).astype(float)
output.to_csv(name[int((method-1)/2)]+"-Error.csv")
