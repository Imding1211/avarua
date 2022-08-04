#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 14:13:15 2021

@author: ding
"""

from ATM_6method_env import AI, ini
from RL_brain import QLearningTable
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

class agent:
    
    def __init__(self, icase, imethod, mesh, state):
        self.nn = ["Pri.",
                   "Con."]
        self.case_name = ['Two interacting blast waves problem',
                          'Detonation Wave Case']
        self.method_name = ['MUSCL_Pri.',
                            'MUSCL_Con.',
                            'THINCEM_Pri.',
                            'THINCEM_Con.',
                            'THINC_Pri.',
                            'THINC_Con.',
                            'ATM1_Pri.',
                            'ATM1_Con.',
                            'ATM2_Pri.',
                            'ATM2_Con.',
                            'ATM3_Pri.',
                            'ATM3_Con.',]
        self.line_style = ["-","-."]
        self.icase = icase
        self.imethod = imethod
        self.mesh = mesh
        self.state = state
        self.env = AI(self.icase, self.mesh, self.imethod)
        self.RL_alpha = QLearningTable(actions=list(range(self.env.n_actions_alpha)))
        self.RL_beta = QLearningTable(actions=list(range(self.env.n_actions_beta)))
        _, self.nT, self.lengh_x = ini(self.icase, self.mesh)
        self.ref = np.loadtxt('12800.csv',delimiter=",",dtype=np.float32,skiprows=1)[:,1]
        self.imax, self.x = self.sett(self.lengh_x, self.mesh)
        _, self.x128 = self.sett(self.lengh_x, 128)

#==============================================================================

    def sett(self, lengh, mesh):
        
        dx = lengh/(100*mesh)
        indexmax = 100*mesh+1
        axis_x = np.arange(0, lengh, dx)
    
        return indexmax, axis_x
        
#==============================================================================
        
    def update(self, epoch):
        
        best_reward = 0.5
        all_reward = []
        best_alpha_lis = []
        best_beta_lis = []
        best_h = np.zeros([self.state,(self.mesh*100)])
        
        for episode in range(epoch):
            
            T_reward = []
            alpha_lis = []
            beta_lis = []
            temp_h = np.zeros([self.state,(self.mesh*100)])   
            observation, nt, _ = self.env.reset()
            
            while True:
            
                if nt % int(self.nT/self.state) == 0:
                    a_alpha = self.RL_alpha.choose_action(str(observation))
                    a_beta = self.RL_beta.choose_action(str(observation))
    
                observation_, reward, done, nt, h, alpha, beta = self.env.step(a_alpha, a_beta)
                
                print (f'\rEpisode = {episode+1}, step = [{nt+1}/{self.nT}]',end='')
                
                T_reward.append(reward)
                
                self.RL_alpha.learn(str(observation), a_alpha, reward, str(observation_))
                self.RL_beta.learn(str(observation), a_beta, reward, str(observation_))
                
                if nt % int(self.nT/self.state) == 0:
                    observation = observation_
                    temp_h[int(nt/int(self.nT/self.state))-1,:] = h[:]
                    alpha_lis.append(alpha)
                    beta_lis.append(beta)
                    
                if done:
                    
                    temp_h[-1,:] = h[:]
                    alpha_lis.append(alpha)
                    beta_lis.append(beta)
                    
                    if best_reward > -T_reward[-1]:
                        best_h[:,:] = temp_h[:,:]
                        best_alpha_lis = alpha_lis
                        best_beta_lis = beta_lis
                        best_reward = -T_reward[-1]
                        all_reward.append(-T_reward[-1])
                        self.plot(temp_h, alpha_lis, beta_lis, all_reward)
                        
                    else:
                        all_reward.append(best_reward)
                        
                    print("")
                    print (f'Error = {round(-T_reward[-1],3)}')
                    break
        print ('Done')        
        return all_reward, best_h, best_alpha_lis, best_beta_lis

#==============================================================================
        
    def plot(self, u, alpha, beta, Error):
        
        msize=3
        file = self.method_name[self.imethod-1]+"_"+str(self.state)
        if not os.path.isdir(file):
            os.mkdir(file)
        if self.state != 1:
            for j in range(self.state):
                plt.title(self.case_name[self.icase-1]+"\nRMSE="+str(round(Error[j],3)))
                
                if self.icase == 1 :
                    plt.axis([0,self.lengh_x,0,30])
                else:
                    plt.axis([0,self.lengh_x,0.001,0.003])
                    plt.yticks([0.001,0.002,0.003])
                    
                plt.xlabel('x')
                plt.ylabel('Density')
                plt.plot(self.x128,self.ref[:],'-',color="black",label='Ref.')
                if self.imethod == 1 or self.imethod == 2:
                    plt.plot(self.x,u[j,:],"o",color="red",markersize=msize,markerfacecolor="white",label=self.method_name[self.imethod-1])
                elif self.imethod == 3 or self.imethod == 4:
                	plt.plot(self.x,u[j,:],"o",color="red",markersize=msize,markerfacecolor="white",label=self.method_name[self.imethod-1]+" Beta="+str(round(beta[j],3)))
                elif self.imethod == 5 or self.imethod == 6:
                    plt.plot(self.x,u[j,:],"o",color="red",markersize=msize,markerfacecolor="white",label=self.method_name[self.imethod-1]+" Beta="+str(round(beta[j],3)))
                elif self.imethod == 11 or self.imethod == 12:
                	plt.plot(self.x,u[j,:],"o",color="red",markersize=msize,markerfacecolor="white",label=self.method_name[self.imethod-1]+" Beta="+str(round(beta[j],3)))
                else:
                    plt.plot(self.x,u[j,:],"o",color="red",markersize=msize,markerfacecolor="white",label=self.method_name[self.imethod-1]+" Alpha="+str(round(alpha[j],1))+" Beta="+str(round(beta[j],3)))
                plt.legend(loc=2)
                plt.savefig(file+'//'+str(j+1)+'-f1_result_'+str(self.state)+'.png',dpi=300)
                plt.clf()
        
        self.plot_paper_size(u, alpha, beta, file, Error)

#==============================================================================

    def plot_paper_size(self, u, alpha, beta, file, Error):
        
        msize=3
        plt.title(self.case_name[self.icase-1]+"\nRMSE="+str(round(Error[-1],3)))
        plt.xlabel('x')
        plt.ylabel('Density')
        plt.plot(self.x128,self.ref[:],'-',color="black",label='Ref.')
        if self.imethod == 1 or self.imethod == 2:
            plt.plot(self.x,u[-1,:],"o",color="red",markersize=msize,markerfacecolor="white",label=self.method_name[self.imethod-1])
        elif self.imethod == 3 or self.imethod == 4:
            plt.plot(self.x,u[-1,:],"o",color="red",markersize=msize,markerfacecolor="white",label=self.method_name[self.imethod-1]+" Beta="+str(round(beta[-1],3)))
        elif self.imethod == 5 or self.imethod == 6:
            plt.plot(self.x,u[-1,:],"o",color="red",markersize=msize,markerfacecolor="white",label=self.method_name[self.imethod-1]+" Beta="+str(round(beta[-1],3)))
        elif self.imethod == 11 or self.imethod == 12:
            plt.plot(self.x,u[-1,:],"o",color="red",markersize=msize,markerfacecolor="white",label=self.method_name[self.imethod-1]+" Beta="+str(round(beta[-1],3)))
        else:
            plt.plot(self.x,u[-1,:],"o",color="red",markersize=msize,markerfacecolor="white",label=self.method_name[self.imethod-1]+" Alpha="+str(round(alpha[-1],1))+" Beta="+str(round(beta[-1],3)))
        plt.legend(loc=2)
        plt.savefig(file+'//'+'f1_result_'+str(self.state)+'_paper_size.png',dpi=300)
        plt.clf()
