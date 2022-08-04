#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 27 19:02:16 2021

@author: ding
"""

from PIL import ImageTk
from solver import solver, ini
from functools import partial
import os
import tkinter as tk
import tkinter.ttk as ttk
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

# =============================================================================

root = tk.Tk()
root.title("ATM")
root.geometry("700x820")
root.config(bg="skyblue")
text_font = ('Courier New','16')

global case_name, method_name, savenum
savenum = 0
case_name = ['Two interacting blast waves problem',
             'Detonation Wave Case']
method_name = ['MUSCL_PRIMITIVE',
               'MUSCL_CONSERVATIVE',
               'THINCEM_PRIMITIVE',
               'THINCEM_CONSERVATIVE',
               'THINC_PRIMITIVE',
               'THINC_CONSERVATIVE',
               'ATM_PRIMITIVE',
               'ATM_CONSERVATIVE',
               'ATM2_PRIMITIVE',
               'ATM2_CONSERVATIVE',
               'ATM3_PRIMITIVE',
               'ATM3_CONSERVATIVE',
               'ATM4_PRIMITIVE',
               'ATM4_CONSERVATIVE',
               'ATM5_PRIMITIVE',
               'ATM5_CONSERVATIVE',
               'ATM6_PRIMITIVE',
               'ATM6_CONSERVATIVE',
               'ATM7_PRIMITIVE',
               'ATM7_CONSERVATIVE',
               'ATM8_PRIMITIVE',
               'ATM8_CONSERVATIVE',
               'ATM(2021)_PRIMITIVE',
               'ATM(2021)_CONSERVATIVE',
               'ATM2(2021)_PRIMITIVE',
               'ATM2(2021)_CONSERVATIVE']

# =============================================================================

def sett(lengh, mesh): 
    
    dx = lengh/(100*mesh)
    indexmax = 100*mesh+1
    axis_x = np.arange(0, lengh, dx)
    
    return indexmax, axis_x

# =============================================================================

def show():
    
    my_img = ImageTk.PhotoImage(file='tk.png')
    if (os.path.isfile('tk.png') == True):
        os.remove('tk.png')
    my_label.config(image=my_img)
    my_label.image = my_img

# =============================================================================

def clear():
    
    plt.clf()
    plt.savefig('tk.png')
    my_img = ImageTk.PhotoImage(file='tk.png')
    if (os.path.isfile('tk.png') == True):
        os.remove('tk.png')
    my_label.config(image=my_img)
    my_label.image = my_img
    
# =============================================================================

def getvalue(event):
    getvalue_()

# =============================================================================

def getvalue_():
    
    global u, imax, x, lengh_x
    global mesh, icase, imethod, line, color
    global alpha, beta

    mesh = int(meshentry.get())
    icase = casecombo.current()+1
    line = linecombo.get()
    color = colorcombo.get()
    
    imethod = methodcombo.current()+1
    
    msize = 3
    
    if imethod == 1 or imethod == 2:
    	alpha = 100
    	beta = 1.2
        
    elif imethod == 3 or imethod == 4:
        alpha = 100
        beta = float(betaentry.get())
        
    elif imethod == 5 or imethod == 6:
        alpha = 100
        beta = float(betaentry.get())
                                
    elif imethod == 23 or imethod == 24:
        alpha = 100
        beta = float(betaentry.get())

    elif imethod == 25 or imethod == 26:
        alpha = 100
        beta = float(betaentry.get())
        
    else:
        alpha = float(alphaentry.get())
        beta = float(betaentry.get())
        
    u, step, lengh_x= ini(icase,mesh)
    imax, x = sett(lengh_x, mesh)
    
    progress["value"] = 0
    for i in range(1,step+1):
        u = solver(icase,imethod,mesh,alpha,beta,u)
        progress["value"] += (100/step)
        root.update_idletasks()
    
    plt.title(case_name[icase-1])
    plt.xlabel('x')
    plt.ylabel('Density')
    if imethod == 1 or imethod == 2:
    	plt.plot(x,u[0:imax-1,0],line,color=color,markerfacecolor="white",markersize=msize,label=method_name[imethod-1]) 
    elif imethod == 3 or imethod == 4:
    	plt.plot(x,u[0:imax-1,0],line,color=color,markerfacecolor="white",markersize=msize,label=method_name[imethod-1]+", Beta="+str(beta))
    elif imethod == 5 or imethod == 6:
    	plt.plot(x,u[0:imax-1,0],line,color=color,markerfacecolor="white",markersize=msize,label=method_name[imethod-1]+", Beta="+str(beta))
    elif imethod == 23 or imethod == 24:
    	plt.plot(x,u[0:imax-1,0],line,color=color,markerfacecolor="white",markersize=msize,label=method_name[imethod-1]+", Beta="+str(beta))
    elif imethod == 25 or imethod == 26:
    	plt.plot(x,u[0:imax-1,0],line,color=color,markerfacecolor="white",markersize=msize,label=method_name[imethod-1]+", Beta="+str(beta))
    else:
        plt.plot(x,u[0:imax-1,0],line,color=color,markerfacecolor="white",markersize=msize,label=method_name[imethod-1]+", Alpha="+str(alpha)+" Beta="+str(beta))
    if icase == 1:
        #plt.axis([0,lengh_x,0,7])
        plt.legend(loc=2)
    elif icase == 2:
        plt.legend(loc=1)
        plt.axis([0,lengh_x,0.001,0.003])
        plt.yticks([0.001,0.002,0.003])
    plt.savefig('tk.png')
    show()
    
# =============================================================================  
    
def save():
    
    global savenum
    savenum += 1
    name = str(nameentry.get())
    if name == "":
        name = "save-"+str(savenum)

    plt.savefig(name+'.png',dpi=300)
    
# =============================================================================    

my_label = tk.Label(root, image=None)
my_label.config(bg="skyblue")
my_label.place(relx=0.5, y=30, anchor="n")

# =============================================================================

progress = ttk.Progressbar(root, orient="horizontal", length=600)
progress.place(relx=0.5, y=520, anchor="n")

# =============================================================================

labelcase = tk.Label(root, text="Case",font=text_font)
labelcase.config(bg="skyblue")
labelcase.place(x=60, y=560, anchor="nw")

icasename = ["Two interacting blast waves problem",
             "Detonation Wave Case"]

casecombo = ttk.Combobox(root,width=45,font=text_font)
casecombo["values"] = icasename
casecombo.current(0)
casecombo.place(x=150, y=560, anchor="nw")

# =============================================================================

labelmethod = tk.Label(root, text="Method",font=text_font)
labelmethod.config(bg="skyblue")
labelmethod.place(x=60, y=600, anchor="nw")

imethodname =  ["MUSCL PRIMITIVE",
                "MUSCL CONSERVATIVE",
                "THINCEM PRIMITIVE",
                "THINCEM CONSERVATIVE",
                "THINC PRIMITIVE",
                "THINC CONSERVATIVE",
                "ATM PRIMITIVE",
                "ATM CONSERVATIVE",
                "ATM2 PRIMITIVE",
                "ATM2 CONSERVATIVE",
                "ATM3 PRIMITIVE",
                "ATM3 CONSERVATIVE",
                "ATM4 PRIMITIVE",
                "ATM4 CONSERVATIVE",
                "ATM5 PRIMITIVE",
                "ATM5 CONSERVATIVE",
                "ATM6 PRIMITIVE",
                "ATM6 CONSERVATIVE",
                "ATM7 PRIMITIVE",
                "ATM7 CONSERVATIVE",
                "ATM8 PRIMITIVE",
                "ATM8 CONSERVATIVE",
                "ATM(2021) PRIMITIVE",
                "ATM(2021) CONSERVATIVE",
                "ATM2(2021) PRIMITIVE",
                "ATM2(2021) CONSERVATIVE"]

methodcombo = ttk.Combobox(root,width=45,font=text_font)
methodcombo["values"] = imethodname
methodcombo.current(0)
methodcombo.place(x=150, y=600, anchor="nw")

# =============================================================================

labelmesh = tk.Label(root, text="Mesh",font=text_font)
labelmesh.config(bg="skyblue")
labelmesh.place(x=60, y=640, anchor="nw")

meshentry = tk.Entry(font=text_font,width=10)
meshentry.place(x=150, y=640, anchor="nw")

# =============================================================================

labelline = tk.Label(root, text="Line",font=text_font)
labelline.config(bg="skyblue")
labelline.place(x=300, y=640, anchor="nw")

linename = ["-",
             "--",
             "-.",
             "o",
             "v"]

linecombo = ttk.Combobox(root,width=10,font=text_font)
linecombo["values"] = linename
linecombo.current(0)
linecombo.place(x=390, y=640, anchor="nw")

# =============================================================================

labelalpha = tk.Label(root, text="Alpha",font=text_font)
labelalpha.config(bg="skyblue")
labelalpha.place(x=60, y=680, anchor="nw")

alphaentry = tk.Entry(font=text_font,width=10)
alphaentry.place(x=150, y=680, anchor="nw")

# =============================================================================

labelcolor = tk.Label(root, text="Color",font=text_font)
labelcolor.config(bg="skyblue")
labelcolor.place(x=300, y=680, anchor="nw")

colorname = ["black",
             "blue",
             "orange",
             "green",
             "red",
             "purple",
             "brown",
             "pink",
             "gray",
             "olive",
             "cyan"]

colorcombo = ttk.Combobox(root,width=10,font=text_font)
colorcombo["values"] = colorname
colorcombo.current(0)
colorcombo.place(x=390, y=680, anchor="nw")

# =============================================================================

labelbeta = tk.Label(root, text="Beta",font=text_font)
labelbeta.config(bg="skyblue")
labelbeta.place(x=60, y=720, anchor="nw")

betaentry = tk.Entry(font=text_font,width=10)
betaentry.place(x=150, y=720, anchor="nw")

# =============================================================================

labelname = tk.Label(root,text="Name",font=text_font)
labelname.config(bg="skyblue")
labelname.place(x=300, y=720, anchor="nw")

nameentry = tk.Entry(font=text_font,width=11)
nameentry.place(x=390, y=720, anchor="nw")

# =============================================================================

Select_button = tk.Button(root, text="Select", command=getvalue_)
Select_button.place(relx=0.25, y=760, anchor="n")

# =============================================================================

save_button = tk.Button(root, text="Save", command=save)
save_button.place(relx=0.75, y=760, anchor="n")

# =============================================================================

clear_button = tk.Button(root, text="clear", command=clear)
clear_button.place(relx=0.5, y=760, anchor="n")

# =============================================================================

root.bind("<Return>", partial(getvalue))
root.mainloop()
