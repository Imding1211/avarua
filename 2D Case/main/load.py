#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 24 13:35:55 2022

@author: ding
"""

def loadsetting(path,head):

    f = open(path)
    lines = []
    for line in f.readlines():
        lines.append(line)
        
    f.close
    par = lines[head:]
    
    return par
