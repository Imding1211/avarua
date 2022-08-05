# -*- coding: utf-8 -*-
"""
Created on Sun Apr 24 00:28:52 2022

@author: ding
"""

import glob
import os
import shutil

def outputmakefile(icase,nohup):
  makeos = checkos()
  path = 'makefile'
  f = open(path, 'w')
  f.write("\n")
  f.write("all: $(OBJECTS)\n")
  if makeos == False:
    f.write("\tf2py -c main/solver/CHAN_solver_IT.f90 -m solver_CHAN_IT\n")
    f.write("\tf2py -c main/solver/CHAN_solver_ITMAX.f -m solver_CHAN_ITMAX\n")
    f.write("\tf2py -c main/solver/HAR_solver_3method.f -m solver_HAR\n")
    f.write("\tpython Ding.py\n")
  if icase == 1:
    if nohup == True:
      f.write("\tnohup python main/CHAN_IT_main.py &\n")
    else:
      f.write("\tpython main/CHAN_IT_main.py\n")
  elif icase == 2:
    if nohup == True:
      f.write("\tnohup python main/CHAN_ITMAX_main.py &\n")
    else:
      f.write("\tpython main/CHAN_ITMAX_main.py\n")
  elif icase == 3:
    if nohup == True:
      f.write("\tnohup python main/HAR_main.py &\n")
    else:
      f.write("\tpython main/HAR_main.py\n")
  f.write("\n")
  f.write(".PHONY:clean\n")
  f.write("clear:\n")
  f.write("\trm -f *.exe\n")
  f.write("\trm -f *.mod\n")
  f.write("\trm -f *.o\n")
  f.write("\trm -f *.dat\n")
  f.write("\trm -f *.out\n")
  f.close()

  return

def checkos():
  
  files = glob.glob("main/*.so")

  if len(files) != 3:
    ans = False
  else:
    ans = True

  return ans

def movesolver():

  path = "main"
  filesos = glob.glob("*.so")
  for file in filesos:
      shutil.move(str(os.path.join(os.getcwd(),file)),str(os.path.join(os.getcwd(),path,file)))

  return

if __name__ == '__main__':
    movesolver()
