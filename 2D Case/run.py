# -*- coding: utf-8 -*-
"""
Created on Sun Apr 24 00:28:52 2022

@author: ding
"""

import subprocess
from Ding import outputmakefile

#================參數設定==========================================================

"""
icase 1 = CHAN_IT
icase 2 = CHAN_ITMAX
icase 3 = HAR
"""

icase = 3
nohup = False

#==================================================================================
outputmakefile(icase,nohup)
make_proc = subprocess.call("make")
#==================================================================================
