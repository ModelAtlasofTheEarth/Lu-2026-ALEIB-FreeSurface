import numpy as np
import os
import glob
import sys
import subprocess

test_timeratio = [1, 2, 4, 8, 16, 32, 64] 
for tRatio in test_timeratio:
    subprocess.run([sys.executable, '/home/nl/env/uw3dev/temtests/Ex_2DTopoRelax_Cartesian_EulerainMC_uw3dev_tratio.py', str(tRatio)])
    for use_FSSA in [True,False]:
        subprocess.run([sys.executable, '/home/nl/env/uw3dev/temtests/Ex_2DTopoRelax_Cartesian_ALE_uw3dev_tratio.py', str(tRatio), str(use_FSSA)])
        subprocess.run([sys.executable, '/home/nl/env/uw3dev/temtests/Ex_2DTopoRelax_Cartesian_ALEIB_uw3dev_tratio.py', str(tRatio), str(use_FSSA)])
