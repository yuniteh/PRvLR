# README:
# This user/mode relies on saving variables to PVD so that they're retrievable on power cycles.
# On a first boot of the user/mode you may see an error stating 'Error while request to save PCE variables:
# Error returned from PCE Variable file open routine'. This error is produced because the PCEVarList.cfg file
# located in /config/mode/caps/*MODENAME*/ is trying to load variables (e.g. ADAPT_FLAG) which don't currently 
# exist, and therefore cannot be loaded. To resolve this issue, run this script to initialise the variables so
# CAPS can find them. You should only need to run this once. If you're seeing the error often, contact CBM. 
#   1: Open PuTTY and navigate to the mode's SLAVE FOLDER:
#       a) cd /
#       b) cd /config/mode/caps/*MODENAME*/SLAVE/ 
#   2: Type the following:
#       a) python initialise.py
#   3: A message should state that initialise is running, and then a second message stating that it has completed.
#   4: Close PuTTY and power cycle the device. Null variables will now exist in the system and CAPS should operate normally.

import pcepy.pce as pce
import numpy as np

# Number of modes/classes.
numClasses = 5
# Number of modes/classes.
numModels = 5
# Number of EMG channels.
numEMG = 6
# Number of features. (10 for '47')
featNum = 4
# Matrix size.
matSize = numEMG * featNum
# Sample threshold
sampThres = 100
# Number of samples in one frame 
window = 200

print('RUNNING INITIALISE...')

# Common variables
pce.set_var('CTRL',1)
pce.set_var('MODE', -1)
pce.set_var('CUR_VAL', 0)
pce.set_var('THRESH_VAL', 0)
pce.set_var('NOISE', 0)
pce.set_var('NOISE_CH', 0)
pce.set_var('NOISE_SCALE', 0)
pce.set_var('COLLECT', 2)
pce.set_var('IN_TRIAL', 0)
pce.set_var('MVC_R', np.zeros((1, numClasses), dtype=float, order='F'))
pce.set_var('MVC_T', np.zeros((1, numClasses), dtype=float, order='F'))
pce.set_var('MVC', np.zeros((numClasses, numEMG), dtype=float, order='F'))
pce.set_var('N_C', np.zeros((1, numClasses), dtype=float, order='F'))
pce.set_var('N_R', np.zeros((1, numClasses), dtype=float, order='F'))
pce.set_var('N_T', np.zeros((1, numClasses), dtype=float, order='F'))
pce.set_var('DAQ_NOISY', np.zeros((numEMG, window), dtype=float, order='F'))
pce.set_var('FEAT_RAW', np.zeros((1, matSize), dtype=float, order='F'))
pce.set_var('FEAT_NOISY', np.zeros((1, matSize), dtype=float, order='F'))


# Variables for PR only
pce.set_var('NEW_CLASS', 0)
pce.set_var('OLD_CLASS', 0)
pce.set_var('VEL', 0)
pce.set_var('RAMP', 1)
pce.set_var('CLASS_ACTIVE', 0)
pce.set_var('CLASS_EST', -1)
pce.set_var('OUT_MAP', np.zeros((1, numClasses), dtype=float, order='F'))
pce.set_var('WG_ADAPT', np.zeros((matSize, numClasses), dtype=float, order='F'))
pce.set_var('CG_ADAPT', np.zeros((1, numClasses), dtype=float, order='F'))
pce.set_var('MID', np.zeros((numClasses, 1), dtype=float, order='F'))

for i in range(0, numClasses):
    pce.set_var('COV' + str(i), np.zeros((matSize, matSize), dtype=float, order='F'))
    pce.set_var('MN' + str(i), np.zeros((1, matSize), dtype=float, order='F'))
    pce.set_var('CLASS_MAV' + str(i), 0)

# Variables for regression only
pce.set_var('Y_LABEL',-1)
pce.set_var('X', np.zeros((sampThres, matSize + 1), dtype=float, order='F'))
pce.set_var('Y', np.zeros((sampThres, 1), dtype=float, order='F'))
pce.set_var('Y_EST', np.zeros((numModels,1), dtype=float, order='F'))
pce.set_var('ZERO_FLAGS', np.zeros((numModels,numModels), dtype=float, order='F'))
pce.set_var('W', np.zeros((matSize + 1, numModels), dtype=float, order='F'))
for i in range(1, numModels):
    pce.set_var('XX' + str(i), np.zeros((matSize + 1, matSize + 1), dtype=float, order='F'))
    pce.set_var('XY' + str(i), np.zeros((matSize + 1, 1), dtype=float, order='F'))
    
print('INITIALISE COMPLETE')