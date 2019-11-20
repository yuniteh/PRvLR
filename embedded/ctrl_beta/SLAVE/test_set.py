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

def run():
    t = pce.get_var('THRESH_VAL')
    print t
    t += 1
    pce.set_var('THRESH_VAL', t)
