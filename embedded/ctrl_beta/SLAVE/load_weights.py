#######################################################################################################################
# Function    : load_weights(args)
# args        : folder directory: (../DATA/foldername), controller type: 'LR' or 'PR'
# Description : This function loads previously calculated weights for the LDA/LR controllers
#######################################################################################################################
# Import all the required modules. These are helper functions that will allow us to get variables from CAPS PC
import sys
import pcepy.pce as pce
import numpy as np

if sys.argv[2] == 'LR':
    print 'Loading linear regression weights...'
    w = np.genfromtxt(str(sys.argv[1]) + 'w.csv', delimiter=',')
    pce.set_var('W', w.astype(float, order='F'))
else:
    numClasses = 5
    out_map = pce.get_var('OUT_MAP').to_np_array()
    print 'Loading LDA weights...'
    cg = pce.get_var('CG_ADAPT').to_np_array()
    wg = np.genfromtxt(str(sys.argv[1]) + 'wg.csv', delimiter=',')
    cg_temp = np.genfromtxt(str(sys.argv[1]) + 'cg.csv', delimiter=',')
    for i in range(0, numClasses):
        cg[0,i] = cg_temp[i]
    
    mid = np.genfromtxt(str(sys.argv[1]) + 'mid.csv', delimiter=',')
    NR = 3 * np.ones((1, numClasses), dtype=float, order='F')
    pce.set_var('N_R', NR)
    # Create vector with just the values of classes trained (for remapping purposes).
    classList = np.nonzero(NR)[1]
    # Update out_map.
    out_map[0,0:len(classList)] = classList
    pce.set_var('OUT_MAP', out_map)
    pce.set_var('WG_ADAPT',wg.astype(float, order='F'))
    pce.set_var('CG_ADAPT',cg.astype(float, order='F'))
    pce.set_var('MID',mid.astype(float, order='F'))

mvc = np.genfromtxt(str(sys.argv[1]) + 'mvc.csv', delimiter=',')
pce.set_var('MVC',mvc.astype(float, order='F'))
print 'COMPLETE'