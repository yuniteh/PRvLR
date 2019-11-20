# Import all the required modules. These are helper functions that will allow us to get variables from CAPS PC
import os
import pcepy.pce as pce
import numpy as np

numEMG = 6
a2d = (2**16 - 1) / 10
eps = .1

def dispose():
    pass

############################################# MAIN FUNCTION LOOP #######################################################
def run():
    #i = 1
    daq = np.array(pce.get_var('DAQ_DATA').to_np_array()[0:numEMG,:], order='F').astype(float)
    ratio = pce.get_var('NOISE_SCALE')/100
    amp = pce.get_var('MVC').to_np_array() 
    dc = 5 * np.ones(daq.shape[1])
    daq = (daq / a2d)  - dc
    noise = ratio * np.random.randn(daq.shape[1])
    ch = int(pce.get_var('NOISE_CH'))
    print ratio
    
    if pce.get_var('NOISE') == 1:
        for i in range(0,numEMG):
            daq[i,:] = daq[i,:] + noise
    elif pce.get_var('NOISE') == 2:
        # SNR dependent noise
        # noise *= np.amax(amp[:,ch])
        # print np.amax(amp[:,ch])
        daq[ch,:] = daq[ch,:] + noise
    elif pce.get_var('NOISE') == 3:
        daq[ch,:] = np.zeros(daq.shape[1]) + eps * np.random.randn(daq.shape[1])
    
    daq = (daq + dc) * a2d
    pce.set_var('DAQ_NOISY', daq.astype(float, order='F'))