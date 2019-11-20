################################################## relaxed_mode ########################################################
################################################### Version 1.0.4 #######################################################
# FILE            : ctrl_beta.py
# VERSION         : 2.0.0
# FUNCTION        : Adapt upper limb data using an LDA classifier and linear regression model. To be used with Unity project.
# DEPENDENCIES    : None
# SLAVE STEP      : Add After Classify
__author__ = 'lhargrove & rwoodward & yteh'
########################################################################################################################

# Import all the required modules. These are helper functions that will allow us to get variables from CAPS PC
import os
import csv
import pcepy.pce as pce
import pcepy.feat as feat
import numpy as np

# Class dictionary
classmap = {0: 'NO MOVEMENT',
            1: 'WRIST SUP.',
            2: 'WRIST PRO.',
            3: 'WRIST FLEX.',
            4: 'WRIST EXT.'}

# Number of classes.
numClasses = 5
# Number of LR models (+1 for empty rest model)
numModels = 5
# Number of EMG channels.
numEMG = 6
# Feature value ('15'; time domain)
featVal = 15
# Number of features. (4 for '15')
featNum = 4
# Matrix size.
matSize = numEMG * featNum
# Threshold multiplier
thresX = 1.2
# Sample threshold
sampThres = 100
# Number of samples in one frame 
window = 200
# Set max ramp length
ramp_max = 10

def dispose():
    pass

############################################# MAIN FUNCTION LOOP #######################################################
def run():
    # Try/catch to see if data already exists.
    try:
        pce.get_var('MODE')
    except:
        # Initialise all variables.
        #initialiseVariables()
        print('no mode')
        
    # Don't do anything if PCE is training.       
    if pce.get_var('TRAIN_STATUS') != 1:
        ctrl = int(pce.get_var('CTRL'))
        mode = int(pce.get_var('MODE'))
        collect = pce.get_var('COLLECT')
        N_T = pce.get_var('N_T').to_np_array()
        N_R = pce.get_var('N_R').to_np_array()
        chan_mav = pce.get_var('CHAN_MAV').to_np_array()[0:numEMG]
        cur_val = np.average(chan_mav)
        pce.set_var('CUR_VAL', cur_val)
        
        # IF MODE == RESET
        if mode == 999:
            # All variables will be initialised back to 0 (or their initial values).
            initialiseVariables()
        elif mode == 888:
            print 'ERROR'
        # IF MODE == START TEST/SAVE PARAMETERS
        elif mode == 500:
            datafolder = 'DATA/' + pce.get_var('DAQ_OUT_FNAME')
            datadir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', datafolder))
            daqdir = datadir + '/DATA/DAQ'
            pvddir = datadir + '/DATA/PVD'
            if not os.path.exists(datadir):
                os.makedirs(daqdir)
                os.makedirs(pvddir)

            saveWeights(datadir)
            f = open(datadir + "/info.txt","w+")
            f.write(pce.get_var('DAQ_OUT_FNAME'))
            f.write('\nctrl = {0:d}; noise = {1:d}; scale = {2:d}'.format(ctrl, int(pce.get_var('NOISE')), int(pce.get_var('NOISE_SCALE'))))
            f.close()
            
            pce.set_var('DAQ_OUT_FNAME', daqdir + '/data.DAQ')
            pce.set_var('PCE_VAR_OUT_FOLDER', pvddir)
            pce.set_var('MODE', -1)
        elif mode == 501:
            # Save weights
            saveWeights(datadir)
            # Update mode
            pce.set_var('MODE', -1)
            
        # UPDATE COLLECT FLAG
        if -5 < collect <= 0:
            collect -= 1
        else:
            collect = 2
            
        # Extract features from raw data
        feat_noisy = feat.extract(featVal, np.array(pce.get_var('DAQ_NOISY').to_np_array()[0:numEMG,:], dtype='uint16', order='F'))
        feat_raw = feat.extract(featVal, np.array(pce.get_var('DAQ_DATA').to_np_array()[0:numEMG,:], order='F'))
        pce.set_var('FEAT_RAW', feat_raw.astype(float,order='F'))
        pce.set_var('FEAT_NOISY', feat_noisy.astype(float,order='F'))
        
        if pce.get_var('NOISE') != 0:
            mvc = pce.get_var('MVC').to_np_array() 
            # if not np.all(mvc):
                # pce.set_var('MODE', 888)
                # print 'COLLECT MVC!!'
            # else:
            feat_data = feat_noisy
        else:
            feat_data = feat_raw

        
        # IF CONTROL == PATTERN REC
        if ctrl == 1:
            print 'PATTERN REC'
        # IF MODE == INITIAL CLASS TRAINER ACTIVATED
            if (pce.get_var('CLASS_ACTIVE') == 0) & (0 <= mode < numClasses):
                # Reset the temp training counter for the specific class.
                N_T[0, mode] = 0
                pce.set_var('N_T', N_T)
                # Toggle the class_active variable to 1.
                pce.set_var('CLASS_ACTIVE', 1)
            # IF MODE == TRAINING
            elif (mode == 0):
                    print('COLLECTING ' + classmap[mode])
                    # Prepare data for the 'no movement' class. Create means and covariance matrices.
                    collect = collectLDA(mode, cur_val, feat_data, 'THRESH_VAL')
            elif (1 <= mode < numClasses) & (N_R[0, 0] >= 1):
                # Compare the current channel MAV against the no movement threshold, if it exceeds then continue.
                if cur_val > (thresX * pce.get_var('THRESH_VAL')):
                    print('COLLECTING ' + classmap[mode])
                    # Prepare data for any movement other than 'no movement'. Create means and covariance matricies.
                    collect = collectLDA(mode, cur_val, feat_data, ('CLASS_MAV' + str(mode)))
               
            # CLASSIFY AND FORWARD PASS
            # To classify the mode must be 0, 'no motion' must be trained (N_R[0,0] >= 1), 
            # and at least one other class needs to be trained (np.sum(N_R[1:numClasses]]) > 1).
            elif (mode == -1) & (N_R[0, 0] >= 1.0):
                out_map = pce.get_var('OUT_MAP').to_np_array()
                # Check that there is a new class to train.
                if pce.get_var('NEW_CLASS') == 1:
                    print('TRAINING')
                    # Create vector with just the values of classes trained (for remapping purposes).
                    classList = np.nonzero(N_R)[1]
                    # Update out_map.
                    out_map[0,0:len(classList)] = classList
                    pce.set_var('OUT_MAP', out_map)
                    # If channel data is poor, the LDA will fail to classify and will throw a singular matrix error. Catch this error.
                    try:
                        # Train using an LDA classifier.
                        (wg_adapt, cg_adapt) = makeLDAClassifier(classList)   
                        # Add weights to WG and CG arrays and set to PCE.
                        updateWgAndCg(wg_adapt, cg_adapt, classList)
                        # Toggle new_class parameter.
                        pce.set_var('NEW_CLASS', 0)
                    except: 
                        print('ERROR: Bad pooled covariance data resulting in singular matrix.')
                # Get weights. Remove non-trained columns.
                wg_adapt = pce.get_var('WG_ADAPT').to_np_array()[:, out_map.tolist()[0]]
                cg_adapt = pce.get_var('CG_ADAPT').to_np_array()[0, out_map.tolist()[0]]
                # Add noise depending on testing phase.
                #if pce.get_var('NOISE') == 1:
                    #feat_data = feat.extract(featVal, np.array(pce.get_var('DAQ_NOISY').to_np_array()[0:numEMG,:], dtype=np.uint16, order='F'))
                # Forward pass to get estimate.
                lda_out = (np.dot(feat_data, wg_adapt) + cg_adapt)
                # Take argmax and remap for class value (i.e. 0 for 'no movement')
                class_dec =  float(out_map[0, (np.argmax(lda_out))])
                # Set estimate to PCE.
                pce.set_var('CLASS_EST', class_dec)
                # Set velocity out.
                velocity = cur_val/(2 * pce.get_var('MID').to_np_array()[class_dec,0])
                
                if pce.get_var('OLD_CLASS') != class_dec:
                    pce.set_var('RAMP', 1)
                    pce.set_var('OLD_CLASS', class_dec)
                # Ramp proportional control
                ramp = pce.get_var('RAMP')
                velocity *= ramp/ramp_max
                pce.set_var('VEL', velocity)
                if ramp != ramp_max:
                    pce.set_var('RAMP', ramp + 1)
                    
                # Print message with class estimation
                print('FORWARD - ' + str(class_dec))
                print('ramp - ' + str(ramp))
        
        # IF CONTROL == REGRESSION
        elif ctrl == 2:
            print 'REGRESSION'
            if (pce.get_var('CLASS_ACTIVE') == 0) & (0 <= mode < numModels):
                # Reset the temp training counter for the specific class.
                N_T[0, mode] = 0
                pce.set_var('N_T', N_T)
                # Toggle the class_active variable to 1.
                pce.set_var('CLASS_ACTIVE', 1)
            elif (0 <= mode < numModels):
                y_label = pce.get_var('Y_LABEL')
                print y_label
                if y_label != -1:
                    print pce.get_var('THRESH_VAL')
                    print("ylabel: ", y_label)
                    print("COLLECTING TRAINING DATA")
                    if y_label == 0:
                        collect = collectLR(mode, feat_data, y_label)
                        # Update threshold
                        updateThresh(cur_val)
                    elif cur_val > (thresX * pce.get_var('THRESH_VAL')):
                        collect = collectLR(mode, feat_data, y_label)                  
            # TRAIN
            elif mode == -100:
                w = pce.get_var('W').to_np_array()
                for i in range(1, numModels):
                    xx = pce.get_var('XX' + str(i)).to_np_array()
                    xy = pce.get_var('XY' + str(i)).to_np_array()
                    w[:,i] = linearModel(xx,xy)
                pce.set_var('W', w.astype(float, order='F'))
                # Set to prediction mode    
                pce.set_var('MODE', -1)
            # PREDICT
            elif mode == -1:
                y_pred = pce.get_var('Y_EST').to_np_array()
                w = pce.get_var('W').to_np_array()
                for i in range(1, numModels):
                    y_pred[i, 0] = predict(feat_data,w[:,i])
                    if y_pred[i, 0] <= 10 or y_pred[i, 0] > 100:
                        y_pred[i, 0] = 0
                    print str(i) + ': ' + str(y_pred[i, 0])
                pce.set_var('Y_EST', y_pred.astype(float, order='F'))
        # IF COLLECTING MVC
        elif ctrl == 3:
            print 'MVC'
            MVC_T = pce.get_var('MVC_T').to_np_array()
            if (pce.get_var('CLASS_ACTIVE') == 0) & (0 <= mode < numClasses):
                # Reset the temp training counter for the specific class.
                MVC_T[0, mode] = 0
                pce.set_var('MVC_T', MVC_T)
                # Toggle the class_active variable to 1.
                pce.set_var('CLASS_ACTIVE', 1)
            else:
                # Collect resting threshold value
                if mode == 0:
                    updateThresh(cur_val)
                    collect = updateMVC(mode)
                elif (numClasses > mode > 0) & (cur_val > (1.5 * pce.get_var('THRESH_VAL'))):
                    mvc = pce.get_var('MVC').to_np_array()
                    MVC_R = pce.get_var('MVC_R').to_np_array()
                    # Update mean MVC for current movement
                    for i in range(0, numEMG):
                        mvc[mode, i] = updateAverage(mvc[mode, i], chan_mav[i], sampThres * MVC_R[0, mode] + MVC_T[0, mode])
                    pce.set_var('MVC', mvc)
                    collect = updateMVC(mode)
        # SET COLLECT FLAG
        print 'mode = ' + str(mode) + '; collect = ' + str(collect)
        pce.set_var('COLLECT', collect)
        # DO NOTHING
        # In the event that no classes have been trained, print a message to the PCE log.
        # This statement is purely for debugging/logging purposes. It can be removed if necessary.
        # else:
            # print('NO ACTION')
      
#######################################################################################################################
# Function    : initialiseVariables(args)
# args        : None.
# Description : This function is used to initialise variables when starting for the first time, or when resetting.
#######################################################################################################################
def initialiseVariables():
    print('RESET')
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
        
#######################################################################################################################
# Function    : saveWeights(args)
# args        : 
# Description : Save weights to csv file.
#######################################################################################################################
def saveWeights(datadir):
    np.savetxt(datadir + '/' + 'mid.csv', pce.get_var('MID').to_np_array(), fmt='%.20f', delimiter=',')
    np.savetxt(datadir + '/' + 'mvc.csv', pce.get_var('MVC').to_np_array(), fmt='%.20f', delimiter=',')
    np.savetxt(datadir + '/' + 'w.csv', pce.get_var('W').to_np_array(), fmt='%.20f', delimiter=',')
    np.savetxt(datadir + '/' + 'wg.csv', pce.get_var('WG_ADAPT').to_np_array(), fmt='%.20f', delimiter=',')
    np.savetxt(datadir + '/' + 'cg.csv', pce.get_var('CG_ADAPT').to_np_array(), fmt='%.20f', delimiter=',')
        
#######################################################################################################################
# Function    : collectLDA(args)
# args        : mode, transmitted value for identification: cur_val, the current data: feat_data, the current feature data: update_var, argument to update.
# Description : This function is used to remove redundant repetition of code between no-movement and all other classes.
#######################################################################################################################
def collectLDA(mode, cur_val, feat_data, update_var):
    # indicate collecting training data
    collect = 1
    # cov and mean are LDA variables.
    cov_C = pce.get_var('COV' + str(mode)).to_np_array()
    mean_C = pce.get_var('MN' + str(mode)).to_np_array()
    # N_C: Total number of windows used for training. This will increment to Inf.
    N_C = pce.get_var('N_C').to_np_array()
    # N_R: Number of training repetitions.
    N_R = pce.get_var('N_R').to_np_array()
    # N_T: Number of windows used for training on the current repetition. 
    N_T = pce.get_var('N_T').to_np_array()

    # Update the running average of training windows.
    update_val = updateAverage(pce.get_var(update_var), cur_val, N_C[0, mode])
    pce.set_var(update_var, update_val)

    # Update the cov and mean for LDA classification.
    (mean_C, cov_C, N_C[0, mode]) = updateMeanAndCov(mean_C, cov_C, N_C[0, mode], feat_data)

    # Update cov, mean, and total
    pce.set_var('COV' + str(mode), cov_C)
    pce.set_var('MN' + str(mode), mean_C)
    pce.set_var('N_C', N_C)
    
    # Increment and set the temp training counter separately.
    N_T[0, mode] = N_T[0, mode] + 1
    pce.set_var('N_T', N_T)
    print(N_T[0, mode])
    
    # Update mean MVC for current movement
    mid = pce.get_var('MID').to_np_array()
    mid[mode, 0] = updateAverage(mid[mode, 0], cur_val, N_T[0, mode])
    pce.set_var('MID', mid)

    # Once the tmp training counter reaches the threshold, stop collecting data for training (As the variable starts at 0 we have to subtract 1 from sampThres).
    if N_T[0, mode] == (sampThres - 1):
        # Update collect flag
        collect = 0
        # Increment the repetition variable to indicate a new training session has been completed.
        N_R[0, mode] += 1        
        pce.set_var('N_R', N_R)
        # Set new_class to 1. This will indicate that a new training session is ready to be trained.
        pce.set_var('NEW_CLASS', 1)
        # Toggle the class_activate variable to 0.
        pce.set_var('CLASS_ACTIVE', 0)
        # Set the mode back to its standby value of -1.
        pce.set_var('MODE', -1)
    
    return collect
        
#######################################################################################################################
# Function    : updateWgAndCg(args)
# args        : wg_adapt, adapted wg weights: cg_adapt, adapted cg weights: classList, list of classes trained.
# Description : This function iteratively updates wg and cg matrices.
#######################################################################################################################            
def updateWgAndCg(wg_adapt, cg_adapt, classList):
    tmp_wg = pce.get_var('WG_ADAPT').to_np_array()
    tmp_cg = pce.get_var('CG_ADAPT').to_np_array()
    for idx, i in enumerate(classList):
        tmp_wg[:, classList[idx]] = wg_adapt[:, idx]
        tmp_cg[0, classList[idx]] = cg_adapt[0, idx]
    pce.set_var('WG_ADAPT', tmp_wg)
    pce.set_var('CG_ADAPT', tmp_cg)

#######################################################################################################################
# Function    : updateMeanAndCov(args)
# args        : meanMat, the previous mean: covMat: the previous covariance: N: the number of points, cur_feat: the current feature vector
# Description : This function iteratively updates means and covariance matrix based on a new feature point.
#######################################################################################################################
def updateMeanAndCov(meanMat, covMat, N, cur_feat):
    ALPHA = N / (N + 1)
    zero_mean_feats_old = cur_feat - meanMat                                    # De-mean based on old mean value
    mean_feats = ALPHA * meanMat + (1 - ALPHA) * cur_feat                       # Update the mean vector
    zero_mean_feats_new = cur_feat - mean_feats                                 # De-mean based on the updated mean value
    point_cov = np.dot(zero_mean_feats_old.transpose(), zero_mean_feats_new)
    point_cov = np.array(point_cov, np.float64, order='F')
    mean_feats = np.array(mean_feats, np.float64, order='F')
    cov_updated = ALPHA * covMat + (1 - ALPHA) * point_cov                      # Update the covariance
    N = N + 1
    
    return (mean_feats, cov_updated, N)

#######################################################################################################################
# Function    : updateMean(args)
# args        : meanMat, the previous mean: N: the number of points, cur_feat: the current feature vector
# Description : This function iteratively updates means based on a new feature point.
#######################################################################################################################
def updateMean(meanMat, N, cur_feat):
    ALPHA = N/(N+1)
    mean_feats = ALPHA * meanMat + (1 - ALPHA) * cur_feat                       # Update the mean vector
    mean_feats = np.array(mean_feats, np.float64,order='F')
    N = N + 1
    
    return (mean_feats, N)

#######################################################################################################################
# Function    : updateAverage(args)
# args        : prevVal, the previous average: N, the number of points: cur_val, the current value
# Description : This function is used to update averages
#######################################################################################################################    
def updateAverage(prevVal, curVal, N):
    ALPHA = N / (N + 1)
    newVal = ALPHA * prevVal + (1 - ALPHA) * curVal
    
    return newVal

#######################################################################################################################
# Function    : makeLDAClassifier(args)
# args        : class_list, the list of class labels in the classifier
# Description : Will compute the LDA weights and biases.
#######################################################################################################################
def makeLDAClassifier(class_list):
    for i in class_list:
        if i == 0:                                                          # Build pooled covariance, assumes that no-movment is always involved
            pooled_cov = pce.get_var('COV' + str(i)).to_np_array();
        else:
            tmpVal = pce.get_var('COV' + str(i)).to_np_array();
            pooled_cov += tmpVal

    num_classes = np.shape(class_list)
    pooled_cov = pooled_cov / num_classes[0]
    inv_pooled_cov = np.linalg.inv(pooled_cov)                              # Find the pooled inverse covariance matrix
    inv_pooled_cov = np.array(inv_pooled_cov, np.float64, order='F')
    pce.set_var("INVPOOL", inv_pooled_cov)

    for i in class_list:
        mVal = pce.get_var('MN' + str(i)).to_np_array();
        tmpWg = np.dot(inv_pooled_cov, mVal.T)
        tmpCg = -0.5 * (mVal.dot(inv_pooled_cov).dot(mVal.T))

        if i == 0:
            Wg = tmpWg;
            Cg = tmpCg;
        else:
            Wg = np.concatenate((Wg,tmpWg), axis=1)
            Cg = np.concatenate((Cg,tmpCg), axis=1)

    Wg = np.array(Wg, np.float64, order='F')                                # Seems like this is needed to make it type compatible with the PCE
    Cg = np.array(Cg, np.float64, order='F')

    return (Wg, Cg)

#######################################################################################################################
# Function    : calcRMS(args)
# args        : data, filtered data
# Description : This function is used to calculate the RMS feature
#######################################################################################################################
def calcRMS(data):
    N = data.shape[0]
    rms = np.sqrt(np.sum(np.square(data), axis = 0)/N)[np.newaxis]
    return rms

#######################################################################################################################
# Function    : collectLR(args)
# args        : mode, transmitted value for identification: feat_data, the current feature data: y_out, the output label
# Description : This function is used to create input and output arrays
#######################################################################################################################
def collectLR(mode, data, y_out):
    # Set collecting data flag to 1
    collect = 1
    
    # Get number of samples, x, and y variables
    N_T = pce.get_var('N_T').to_np_array()
    x = pce.get_var('X').to_np_array()
    y = pce.get_var('Y').to_np_array()
    
    # Add bias ones to x data
    temp = np.c_[data,1]
    
    # Save x and y data
    x[N_T[0, mode],:] = temp
    y[N_T[0, mode],:] = y_out
    pce.set_var('X',x.astype(float, order='F'))
    pce.set_var('Y',y.astype(float, order='F'))
    
    # Increment the number of x points
    N_T[0, mode] += 1
    
    if N_T[0, mode] == (sampThres - 1):
        # set collecting data flag to 0
        collect = 0
        # Increment the repetition variable to indicate a new training session has been completed.
        N_R = pce.get_var('N_R').to_np_array()
        N_R[0, mode] += 1        
        pce.set_var('N_R', N_R)
        # Set the mode back to its standby value of -1.
        pce.set_var('MODE', -1)
        # Set the y_label back to its standby value of -1.
        pce.set_var('Y_LABEL', -1)
        # Toggle the class_activate variable to 0.
        pce.set_var('CLASS_ACTIVE', 0)
        # Zero y matrix for all other models
        y_zeros = np.zeros((sampThres, 1), dtype=float, order='F')
        # Randomly select samples
        ind = np.random.choice(sampThres - 1, sampThres/(2*6), replace=False)
        # Update X'X and X'Y matrices
        for i in range(1, numModels):
            xx = pce.get_var('XX' + str(i)).to_np_array()
            xy = pce.get_var('XY' + str(i)).to_np_array()
            if i == mode:
                #if np.any(y):
                (xx, xy) = updateMat(xx, xy, x, y)
                #else:
                    #(xx, xy) = updateMat(xx, xy, x[0:sampThres-1:2,:], y[0:sampThres-1:2,:])
            else: #if np.any(y):
                (xx, xy) = updateMat(xx, xy, x, y_zeros)
                #(xx, xy) = updateMat(xx, xy, x[ind,:], y_zeros[ind,:])
            pce.set_var('XX' + str(i),xx.astype(float, order='F'))
            pce.set_var('XY' + str(i),xy.astype(float, order='F'))

    pce.set_var('N_T', N_T)
    
    return collect
#######################################################################################################################
# Function    : linearModel(args)
# args        : xx, X'X matrix: xy, X'Y matrix
# Description : This function is used to create the linear regression model
#######################################################################################################################
def linearModel(xx, xy):
    w = np.dot(np.linalg.inv(xx),xy)
    return w[:,0]

#######################################################################################################################
# Function    : predict(args)
# args        : x, input array: w, model weights
# Description : This function is used to predict a y value
#######################################################################################################################    
def predict(data, w):
    x = np.c_[data,1] # CHANGE THIS
    y = np.dot(x,w)
    return y
    
#######################################################################################################################
# Function    : calcMat(args)
# args        : x, input array: y, output array
# Description : This function is used to calculate X'X and X'Y
#######################################################################################################################
def calcMat(x, y):
    xx = np.dot(np.transpose(x),x)
    xy = np.dot(np.transpose(x),y)
    return (xx, xy)
    
#######################################################################################################################
# Function    : updateMat(args)
# args        : xx_old, old X'X matrix: xy_old, old X'Y matrix: x, input array: y, output array
# Description : This function is used to update X'X and X'Y
#######################################################################################################################
def updateMat(xx_old, xy_old, x, y):
    (xx, xy) = calcMat(x, y)
    xx_new = xx_old + xx
    xy_new = xy_old + xy
    return (xx_new, xy_new)

#######################################################################################################################
# Function    : udpateMVC(args)
# args        : mode, transmitted value for identification: cur_val, the current data: feat_data, the current feature data: update_var, argument to update.
# Description : This function is used to remove redundant repetition of code between no-movement and all other classes.
#######################################################################################################################
def updateMVC(mode):
    collect = 1
    # N_MVC: Total number of windows used
    MVC_T = pce.get_var('MVC_T').to_np_array()
    # Update number of samples
    MVC_T[0, mode] += 1
    # Set MVC samples
    pce.set_var('MVC_T', MVC_T)    
    # If collected max samples
    if MVC_T[0, mode] == (sampThres - 1):
        pce.set_var('CLASS_ACTIVE', 0)
        pce.set_var('MODE', -1)
        MVC_R = pce.get_var('MVC_R').to_np_array()
        MVC_R[0, mode] += 1
        pce.set_var('MVC_R', MVC_R)
        collect = 0
    
    return collect

#######################################################################################################################
# Function    : updateThresh(args)
# args        : 
# Description : This function is used to update the no movement threshold average
#######################################################################################################################
def updateThresh(cur_val):
    N_C = pce.get_var('N_C').to_np_array()
    thres_new = updateAverage(pce.get_var('THRESH_VAL'), cur_val, N_C[0, 0])
    N_C[0, 0] = N_C[0, 0] + 1
    pce.set_var('THRESH_VAL', thres_new)
    pce.set_var('N_C', N_C)

