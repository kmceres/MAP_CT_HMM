import cthmm_gamma as ghmm
import pandas as pd
import time
import numpy as np


#generate seed
seed = 1585510210#int(time.time()) % (2**32-1)  1585510210
print(seed)
#set seed
np.random.seed(seed)

#generate samp_size # of random uniform samples of each variable
samp_size = 250

for n in range(2,4): #1585341561
    #set number of model states
    nstates = n

    #initialize data structures
    qmat = np.zeros((samp_size, nstates, nstates))
    ip = np.zeros((samp_size, nstates, 1))
    shape = np.zeros((samp_size, nstates,1))
    scale = np.zeros((samp_size, nstates,1))
    
    for i in range(samp_size):
        # get random Q matrix values, keep condition that rows must add to 0
        # want values of Q matrix to be small enough so that the matrix exponential changes with time interval length
        qmat[i] = np.random.uniform(low=0.01, high = 5, size=(nstates, nstates))
        for k in range(nstates):
            qmat[i][k][k] = 0
            qsum = 0
            for j in range(nstates):
                qsum += qmat[i][k][j]
            qmat[i][k][k] = -qsum
    
        #get random initial probabilities, must sum to 1
        ip[i] = (np.random.uniform(low=0, high=1, size=nstates)).reshape(nstates,1)
        denom = sum(ip[i])
        for k in range(len(ip[i])):
            ip[i][k] = ip[i][k] / denom


        #get random gamma shape parameters
        shape[i] = np.random.uniform(low=1, high = 20, size=(nstates, 1))
        shape = np.sort(shape, axis = 1)
    
        #get random gamma scale parameters
        scale[i] = np.random.uniform(low=0.1, high=1, size=(nstates, 1))


    sequence, tau = ghmm.read_file('/Users/kristinaceres/Desktop/CT_HMM/hmm_ready_data.csv')
    
    #initialize data structures
    Qmat = np.zeros((samp_size, nstates,nstates))
    AIC = [0]*samp_size
    g_shape = np.zeros((samp_size, nstates,1))
    g_scale = np.zeros((samp_size, nstates,1))
    Ip = np.zeros((samp_size, nstates, 1))
    #want most cows in low states
    #ip = -np.sort(-ip, axis = 1)
    
    #initialize dataframes for storage and export
    qmat_df = pd.DataFrame()
    scale_df = pd.DataFrame()
    shape_df = pd.DataFrame()
    ip_df = pd.DataFrame()
    aic_df = pd.DataFrame()
    qmat_df_out = pd.DataFrame()
    scale_df_out = pd.DataFrame()
    shape_df_out = pd.DataFrame()
    ip_df_out = pd.DataFrame()
    
    #initialize shape of dataframes
    qmat_df[0] = [0]*(nstates*nstates)
    scale_df[0] = [0]*nstates
    shape_df[0] = [0]*nstates
    ip_df[0] = [0]*nstates
    aic_df[0] = [0]
    qmat_df_out[0] = [0]*(nstates*nstates)
    scale_df_out[0] = [0]*nstates
    shape_df_out[0] = [0]*nstates
    ip_df_out[0] = [0]*nstates


    for sim in range(samp_size):
        try:
            # save inputs
            qmat_df[sim] = qmat[sim].reshape(nstates * nstates,1)
            scale_df[sim] = scale[sim]
            shape_df[sim] = shape[sim]
            ip_df[sim] = ip[sim]
    
            print("shape :" + str(shape[sim]))
            print("scale :" + str(scale[sim]))
            print("ip: "+ str(ip[sim]))
            print(qmat[sim])
    
            Qmat[sim], g_shape[sim], g_scale[sim], AIC[sim], Ip[sim]= ghmm.train(sequence, tau, qmat[sim], ip[sim], shape[sim], scale[sim])


        except:
            print(" EXCEPTION __________________________________________________")
            print(str(sim))
            print("qmat ")
            print(qmat[sim])
            print("scale")
            print(scale[sim])
            print("shape")
            print(shape[sim])
            print("ip")
            print(ip[sim])
            print("____________________________________________________________")
            pass
    
        #save outputs
        aic_df[sim] = AIC[sim]
        qmat_df_out[sim]= Qmat[sim].reshape(nstates * nstates, 1)
        scale_df_out[sim]= g_scale[sim]
        shape_df_out[sim]=g_shape[sim]
        ip_df_out[sim] = Ip[sim]
    
        print(sim)
    
    print("seed: " + str(seed))
    #write to csv
    qmat_df.to_csv(path_or_buf="qmat_"+str(n)+str(seed)+".csv", index=True)
    scale_df.to_csv(path_or_buf="scale_"+str(n)+str(seed)+".csv",index=True)
    shape_df.to_csv(path_or_buf="shape_"+str(n)+str(seed)+".csv",index=True)
    ip_df.to_csv(path_or_buf="ip_"+str(n)+str(seed)+".csv",index=True)
    aic_df.to_csv(path_or_buf="aic_"+str(n)+str(seed)+".csv",index=True)
    qmat_df_out.to_csv(path_or_buf="qmat_out"+str(n)+str(seed)+".csv",index=True)
    scale_df_out.to_csv(path_or_buf="scale_out"+str(n)+str(seed)+".csv",index=True)
    shape_df_out.to_csv(path_or_buf="shape_out"+str(n)+str(seed)+".csv",index=True)
    ip_df_out.to_csv(path_or_buf="ip_out"+str(n)+str(seed)+".csv", index=True)


