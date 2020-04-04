#!/usr/bin/env python3

''' calculate bootstrap confidence intervals for model parameters
Arguments:
    - can edit the filename and the number of iterations
Outputs:
   95% confidence intervals for tp and ep parameters
Usage:
   python3 bootstrap_CI.py
'''

import numpy as np
import pandas as pd
from sklearn.utils import resample
import cthmm_gamma as ghmm
import time

#adapted from https://machinelearningmastery.com/calculate-bootstrap-confidence-intervals-machine-learning-results-python/
seed = int(time.time()) % (2**32-1)
print(seed)
#set seed
np.random.seed(seed)

#load data
sequence, tau  =ghmm.read_file('/Users/kristinaceres/Development/HMM_project/Post_review_1/data/hmm_ready_data_with_categories_test_4.csv')

#initializaze bootstrap
n_iterations = 1000
n_size = int(len(sequence) * .5)

#run bootstrap 2 states ############################################################################
#creating lists to store transition probability parameters
qmat0_2 = list()
qmat1_2 = list()
qmat2_2 = list()
qmat3_2 = list()

qmat0_3 = list()
qmat1_3 = list()
qmat2_3 = list()
qmat3_3 = list()
qmat4_3 = list()
qmat5_3 = list()
qmat6_3 = list()
qmat7_3 = list()
qmat8_3 = list()

#creating lists to store emission probability lists
shape0_2 = list()
shape1_2 = list()
scale0_2 = list()
scale1_2 = list()

mean0_2 = list()
mean1_2 = list()

shape0_3 = list()
shape1_3 = list()
shape2_3 = list()
scale0_3 = list()
scale1_3 = list()
scale2_3 = list()

mean0_3 = list()
mean1_3 = list()
mean2_3 = list()

# initial probabilities
ip0_2 = list()
ip1_2 = list()

ip0_3 = list()
ip1_3 = list()
ip2_3 = list()

for i in range(n_iterations):
    try:
    #prepare train and test sets
        cowid_sample = resample(list(sequence.keys()), n_samples=n_size)
        train_set = dict()
        tau_train = dict()
        for key in cowid_sample:
            train_set[key] = sequence[key]
            tau_train[key] = tau[key]
        #fit model
        Qmat, shape, scale, AIC, Ip= ghmm.train(train_set,
                                                    tau_train,
                                                    [[-4.35400581,4.35400581],[0.38016362,-0.38016362]],
                                                    [0.44420633,0.555793673],
                                                    [4.419186,10.179886],
                                                    [0.8357585,	0.4689414])
        qmat0_2.append(Qmat[0][0])
        qmat1_2.append(Qmat[0][1])
        qmat2_2.append(Qmat[1][0])
        qmat3_2.append(Qmat[1][1])

        scale0_2.append(scale[0])
        scale1_2.append(scale[1])
        shape0_2.append(shape[0])
        shape1_2.append(shape[1])

        mean0_2.append(shape[0]*scale[0])
        mean1_2.append(shape[1]*scale[1])

        ip0_2.append(Ip[0])
        ip1_2.append(Ip[1])
        print(i)
    except:
        print(i)


#confidence intervals
alpha = 0.05
p1 = ((alpha)/2) * 100
p2 = (1-alpha + ((alpha)/2)) * 100

trans_param_list2 = [qmat0_2, qmat1_2, qmat2_2, qmat3_2]
trans_lower2 = list()
trans_upper2 = list()

for k in range(len(trans_param_list2)):
    trans_lower2.append(np.percentile(trans_param_list2[k], p1))
    trans_upper2.append(np.percentile(trans_param_list2[k], p2))

shape_param_list2 = [shape0_2, shape1_2]
scale_param_list2 = [scale0_2, scale1_2]
shape_lower2 = list()
shape_upper2 = list()
scale_lower2 = list()
scale_upper2 = list()
ip_param_list2 = [ip0_2, ip1_2]
ip_lower2 = list()
ip_upper2 = list()
mean_param_list2 = [mean0_2, mean1_2]
mean_lower2 = list()
mean_upper2 = list()


for k in range(len(shape_param_list2)):
    shape_lower2.append(np.percentile(shape_param_list2[k], p1))
    shape_upper2.append(np.percentile(shape_param_list2[k], p2))
    scale_lower2.append(np.percentile(scale_param_list2[k], p1))
    scale_upper2.append(np.percentile(scale_param_list2[k], p2))
    ip_lower2.append(max(0, np.percentile(ip_param_list2[k], p1)))
    ip_upper2.append(min(1, np.percentile(ip_param_list2[k], p2)))
    mean_lower2.append(np.percentile(mean_param_list2[k], p1))
    mean_upper2.append(np.percentile(mean_param_list2[k], p2))


trans_parameters2 = ['qmat02', 'qmat12', 'qmat22', 'qmat32']
trans_ci_df2 = pd.DataFrame({'parameters': trans_parameters2, 'lower': trans_lower2, 'upper': trans_upper2})
print(trans_ci_df2)

shape_parameters2 = ['shape0_2', 'shape1_2']
scale_parameters2 = ['scale0_2', 'scale1_2']
shape_ci_df2 = pd.DataFrame({'parameters': shape_parameters2, 'lower': shape_lower2, 'upper': shape_upper2})
print(shape_ci_df2)

scale_ci_df2 = pd.DataFrame({'parameters': scale_parameters2, 'lower': scale_lower2, 'upper': scale_upper2})
print(scale_ci_df2)

ip_parameters2 = ['ip0_2', 'ip1_2']
ip_ci_df2 = pd.DataFrame({'parameters': ip_parameters2, 'lower': ip_lower2, 'upper': ip_upper2})
print(ip_ci_df2)

mean_parameters2 = ['mean0_2', 'mean1_2']
mean_ci_df2 = pd.DataFrame({'parameters': mean_parameters2, 'lower': mean_lower2, 'upper': mean_upper2})
print(mean_ci_df2)

trans_ci_df2.to_csv(path_or_buf = "q_ci2_"+str(seed)+".csv", index = True)
shape_ci_df2.to_csv(path_or_buf = "shape_ci2_"+str(seed)+".csv", index = True)
scale_ci_df2.to_csv(path_or_buf = "scale_ci2_"+str(seed)+".csv", index = True)
ip_ci_df2.to_csv(path_or_buf= "ip_ci2_"+str(seed)+".csv", index = True)
mean_ci_df2.to_csv(path_or_buf= "mean_ci2_"+str(seed)+".csv", index = True)
#run bootstrap 3 states ############################################################################

for i in range(n_iterations):
    # prepare train and test sets
    try:
        cowid_sample = resample(list(sequence.keys()), n_samples=n_size)
        train_set = dict()
        tau_train = dict()
        for key in cowid_sample:
            train_set[key] = sequence[key]
            tau_train[key] = tau[key]
        # fit model
        Qmat, shape, scale, AIC, Ip = ghmm.train(train_set,
                                                 tau_train,
                                                 [[-0.7976487, 0.59686338, 0.20078534],
                                                  [2.26019584, -6.86315115, 4.60295532],
                                                  [0.99030801, 0.66373244, -1.6540404]],
                                                 [0.281243205,	0.3950746154,0.323682180],
                                                 [1.088811,9.882592,19.807027],
                                                 [0.4182357,	0.7552176,	0.4052078])
        qmat0_3.append(Qmat[0][0])
        qmat1_3.append(Qmat[0][1])
        qmat2_3.append(Qmat[0][2])
        qmat3_3.append(Qmat[1][0])
        qmat4_3.append(Qmat[1][1])
        qmat5_3.append(Qmat[1][2])
        qmat6_3.append(Qmat[2][0])
        qmat7_3.append(Qmat[2][1])
        qmat8_3.append(Qmat[2][2])

        scale0_3.append(scale[0])
        scale1_3.append(scale[1])
        scale2_3.append(scale[2])
        shape0_3.append(shape[0])
        shape1_3.append(shape[1])
        shape2_3.append(shape[2])

        mean0_3.append(shape[0] * scale[0])
        mean1_3.append(shape[1] * scale[1])
        mean2_3.append(shape[2] * scale[2])

        ip0_3.append(Ip[0])
        ip1_3.append(Ip[1])
        ip2_3.append(Ip[2])
        print(i)
    except:
        print(i)


trans_param_list3 = [qmat0_3, qmat1_3, qmat2_3,
                     qmat3_3, qmat4_3, qmat5_3,
                     qmat6_3, qmat7_3, qmat8_3]
trans_lower3 = list()
trans_upper3 = list()

for k in range(len(trans_param_list3)):
    trans_lower3.append(np.percentile(trans_param_list3[k], p1))
    trans_upper3.append(np.percentile(trans_param_list3[k], p2))

shape_param_list3 = [shape0_3, shape1_3, shape2_3]
scale_param_list3 = [scale0_3, scale1_3, scale2_3]
shape_lower3 = list()
shape_upper3 = list()
scale_lower3 = list()
scale_upper3 = list()

ip_param_list3 = [ip0_3, ip1_3, ip2_3]
ip_lower3 = list()
ip_upper3 = list()

mean_param_list3 = [mean0_3, mean1_3, mean2_3]
mean_upper3 = list()
mean_lower3 = list()

for k in range(len(shape_param_list3)):
    shape_lower3.append(np.percentile(shape_param_list3[k], p1))
    shape_upper3.append(np.percentile(shape_param_list3[k], p2))
    scale_lower3.append(np.percentile(scale_param_list3[k], p1))
    scale_upper3.append(np.percentile(scale_param_list3[k], p2))
    ip_lower3.append(max(0, np.percentile(ip_param_list3[k], p1)))
    ip_upper3.append(min(1, np.percentile(ip_param_list3[k], p2)))
    mean_lower3.append(np.percentile(mean_param_list3[k], p1))
    mean_upper3.append(np.percentile(mean_param_list3[k], p2))

trans_parameters3 = ['qmat03', 'qmat13', 'qmat23', 'qmat33', 'qmat43', 'qmat53', 'qmat63', 'qmat73', 'qmat83']
trans_ci_df3 = pd.DataFrame({'parameters': trans_parameters3, 'lower': trans_lower3, 'upper': trans_upper3})
print(trans_ci_df3)

shape_parameters3 = ['shape03', 'shape13', 'shape23']
scale_parameters3 = ['scale03', 'scale13', 'scale23']
shape_ci_df3 = pd.DataFrame({'parameters': shape_parameters3, 'lower': shape_lower3, 'upper': shape_upper3})
print(shape_ci_df3)
scale_ci_df3 = pd.DataFrame({'parameters': scale_parameters3, 'lower': scale_lower3, 'upper': scale_upper3})
print(scale_ci_df3)
ip_parameters3 = ['ip03', 'ip13', 'ip23']
ip_ci_df3 = pd.DataFrame({'parameters': ip_parameters3, 'lower': ip_lower3, 'upper': ip_upper3})
print(ip_ci_df3)

mean_parameters3 = ['mean03', 'mean13', 'mean23']
mean_ci_df3 = pd.DataFrame({'parameters': mean_parameters3, 'lower': mean_lower3, 'upper': mean_upper3})
print(mean_ci_df3)

trans_ci_df3.to_csv(path_or_buf = "q_ci3_"+str(seed)+".csv", index = True)
shape_ci_df3.to_csv(path_or_buf = "shape_ci3_"+str(seed)+".csv", index = True)
scale_ci_df3.to_csv(path_or_buf = "scale_ci3_"+str(seed)+".csv", index = True)
ip_ci_df3.to_csv(path_or_buf = "ip_ci3_"+str(seed)+".csv", index = True)
mean_ci_df3.to_csv(path_or_buf = "ip_ci3_"+str(seed)+".csv", index = True)