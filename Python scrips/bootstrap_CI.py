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
sequence, tau  =ghmm.read_file('/Users/kristinaceres/PhD/HMM_project/CT_HMM/data/hmm_ready_data.csv')

#initializaze bootstrap
n_iterations = 1000
n_size = int(len(sequence) * .5)


#creating lists to store transition probability parameters
# 2 states
qmat0_2 = list()
qmat1_2 = list()
qmat2_2 = list()
qmat3_2 = list()

# 3 states
qmat0_3 = list()
qmat1_3 = list()
qmat2_3 = list()
qmat3_3 = list()
qmat4_3 = list()
qmat5_3 = list()
qmat6_3 = list()
qmat7_3 = list()
qmat8_3 = list()


# 4 states
qmat0_4 = list()
qmat1_4 = list()
qmat2_4 = list()
qmat3_4 = list()
qmat4_4 = list()
qmat5_4 = list()
qmat6_4 = list()
qmat7_4 = list()
qmat8_4 = list()
qmat9_4 = list()
qmat10_4 = list()
qmat11_4 = list()
qmat12_4 = list()
qmat13_4 = list()
qmat14_4 = list()
qmat15_4 = list()

# 5 states
qmat0_5 = list()
qmat1_5 = list()
qmat2_5 = list()
qmat3_5 = list()
qmat4_5 = list()
qmat5_5 = list()
qmat6_5 = list()
qmat7_5 = list()
qmat8_5 = list()
qmat9_5 = list()
qmat10_5 = list()
qmat11_5 = list()
qmat12_5 = list()
qmat13_5 = list()
qmat14_5 = list()
qmat15_5 = list()
qmat16_5 = list()
qmat17_5 = list()
qmat18_5 = list()
qmat19_5 = list()
qmat20_5 = list()
qmat21_5 = list()
qmat22_5 = list()
qmat23_5 = list()
qmat24_5 = list()

#creating lists to store emission probability lists
# 2 states
shape0_2 = list()
shape1_2 = list()
scale0_2 = list()
scale1_2 = list()

mean0_2 = list()
mean1_2 = list()

# 3 states
shape0_3 = list()
shape1_3 = list()
shape2_3 = list()
scale0_3 = list()
scale1_3 = list()
scale2_3 = list()

mean0_3 = list()
mean1_3 = list()
mean2_3 = list()


#4 states
shape0_4 = list()
shape1_4 = list()
shape2_4 = list()
shape3_4 = list()
scale0_4 = list()
scale1_4 = list()
scale2_4 = list()
scale3_4 = list()

mean0_4 = list()
mean1_4 = list()
mean2_4 = list()
mean3_4 = list()

# 5 states
shape0_5 = list()
shape1_5 = list()
shape2_5 = list()
shape3_5 = list()
shape4_5 = list()
scale0_5 = list()
scale1_5 = list()
scale2_5 = list()
scale3_5 = list()
scale4_5 = list()

mean0_5 = list()
mean1_5 = list()
mean2_5 = list()
mean3_5 = list()
mean4_5 = list()

# initial probabilities
ip0_2 = list()
ip1_2 = list()

ip0_3 = list()
ip1_3 = list()
ip2_3 = list()


# 4 states
ip0_4 = list()
ip1_4 = list()
ip2_4 = list()
ip3_4 = list()

# 5 states
ip0_5 = list()
ip1_5 = list()
ip2_5 = list()
ip3_5 = list()
ip4_5 = list()

#confidence intervals
alpha = 0.05
p1 = ((alpha)/2) * 100
p2 = (1-alpha + ((alpha)/2)) * 100
'''
#run bootstrap 2 states ############################################################################
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
'''
# run bootstrap_CI for 4 state model ##########################################################################
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
                                                 #qmat
                                                 [[-5.733124, 2.13314049, 1.15067234, 2.44931111],
                                                  [1.71098830,	-4.177774,	1.10838173,	1.35840404],
                                                  [4.12366395,	1.02293143,	-8.971035,	3.82443944],
                                                  [4.8377172,	4.01029238,	1.26020116,	-10.108211]],
                                                 #ip
                                                 [0.232910679,	0.285337771,	0.4534965523,	0.0282549975],
                                                 #shape
                                                 [2.830972,	8.648022,	13.369714,	18.879002],
                                                 #scale
                                                 [0.6226340,	0.9698231,	0.9578032,	0.4693736])

        qmat0_4.append(Qmat[0][0])
        qmat1_4.append(Qmat[0][1])
        qmat2_4.append(Qmat[0][2])
        qmat3_4.append(Qmat[0][3])

        qmat4_4.append(Qmat[1][0])
        qmat5_4.append(Qmat[1][1])
        qmat6_4.append(Qmat[1][2])
        qmat7_4.append(Qmat[1][3])

        qmat8_4.append(Qmat[2][0])
        qmat9_4.append(Qmat[2][1])
        qmat10_4.append(Qmat[2][2])
        qmat11_4.append(Qmat[2][3])

        qmat12_4.append(Qmat[3][0])
        qmat13_4.append(Qmat[3][1])
        qmat14_4.append(Qmat[3][2])
        qmat15_4.append(Qmat[3][3])

        scale0_4.append(scale[0])
        scale1_4.append(scale[1])
        scale2_4.append(scale[2])
        scale3_4.append(scale[3])
        shape0_4.append(shape[0])
        shape1_4.append(shape[1])
        shape2_4.append(shape[2])
        shape3_4.append(shape[3])

        mean0_4.append(shape[0] * scale[0])
        mean1_4.append(shape[1] * scale[1])
        mean2_4.append(shape[2] * scale[2])
        mean3_4.append(shape[3] * scale[3])

        ip0_4.append(Ip[0])
        ip1_4.append(Ip[1])
        ip2_4.append(Ip[2])
        ip3_4.append(Ip[3])
        print(i)

    except:
        print(i)


trans_param_list4 = [qmat0_4, qmat1_4, qmat2_4, qmat3_4,
                     qmat4_4, qmat5_4, qmat6_4, qmat7_4,
                     qmat8_4, qmat9_4, qmat10_4, qmat11_4,
                     qmat12_4, qmat13_4, qmat14_4, qmat15_4]
trans_lower4 = list()
trans_upper4 = list()

for k in range(len(trans_param_list4)):
    trans_lower4.append(np.percentile(trans_param_list4[k], p1))
    trans_upper4.append(np.percentile(trans_param_list4[k], p2))

shape_param_list4 = [shape0_4, shape1_4, shape2_4, shape3_4]
scale_param_list4 = [scale0_4, scale1_4, scale2_4, scale3_4]
shape_lower4 = list()
shape_upper4 = list()
scale_lower4 = list()
scale_upper4 = list()

ip_param_list4 = [ip0_4, ip1_4, ip2_4, ip3_4]
ip_lower4 = list()
ip_upper4 = list()

mean_param_list4 = [mean0_4, mean1_4, mean2_4, mean3_4]
mean_upper4 = list()
mean_lower4 = list()

for k in range(len(shape_param_list4)):
    shape_lower4.append(np.percentile(shape_param_list4[k], p1))
    shape_upper4.append(np.percentile(shape_param_list4[k], p2))
    scale_lower4.append(np.percentile(scale_param_list4[k], p1))
    scale_upper4.append(np.percentile(scale_param_list4[k], p2))
    ip_lower4.append(max(0, np.percentile(ip_param_list4[k], p1)))
    ip_upper4.append(min(1, np.percentile(ip_param_list4[k], p2)))
    mean_lower4.append(np.percentile(mean_param_list4[k], p1))
    mean_upper4.append(np.percentile(mean_param_list4[k], p2))

trans_parameters4 = ['qmat04', 'qmat14', 'qmat24', 'qmat34',
                     'qmat44', 'qmat54', 'qmat64', 'qmat74',
                     'qmat84', 'qmat94', 'qmat104', 'qmat114',
                     'qmat124', 'qmat134', 'qmat144', 'qmat154']
trans_ci_df4 = pd.DataFrame({'parameters': trans_parameters4, 'lower': trans_lower4, 'upper': trans_upper4})
print(trans_ci_df4)

shape_parameters4 = ['shape04', 'shape14', 'shape24', 'shape34']
scale_parameters4 = ['scale04', 'scale14', 'scale24', 'scale34']
shape_ci_df4 = pd.DataFrame({'parameters': shape_parameters4, 'lower': shape_lower4, 'upper': shape_upper4})
print(shape_ci_df4)
scale_ci_df4 = pd.DataFrame({'parameters': scale_parameters4, 'lower': scale_lower4, 'upper': scale_upper4})
print(scale_ci_df4)
ip_parameters4 = ['ip04', 'ip14', 'ip24', 'ip34']
ip_ci_df4 = pd.DataFrame({'parameters': ip_parameters4, 'lower': ip_lower4, 'upper': ip_upper4})
print(ip_ci_df4)

mean_parameters4 = ['mean04', 'mean14', 'mean24', 'mean34']
mean_ci_df4 = pd.DataFrame({'parameters': mean_parameters4, 'lower': mean_lower4, 'upper': mean_upper4})
print(mean_ci_df4)

trans_ci_df4.to_csv(path_or_buf = "q_ci4_"+str(seed)+".csv", index = True)
shape_ci_df4.to_csv(path_or_buf = "shape_ci4_"+str(seed)+".csv", index = True)
scale_ci_df4.to_csv(path_or_buf = "scale_ci4_"+str(seed)+".csv", index = True)
ip_ci_df4.to_csv(path_or_buf = "ip_ci4_"+str(seed)+".csv", index = True)
mean_ci_df4.to_csv(path_or_buf = "ip_ci4_"+str(seed)+".csv", index = True)

# run bootstrap_CI for 5 state model ##########################################################################
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
                                                 #qmat
                                                 [[-13.788765, 4.7907883, 3.65662112, 3.28642967, 2.05492582],
                                                  [0.1150887,	-6.355558,	3.5294438,	1.57641998,	1.13460587],
                                                  [2.87737005,	0.84359171,	-10.467957,	4.80898094,	1.9380146],
                                                  [0.68257322,	2.5453775,	0.59526119,	-7.545565,	3.7223532],
                                                  [1.4551189,	4.3630441,	2.57657329,	4.10989615,	-12.504632]],
                                                 #ip
                                                 [0.047094735,0.140309388,0.31488274,0.125839746,0.371873390],
                                                 #shape
                                                 [9.688221, 11.377515, 14.952418, 15.012799, 18.503526],
                                                 #scale
                                                 [0.7431691,	0.6044293,	0.3611778,	0.2127581,	0.8774134])

        qmat0_5.append(Qmat[0][0])
        qmat1_5.append(Qmat[0][1])
        qmat2_5.append(Qmat[0][2])
        qmat3_5.append(Qmat[0][3])
        qmat4_4.append(Qmat[0][4])

        qmat5_5.append(Qmat[1][0])
        qmat6_5.append(Qmat[1][1])
        qmat7_5.append(Qmat[1][2])
        qmat8_5.append(Qmat[1][3])
        qmat9_5.append(Qmat[1][4])

        qmat10_5.append(Qmat[2][0])
        qmat11_5.append(Qmat[2][1])
        qmat12_5.append(Qmat[2][2])
        qmat13_5.append(Qmat[2][3])
        qmat14_5.append(Qmat[2][4])

        qmat15_5.append(Qmat[3][0])
        qmat16_5.append(Qmat[3][1])
        qmat17_5.append(Qmat[3][2])
        qmat18_5.append(Qmat[3][3])
        qmat19_5.append(Qmat[3][4])

        qmat20_5.append(Qmat[4][0])
        qmat21_5.append(Qmat[4][1])
        qmat22_5.append(Qmat[4][2])
        qmat23_5.append(Qmat[4][3])
        qmat24_5.append(Qmat[4][4])


        scale0_5.append(scale[0])
        scale1_5.append(scale[1])
        scale2_5.append(scale[2])
        scale3_5.append(scale[3])
        scale4_5.append(scale[4])
        shape0_5.append(shape[0])
        shape1_5.append(shape[1])
        shape2_5.append(shape[2])
        shape3_5.append(shape[3])
        shape4_5.append(shape[4])

        mean0_5.append(shape[0] * scale[0])
        mean1_5.append(shape[1] * scale[1])
        mean2_5.append(shape[2] * scale[2])
        mean3_5.append(shape[3] * scale[3])
        mean4_5.append(shape[4] * scale[4])

        ip0_5.append(Ip[0])
        ip1_5.append(Ip[1])
        ip2_5.append(Ip[2])
        ip3_5.append(Ip[3])
        ip4_5.append(Ip[4])
        print(i)

    except:
        print(i)

trans_param_list5 = [qmat0_5, qmat1_5, qmat2_5, qmat3_5, qmat4_5,
                     qmat5_5, qmat6_5, qmat7_5, qmat8_5, qmat9_5,
                     qmat10_5, qmat11_5, qmat12_5, qmat13_5, qmat14_5,
                     qmat15_5, qmat16_5, qmat17_5, qmat18_5, qmat19_5,
                     qmat20_5, qmat21_5, qmat22_5, qmat23_5, qmat24_5]
trans_lower5 = list()
trans_upper5 = list()

for k in range(len(trans_param_list5)):
    trans_lower5.append(np.percentile(trans_param_list5[k], p1))
    trans_upper5.append(np.percentile(trans_param_list5[k], p2))

shape_param_list5 = [shape0_5, shape1_5, shape2_5, shape3_5, shape4_5]
scale_param_list5 = [scale0_5, scale1_5, scale2_5, scale3_5, scale4_5]
shape_lower5 = list()
shape_upper5 = list()
scale_lower5 = list()
scale_upper5 = list()


ip_param_list5 = [ip0_5, ip1_5, ip2_5, ip3_5, ip4_5]
ip_lower5 = list()
ip_upper5 = list()

mean_param_list5 = [mean0_5, mean1_5, mean2_5, mean3_5, mean4_5]
mean_upper5 = list()
mean_lower5 = list()

for k in range(len(shape_param_list5)):
    shape_lower5.append(np.percentile(shape_param_list5[k], p1))
    shape_upper5.append(np.percentile(shape_param_list5[k], p2))
    scale_lower5.append(np.percentile(scale_param_list5[k], p1))
    scale_upper5.append(np.percentile(scale_param_list5[k], p2))
    ip_lower5.append(max(0, np.percentile(ip_param_list5[k], p1)))
    ip_upper5.append(min(1, np.percentile(ip_param_list5[k], p2)))
    mean_lower5.append(np.percentile(mean_param_list5[k], p1))
    mean_upper5.append(np.percentile(mean_param_list5[k], p2))

trans_parameters5 = ['qmat05', 'qmat15', 'qmat25', 'qmat35', 'qmat45',
                     'qmat55', 'qmat65', 'qmat75', 'qmat85', 'qmat95',
                     'qmat105', 'qmat115', 'qmat125', 'qmat135', 'qmat145',
                     'qmat155', 'qmat165' 'qmat175', 'qmat185', 'qmat195',
                     'qmat205', 'qmat215', 'qmat225', 'qmat235', 'qmat245']
trans_ci_df5 = pd.DataFrame({'parameters': trans_parameters5, 'lower': trans_lower5, 'upper': trans_upper5})
print(trans_ci_df5)

shape_parameters5 = ['shape05', 'shape15', 'shape25', 'shape35', 'shape45']
scale_parameters5 = ['scale05', 'scale15', 'scale25', 'scale35' 'scale45']
shape_ci_df5 = pd.DataFrame({'parameters': shape_parameters5, 'lower': shape_lower5, 'upper': shape_upper5})
print(shape_ci_df5)
scale_ci_df5 = pd.DataFrame({'parameters': scale_parameters5, 'lower': scale_lower5, 'upper': scale_upper5})
print(scale_ci_df5)
ip_parameters4 = ['ip05', 'ip15', 'ip25', 'ip35', 'ip45']
ip_ci_df5 = pd.DataFrame({'parameters': ip_parameters5, 'lower': ip_lower5, 'upper': ip_upper5})
print(ip_ci_df5)

mean_parameters5 = ['mean05', 'mean15', 'mean25', 'mean35', 'mean45']
mean_ci_df5 = pd.DataFrame({'parameters': mean_parameters5, 'lower': mean_lower5, 'upper': mean_upper5})
print(mean_ci_df5)

trans_ci_df5.to_csv(path_or_buf = "q_ci5_"+str(seed)+".csv", index = True)
shape_ci_df5.to_csv(path_or_buf = "shape_ci5_"+str(seed)+".csv", index = True)
scale_ci_df5.to_csv(path_or_buf = "scale_ci5_"+str(seed)+".csv", index = True)
ip_ci_df5.to_csv(path_or_buf = "ip_ci5_"+str(seed)+".csv", index = True)
mean_ci_df5.to_csv(path_or_buf = "ip_ci5_"+str(seed)+".csv", index = True)
