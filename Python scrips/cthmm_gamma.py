#!/usr/bin/env python3

"""
Estimate hidden state transition rate matrix and log normal emission probabilities via EM.
Uses the Expm method of EM described in:
    Liu, Yu-Ying, et al. "Efficient learning of continuous-time hidden markov models
    for disease progression." Advances in neural information processing systems. 2015.

Usage:
   python3 cthmm_gamma.py
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import csv
from collections import defaultdict
import random
import pandas as pd
from scipy import linalg as la
from scipy.stats import gamma
import math
import gammaMLE

from sys import exit


# create output classes for fb and em functions
class FbOutput:
    fwd = dict()
    back = dict()
    Fmat = dict()
    likelihood = list()
    posterior = dict()
    max_posterior = dict()
    Ptrans = dict()
    final_likelihood = list()

class EmOutput:
    g_shape = list()
    g_scale = list()
    Qmat = np.ndarray
    ip = list()


def read_file(filename):
    """
    Reads csv and outputs a dictionary with keys cowID and elements
        CFU sequence, removes keys that only have one CFU observation
    :return:
    """

    with open(filename) as f:
        reader = csv.DictReader(f)
        cfu_dict = defaultdict(list)
        tau_dict = defaultdict(list)
        for row in reader:
            cfu_dict[row['CombinedID']].append(str(round(float(row['cor_totCFU']),2)))
            tau_dict[row['CombinedID']].append(row['time'])
        print(tau_dict)
    
        return cfu_dict, tau_dict


def get_probabilities(obs, ip, g_shape, g_scale):
    """
    Generates the transition and emission probabilities table given the parameters
    :param obs:  dictionary of CFU observations
    :param ip: initial state probabilities input by user
    :param g_shape: current estimate of gamma emission shape parameters
    :param g_scale: current estimate of gamma emission scale parameters

    :return:
    emission probabilities: emission pdf calculated for observed data using estimated gamma parameters
    initial probabilities: ip as a dictionary
    """
    
    nstates = len(ip)
    emission_probabilities = dict()
    for state in range(nstates):
        emission_probabilities[state] = dict()
    
    obs_list = []
    for seq in obs.values():
        for val in seq:
            obs_list.append(float(val))
    
    obs_list = np.unique(obs_list)
    
    max_obs = max(obs_list)
   # discrete_gamma = dict()
   # bin_vals = np.array(range(0,int(round(max_obs+1,2)*100)))/100
    for s in obs:
        for i in range(len(obs[s])):
            for state in range(nstates):
                emission_probabilities[state][(obs[s][i])] = gamma.pdf(float(obs[s][i]), a=g_shape[state], scale=g_scale[state])
       # sumbin = 0
       #discrete_gamma[state] = dict() #dictionary of dictionaries to hold gamma pmf
       # for bin in range(len(bin_vals)-1):
    
         #   discrete_gamma[state][bin_vals[bin]] = gamma.cdf(x=(bin_vals[bin+1]), a=g_shape[state], scale=g_scale[state]) \
         #    - gamma.cdf(x=bin_vals[bin], a= g_shape[state], scale= g_scale[state])
          #  sumbin += discrete_gamma[state][bin_vals[bin]]
        #for x in obs_list:
         #   for bin in range(len(bin_vals)-1):
          #      if x >= bin_vals[bin] and x < bin_vals[bin+1]:
           #         emission_probabilities[state][str(x)] = discrete_gamma[state][bin_vals[bin]]
##### issue is that if obs are integers (like 3) the obs will be '3' but the emission will be '3.0'
    init_probs = dict()
    for i in range(len(ip)):
        init_probs[i] = ip[i]
    
    return emission_probabilities, init_probs


def checkQ(Qmat, nstates):
    """
    Checks if rows of Q matrix sum to 0
    :param Qmat: transition rate matrix
    :param nstates: number of states
    :raises:
    Value error if rows don't sum to 0
    """

    for row in range(nstates):
        if round(sum(Qmat[row]), 3) != 0:
            print(Qmat)
            raise ValueError("Row of Q matrix doesn't sum to 0")


def setup_fb(Qmat, tau):
    """
    Setup parameters for the forward backward algorithm
    :param Qmat: transition rate matrix
    :param tau: dictionary of time intervals
    :return:
    P_trans: dictonary that stores a state transition probability matrix
             for each time interval
    F: dictionary, holds probability of transitioning from state k to state j
                at any time interval
    """

    Qmat = np.array(Qmat)
    
    nstates = np.shape(Qmat)[0]  # number of hidden states
    
    # throws error if rows of Qmat don't sum to 0
    checkQ(Qmat, nstates)
    
    # setup P_trans
    P_trans = dict()
    
    for seq in tau.values():
        for time in range(len(seq)):
            if seq[time] != 'NA':
                P_trans[seq[time]] = la.expm(Qmat * float(seq[time]))
    
    # setup F
    F = dict()
    for val in tau.values():
        for times in range(len(val)):
            if val[times] != 'NA':
                F[val[times]] = dict()
                for state in range(nstates):
                    F[val[times]][state] = dict()
    for t in F:
        F[t] = np.zeros((nstates, nstates))
    
    # setup Ftemp which holds Fmat values from each sequence.
    # at the end of iterating through all states, Ftemp[v] will be added to the
    # previous Ftemp[v]'s from all prior sequences. That sum is stored in F[v]
    Ftemp = dict()
    for q in tau.values():
        for tm in range(len(q)):
            if q[tm] != 'NA':
                Ftemp[q[tm]] = dict()
    
    return nstates, P_trans, F, Ftemp


def forward_backward(obs, tau, Qmat, ip, g_shape, g_scale):
    """
    Calculates the forward and backward probabilities of a given observation
    :param obs: sequence data to train on. Map from sequence name
                to list representing the sequence
    :param tau: dictionary of time intervals between consecutive samples from the same cow
    :param Qmat: transition rate matrix
    :param ip: dictionary of initial probabilities
    :param g_shape: current estimate of gamma emission shape parameters
    :param g_scale: current estimate of gamma emission scale parameters

    :return:
    fwd: dictionary of forward probabilities
             likelihood_f: P(obs) calculated using the forward algorithm
             back: dictionary of backward probabilities
             likelihood_b: P(obs) calculated using the backward algorithm
             posterior: dictionary of posterior probabilities (fwd*back/likelihood_f)
             max_posterior: state with highest posterior probability for each time interval
             P_trans: state transition probability matrix for each time interval
             F: dictionary, holds probability of transitioning from state k to state j at any time interval
    """
    
    # Forward probability calculation
    # initialization
    fwd = dict()
    back = dict()
    posterior = dict()
    max_posterior = dict()
    seq_likelihood = dict()
    ind_likelihood = dict()
    
    # get other setup parameters
    nstates, P_trans, F, Ftemp = setup_fb(Qmat, tau)
    emiss_probs, init_probs = get_probabilities(obs, ip, g_shape, g_scale)
    
    # For each sequence, calculate a separate set of
    # forward/backward probabilities
    for s in tau:
        # get transition probability from np.array Qmat
        # P_trans is the probability of state transmission, and this probabiliity
        # depends on the length of the time interval tau_v, so there is a different
        # Probability matrix for each time interval v
    
        fwd[s] = dict()
        # Initialize the forward probabilities.
        # Each state will have a list of probabilities for every observation in
        # the sequence
        for state in range(nstates):
            fwd[s][state] = [0] * len(obs[s])
            fwd[s][state][0] = init_probs[state] * emiss_probs[state][obs[s][0]]
    
        # For each observation in the sequence (starting with the second),
        # calculate the probability of transitioning from any state at the
        # previous observation to every state at the current observation
        i = 1  # i is the counter for CFU observation in the sequence
        for t in tau[s]:  # t is the actual time interval
            if t != 'NA': # (t == 'NA') == False:
                for k in range(nstates):
                    cum_sum = 0
                    # Sum probabilities of transitioning to state k from
                    # all possible states j at the previous observation
                    for j in range(nstates):
                        fwd_k = fwd[s][j][i - 1] * P_trans[t][j][k]
                        cum_sum += fwd_k
                    fwd[s][k][i] = emiss_probs[k][obs[s][i]] * cum_sum
                i += 1
    
        # Backward probability calculation
        # Initialize backward probabilities
        back[s] = dict()
        for k in range(nstates):
            # Initialize last backwards probability to some value
            back[s][k] = [0] * len(obs[s])
            back[s][k][len(obs[s]) - 1] = 1  # initialize to 1
    
        # i counts observations in an obs sequence
        # t is the acutal time interval length
        i = len(obs[s]) - 2
        for time in reversed(tau[s]):
            if time != 'NA':
                for k in range(nstates):
                    back_k = 0
                    for j in range(nstates):
                        back_k += (P_trans[time][k][j] * emiss_probs[j][obs[s][i + 1]] * back[s][j][i + 1])
                    back[s][k][i] = back_k
            i -= 1
    
        # log likelihood of the sequence = sum of all individual likelihoods for that sequence
        seq_likelihood[s] = 0
    
        # likelihood for each observation
        ind_likelihood[s] = [0] * len(obs[s])
        for i in range(len(obs[s])):
            sum_k = 0
            for k in range(nstates):
                sum_k += (fwd[s][k][i] * back[s][k][i])
            ind_likelihood[s][i] = sum_k
    
        like = 0
        for state in range(nstates):
            like += fwd[s][state][-1]
        seq_likelihood[s] = np.log10(like)
    
        # posterior probabilities
        posterior[s] = dict()
        for k in range(nstates):
            posterior[s][k] = [0] * len(obs[s])
            for i in range(len(obs[s])):
                posterior[s][k][i] = (fwd[s][k][i]
                                      * back[s][k][i]) / ind_likelihood[s][i]
    
        # posterior decoding: maximum posterior probability at each state
        max_posterior[s] = dict()
        max_posterior[s] = [0] * len(obs[s])
        for i in range(len(obs[s])):
            key_max = max(posterior[s].keys(), key=(lambda w: posterior[s][w][i]))
            max_posterior[s][i] = key_max
    
        # calculate probability of transitioning from state k to state j at any time interval
        # sum over time intervals
        # i is each element of obs
    
        i = 0  # i is the counter for CFU observation in the sequence
        for interval in tau[s]:  # t is the actual time interval
            if interval != 'NA':
                denom = 0
                Ftemp[interval] = np.zeros((nstates, nstates))
                for k in range(nstates):
                    for j in range(nstates):
                        temp = (fwd[s][k][i] * P_trans[interval][k][j] * emiss_probs[j][obs[s][i + 1]] *
                                back[s][j][i + 1])
                        Ftemp[interval][k][j] += temp
                        denom += temp
                # forget about making count table for now, see if it works without
                for k in range(nstates):
                    for j in range(nstates):
                        Ftemp[interval][k][j] = Ftemp[interval][k][j] / denom
                        F[interval][k][j] += Ftemp[interval][k][j]
    
                i += 1
    
    # calculating likelihood of all sequences, summing the log likelihood for all sequences
    final_likelihood = 0
    for s in seq_likelihood:
        final_likelihood += seq_likelihood[s]
    
    print("log likelihood" + str(final_likelihood))
    
    out = FbOutput()
    out.fwd = fwd
    out.back = back
    out.Fmat = F
    out.likelihood = ind_likelihood
    out.posterior = posterior
    out.max_posterior = max_posterior
    out.Ptrans = P_trans
    out.final_likelihood = final_likelihood
    return out


def setupEM(F, P, Qmat, nstates):
    """
    Sets up parameters for Expm method of EM
    :param F: counts of probabilities of transitioning from state k to j
    :param P: state transition probability matrix for each time interval
    :param Qmat: state transition rate matrix for each time interval
    :param nstates: number of hidden states
    :return:
    Exp_tau: Total amount of time chain remains in state k
    Exp_n: Expected number of transitions from state k to j
    """

    # Initialize helper arrays
    Z = np.zeros((nstates, nstates))  # component of matrix A, doesnt get updated
    D = np.zeros((nstates, nstates, nstates))  # one nxn matrix for each state
    
    # Initialize Exp_tau, a dictionary that stores expected number of transitions
    # in a time period, and
    # Initialize Exp_n, a dictionary that stores expected numbers of transtitions from
    # state k to state j
    Exp_tau = dict()
    for k in range(nstates):
        Exp_tau[k] = 0
    
    Exp_n = dict()
    for k in range(nstates):
        Exp_n[k] = dict()
        for j in range(nstates):
            Exp_n[k][j] = 0
    
    # perform Exmp method for EM
    # create matrix from F dictionary for easier computation
    Fmat = dict()
    for v in F:
        m = [0] * nstates
        m = np.array(m)
        m = m.reshape(1, nstates)
        for k in range(nstates):
            Ftemp = np.array(F[v][k])
            Ftemp = Ftemp.reshape(1, nstates)
            m = np.append(m, Ftemp, axis=0)
    
        m = np.delete(m, list(range(0, nstates)))
        m = m.reshape(nstates, nstates)
        Fmat[v] = m
    
    # get expected value of total time chain spent in state k
    for i in range(nstates):
        for time in F:
            B = np.zeros((nstates, nstates))
            B[i][i] = 1
            A = np.bmat([[Qmat, B], [Z, Qmat]])
    
            # keep only upper right quadrant of D
            # convert P back to non-log scale for ease of operations
            D[i] = ((la.expm(A * float(time)))[0:nstates, nstates:2 * nstates]) / (P[time])
            Exp_tau[i] += (Fmat[time] * D[i]).sum().sum()  # element wise multiplication of Fmat and
    
    # get expected value of number of transitions from i to j
    for i in range(nstates):
        for j in range(nstates):
    
            for interval in F:
                B = np.zeros((nstates, nstates))
                B[i][j] = 1
                A = np.bmat([[Qmat, B], [Z, Qmat]])
                N = Qmat[i][j] * (la.expm(A * float(interval))[0:nstates, nstates:2 * nstates]) / (P[interval])
    
                Exp_n[i][j] += (Fmat[interval] * N).sum().sum()
    
    return Exp_tau, Exp_n


def em(fb, obs, Qmat, g_shape, g_scale, ip):
    """
    Performs one EM stem
    :param fb: relevant output variables from forward-backward
    :param obs: the sequence under analysis
    :param Qmat: transition rate matrix
    :param g_shape: current estimate of gamma emission shape parameters
    :param g_scale: current estimate of gamma emission scale parameters
    :return:
    Qmat: updated transition rate matrix
    :param g_shape: updated estimate of gamma emission shape parameters
    :param g_scale: updated estimate of gamma emission scale parameters
    """

    fwd = fb.fwd
    back = fb.back
    Fmat = fb.Fmat
    Ptrans = fb.Ptrans
    Px = fb.likelihood
    posterior = fb.posterior
    nstates = np.shape(Qmat)[0]
    
    Exp_tau, Exp_n = setupEM(Fmat, Ptrans, Qmat, nstates)
    
    # update Q = exp number of transitions from j to j / exp total number
    # of transitions in time period
    # Qmat = np.zeros((nstates, nstates))
    
    for k in range(nstates):
        qsum = 0
        for j in range(nstates):
            if k != j:
                Qmat[k][j] = Exp_n[k][j] / Exp_tau[k]
                qsum += Qmat[k][j]
        Qmat[k][k] = -qsum
    
    # throw error if rows of Q don't sum to 0
    checkQ(Qmat, nstates)
    # quit(0)
    # Initialize E, a dictionary that stores expected counts of
    # emissions from states
    
    # For each sequence
    xbar = [0] * nstates
    logxbar = [0] * nstates
    cum_sum_numx  = [0] * nstates
    cum_sum_denom = [0] * nstates
    cum_sum_numlogx = [0] * nstates
    
    for s in obs:
        # Calculate E, a dictionary that stores expected number of times b
        # is observed while in state k
        for i in range(len(obs[s])):
            for k in range(nstates):
                #calculate xbar
                cum_sum_numx[k] += posterior[s][k][i] * float(obs[s][i])
                cum_sum_denom[k] += posterior[s][k][i]
    
                #calculate logxbar
                cum_sum_numlogx[k] +=  posterior[s][k][i] * np.log(float(obs[s][i]))
    # calculate new initial probs
    ip_temp = np.zeros((nstates,1))
    ip_denom = 0
    for k in range(nstates):
        for s in obs:
            ip_temp[k] += posterior[s][k][1]
    ip_denom = sum(ip_temp)
    for k in range(nstates):
        ip[k] = ip_temp[k]/ip_denom


    for k in range(nstates):
        xbar[k] = cum_sum_numx[k] / cum_sum_denom[k]
        logxbar[k] = cum_sum_numlogx[k] / cum_sum_denom[k]
        g_shape[k], g_scale[k] = gammaMLE.gamma_MLE(xbar[k], logxbar[k], obs, g_shape[k], g_scale[k])


    out = EmOutput()
    out.g_shape = g_shape
    out.g_scale = g_scale
    out.Qmat = Qmat
    out.ip = ip
    return out


def saveplot(log_likelihoods):
    """
    Helper function to save plot of log likelihoods over iterations to file for visualization.
    :param log_likelihoods: list of log likelihoods over iterations
    :return:
    plot of likelihoods to file
    """

    plt.title(" ")
    plt.xlabel("Iteration")
    plt.ylabel("Log likelihood")
    plt.plot(range(len(log_likelihoods)), log_likelihoods, 'r-')
    plt.savefig("3.3em_.svg")


def train(sequence, tau, Qmat, ip, g_shape, g_scale, stop_diff=0.0001):
    """
    Uses Expm method of EM to infer the parameters, iterating until
        a valid stopping condition is reached.
    :param sequence: sequence data to train on. Map from sequence name
                  to list representing the sequence
    :param tau: sequence of time intervals between observations
    :param Qmat: state transition rate matrix
    :param ip: initial state probabilities
    :param g_shape: current estimate of gamma emission shape parameters
    :param g_scale: current estimate of gamma emission scale parameters
    :param stop_diff: acceptable difference between likelihoods for convergence
    :return:
    Qmat: EM estimates of transition rate matrix
    ep: EM estimates of emission distribution means
    sigmas: EM estimates of emission distribution standard deviations
    AIC
    """

    # Two iterations of forward backward and em
    fb = forward_backward(sequence, tau, Qmat, ip, g_shape, g_scale)
    update = em(fb, sequence, Qmat, g_shape, g_scale, ip)
    log_likelihoods = [fb.final_likelihood]
    
    fb = forward_backward(sequence, tau, Qmat, ip, g_shape, g_scale)
    update = em(fb, sequence, Qmat, g_shape, g_scale, ip)
    log_likelihoods.append(fb.final_likelihood)
    
    i = 1
    while abs(log_likelihoods[i] - log_likelihoods[i - 1]) > stop_diff:
        fb = forward_backward(sequence, tau, Qmat, ip, g_shape, g_scale)
        update = em(fb, sequence, Qmat, g_shape, g_scale, ip)
    
        i = i + 1
        log_likelihoods.append(fb.final_likelihood)
    
    posterior = fb.posterior
    for s in posterior:
        for k in posterior[s]:
            for i in range(len(posterior[s][k])):
                posterior[s][k][i] = posterior[s][k][i]
    
    max_posterior = fb.max_posterior
    
    max_posterior_plot = dict()
    for s in max_posterior:
        max_posterior_plot[s] = [0] * len(max_posterior[s])
        for i in range(len(max_posterior[s])):
            if max_posterior[s][i] == 0:
                max_posterior_plot[s][i] = 0
            else:
                max_posterior_plot[s][i] = 1

#    max_posterior_df = pd.DataFrame.from_dict(max_posterior_plot, orient='index')
 #  max_posterior_df.to_csv(path_or_buf = "max_posterior_df3.csv", index = True)


    # split dictionary of dictionary to individual dictionaries of posterior
    # probabilities for each hidden state
    posterior0 = dict()
    for s in posterior:
        posterior0[s] = posterior[s][0]
    
    posterior1 = dict()
    for s in posterior:
       posterior1[s] = posterior[s][1]

   # posterior2 = dict()
   # for s in posterior:
   #     posterior2[s] = posterior[s][2]

   # posterior3 = dict()
   # for s in posterior:
   #     posterior3[s] = posterior[s][3]

   # posterior4 = dict()
   # for s in posterior:
   #     posterior4[s] = posterior[s][4]

    # create dataframes for posterior probabilities, and export to csv
   # posterior0_df_3 = pd.DataFrame.from_dict(posterior0, orient = 'index')
   # posterior1_df_3 = pd.DataFrame.from_dict(posterior1, orient = 'index')
   # posterior2_df_3 = pd.DataFrame.from_dict(posterior2, orient = 'index')
   # posterior3_df_5 = pd.DataFrame.from_dict(posterior3, orient= 'index')
    #posterior4_df_5 = pd.DataFrame.from_dict(posterior4, orient='index')
    
    #posterior0_df_3.to_csv(path_or_buf = "posterior0_df_3.csv", index = True)
    #posterior1_df_3.to_csv(path_or_buf = "posterior1_df_3.csv", index = True)
    #posterior2_df_3.to_csv(path_or_buf = "posterior2_df_3.csv", index = True)
    #posterior3_df_5.to_csv(path_or_buf = "posterior3_df_5.csv", index = True)
    #posterior4_df_5.to_csv(path_or_buf="posterior4_df_5.csv", index=True)
    
    print(log_likelihoods[-1])
    p = 2 * len(g_shape) + len(g_shape)*len(g_shape) - 1 #  #ip params - 1 + (nstates * nstates - nstates) + #shape + #scale
    AIC = -2 * log_likelihoods[-1] + 2 * p
    
    n = sum(len(x) for x in sequence.values())
    print(n)
    BIC = np.log(n)*p - 2 * log_likelihoods[-1]
    
    print("AIC: " + str(AIC))
    print("BIC: " + str(BIC))
    print(max_posterior)
    
    return Qmat, g_shape, g_scale, AIC, ip

def main():
    parser = argparse.ArgumentParser(description='Compute Qmat, and emiss probs via EM.')
    parser.add_argument('-f',
                        action="store",
                        dest="f", type=str,
                        default='hmm_ready_data.csv')
    parser.add_argument('-Qmat',
                        action="store",
                        dest="Qmat",
                        type=list,
                        default=[[-0.7976487,	0.59686338,	0.20078534],
                                 [2.26019584,	-6.86315115,	4.60295532],
                                 [0.99030801,	0.66373244,	-1.6540404]])

                                #[[-4.35400581,4.35400581],[0.38016362,-0.38016362]])
    
                                #
    
                                #[[-5.733124,	2.13314049,	1.15067234,	2.44931111],
                                #[1.71098830,	-4.177774,	1.10838173,	1.35840404],
                                #[4.12366395,	1.02293143,	-8.971035,	3.82443944],
                                #[4.8377172,	4.01029238,	1.26020116,	-10.108211]])
    
                                #[[-13.788765,	4.7907883,	3.65662112,	3.28642967,	2.05492582],
                                #[0.1150887,	-6.355558,	3.5294438,	1.57641998,	1.13460587],
                                #[2.87737005,	0.84359171,	-10.467957,	4.80898094,	1.9380146],
                                #[0.68257322,	2.5453775,	0.59526119,	-7.545565,	3.7223532],
                                #[1.4551189,	4.3630441,	2.57657329,	4.10989615,	-12.504632]])
    parser.add_argument('-ip',
                        action="store",
                        dest="ip",
                        type=list,
                        default= [0.281243205,	0.3950746154,0.323682180])
                            #[0.44420633,0.555793673])
                            #
                            #
                            #[0.232910679,	0.285337771,	0.4534965523,	0.0282549975])
                            #[0.047094735,0.140309388,0.31488274,0.125839746,0.371873390])
    parser.add_argument('-g_shape',
                        action="store",
                        dest="g_shape",
                        type=list,
                        default=[1.088811,9.882592,19.807027])
                                #[4.419186,10.179886])
                                #
                                #[2.830972,	8.648022,	13.369714,	18.879002])
                                #[9.688221,	11.377515,	14.952418,	15.012799,	18.503526])
    parser.add_argument('-g_scale',
                        action="store",
                        dest="g_scale",
                        type=list,
                        default=[0.4182357,	0.7552176,	0.4052078])
                                #[0.8357585,	0.4689414])
                                #
                                #[0.6226340,	0.9698231,	0.9578032,	0.4693736])
                                #[0.7431691,	0.6044293,	0.3611778,	0.2127581,	0.8774134])
    
    args = parser.parse_args()
    sequence, tau = read_file('/Users/kristinaceres/Development/HMM_project/Post_review_1/data/hmm_ready_data_with_categories_test_4.csv') #fake_data.csv')
    Qmat = args.Qmat
    ip = args.ip
    g_shape = args.g_shape
    g_scale = args.g_scale
    
    Qmat, g_shape, g_scale, AIC, ip = train(sequence, tau, Qmat, ip, g_shape, g_scale)
    print("ip")
    print(ip)
    
    print("Qmat")
    print(Qmat)


    print("gamma shape parameter")
    print(g_shape)
    
    print("gamma scale parameter")
    print(g_scale)


if __name__ == '__main__':
    main()
