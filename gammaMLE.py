import numpy as np
from scipy.special import digamma
from math import gamma
from scipy.special import polygamma
import random
import math

def newton_update(a_old, xbar, logxbar):
    """
    one iteration of the newton update
    :param a_old: old shape parameter value to be updated
    :param xbar: calculated value of xbar from EM
    :param logxbar: calculated value of logxbar from E<
    :return: k_new
    """
    #d1_l = np.log(1/theta) + logxbar - digamma(a_old)
  #  d2_l = -polygamma(1,a_old)

   # a_new = a_old - d1_l / d2_l
    update_val =  (np.log(a_old) - digamma(a_old) - np.log(xbar) + logxbar) \
                  / (1 / a_old - polygamma(1, a_old))


    a_new = a_old - update_val

    return(a_new)


def gamma_MLE(xbar, logxbar, obs, g_shape, g_scale):
    """
    Calculates approximate maximum likelihood estimates of gamma distribution parameters:
    k (shape), and theta (scale)
    :param xbar:  calculated value of xbar from EM
    :param logxbar:  calculated value of logxbar from EM
    :return:
    """

    # initialize with MOM estimators
    denom = 0
    num_obs = 0
    for cow in obs.values():
        for val in cow:
            denom += (float(val) - xbar)**2
            num_obs += 1

    a_mom = (num_obs * xbar**2) / (denom)
    a = [a_mom]


    # update
    a.append(newton_update(a[0], xbar, logxbar))

    i = 1
    while abs(a[i] - a[i - 1]) > 0.0001:
        if i > 10000:
            raise ValueError("too many NR iterations")
        if a[i] > 1000000000:
            raise ValueError("a too large, will not converge")
        a.append(newton_update(a[i], xbar, logxbar))
        i += 1

    theta = xbar / a[-1]

    return a[-1], theta
