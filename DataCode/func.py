# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 21:51:18 2020

@author: Tingting
"""

import numpy as np
import scipy.special as sc

def mulf(mat, vec):
    return 1 - ((1 - vec)[None,:]**mat).prod(axis = 1)

def lambda_eff_weibull(_alpha_inf, _beta_inf, _alpha_rem, _beta_rem):
    return ((_beta_rem/_beta_inf)**_alpha_inf)*(_alpha_inf/_alpha_rem)*sc.gamma(_alpha_inf/_alpha_rem)

def lambda_eff_inc_weibull(_alpha_inf, _beta_inf, _alpha_rem, _beta_rem, _tau):
    return ((_beta_rem / _beta_inf) ** _alpha_inf) * (_alpha_inf / _alpha_rem) \
           * (1 - sc.gammainc(_alpha_inf / _alpha_rem, (_tau / _beta_rem)**_alpha_rem)) \
           * sc.gamma(_alpha_inf / _alpha_rem) / np.exp(-(_tau / _beta_rem) ** _alpha_rem)


def srv_exponent(_tau, _lambda):
    return np.exp(- _lambda * _tau)

def srv_weibull_scale(_tau, _alpha, _beta):
    return np.exp(- (_tau / _beta) ** _alpha)

def srv_weibull_rate(_tau, _alpha, _theta):
    return np.exp(- (_tau * _theta) ** _alpha)

def srv_gamma_scale(_tau, _alpha, _beta):
    return 1 - sc.gammainc(_alpha, _tau / _beta)

def srv_gamma_rate(_tau, _alpha, _theta):
    return 1 - sc.gammainc(_alpha, _theta * _tau)

def srv_loglogistic_scale(_tau, _alpha, _beta):
    return 1 / ( 1 + (_tau / _beta) ** _alpha)

def srv_lognormal(_tau, _alpha, _beta):
    if _tau == 0:
        return 1
    else:
        return 0.5 - 0.5 * sc.erf((np.log(_tau) - _beta) / (_alpha * (2 ** 0.5)))

srv_func_dict = {'exponent': srv_exponent, 
                 'weibull_scale': srv_weibull_scale, 'weibull_rate': srv_weibull_rate, 
                 'gamma_scale': srv_gamma_scale, 'gamma_rate': srv_gamma_rate, 'loglogistic_scale': srv_loglogistic_scale}

def binary_search(func, y, xl, xr, lbound = -np.inf, rbound = np.inf, tol = 1e-6):
    def overbound(_x):
            if _x < lbound:
                return -1
            elif _x > rbound:
                return 1
            else:
                return 0    
            
    res = {'message': None, 'fun': None, 'x': None, 'success': False}
    if xl >= xr or lbound >= rbound or overbound(xl) or overbound(xr):
        res['message'] = 'the parameters are illegal.'
        return res
    yl = func(xl)
    yr = func(xr)
    if yl == yr:
        res['message'] = 'func(xl) is equal to func(xr).'
        return res
    elif yl < yr:
        xlower  = xl
        xhigher = xr
        ylower = yl
        yhigher = yr
    else:
        xlower  = xr
        xhigher = xl
        ylower = yr
        yhigher = yl
    overboundmark = False
    while True:
        if overboundmark and (y < ylower or y > yhigher):
            res['message'] = 'beyound bounds.'
            return res
        if y < ylower:
            delta_x = xlower - xhigher
            xhigher = xlower
            xlower = xlower + 2 * delta_x
            if overbound(xlower) == 1:
                xlower = rbound
                overboundmark = True
            elif overbound(xlower) == -1:
                xlower = lbound
                overboundmark = True
            else:
                pass
        elif y > yhigher:
            delta_x = xhigher - xlower
            xlower = xhigher
            xhigher = xhigher + 2 * delta_x
            if overbound(xhigher) == 1:
                xhigher = rbound
                overboundmark = True
            elif overbound(xhigher) == -1:
                xhigher = lbound
                overboundmark = True
            else:
                pass
        else:
            break
        ylower = func(xlower)
        yhigher = func(xhigher)   
    if abs((ylower - y) / y) < tol:
        res['message'] = 'the x has been found'
        res['success'] = True
        res['fun'] = ylower
        res['x'] = xlower
        return res
    if abs((yhigher - y) / y) < tol:
        res['message'] = 'the x has been found'
        res['success'] = True
        res['fun'] = yhigher
        res['x'] = xhigher
        return res
    while True:
        ytmp = func((xlower + xhigher) / 2)
        if abs((ytmp - y) / y) < tol:
            res['message'] = 'the x has been found'
            res['success'] = True
            res['fun'] = ytmp
            res['x'] = (xlower + xhigher) / 2
            return res
        else:
            if ytmp > y:
                xhigher = (xlower + xhigher) / 2
            else:
                xlower = (xlower + xhigher) / 2
        



































        
        
    