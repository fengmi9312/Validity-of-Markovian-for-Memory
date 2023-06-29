# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 22:28:54 2023

@author: fengmi9312
"""

import numpy as np
from scipy.special import gamma, erf, gammainc

#Weibull survival function
def weibull_srv_func(_tau, _alpha, _beta):
    return np.exp(-(_tau / _beta) ** _alpha)

#gamma survival function
def gamma_srv_func(_tau, _alpha, _beta):
    return 1 - gammainc(_alpha, _tau / _beta)

#log-normal survival function
def lognormal_srv_func(_tau, _alpha, _beta):
    return 0.5 - 0.5 * erf((np.log(_tau) - _alpha)/(_beta * (2**0.5)))

#calculate the mean value of Weibull distribution
def weibull_mean(_alpha, _beta):
    return _beta * gamma(1 + 1 / _alpha)

#calculate the mean value of gamma distribution
def gamma_mean(_alpha, _beta):
    return _alpha * _beta

#calculate the mean value of log-normal distribution
def lognormal_mean(_alpha, _beta):
    return np.exp(_alpha + (_beta ** 2) / 2)

#generating survival function
def srv(_srv_func, _alpha, _beta, _length, _step):
    return _srv_func(np.arange(_length + 1) * _step, _alpha, _beta)

#generating distribution
def dist(_srv_func, _alpha, _beta, _length, _step):
    return (_srv_func(np.arange(_length) * _step, _alpha, _beta) - 
           _srv_func(np.arange(1,_length + 1) * _step, _alpha, _beta)) / _step

#generating hazard function
def haz(_srv_func, _alpha, _beta, _length, _step):
    return (1 - 
           _srv_func(np.arange(1,_length + 1) * _step, _alpha, _beta) / 
           _srv_func(np.arange(_length) * _step, _alpha, _beta)) / _step

#generating the generation time distribution, as well as its mean value and the value of lambda_eff
def gen(_srv_func_inf, _alpha_inf, _beta_inf, _srv_func_rem, _alpha_rem, _beta_rem,
        _length, _step):
    _tol = 1e-6
    _res = []
    _srv_inf = 0
    _i = 0
    _tau = 0
    _lam = 0
    _term = 0
    _mean = 0
    while True:
        _tau = _i * _step
        _srv_inf = _srv_func_inf(_tau, _alpha_inf, _beta_inf)
        if _i < _length:
            if _srv_inf == 0:
                break
        else:
            if _srv_inf < _tol:
                break
        _term = (1 - _srv_func_inf(_tau+_step, _alpha_inf, _beta_inf) / 
                     _srv_inf) * _srv_func_rem(_tau, _alpha_rem, _beta_rem)
        _res.append(_term)
        _lam += _term
        _mean += _tau * _term
        _i += 1
    return np.array(_res) / _lam / _step, _mean / _lam, _lam
    
#g_func(g) = 0 is the Euler-Lokta equation
def g_func(_g, _r0, _psi_gen, _step):
    _tau = np.arange(len(_psi_gen)) * _step
    return _r0 * (np.exp(-_g * _tau) * _psi_gen * _step).sum() - 1

#find a appropriate duration containing the solution of g_func(g) = 0
def findDuration(_r0, _psi_gen, _step):
    _left, _right = -1, 1
    while g_func(_left, _r0, _psi_gen, _step) <= 0:
        _left *= 2
    while g_func(_right, _r0, _psi_gen, _step) >= 0:
        _right *= 2
    return _left, _right
    
#using bisection method to find the solution of g_func(g) = 0
def findGrowth(_r0, _psi_gen, _step):
    _left, _right = findDuration(_r0, _psi_gen, _step)
    _tol = 1e-6
    _f_left, _f_right = g_func(_left, _r0, _psi_gen, _step), g_func(_right, _r0, _psi_gen, _step)
    if _f_left * _f_right > 0:
        return 
    while abs(_right - _left) > _tol:
        _mid = (_left + _right) / 2
        _f_mid = g_func(_mid, _r0, _psi_gen, _step)    
        if _f_mid == 0:
            break
        elif _f_mid> 0:
            _left = _mid
        else:
            _right = _mid
    return _mid




srv_funcs = {'Weibull': weibull_srv_func,
             'gamma': gamma_srv_func,
             'log-normal': lognormal_srv_func} # a dictionary maping the survival functions

means = {'Weibull': weibull_mean, 'gamma': gamma_mean, 'log-normal': lognormal_mean} # a dictionary maping the functions of calulating the distribution mean values


##################################### the parameters in this area can be customized by users ####################################################################
psi_length = 400                                                                # the time point number
step = 0.01                                                                     # the time step
alpha = {'inf': 1.5, 'rem':2}                                                     # the parameter alpha of infection and removal time distributions
beta = {'inf': 3, 'rem':2}                                                      # the parameter beta of infection and removal time distributions
lambda_max = 1                                                                  # the maximum eigenvalue 
func_forms = {'inf': 'Weibull', 'rem': 'Weibull'}                               # the distribution form of infection and removal processes

#################################################################################################################################################


#################################### the codes in this area generate results###########################################################################
# generate the distributions, hazard and survival functions and their mean values
distributions, hazards, survivals, mean_vals = {}, {}, {}, {}
for i in ['inf', 'rem']:
    distributions[i] = dist(srv_funcs[func_forms[i]], alpha[i], beta[i],psi_length, step)
    hazards[i] = haz(srv_funcs[func_forms[i]], alpha[i], beta[i],psi_length, step)
    survivals[i] = srv(srv_funcs[func_forms[i]], alpha[i], beta[i],psi_length, step)
    mean_vals[i] = means[func_forms[i]](alpha[i], beta[i])
     
# generate the genration time distribution and its mean value, and the value of lambda_eff
dist_gen, mean_gen, lambda_eff = gen(srv_funcs[func_forms['inf']], alpha['inf'], beta['inf'],
               srv_funcs[func_forms['rem']], alpha['rem'], beta['rem'], psi_length, step)

r0 = lambda_eff * lambda_max                                                    # the value of basic reproduction number R_0
g = findGrowth(r0, dist_gen, step)                                              # the value of growth rate g
mar_gamma, mar_mu = g*r0/(lambda_max*(r0-1)), g/(r0-1)                          # the infection and removal rates of the Markovian dynamics in transient-state equivalence

##########################################################################################################################################################################






# print some results in the console
print('average infection time T_inf: ', mean_vals['inf'])
print('average removal time T_rem: ', mean_vals['rem'])
print('average generation time T_gen: ', mean_gen)
print('basic reproduction number R_0: ', r0)
print('growth rate g: ', g)


#################################### the codes in this area plot figures###########################################################################
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

fig = plt.figure(figsize=(12, 9))    
axes = []
gs = GridSpec(6, 6, figure=fig)
for i in range(2):
    axes.append([])
    for j in range(3):
        axes[-1].append(fig.add_subplot(gs[j*2:(j+1)*2, i*2:(i+1)*2]))
axes.append([fig.add_subplot(gs[1:3, 4:6])])
axes.append([fig.add_subplot(gs[3:5, 4:6])])

colors = {'inf':'tab:red', 'rem':'tab:blue'}
titles = {'inf':'Infection', 'rem': 'Removal'}
for i, idx in enumerate(['inf', 'rem']):
    plt.sca(axes[i][0])
    plt.plot(np.arange(len(distributions[idx]))*step, distributions[idx], color = colors[idx])
    plt.title(titles[idx] + ' Time Distribution: ' + r'$\psi_{\mathrm{' + idx + r'}}(\tau)$')
    plt.axvline(mean_vals[idx], linestyle = '--', color = 'tab:gray')
    plt.sca(axes[i][1])
    plt.plot(np.arange(len(hazards[idx]))*step, hazards[idx], color = colors[idx])
    plt.title(titles[idx] + ' Hazard Function: ' + r'$\omega_{\mathrm{' + idx + r'}}(\tau)$')
    plt.sca(axes[i][2])
    plt.plot(np.arange(len(survivals[idx]))*step, survivals[idx], color = colors[idx])
    plt.title(titles[idx] + ' Survival Function: ' + r'$\Psi_{\mathrm{' + idx + r'}}(\tau)$')
        
plt.sca(axes[2][0])
plt.plot(np.arange(len(dist_gen))*step, dist_gen, color = 'tab:orange')
plt.axvline(mean_gen, linestyle = '--', color = 'tab:gray')
plt.title('Genration Time Distribution: ' + r'$\Psi_{\mathrm{gen}}(\tau)$')
    
plt.sca(axes[3][0])
plt.bar([r'$\gamma$', r'$\mu$'], [mar_gamma, mar_mu], color = ['tab:blue', 'tab:orange'], width = 0.6)
plt.title('Parameters of Markovian dynamics \n in Transient-state Equivalence')
plt.ylabel('Value')

for i in range(3):
    for ax in axes[i]:
        plt.sca(ax)
        plt.xlim(0 - psi_length * step / 25, psi_length * step)
        plt.xlabel(r'$\tau$')
        plt.ylabel('Value')
    
plt.subplots_adjust(wspace=2, hspace=3)
#########################################################################################################################
    
    
    
    
    
    









