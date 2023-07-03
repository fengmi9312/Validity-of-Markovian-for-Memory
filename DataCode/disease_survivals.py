# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 14:48:21 2022

@author: 20481756
"""
import numpy as np
from scipy.special import gammainc, erf

def from_survival_to_hazard(_survivals):
    _hazards = []
    for i in range(1, len(_survivals)):
        if _survivals[i] == 0:
            _hazards.append(1)
        else:
            _hazards.append(1 - _survivals[i] / _survivals[i - 1])
    return _hazards

def from_survival_to_dist(_survival):
    return _survival[:-1] - _survival[1:]

def from_dist_to_survival(_dist):
    return np.append(1, 1 - _dist.cumsum())

def from_hazard_to_survival(_hazard):
    return np.append(1, (1 - _hazard).cumprod())

def sum_survival(_survival0, _survival1):
    return from_dist_to_survival(np.convolve(from_survival_to_dist( _survival0), 
                                             from_survival_to_dist(_survival1), 'full')[:min(len(_survival0), len(_survival1)) - 1])

def get_gentime(_survivals_inf, _survivals_rem, _time_step):
    _hazards_inf = from_survival_to_hazard(_survivals_inf)
    _res = 0
    _cum = 0
    for i in range(len(_survivals_rem) - 1):
        _res += _hazards_inf[i] * _survivals_rem[i] * i * _time_step
        _cum += _hazards_inf[i] * _survivals_rem[i]
    return _res / _cum

def get_r0(_survivals_inf, _survivals_rem, _time_step):
    _hazards_inf = from_survival_to_hazard(_survivals_inf)
    _cum = 0
    for i in range(len(_survivals_rem) - 1):
        _cum += _hazards_inf[i] * _survivals_rem[i]
    return _cum

def get_meantime(_survivals, _time_step):
    _res = 0
    _tot = 0
    _dist = from_survival_to_dist(_survivals)
    for i in range(len(_dist)):
        _res += _dist[i] * i * _time_step
        _tot += _dist[i]
    return _res / _tot

def get_all_meantime(_survival_arrs, _time_step):
    return get_meantime(_survival_arrs[0], _time_step), get_meantime(_survival_arrs[1], _time_step), get_gentime(_survival_arrs[0], _survival_arrs[1],_time_step)


def covid_19(_length, _step, _beta0 = 5.665):
    _tau = np.arange(_length) * _step
    
    #_mean_incu, _sigma_incu = np.exp(1.621), np.exp(0.418)
    _shape_incu, _scale_incu = 5.807, 0.948
    _shape_asym, _scale_asym = 5.0, 1.0
    #_shape_pres, _scale_pres = 1.058, 2.174
    _shape_symp, _scale_symp = 2.768, 1.1563
    #_miu_incu = np.log(_mean_incu) - (_sigma_incu ** 2) / 2
    #_survivals_incu = np.append(1, 0.5 - 0.5 * erf((np.log(_tau[1:]) - _miu_incu) / (_sigma_incu * (2 ** 0.5))))
    _survivals_incu = 1 - gammainc(_shape_incu, _tau / _scale_incu)
    _survivals_asym = 1 - gammainc(_shape_asym, _tau / _scale_asym)
    #_survivals_pres = 1 - gammainc(_shape_pres, _tau / _scale_pres)
    _survivals_symp = 1 - gammainc(_shape_symp, _tau / _scale_symp)
    _survivals_rem = sum_survival(_survivals_symp, _survivals_incu) * 0.8 + _survivals_asym * 0.2
    
    _alpha, _beta, _r0 = 2.826, _beta0, 2.0
    _hazards_inf = _r0 * (np.exp(- (_tau / _beta) ** _alpha) - np.exp(- ((_tau + _step) / _beta) ** _alpha)) / _survivals_rem
    _survivals_inf = from_hazard_to_survival(_hazards_inf)
    return _survivals_inf, _survivals_rem

def sars(_length, _step, _miu = 0.25):
    _tau = np.arange(_length) * _step
    _r0 = 3.0
    _shape_incu, _scale_incu = 2.4, 2.6
    _shape_symp0, _scale_symp0 = 8.9, 2.6
    _shape_symp1, _scale_symp1 = 1.9, 2.5
    _survivals_incu = 1 - gammainc(_shape_incu, _tau / _scale_incu)
    _survivals_symp0 = 1 - gammainc(_shape_symp0, _tau / _scale_symp0)
    _survivals_symp1 = 1 - gammainc(_shape_symp1, _tau / _scale_symp1)
    _survivals_symp = sum_survival(_survivals_symp0, _survivals_symp1)
    _survivals_rem = sum_survival(_survivals_incu, _survivals_symp)
    _c_inf = _r0 / ((_survivals_incu * 0.1 + (_survivals_rem - _survivals_incu)) * 0.25 * _step).sum()
    _hazards_inf = _c_inf * (_survivals_incu * 0.1 + (_survivals_rem - _survivals_incu)) * _miu * _step / _survivals_rem
    _survivals_inf = from_hazard_to_survival(_hazards_inf)[:-1]
    return _survivals_inf, _survivals_rem

def h1n1(_length, _step, _r0 = 1.31):
    _tau = np.arange(_length) * _step
    _median_incu, _mean_incu = 4, 4.3
    _median_symp, _mean_symp = 7, 9.3
    _miu_incu = np.log(_median_incu)
    _sigma_incu = ((np.log(_mean_incu) - _miu_incu) * 2) ** 0.5
    _miu_symp = np.log(_median_symp)
    _sigma_symp = ((np.log(_mean_symp) - _miu_symp) * 2) ** 0.5
    _survivals_incu = np.append(1, 0.5 - 0.5 * erf((np.log(_tau[1:]) - _miu_incu) / (_sigma_incu * (2 ** 0.5))))
    _survivals_symp = np.append(1, 0.5 - 0.5 * erf((np.log(_tau[1:]) - _miu_symp) / (_sigma_symp * (2 ** 0.5))))
    _survivals_rem = sum_survival(_survivals_incu, _survivals_symp)
    _lambda_c = 0.25184078285552103
    x = (_survivals_rem - _survivals_incu).sum() / _lambda_c
    _c_inf = _r0 / x
    _hazards_inf = _c_inf * (_survivals_rem - _survivals_incu) / (_lambda_c * _survivals_rem)
    _survivals_inf = from_hazard_to_survival(_hazards_inf)[:-1]
    return _survivals_inf, _survivals_rem

def smallpox(_length, _step, _b = 0.157):
    _tau = np.arange(_length) * _step
    _mean_incu, _sd_incu = 11.6, 1.9
    _mean_fever, _sd_fever = 2.49, 0.88
    _mean_rash, _sd_rash = 16, 2.83
    
    _shape_incu, _scale_incu = (_mean_incu ** 2) / (_sd_incu ** 2), (_sd_incu ** 2) / _mean_incu
    _shape_fever, _scale_fever = (_mean_fever ** 2) / (_sd_fever ** 2), (_sd_fever ** 2) / _mean_fever
    _shape_rash, _scale_rash = (_mean_rash ** 2) / (_sd_rash ** 2), (_sd_rash ** 2) / _mean_rash
    _survivals_incu = 1 - gammainc(_shape_incu, _tau / _scale_incu)
    _survivals_fever = 1 - gammainc(_shape_fever, _tau / _scale_fever)
    _survivals_rash = 1 - gammainc(_shape_rash, _tau / _scale_rash)
    _survivals_incu_fever = sum_survival(_survivals_incu, _survivals_fever)
    _survivals_rem = sum_survival(_survivals_incu_fever, _survivals_rash)
    _lambda_c = 0.25184078285552103
    _r0 = 6.87
    x = (_b * (_survivals_incu_fever - _survivals_incu) + (_survivals_rem - _survivals_incu_fever)).sum() / _lambda_c
    _c_inf = _r0 / x
    _hazards_inf = _c_inf * (_b * (_survivals_incu_fever - _survivals_incu) + (_survivals_rem - _survivals_incu_fever)) / (_lambda_c * _survivals_rem)
    _survivals_inf = from_hazard_to_survival(_hazards_inf)[:-1]
    return _survivals_inf, _survivals_rem
'''
import matplotlib.pyplot as plt


step = 1 / 24
x = sars(4000, step)
print(get_all_meantime(x, step))
tau = np.arange(4000 - 1) * step
alpha = 2.5
beta = 4.128083910048251 / gamma(1+1/alpha)
alpha = 4.655165822
beta = 38.01153843
plt.plot(tau, (np.exp(-(tau/beta)**alpha) - np.exp(-((tau+step)/beta)**alpha)) / step, color = 'red')
#plt.plot(tau, from_survival_to_dist(x[0]) / step, color = 'blue', linestyle = '-')
plt.plot(tau, from_survival_to_dist(x[1]) / step, color = 'blue', linestyle = '--')
plt.xlim(0, 80)



print(get_all_meantime(x, step)[2] / get_all_meantime(x, step)[1])
x = sars(4000, step)
print(get_all_meantime(x, step))
plt.plot(np.arange(4000 - 1) * step, from_survival_to_dist(x[0]), color = 'red', linestyle = '-')
plt.plot(np.arange(4000 - 1) * step, from_survival_to_dist(x[1]), color = 'red', linestyle = '--')
print(get_all_meantime(x, step)[2] / get_all_meantime(x, step)[1])
x = h1n1(4000, step)
print(get_all_meantime(x, step))
plt.plot(np.arange(4000 - 1) * step, from_survival_to_dist(x[0]), color = 'orange', linestyle = '-')
plt.plot(np.arange(4000 - 1) * step, from_survival_to_dist(x[1]), color = 'orange', linestyle = '--')
print(get_all_meantime(x, step)[2] / get_all_meantime(x, step)[1])
x = smallpox(4000, step)
print(get_all_meantime(x, step))
plt.plot(np.arange(4000 - 1) * step, from_survival_to_dist(x[0]), color = 'green', linestyle = '-')
plt.plot(np.arange(4000 - 1) * step, from_survival_to_dist(x[1]), color = 'green', linestyle = '--')
print(get_all_meantime(x, step)[2] / get_all_meantime(x, step)[1])
'''
'''

step = 1 / 24
x = covid_19(4000, step)
print(from_survival_to_dist(x[0]).sum())
x = sars(4000, step)
print(from_survival_to_dist(x[0]).sum())
x = h1n1(4000, step)
print(from_survival_to_dist(x[0]).sum())
x = smallpox(4000, step)
print(from_survival_to_dist(x[0]).sum())
'''






