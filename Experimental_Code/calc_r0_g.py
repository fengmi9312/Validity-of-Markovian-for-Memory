# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 12:51:16 2023

@author: admin
"""

import numpy as np
import func
from occur import occur, gnr, calc_lambda
from model_n import calc_general_k, calc_k
import data_import as di
import metapopulation_simulation as ms
from scipy.special import gamma, erf, gammainc
from scipy.optimize import minimize, Bounds

def get_mean_from_weibull(alpha, beta):
    return beta * gamma(1 + 1.0 / alpha)

def get_beta_from_weibull(alpha, mean_value):
    return mean_value / gamma(1 + 1.0 / alpha)

def get_mean_from_gamma(alpha, beta):
    return alpha * beta

def get_beta_from_gamma(alpha, mean_value):
    return mean_value / alpha

def get_mean_from_lognormal(alpha, beta):
    return np.exp(beta + (alpha ** 2) / 2)

def get_beta_from_lognormal(alpha, mean_value):
    return np.log(mean_value) - (alpha ** 2) / 2

def lognormal_survivals(_alpha, _beta, _length, _step):
    _tau = np.arange(1, _length) * _step
    return np.append(1, 0.5 - 0.5 * erf((np.log(_tau) - _beta) / (_alpha * (2 ** 0.5))))

def gamma_survivals(_alpha, _beta, _length, _step):
    _tau = _tau = np.arange(_length) * _step
    return 1 - gammainc(_alpha, _tau / _beta)


def get_lambda_from_cum(_alpha_inf, _beta_inf, _alpha_rem, _beta_rem):
    return (_alpha_inf / _alpha_rem) * ((_beta_rem / _beta_inf) ** _alpha_inf) * gamma(_alpha_inf / _alpha_rem)



def find_r(_occur_inf, _occur_rem, _r0):
    def _loss(_x):
        return (1 / _r0 - (np.exp(- _x[0] * np.arange(_occur_inf.get_len()) * _occur_inf.get_step()) * 
                       _occur_inf.get_rate() * _occur_rem.get_srv()[:-1]).sum() / calc_lambda(_occur_inf, _occur_rem)) ** 2
    return minimize(_loss, [0,], bounds =  Bounds(-np.ones(1) * np.inf, np.ones(1) * np.inf), method = 'SLSQP', tol = 1e-16)


day_div = 24
delay = 14
eta = 0.95
step = 1.0 / day_div
i0 = 0.01
time0 = 0
group_div = np.array([2,] * 8)
group_amount = len(group_div)
simu_node_amount = 12000
occur_length = 4000
steady_prdt = 0.6
amount_seeds = (np.ones(len(group_div)) * 30).astype(int)
flevel = 0.2
params_fitting = {'method': 'SLSQP', 'tol': 1e-16}

# obtain the parameters for simulations and calculations
params = {'contacts': di.import_contacts('United States of America'), 
          'populations': di.import_populations('United States of America')[0], 
          'ifrs': di.import_ifrs(), 
          'ylls': di.import_ylls('United States of America')}
params_adj = di.rescale_params(params, group_div)
contacts_dir = params_adj['contacts']['home'] + params_adj['contacts']['school'] \
             + params_adj['contacts']['work'] + params_adj['contacts']['other_locations']
contacts = contacts_dir + contacts_dir.T
simu_structure_params = {'node_amount':simu_node_amount, 'group_amount': group_amount, 'populations':params_adj['populations'], 
                       'contacts': contacts, 'ifrs': params_adj['ifrs'], 'ylls': params_adj['ylls'], 'k': 1, 'step': step, 'time': 0}
simu_vaccination_params = {'delay': delay, 'eta': eta}

simu_test = ms.simu(**simu_structure_params)
simu_test.set_spreading_func('weibull', 'weibull')
simu_test.set_total_spreading_params([1, 1], [1, 1])
simu_test.set_amount_seeds(amount_seeds)
populations_from_simu = np.array(simu_test.get_populations())
i0_from_simu = np.array(simu_test.get_i_amount_arr()) / np.array(simu_test.get_population_amounts())
calc_params = {'i0': i0_from_simu, 'populations': populations_from_simu, 'contacts': params_adj['contacts'], 
               'ifrs': params_adj['ifrs'], 'ylls': params_adj['ylls'], 
               'contacts_prop': None, 'k': simu_structure_params['k'], 'occur_inf': None, 'occur_rem': None,
               'delay': delay, 'eta': eta, 'step': step, 'time': time0} 
vgnr_inf = gnr(func.srv_weibull_scale, 'func', [1, 1])
vgnr_rem = gnr(func.srv_weibull_scale, 'func', [1, 1])
calc_params['occur_inf'] = occur(vgnr_inf, vgnr_type = 'srv', length = occur_length - 1, step = calc_params['step'])
calc_params['occur_rem'] = occur(vgnr_rem, vgnr_type = 'srv', length = occur_length - 1, step = calc_params['step'])
population_amounts = simu_test.get_population_amounts()
del simu_test

calc_selected = 7
print("calc_selected: " + str(calc_selected))

trans_list = [['weibull', 5, 7, 6], 
              ['weibull', 5, 5, 6],
              ['weibull', 7, 7, 6],
              ['lognormal', 5, 7, 6],
              ['gamma', 5, 7, 6],
              ['weibull', 5, 7, 4],
              ['weibull', 5, 7, 8],]
param_range = {'weibull': np.arange(-6, 25, 1), 'gamma': np.arange(-6, 25, 1), 'lognormal': np.arange(6, 37, 1)}
srv_funcs = {'weibull': func.srv_weibull_scale, 'gamma': func.srv_gamma_scale, 'lognormal': func.srv_lognormal}


func_type, mean_inf, mean_rem, steady_level = trans_list[calc_selected]


ks = []
get_beta = {'weibull': get_beta_from_weibull, 'gamma': get_beta_from_gamma, 'lognormal': get_beta_from_lognormal}
for i in param_range[func_type]:
    for j in param_range[func_type]:
        if func_type == 'weibull':
            alpha_inf = np.exp(i*0.05)
            alpha_rem = np.exp(j*0.05)
        elif func_type == 'lognormal':
            alpha_inf = i*0.05
            alpha_rem = j*0.05
        elif func_type == 'gamma':
            alpha_inf = np.exp(i*0.1)
            alpha_rem = np.exp(j*0.1)
        beta_inf = get_beta[func_type](alpha_inf, mean_inf)
        beta_rem = get_beta[func_type](alpha_rem, mean_rem)
        if func_type == 'weibull':
            ks.append(calc_k(alpha_inf, beta_inf, alpha_rem, beta_rem, 4000, 0.6, calc_params))
        else:
            vgnr_inf = gnr(srv_funcs[func_type], 'func', [alpha_inf, beta_inf])
            vgnr_rem = gnr(srv_funcs[func_type], 'func', [alpha_rem, beta_rem])
            calc_params['occur_inf'] = occur(vgnr_inf, vgnr_type = 'srv', length = occur_length - 1, step = calc_params['step'])
            calc_params['occur_rem'] = occur(vgnr_rem, vgnr_type = 'srv', length = occur_length - 1, step = calc_params['step'])
            ks.append(calc_general_k(calc_params['occur_inf'].get_srv(), calc_params['occur_rem'].get_srv(), steady_prdt, calc_params))


r0_g_data = {'r0': [], 'g': []}
if steady_level == 6:
    key = func_type + '_i' + str(mean_inf) + 'r' + str(mean_rem)
else:
    key = func_type + '_i' + str(mean_inf) + 'r' + str(mean_rem) + '_level_' + str(steady_level)
idx = 0
for i in param_range[func_type]:
    for j in param_range[func_type]:
        r = ks[idx] * 3.9707627520110265
        if func_type == 'weibull':
            alpha_inf = np.exp(i*0.05)
            alpha_rem = np.exp(j*0.05)
        elif func_type == 'lognormal':
            alpha_inf = i*0.05
            alpha_rem = j*0.05
        elif func_type == 'gamma':
            alpha_inf = np.exp(i*0.1)
            alpha_rem = np.exp(j*0.1)
        else:
            pass
        beta_inf = get_beta[func_type](alpha_inf, mean_inf)
        beta_rem = get_beta[func_type](alpha_rem, mean_rem)
        vgnr_inf = gnr(srv_funcs[func_type], 'func', [alpha_inf, beta_inf])
        vgnr_rem = gnr(srv_funcs[func_type], 'func', [alpha_rem, beta_rem])
        occur_inf = occur(vgnr_inf, vgnr_type = 'srv', length = occur_length - 1, step = calc_params['step'])
        occur_rem = occur(vgnr_rem, vgnr_type = 'srv', length = occur_length - 1, step = calc_params['step'])
        if func_type == 'weibull':
            lambda_real = get_lambda_from_cum(alpha_inf, beta_inf, alpha_rem, beta_rem)
        else:
            lambda_real = calc_lambda(occur_inf, occur_rem)
        growth_real = find_r(occur_inf, occur_rem,  r * lambda_real).x[0]
        r0_g_data['r0'].append(r * lambda_real)
        r0_g_data['g'].append(growth_real)
        print(idx)
        idx += 1
        
import pandas as pd   
writer = pd.ExcelWriter("../Experimental_Data/r0_g_data/r0_g_calc_"+ key + ".xlsx")     
pd.DataFrame(r0_g_data).to_excel(writer, sheet_name = 'r0_g')
writer.close()
print("calc_selected: " + str(calc_selected))