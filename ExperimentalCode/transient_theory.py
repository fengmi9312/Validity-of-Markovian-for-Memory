# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 13:47:17 2023

@author: admin
"""



import data_import as di
from model_n import calc_general_k, sir_model
import metapopulation_simulation as ms
from scipy.optimize import minimize, Bounds
from scipy.special import gamma
import func
from occur import occur, gnr, calc_lambda
from copy import deepcopy
import pandas as pd
import numpy as np

def get_mean_from_weibull(alpha, beta):
    return beta * gamma(1 + 1.0 / alpha)

def get_beta_from_weibull(alpha, mean_value):
    return mean_value / gamma(1 + 1.0 / alpha)

def lambda_eff_weibull(_alpha_inf, _beta_inf, _alpha_rem, _beta_rem):
    return ((_beta_rem/_beta_inf)**_alpha_inf)*(_alpha_inf/_alpha_rem)*gamma(_alpha_inf/_alpha_rem)

def get_mean_from_cum(_alpha_inf, _beta_inf, _alpha_rem, _beta_rem):
    return _beta_rem * gamma((_alpha_inf + 1) / _alpha_rem) / gamma(_alpha_inf / _alpha_rem)

def split_int(x, a):
    tmp = [a] * (x // a)
    sum_tmp = sum(tmp)
    if sum_tmp != x:
        tmp.append(x - sum_tmp)
    return np.array(tmp)

def get_mean_from_gamma(alpha, beta):
    return alpha * beta

def get_beta_from_gamma(alpha, mean_value):
    return mean_value / alpha

def get_mean_from_lognormal(alpha, beta):
    return np.exp(beta + (alpha ** 2) / 2)

def get_beta_from_lognormal(alpha, mean_value):
    return np.log(mean_value) - (alpha ** 2) / 2

def calc_generation_time(_occur_inf, _occur_rem):
    if _occur_inf.get_step() != _occur_rem.get_step():
        print('The steps of a and b must be equal.')
        return None
    __length = max(_occur_inf.get_actl_len(), _occur_rem.get_actl_len()) * 2
    x = occur(_occur_inf.get_vgnr(), _occur_inf.get_vgnr_type(), __length, _occur_inf.get_step())
    y = occur(_occur_rem.get_vgnr(), _occur_rem.get_vgnr_type(), __length, _occur_rem.get_step())
    z = x.get_rate() * y.get_srv()[:-1]
    return (np.arange(len(z)) * x.get_step() * z).sum() / z.sum()

def calc_mean_time(_occur):
    return (np.arange(_occur.get_actl_len()) * _occur.get_step() * _occur.get_dist()[:_occur.get_actl_len()]).sum()
##########################################################

# common parameters
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
amount_seeds = (np.ones(len(group_div)) * 30).astype(int)

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

def find_r(_occur_inf, _occur_rem, _r0):
    def _loss(_x):
        return (1 / _r0 - (np.exp(- _x[0] * np.arange(_occur_inf.get_len()) * _occur_inf.get_step()) * 
                       _occur_inf.get_rate() * _occur_rem.get_srv()[:-1]).sum() / calc_lambda(_occur_inf, _occur_rem)) ** 2
    return minimize(_loss, [0,], bounds =  Bounds(-np.ones(1) * np.inf, np.ones(1) * np.inf), method = 'SLSQP', tol = 1e-16)

def get_m(_gamma_inf, _mu_rem, _calc_params, _spread_len):
    m_calc_params = deepcopy(_calc_params)
    m_vgnr_inf = gnr(func.srv_exponent, 'func', [_gamma_inf,])
    m_vgnr_rem = gnr(func.srv_exponent, 'func', [_mu_rem,])
    m_calc_params['occur_inf'] = occur(m_vgnr_inf, vgnr_type = 'srv', length = occur_length - 1, step = m_calc_params['step'])
    m_calc_params['occur_rem'] = occur(m_vgnr_rem, vgnr_type = 'srv', length = occur_length - 1, step = m_calc_params['step'])
    m_calc_once = sir_model(**m_calc_params)
    for i in range(_spread_len):
        m_calc_once.spread_once()
    return m_calc_once.get_s_tot(), m_calc_once.get_i_tot(), m_calc_once.get_r_tot()



def loss(_params, _calc_params, _calc_once):
    _res = get_m(_params[0], _params[1], _calc_params, _calc_once.get_current_time_int())
    return ((_res[0] - _calc_once.get_s_tot()) ** 2 + (_res[1] - _calc_once.get_i_tot()) ** 2 + (_res[2] - _calc_once.get_r_tot()) ** 2).sum() / ((_calc_once.get_current_time_int() + 1) * step)

param_range = {'weibull': np.arange(-6, 25, 1), 'gamma': np.arange(-6, 25, 1), 'lognormal': np.arange(6, 37, 1)}
get_mean = {'weibull': get_mean_from_weibull, 'gamma': get_mean_from_gamma, 'lognormal': get_mean_from_lognormal}
get_beta = {'weibull': get_beta_from_weibull, 'gamma': get_beta_from_gamma, 'lognormal': get_beta_from_lognormal}
srv_funcs = {'weibull': func.srv_weibull_scale, 'gamma': func.srv_gamma_scale, 'lognormal': func.srv_lognormal}

res= {'ratio': []}

trans_list = [['weibull', 5, 7, 6], 
              ['weibull', 5, 5, 6],
              ['weibull', 7, 7, 6],
              ['gamma', 5, 7, 6],
              ['lognormal', 5, 7, 6],
              ['weibull', 5, 7, 4],
              ['weibull', 5, 7, 8],]

if __name__ == '__main__':
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    func_type, mean_inf, mean_rem, steady_level = trans_list[rank]
    steady_prdt = steady_level / 10
    
    writer = pd.ExcelWriter("trans_equivalence_data/"+ func_type + "_" + str(mean_inf) + "_" + str(mean_rem) + "_" + str(steady_level) + ".xlsx")
    for i in param_range[func_type]:
        for j in param_range[func_type]:
            #print(i, j)
            if func_type == 'lognormal':
                alpha_inf = i * 0.05
                alpha_rem = j * 0.05
            else:
                alpha_inf = np.exp(i * 0.05)
                alpha_rem = np.exp(j * 0.05)
            beta_inf = get_beta[func_type](alpha_inf, mean_inf)
            beta_rem = get_beta[func_type](alpha_rem, mean_rem)
            vgnr_inf = gnr(srv_funcs[func_type], 'func', [alpha_inf, beta_inf])
            vgnr_rem = gnr(srv_funcs[func_type], 'func', [alpha_rem, beta_rem])
            calc_params['occur_inf'] = occur(vgnr_inf, vgnr_type = 'srv', length = occur_length - 1, step = calc_params['step'])
            calc_params['occur_rem'] = occur(vgnr_rem, vgnr_type = 'srv', length = occur_length - 1, step = calc_params['step'])
            calc_params['k'] = calc_general_k(calc_params['occur_inf'].get_srv(), calc_params['occur_rem'].get_srv(), steady_prdt, calc_params)
            calc_once = sir_model(**calc_params)
            
            lambda_eff = calc_lambda(calc_params['occur_inf'], calc_params['occur_rem'])
            r0 = 3.9707627520110265 * calc_params['k'] * lambda_eff
            r = find_r(calc_params['occur_inf'], calc_params['occur_rem'], r0).x[0]
            m_calc_params = deepcopy(calc_params)
            m_vgnr_inf = gnr(func.srv_exponent, 'func', [(r * lambda_eff) / (r0 - 1),])
            m_vgnr_rem = gnr(func.srv_exponent, 'func', [r / (r0 - 1),])
            m_calc_params['occur_inf'] = occur(m_vgnr_inf, vgnr_type = 'srv', length = occur_length - 1, step = m_calc_params['step'])
            m_calc_params['occur_rem'] = occur(m_vgnr_rem, vgnr_type = 'srv', length = occur_length - 1, step = m_calc_params['step'])
            m_calc_once = sir_model(**m_calc_params)
            while True:
                calc_once.spread_once()
                m_calc_once.spread_once() 
                if calc_once.getc_i_tot() < 1e-4:
                    break
            if func_type == 'lognormal':
                alpha_inf, alpha_rem = i * 0.05, j * 0.05
            else:
                alpha_inf = np.exp(i * 0.05)
                alpha_rem = np.exp(j * 0.05)
            beta_inf = get_beta[func_type](alpha_inf, mean_inf)
            beta_rem = get_beta[func_type](alpha_rem, mean_rem)
            if func_type == 'weibull':
                res['ratio'].append(get_mean_from_cum(alpha_inf, beta_inf, alpha_rem, beta_rem) / get_mean_from_weibull(alpha_rem, beta_rem))
            else:
                vgnr_inf = gnr(srv_funcs[func_type], 'func', [alpha_inf, beta_inf])
                vgnr_rem = gnr(srv_funcs[func_type], 'func', [alpha_rem, beta_rem])
                occur_inf = occur(vgnr_inf, vgnr_type = 'srv', length = occur_length - 1, step = calc_params['step'])
                occur_rem = occur(vgnr_rem, vgnr_type = 'srv', length = occur_length - 1, step = calc_params['step'])
                res['ratio'].append(calc_generation_time(occur_inf, occur_rem)/calc_mean_time(occur_rem)) 
            curve_data = {}
            curve_data['non_s'] = calc_once.get_s_tot()
            curve_data['non_i'] = calc_once.get_i_tot()
            curve_data['non_r'] = calc_once.get_r_tot()
            curve_data['m_s'] = m_calc_once.get_s_tot()
            curve_data['m_i'] = m_calc_once.get_i_tot()
            curve_data['m_r'] = m_calc_once.get_r_tot()
            pd.DataFrame(curve_data).to_excel(writer, sheet_name = 'curves_' + str(i) + '_' + str(j))
     
    pd.DataFrame(res).to_excel(writer, sheet_name = 'result')
    writer.close() 
                
'''                
    
    
    
colors =  {'weibull': 'tab:red', 'gamma': 'tab:blue', 'lognormal': 'tab:orange'}
for func_type in func_types:
    plt.plot(x[func_type], loss_calc[func_type], '.', color = colors[func_type])
                
             
                
                
                
                
                
            res = minimize(loss, x0 = np.array([1,1]), args = (calc_params, calc_once), bounds =  Bounds(np.zeros(2), np.ones(2) * np.inf), 
                         jac = None, options={'ftol': 1e-16, 'disp': True, 'maxiter':100000}, **params_fitting)
            results[func_type].append(res.fun)
            param_res[func_type].append(res.x)
            m_res = get_m(res.x[0], res.x[1], calc_params, calc_once.get_current_time_int())
            re[func_type].append((abs(m_res[0] - calc_once.get_s_tot()) + abs(m_res[1] - calc_once.get_i_tot()) + abs(m_res[2] - calc_once.get_r_tot())).sum() / 
                                      ((calc_once.get_current_time_int() + 1) * step))
            n_calc[func_type].append([calc_once.get_time_line(), calc_once.get_s_tot(), calc_once.get_i_tot(), calc_once.get_r_tot()])
            m_calc[func_type].append([np.arange(len(m_res[0])) * step, m_res[0], m_res[1], m_res[2]])
    
for i in n_calc['weibull']:
    plt.plot(i[0], i[2], color = 'tab:red')
    plt.plot(i[0], i[3], color = 'tab:green')
        
for i in m_calc['weibull']:
    plt.plot(i[0], i[2], color = 'tab:blue')
    plt.plot(i[0], i[3], color = 'tab:orange')
        

x = []
func_type = 'lognormal'
for i in param_range['lognormal']:
    for j in param_range['lognormal']:
        alpha_inf = i * 0.05
        alpha_rem = j * 0.05
        beta_inf = get_beta['lognormal'](alpha_inf, mean_inf)
        beta_rem = get_beta['lognormal'](alpha_rem, mean_rem)
        vgnr_inf = gnr(srv_funcs[func_type], 'func', [alpha_inf, beta_inf])
        vgnr_rem = gnr(srv_funcs[func_type], 'func', [alpha_rem, beta_rem])
        occur_inf = occur(vgnr_inf, vgnr_type = 'srv', length = occur_length - 1, step = calc_params['step'])
        occur_rem = occur(vgnr_rem, vgnr_type = 'srv', length = occur_length - 1, step = calc_params['step'])
        x.append(calc_generation_time(occur_inf, occur_rem)/calc_mean_time(occur_rem))

'''













