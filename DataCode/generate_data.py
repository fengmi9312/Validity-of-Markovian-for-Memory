# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 15:18:13 2023

@author: admin
"""

import numpy as np
from disease_survivals import get_all_meantime, covid_19, sars, h1n1, smallpox
from scipy.optimize import minimize, Bounds
from scipy.special import gamma
import func
from occur import occur, gnr, calc_lambda
import string
from copy import deepcopy
from model_n import calc_general_k, calc_k, sir_model
import data_import as di
import metapopulation_simulation as ms
from scipy.stats import spearmanr


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





def arr_extent(_arr, _length):
    if len(_arr) < _length:
        return np.append(_arr, np.ones(_length - len(_arr)) * _arr[-1])
    else:
        return _arr


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

def lambda_eff_weibull(_alpha_inf, _beta_inf, _alpha_rem, _beta_rem):
    return ((_beta_rem/_beta_inf)**_alpha_inf)*(_alpha_inf/_alpha_rem)*gamma(_alpha_inf/_alpha_rem)

def get_mean_from_cum(_alpha_inf, _beta_inf, _alpha_rem, _beta_rem):
    return _beta_rem * gamma((_alpha_inf + 1) / _alpha_rem) / gamma(_alpha_inf / _alpha_rem)

def get_lambda_from_cum(_alpha_inf, _beta_inf, _alpha_rem, _beta_rem):
    return (_alpha_inf / _alpha_rem) * ((_beta_rem / _beta_inf) ** _alpha_inf) * gamma(_alpha_inf / _alpha_rem)

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

def find_r(_occur_inf, _occur_rem, _r0):
    def _loss(_x):
        return (1 / _r0 - (np.exp(- _x[0] * np.arange(_occur_inf.get_len()) * _occur_inf.get_step()) * 
                       _occur_inf.get_rate() * _occur_rem.get_srv()[:-1]).sum() / calc_lambda(_occur_inf, _occur_rem)) ** 2
    return minimize(_loss, [0,], bounds =  Bounds(-np.ones(1) * np.inf, np.ones(1) * np.inf), method = 'SLSQP', tol = 1e-16)


def get_func_ir_name(_func_type, _mean_inf, _mean_rem, _steady_level):
    if _steady_level == 6:
        return _func_type + '_i' + str(_mean_inf) + 'r' + str(_mean_rem)
    else:
        return _func_type + '_i' + str(_mean_inf) + 'r' + str(_mean_rem) + '_level' + str(_steady_level)



trans_list = [['weibull', 5, 7, 6], 
              ['weibull', 5, 5, 6],
              ['weibull', 7, 7, 6],
              ['lognormal', 5, 7, 6],
              ['gamma', 5, 7, 6],
              ['weibull', 5, 7, 4],
              ['weibull', 5, 7, 8],]
param_range = {'weibull': np.arange(-6, 25, 1), 'gamma': np.arange(-6, 25, 1), 'lognormal': np.arange(6, 37, 1)}
srv_funcs = {'weibull': func.srv_weibull_scale, 'gamma': func.srv_gamma_scale, 'lognormal': func.srv_lognormal}

get_mean = {'weibull': get_mean_from_weibull, 'gamma': get_mean_from_gamma, 'lognormal': get_mean_from_lognormal}
get_beta = {'weibull': get_beta_from_weibull, 'gamma': get_beta_from_gamma, 'lognormal': get_beta_from_lognormal}


def get_ks():
    ks = {}
    for func_type, mean_inf, mean_rem, steady_level in trans_list[:5]:
        func_ir_name = func_type + '_i' + str(mean_inf) + 'r' + str(mean_rem)
        ks[func_ir_name] = []
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
                    ks[func_type + '_i' + str(mean_inf) + 'r' + str(mean_rem)].append(calc_k(alpha_inf, beta_inf, alpha_rem, beta_rem, 4000, 0.6, calc_params))
                else:
                    vgnr_inf = gnr(srv_funcs[func_type], 'func', [alpha_inf, beta_inf])
                    vgnr_rem = gnr(srv_funcs[func_type], 'func', [alpha_rem, beta_rem])
                    occur_inf = occur(vgnr_inf, vgnr_type = 'srv', length = occur_length - 1, step = calc_params['step'])
                    occur_rem = occur(vgnr_rem, vgnr_type = 'srv', length = occur_length - 1, step = calc_params['step'])
                    ks[func_type + '_i' + str(mean_inf) + 'r' + str(mean_rem)].append(calc_general_k(occur_inf.get_srv(), occur_rem.get_srv(), 0.6, calc_params))
    print('ks calculated')
    return ks

def generate_figure_theories_data(_simu_data, _calc_data, _simu_data_m, _simu_data_nm, 
                                  _steady_simu_data, _steady_calc_data, _steady_simu_mar_data, _steady_calc_mar_data, _transient_data):    
    _figure_theories_data = {'s_curves_simu(subplot_a)': None, 'i_curves_simu(subplot_b)': None, 'r_curves_simu(subplot_c)': None, 
                            's_curve_calc(subplot_a)': None, 'i_curve_calc(subplot_b)': None, 'r_curve_calc(subplot_c)': None,
                            'r_curves_simu(subplot_d)': None,
                            'r_steady_state(subplot_e)': None,
                            'curves_calc(subplot_f)': None, 'curves_calc(subplot_g)': None, 'curves_calc(subplot_h)': None, 'curves_calc(subplot_i)': None,
                            'weibull_errors_i5r7(subplot_j)': None, 'weibull_errors_i5r5(subplot_j)': None, 'weibull_errors_i7r7(subplot_j)': None, 
                            'lognormal_errors_i5r7(subplot_j)': None, 'gamma_errors_i5r7(subplot_j)': None}
    
    # generate data of subplots a--c
    for i, target in enumerate(['s', 'i', 'r']):
        max_len = 0
        for j in range(100):
            arr_len = len(_simu_data['data_' + str(j)][target].to_numpy())
            if arr_len > max_len:
                max_len = arr_len
        key_name = target + '_curves_simu(subplot_' + string.ascii_lowercase[i] + ')'
        _figure_theories_data[key_name] = {'time': np.arange(max_len) * step}
        for j in range(100):
            _figure_theories_data[key_name]['curve_simu_' + str(j)] = \
            arr_extent(_simu_data['data_' + str(j)][target].to_numpy(), max_len)
        key_name = target + '_curve_calc(subplot_' + string.ascii_lowercase[i] + ')'
        _figure_theories_data[key_name] = {'time': np.arange(len(_calc_data['data_calc'][target])) * step}
        _figure_theories_data[key_name]['curve_calc'] = _calc_data['data_calc'][target]
    
    # generate data of subplot d
    max_len = 0
    for i in range(100):
        arr_len = max(len(_simu_data_nm['data_' + str(i)]['r'].to_numpy()), len(_simu_data_m['data_' + str(i)]['r'].to_numpy()))
        if arr_len > max_len:
            max_len = arr_len
    key_name = 'r_curves_simu(subplot_d)'
    _figure_theories_data[key_name] = {'time': np.arange(max_len) * step}
    for i in range(100):
        _figure_theories_data[key_name]['nm_' + str(i)] = arr_extent(_simu_data_nm['data_' + str(i)]['r'].to_numpy(), max_len)
        _figure_theories_data[key_name]['m_' + str(i)] = arr_extent(_simu_data_m['data_' + str(i)]['r'].to_numpy(), max_len)
    
    # generate data of subplot e
    eig_max = 2.3135629562408058
    key_name = 'r_steady_state(subplot_e)'
    _figure_theories_data[key_name] = {'r0': np.linspace(0.01, 2, 40) * eig_max}
    _steady_res_nm = []
    _steady_res_m = []
    for i in range(100):
        _steady_res_nm.append([])
        _steady_res_m.append([])
        for j in range(40):
            _steady_res_nm[-1].append(_steady_simu_data['data_' + str(j)]['c'][i])
            _steady_res_m[-1].append(_steady_simu_mar_data['data_' + str(j)]['c'][i])
    _figure_theories_data[key_name]['steady_state_simu_non_m'] = np.array(_steady_res_nm).mean(axis = 0)
    _figure_theories_data[key_name]['steady_state_simu_m'] = np.array(_steady_res_m).mean(axis = 0)
    _figure_theories_data[key_name]['steady_state_calc'] = _steady_calc_data['data']['c']
    
    # generate data of subplots f--i
    mean_inf, mean_rem = 5, 7
    example_params = [[np.exp(0.45), np.exp(0.45)], [np.exp(-0.3), np.exp(1.2)]]
    for i in range(2):
        alpha_inf, alpha_rem = example_params[i]
        beta_inf = get_beta_from_weibull(alpha_inf, mean_inf)
        beta_rem = get_beta_from_weibull(alpha_rem, mean_rem)
        vgnr_inf = gnr(func.srv_weibull_scale, 'func', [alpha_inf, beta_inf])
        vgnr_rem = gnr(func.srv_weibull_scale, 'func', [alpha_rem, beta_rem])
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
        key_name = 'curves_calc(subplot_' + string.ascii_lowercase[5 + i * 2] + ')'
        _figure_theories_data[key_name] = {'time': calc_once.get_time_line()}
        _figure_theories_data[key_name]['s_curve'] = calc_once.get_s_tot()
        _figure_theories_data[key_name]['r_curve'] = calc_once.get_r_tot()
        _figure_theories_data[key_name]['inferred_s_curve'] = \
        calc_params['populations'] @ (calc_once.get_s()[0][:, None] * \
                                      np.exp(- lambda_eff * calc_params['k'] * contacts @ (calc_params['populations'] * (np.array(calc_once.get_r()) - calc_once.get_r()[0])).T))
        key_name = 'curves_calc(subplot_' + string.ascii_lowercase[6 + i * 2] + ')'
        _figure_theories_data[key_name] = {'time': calc_once.get_time_line()}
        _figure_theories_data[key_name]['i_curve_nm'] = calc_once.get_i_tot()
        _figure_theories_data[key_name]['r_curve_nm'] = calc_once.get_r_tot()
        _figure_theories_data[key_name]['i_curve_m'] = m_calc_once.get_i_tot()
        _figure_theories_data[key_name]['r_curve_m'] = m_calc_once.get_r_tot()
    
    # generate data of subplot j
    for func_type, mean_inf, mean_rem, steady_level in trans_list[:5]:
        func_ir_name = get_func_ir_name(func_type, mean_inf, mean_rem, steady_level)
        ratio_list= _transient_data[func_ir_name]['result']['ratio'].to_numpy()
        key_name = 'errors_'+func_ir_name+'(subplot_j)'
        _figure_theories_data[key_name] = {'log_ratio': np.log(ratio_list)}
        _figure_theories_data[key_name]['loss'] = []
        for m in param_range[func_type]:
            for n in param_range[func_type]:
                res = _transient_data[func_ir_name]['curves_' + str(m) + '_' + str(n)]
                level_len = 0
                c_tot = res['non_i'].to_numpy() + res['non_r'].to_numpy()
                for j in c_tot:
                    level_len += 1
                    if j > (c_tot[-1] - c_tot[0]) * 0.5 + c_tot[0]:
                        break
                non_s, non_i, non_r, m_s, m_i, m_r = res['non_s'].to_numpy()[:level_len], res['non_i'].to_numpy()[:level_len], \
                res['non_r'].to_numpy()[:level_len],res['m_s'].to_numpy()[:level_len], res['m_i'].to_numpy()[:level_len], res['m_r'].to_numpy()[:level_len]
                xxx = 2
                yyy = 1 / xxx
                loss = ((abs(non_s - m_s) ** xxx).sum() / (non_s ** xxx).sum() + \
                       (abs(non_i - m_i) ** xxx).sum() / (non_i ** xxx).sum() + \
                       (abs(non_r - m_r) ** xxx).sum() / (non_r ** xxx).sum()) ** yyy
                _figure_theories_data[key_name]['loss'].append(loss)
                
    return _figure_theories_data

def generate_si_figure_percentile_data(_transient_data):
    _si_figure_percentile_data = {'weibull_errors_i5r7_p20(subplot_a)': None, 'weibull_errors_i5r5_p20(subplot_a)': None, 'weibull_errors_i7r7_p20(subplot_a)': None, 
                                  'lognormal_errors_i5r7_p20(subplot_a)': None, 'gamma_errors_i5r7_p20(subplot_a)': None,
                                  'weibull_errors_i5r7_p80(subplot_b)': None, 'weibull_errors_i5r5_p80(subplot_b)': None, 'weibull_errors_i7r7_p80(subplot_b)': None, 
                                  'lognormal_errors_i5r7_p80(subplot_b)': None, 'gamma_errors_i5r7_p80(subplot_b)': None,
                                  'weibull_errors_i5r7_s60_p50(subplot_c)': None, 'weibull_errors_i5r7_s30_p50(subplot_c)': None,'weibull_errors_i5r7_s80_p50(subplot_c)': None,
                                  'weibull_errors_i5r7_s60_p5(subplot_d)': None, 'weibull_errors_i5r7_s30_p5(subplot_d)': None, 'weibull_errors_i5r7_s80_p5(subplot_d)': None}
    param_range = {'weibull': np.arange(-6, 25, 1), 'gamma': np.arange(-6, 25, 1), 'lognormal': np.arange(6, 37, 1)}
    #generate data of suplots a--b
    for idx, percentile in enumerate([10, 90]):
        for i in range(5):
            func_type, mean_inf, mean_rem, steady_level = trans_list[i]
            func_ir_name = get_func_ir_name(func_type, mean_inf, mean_rem, steady_level)
            ratio_list= _transient_data[func_ir_name]['result']['ratio'].to_numpy()
            key_name = 'errors_'+ func_ir_name + '_p' + str(percentile) + '(subplot_'+string.ascii_lowercase[idx]+')'
            _si_figure_percentile_data[key_name] = {'log_ratio': np.log(ratio_list)}
            _si_figure_percentile_data[key_name]['loss'] = []
            for m in param_range[func_type]:
                for n in param_range[func_type]:
                    res = _transient_data[func_ir_name]['curves_' + str(m) + '_' + str(n)]
                    level_len = 0
                    c_tot = res['non_i'].to_numpy() + res['non_r'].to_numpy()
                    for j in c_tot:
                        level_len += 1
                        if j > (c_tot[-1] - c_tot[0]) * percentile / 100 + c_tot[0]:
                            break
                    non_s, non_i, non_r, m_s, m_i, m_r = res['non_s'].to_numpy()[:level_len], res['non_i'].to_numpy()[:level_len], \
                    res['non_r'].to_numpy()[:level_len],res['m_s'].to_numpy()[:level_len], res['m_i'].to_numpy()[:level_len], res['m_r'].to_numpy()[:level_len]
                    xxx = 2
                    yyy = 1 / xxx
                    loss = ((abs(non_s - m_s) ** xxx).sum() / (non_s ** xxx).sum() + \
                           (abs(non_i - m_i) ** xxx).sum() / (non_i ** xxx).sum() + \
                           (abs(non_r - m_r) ** xxx).sum() / (non_r ** xxx).sum()) ** yyy
                    _si_figure_percentile_data[key_name]['loss'].append(loss)
       
    #generate data of suplots c--d
    for idx, percentile in enumerate([50, 5]):
        for i in [0, 5, 6]:
            func_type, mean_inf, mean_rem, steady_level = trans_list[i]
            func_ir_name = get_func_ir_name(func_type, mean_inf, mean_rem, steady_level)
            ratio_list= _transient_data[func_ir_name]['result']['ratio'].to_numpy()
            key_name = 'errors_'+ func_ir_name + '_p' + str(percentile) + '(subplot_'+string.ascii_lowercase[idx + 2]+')'
            _si_figure_percentile_data[key_name] = {'log_ratio': np.log(ratio_list)}
            _si_figure_percentile_data[key_name]['loss'] = []
            for m in param_range[func_type]:
                for n in param_range[func_type]:
                    res = _transient_data[func_ir_name]['curves_' + str(m) + '_' + str(n)]
                    level_len = 0
                    c_tot = res['non_i'].to_numpy() + res['non_r'].to_numpy()
                    for j in c_tot:
                        level_len += 1
                        if j > (c_tot[-1] - c_tot[0]) * percentile / 100 + c_tot[0]:
                            break
                    non_s, non_i, non_r, m_s, m_i, m_r = res['non_s'].to_numpy()[:level_len], res['non_i'].to_numpy()[:level_len], \
                    res['non_r'].to_numpy()[:level_len],res['m_s'].to_numpy()[:level_len], res['m_i'].to_numpy()[:level_len], res['m_r'].to_numpy()[:level_len]
                    xxx = 2
                    yyy = 1 / xxx
                    loss = ((abs(non_s - m_s) ** xxx).sum() / (non_s ** xxx).sum() + \
                           (abs(non_i - m_i) ** xxx).sum() / (non_i ** xxx).sum() + \
                           (abs(non_r - m_r) ** xxx).sum() / (non_r ** xxx).sum()) ** yyy
                    _si_figure_percentile_data[key_name]['loss'].append(loss)
    return _si_figure_percentile_data

def generate_si_figure_markovian_time_data():
    mean_inf, mean_rem = 5, 7
    _si_figure_markovian_time = {'log_alpha': np.arange(-6, 25, 1) * 0.05, 'm_mean_inf': [], 'm_mean_rem': []}
    for i in np.arange(-6, 25, 1):
        alpha_inf, alpha_rem = np.exp(i*0.05), np.exp(i*0.05)
        beta_inf = get_beta_from_weibull(alpha_inf, mean_inf)
        beta_rem = get_beta_from_weibull(alpha_rem, mean_rem)
        vgnr_inf = gnr(func.srv_weibull_scale, 'func', [alpha_inf, beta_inf])
        vgnr_rem = gnr(func.srv_weibull_scale, 'func', [alpha_rem, beta_rem])
        occur_inf = occur(vgnr_inf, vgnr_type = 'srv', length = occur_length - 1, step = step)
        occur_rem = occur(vgnr_rem, vgnr_type = 'srv', length = occur_length - 1, step = step)
        lambda_eff = calc_lambda(occur_inf, occur_rem)
        r0 = 3.9707627520110265 * calc_k(alpha_inf, beta_inf, alpha_rem, beta_rem, 4000, 0.6, calc_params) * lambda_eff
        r = find_r(occur_inf, occur_rem, r0).x[0]
        _si_figure_markovian_time['m_mean_inf'].append((r0 - 1) / (r * lambda_eff))
        _si_figure_markovian_time['m_mean_rem'].append((r0 - 1) / r)
    return _si_figure_markovian_time

def generate_si_figure_error_percentile_data(_transient_data, _log_ratios):
    percentile_range = np.arange(5, 96)
    _si_figure_error_percentile_data = {'weibull(subplot_a)': {'percentile': percentile_range}, 
                                        'lognormal(subplot_b)': {'percentile': percentile_range}, 
                                        'gamma(subplot_c)':{'percentile': percentile_range}}
    
    for i, idx in enumerate([0, 3, 4]):
        func_type, mean_inf, mean_rem, steady_level = trans_list[idx]
        func_ir_name = get_func_ir_name(func_type, mean_inf, mean_rem, steady_level)
        param_idx = 0
        key_name = func_type+'(subplot_'+string.ascii_lowercase[i]+')'
        for i in param_range[func_type]:
            for j in param_range[func_type]:
                log_ratio = _log_ratios[func_ir_name][param_idx]
                param_idx += 1
                if abs(log_ratio) > 0.05:
                    continue
                percentile = _si_figure_error_percentile_data[key_name]['percentile'][0]
                _si_figure_error_percentile_data[key_name]['loss_' + str(i) + '_' + str(j)] = []
                res = _transient_data[func_ir_name]['curves_' + str(i) + '_' + str(j)]
                level_len = 0
                c_tot = res['non_i'].to_numpy() + res['non_r'].to_numpy()
                for m in c_tot:
                    level_len += 1
                    if m > (c_tot[-1] - c_tot[0]) * percentile / 100 + c_tot[0]:
                        non_s, non_i, non_r, m_s, m_i, m_r = res['non_s'].to_numpy()[:level_len], res['non_i'].to_numpy()[:level_len], \
                        res['non_r'].to_numpy()[:level_len],res['m_s'].to_numpy()[:level_len], res['m_i'].to_numpy()[:level_len], res['m_r'].to_numpy()[:level_len]
                        xxx = 2
                        yyy = 1 / xxx
                        loss = ((abs(non_s - m_s) ** xxx).sum() / (non_s ** xxx).sum() + \
                               (abs(non_i - m_i) ** xxx).sum() / (non_i ** xxx).sum() + \
                               (abs(non_r - m_r) ** xxx).sum() / (non_r ** xxx).sum()) ** yyy
                        _si_figure_error_percentile_data[key_name]['loss_' + str(i) + '_' + str(j)].append(loss)
                        if percentile == _si_figure_error_percentile_data[key_name]['percentile'][-1]:
                            break
                        percentile += 1
    return _si_figure_error_percentile_data
      



def generate_figure_analysis_data(_explain_data, _forecasting_data, _r0_g_data, _log_ratios, _ks):
    _figure_analysis_data = {'fitting_period(subplot_a)': {}}
    simu_curve_data = [[], [], []]
    calc_curve_data = [[], [], []]
    for i in range(3):
        simu_curve_data[i] = [_explain_data['data_curves_'+ str(j)]['simu_curve_' + str(i)]['c'].to_numpy() for j in range(100)]
        calc_curve_data[i] = [_explain_data['data_curves_'+ str(j)]['calc_curve_' + str(i)]['c_exponent'].to_numpy() for j in range(100)]
    for i in range(3):
        simu_len = max([len(simu_curve_data[i][j]) for j in range(100)])
        key_name = 'explain_data_simu_' + str(i) + '(subplot_a)'
        _figure_analysis_data[key_name] = {'time': np.arange(simu_len) * step}
        _figure_analysis_data['fitting_period(subplot_a)']['fitting_period_'+str(i)] = []
        _figure_analysis_data['fitting_period(subplot_a)']['c_val_'+str(i)] = []
        for j in range(100):
            t_fitting = 0
            _c_val = 0
            for c_t in simu_curve_data[i][j]:
                if c_t > (simu_curve_data[i][j][-1] - simu_curve_data[i][j][0]) * 0.2 + simu_curve_data[i][j][0]:
                    _c_val = c_t
                    break
                t_fitting += step
            _figure_analysis_data['fitting_period(subplot_a)']['fitting_period_'+str(i)].append(t_fitting)
            _figure_analysis_data['fitting_period(subplot_a)']['c_val_'+str(i)].append(_c_val)
            _figure_analysis_data[key_name]['curve_' + str(j)] = arr_extent(simu_curve_data[i][j], simu_len)
        key_name = 'explain_data_calc_' + str(i) + '(subplot_a)'
        calc_len = max([len(calc_curve_data[i][j]) for j in range(100)])
        _figure_analysis_data[key_name] = {'time': np.arange(calc_len) * step}
        _figure_analysis_data[key_name]['curve'] = np.array([arr_extent(calc_curve_data[i][j], calc_len) for j in range(100)]).mean(axis = 0)
    
    fitted_x, fitted_y = [], []
    for func_type, mean_inf, mean_rem, steady_level in trans_list[:5]:
        idx = 0
        func_ir_name = get_func_ir_name(func_type, mean_inf, mean_rem, steady_level)
        key_name = func_ir_name + '_r0'+'(subplot_b)'
        _figure_analysis_data[key_name] = {'log_ratio': _log_ratios[func_ir_name], 'hat_r0': [], 'r0': _r0_g_data[func_ir_name]['r0_g']['r0'].to_numpy()}
        for i in param_range[func_type]:
            for j in param_range[func_type]:
                r = _ks[func_ir_name][idx] * 3.9707627520110265
                exp_inf = _forecasting_data[func_ir_name]['data_' + str(i) + '_' + str(j)]['params']['exponent_inf'].to_numpy()[:100]
                exp_rem = _forecasting_data[func_ir_name]['data_' + str(i) + '_' + str(j)]['params']['exponent_rem'].to_numpy()[:100]
                lambda_fitting = (exp_inf / exp_rem).mean()
                _figure_analysis_data[key_name]['hat_r0'].append(r * lambda_fitting)
                fitted_x.append(_log_ratios[func_ir_name][idx])
                fitted_y.append(- np.log(np.log(_figure_analysis_data[key_name]['hat_r0'][-1]) / np.log(_figure_analysis_data[key_name]['r0'][-1])))
                idx += 1
                
        
    def fitted_func(_x, _a):
        return _a * _x
    from scipy.optimize import curve_fit
    key_name = 'curve_fitted'
    _figure_analysis_data[key_name] = {}
    _figure_analysis_data[key_name]['log_ratio']= np.linspace(min(fitted_x), max(fitted_x), 100)
    fitted_a = curve_fit(fitted_func, fitted_x, fitted_y)
    _figure_analysis_data[key_name]['res']= fitted_func(_figure_analysis_data[key_name]['log_ratio'], fitted_a[0][0])
    #print(fitted_a)
    return _figure_analysis_data

def generate_figure_forecasting_data(_curve_data, _forecasting_data, _real_fitting_data, _log_ratios):
    _figure_forecasting_data = {}
    _fitting_period_data = {'fitting_period':{}}
    data_types = ['simu', 'weibull', 'exponent']
    curve_names = {'simu': 'simulation', 'weibull': 'non-Markovian', 'exponent': 'Markovian'}
    for idx, i in enumerate(['-6_24', '9_9', '24_-6']):
        data_len = min([len(_curve_data[i]['simu_curve_' + str(j)]) for j in range(100)])
        key_name = 'curves_' + i + '(subplot_' + string.ascii_lowercase[idx] + ')'
        _figure_forecasting_data[key_name] = {'time': np.arange(data_len) * step}
        _fitting_period_data['fitting_period']['fitting_period_' + i] = []
        _fitting_period_data['fitting_period']['c_val_' + i] = []
        for data_type in data_types:
            if data_type == 'simu':
                _figure_forecasting_data[key_name][curve_names[data_type]] = \
                    np.array([_curve_data[i]['simu_curve_' + str(j)]['c'].to_numpy()[:data_len] for j in range(100)]).mean(axis = 0)
                _figure_forecasting_data[key_name][curve_names[data_type] + '_std'] = \
                    np.array([_curve_data[i]['simu_curve_' + str(j)]['c'].to_numpy()[:data_len] for j in range(100)]).std(axis = 0)
                for j in range(100):
                    t_fitting = 0
                    _c_val = 0
                    curve_init, curve_end = _curve_data[i]['simu_curve_' + str(j)]['c'].to_numpy()[0], _curve_data[i]['simu_curve_' + str(j)]['c'].to_numpy()[-1]
                    for c_t in _curve_data[i]['simu_curve_' + str(j)]['c'].to_numpy():
                        if c_t > (curve_end - curve_init) * 0.2 + curve_init:
                            _c_val = c_t
                            break
                        t_fitting += step
                    _fitting_period_data['fitting_period']['fitting_period_' + i].append(t_fitting)
                    _fitting_period_data['fitting_period']['c_val_' + i].append(_c_val)
            else:
                _figure_forecasting_data[key_name][curve_names[data_type]] = \
                    np.array([_curve_data[i]['calc_curve_' + str(j)]['c_' + data_type].to_numpy()[:data_len] for j in range(100)]).mean(axis = 0)
                _figure_forecasting_data[key_name][curve_names[data_type] + '_std'] = \
                    np.array([_curve_data[i]['calc_curve_' + str(j)]['c_' + data_type].to_numpy()[:data_len] for j in range(100)]).std(axis = 0)
    real_survivals = {'covid_19': covid_19, 'sars': sars, 'h1n1': h1n1, 'smallpox': smallpox}
    
    key_name = 'real_fitting(subplot_f)'
    _figure_forecasting_data[key_name] = {}
    for fitting_type in ['exponent', 'weibull']:
        for disease in ['covid_19', 'sars', 'h1n1', 'smallpox']: 
            x = get_all_meantime(real_survivals[disease](4000, 1 / 24), 1/ 24) 
            log_ratio = np.log(x[2] / x[1])
            c_tmp = _real_fitting_data[disease]['res']['c_' + fitting_type+'_steady'].to_numpy()
            c_tmp_tot = _real_fitting_data[disease]['res']['c_tot_steady'].to_numpy()
            res = (c_tmp - c_tmp_tot) / c_tmp_tot
            _figure_forecasting_data[key_name][disease + '_' + curve_names[fitting_type]] = [log_ratio, res.mean(), res.std()]
            
    
    for func_type, mean_inf, mean_rem, steady_level in trans_list[:5]:
        for idx, fitting_type in enumerate(['exponent', 'weibull']): 
            func_ir_name = get_func_ir_name(func_type, mean_inf, mean_rem, steady_level)
            key_name = func_ir_name + '_' + curve_names[fitting_type] + '(subplot_f)'
            _figure_forecasting_data[key_name] = {'log_ratio': _log_ratios[func_ir_name], 'res': []}
            if func_ir_name == 'weibull_i5r7':
                key_name_grid = 'grid_' + curve_names[fitting_type] + '(subplot_' + string.ascii_lowercase[idx + 3] + ')'
                _figure_forecasting_data[key_name_grid] = {}
            for i in param_range[func_type]:
                if func_ir_name == 'weibull_i5r7':
                    _figure_forecasting_data[key_name_grid]['inf_'+str(i)] = []
                for j in param_range[func_type]:
                    c_tmp = _forecasting_data[func_ir_name]['data_' + str(i) + '_' + str(j)]['res']['c_' + fitting_type+'_steady'].to_numpy()[:100]
                    c_tmp_tot = _forecasting_data[func_ir_name]['data_' + str(i) + '_' + str(j)]['res']['c_tot_steady'].to_numpy()[:100]
                    if func_ir_name == 'weibull_i5r7':
                        _figure_forecasting_data[key_name_grid]['inf_'+str(i)].append(((c_tmp - c_tmp_tot) / c_tmp_tot).mean())
                    _figure_forecasting_data[key_name]['res'].append(((c_tmp - c_tmp_tot) / c_tmp_tot).mean())
    return _figure_forecasting_data, _fitting_period_data
  
def generate_figure_prevention_data(_vac_curve_data, _vac_data, _real_vac_fitting_data, _log_ratios):
    _figure_prevention_data = {}
    curve_names = {'simu': 'simulation', 'weibull': 'non-Markovian', 'exponent': 'Markovian'}
    ages = ['under_20', '20-49', '20+', '60+', 'all_ages']
    target = 'c'
    file_mark = '9_9'
    key_name = 'simulation_curves(subplot_a)'
    _figure_prevention_data[key_name] = {}
    for i, age in enumerate(['no_vaccine'] + ages):
        min_len = min([len(_vac_curve_data[file_mark]['simu_curve_' + age + '_' + str(j)][target]) for j in range(100)])
        _figure_prevention_data[key_name]['time'] = np.arange(min_len) * step
        _figure_prevention_data[key_name][age] = \
        np.array([_vac_curve_data[file_mark]['simu_curve_' + age + '_' + str(j)][target].to_numpy()[:min_len] for j in range(100)]).mean(axis = 0)
        _figure_prevention_data[key_name][age + '_std'] = \
        np.array([_vac_curve_data[file_mark]['simu_curve_' + age + '_' + str(j)][target].to_numpy()[:min_len] for j in range(100)]).std(axis = 0)
    for idx, fitting_type in enumerate(['exponent', 'weibull']):
        key_name = curve_names[fitting_type] + '_curves(subplot_'+string.ascii_lowercase[idx + 1]+')' 
        _figure_prevention_data[key_name] = {'time': np.arange(min_len) * step}
        for i, age in enumerate(['no_vaccine'] + ages):
            min_len = min([len(_vac_curve_data[file_mark]['calc_curve_' + age + '_' + str(j)][fitting_type + '_' + target]) for j in range(100)])
            _figure_prevention_data[key_name][age] = \
            np.array([_vac_curve_data[file_mark]['calc_curve_' + age + '_' + str(j)][fitting_type + '_' + target].to_numpy()[:min_len] for j in range(100)]).mean(axis = 0)
            _figure_prevention_data[key_name][age + '_std'] = \
            np.array([_vac_curve_data[file_mark]['calc_curve_' + age + '_' + str(j)][fitting_type + '_' + target].to_numpy()[:min_len] for j in range(100)]).std(axis = 0)
    
    real_survivals = {'covid_19': covid_19, 'sars': sars, 'h1n1': h1n1, 'smallpox': smallpox}
    key_name = 'real_fitting(subplot_i)'
    key_name1 = 'real_fitting(subplot_j)'
    _figure_prevention_data[key_name] = {}
    _figure_prevention_data[key_name1] = {}
    for fitting_type in ['exponent', 'weibull']:
        for disease in ['covid_19', 'sars', 'h1n1', 'smallpox']: 
            x = get_all_meantime(real_survivals[disease](4000, 1 / 24), 1/ 24) 
            log_ratio = np.log(x[2] / x[1])
            _res = []
            _opt_res = []
            for m in range(100):
                res_simu = []
                res_calc = []
                _res_once = 0
                _res_tot = 0
                for age in ['no_vaccine'] + ages:
                    if age != 'no_vaccine':
                        res_simu.append(_real_vac_fitting_data[disease]['simu_data_' + age][target][m])
                        res_calc.append(_real_vac_fitting_data[disease]['calc_data_' + age][fitting_type + '_' + target][m])
                    _x_simu = _real_vac_fitting_data[disease]['simu_data_' + age][target][m]
                    _x_target = _real_vac_fitting_data[disease]['calc_data_' + age][fitting_type + '_' + target][m]
                    _res_once += (_x_simu - _x_target) ** 2
                    _res_tot += _x_simu ** 2
                _res.append((_res_once ** 0.5) / (_res_tot ** 0.5))
                if np.argmin(res_simu) == np.argmin(res_calc):
                    _opt_res.append(0)
                else:
                    _opt_res.append(1)
            _figure_prevention_data[key_name][disease + '_' + curve_names[fitting_type]] = [log_ratio, np.array(_res).mean(), np.array(_res).std()]
            _figure_prevention_data[key_name1][disease + '_' + curve_names[fitting_type]] = [log_ratio, np.array(_opt_res).mean()]
    
    
    for func_type, mean_inf, mean_rem, steady_level in trans_list[:5]:
        for idx, fitting_type in enumerate(['exponent', 'weibull']):
            func_ir_name = get_func_ir_name(func_type, mean_inf, mean_rem, steady_level)
            key_name = func_ir_name + '_' + fitting_type + '(subplot_i)'
            key_name1 = func_ir_name + '_' + fitting_type + '(subplot_j)'
            _figure_prevention_data[key_name] = {'log_ratio': _log_ratios[func_ir_name], 'res': []}
            _figure_prevention_data[key_name1] = {'log_ratio': _log_ratios[func_ir_name], 'res': []}
            if func_ir_name == 'weibull_i5r7':
                key_name_grid = 'grid_' + curve_names[fitting_type] +'(subplot_'+string.ascii_lowercase[idx + 4]+')'
                key_name_grid1 = 'grid_' + curve_names[fitting_type] +'(subplot_'+string.ascii_lowercase[idx + 6]+')'
                _figure_prevention_data[key_name_grid] = {}
                _figure_prevention_data[key_name_grid1] = {}
            for i, inf_pow in enumerate(param_range[func_type]):
                if func_ir_name == 'weibull_i5r7':
                    _figure_prevention_data[key_name_grid]['inf_' + str(inf_pow)] = []
                    _figure_prevention_data[key_name_grid1]['inf_' + str(inf_pow)] = []
                for j, rem_pow in enumerate(param_range[func_type]):
                    _data_one = _vac_data[func_ir_name]['data_'+str(inf_pow)+'_'+str(rem_pow)]
                    _res = 0
                    _opt_res = 0
                    for m in range(100):
                        res_simu = []
                        res_calc = []
                        _res_once = 0
                        _res_tot = 0
                        for age in ['no_vaccine'] + ages:
                            if age != 'no_vaccine':
                                res_simu.append(_data_one['simu_data_' + age][target][m])
                                res_calc.append(_data_one['calc_data_' + age][fitting_type + '_' + target][m])
                            _x_simu = _data_one['simu_data_' + age][target][m]
                            _x_target = _data_one['calc_data_' + age][fitting_type + '_' + target][m]
                            _res_once += (_x_simu - _x_target) ** 2
                            _res_tot += _x_simu ** 2
                        _res += (_res_once ** 0.5) / (_res_tot ** 0.5)
                        if np.argmin(res_simu) == np.argmin(res_calc):
                            _opt_res += 0
                        else:
                            _opt_res += 1
                    _figure_prevention_data[key_name]['res'].append(_res / 100)
                    _figure_prevention_data[key_name1]['res'].append(_opt_res / 100)
                    if func_ir_name == 'weibull_i5r7':
                        _figure_prevention_data[key_name_grid]['inf_' + str(inf_pow)].append(_res / 100)
                        _figure_prevention_data[key_name_grid1]['inf_' + str(inf_pow)].append(_opt_res / 100)
                        
        key_name = 'vac_eff(subplot_d)'
        _figure_prevention_data[key_name] = {}
        _data_one = _vac_data['weibull_i5r7']['data_24_24']
        _figure_prevention_data[key_name]['simulation'] = \
            np.array([[_data_one['simu_data_' + age][target][m] for age in (['no_vaccine'] + ages)] for m in range(100)]).mean(axis = 0)
        _figure_prevention_data[key_name]['simulation_std'] = \
            np.array([[_data_one['simu_data_' + age][target][m] for age in (['no_vaccine'] + ages)] for m in range(100)]).std(axis = 0)
        _figure_prevention_data[key_name]['Markovian'] = \
            np.array([[_data_one['calc_data_' + age]['exponent_' + target][m] for age in (['no_vaccine'] + ages)] for m in range(100)]).mean(axis = 0)
        _figure_prevention_data[key_name]['Markovian_std'] = \
            np.array([[_data_one['calc_data_' + age]['exponent_' + target][m] for age in (['no_vaccine'] + ages)] for m in range(100)]).std(axis = 0)
        _figure_prevention_data[key_name]['non-Markovian'] = \
            np.array([[_data_one['calc_data_' + age]['weibull_' + target][m] for age in (['no_vaccine'] + ages)] for m in range(100)]).mean(axis = 0)
        _figure_prevention_data[key_name]['non-Markovian_std'] = \
            np.array([[_data_one['calc_data_' + age]['weibull_' + target][m] for age in (['no_vaccine'] + ages)] for m in range(100)]).std(axis = 0)
    
    return _figure_prevention_data
   
def generate_figure_prevention_opt_data(_vac_curve_data, _vac_data, _real_vac_fitting_data, _log_ratios):
    _figure_prevention_opt_data = {}
    curve_names = {'simu': 'simulation', 'weibull': 'non-Markovian', 'exponent': 'Markovian'}
    ages = ['under_20', '20-49', '20+', '60+', 'all_ages']
    target = 'c'
    real_survivals = {'covid_19': covid_19, 'sars': sars, 'h1n1': h1n1, 'smallpox': smallpox}
    for age_idx in range(5):
        key_name = 'real_fitting(subplot_'+string.ascii_lowercase[age_idx]+')'
        _figure_prevention_opt_data[key_name] = {}
        for fitting_type in ['exponent', 'weibull']:
            for disease in ['covid_19', 'sars', 'h1n1', 'smallpox']: 
                x = get_all_meantime(real_survivals[disease](4000, 1 / 24), 1/ 24) 
                log_ratio = np.log(x[2] / x[1])
                _res = 0
                for m in range(100):
                    _res_simu = []
                    _res_calc = []
                    for age in ages:
                        if ages[age_idx] != age:
                            _res_simu.append(_real_vac_fitting_data[disease]['simu_data_' + age][target][m])
                            _res_calc.append(_real_vac_fitting_data[disease]['calc_data_' + age][fitting_type + '_' + target][m])
                    if np.argmin(_res_simu) != np.argmin(_res_calc):
                        _res += 1
                _figure_prevention_opt_data[key_name][disease + '_' + curve_names[fitting_type]] = [log_ratio, _res / 100]
        for func_type, mean_inf, mean_rem, steady_level in trans_list[:5]:
            for idx, fitting_type in enumerate(['exponent', 'weibull']):
                func_ir_name = get_func_ir_name(func_type, mean_inf, mean_rem, steady_level)
                key_name = func_ir_name + '_' + fitting_type + '(subplot_'+string.ascii_lowercase[age_idx]+')'
                _figure_prevention_opt_data[key_name] = {'log_ratio': _log_ratios[func_ir_name], 'res': []}
                for i, inf_pow in enumerate(param_range[func_type]):
                    for j, rem_pow in enumerate(param_range[func_type]):
                        _data_one = _vac_data[func_ir_name]['data_'+str(inf_pow)+'_'+str(rem_pow)]
                        _res = 0
                        for m in range(100):
                            _res_simu = []
                            _res_calc = []
                            for age in ages:
                                if ages[age_idx] != age:
                                    _res_simu.append(_data_one['simu_data_' + age][target][m])
                                    _res_calc.append(_data_one['calc_data_' + age][fitting_type + '_' + target][m])
                            if np.argmin(_res_simu) != np.argmin(_res_calc):
                                _res += 1
                        _figure_prevention_opt_data[key_name]['res'].append(_res / 100)
    return _figure_prevention_opt_data    
   
