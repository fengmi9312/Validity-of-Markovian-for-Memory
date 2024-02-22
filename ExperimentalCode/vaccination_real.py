# -*- coding: utf-8 -*-
"""
Created on Fri May 20 14:30:47 2022

@author: 20481756
"""
import sys
sys.path.append('../Dependencies/CodeDependencies')

# If you are running the code on Linux, you can utilize the following code that is compatible with Slurm and mpi4py.
sys.path.append('../Dependencies/ModuleDependencies(Linux)')
'''
from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
'''
# If you are running the code on Windows, you can utilize the following 2 lines of code.
sys.path.append('../Dependencies/ModuleDependencies(Windowns)')
rank = 0 # The rank should be limited to the range(4)


if rank >= 4:
    sys.exit()
from scipy.optimize import Bounds
import numpy as np
import data_import as di
from model_n import sir_model
import metapopulation_simulation as ms
from scipy.special import gamma
import pandas as pd
from occur import occur, gnr
import disease_survivals as ds
from copy import deepcopy
import func

def get_fitting_with_real(_idx, _spreading_params, _k, _amount_seeds, 
             simu_structure_params, calc_params, params_fitting, 
             fitting_level= 0.1, occur_length = 4000):
    
    ################ initialize the simulation model ##########################
    _simu_structure_params = deepcopy(simu_structure_params)
    _simu_structure_params['k'] = _k
    simu_once = ms.simu(**_simu_structure_params)
    simu_once.set_generator_seed(_idx)
    simu_once.set_spreading_func('general', 'general')
    simu_once.set_total_spreading_params([_spreading_params['step_inf'],], [_spreading_params['step_rem'],])
    simu_once.set_total_spreading_survivals(_spreading_params['srv_inf'], _spreading_params['srv_rem'])
    simu_once.set_amount_seeds(_amount_seeds)
    ###########################################################################
    
    #################### get the simulation data ##############################
    _simu_once_data = {'s': [], 'i': [], 'r': [], 'w': [], 'c': [], 'd': [], 'y': []}
    _simu_once_arr_data = {'c': []}
    for _target in _simu_once_data.keys():
        _simu_once_data[_target].append(simu_once.get_x_amount(_target) / simu_once.get_node_amount())
    _simu_once_arr_data['c'].append(np.array(simu_once.get_x_amount_arr('c')) / np.array(simu_once.get_population_amounts()))
    while True:
        if simu_once.get_i_amount() == 0:
            break
        simu_once.spread_once()
        for _target in _simu_once_data.keys():
            _simu_once_data[_target].append(simu_once.get_x_amount(_target) / simu_once.get_node_amount())
        _simu_once_arr_data['c'].append(np.array(simu_once.get_x_amount_arr('c')) / np.array(simu_once.get_population_amounts()))
    for _target in _simu_once_data.keys():
        _simu_once_data[_target] = np.array(_simu_once_data[_target])
    _simu_once_arr_data['c'] = np.array(_simu_once_arr_data['c'])
    del simu_once
    ###########################################################################
    
    ##################### get the simulation data for fitting #################
    _fitting_data = {'confirmed': [], 'removal': [], 'confirmed_arr': []}
    _fitting_end_c = (_simu_once_data['c'][-1] - _simu_once_data['c'][0]) * fitting_level + _simu_once_data['c'][0]
    for _i in range(len(_simu_once_data['c'])):
        _fitting_data['confirmed'].append(_simu_once_data['c'][_i])
        _fitting_data['removal'].append(_simu_once_data['r'][_i])
        _fitting_data['confirmed_arr'].append(_simu_once_arr_data['c'][_i])
        if _fitting_data['confirmed'][-1] > _fitting_end_c:
            break
    _fitting_data['confirmed'] = np.array(_fitting_data['confirmed'])
    _fitting_data['removal'] = np.array(_fitting_data['removal'])
    _fitting_data['confirmed_arr'] = np.array(_fitting_data['confirmed_arr'])
    _fitting_length = len(_fitting_data['confirmed'])
    _i_in_data = np.append(_fitting_data['confirmed_arr'][0:1], 
                           np.diff(_fitting_data['confirmed_arr'], axis = 0), axis = 0)
    ###########################################################################
    
    ################## initialize fitting model ###############################
    _calc_params = deepcopy(calc_params)
    _calc_params['real_data'] = _fitting_data
    _calc_params['occur_inf'] = None
    _calc_params['occur_rem'] = None
    _calc_params['k'] = _k
    fitting_once = [sir_model(**_calc_params) for i in range(2)]
    _init_params_rem = [np.array([1.0, ]), np.array([1.0, 1.0])]
    _init_params_inf = [np.array([1.0, ]), np.array([1.0, 1.0])]
    _fitting_func = ['exponent', 'weibull_rate']
    _bounds = [Bounds(np.ones(1) * 0.001, np.ones(1) * np.inf),
               Bounds(np.ones(2) * 0.001, np.ones(2) * np.inf)]
    _rem_fitting_res = [None, None]
    _inf_fitting_res = [None, None]
    ###########################################################################
    
    ######################## fitting ##########################################
    for i in range(2):
        _init_params = _init_params_rem[i]
        while True:
            _rem_fitting_res[i] =\
            fitting_once[i].fit_rem_from_data(init_params = _init_params, 
                                              occur_length = occur_length,  
                                              fitting_func = _fitting_func[i], 
                                              bounds = _bounds[i], **params_fitting)
            if _rem_fitting_res[i].success:
                break
            else:
                _init_params = _init_params_rem[i] * np.random.uniform(0.7, 1.3) 
        _init_params = _init_params_inf[i]
        while True:     
            _inf_fitting_res[i] = \
            fitting_once[i].fit_inf_from_data(init_params = _init_params, 
                                              occur_length = occur_length, 
                                              fitting_func = _fitting_func[i], 
                                              bounds = _bounds[i],  **params_fitting)
            if _inf_fitting_res[i].success:
                break
            else:
                _init_params = _init_params_inf[i] * np.random.uniform(0.7, 1.3)
    ###########################################################################
    
    ###################### return data ###############################        
    _params = {'exponent_inf': _inf_fitting_res[0].x[0], 'exponent_rem': _rem_fitting_res[0].x[0], 
               'weibull_alpha_inf':_inf_fitting_res[1].x[0], 'weibull_beta_inf':1.0 / _inf_fitting_res[1].x[1], 
               'weibull_alpha_rem':_rem_fitting_res[1].x[0], 'weibull_beta_rem':1.0 / _rem_fitting_res[1].x[1],
               'fitting_length': _fitting_length}
    return {'params':_params, 'i_in_data': _i_in_data, 'simu_curves': _simu_once_data}

def get_v_with_real(_idx, _simu_spreading_params, _fitting_spreading_params, spreading_len,
               _i_in_data, _k, _amount_seeds, 
               simu_structure_params, calc_params, 
               daily_vac_amount, vac_times = 7, vac_duration = 24, occur_length = 4000):
    
    def _get_alloc(s_arr, vac_amount, vac_group = None, adjust = True):
        if vac_group is None or len(vac_group) == 0:
            return None
        if adjust:
            vac_alloc = np.zeros(len(s_arr), dtype = np.int64)
        else:
            vac_alloc = np.zeros(len(s_arr))
        if s_arr[vac_group].sum() == 0:
            proportion = np.ones(len(vac_group)) / len(vac_group)
        else:
            proportion = s_arr[vac_group] / s_arr[vac_group].sum()
        if adjust:
            vac_alloc_tmp = (proportion * vac_amount).astype(np.int32)
            for i in np.argsort(proportion * vac_amount - vac_alloc_tmp)[::-1]:
                vac_alloc_tmp[i] += 1
                if vac_alloc_tmp.sum() == vac_amount:
                    break
            vac_alloc[vac_group] = vac_alloc_tmp
        else:
            vac_alloc[vac_group] = proportion * vac_amount
        return vac_alloc
    
    vac_groups = {'no_vaccine': None, 'under_20': np.arange(0, 2).tolist(), 
                  '20-49': np.arange(2, 5).tolist(), 
                  '20+': np.arange(2, 8).tolist(), '60+': np.arange(6, 8).tolist(), 
                  'all_ages': np.arange(0, 8).tolist()}
    vac_dates = [len(_i_in_data) + i * vac_duration for i in range(vac_times)]
    spreading_length = max(vac_dates[-1] + 30 * 24, spreading_len)
    
    _simu_structure_params = deepcopy(simu_structure_params)
    _simu_structure_params['k'] = _k
    simu_once = ms.simu(**_simu_structure_params)
    _node_amount = simu_once.get_node_amount()
    del simu_once
    simu_data = {'no_vaccine': None, 'under_20': None, '20-49': None, '20+': None,'60+': None, 'all_ages': None}
    for vac_group_key in vac_groups.keys():
        simu_data[vac_group_key] = {'s': [], 'i': [], 'r': [], 'w': [], 'c': [], 'd': [], 'y': []}
        _simu_structure_params = deepcopy(simu_structure_params)
        _simu_structure_params['k'] = _k
        simu_once = ms.simu(**_simu_structure_params)
        simu_once.set_generator_seed(_idx)
        simu_once.set_spreading_func('general', 'general')
        simu_once.set_total_spreading_params([_simu_spreading_params['step_inf'],], [_simu_spreading_params['step_rem'],])
        simu_once.set_total_spreading_survivals(_simu_spreading_params['srv_inf'], _simu_spreading_params['srv_rem'])
        simu_once.set_amount_seeds(_amount_seeds)
        for _target in simu_data[vac_group_key].keys():
            simu_data[vac_group_key][_target].append(simu_once.get_x_amount(_target) / simu_once.get_node_amount())
        for i in range(spreading_length - 1):
            if i in vac_dates:
                vac_alloc_simu = func.get_alloc(np.array(simu_once.get_s_amount_arr()), _simu_structure_params['populations'], daily_vac_amount, vac_groups[vac_group_key])
                if vac_alloc_simu is not None:
                    simu_once.add_vaccine(vac_alloc_simu)
            simu_once.spread_once()
            for _target in simu_data[vac_group_key].keys():
                simu_data[vac_group_key][_target].append(simu_once.get_x_amount(_target) / simu_once.get_node_amount())
        for _target in simu_data[vac_group_key].keys():
            simu_data[vac_group_key][_target] = np.array(simu_data[vac_group_key][_target])
        del simu_once
    
    calc_data = {'no_vaccine': {}, 'under_20': {}, '20-49': {}, '20+': {}, '60+': {}, 'all_ages': {}}
    for vac_group_key in vac_groups.keys():
        calc_once = []
        _srv_func = [func.srv_exponent, func.srv_weibull_scale, func.srv_gamma_scale]
        _inf_params = [(_fitting_spreading_params['exponent_inf'],),
                      (_fitting_spreading_params['weibull_alpha_inf'], _fitting_spreading_params['weibull_beta_inf'])]
        _rem_params = [(_fitting_spreading_params['exponent_rem'],),
                      (_fitting_spreading_params['weibull_alpha_rem'], _fitting_spreading_params['weibull_beta_rem'])]
        for i in range(2):
            _calc_params = deepcopy(calc_params)
            vgnr_inf = gnr(_srv_func[i], 'func', _inf_params[i])
            vgnr_rem = gnr(_srv_func[i], 'func', _rem_params[i])
            _calc_params['occur_inf'] = occur(vgnr_inf, vgnr_type = 'srv', length = occur_length, step = _calc_params['step'])
            _calc_params['occur_rem'] = occur(vgnr_rem, vgnr_type = 'srv', length = occur_length, step = _calc_params['step'])
            _calc_params['k'] = _k
            calc_once.append(sir_model(**_calc_params))
            calc_once[i].set_i_in_init(_i_in_data)
            for j in range(len(_i_in_data) - 1, spreading_length - 1):
                if j in vac_dates:
                    vac_alloc_calc = func.get_alloc(np.array(calc_once[i].getc_s_eff()), _calc_params['populations'], daily_vac_amount / _node_amount, vac_groups[vac_group_key], 
                                                adjust = False)
                    if vac_alloc_calc is not None:
                        calc_once[i].add_vaccine(vac_alloc_calc / calc_once[i].get_populations())
                calc_once[i].spread_once()
        for _target in ['s', 'i', 'r', 'c', 'd', 'y']:
            calc_data[vac_group_key]['exponent_' + _target] = calc_once[0].get_x_tot(_target)
            calc_data[vac_group_key]['weibull_' + _target] = calc_once[1].get_x_tot(_target)
        calc_data[vac_group_key]['exponent_w'] = calc_data[vac_group_key]['exponent_r'] - calc_data[vac_group_key]['exponent_d']
        calc_data[vac_group_key]['weibull_w'] = calc_data[vac_group_key]['weibull_r'] - calc_data[vac_group_key]['weibull_d']
        del calc_once
    
    return {'simu_data': simu_data, 'calc_data': calc_data}



np.random.seed(0)

# self-defined important function
def get_mean_from_weibull(alpha, beta):
    return beta * gamma(1 + 1.0 / alpha)

def get_beta_from_weibull(alpha, mean_value):
    return mean_value / gamma(1 + 1.0 / alpha)


def split_int(x, a):
    tmp = [a] * (x // a)
    sum_tmp = sum(tmp)
    if sum_tmp != x:
        tmp.append(x - sum_tmp)
    return np.array(tmp)

##########################################################

# common parameters
day_div = 24
delay = day_div * 3
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
params_fitting = {'method': 'L-BFGS-B', 'tol': 1e-16}
survival_funcs = {'covid_19': ds.covid_19, 'sars': ds.sars, 'h1n1': ds.h1n1, 'smallpox': ds.smallpox}
disease = ['covid_19', 'sars', 'h1n1', 'smallpox']


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
                       'contacts': contacts, 'ifrs': params_adj['ifrs'], 'ylls': params_adj['ylls'], 'k': 0.3, 'step': step, 'time': 0}
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
del simu_test


if __name__ == '__main__':
    res = {}
    param_list = {}
    fitting_data = {}
    real_survivals = survival_funcs[disease[rank]](occur_length, step)
    simu_once = ms.simu(**simu_structure_params)
    simu_once.set_spreading_func('general', 'general')
    simu_once.set_total_spreading_params([step,], [step,])
    simu_once.set_total_spreading_survivals(real_survivals[0], real_survivals[1])
    simu_once.set_amount_seeds(amount_seeds)
    populations_from_simu = np.array(simu_once.get_populations())
    i0_from_simu = np.array(simu_once.get_i_amount_arr()) / np.array(simu_once.get_population_amounts())
    del simu_once
    calc_params = {'i0': i0_from_simu, 'populations': populations_from_simu, 'contacts': params_adj['contacts'], 
                   'ifrs': params_adj['ifrs'], 'ylls': params_adj['ylls'], 
                   'contacts_prop': None, 'k': 1, 'occur_inf': None, 'occur_rem': None,
                   'delay': delay, 'eta': eta, 'step': step, 'time': time0} 
    _vgnr_inf = gnr(real_survivals[0], 'arr', None)
    _vgnr_rem = gnr(real_survivals[1], 'arr', None)
    calc_params['occur_inf'] = occur(_vgnr_inf, vgnr_type = 'srv', length = occur_length - 1, step = calc_params['step'])
    calc_params['occur_rem'] = occur(_vgnr_rem, vgnr_type = 'srv', length = occur_length - 1, step = calc_params['step'])
    k_calc_obtain = sir_model(**calc_params)
    k = k_calc_obtain.get_k_from_steady(steady_prdt, target = 'c', prdt_method = 'calc', tol = 1e-10)
    steady_data = {}
    writer = pd.ExcelWriter("../ExperimentalData/vaccination_data_real/data_"+ disease[rank] + ".xlsx")
    for i in range(100):
        simu_once = ms.simu(**simu_structure_params)  
        spreading_params = {'step_inf': step, 'step_rem': step, 'srv_inf': real_survivals[0], 'srv_rem': real_survivals[1]}
        fitting_res = get_fitting_with_real(i + 101, spreading_params, k, amount_seeds, 
                 simu_structure_params, calc_params, params_fitting, flevel)
        calc_res = get_v_with_real(i + 101, spreading_params, fitting_res['params'], 
                              len(fitting_res['simu_curves']['c']), fitting_res['i_in_data'], k, 
                              amount_seeds, simu_structure_params, calc_params, 300)
        for _key in calc_res.keys():
            for _key1 in calc_res[_key].keys():
                if i == 0:
                    steady_data[_key + '_' + _key1] = {}
                for _key2 in calc_res[_key][_key1]:
                    if i == 0:
                        steady_data[_key + '_' + _key1][_key2] =  calc_res[_key][_key1][_key2][-1:].tolist()
                    else:
                        steady_data[_key + '_' + _key1][_key2].append(calc_res[_key][_key1][_key2][-1])
    for _key in calc_res.keys():
        for _key1 in calc_res[_key].keys():
                pd.DataFrame(steady_data[_key + '_' + _key1]).to_excel(writer, sheet_name = _key + '_' + _key1)
    writer.close()         
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

