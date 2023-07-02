# -*- coding: utf-8 -*-
"""
Created on Fri May 20 14:30:47 2022

@author: 20481756
"""

import numpy as np
import data_import as di
import metapopulation_simulation as ms
from model_n import get_fitting_data_with_real, get_calc_data, sir_model
from occur import occur, gnr
import disease_survivals as ds


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
params_fitting = {'method': 'L-BFGS-B', 'tol': 1e-16}
survival_funcs = {'covid_19': ds.covid_19, 'sars': ds.sars, 'h1n1': ds.h1n1, 'smallpox': ds.smallpox}
disease = ['covid_19', 'sars', 'h1n1', 'smallpox']

##########################################################

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
                       'contacts': contacts, 'ifrs': params_adj['ifrs'], 'ylls': params_adj['ylls'], 'k': 0.6, 'step': step, 'time': 0}

if __name__ == '__main__':
    res = {}
    param_list = {}
    import pandas as pd
    rank = 0
    
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
    writer = pd.ExcelWriter("../Experimental_Data/real_fitting_data/data_"+ disease[rank] + ".xlsx")
    calc_params = {'i0': i0_from_simu, 'populations': populations_from_simu, 'contacts': params_adj['contacts'], 
                   'ifrs': params_adj['ifrs'], 'ylls': params_adj['ylls'], 
                   'contacts_prop': None, 'k': 1, 'occur_inf': None, 'occur_rem': None,
                   'delay': delay, 'eta': eta, 'step': step, 'time': time0} 
    _vgnr_inf = gnr(real_survivals[0], 'arr', None)
    _vgnr_rem = gnr(real_survivals[1], 'arr', None)
    calc_params['occur_inf'] = occur(_vgnr_inf, vgnr_type = 'srv', length = occur_length - 1, step = calc_params['step'])
    calc_params['occur_rem'] = occur(_vgnr_rem, vgnr_type = 'srv', length = occur_length - 1, step = calc_params['step'])
    k_calc_obtain = sir_model(**calc_params)
    simu_structure_params['k'] = k_calc_obtain.get_k_from_steady(steady_prdt, target = 'c', prdt_method = 'calc', tol = 1e-10)
    calc_params['k'] = simu_structure_params['k']
    for i in range(100):
        print(i)
        simu_once = ms.simu(**simu_structure_params)
        simu_once.set_generator_seed(i + 1)
        simu_once.set_spreading_func('general', 'general')
        simu_once.set_total_spreading_params([step,], [step,])
        simu_once.set_total_spreading_survivals(real_survivals[0], real_survivals[1])
        simu_once.set_amount_seeds(amount_seeds)
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
        _fitting_end_c = (_simu_once_data['c'][-1] - _simu_once_data['c'][0]) * flevel + _simu_once_data['c'][0]
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
        fitting_res = get_fitting_data_with_real(_fitting_data, calc_params, params_fitting)
        if i == 0:
            for key in fitting_res.keys():
                param_list[key] = [fitting_res[key]]
        else:
            for key in fitting_res.keys():
                param_list[key].append(fitting_res[key])
        calc_res = get_calc_data(fitting_res, _i_in_data, simu_structure_params['k'], _simu_once_data, calc_params)
        if i == 0:
            for res_key in calc_res['res'].keys():
                res[res_key] = [calc_res['res'][res_key], ]
        else:
            for res_key in calc_res['res'].keys():
                res[res_key].append(calc_res['res'][res_key])
    pd.DataFrame(param_list).to_excel(writer, sheet_name = 'params')   
    pd.DataFrame(res).to_excel(writer, sheet_name = 'res')    
    writer.close()
    


'''
params_res = get_fitting_data_with_real(fitting_data, calc_params, params_fitting)

inf_params = [(params_res['exponent_inf'],),
              (params_res['weibull_alpha_inf'], params_res['weibull_beta_inf']),
              (params_res['gamma_alpha_inf'], params_res['gamma_beta_inf'])]
rem_params = [(params_res['exponent_rem'],),
              (params_res['weibull_alpha_rem'], params_res['weibull_beta_rem']),
              (params_res['gamma_alpha_rem'], params_res['gamma_beta_rem'])]
srv_func = [func.srv_exponent, func.srv_weibull_scale, func.srv_gamma_scale]
nmfuncs = ['exponent', 'weibull', 'gamma']
rem_curves = {}
for idx, nmfunc in enumerate(nmfuncs):
    i_in = np.append(np.array(real_curves['USA'][0][0]), np.diff(np.array(real_curves['USA'][0])))
    _data_len = len(i_in)
    rem_srv = np.array([1 - srv_func[idx](j * calc_params['step'], *rem_params[idx]) for j in range(_data_len)])
    rem_curves[nmfunc] = np.convolve(rem_srv, i_in, 'full')[:_data_len] + (np.array(real_curves['USA'][2]) + np.array(real_curves['USA'][3]))[0]


fitting_once = []
for i in range(3):
    vgnr_inf = gnr(srv_func[i], 'func', inf_params[i])
    vgnr_rem = gnr(srv_func[i], 'func', rem_params[i])
    calc_params['occur_inf'] = occur(vgnr_inf, vgnr_type = 'srv', length = occur_length, step = calc_params['step'])
    calc_params['occur_rem'] = occur(vgnr_rem, vgnr_type = 'srv', length = occur_length, step = calc_params['step'])
    fitting_once.append(sir_model(**calc_params))
    fitting_once[i].set_i_in_init(np.tile(np.append(np.array(real_curves['USA'][0][0]), np.diff(np.array(real_curves['USA'][0])))[:inf_data_len], [8, 1]).T)
colors = ['tab:red', 'tab:blue', 'tab:green']
for i in range(3):
    for j in range(len(real_curves['USA'][0]) - 1):
        fitting_once[i].spread_once()
    plt.plot(fitting_once[i].get_c_tot(), color = colors[i])
plt.plot(real_curves['USA'][0], color = 'black')
'''






