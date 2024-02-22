# -*- coding: utf-8 -*-
"""
Created on Fri May 20 14:30:47 2022

@author: 20481756
"""
import sys
sys.path.append('../Dependencies/CodeDependencies')

# If you are running the code on Linux, you can utilize the following code that is compatible with Slurm and mpi4py.
'''
sys.path.append('../Dependencies/ModuleDependencies(Linux)')
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
    import pandas as pd
    
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
    writer = pd.ExcelWriter("../ExperimentalData/forecasting_data_real/data_"+ disease[rank] + ".xlsx")
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
        calc_res = get_calc_data(fitting_res, _i_in_data, simu_structure_params['k'], _simu_once_data, calc_params)
        if i == 0:
            for res_key in calc_res['res'].keys():
                res[res_key] = [calc_res['res'][res_key], ]
        else:
            for res_key in calc_res['res'].keys():
                res[res_key].append(calc_res['res'][res_key])
    pd.DataFrame(res).to_excel(writer, sheet_name = 'res')    
    writer.close()
    




