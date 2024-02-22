# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 20:46:44 2022

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
from model_n import calc_k, sir_model
import metapopulation_simulation as ms
from scipy.special import gamma
import func
from occur import gnr, occur
import pandas as pd

import scipy.special as sc
from scipy.optimize import fsolve

def lambda_eff_weibull(_alpha_inf, _beta_inf, _alpha_rem, _beta_rem):
    return ((_beta_rem/_beta_inf)**_alpha_inf)*(_alpha_inf/_alpha_rem)*sc.gamma(_alpha_inf/_alpha_rem)

def get_beta_inf(_alpha_inf, _alpha_rem, _beta_rem, _lambda_eff):
    def get_lambda_sub(_beta_inf):
        return lambda_eff_weibull(_alpha_inf,_beta_inf, _alpha_rem, _beta_rem) - _lambda_eff
    return fsolve(get_lambda_sub, 1.0, xtol = 1e-6)[0]

def get_mean_from_weibull(alpha, beta):
    return beta * gamma(1 + 1.0 / alpha)


np.random.seed(0)

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
steady_prdt = 0.6
amount_seeds = (np.ones(len(group_div)) * 5).astype(int)
simu_times = 100
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
                       'contacts': contacts, 'ifrs': params_adj['ifrs'], 'ylls': params_adj['ylls'], 'k': 1, 'step': step, 'time': 0}
simu_vaccination_params = {'delay': delay, 'eta': eta}
param_simu_obtain = ms.simu(**simu_structure_params)
param_simu_obtain.set_spreading_func('weibull', 'weibull')
param_simu_obtain.set_total_spreading_params([1, 1], [1, 1])
param_simu_obtain.set_amount_seeds(amount_seeds)
populations_from_simu = np.array(param_simu_obtain.get_populations())
population_amounts_from_simu = np.array(param_simu_obtain.get_population_amounts())
i0_from_simu = np.array(param_simu_obtain.get_i_amount_arr()) / np.array(param_simu_obtain.get_population_amounts())
del param_simu_obtain

alpha_inf, beta_inf, alpha_rem, beta_rem = 3, 4, 2.5, 3.6
mean_inf = get_mean_from_weibull(alpha_inf, beta_inf)
mean_rem = get_mean_from_weibull(alpha_rem, beta_rem)
lambda_eff = lambda_eff_weibull(alpha_inf, beta_inf, alpha_rem, beta_rem)

calc_params = {'i0': i0_from_simu, 'populations': populations_from_simu, 'contacts': params_adj['contacts'], 
               'ifrs': params_adj['ifrs'], 'ylls': params_adj['ylls'], 
               'contacts_prop': None, 'k': 1, 'occur_inf': None, 'occur_rem': None,
               'delay': delay, 'eta': eta, 'step': step, 'time': time0}

vgnr_inf = gnr(func.srv_weibull_scale, 'func', [alpha_inf, beta_inf])
vgnr_rem = gnr(func.srv_weibull_scale, 'func', [alpha_rem, beta_rem])
calc_params['occur_inf'] = occur(vgnr_inf, vgnr_type = 'srv', length = occur_length, step = calc_params['step'])
calc_params['occur_rem'] = occur(vgnr_rem, vgnr_type = 'srv', length = occur_length, step = calc_params['step'])
k = calc_k(alpha_inf, beta_inf, alpha_rem, beta_rem, occur_length, steady_prdt, calc_params)

calc_params['k'] = k
simu_structure_params['k'] = k

if rank == 0:
    writer_simu = pd.ExcelWriter("../ExperimentalData/theory_compare_ijk_data/simu_data.xlsx")
    writer_calc = pd.ExcelWriter("../ExperimentalData/theory_compare_ijk_data/calc_data.xlsx")
    simu_curve_data = []
    max_spread_len = 0
    for i in range(simu_times):
        simu_curve_data.append({'s': [], 'i': [], 'r': [], 'c': [], 'd': [], 'y': []})
        simu_once = ms.simu(**simu_structure_params)
        simu_once.set_generator_seed(i + 1)
        simu_once.set_spreading_func('weibull', 'weibull')
        simu_once.set_total_spreading_params([alpha_inf, beta_inf],
                                             [alpha_rem, beta_rem])
        simu_once.set_amount_seeds(amount_seeds)
        spread_len = 0
        for data_key in simu_curve_data[-1].keys():
            simu_curve_data[-1][data_key].append(simu_once.get_x_amount(data_key) / simu_once.get_node_amount())
        while simu_once.get_i_amount() > 0:
            simu_once.spread_once()
            for data_key in simu_curve_data[-1].keys():
                simu_curve_data[-1][data_key].append(simu_once.get_x_amount(data_key) / simu_once.get_node_amount())
            spread_len += 1
        if spread_len > max_spread_len:
            max_spread_len = spread_len
    for i in range(simu_times): 
        pd.DataFrame(simu_curve_data[i]).to_excel(writer_simu, sheet_name = 'data_' + str(i))  
    writer_simu.close()
    
    calc_once = sir_model(**calc_params)
    for i in range(max_spread_len):
        calc_once.spread_once()
    calc_curve_data = {'s': None, 'i': None, 'r': None, 'c': None, 'd': None, 'y': None}
    for data_key in calc_curve_data.keys():
        calc_curve_data[data_key] = calc_once.get_x_tot(data_key)
        
    pd.DataFrame(calc_curve_data).to_excel(writer_calc, sheet_name = 'data_calc')  
    writer_calc.close()


elif rank == 1:
    simu_curve_data_n = []
    writer_simu_m = pd.ExcelWriter("../ExperimentalData/theory_compare_ijk_data/simu_data_m.xlsx")
    writer_simu_nm = pd.ExcelWriter("../ExperimentalData/theory_compare_ijk_data/simu_data_nm.xlsx")
    for alpha_inf_n, alpha_rem_n, beta_rem_n in ((1, 1, 7.2), (1.5, 2, 2.4)):
        beta_inf_n = get_beta_inf(alpha_inf_n, alpha_rem_n, beta_rem_n, lambda_eff)
        simu_curve_data_n.append([])
        for i in range(simu_times):
            simu_curve_data_n[-1].append({'s': [], 'i': [], 'r': [], 'c': [], 'd': [], 'y': []})
            simu_once = ms.simu(**simu_structure_params)
            simu_once.set_generator_seed(i + 1)
            simu_once.set_spreading_func('weibull', 'weibull')
            simu_once.set_total_spreading_params([alpha_inf_n, beta_inf_n],
                                                 [alpha_rem_n, beta_rem_n])
            simu_once.set_amount_seeds(amount_seeds)
            for data_key in simu_curve_data_n[-1][-1].keys():
                simu_curve_data_n[-1][-1][data_key].append(simu_once.get_x_amount(data_key) / simu_once.get_node_amount())
            while simu_once.get_i_amount() > 0:
                simu_once.spread_once()
                for data_key in simu_curve_data_n[-1][-1].keys():
                    simu_curve_data_n[-1][-1][data_key].append(simu_once.get_x_amount(data_key) / simu_once.get_node_amount())
            del simu_once
    for i in range(simu_times): 
        pd.DataFrame(simu_curve_data_n[0][i]).to_excel(writer_simu_m, sheet_name = 'data_' + str(i))  
    for i in range(simu_times): 
        pd.DataFrame(simu_curve_data_n[1][i]).to_excel(writer_simu_nm, sheet_name = 'data_' + str(i))  
    writer_simu_m.close()
    writer_simu_nm.close()

elif rank == 2:
    steady_point_len = 40
    writer_steady_simu = pd.ExcelWriter("../ExperimentalData/theory_compare_ijk_data/steady_simu_data.xlsx")
    writer_steady_calc = pd.ExcelWriter("../ExperimentalData/theory_compare_ijk_data/steady_calc_data.xlsx")
    from copy import deepcopy
    simu_steady_data = []
    calc_steady_data = {'c': [], 'd': [], 'y': []}
    for lambda_eff_n in np.linspace(0.01, 2, steady_point_len):
        simu_steady_data.append({'c': [], 'd': [], 'y': []})
        beta_inf_n = get_beta_inf(alpha_inf, alpha_rem, beta_rem, lambda_eff_n)
        for i in range(simu_times):
            simu_once = ms.simu(**simu_structure_params)
            simu_once.set_generator_seed(i + 1)
            simu_once.set_spreading_func('weibull', 'weibull')
            simu_once.set_total_spreading_params([alpha_inf, beta_inf_n],
                                                 [alpha_rem, beta_rem])
            simu_once.set_amount_seeds(amount_seeds)
            while simu_once.get_i_amount() > 0:
                simu_once.spread_once()
            for data_key in simu_steady_data[-1].keys():
                simu_steady_data[-1][data_key].append(simu_once.get_x_amount(data_key) / simu_once.get_node_amount())
            del simu_once
        _vgnr_inf = gnr(func.srv_weibull_scale, 'func', (alpha_inf, beta_inf_n))
        _vgnr_rem = gnr(func.srv_weibull_scale, 'func', (alpha_rem, beta_rem))
        _calc_params = deepcopy(calc_params)
        _calc_params['occur_inf'] = occur(_vgnr_inf, vgnr_type = 'srv', length = occur_length, step = _calc_params['step'])
        _calc_params['occur_rem'] = occur(_vgnr_rem, vgnr_type = 'srv', length = occur_length, step = _calc_params['step'])
        calc_once = sir_model(**_calc_params)
        for data_key in calc_steady_data.keys():
            calc_steady_data[data_key].append(calc_once.get_prdt_from_init(target = data_key))
    for i in range(steady_point_len): 
        pd.DataFrame(simu_steady_data[i]).to_excel(writer_steady_simu, sheet_name = 'data_' + str(i))  
    pd.DataFrame(calc_steady_data).to_excel(writer_steady_calc, sheet_name = 'data')
    writer_steady_simu.close()
    writer_steady_calc.close()
elif rank == 3:
    alpha_inf, beta_inf, alpha_rem, beta_rem = 1, 4, 1, 4
    steady_point_len = 40
    writer_steady_simu_mar = pd.ExcelWriter("../ExperimentalData/theory_compare_ijk_data/steady_simu_mar_data.xlsx")
    writer_steady_calc_mar = pd.ExcelWriter("../ExperimentalData/theory_compare_ijk_data/steady_calc_mar_data.xlsx")
    from copy import deepcopy
    simu_steady_data = []
    calc_steady_data = {'c': [], 'd': [], 'y': []}
    for lambda_eff_n in np.linspace(0.01, 2, steady_point_len):
        simu_steady_data.append({'c': [], 'd': [], 'y': []})
        beta_inf_n = get_beta_inf(alpha_inf, alpha_rem, beta_rem, lambda_eff_n)
        for i in range(simu_times):
            simu_once = ms.simu(**simu_structure_params)
            simu_once.set_generator_seed(i + 1)
            simu_once.set_spreading_func('weibull', 'weibull')
            simu_once.set_total_spreading_params([alpha_inf, beta_inf_n],
                                                 [alpha_rem, beta_rem])
            simu_once.set_amount_seeds(amount_seeds)
            while simu_once.get_i_amount() > 0:
                simu_once.spread_once()
            for data_key in simu_steady_data[-1].keys():
                simu_steady_data[-1][data_key].append(simu_once.get_x_amount(data_key) / simu_once.get_node_amount())
            del simu_once
        _vgnr_inf = gnr(func.srv_weibull_scale, 'func', (alpha_inf, beta_inf_n))
        _vgnr_rem = gnr(func.srv_weibull_scale, 'func', (alpha_rem, beta_rem))
        _calc_params = deepcopy(calc_params)
        _calc_params['occur_inf'] = occur(_vgnr_inf, vgnr_type = 'srv', length = occur_length, step = _calc_params['step'])
        _calc_params['occur_rem'] = occur(_vgnr_rem, vgnr_type = 'srv', length = occur_length, step = _calc_params['step'])
        calc_once = sir_model(**_calc_params)
        for data_key in calc_steady_data.keys():
            calc_steady_data[data_key].append(calc_once.get_prdt_from_init(target = data_key))
    for i in range(steady_point_len): 
        pd.DataFrame(simu_steady_data[i]).to_excel(writer_steady_simu_mar, sheet_name = 'data_' + str(i))  
    pd.DataFrame(calc_steady_data).to_excel(writer_steady_calc_mar, sheet_name = 'data')
    writer_steady_simu_mar.close()
    writer_steady_calc_mar.close()
else:
    pass







