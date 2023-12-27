# -*- coding: utf-8 -*-
"""
Created on Fri May 20 14:30:47 2022

@author: 20481756
"""
import sys
sys.path.append('../Dependencies/CodeDependencies')
sys.path.append('../Dependencies/ModuleDependencies(Windows)')
import numpy as np
import data_import as di
from model_n import get_general_fitting_data, calc_general_k, get_general_v_data
import metapopulation_simulation as ms
from scipy.special import gamma, gammainc, erf
import pandas as pd

np.random.seed(0)

# self-defined important function
def get_mean_from_weibull(alpha, beta):
    return beta * gamma(1 + 1.0 / alpha)

def get_beta_from_weibull(alpha, mean_value):
    return mean_value / gamma(1 + 1.0 / alpha)

def gamma_survivals(_alpha, _beta, _length, _step):
    _tau = np.arange(_length) * _step
    return 1 - gammainc(_alpha, _tau / _beta)

def lognormal_survivals(_alpha, _beta, _length, _step):
    _tau = np.arange(1, _length) * _step
    return np.append(1, 0.5 - 0.5 * erf((np.log(_tau) - _beta) / (_alpha * (2 ** 0.5))))

def split_int(x, a):
    tmp = [a] * (x // a)
    sum_tmp = sum(tmp)
    if sum_tmp != x:
        tmp.append(x - sum_tmp)
    return np.array(tmp)

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
steady_prdt = 0.6
amount_seeds = (np.ones(len(group_div)) * 30).astype(int)
flevel = 0.2
params_fitting = {'method': 'L-BFGS-B', 'tol': 1e-16}


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


_distribution = 'gamma'  # two options: 'gamma' and 'lognormal'
mean_inf = 5 # Here, you can set the value of the average infection time for the gamma/log-normal infection time distribution
mean_rem = 7 # Here, you can set the value of the average removal time for the gamma/log-normal removal time distribution

rank = 0
p_num = 100
t_num = 100 # To accommodate the huge calculations, we have divided the calculation parameters into t_num parts. 
            # Each calculation is assigned a specific parameter part.

'''
import pandas as pd
from mpi4py import MPI
comm = MPI.COMM_WORLD
t_num = comm.Get_size()
rank = comm.Get_rank()
'''          # These five lines of code ensure that the code runs on the server utilizing Slurm to allocate computational resources to each calculation.

pair_list = []
for i in np.arange(24, -7, -1):
    for j in np.arange(-6, 25):
        pair_list.append((i, j))
group_param_list = [[] for i in range(t_num)]
for i in range(t_num):
    for j in range(len(pair_list) // t_num + 1):
        if j * t_num + i < len(pair_list):
            group_param_list[i].append(pair_list[j * t_num + i])
for inf_pow, rem_pow in group_param_list[rank]:
    add_mark = False
    if (inf_pow, rem_pow) == (-6, -6) or (inf_pow, rem_pow) == (-6, 24) or (inf_pow, rem_pow) == (24, -6) or \
    (inf_pow, rem_pow) == (24, 24) or (inf_pow, rem_pow) == (9, 9):
        add_mark = True
    file_mark = str(inf_pow) + "_" + str(rem_pow)
    steady_data = {}
    writer = pd.ExcelWriter("../ExperimentalDataTmp/fitting_vac_data_"+_distribution+"_"+"i"+str(mean_inf)+"r"+str(mean_rem)+"/data_vac_"+ file_mark + ".xlsx")
    if add_mark:
        writer_add = pd.ExcelWriter("../ExperimentalDataTmp/fitting_vac_data_"+_distribution+"_"+"i"+str(mean_inf)+"r"+str(mean_rem)+"/data_vac_curves_"+ file_mark + ".xlsx")
    for i in range(p_num):
        alpha_inf = np.e ** (inf_pow * 0.1)
        alpha_rem = np.e ** (rem_pow * 0.1)
        beta_inf = mean_inf / alpha_inf
        beta_rem = mean_rem / alpha_rem
        if _distribution == 'gamma':
            survivals_inf = gamma_survivals(alpha_inf, beta_inf, occur_length, step) # using gamma distribution
            survivals_rem = gamma_survivals(alpha_rem, beta_rem, occur_length, step) # using gamma distribution
        elif _distribution == 'lognormal':
            survivals_inf = lognormal_survivals(alpha_inf, beta_inf, occur_length, step)  # using lognormal distribution
            survivals_rem = lognormal_survivals(alpha_rem, beta_rem, occur_length, step)  # using lognormal distribution
        else:
            pass
        k = calc_general_k(survivals_inf, survivals_rem, steady_prdt, calc_params)
        fitting_res = get_general_fitting_data(i + 1, (survivals_inf, survivals_rem, step), k, amount_seeds, 
                                       simu_structure_params, calc_params, params_fitting, fitting_level = flevel)
        calc_res = get_general_v_data(i + 1, (survivals_inf, survivals_rem, step), fitting_res['params'], 
                              len(fitting_res['simu_curves']['c']), fitting_res['i_in_data'], k, 
                              amount_seeds, simu_structure_params, calc_params, 300)
        if add_mark:
            for vac_key in calc_res['simu_data'].keys():
                pd.DataFrame(calc_res['simu_data'][vac_key]).to_excel(writer_add, sheet_name = 'simu_curve_' + vac_key + '_' + str(i)) 
            for vac_key in calc_res['calc_data'].keys():
                pd.DataFrame(calc_res['calc_data'][vac_key]).to_excel(writer_add, sheet_name = 'calc_curve_' + vac_key + '_' + str(i)) 
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
    if add_mark:
        writer_add.close()            

