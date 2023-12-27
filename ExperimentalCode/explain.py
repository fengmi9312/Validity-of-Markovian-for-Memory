# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 20:46:44 2022

@author: 20481756
"""

import numpy as np
import data_import as di
import metapopulation_simulation as ms
from scipy.special import gamma
import pandas as pd
from model_n import get_fitting_data, calc_k, get_calc_data

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
day_div = 100
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

alpha_inf, beta_inf, alpha_rem, beta_rem = 2.5, 5, 2.5, 5
mean_inf = get_mean_from_weibull(alpha_inf, beta_inf)
mean_rem = get_mean_from_weibull(alpha_rem, beta_rem)
lambda_eff = lambda_eff_weibull(alpha_inf, beta_inf, alpha_rem, beta_rem)

calc_params = {'i0': i0_from_simu, 'populations': populations_from_simu, 'contacts': params_adj['contacts'], 
               'ifrs': params_adj['ifrs'], 'ylls': params_adj['ylls'], 
               'contacts_prop': None, 'k': 1, 'occur_inf': None, 'occur_rem': None,
               'delay': delay, 'eta': eta, 'step': step, 'time': time0}

k = calc_k(alpha_inf, beta_inf, alpha_rem, beta_rem, occur_length, steady_prdt, calc_params)

simu_structure_params['k'] = k
flevel = 0.2
params_fitting = {'method': 'L-BFGS-B', 'tol': 1e-12}
from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
writer = pd.ExcelWriter("../ExperimentalDataTmp/explain_data/data_"+ str(rank) + ".xlsx")
writer_add = pd.ExcelWriter("../ExperimentalDataTmp/explain_data/data_curves_"+ str(rank) + ".xlsx")
param_list = {}
res_list = {}
rem_res_list = {}
for idx, alpha_inf_n in enumerate((2.5, 1.5, 3.5)):
    fitting_res = get_fitting_data(rank + 1, (alpha_inf_n, get_beta_inf(alpha_inf_n, alpha_rem, beta_rem, lambda_eff), alpha_rem, beta_rem), 
                                   k, amount_seeds, 
                                   simu_structure_params, calc_params, params_fitting, fitting_level = flevel)
    pd.DataFrame(fitting_res['simu_curves']).to_excel(writer_add, sheet_name = 'simu_curve_'+str(idx)) 
    if idx == 0:
        for key in fitting_res['params'].keys():
            param_list[key] = [fitting_res['params'][key]]
    else:
        for key in fitting_res['params'].keys():
            param_list[key].append(fitting_res['params'][key])
    calc_res = get_calc_data(fitting_res['params'], fitting_res['i_in_data'], k, fitting_res['simu_curves'], calc_params)
    pd.DataFrame(calc_res['calc_curves']).to_excel(writer_add, sheet_name = 'calc_curve_'+str(idx)) 
    pd.DataFrame(calc_res['rem_curves']).to_excel(writer_add, sheet_name = 'rem_curve'+str(idx)) 
    if idx == 0:
        for key in calc_res['res'].keys():
            res_list[key] = [calc_res['res'][key]]
        for key in calc_res['rem_res'].keys():
            rem_res_list[key] = [calc_res['rem_res'][key]]
    else:
        for key in calc_res['res'].keys():
            res_list[key].append(calc_res['res'][key])
        for key in calc_res['rem_res'].keys():
            rem_res_list[key].append(calc_res['rem_res'][key])
pd.DataFrame(param_list).to_excel(writer, sheet_name = 'params')    
pd.DataFrame(res_list).to_excel(writer, sheet_name = 'res')  
pd.DataFrame(rem_res_list).to_excel(writer, sheet_name = 'rem_res')
writer.close()
writer_add.close()


