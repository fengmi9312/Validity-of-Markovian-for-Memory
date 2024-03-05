# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 15:22:30 2023

@author: admin
"""
import sys
sys.path.append('../Dependencies/CodeDependencies')
import pandas as pd
import numpy as np
import func



trans_list = [['weibull', 5, 7, 6], 
              ['weibull', 5, 5, 6],
              ['weibull', 7, 7, 6],
              ['lognormal', 5, 7, 6],
              ['gamma', 5, 7, 6], 
              ['weibull', 7, 5, 6],
              ['weibull', 6, 6, 6],
              ['weibull', 5, 7, 4],
              ['weibull', 5, 7, 8],]
param_range = {'weibull': np.arange(-6, 25, 1), 'gamma': np.arange(-6, 25, 1), 'lognormal': np.arange(6, 37, 1)}
srv_funcs = {'weibull': func.srv_weibull_scale, 'gamma': func.srv_gamma_scale, 'lognormal': func.srv_lognormal}

def get_func_ir_name(_func_type, _mean_inf, _mean_rem, _steady_level):
    if _steady_level == 6:
        return _func_type + '_i' + str(_mean_inf) + 'r' + str(_mean_rem)
    else:
        return _func_type + '_i' + str(_mean_inf) + 'r' + str(_mean_rem) + '_level' + str(_steady_level)






    
def import_transient_data():
    transient_data = {}
    for func_type, mean_inf, mean_rem, steady_level in trans_list[:5] + trans_list[7:9]:
        func_ir_name = get_func_ir_name(func_type, mean_inf, mean_rem, steady_level)
        transient_data[func_ir_name] = \
        pd.read_excel("../ExperimentalData/transient_equivalence_data/"+ 
                      func_type + "_" + str(mean_inf) + "_" + str(mean_rem) + "_" + str(steady_level) + ".xlsx", index_col = 0, sheet_name = None)
    return transient_data

def import_r0_g_data():
    r0_g_data = {}
    for func_type, mean_inf, mean_rem, steady_level in trans_list:
        func_ir_name = get_func_ir_name(func_type, mean_inf, mean_rem, steady_level)
        if steady_level == 6:
            file_name = func_type + '_i' + str(mean_inf) + 'r' + str(mean_rem)
        else:
            file_name = func_type + '_i' + str(mean_inf) + 'r' + str(mean_rem) + '_level_' + str(steady_level)
        r0_g_data[func_ir_name] = pd.read_excel("../ExperimentalData/r0_g_data/r0_g_calc_"+ file_name + ".xlsx", index_col = 0, sheet_name = None)
    return r0_g_data
        
def import_simu_data():
    _data = pd.read_excel('../ExperimentalData/theory_compare_ijk_data/simu_data.xlsx', index_col = 0, sheet_name = None)
    print('simu data imported')
    return _data
    
def import_calc_data():
    _data = pd.read_excel('../ExperimentalData/theory_compare_ijk_data/calc_data.xlsx', index_col = 0, sheet_name = None)
    print('calc data imported')
    return _data
    
def import_simu_data_m():
    _data = pd.read_excel('../ExperimentalData/theory_compare_ijk_data/simu_data_m.xlsx', index_col = 0, sheet_name = None)
    print('m-simu data imported')
    return _data
    
def import_simu_data_nm():
    _data = pd.read_excel('../ExperimentalData/theory_compare_ijk_data/simu_data_nm.xlsx', index_col = 0, sheet_name = None)
    print('nm-simu data imported')
    return _data
    
def import_steady_simu_data():
    _data = pd.read_excel('../ExperimentalData/theory_compare_ijk_data/steady_simu_data.xlsx', index_col = 0, sheet_name = None)
    print('steady-simu data imported')
    return _data
    
def import_steady_calc_data():
    _data = pd.read_excel('../ExperimentalData/theory_compare_ijk_data/steady_calc_data.xlsx', index_col = 0, sheet_name = None)
    print('steady-calc data imported')
    return _data

def import_steady_simu_mar_data():
    _data = pd.read_excel('../ExperimentalData/theory_compare_ijk_data/steady_simu_mar_data.xlsx', index_col = 0, sheet_name = None)
    print('steady-simu-mar data imported')
    return _data
    
def import_steady_calc_mar_data():
    _data = pd.read_excel('../ExperimentalData/theory_compare_ijk_data/steady_calc_mar_data.xlsx', index_col = 0, sheet_name = None)
    print('steady-calc-mar data imported')
    return _data
    
def import_explain_data():
    explain_data = {}
    for i in range(100):
        explain_data['data_curves_' + str(i)] = pd.read_excel('../ExperimentalData/explain_data/data_curves_'+ str(i) + '.xlsx', index_col = 0, sheet_name = None)
    print('explain data imported')
    return explain_data

def import_forecasting_curve_data():
    forecasting_curve_data = {}
    for file_mark in ['-6_24', '9_9', '24_-6']:
        forecasting_curve_data[file_mark] = pd.read_excel('../ExperimentalData/forecasting_data_i5r7/data_curves_'+ file_mark + '.xlsx', index_col = 0, sheet_name = None)
    print('forecasting-curve data imported')
    return forecasting_curve_data
 
def import_forecasting_data():   
    folder_name = {'weibull_i5r7': 'forecasting_data_i5r7', 
                   'weibull_i5r5': 'forecasting_data_i5r5', 
                   'weibull_i7r7': 'forecasting_data_i7r7', 
                   'weibull_i7r5': 'forecasting_data_i7r5', 
                   'weibull_i6r6': 'forecasting_data_i6r6', 
                   'lognormal_i5r7': 'forecasting_data_lognormal', 
                   'gamma_i5r7': 'forecasting_data_gamma'}
    forecasting_data = {}
    for func_type, mean_inf, mean_rem, steady_level in trans_list:
        forecasting_data[func_type + '_i' + str(mean_inf) + 'r' +  str(mean_rem)] = {}
        for inf_pow in param_range[func_type]:
            for rem_pow in param_range[func_type]:
                forecasting_data[func_type + '_i' + str(mean_inf) + 'r' +  str(mean_rem)]['data_'+str(inf_pow)+'_'+str(rem_pow)] = \
                pd.read_excel('../ExperimentalData/' + folder_name[func_type + '_i' + str(mean_inf) + 'r' +  str(mean_rem)] + 
                              '/data_' + str(inf_pow) + '_' + str(rem_pow) + '.xlsx', index_col = 0, sheet_name = None)
    print('forecasting data imported')
    return forecasting_data

def import_real_forecasting_data():
    real_forecasting_data = {}
    for disease in ['covid_19', 'sars', 'h1n1', 'smallpox']:
        real_forecasting_data[disease] = pd.read_excel('../ExperimentalData/forecasting_data_real/data_' + disease + '.xlsx', index_col = 0, sheet_name = None)
    print('real-forecasting data imported')
    return real_forecasting_data

def import_vac_curve_data():
    vac_curve_data = {}
    for file_mark in ['9_9']:
        vac_curve_data[file_mark] = pd.read_excel('../ExperimentalData/vaccination_data_i5r7/data_vac_curves_'+ file_mark + '.xlsx', index_col = 0, sheet_name = None)
    print('vac-curve data imported')
    return vac_curve_data

def import_vac_data():
    folder_name = {'weibull_i5r7': 'vaccination_data_i5r7', 'weibull_i5r5': 'vaccination_data_i5r5', 
                   'weibull_i7r7': 'vaccination_data_i7r7', 
                   'lognormal_i5r7': 'vaccination_data_lognormal', 'gamma_i5r7': 'vaccination_data_gamma'}
    vac_data = {}
    for func_type, mean_inf, mean_rem, steady_level in trans_list[:5]:
        vac_data[func_type + '_i' + str(mean_inf) + 'r' +  str(mean_rem)] = {}
        for inf_pow in param_range[func_type] :
            for rem_pow in param_range[func_type]:
                vac_data[func_type + '_i' + str(mean_inf) + 'r' +  str(mean_rem)]['data_'+str(inf_pow)+'_'+str(rem_pow)] = \
                pd.read_excel('../ExperimentalData/'+folder_name[func_type + '_i' + str(mean_inf) + 'r' +  str(mean_rem)] +
                               '/data_vac_' + str(inf_pow) + '_' + str(rem_pow) + '.xlsx', 
                              index_col = 0, sheet_name = None)    
    print('vac data imported')
    return vac_data
   
def import_real_vac_data():         
    real_vac_data = {}
    for disease in ['covid_19', 'sars', 'h1n1', 'smallpox']:
        real_vac_data[disease] = pd.read_excel('../ExperimentalData/vaccination_data_real/data_' + disease + '.xlsx', index_col = 0, sheet_name = None)
    print('real-vac data imported')
    return real_vac_data

'''
def import_vac_curve_data_revision():
    vac_curve_data = {}
    for file_mark in ['9_9']:
        vac_curve_data[file_mark] = pd.read_excel('../ExperimentalData/validity_vaccine_x/simu_vac_fitting_i5r7/fitting_vac_data/data_vac_curves_'+ file_mark + '.xlsx', index_col = 0, sheet_name = None)
    print('vac-curve-revision data imported')
    return vac_curve_data

def import_vac_data_revision():
    folder_name = {'weibull_i5r7': 'simu_vac_fitting_i5r7/fitting_vac_data', 'weibull_i5r5': 'simu_vac_fitting_i5r5/fitting_vac_data', 
                   'weibull_i7r7': 'simu_vac_fitting_i7r7/fitting_vac_data', 
                   'lognormal_i5r7': 'simu_vac_fitting_lognormal/fitting_vac_data', 'gamma_i5r7': 'simu_vac_fitting_gamma/fitting_vac_data'}
    vac_data = {}
    for func_type, mean_inf, mean_rem, steady_level in trans_list[:5]:
        vac_data[func_type + '_i' + str(mean_inf) + 'r' +  str(mean_rem)] = {}
        for inf_pow in param_range[func_type] :
            for rem_pow in param_range[func_type]:
                vac_data[func_type + '_i' + str(mean_inf) + 'r' +  str(mean_rem)]['data_'+str(inf_pow)+'_'+str(rem_pow)] = \
                pd.read_excel('../ExperimentalData/validity_vaccine_x/'+folder_name[func_type + '_i' + str(mean_inf) + 'r' +  str(mean_rem)] +
                               '/data_vac_' + str(inf_pow) + '_' + str(rem_pow) + '.xlsx', 
                              index_col = 0, sheet_name = None)    
    print('vac data-revision imported')
    return vac_data
   
def import_real_vac_data_revision():         
    real_vac_data = {}
    for disease in ['covid_19', 'sars', 'h1n1', 'smallpox']:
        real_vac_data[disease] = pd.read_excel('../ExperimentalData/validity_vaccine_x/simu_vac_fitting_real/real_vac_fitting_data/data_' + disease + '.xlsx', index_col = 0, sheet_name = None)
    print('real-vac data-revision imported')
    return real_vac_data
'''
























    
