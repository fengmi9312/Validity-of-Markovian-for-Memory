# -*- coding: utf-8 -*-
"""
Created on Mon Jan 24 15:12:48 2022

@author: 20481756
"""

import pandas as pd
import numpy as np
from scipy.interpolate import interp1d

def import_contacts(country):
    sheet_names_1 = pd.ExcelFile('data/contact_matrices_152_countries/MUestimates_all_locations_1.xlsx').sheet_names
    sheet_names_2 = pd.ExcelFile('data/contact_matrices_152_countries/MUestimates_all_locations_2.xlsx').sheet_names
    idx = 0
    if country in sheet_names_1:
        idx = 1
    elif country in sheet_names_2:
        idx = 2
    else:
        pass
    if idx == 0:
        return None
    home = pd.read_excel('data/contact_matrices_152_countries/MUestimates_home_' + str(idx) + '.xlsx', 
                      sheet_name = country, header = None).values
    school = pd.read_excel('data/contact_matrices_152_countries/MUestimates_school_' + str(idx) + '.xlsx', 
                      sheet_name = country, header = None).values
    work = pd.read_excel('data/contact_matrices_152_countries/MUestimates_work_' + str(idx) + '.xlsx', 
                      sheet_name = country, header = None).values
    other_locations = pd.read_excel('data/contact_matrices_152_countries/MUestimates_other_locations_' + str(idx) + '.xlsx', 
                      sheet_name = country, header = None).values
    return  {'home': home,
             'school': school,
             'work': work ,
             'other_locations': other_locations}

               
def import_populations(country):
    data = pd.read_excel('data/populations/WPP2019_POP_F07_1_POPULATION_BY_AGE_BOTH_SEXES.xlsx', sheet_name = 'ESTIMATES')
    data_line = data.loc[(data['Unnamed: 2'] == country)\
                    & (data['Unnamed: 7'] == 2020)]
    data_array = data_line.values[0][8:]*1000
    return (np.append(data_array[:15], sum(data_array[15:]))/sum(data_array)).astype('float64'), sum(data_array)

def import_ifrs():
    original_data_x = np.linspace(4.5,94.5,10)
    original_data_y = np.array([2/858, 5/1591, 23/13304, 61/22423, 198/34793, \
                                607/42515, 1669/34181, 4544/32323, 7728/37268, 3981/18142])
    func = interp1d(original_data_x, original_data_y, kind = 'linear',
                            fill_value =  'extrapolate')
    x_list = [i * 5 + 2 for i in range(15)]
    x_list.append(87.5)
    return func(x_list)


def import_ylls(country):
    le = {'United States of America':78.5, 'Germany':81.7, 'Brazil':75.9, 'Israel':82.6}
    data = pd.read_excel('data/populations/WPP2019_POP_F07_1_POPULATION_BY_AGE_BOTH_SEXES.xlsx', sheet_name = 'ESTIMATES')
    data_line = data.loc[(data['Unnamed: 2'] == country)\
                    & (data['Unnamed: 7'] == 2020)]
    data_array = data_line.values[0][8:]*1000
    data_array = data_array / sum(data_array)
    mid_age = np.append(np.linspace(2,97,20),100)
    rl = []
    for i in range(16):
        if i < 15:
            mean_age = mid_age[i]
        else:
            mean_age = sum(mid_age[j] * data_array[j] for j in range(15, 21)) / sum(data_array[15:])
        rl.append(max(0, le[country] - mean_age))
    return np.array(rl)



def rescale_params(param_dict, scale_arr):
    #check the params_dict
    locations = ['home', 'school', 'work', 'other_locations']
    populations = param_dict['populations']
    contacts_tot = {}
    for location in locations:
        contacts_tot[location] = param_dict['contacts'][location] * populations[:, None]
    ifrs_tot = param_dict['ifrs'] * populations
    ylls_tot = param_dict['ylls'] * populations
    scale_cum_arr = np.append(0, scale_arr.cumsum())
    contacts_n = {'home':[], 'school':[], 'work':[], 'other_locations':[]}
    populations_n = []
    ifrs_n = []
    ylls_n = []
    for i in range(1, len(scale_cum_arr)):
        populations_n.append(populations[scale_cum_arr[i - 1]:scale_cum_arr[i]].sum())
        ifrs_n.append(ifrs_tot[scale_cum_arr[i - 1]:scale_cum_arr[i]].sum())
        ylls_n.append(ylls_tot[scale_cum_arr[i - 1]:scale_cum_arr[i]].sum())
        for location in locations:
            contacts_n[location].append([])
            for j in range(1, len(scale_cum_arr)):
                contacts_n[location][i - 1].append(contacts_tot[location][scale_cum_arr[i - 1]:scale_cum_arr[i], 
                                                            scale_cum_arr[j - 1]:scale_cum_arr[j]].sum())
    populations_n = np.array(populations_n)
    ifrs_n = np.array(ifrs_n / populations_n)
    ylls_n = np.array(ylls_n / populations_n)
    for location in locations:
        contacts_n[location] = np.array(contacts_n[location] / populations_n[:, None])
    return {'contacts':contacts_n, 'populations':populations_n, 'ifrs': ifrs_n, 'ylls':ylls_n}








