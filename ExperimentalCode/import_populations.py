# -*- coding: utf-8 -*-
"""
Created on Sun Dec 27 15:37:39 2020

@author: Tingting
"""

import numpy as np
import pandas as pd

def import_populations(country):
    data = pd.read_excel('data/populations/WPP2019_POP_F07_1_POPULATION_BY_AGE_BOTH_SEXES.xlsx', sheet_name = 'ESTIMATES')
    if country == 'USA':
        country = 'United States of America'
    elif country == 'South Korea':
        country = 'Republic of Korea'
    data_line = data.loc[(data['Unnamed: 2'] == country)\
                    & (data['Unnamed: 7'] == 2020)]
    data_array = data_line.values[0][8:]*1000
    return (np.append(data_array[:15], sum(data_array[15:]))/sum(data_array)).astype('float64'), sum(data_array)
