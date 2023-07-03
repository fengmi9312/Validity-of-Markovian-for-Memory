# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 13:35:31 2020

@author: Tingting
"""

import pandas as pd
from import_populations import import_populations
from math import log
import numpy as np

def import_curves(country, start = None, end = None):
    tmp, population = import_populations(country)
    if country == 'Britain':
        country = 'United Kingdom'
    elif country == 'USA':
        country = 'US'
    elif country == 'South Korea':
        country = 'Korea, South'
    else:
        pass
    c_data = pd.read_csv("data/curves/time_series_covid19_confirmed_global.csv")
    rec_data = pd.read_csv("data/curves/time_series_covid19_recovered_global.csv")
    death_data = pd.read_csv("data/curves/time_series_covid19_deaths_global.csv")
    c_actual = c_data.loc[(c_data['Country/Region'] == country)&\
               (c_data['Province/State'].isnull()),\
                '1/22/20':].values[0,:]/population
    rec_actual = rec_data.loc[(rec_data['Country/Region'] == country)&\
               (rec_data['Province/State'].isnull()),\
                '1/22/20':].values[0,:]/population
    death_actual = death_data.loc[(death_data['Country/Region'] == country)&\
               (death_data['Province/State'].isnull()),\
                '1/22/20':].values[0,:]/population
    inf_actual = c_actual - rec_actual - death_actual
    if start == None:
        start = np.argwhere(c_actual > 0)[0][0]
    return c_actual[start:end], inf_actual[start:end], \
            rec_actual[start:end], death_actual[start:end]

'''
import matplotlib.pyplot as plt
country = 'South Korea'
plt.figure()
c, i , rec, d = import_curves(country)
plt.plot(np.arange(len(c)), c, color = 'blue')
plt.plot(np.arange(len(i)), i, color = 'red')
plt.plot(np.arange(len(rec)), rec, color = 'green')
plt.plot(np.arange(len(d)), d, color = 'grey')
plt.title(country)
'''

def judge_init(s, r):
    length = len(s)
    x = []    
    for i in np.arange(length - 100):
        tmp = []
        for j in np.arange(i + 1, length):
            if s[j] == 0 or s[i] == 0 or r[i] >= r[j]:
                x_tmp = np.inf
            else:
                x_tmp = log(s[j] / s[i]) / (r[i] - r[j])
            tmp.append(x_tmp)
        if np.inf in tmp:
            var = np.inf
        else:
            var = np.var(tmp)
        x.append(var)
    return x.index(min(x))
    

'''
c_g, i_g, r_g, d_g = import_curves('Germany')
c_f, i_f, r_f, d_f = import_curves('Britain')
c_b, i_b, r_b, d_b = import_curves('Italy')
c_u, i_u, r_u, d_u = import_curves('USA')
plt.subplot(221)
plt.plot(np.arange(len(c_g))*0.01, c_g)
plt.plot(np.arange(len(i_g))*0.01, i_g)
plt.plot(np.arange(len(r_g))*0.01, r_g)
plt.plot(np.arange(len(d_g))*0.01, d_g)
plt.subplot(222)
plt.plot(np.arange(len(c_f))*0.01, c_f)
plt.plot(np.arange(len(i_f))*0.01, i_f)
plt.plot(np.arange(len(r_f))*0.01, r_f)
plt.plot(np.arange(len(d_f))*0.01, d_f)
plt.subplot(223)
plt.plot(np.arange(len(c_b))*0.01, c_b)
plt.plot(np.arange(len(i_b))*0.01, i_b)
plt.plot(np.arange(len(r_b))*0.01, r_b)
plt.plot(np.arange(len(d_b))*0.01, d_b)
plt.subplot(224)
plt.plot(np.arange(len(c_u))*0.01, c_u)
plt.plot(np.arange(len(i_u))*0.01, i_u)
plt.plot(np.arange(len(r_u))*0.01, r_u)
plt.plot(np.arange(len(d_u))*0.01, d_u)
'''