# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 10:17:39 2023

@author: admin
"""

import pandas as pd

def import_from_file(_file_dir):
    _data = pd.read_excel(_file_dir, index_col = 0, sheet_name = None)
    for _key in _data.keys():
        for _key1 in _data[_key].keys():
            _data[_key][_key1] = _data[_key][_key1].to_numpy()
    return _data

figure_theories_data_imported = import_from_file('./generated_figure_data/theories_data(Figure_2).xlsx')
si_figure_percentile_data_imported = import_from_file('./generated_figure_data/percentile_data(Figure_S3).xlsx')
si_figure_markovian_time_data_imported = import_from_file('./generated_figure_data/markovian_time_data(Figure_S1).xlsx')
si_figure_error_percentile_data_imported = import_from_file('./generated_figure_data/error_percentile_data(Figure_S2).xlsx')
figure_analysis_data_imported = import_from_file('./generated_figure_data/analysis_data(Figure_4).xlsx')
figure_forecasting_data_imported = import_from_file('./generated_figure_data/forecasting_data(Figure_5).xlsx')
figure_prevention_data_imported = import_from_file('./generated_figure_data/prevention_data(Figure_6).xlsx')
fitting_period_data_imported = import_from_file('./generated_figure_data/fitting_period_data(Figure_4_and_5).xlsx')

import figure_func
import matplotlib.pyplot as plt

figure_exported = True
figure_closed = True

fig = figure_func.figure_theories(figure_theories_data_imported)[0]
if figure_exported:
    fig.savefig('./generated_figures/theories.pdf', format = 'pdf')
if figure_closed:
    plt.close(fig)

fig = figure_func.figure_analysis(figure_analysis_data_imported)[0]
if figure_exported:
    fig.savefig('./generated_figures/analysis.pdf', format = 'pdf')
if figure_closed:
    plt.close(fig)

fig = figure_func.figure_forecasting(figure_forecasting_data_imported, fitting_period_data_imported)[0]
if figure_exported:
    fig.savefig('./generated_figures/forecasting.pdf', format = 'pdf')
if figure_closed:
    plt.close(fig)

fig = figure_func.figure_prevention(figure_prevention_data_imported, fitting_period_data_imported)[0]
if figure_exported:
    fig.savefig('./generated_figures/prevention.pdf', format = 'pdf')
if figure_closed:
    plt.close(fig)

fig = figure_func.si_figure_percentile(si_figure_percentile_data_imported)[0]
if figure_exported:
    fig.savefig('./generated_figures/si_percentile.pdf', format = 'pdf')
if figure_closed:
    plt.close(fig)

fig = figure_func.si_figure_markovian_time(si_figure_markovian_time_data_imported)[0]
if figure_exported:
    fig.savefig('./generated_figures/si_markovian_time.pdf', format = 'pdf')
if figure_closed:
    plt.close(fig)

fig = figure_func.si_figure_error_percentile(si_figure_error_percentile_data_imported)[0]
if figure_exported:
    fig.savefig('./generated_figures/si_error_percentile.pdf', format = 'pdf')
if figure_closed:
    plt.close(fig)
