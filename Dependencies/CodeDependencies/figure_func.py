# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 15:37:31 2023

@author: admin
"""

import matplotlib.pyplot as plt
import cmocean
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid.inset_locator import inset_axes
import matplotlib.patches as patches
import string
import numpy as np

trans_list = [['weibull', 5, 7, 6], 
              ['weibull', 5, 5, 6],
              ['weibull', 7, 7, 6],
              ['lognormal', 5, 7, 6],
              ['gamma', 5, 7, 6],
              ['weibull', 7, 5, 6],
              ['weibull', 6, 6, 6],
              ['weibull', 5, 7, 4],
              ['weibull', 5, 7, 8],]
def get_func_ir_name(_func_type, _mean_inf, _mean_rem, _steady_level):
    if _steady_level == 6:
        return _func_type + '_i' + str(_mean_inf) + 'r' + str(_mean_rem)
    else:
        return _func_type + '_i' + str(_mean_inf) + 'r' + str(_mean_rem) + '_level' + str(_steady_level)

def set_ax_format(_ax, _xlabel, _ylabel, _labelsize, _ticksize, _title_size, _spinetickwidth, _spinewidth, 
                  _xtick_amount, _ytick_amount, _xlims, _ylims, _frame_type, 
                  _title, _idx_title, _offset,
                  _xtick_int = False, _gap_prop = 1 / 25, _show_prop = False, _corner = True, _x_ext = False, _y_ext = False, _set_lims = True):
    _ax.set_xlabel(_xlabel, fontsize = _labelsize)
    _ax.set_ylabel(_ylabel, fontsize = _labelsize)
    _ax.spines['left'].set_linewidth(_spinewidth)
    _ax.spines['right'].set_linewidth(_spinewidth)
    _ax.spines['bottom'].set_linewidth(_spinewidth)
    _ax.spines['top'].set_linewidth(_spinewidth)
    _ax.tick_params('both', width = _spinetickwidth)
    
    if _corner:
        _x_left , _y_lower = _xlims[0] - _gap_prop * (_xlims[1] - _xlims[0]), _ylims[0]  - _gap_prop * (_ylims[1] - _ylims[0])
    else:
        _x_left , _y_lower = _xlims[0], _ylims[0]
    if _x_ext:
        _x_right = _xlims[1] + _gap_prop * (_xlims[1] - _xlims[0])
    else:
        _x_right = _xlims[1]
    if _y_ext:
        _y_upper = _ylims[1] + _gap_prop * (_ylims[1] - _ylims[0])
    else:
        _y_upper = _ylims[1]
    if _set_lims:
        _ax.set_xlim(left = _x_left, right = _x_right)
        _ax.set_ylim(bottom = _y_lower, top = _y_upper)
    if _frame_type == 0:
        _ax.spines['right'].set_visible(False)
        _ax.spines['top'].set_visible(False)
        _ax.spines['bottom'].set_bounds(low = _xlims[0], high = _xlims[1])
        _ax.spines['left'].set_bounds(low = _ylims[0], high = _ylims[1])
    elif _frame_type == 1:
        _ax.spines['right'].set_visible(False)
        _ax.spines['top'].set_visible(False)
    else:
        pass
    if _show_prop:
        _ax.set_yticks(np.linspace(_ylims[0], _ylims[1], _ytick_amount), list(map("{:.0f}".format, np.linspace(_ylims[0], _ylims[1], _ytick_amount) * 100)), fontsize = _ticksize)
    else:
        _ax.set_yticks(np.linspace(_ylims[0], _ylims[1], _ytick_amount), np.linspace(_ylims[0], _ylims[1], _ytick_amount), fontsize = _ticksize)
    if _xtick_int:
        _ax.set_xticks(np.linspace(_xlims[0], _xlims[1], _xtick_amount, dtype = np.int), np.linspace(_xlims[0], _xlims[1], _xtick_amount, dtype = np.int), fontsize = _ticksize)
    else:
        _ax.set_xticks(np.linspace(_xlims[0], _xlims[1], _xtick_amount), np.linspace(_xlims[0], _xlims[1], _xtick_amount), fontsize = _ticksize)
    _ax.set_title(_title, fontsize = _title_size)
    _ax.text(-0.1, 1.03, _idx_title, fontsize = _title_size, fontweight = 'bold', horizontalalignment='right', verticalalignment='bottom', transform=_ax.transAxes)
    return _ax
    
spinetickwidth = 1.5
offset = np.array([-0.15, 0])
labelsize = 12
spinewidth = 1.5
ticksize = 12
gap_prop = 1 / 25

def figure_theories(_figure_theories_data):
    _fig = plt.figure(figsize = [12, 9])
    _axes = []
    gs = GridSpec(31, 12, figure = _fig)
    _axes.append([])
    _axes[-1].append(_fig.add_subplot(gs[0:6, 0:3]))
    _axes[-1].append(_fig.add_subplot(gs[0:6, 4:7]))
    _axes[-1].append(_fig.add_subplot(gs[8:14, 0:3]))
    _axes.append([])
    _axes[-1].append(_fig.add_subplot(gs[8:14, 4:7]))
    _axes[-1].append(_fig.add_subplot(gs[2:12, 8:12]))
    _axes.append([])
    _axes[-1].append(_fig.add_subplot(gs[17:23, 0:3]))
    _axes[-1].append(_fig.add_subplot(gs[17:23, 4:7]))
    _axes[-1].append(_fig.add_subplot(gs[25:31, 0:3]))
    _axes[-1].append(_fig.add_subplot(gs[25:31, 4:7]))
    _axes[-1].append(_fig.add_subplot(gs[19:29, 8:12]))
    
    legend_title = ['Susceptible:', 'Infected:', 'Removed:']
    legend_y = [0.6,0.6,0.1]
    colors = ['tab:orange', 'tab:red', 'tab:purple']
    colors_calc = ['tab:brown', 'tab:blue', 'tab:green']
    for i, target in enumerate(['s', 'i', 'r']):
        ax = _axes[0][i]
        plt.sca(ax)
        key_name = target + '_curves_simu(subplot_' + string.ascii_lowercase[i] + ')'
        for j in range(100):
            if j == 99:
                plt.plot(_figure_theories_data[key_name]['time'], _figure_theories_data[key_name]['curve_simu_' + str(j)], 
                         linewidth = 1, alpha = 0.1, color = colors[i], label = 'Simulation')
            else:
                plt.plot(_figure_theories_data[key_name]['time'], _figure_theories_data[key_name]['curve_simu_' + str(j)], 
                         linewidth = 1, alpha = 0.1, color = colors[i])
        key_name = target + '_curve_calc(subplot_' + string.ascii_lowercase[i] + ')'
        plt.plot(_figure_theories_data[key_name]['time'], _figure_theories_data[key_name]['curve_calc'],
                 linewidth = 3, color = colors_calc[i], linestyle = '-', label = 'Calculation')
        plt.text(0.48, legend_y[i] + 0.4, legend_title[i], fontsize = 10, fontweight = 'bold', horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
        leg = plt.legend(fontsize = 10, bbox_to_anchor=(0.45, legend_y[i]), loc = 'lower left', frameon = False)
        leg.get_lines()[0].set_alpha(1)
        
    ax = _axes[1][0]
    plt.sca(ax)
    key_name = 'r_curves_simu(subplot_d)'
    for i in range(100):
        if i == 99:
            plt.plot(_figure_theories_data[key_name]['time'], _figure_theories_data[key_name]['nm_' + str(i)], 
                     linewidth = 1, alpha = 0.1, color = 'tab:orange', label = 'Memory-dependent')
            plt.plot(_figure_theories_data[key_name]['time'],  _figure_theories_data[key_name]['m_' + str(i)], 
                     linewidth = 1, alpha = 0.1, color = 'tab:green', label = 'Memoryless')
        else:
            plt.plot(_figure_theories_data[key_name]['time'], _figure_theories_data[key_name]['nm_' + str(i)],
                     linewidth = 1, alpha = 0.1, color = 'tab:orange')
            plt.plot(_figure_theories_data[key_name]['time'], _figure_theories_data[key_name]['m_' + str(i)], 
                     linewidth = 1, alpha = 0.1, color = 'tab:green')
    
    leg = plt.legend(fontsize = 8, loc = 'lower right')    
    for lh in leg.legendHandles: 
        lh.set_alpha(1)

    ax = _axes[1][1]
    eig_max = 2.3131550279038704
    plt.sca(ax)
    key_name = 'r_steady_state(subplot_e)'
    
    plt.plot(_figure_theories_data[key_name]['r0'], _figure_theories_data[key_name]['steady_state_calc'], 
             linewidth = 3, color = 'tab:orange', label = 'Calculation')    
    plt.plot(_figure_theories_data[key_name]['r0'], _figure_theories_data[key_name]['steady_state_simu_non_m'], 
                     marker = '+', markerfacecolor='none',markeredgecolor = 'tab:red', linestyle = 'none', label = 'Memory-dependent', markersize = 8, markeredgewidth = 1.5)
    plt.plot(_figure_theories_data[key_name]['r0'], _figure_theories_data[key_name]['steady_state_simu_m'], 
                     marker = 'x', markerfacecolor='none',markeredgecolor = 'tab:blue', linestyle = 'none', label = 'Memoryless', markersize = 8, markeredgewidth = 1.5)
    plt.axvline(1, color = 'tab:gray', linewidth = 1.5, linestyle = '--')
    leg = plt.legend(fontsize = 10, loc = 'lower right')
    leg.get_lines()[0].set_alpha(1)
    ax.text(0.21  * eig_max, 0.4, 'Absorbing\nPhase', fontsize=10, horizontalalignment= 'center')
    ax.text(1.2 * eig_max, 0.45, 'Active Phase', fontsize=10, horizontalalignment= 'center')
    ax.axvspan(0, 1, ymin = -1/25, alpha = 0.2, color='tab:green', linewidth = 0)
    ax.axvspan(1, 5 * (24 / 25), ymin = -1/25, alpha = 0.2, color='tab:red', linewidth = 0)
    
    markevery_list = [150, 75]
    for i in range(2):
        ax = _axes[2][i*2]
        plt.sca(ax)
        key_name = 'curves_calc(subplot_' + string.ascii_lowercase[5 + i * 2] + ')'
        plt.plot(_figure_theories_data[key_name]['time'], _figure_theories_data[key_name]['s_curve'], 
                 color = 'tab:blue', linewidth = 3, linestyle = '-', label = 'Susceptible')
        plt.plot(_figure_theories_data[key_name]['time'], _figure_theories_data[key_name]['r_curve'], 
                 color = 'tab:green', linewidth = 3, linestyle = '-', label = 'Removed')
        plt.plot(_figure_theories_data[key_name]['time'], _figure_theories_data[key_name]['inferred_s_curve'], 
                 'x' ,color = 'black', markevery = markevery_list[i], markeredgewidth = 2, label = 'Inferred Susceptible')
        if i == 0:
            plt.legend(fontsize = 8, loc = 'upper right')
        ax = _axes[2][i*2 + 1]
        plt.sca(ax)
        key_name = 'curves_calc(subplot_' + string.ascii_lowercase[6 + i * 2] + ')'
        plt.plot(_figure_theories_data[key_name]['time'], _figure_theories_data[key_name]['i_curve_nm'], 
                 color = 'tab:red', linewidth = 3, linestyle = '-', label = 'Infected')
        plt.plot(_figure_theories_data[key_name]['time'], _figure_theories_data[key_name]['r_curve_nm'], 
                 color = 'tab:green', linewidth = 3, linestyle = '-', label = 'Removed')
        plt.plot(_figure_theories_data[key_name]['time'], _figure_theories_data[key_name]['i_curve_m'], 
                 color = 'tab:orange', linewidth = 3, linestyle = ':', label = 'Infected')
        plt.plot(_figure_theories_data[key_name]['time'], _figure_theories_data[key_name]['r_curve_m'], 
                 color = 'tab:purple', linewidth = 3, linestyle = ':', label = 'Removed ')
        if i == 0:
            handles, handle_labels = ax.get_legend_handles_labels()
            leg = ax.legend(title_fontsize = 10, handles = handles[:2], labels = handle_labels[:2], fontsize = 8, 
                            bbox_to_anchor=(0.45, 0.5, 1, 0.1), frameon = False, loc = 'lower left')
            plt.text(0.48, 0.8, 'non-Markovian:', fontsize = 10, fontweight = 'bold', horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
            ax.legend(title_fontsize = 10, handles = handles[2:], labels = handle_labels[2:], fontsize = 8, 
                      bbox_to_anchor=(0.45, 0.15, 1, 0.1), frameon = False, loc = 'lower left')
            plt.text(0.48, 0.45, 'Markovian:', fontsize = 10, fontweight = 'bold', horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
            ax.add_artist(leg)
    
    ax = _axes[2][4]
    plt.sca(ax)
    colors = ['tab:blue', 'tab:red', 'tab:purple', 'tab:green', 'tab:orange', 'tab:red', 'tab:orange']
    labels = [r'Weibull, $T_{\mathrm{inf}} = 5$, $T_{\mathrm{rem}} = 7$', r'Weibull, $T_{\mathrm{inf}} = 5$, $T_{\mathrm{rem}} = 5$', r'Weibull, $T_{\mathrm{inf}} = 7$, $T_{\mathrm{rem}} = 7$', 
              r'Log-normal, $T_{\mathrm{inf}} = 5$, $T_{\mathrm{rem}} = 7$', r'Gamma, $T_{\mathrm{inf}} = 5$, $T_{\mathrm{rem}} = 7$']

    markers = ['+', 'x', 's', '^', 'D', 'o', 'v']
    for i in range(5):
        func_type, mean_inf, mean_rem, steady_level = trans_list[i]
        func_ir_name = get_func_ir_name(func_type, mean_inf, mean_rem, steady_level)
        key_name = 'errors_'+func_ir_name+'(subplot_j)'
        plt.plot(_figure_theories_data[key_name]['log_ratio'], 
                 _figure_theories_data[key_name]['loss'], 
                 markers[i], color = colors[i], fillstyle='none', markersize = 3, label = labels[i])
    leg = plt.legend(fontsize = 10, ncol = 2, loc = 'upper center', bbox_to_anchor=(0.5, 1.45), frameon = False)
    plt.axvline(0, color = 'tab:gray', linewidth = 1.5, linestyle = '--')
    for lh in leg.legendHandles: 
        lh.set_markersize(6)
    
    titles = [[None,None, None], 
              [None, None], [None, None, None, None, None]]
    idx_titles = [['a', 'b', 'c'], ['d', 'e'], ['f', 'g', 'h', 'i', 'j']]
    xlabels = [['Time (d)', 'Time (d)', 'Time (d)'],['Time (d)', r'$R_0$'],['Time (d)', 'Time (d)', 'Time (d)', 'Time (d)', r'$\ln{\eta}$']]
    ylabels = [['Fraction (%)', 'Fraction (%)', 'Fraction (%)'],['Fraction (%)', 'Fraction (%)'],
               ['Fraction (%)', 'Fraction (%)', 'Fraction (%)', 'Fraction (%)', r'$\varepsilon$']]
    xtick_amounts = [[3, 3, 3],[3, 6],[4, 4, 4, 4, 4]]
    ytick_amounts = [[3, 3, 3],[3, 5],[3, 3, 3, 3, 5]]
    xlims = [[(0, 100),(0, 100),(0, 100)],[(0, 120), (0, 5)],[(0, 150), (0, 150), (0, 75), (0, 75), (-1, 2)]]
    ylims = [[(0, 1), (0, 0.16), (0, 0.7)],[(0, 0.7), (0, 1)],[(0, 1), (0, 0.6), (0, 1), (0, 0.6), (0, 4)]]
    y_exts = [[False, False, False], [False, False] ,[True, True, True, True, False]]
    frame_types = [[0, 0, 0],[2, 2],[2, 2, 2, 2, 1]] 
    xtick_ints = [[True, True, True],[True, False],[True, True, True, True, False]]
    show_props = [[True, True, True],[True, True],[True, True, True, True, False]]
    for i, ax_arr in enumerate(_axes):
        for j, ax in enumerate(ax_arr):
            set_ax_format(ax, xlabels[i][j], ylabels[i][j], labelsize, ticksize, 18, spinetickwidth, spinewidth,
                          xtick_amounts[i][j], ytick_amounts[i][j], xlims[i][j], ylims[i][j], frame_types[i][j], 
                          titles[i][j], idx_titles[i][j], offset, _xtick_int = xtick_ints[i][j], _show_prop = show_props[i][j], _y_ext = y_exts[i][j])
    
            
    plt.subplots_adjust(left=0.06, bottom=0.06, right=0.94, top=0.96, wspace=0.25, hspace=1)
    return _fig, _axes
 
def figure_analysis(_figure_analysis_data):
    _fig = plt.figure(figsize = [2 * 8, 2 * 3])
    _axes = []
    gs = GridSpec(1, 2, figure = _fig)
    _axes.append([])
    _axes[-1].append(_fig.add_subplot(gs[0:1, 0:1]))
    _axes[-1].append(_fig.add_subplot(gs[0:1, 1:2]))

    ax = _axes[0][0]
    plt.sca(ax)
    colors = ['tab:red', 'tab:blue', 'tab:green']
    colors_marker = ['darkred', 'darkblue', 'darkgreen']
    labels = [r'$T_{\mathrm{gen}} = T_{\mathrm{rem}}$', r'$T_{\mathrm{gen}} < T_{\mathrm{rem}}$', r'$T_{\mathrm{gen}} > T_{\mathrm{rem}}$']
    markers = ['x', '+', '^']
    for i in [1, 0, 2]:
        key_name = 'explain_data_simu_' + str(i) + '(subplot_a)'
        for j in range(100):
            if j == 99:
                plt.plot(_figure_analysis_data[key_name]['time'], _figure_analysis_data[key_name]['curve_' + str(j)], 
                         linewidth = 1, alpha = 0.1, color = colors[i], label = 'Simulations, ' + labels[i])
            else:
                plt.plot(_figure_analysis_data[key_name]['time'], _figure_analysis_data[key_name]['curve_' + str(j)], 
                         linewidth = 1, alpha = 0.1, color = colors[i])
    for i in [1, 0, 2]:
        key_name = 'explain_data_calc_' + str(i) + '(subplot_a)'
        plt.plot(_figure_analysis_data[key_name]['time'], _figure_analysis_data[key_name]['curve'], 
                 markers[i], fillstyle='none', color = colors_marker[i], markersize = 12,
                 markevery = 360, markeredgewidth = 2, label = 'Fitted(average), ' + labels[i])
    fitting_c_val = np.array([_figure_analysis_data['fitting_period(subplot_a)']['c_val_'+str(i)].mean() for i in range(i)]).mean()
    fitting_c_val_init = np.array([ _figure_analysis_data['explain_data_calc_' + str(i) + '(subplot_a)']['curve'][0] for i in range(i)]).mean()
    ax.axhspan(fitting_c_val_init, fitting_c_val, 
               xmin = 0, xmax = 1, alpha = 0.2, color='tab:gray', linewidth = 0)
    
    leg = plt.legend(fontsize = 12, loc = 'lower right')
    for i in range(3):
        leg.get_lines()[i].set_alpha(1)

    ax = _axes[0][1]
    plt.sca(ax)
    axin = ax.inset_axes([0.15, 0.6, 0.4, 0.4])
    markers = {'exponent':{'weibull_i5r7': '*', 'weibull_i5r5': 'o', 'weibull_i7r7': 's', 'weibull_i7r5': '^', 'weibull_i6r6': 'D'},
               'weibull':{'weibull_i5r7': '+', 'weibull_i5r5': 'x', 'weibull_i7r7': '1', 'weibull_i7r5': '2', 'weibull_i6r6': '3'}}
    colors = {'weibull_i5r7': 'tab:blue', 'weibull_i5r5': 'tab:red', 'weibull_i7r7': 'tab:purple', 'weibull_i7r5': 'tab:green', 'weibull_i6r6': 'tab:orange'}
    labels = {'weibull_i5r7': r'Weibull, $T_{\mathrm{inf}} = 5$, $T_{\mathrm{rem}} = 7$', 
              'weibull_i5r5': r'Weibull, $T_{\mathrm{inf}} = 5$, $T_{\mathrm{rem}} = 5$', 
              'weibull_i7r7': r'Weibull, $T_{\mathrm{inf}} = 7$, $T_{\mathrm{rem}} = 7$', 
              'weibull_i7r5': r'Weibull, $T_{\mathrm{inf}} = 7$, $T_{\mathrm{rem}} = 5$',
              'weibull_i6r6': r'Weibull, $T_{\mathrm{inf}} = 6$, $T_{\mathrm{rem}} = 6$'}
    for func_type, mean_inf, mean_rem, steady_level in trans_list[:3]+trans_list[5:7]:
        func_ir_name = get_func_ir_name(func_type, mean_inf, mean_rem, steady_level)
        key_name = func_ir_name + '_r0'+'(subplot_b)'
        x = _figure_analysis_data[key_name]['log_ratio']
        y = _figure_analysis_data[key_name]['hat_r0']
        z = _figure_analysis_data[key_name]['r0']
        plt.plot(x, np.log(np.log(z) / np.log(y)), 
                 markers['exponent'][func_ir_name], color = colors[func_ir_name], fillstyle='none', label = labels[func_ir_name])
        axin.plot(np.exp(x), np.log(y), markers['exponent'][func_ir_name], color = colors[func_ir_name], fillstyle='none', markersize = 4)
    plt.hlines(0, - 1.5 - 3.5 /25, xmax = 0, color = 'tab:gray', linewidth = 2, linestyle = ':')
    plt.vlines(0, - 2 - 6 /25, ymax = 0, color = 'tab:gray', linewidth = 2, linestyle = ':') 
    plt.plot(_figure_analysis_data['curve_fitted']['log_ratio'], _figure_analysis_data['curve_fitted']['res'], color = 'black', linewidth = 3, linestyle = '--')
    plt.legend(fontsize = 12, loc = 4)
    
    set_ax_format(axin, r'$\eta$', r'$\ln{\hat{R}_0}$', 18, 18, 25, spinetickwidth, spinewidth,
                  4, 4, (0, 6), (0, 3), 2, 
                  None, None, offset, _xtick_int = False, _show_prop = False)
    
    titles = [[None, None]]
    idx_titles = [['d', 'e']]
    xlabels = [['Time (d)', r'$\ln\eta$']]
    ylabels = [['Fraction (%)', r'$\ln(\frac{\ln{R_0}}{\ln{\hat{R}_0}})$']]
    xtick_amounts = [[5, 8]]
    ytick_amounts = [[4, 4]]
    xlims = [[(0, 400),(-1.5, 2)]]
    ylims = [[(0, 0.75), (-2, 4)]]
    y_exts = [[False, False]]
    frame_types = [[0, 0]] 
    xtick_ints = [[True, False]]
    show_props = [[True, False]]
    for i, ax_arr in enumerate(_axes):
        for j, ax in enumerate(ax_arr):
            set_ax_format(ax, xlabels[i][j], ylabels[i][j], 18, 18, 25, spinetickwidth, spinewidth,
                          xtick_amounts[i][j], ytick_amounts[i][j], xlims[i][j], ylims[i][j], frame_types[i][j], 
                          titles[i][j], idx_titles[i][j], offset, _xtick_int = xtick_ints[i][j], _show_prop = show_props[i][j], _y_ext = y_exts[i][j])
          
    plt.subplots_adjust(left=0.065, bottom=0.12, right=0.96, top=0.9, wspace=0.3, hspace=0.4)
    return _fig, _axes

def figure_analysis_addition(_figure_analysis_data):
    _fig = plt.figure(figsize = [1 * 8, 2 * 3])
    _axes = []
    gs = GridSpec(1, 1, figure = _fig)
    _axes.append([])
    _axes[-1].append(_fig.add_subplot(gs[0:1, 0:1]))

    ax = _axes[0][0]
    plt.sca(ax)
    axin = ax.inset_axes([0.15, 0.6, 0.4, 0.4])
    markers = {'exponent':{'weibull_i5r7': '*', 'weibull_i5r5': 'o', 'weibull_i7r7': 's', 'lognormal_i5r7': '^', 'gamma_i5r7': 'D'},
               'weibull':{'weibull_i5r7': '+', 'weibull_i5r5': 'x', 'weibull_i7r7': '1', 'lognormal_i5r7': '2', 'gamma_i5r7': '3'}}
    colors = {'weibull_i5r7': 'tab:blue', 'weibull_i5r5': 'tab:red', 'weibull_i7r7': 'tab:purple', 'lognormal_i5r7': 'tab:green', 'gamma_i5r7': 'tab:orange'}
    labels = {'weibull_i5r7': r'Weibull, $T_{\mathrm{inf}} = 5$, $T_{\mathrm{rem}} = 7$', 
              'weibull_i5r5': r'Weibull, $T_{\mathrm{inf}} = 5$, $T_{\mathrm{rem}} = 5$', 
              'weibull_i7r7': r'Weibull, $T_{\mathrm{inf}} = 7$, $T_{\mathrm{rem}} = 7$', 
              'lognormal_i5r7': r'Log-normal, $T_{\mathrm{inf}} = 5$, $T_{\mathrm{rem}} = 7$',
              'gamma_i5r7': r'Gamma, $T_{\mathrm{inf}} = 5$, $T_{\mathrm{rem}} = 7$'}
    for func_type, mean_inf, mean_rem, steady_level in trans_list[:5]:
        func_ir_name = get_func_ir_name(func_type, mean_inf, mean_rem, steady_level)
        if func_type != 'weibull':
            key_name = func_ir_name + '_r0'+'(si_subplot_a)'
        else:
            key_name = func_ir_name + '_r0'+'(subplot_b)'
        x = _figure_analysis_data[key_name]['log_ratio']
        y = _figure_analysis_data[key_name]['hat_r0']
        z = _figure_analysis_data[key_name]['r0']
        plt.plot(x, np.log(np.log(z) / np.log(y)), 
                 markers['exponent'][func_ir_name], color = colors[func_ir_name], fillstyle='none', label = labels[func_ir_name])
        axin.plot(np.exp(x), np.log(y), markers['exponent'][func_ir_name], color = colors[func_ir_name], fillstyle='none', markersize = 4)
    plt.hlines(0, - 1.5 - 3.5 /25, xmax = 0, color = 'tab:gray', linewidth = 2, linestyle = ':')
    plt.vlines(0, - 2 - 6 /25, ymax = 0, color = 'tab:gray', linewidth = 2, linestyle = ':') 
    #plt.plot(_figure_analysis_data['curve_fitted']['log_ratio'], _figure_analysis_data['curve_fitted']['res'], color = 'black', linewidth = 3, linestyle = '--')
    plt.legend(fontsize = 12, loc = 4)
    
    set_ax_format(axin, r'$\eta$', r'$\ln{\hat{R}_0}$', 18, 18, 25, spinetickwidth, spinewidth,
                  4, 4, (0, 6), (0, 3), 2, 
                  None, None, offset, _xtick_int = False, _show_prop = False)
    
    titles = [[None, None]]
    idx_titles = [[None]]
    xlabels = [[r'$\ln\eta$']]
    ylabels = [[r'$\ln(\frac{\ln{R_0}}{\ln{\hat{R}_0}})$']]
    xtick_amounts = [[8]]
    ytick_amounts = [[4]]
    xlims = [[(-1.5, 2)]]
    ylims = [[(-2, 4)]]
    y_exts = [[False]]
    frame_types = [[0]] 
    xtick_ints = [[False]]
    show_props = [[False]]
    for i, ax_arr in enumerate(_axes):
        for j, ax in enumerate(ax_arr):
            set_ax_format(ax, xlabels[i][j], ylabels[i][j], 18, 18, 25, spinetickwidth, spinewidth,
                          xtick_amounts[i][j], ytick_amounts[i][j], xlims[i][j], ylims[i][j], frame_types[i][j], 
                          titles[i][j], idx_titles[i][j], offset, _xtick_int = xtick_ints[i][j], _show_prop = show_props[i][j], _y_ext = y_exts[i][j])
          
    plt.subplots_adjust(left=0.16, bottom=0.125, right=0.95, top=0.955, wspace=0.3, hspace=0.4)
    return _fig, _axes
   
def figure_forecasting(_figure_forecasting_data, _fitting_period_data):
    _fig = plt.figure(figsize = [16, 9])
    _axes = []
    gs = GridSpec(50, 41, figure=_fig)
    _axes.append([])
    _axes[-1].append(_fig.add_subplot(gs[0:14, 0:10]))
    _axes[-1].append(_fig.add_subplot(gs[18:32, 0:10]))
    _axes[-1].append(_fig.add_subplot(gs[36:50, 0:10]))
    _axes.append([])
    _axes[-1].append(_fig.add_subplot(gs[1:19, 12:22]))
    _axes[-1].append(_fig.add_subplot(gs[25:43, 12:22]))
    cbar_ax = _fig.add_subplot(gs[49:50, 12:22])
    _axes.append([])
    _axes[-1].append(_fig.add_subplot(gs[10:44, 24:41]))

    colors = {'simu': 'tab:blue', 'weibull': 'tab:red', 'gamma': 'tab:green', 'exponent': 'tab:orange'}
    linestyles = {'simu': '-', 'weibull': '--', 'gamma': '-.', 'exponent': ':'}
    labels = {'simu': 'Simulation', 'weibull': 'non-Markovian', 'exponent': 'Markovian'}
    
    param_settings = [r'$\ln{\alpha_{inf}} = -0.3, \ln{\alpha_{rem}} = 1.2$', 
              r'$\ln{\alpha_{inf}} = 0.45, \ln{\alpha_{rem}} = 0.45$', 
              r'$\ln{\alpha_{inf}} = 1.2, \ln{\alpha_{rem}} = -0.3$']
    data_types = ['simu', 'weibull', 'exponent']
    curve_names = {'simu': 'simulation', 'weibull': 'non-Markovian', 'exponent': 'Markovian'}
    for idx, i in enumerate(['-6_24', '9_9', '24_-6']):
        ax = _axes[0][idx]
        plt.sca(ax)
        key_name = 'curves_' + i + '(subplot_' + string.ascii_lowercase[idx] + ')'
        for data_type in data_types:
            time_line, curve_mean, curve_std = _figure_forecasting_data[key_name]['time'], _figure_forecasting_data[key_name][curve_names[data_type]],\
                                               _figure_forecasting_data[key_name][curve_names[data_type] + '_std']
            plt.plot(time_line, curve_mean, color = colors[data_type], linewidth = 2.5,
                     linestyle = linestyles[data_type], label = labels[data_type])
            plt.fill_between(time_line, curve_mean - curve_std, curve_mean + curve_std, 
                             color = colors[data_type], alpha = 0.2, linewidth = 0)
            #plt.axvline(np.array(_fitting_period_data['fitting_period']['fitting_period_' + i]).mean(), color = 'tab:gray', linewidth = 2, linestyle = ':')  
        ax.axhspan(_figure_forecasting_data[key_name]['simulation'][0], 
                   np.array(_fitting_period_data['fitting_period']['c_val_' + i]).mean(), 
                   xmin = 0, xmax = 1, alpha = 0.2, color='tab:gray', linewidth = 0)
        #ax.add_patch(plt.Rectangle((0, 0), np.array(_fitting_period_data['fitting_period']['fitting_period_' + i]).mean(), 
                      #np.array(_fitting_period_data['fitting_period']['c_val_' + i]).mean(), fill=True, facecolor = 'tab:gray', alpha = 0.3, transform=ax.transData))
        plt.text(0.5, 0.98, param_settings[idx], fontsize = 14, fontweight = 'bold', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
        if idx == 0:
            plt.legend(loc = 'lower right', fontsize = 12)
    
    grid_data= {}
    for i, fitting_type in enumerate(['exponent', 'weibull']):
        grid_data[curve_names[fitting_type]] = np.array([_figure_forecasting_data['grid_'+curve_names[fitting_type] + '(subplot_'+string.ascii_lowercase[i + 3]+')']['inf_'+str(inf_pow)] 
                                                 for inf_pow in range(-6, 25)])
    imin = min(np.min(grid_data['Markovian']), np.min(grid_data['non-Markovian']))
    imax = max(np.max(grid_data['Markovian']), np.max(grid_data['non-Markovian']))
    titles = ['Markovian', 'non-Markovian']
    for idx, fitting_type in enumerate(['exponent', 'weibull']):
        ax = _axes[1][idx]
        plt.sca(ax)
        key_name = curve_names[fitting_type]
        im = plt.imshow(grid_data[key_name].T, cmap = cmocean.cm.curl, vmin = imin, vmax = imax)
        plt.title(titles[idx], fontsize = 15)
        if idx == 1:
            plt.xlabel(r'$\ln{\alpha_{inf}}$', fontsize = 15)
        plt.ylabel(r'$\ln{\alpha_{rem}}$', fontsize = 15)
        ax.invert_yaxis()
        ax.add_patch(patches.Rectangle((0 - 0.5, 30 - 0.5), 1, 1, linewidth = 1.5, edgecolor='tab:green', facecolor='none'))
        ax.add_patch(patches.Rectangle((15 - 0.5, 15 - 0.5), 1, 1, linewidth = 1.5, edgecolor='tab:blue', facecolor='none'))
        ax.add_patch(patches.Rectangle((30 - 0.5, 0 - 0.5), 1, 1, linewidth = 1.5, edgecolor='tab:red', facecolor='none'))
        plt.xticks(np.arange(32)[::10], np.arange(-6, 25)[::10] / 20, fontsize = 12)
        plt.yticks(np.arange(32)[::10], np.arange(-6, 25)[::10] / 20, fontsize = 12)
        ax.yaxis.set_label_coords(-0.18, 0.5)

    cb = _fig.colorbar(im, cax = cbar_ax, orientation = 'horizontal')
    cb.ax.tick_params(labelsize=18, width = 1.5)
    cb.outline.set_linewidth(1.5)
    ax = _axes[2][0]
    plt.sca(ax)

    markers = {'exponent':{'weibull_i5r7': '*', 'weibull_i5r5': 'o', 'weibull_i7r7': 's', 'lognormal_i5r7': '^', 'gamma_i5r7': 'D'},
               'weibull':{'weibull_i5r7': '+', 'weibull_i5r5': 'x', 'weibull_i7r7': '1', 'lognormal_i5r7': '2', 'gamma_i5r7': '3'}}
    colors = {'exponent':{'weibull_i5r7': 'darkslateblue', 'weibull_i5r5': 'blueviolet', 'weibull_i7r7': 'cornflowerblue', 
                          'lognormal_i5r7': 'tab:blue', 'gamma_i5r7': 'midnightblue'},
              'weibull':{'weibull_i5r7': 'saddlebrown', 'weibull_i5r5': 'chocolate', 'weibull_i7r7': 'lightcoral', 
                         'lognormal_i5r7': 'tab:red', 'gamma_i5r7': 'darkred'}}
    labels = {'weibull_i5r7': r'Weibull, $T_{\mathrm{inf}} = 5$, $T_{\mathrm{rem}} = 7$', 
              'weibull_i5r5': r'Weibull, $T_{\mathrm{inf}} = 5$, $T_{\mathrm{rem}} = 5$', 
              'weibull_i7r7':r'Weibull, $T_{\mathrm{inf}} = 7$, $T_{\mathrm{rem}} = 7$', 
              'lognormal_i5r7':r'Log-normal, $T_{\mathrm{inf}} = 5$, $T_{\mathrm{rem}} = 7$',
              'gamma_i5r7':r'Gamma, $T_{\mathrm{inf}} = 5$, $T_{\mathrm{rem}} = 7$'}
    alpha_types = {'exponent': 1, 'weibull': 1}
    for fitting_type in ['exponent', 'weibull']:
        for func_type, mean_inf, mean_rem, steady_level in trans_list[:5]:
            func_ir_name = get_func_ir_name(func_type, mean_inf, mean_rem, steady_level)
            key_name = func_ir_name + '_' + curve_names[fitting_type] + '(subplot_f)'
            plt.plot(_figure_forecasting_data[key_name]['log_ratio'], _figure_forecasting_data[key_name]['res'], linestyle = 'None',
                        marker = markers[fitting_type][func_ir_name], fillstyle = 'none', color = colors[fitting_type][func_ir_name], 
                        alpha= alpha_types[fitting_type], markersize = 6, markeredgewidth = 1, label = labels[func_ir_name])
    leg = ax.legend(fontsize = 12, bbox_to_anchor=(0, 1, 1, 0.1), frameon = False, loc = 'lower left', ncol = 2)
    plt.text(0.035, 1.3, 'Markovian:', fontsize = 14, fontweight = 'bold', horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
    plt.text(0.56, 1.3, 'non-Markovian:', fontsize = 14, fontweight = 'bold', horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
    for lh in leg.legendHandles: 
        lh.set_markersize(8)
        lh.set_markeredgewidth(1.6)
        lh.set_alpha(1)
    fmts = {'covid_19': 'D', 'sars': 'o', 'h1n1': 's', 'smallpox': 'v'}
    colors = {'exponent':{'covid_19': 'tab:purple', 'sars': 'dimgray', 'h1n1': 'tab:orange', 'smallpox': 'tab:green'},
              'weibull':{'covid_19': 'purple', 'sars': 'black', 'h1n1': 'goldenrod', 'smallpox': 'darkgreen'}}
    labels = {'covid_19': {'exponent': 'Markovian, Covid-19', 'weibull': 'non-Markovian, Covid-19'},
              'sars': {'exponent': 'Markovian, SARS', 'weibull': 'non-Markovian, SARS'}, 
              'h1n1': {'exponent': 'Markovian, H1N1', 'weibull': 'non-Markovian, H1N1'}, 
              'smallpox': {'exponent': 'Markovian, Smallpox', 'weibull': 'non-Markovian, Smallpox'}}
    lines_real = []
    key_name = 'real_fitting(subplot_f)'
    for data_type in ['exponent', 'weibull']:
        for disease in ['covid_19', 'sars', 'h1n1',  'smallpox']:
            if data_type == 'exponent':
                fillcolor = 'white'
            else:
                fillcolor = colors[data_type][disease]
            real_res = _figure_forecasting_data[key_name][disease + '_' + curve_names[data_type]]
            lines_real.append(plt.errorbar(real_res[0], real_res[1], 
                         yerr = real_res[2], markersize = 8, markeredgewidth = 2, elinewidth = 2,
                         color = colors[data_type][disease], fmt=fmts[disease], capsize=6, label = labels[disease][data_type], mfc = fillcolor)[0])
    ax.yaxis.set_label_coords(-0.07, 0.5)
    ph = [plt.errorbar([],[],marker="", ls="")[0]]*2
    handles = ph[:1] + lines_real[:4] + ph[1:] + lines_real[4:]
    labels = [None, 'COVID-19', 'SARS', 'H1N1', 'Smallpox', None, 'COVID-19', 'SARS', 'H1N1', 'Smallpox']
    leg1 = ax.legend(handles = handles, labels = labels, fontsize = 12, bbox_to_anchor=(0.48, 0.65, 1, 0.1), loc = 'lower left', ncol = 2)
    plt.text(0.51, 0.87, 'Markovian:', fontsize = 12, fontweight = 'bold', horizontalalignment='left', verticalalignment='center', transform=ax.transAxes, zorder = 100)
    plt.text(0.74, 0.87, 'non-Markovian:', fontsize = 12, fontweight = 'bold', horizontalalignment='left', verticalalignment='center', transform=ax.transAxes, zorder = 100)
    ax.add_artist(leg)
    for lh in leg1.legendHandles: 
        lh.set_markersize(6)
    
    plt.axhline(0, color = 'tab:gray', linestyle = ':', linewidth = 2)
    plt.axvline(0, color = 'tab:gray', linestyle = '--', linewidth = 2)
    xlim_left, xlim_right = -1 - 3 / 25, 2.0
    ylim_lower, ylim_upper = -0.8 - 1.6 / 25, 0.8
    ax.axvspan(xlim_left, 0, ymin=0, ymax=(0 - ylim_lower) / (ylim_upper - ylim_lower), alpha=0.2, color='tab:blue')
    ax.axvspan(0, xlim_right, ymin=(0 - ylim_lower) / (ylim_upper - ylim_lower), ymax=1, alpha=0.2, color='tab:orange')
    ax.axvspan(xlim_left, 0, ymin=(0 - ylim_lower) / (ylim_upper - ylim_lower), ymax=1, alpha=0.2, color='tab:red')
    ax.axvspan(0, xlim_right, ymin=0, ymax=(0 - ylim_lower) / (ylim_upper - ylim_lower), alpha=0.2, color='tab:green')
    plt.text(0.12, 0.92, r'(i). $T_{\mathrm{gen}} < T_{\mathrm{rem}}$' + '\n' + 'overestimate', 
             color = 'tab:red', weight = 'bold', transform = ax.transAxes)
    plt.text(0.42, 0.92, r'(ii). $T_{\mathrm{gen}} > T_{\mathrm{rem}}$' + '\n' + 'overestimate', 
             color = 'tab:orange', weight = 'bold', transform = ax.transAxes)
    plt.text(0.12, 0.02, r'(iii). $T_{\mathrm{gen}} < T_{\mathrm{rem}}$' + '\n' + 'underestimate', 
             color = 'tab:blue', weight = 'bold', transform = ax.transAxes)
    plt.text(0.42, 0.02, r'(iv). $T_{\mathrm{gen}} > T_{\mathrm{rem}}$' + '\n' + 'underestimate', 
             color = 'tab:green', weight = 'bold', transform = ax.transAxes)

    titles = [[None, None, None], ['Markovian', 'non-Markovian'], [None]]
    idx_titles = [['a', 'b', 'c'], ['d', 'e'], ['f']]
    xlabels = [[None, None, 'Time (d)'], [None, r'$\ln{\alpha_{inf}}$'], [r'$\ln{\eta}$']]
    ylabels = [[None, 'Fraction (%)', None],[r'$\ln{\alpha_{rem}}$', r'$\ln{\alpha_{rem}}$'], [r'$\varepsilon^{+}$']]
    xtick_amounts = [[4, 4, 4], [4, 4], [4]]
    ytick_amounts = [[4, 4, 4], [4, 4], [5]]
    xlims = [[(0, 54),(0, 90), (0, 270)], [(0, 30), (0, 30)], [(-1, 2)]]
    ylims = [[(0, 1),(0, 0.75), (0, 0.75)], [(0, 30), (0, 30)], [(-0.8, 0.8)]]
    y_exts = [[False, False, False], [False, False], [False]]
    frame_types = [[0, 0, 0], [2, 2], [2]] 
    xtick_ints = [[True, True, True], [False, False], [False]]
    show_props = [[True, True, True], [False, False], [False]]
    set_lims = [[True, True, True], [False, False], [True]]
    for i, ax_arr in enumerate(_axes):
        for j, ax in enumerate(ax_arr):
            set_ax_format(ax, xlabels[i][j], ylabels[i][j], 18, 18, 20, spinetickwidth, spinewidth,
                          xtick_amounts[i][j], ytick_amounts[i][j], xlims[i][j], ylims[i][j], frame_types[i][j], 
                          titles[i][j], idx_titles[i][j], offset, _xtick_int = xtick_ints[i][j], _show_prop = show_props[i][j], _y_ext = y_exts[i][j], _set_lims = set_lims[i][j])
    for i in range(2):
        ax = _axes[1][i]
        ax.set_yticks(np.linspace(0, 30, 4), np.linspace(-0.3, 1.2, 4), fontsize = 18)
        ax.set_xticks(np.linspace(0, 30, 4), np.linspace(-0.3, 1.2, 4), fontsize = 18)
    _axes[-1][-1].set_yticks([-0.8, -0.4, 0.0, 0.4, 0.8], [-0.8, -0.4, 0.0, 0.4, 0.8], fontsize = 18)
    plt.subplots_adjust(left=0.050, bottom=0.08, right=0.98, top=0.96, wspace=1, hspace=0)
    return _fig, _axes

def figure_prevention(_figure_prevention_data, _fitting_period_data):
    fig = plt.figure(figsize = [15, 9])
    axes = []
    gs = GridSpec(60, 40, figure = fig)
    axes.append([])
    cbar_axes = []
    axes[-1].append(fig.add_subplot(gs[0:16, 0:8]))
    axes[-1].append(fig.add_subplot(gs[22:38, 0:8]))
    axes[-1].append(fig.add_subplot(gs[44:60, 0:8]))
    axes.append([])
    axes[-1].append(fig.add_subplot(gs[8:16, 12:22]))
    axes.append([])
    axes[-1].append(fig.add_subplot(gs[20:36, 11:16]))
    axes[-1].append(fig.add_subplot(gs[39:55, 11:16]))
    cbar_axes.append(fig.add_subplot(gs[59:60, 11:16]))
    axes.append([])
    axes[-1].append(fig.add_subplot(gs[20:36, 18:23]))
    axes[-1].append(fig.add_subplot(gs[39:55, 18:23]))
    cbar_axes.append(fig.add_subplot(gs[59:60, 18:23]))
    axes.append([])
    axes[-1].append(fig.add_subplot(gs[12:33, 26:40]))
    axes.append([])
    axes[-1].append(fig.add_subplot(gs[39:60, 26:40]))
    colors = {'no_vaccine': 'tab:blue', 'under_20': 'tab:orange', '20-49': 'tab:green', '20+': 'tab:red', '60+': 'tab:purple', 'all_ages': 'tab:brown'}
    data_types = ['simulation',  'Markovian', 'non-Markovian']
    param_settings = ['Simulation',  'Markovian', 'non-Markovian']
    ylims = [0.64, 0.64, 0.64]
    xlims = [60, 60, 60]
    titles = ['Simulation', 'Markovian','non-Markovian']    
    ages = ['under_20', '20-49', '20+', '60+', 'all_ages']
    age_names = {'no_vaccine': 'no vaccine', 'under_20': 'under 20', '20-49': '20-49', '20+': '20+', '60+': '60+', 'all_ages': 'all ages'}
    curve_names = {'simu': 'simulation', 'weibull': 'non-Markovian', 'exponent': 'Markovian'}
    with plt.style.context('fivethirtyeight'):
        for i, way in enumerate(['simu', 'exponent', 'weibull']):
            ax = axes[0][i]
            plt.sca(ax)
            plt.title(titles[i], x = 0.5, y = 0.9, fontsize = 10)
            key_name = curve_names[way] + '_curves(subplot_'+string.ascii_lowercase[i]+')' 
            for age in ages:
                age_curve = _figure_prevention_data[key_name][age]
                age_std = _figure_prevention_data[key_name][age + '_std']
                plt.plot(_figure_prevention_data[key_name]['time'], age_curve, 
                         color = colors[age], linewidth = 1.5, label = age_names[age])
                plt.fill_between(_figure_prevention_data[key_name]['time'], age_curve - age_std, age_curve + age_std, 
                         color = colors[age], linewidth = 0, alpha = 0.2)
            ax.axhspan(_figure_prevention_data[key_name]['no_vaccine'][0], 
                       np.array(_fitting_period_data['fitting_period']['c_val_9_9']).mean(), 
                       xmin = 0, xmax = 1, alpha = 0.2, color='tab:gray', linewidth = 0)
            #plt.axvline(np.array(_fitting_period_data['fitting_period']['fitting_period_9_9']).mean(), color = 'tab:gray', linewidth = 2, linestyle = ':') 
            plt.text(0.5, 0.98, param_settings[i], fontsize = 14, fontweight = 'bold', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
            if i == 0:
                plt.legend(loc = 'upper left', fontsize = 11, ncol = 1, frameon=False)
            
    ax = axes[1][0]
    plt.sca(ax)     
    vac_eff = {}   
    vac_eff_std = {}
    key_name = 'vac_eff(subplot_d)'
    for i, age in enumerate(['no_vaccine'] + ages):
        vac_eff[age] = []
        vac_eff_std[age] = []
        for data_type in data_types:
            vac_eff[age].append(_figure_prevention_data[key_name][data_type][i])
            vac_eff_std[age].append(_figure_prevention_data[key_name][data_type + '_std'][i])
   
    width = 0.1  # the width of the bars
    multiplier = 0
    for age in ['no_vaccine'] + ages:
        offset = width * multiplier
        ax.bar(np.arange(len(data_types)) + offset, vac_eff[age], width, label=age_names[age], yerr = vac_eff_std[age], color = colors[age])
        multiplier += 1
    plt.legend(fontsize = 12, ncol = 2, bbox_to_anchor=(0, 1, 1, 0.1), frameon = False, loc='lower left')
    
    for type_idx in range(2):
        grid_data= {}
        for i, way in enumerate(['exponent', 'weibull']):
            grid_data[curve_names[way]] = np.array([_figure_prevention_data['grid_'+curve_names[way] + '(subplot_'+string.ascii_lowercase[i + 4 + type_idx * 2]+')']['inf_'+str(inf_pow)] 
                                                     for inf_pow in range(-6, 25)])
        imin = min(np.min(grid_data['Markovian']), np.min(grid_data['non-Markovian']))
        imax = max(np.max(grid_data['Markovian']), np.max(grid_data['non-Markovian']))
        for i, way in enumerate(['exponent', 'weibull']):
            ax = axes[2 + type_idx][i]
            plt.sca(ax)
            plt.title(titles[i + 1], fontsize = 10)
            if type_idx == 0:
                im = plt.imshow(grid_data[curve_names[way]].T, cmap = cmocean.cm.amp, vmin = imin, vmax = imax)
            else:
                im = plt.imshow(grid_data[curve_names[way]].T, cmap = cmocean.cm.ice, vmin = imin, vmax = imax)
            plt.ylabel(r'$\ln{\alpha_{rem}}$', fontsize = 12)
            if i == 1:
                plt.xlabel(r'$\ln{\alpha_{inf}}$', fontsize = 12)
            ax.invert_yaxis()
            ax.add_patch(patches.Rectangle((15 - 0.5, 15 - 0.5), 1, 1, linewidth = 1.5, edgecolor='tab:green', facecolor='none'))
            
        fig.colorbar(im, cax = cbar_axes[type_idx], orientation = 'horizontal')

    for type_idx in range(2):
        ax = axes[4 + type_idx][0]
        plt.sca(ax)
         
        markers = {'exponent':{'weibull_i5r7': '*', 'weibull_i5r5': 'o', 'weibull_i7r7': 's', 'lognormal_i5r7': '^', 'gamma_i5r7': 'D'},
                   'weibull':{'weibull_i5r7': '+', 'weibull_i5r5': 'x', 'weibull_i7r7': '1', 'lognormal_i5r7': '2', 'gamma_i5r7': '3'}}
        colors = {'exponent':{'weibull_i5r7': 'darkslateblue', 'weibull_i5r5': 'blueviolet', 'weibull_i7r7': 'cornflowerblue', 
                              'lognormal_i5r7': 'tab:blue', 'gamma_i5r7': 'midnightblue'},
                  'weibull':{'weibull_i5r7': 'saddlebrown', 'weibull_i5r5': 'chocolate', 'weibull_i7r7': 'lightcoral', 
                             'lognormal_i5r7': 'tab:red', 'gamma_i5r7': 'darkred'}}
        labels = {'weibull_i5r7': r'Weibull, $T_{\mathrm{inf}} = 5$, $T_{\mathrm{rem}} = 7$', 
                  'weibull_i5r5': r'Weibull, $T_{\mathrm{inf}} = 5$, $T_{\mathrm{rem}} = 5$', 
                  'weibull_i7r7':r'Weibull, $T_{\mathrm{inf}} = 7$, $T_{\mathrm{rem}} = 7$', 
                  'lognormal_i5r7':r'Log-normal, $T_{\mathrm{inf}} = 5$, $T_{\mathrm{rem}} = 7$',
                  'gamma_i5r7':r'Gamma, $T_{\mathrm{inf}} = 5$, $T_{\mathrm{rem}} = 7$'}
        alpha_types = {'exponent': 1, 'weibull': 1}
        for fitting_type in ['exponent', 'weibull']:
            for func_type, mean_inf, mean_rem, steady_level in trans_list[:5]:
                func_ir_name = get_func_ir_name(func_type, mean_inf, mean_rem, steady_level)
                key_name = func_ir_name + '_' + fitting_type + '(subplot_'+string.ascii_lowercase[8 + type_idx]+')'
                plt.plot(_figure_prevention_data[key_name]['log_ratio'], _figure_prevention_data[key_name]['res'], linestyle = 'None',
                            marker = markers[fitting_type][func_ir_name], fillstyle = 'none', color = colors[fitting_type][func_ir_name], alpha= alpha_types[fitting_type],
                            markersize = 6, markeredgewidth = 1, label = labels[func_ir_name])
        if type_idx == 0:
            plt.text(0.035, 1.47, 'Markovian:', fontsize = 12, fontweight = 'bold', horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
            plt.text(0.54, 1.47, 'non-Markovian:', fontsize = 12, fontweight = 'bold', horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
            leg = ax.legend(fontsize = 9, bbox_to_anchor=(0, 1, 1, 0.1), frameon = False, loc = 'lower left', ncol = 2)
            for lh in leg.legendHandles: 
                lh.set_markersize(8)
                lh.set_markeredgewidth(1.6)
                lh.set_alpha(1)
        fmts = {'covid_19': 'D', 'sars': 'o', 'h1n1': 's', 'smallpox': 'v'}
        colors = {'exponent':{'covid_19': 'tab:purple', 'sars': 'dimgray', 'h1n1': 'tab:orange', 'smallpox': 'tab:green'},
                  'weibull':{'covid_19': 'purple', 'sars': 'black', 'h1n1': 'goldenrod', 'smallpox': 'darkgreen'}}
        lines_real = []
        for data_type in ['exponent', 'weibull']:
            for disease in ['covid_19', 'sars', 'h1n1',  'smallpox']:
                if data_type == 'exponent':
                    fillcolor = 'white'
                else:
                    fillcolor = colors[data_type][disease]
                key_name = 'real_fitting'+ '(subplot_'+string.ascii_lowercase[8 + type_idx]+')'
                real_res = _figure_prevention_data[key_name][disease + '_' + curve_names[data_type]]
                if type_idx == 0:
                    lines_real.append(plt.errorbar(real_res[0], real_res[1], 
                                 yerr = real_res[2], markersize = 8, markeredgewidth = 2, elinewidth = 2,
                                 color = colors[data_type][disease], fmt=fmts[disease], capsize=6, mfc = fillcolor)[0])
                else:
                    lines_real.append(plt.plot(real_res[0], real_res[1], 
                                 markersize = 8, markeredgewidth = 2, linestyle = 'None',
                                 color = colors[data_type][disease], marker=fmts[disease], markerfacecolor = fillcolor)[0])
        
        ax.yaxis.set_label_coords(-0.125, 0.5)
        ax.xaxis.set_label_coords(0.5, -0.1)
        if type_idx == 0:
            ph = [plt.errorbar([],[],marker="", ls="")[0]]*2
            handles = ph[:1] + lines_real[:4] + ph[1:] + lines_real[4:]
            labels = [None, 'COVID-19', 'SARS', 'H1N1', 'Smallpox', None, 'COVID-19', 'SARS', 'H1N1', 'Smallpox']
            leg1 = ax.legend(handles = handles, labels = labels, fontsize = 9, bbox_to_anchor=(0.36, 0.58, 1, 2), loc = 'lower left', ncol = 2)
            plt.text(0.39, 0.93, 'Markovian:', fontsize = 9, fontweight = 'bold', horizontalalignment='left', verticalalignment='center', transform=ax.transAxes, zorder = 100)
            plt.text(0.61, 0.93, 'non-Markovian:', fontsize = 9, fontweight = 'bold', horizontalalignment='left', verticalalignment='center', transform=ax.transAxes, zorder = 100)
            ax.add_artist(leg)
            for lh in leg1.legendHandles: 
                lh.set_markersize(6)
        plt.axvline(0, color = 'tab:gray', linestyle = '--', linewidth = 2)
        xlim_left, xlim_right = -1.0 - 3 / 25, 2.0
        plt.xlim(xlim_left, xlim_right)
        plt.ylim(ymin = 0)
        ax.axvspan(xlim_left, 0, ymin=0, ymax=1, alpha=0.2, color='tab:blue')
        ax.axvspan(0, xlim_right, ymin=0, ymax=1, alpha=0.2, color='tab:orange')
        if type_idx == 0:
            plt.text(0.12, 0.03, r'(i). $T_{\mathrm{gen}} < T_{\mathrm{rem}}$', 
                     color = 'tab:blue', weight = 'bold', transform = ax.transAxes)
            plt.text(0.38, 0.03, r'(ii). $T_{\mathrm{gen}} > T_{\mathrm{rem}}$', 
                     color = 'tab:orange', weight = 'bold', transform = ax.transAxes)
        else:
            plt.text(0.12, 0.92, r'(i). $T_{\mathrm{gen}} < T_{\mathrm{rem}}$', 
                     color = 'tab:blue', weight = 'bold', transform = ax.transAxes)
            plt.text(0.38, 0.92, r'(ii). $T_{\mathrm{gen}} > T_{\mathrm{rem}}$', 
                     color = 'tab:orange', weight = 'bold', transform = ax.transAxes)

    titles = [[None, None, None],[None], ['Markovian', 'non-Markovian'], ['Markovian', 'non-Markovian'], [None], [None]]
    idx_titles = [['a', 'b', 'c'], ['d'], ['e', 'f'], ['g', 'h'], [None], [None]]
    xlabels = [[None, None, 'Time (d)'], [None], [None, r'$\ln{\alpha_{inf}}$'], [None, r'$\ln{\alpha_{inf}}$'], [None], [r'$\ln{\eta}$']]
    ylabels = [[None, 'Fraction (%)', None], ['Fraction (%)'], [r'$\ln{\alpha_{rem}}$', r'$\ln{\alpha_{rem}}$'], 
               [None, None], [r'$\varepsilon^{*}$'], [r'$\hat{\varepsilon}$']]
    xtick_amounts = [[4, 4, 4], [3], [4, 4], [4, 4], [4], [4]]
    ytick_amounts = [[4, 4, 4], [4], [4, 4], [4, 4], [5], [5]]
    xlims = [[(0, 60),(0, 60), (0, 60)], [(0, 2)], [(0, 30), (0, 30)], [(0, 30), (0, 30)], [(-1, 2)], [(-1, 2)]]
    ylims = [[(0, 0.63),(0, 0.63), (0, 0.63)],[(0, 0.63)], [(0, 30), (0, 30)], [(0, 30), (0, 30)], [(-0.1, 1)], [(-0.1, 1)]]
    y_exts = [[False, False, False],[False], [False, False], [False, False], [True], [True]]
    frame_types = [[0, 0, 0], [2], [2, 2], [2, 2], [2], [2]] 
    xtick_ints = [[True, True, True],[False], [False, False], [False, False], [False], [False]]
    show_props = [[True, True, True],[True], [False, False], [False, False], [False], [False]]
    set_lims = [[True, True, True],[False], [False, False], [False, False], [True], [True]]
    for i, ax_arr in enumerate(axes):
        for j, ax in enumerate(ax_arr):
            set_ax_format(ax, xlabels[i][j], ylabels[i][j], 18, 18, 16, spinetickwidth, spinewidth,
                          xtick_amounts[i][j], ytick_amounts[i][j], xlims[i][j], ylims[i][j], frame_types[i][j], 
                          titles[i][j], idx_titles[i][j], offset, _xtick_int = xtick_ints[i][j], _show_prop = show_props[i][j], _y_ext = y_exts[i][j], _set_lims = set_lims[i][j])
    for type_idx in range(2):
        for i in range(2):
            ax = axes[2 + type_idx][i]
            ax.set_yticks(np.linspace(0, 30, 4), np.linspace(-0.3, 1.2, 4), fontsize = 12)
            ax.set_xticks(np.linspace(0, 30, 4), np.linspace(-0.3, 1.2, 4), fontsize = 12)
    axes[-2][-1].set_yticks([0, 0.25, 0.5, 0.75, 1.0], [0.0, 0.25, 0.5, 0.75, 1.0], fontsize = 18)
    axes[-1][-1].set_yticks([0, 0.25, 0.5, 0.75, 1.0], [0.0, 0.25, 0.5, 0.75, 1.0], fontsize = 18)
    plt.subplots_adjust(left=0.055, bottom=0.08, right=0.98, top=0.96, wspace=1, hspace=0.5)
    axes[1][0].yaxis.label.set_size(14)
    axes[1][0].set_yticks(np.linspace(0, 0.63, 4), list(map("{:.0f}".format, np.linspace(0, 0.63, 4) * 100)), fontsize = 14)
    axes[1][0].set_xticks(np.arange(3) + width * 2.5,  ['simulation',  'Markovian', 'non-Markovian'], fontsize = 10)
    axes[-2][-1].text(-0.03, 1.02, 'i', fontsize = 20, fontweight = 'bold', horizontalalignment='right', verticalalignment='bottom', transform=axes[-2][-1].transAxes)
    axes[-1][-1].text(-0.03, 1.02, 'j', fontsize = 20, fontweight = 'bold', horizontalalignment='right', verticalalignment='bottom', transform=axes[-1][-1].transAxes)
    return fig, axes
    
def si_figure_percentile(_si_figure_percentile_data):
    _fig = plt.figure(figsize = [2 * 6, 2 * 4.8])
    _axes = []
    gs = GridSpec(2, 2, figure = _fig)
    _axes.append([])
    _axes[-1].append(_fig.add_subplot(gs[0:1, 0:1]))
    _axes[-1].append(_fig.add_subplot(gs[0:1, 1:2]))
    _axes.append([])
    _axes[-1].append(_fig.add_subplot(gs[1:2, 0:1]))
    _axes[-1].append(_fig.add_subplot(gs[1:2, 1:2]))
    
    colors = ['tab:blue', 'tab:red', 'tab:purple', 'tab:green', 'tab:orange', None, None, 'tab:red', 'tab:orange']
    labels = [r'Weibull, $T_{\mathrm{inf}} = 5$, $T_{\mathrm{rem}} = 7$', r'Weibull, $T_{\mathrm{inf}} = 5$, $T_{\mathrm{rem}} = 5$', r'Weibull, $T_{\mathrm{inf}} = 7$, $T_{\mathrm{rem}} = 7$', 
              r'Log-normal, $T_{\mathrm{inf}} = 5$, $T_{\mathrm{rem}} = 7$', r'Gamma, $T_{\mathrm{inf}} = 5$, $T_{\mathrm{rem}} = 7$']
    markers = ['+', 'x', 's', '^', 'D', None, None, 'o', 'v']
    for i, percentile in enumerate([10, 90]):
        ax = _axes[0][i]
        plt.sca(ax)
        plt.axvline(0, color = 'tab:gray', linewidth = 1.5, linestyle = '--')
        for j in range(5):
            func_type, mean_inf, mean_rem, steady_level = trans_list[j]
            func_ir_name = func_ir_name = get_func_ir_name(func_type, mean_inf, mean_rem, steady_level)
            key_name = 'errors_'+ func_ir_name + '_p' + str(percentile) + '(subplot_'+string.ascii_lowercase[i]+')'
            res = _si_figure_percentile_data[key_name]
            plt.plot(res['log_ratio'], res['loss'], markers[j], color = colors[j], fillstyle='none', markersize = 6, label = labels[j])
        if i == 0:
            leg = plt.legend(fontsize = 12, ncol = 1, loc = 'upper right', frameon = False)
            for lh in leg.legendHandles: 
                lh.set_markersize(8)
    labels = [r'$\tilde{c} = 0.6$', None, None, None, None, None, None, r'$\tilde{c} = 0.4$', r'$\tilde{c} = 0.8$']        
    for i, percentile in enumerate([50, 5]):
        ax = _axes[1][i]
        plt.sca(ax)
        plt.axvline(0, color = 'tab:gray', linewidth = 1.5, linestyle = '--')
        for j in [0, 7, 8]:
            func_type, mean_inf, mean_rem, steady_level = trans_list[j]
            func_ir_name = get_func_ir_name(func_type, mean_inf, mean_rem, steady_level)
            key_name = 'errors_'+ func_ir_name + '_p' + str(percentile) + '(subplot_'+string.ascii_lowercase[i + 2]+')'
            res = _si_figure_percentile_data[key_name]
            plt.plot(res['log_ratio'], res['loss'], markers[j], color = colors[j], fillstyle='none', markersize = 6, label = labels[j])
        if i == 0:
            leg = plt.legend(fontsize = 18, loc = 'upper center', frameon = False)
            for lh in leg.legendHandles: 
                lh.set_markersize(8)
    
    titles = [['10 pencentile','90 percentile'], [r'50 percentile with different $\tilde{c}$', r'5 percentile with different $\tilde{c}$']]        
    idx_titles = [['a', 'b'], ['c', 'd']]
    xlabels = [[r'$\ln{\eta}$', r'$\ln{\eta}$'], [r'$\ln{\eta}$', r'$\ln{\eta}$']]
    ylabels = [[r'$\varepsilon$', r'$\varepsilon$'], [r'$\varepsilon$', r'$\varepsilon$']]
    xtick_amounts = [[4, 4], [4, 4]]
    ytick_amounts = [[4, 5], [5, 5]]
    xlims = [[(-1, 2),(-1, 2)], [(-1, 2), (-1, 2)]]
    ylims = [[(0, 30),(0, 4)], [(0, 4), (0, 20)]]
    y_exts = [[False, False], [False, False]]
    frame_types = [[1,1], [1, 1]] 
    xtick_ints = [[False, False], [False, False]]
    show_props = [[False, False], [False, False]]
    set_lims = [[True, True], [True, True]]
    for i, ax_arr in enumerate(_axes):
        for j, ax in enumerate(ax_arr):
            set_ax_format(ax, xlabels[i][j], ylabels[i][j], 18, 16, 18, spinetickwidth, spinewidth,
                          xtick_amounts[i][j], ytick_amounts[i][j], xlims[i][j], ylims[i][j], frame_types[i][j], 
                          titles[i][j], idx_titles[i][j], offset, _xtick_int = xtick_ints[i][j], _show_prop = show_props[i][j], _y_ext = y_exts[i][j], _set_lims = set_lims[i][j])
    plt.subplots_adjust(left=0.075, bottom=0.075, right=0.97, top=0.95, wspace=0.25, hspace=0.3)
    return _fig, _axes
                
def si_figure_markovian_time(_si_figure_markovian_time_data):
    _fig = plt.figure(figsize = [1 * 6, 1 * 4.8])
    _axes = []
    gs = GridSpec(1, 1, figure = _fig)
    _axes.append([])
    _axes[-1].append(_fig.add_subplot(gs[0:1, 0:1]))
    ax = _axes[0][0]
    plt.sca(ax)
    plt.plot(_si_figure_markovian_time_data['log_alpha'], _si_figure_markovian_time_data['m_mean_rem'], color = 'tab:blue', linewidth = 3)
    plt.hlines(7, -0.3 - 2.6 /25, xmax = 0, color = 'tab:gray', linewidth = 2, linestyle = ':')
    plt.vlines(0, 5 / 5 / 25, ymax = 7, color = 'tab:gray', linewidth = 2, linestyle = ':') 
    set_ax_format(ax, r'$\ln{\alpha_{\mathrm{gen}}}$ ($\ln{\alpha_{\mathrm{rem}}}$)',r'$T^\dagger_{\mathrm{gen}}$ ($T^\dagger_{\mathrm{rem}})$', 14, 14, 16, spinetickwidth, spinewidth,
                  2, 2, (-0.3, 1.3), (5, 10), 2, 
                  None, None, offset, _xtick_int = False, _show_prop = False, _y_ext = True, _set_lims = True)
    plt.subplots_adjust(left=0.12, bottom=0.14, right=0.92, top=0.92, wspace=0.25, hspace=0.3)
    ax.set_xticks([-0.3, 0.0, 0.3, 0.6, 0.9, 1.2],  [-0.3, 0.0, 0.3, 0.6, 0.9, 1.2], fontsize = 14)
    ax.yaxis.set_label_coords(-0.05, 0.5)
    ax.set_yticks(np.arange(5,11), np.arange(5,11), fontsize = 14)
    return _fig, _axes

def si_figure_error_percentile(_si_figure_error_percentile_data):
    _fig = plt.figure(figsize = [3 * 6, 1 * 4.8])
    _axes = []
    gs = GridSpec(1, 3, figure = _fig)
    _axes.append([])
    _axes[-1].append(_fig.add_subplot(gs[0:1, 0:1]))
    _axes[-1].append(_fig.add_subplot(gs[0:1, 1:2]))
    _axes[-1].append(_fig.add_subplot(gs[0:1, 2:3]))
    ax_idx = 0
    colors = {'weibull': 'tab:red', 'lognormal': 'tab:blue', 'gamma': 'orange'}
    for i, idx in enumerate([0, 3, 4]):
        func_type, mean_inf, mean_rem, steady_level = trans_list[idx]
        key_name = func_type+'(subplot_'+string.ascii_lowercase[i]+')'
        ax = _axes[0][ax_idx]
        ax_idx += 1
        plt.sca(ax)
        data_keys = _si_figure_error_percentile_data[key_name].keys()
        for j in data_keys:
            if j != 'percentile':
                plt.plot(_si_figure_error_percentile_data[key_name]['percentile'], _si_figure_error_percentile_data[key_name][j], color = colors[func_type], linewidth = 0.5)
    titles = [['Weibull','Log-normal', 'Gamma'], ]       
    idx_titles = [['a', 'b', 'c']]
    xlabels = [[r'Percentile $\theta$', r'Percentile $\theta$', r'Percentile $\theta$']]
    ylabels = [[r'$\varepsilon$', r'$\varepsilon$', r'$\varepsilon$']]
    xtick_amounts = [[5, 5, 5]]
    ytick_amounts = [[5, 5, 5]]
    xlims = [[(0, 100),(0, 100), (0, 100)]]
    ylims = [[(0, 1),(0, 1), (0, 1)]]
    y_exts = [[False, False, False]]
    frame_types = [[2, 2, 2]] 
    xtick_ints = [[True, True, True]]
    show_props = [[False, False, False]]
    set_lims = [[True, True, True]]
    for i, ax_arr in enumerate(_axes):
        for j, ax in enumerate(ax_arr):
            set_ax_format(ax, xlabels[i][j], ylabels[i][j], 18, 16, 20, spinetickwidth, spinewidth,
                          xtick_amounts[i][j], ytick_amounts[i][j], xlims[i][j], ylims[i][j], frame_types[i][j], 
                          titles[i][j], idx_titles[i][j], offset, _xtick_int = xtick_ints[i][j], _show_prop = show_props[i][j], _y_ext = y_exts[i][j], _set_lims = set_lims[i][j])
    plt.subplots_adjust(left=0.055, bottom=0.14, right=0.98, top=0.9, wspace=0.25, hspace=0.3)
    return _fig, _axes









def figure_prevention_opt(_figure_prevention_opt_data):
    fig = plt.figure(figsize = [15, 9])
    axes = []
    gs = GridSpec(2, 6, figure = fig)
    axes.append([])
    axes[-1].append(fig.add_subplot(gs[0:1, 0:2]))
    axes[-1].append(fig.add_subplot(gs[0:1, 2:4]))
    axes[-1].append(fig.add_subplot(gs[0:1, 4:6]))
    axes[-1].append(fig.add_subplot(gs[1:2, 1:3]))
    axes[-1].append(fig.add_subplot(gs[1:2, 3:5]))
    ylims = [0.64, 0.64, 0.64]
    xlims = [60, 60, 60]
    titles = ['Simulation', 'Markovian','non-Markovian']    
    curve_names = {'simu': 'simulation', 'weibull': 'non-Markovian', 'exponent': 'Markovian'}
    
    
    
    
    for type_idx in range(5):
        ax = axes[0][type_idx]
        plt.sca(ax)
        markers = {'exponent':{'weibull_i5r7': '*', 'weibull_i5r5': 'o', 'weibull_i7r7': 's', 'lognormal_i5r7': '^', 'gamma_i5r7': 'D'},
                   'weibull':{'weibull_i5r7': '+', 'weibull_i5r5': 'x', 'weibull_i7r7': '1', 'lognormal_i5r7': '2', 'gamma_i5r7': '3'}}
        colors = {'weibull_i5r7': 'tab:blue', 'weibull_i5r5': 'tab:red', 'weibull_i7r7': 'tab:purple', 'lognormal_i5r7': 'tab:green', 'gamma_i5r7': 'tab:orange'}
        labels = {'weibull_i5r7': r'Weibull, $T_{\mathrm{inf}} = 5$, $T_{\mathrm{rem}} = 7$', 
                  'weibull_i5r5': r'Weibull, $T_{\mathrm{inf}} = 5$, $T_{\mathrm{rem}} = 5$', 
                  'weibull_i7r7':r'Weibull, $T_{\mathrm{inf}} = 7$, $T_{\mathrm{rem}} = 7$', 
                  'lognormal_i5r7':r'Log-normal, $T_{\mathrm{inf}} = 5$, $T_{\mathrm{rem}} = 7$',
                  'gamma_i5r7':r'Gamma, $T_{\mathrm{inf}} = 5$, $T_{\mathrm{rem}} = 7$'}
        alpha_types = {'exponent': 1, 'weibull': 0.6}
        for fitting_type in ['exponent', 'weibull']:
            for func_type, mean_inf, mean_rem, steady_level in trans_list[:5]:
                func_ir_name = get_func_ir_name(func_type, mean_inf, mean_rem, steady_level)
                key_name = func_ir_name + '_' + fitting_type + '(subplot_'+string.ascii_lowercase[type_idx]+')'
                plt.plot(_figure_prevention_opt_data[key_name]['log_ratio'], _figure_prevention_opt_data[key_name]['res'], linestyle = 'None',
                            marker = markers[fitting_type][func_ir_name], fillstyle = 'none', color = colors[func_ir_name], alpha= alpha_types[fitting_type],
                            markersize = 6, markeredgewidth = 0.8, label = labels[func_ir_name])
        if type_idx == 0:
            plt.text(0.04, 1.5, 'Markovian:', fontsize = 12, fontweight = 'bold', horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
            plt.text(0.54, 1.5, 'non-Markovian:', fontsize = 12, fontweight = 'bold', horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
            leg = ax.legend(fontsize = 9, bbox_to_anchor=(0, 1, 1, 0.1), frameon = False, loc = 'lower left', ncol = 2)
            for lh in leg.legendHandles: 
                lh.set_markersize(8)
                lh.set_markeredgewidth(1)
                lh.set_alpha(1)
        fmts = {'covid_19': 'D', 'sars': 'o', 'h1n1': 's', 'smallpox': 'v'}
        colors = {'exponent':{'covid_19': 'tab:purple', 'sars': 'dimgray', 'h1n1': 'red', 'smallpox': 'blue'},
                  'weibull':{'covid_19': 'purple', 'sars': 'black', 'h1n1': 'darkred', 'smallpox': 'darkblue'}}
        lines_real = []
        for data_type in ['exponent', 'weibull']:
            for disease in ['covid_19', 'sars', 'h1n1',  'smallpox']:
                if data_type == 'exponent':
                    fillcolor = 'white'
                else:
                    fillcolor = colors[data_type][disease]
                key_name = 'real_fitting'+ '(subplot_'+string.ascii_lowercase[type_idx]+')'
                real_res = _figure_prevention_opt_data[key_name][disease + '_' + curve_names[data_type]]
                lines_real.append(plt.plot(real_res[0], real_res[1], 
                             markersize = 8, markeredgewidth = 2, linestyle = 'None',
                             color = colors[data_type][disease], marker=fmts[disease], markerfacecolor = fillcolor)[0])
        
        ax.yaxis.set_label_coords(-0.125, 0.5)
        ax.xaxis.set_label_coords(0.5, -0.1)
        ph = [plt.errorbar([],[],marker="", ls="")[0]]*2
        handles = ph[:1] + lines_real[:4] + ph[1:] + lines_real[4:]
        labels = [None, 'COVID-19', 'SARS', 'H1N1', 'Smallpox', None, 'COVID-19', 'SARS', 'H1N1', 'Smallpox']
        if type_idx == 0:
            leg1 = ax.legend(handles = handles, labels = labels, fontsize = 9, bbox_to_anchor=(0.36, 0.58, 1, 2), loc = 'lower left', ncol = 2)
            plt.text(0.4, 0.93, 'Markovian:', fontsize = 9, fontweight = 'bold', horizontalalignment='left', verticalalignment='center', transform=ax.transAxes, zorder = 100)
            plt.text(0.62, 0.93, 'non-Markovian:', fontsize = 9, fontweight = 'bold', horizontalalignment='left', verticalalignment='center', transform=ax.transAxes, zorder = 100)
            ax.add_artist(leg)
            for lh in leg1.legendHandles: 
                lh.set_markersize(6)
        plt.axvline(0, color = 'tab:gray', linestyle = '--', linewidth = 2)
        xlim_left, xlim_right = -1.0 - 3 / 25, 2.0
        plt.xlim(xlim_left, xlim_right)
        plt.ylim(ymin = 0)
        ax.axvspan(xlim_left, 0, ymin=0, ymax=1, alpha=0.2, color='tab:blue')
        ax.axvspan(0, xlim_right, ymin=0, ymax=1, alpha=0.2, color='tab:orange')
        plt.text(0.12, 0.92, r'(i). $T_{\mathrm{gen}} < T_{\mathrm{rem}}$', 
                 color = 'tab:blue', weight = 'bold', transform = ax.transAxes)
        plt.text(0.38, 0.92, r'(ii). $T_{\mathrm{gen}} > T_{\mathrm{rem}}$', 
                 color = 'tab:orange', weight = 'bold', transform = ax.transAxes)

    titles = [[None, None, None, None, None]]
    idx_titles = [['a', 'b', 'c', 'd', 'e']]
    xlabels = [[r'$\ln{\eta}$', r'$\ln{\eta}$', r'$\ln{\eta}$', r'$\ln{\eta}$', r'$\ln{\eta}$']]
    ylabels = [[r'$\hat{\varepsilon}$', r'$\hat{\varepsilon}$', r'$\hat{\varepsilon}$', r'$\hat{\varepsilon}$',r'$\hat{\varepsilon}$']]
    xtick_amounts = [[4, 4, 4, 4, 4]]
    ytick_amounts = [[5, 5, 5, 5, 5]]
    xlims = [[(-1, 2), (-1, 2), (-1, 2), (-1, 2), (-1, 2)]]
    ylims = [[(0, 1), (0, 1), (0, 1), (0, 1), (0, 1)]]
    y_exts = [[True, True, True, True, True]]
    frame_types = [[2, 2, 2, 2, 2]] 
    xtick_ints = [[False, False, False, False, False]]
    show_props = [[False, False, False, False, False]]
    set_lims = [[True, True, True, True, True]]
    for i, ax_arr in enumerate(axes):
        for j, ax in enumerate(ax_arr):
            set_ax_format(ax, xlabels[i][j], ylabels[i][j], 18, 18, 16, spinetickwidth, spinewidth,
                          xtick_amounts[i][j], ytick_amounts[i][j], xlims[i][j], ylims[i][j], frame_types[i][j], 
                          titles[i][j], idx_titles[i][j], offset, _xtick_int = xtick_ints[i][j], _show_prop = show_props[i][j], _y_ext = y_exts[i][j], _set_lims = set_lims[i][j])
    
            ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0], [0.0, 0.25, 0.5, 0.75, 1.0], fontsize = 18)
            ax.text(-0.03, 1.02, string.ascii_lowercase[j], fontsize = 20, fontweight = 'bold', horizontalalignment='right', verticalalignment='bottom', transform=axes[-1][-1].transAxes)
    return fig, axes


def si_figure_sensitivity(_si_figure_sensitivity_data):

    def set_colorbar(_fig, _ax, _im, _fontsize):
        cax = inset_axes(_ax, width="5%",   height="100%", loc='lower left',  bbox_to_anchor=(1.05, 0., 1, 1),
                         bbox_transform=_ax.transAxes, borderpad=0)
        _cbar = _fig.colorbar(_im, cax=cax)
        _cbar.ax.tick_params(labelsize=_fontsize)
        return _cbar

    fig = plt.figure(figsize = [15, 7])
    axes = []
    gs = GridSpec(7, 15, figure = fig)
    axes.append([])
    axes[-1].append(fig.add_subplot(gs[1:5, 0:4]))
    axes.append([])
    axes[-1].append(fig.add_subplot(gs[0:3, 6:10]))
    axes[-1].append(fig.add_subplot(gs[0:3, 11:15]))
    axes[-1].append(fig.add_subplot(gs[4:7, 6:10]))
    axes[-1].append(fig.add_subplot(gs[4:7, 11:15]))

    ax = axes[0][0]
    plt.sca(ax)
    original_eta_data=_si_figure_sensitivity_data['eta(subplot_a)']
    eta_data = np.array([original_eta_data['inf_pow_' + str(inf_pow)] for inf_pow in np.arange(-6, 25)])
    im = plt.imshow(eta_data.T)
    plt.gca().invert_yaxis() 
    

    locs = [0.55,0.2,0.2,0.55]
    colors = ['tab:red', 'tab:green', 'tab:blue', 'tab:orange']
    markers = ['1', 'x', '+', '2']
    plot_name = ['prime_alpha_inf(subplot_b)', 'prime_beta_inf(subplot_c)', 
                 'prime_alpha_rem(subplot_d)', 'prime_beta_rem(subplot_e)']
    inset_name = ['prime_alpha_inf_inset(subplot_b)', 'prime_beta_inf_inset(subplot_c)', 
                 'prime_alpha_rem_inset(subplot_d)', 'prime_beta_rem_inset(subplot_e)']
    curve_name = ['prime_alpha_inf', 'prime_beta_inf', 
                 'prime_alpha_rem', 'prime_beta_rem']
    for i in range(4):
        ax = axes[1][i]
        plt.sca(ax)
        plot_data = _si_figure_sensitivity_data[plot_name[i]]
        plt.plot(plot_data['eta'], plot_data[curve_name[i]], markers[i], color = colors[i], markersize = 7, markeredgewidth = 1.2)
        inset_ax = axes[1][i].inset_axes([0.12, locs[i], 0.4, 0.4])
        #plt.sca(inset_ax)
        original_inset_data=_si_figure_sensitivity_data[inset_name[i]]
        inset_data = np.array([original_inset_data['inf_pow_' + str(inf_pow)] for inf_pow in np.arange(-6, 25)])
        im = inset_ax.imshow(inset_data.T)
        inset_ax.invert_yaxis()
        inset_ax.set_xlabel(r'$\ln{\alpha_{\mathrm{inf}}}$', fontsize = 12)
        inset_ax.set_ylabel(r'$\ln{\alpha_{\mathrm{rem}}}$', fontsize = 12) 
        inset_ax.set_xticks(np.linspace(0, 30, 4), np.linspace(-0.3, 1.2, 4), fontsize = 9)
        inset_ax.set_yticks(np.linspace(0, 30, 4), np.linspace(-0.3, 1.2, 4), fontsize = 9)
        set_colorbar(fig, inset_ax, im, 9)
    titles = [[None,], [None, None, None, None,]]       
    idx_titles = [['a'], ['b', 'c', 'd', 'e']]
    xlabels = [[r'$\ln{\alpha_{\mathrm{inf}}}$'], [r'$\eta$', r'$\eta$', r'$\eta$', r'$\eta$']]
    ylabels = [[r'$\ln{\alpha_{\mathrm{rem}}}$'], 
               [r'$\partial R_0 / \partial \alpha_{\mathrm{inf}}$', r'$\partial R_0 / \partial \beta_{\mathrm{inf}}$',
                r'$\partial R_0 / \partial \alpha_{\mathrm{rem}}$', r'$\partial R_0 / \partial \beta_{\mathrm{rem}}$']]
    xtick_amounts = [[4], [8, 8, 8, 8]]
    ytick_amounts = [[4], [6, 7, 8, 6]]
    xlims = [[(0, 30)], [(0, 7), (0, 7), (0, 7), (0, 7)]]
    ylims = [[(0, 30)], [(0, 140), (-12, 0), (-175, 0), (0, 10)]]
    x_exts = [[False], [True, True, True, True]]
    y_exts = [[False], [True, True, True, True]]
    frame_types = [[2], [2, 2, 2, 2]] 
    xtick_ints = [[False], [True, True, True, True]]
    show_props = [[False], [False, False, False, False]]
    set_lims = [[False], [True, True, True, True]]
    for i, ax_arr in enumerate(axes):
        for j, ax in enumerate(ax_arr):
            set_ax_format(ax, xlabels[i][j], ylabels[i][j], 18, 16, 20, spinetickwidth, spinewidth,
                          xtick_amounts[i][j], ytick_amounts[i][j], xlims[i][j], ylims[i][j], frame_types[i][j], 
                          titles[i][j], idx_titles[i][j], offset, _xtick_int = xtick_ints[i][j], _show_prop = show_props[i][j], 
                          _x_ext = x_exts[i][j], _y_ext = y_exts[i][j], _set_lims = set_lims[i][j])
            if i == 1:
                ax.set_yticks(np.linspace(ylims[i][j][0], ylims[i][j][1], ytick_amounts[i][j]).astype(np.int), 
                              np.linspace(ylims[i][j][0], ylims[i][j][1], ytick_amounts[i][j]).astype(np.int))
    ax = axes[0][0]
    ax.set_xticks(np.linspace(0, 30, 4), np.linspace(-0.3, 1.2, 4))
    ax.set_yticks(np.linspace(0, 30, 4), np.linspace(-0.3, 1.2, 4))
    cbar = set_colorbar(fig, axes[0][0], im, 12)
    cbar.outline.set_linewidth(spinewidth)
    cbar.ax.tick_params(labelsize=18, width = spinewidth)
    plt.subplots_adjust(left=0.055, bottom=0.1, right=0.99, top=0.94, wspace=0.4, hspace=0.05)    
    return fig, axes




def si_figure_sensitivity_direct():
    from scipy.special import gamma, digamma
    def _get_mean_from_weibull(_alpha, _beta):
        return _beta * gamma(1 + 1.0 / _alpha)

    def _get_beta_from_weibull(_alpha, _mean_value):
        return _mean_value / gamma(1 + 1.0 / _alpha)


    def _calc_eta(_alpha_inf, _alpha_rem):
        return gamma((1 + _alpha_inf) / _alpha_rem) / (gamma(_alpha_inf / _alpha_rem) * gamma(1 + 1 / _alpha_rem))

    def _gamma_prime(x):
        return digamma(x) * gamma(x)

    def _gamma_prime_y(x, step = 0.001):
        return (gamma(x + step) - gamma(x)) / step

    def _derivative_func(_alpha_inf, _alpha_rem, _beta_inf = 1, _beta_rem = 1, lambda_max = 1):
        p = _alpha_inf / _alpha_rem + 1
        q = _beta_rem / _beta_inf
        gamma_val = gamma(p)
        gamma_prime_val = _gamma_prime(p)
        return [lambda_max * (q ** _alpha_inf) * gamma_prime_val / _alpha_rem + lambda_max * gamma_val * (q ** _alpha_inf) * np.log(q),
                - lambda_max * gamma_val * (q ** _alpha_inf) / _beta_inf,
                - lambda_max * (q ** _alpha_inf)* gamma_prime_val / (_alpha_rem ** 2),
                  lambda_max * gamma_val * (q ** _alpha_inf) / _beta_rem]

    alpha_inf = 3
    mean_inf = 5
    beta_inf = _get_beta_from_weibull(alpha_inf, mean_inf)
    alpha_rem = 2
    mean_rem = 7
    beta_rem = _get_beta_from_weibull(alpha_rem, mean_rem)

    eta = []
    prime_alpha_inf = []
    prime_beta_inf = []
    prime_alpha_rem = []
    prime_beta_rem = []


    for inf_pow in np.arange(-6, 25):
        eta.append([])
        prime_alpha_inf.append([])
        prime_beta_inf.append([])
        prime_alpha_rem.append([])
        prime_beta_rem.append([])
        for rem_pow in np.arange(-6, 25):
            alpha_inf = np.e ** (inf_pow * 0.05)
            alpha_rem = np.e ** (rem_pow * 0.05)
            beta_inf = _get_beta_from_weibull(alpha_inf, mean_inf)
            beta_rem = _get_beta_from_weibull(alpha_rem, mean_rem)
            eta[-1].append(_calc_eta(alpha_inf, alpha_rem))
            primes = _derivative_func(alpha_inf, alpha_rem, beta_inf, beta_rem)
            prime_alpha_inf[-1].append(primes[0])
            prime_beta_inf[-1].append(primes[1])
            prime_alpha_rem[-1].append(primes[2])
            prime_beta_rem[-1].append(primes[3])
            
    eta = np.array(eta)
    prime_alpha_inf = np.array(prime_alpha_inf)
    prime_beta_inf = np.array(prime_beta_inf)
    prime_alpha_rem = np.array(prime_alpha_rem)
    prime_beta_rem = np.array(prime_beta_rem)
    prime_params = [prime_alpha_inf, prime_beta_inf, prime_alpha_rem, prime_beta_rem]

    def set_colorbar(_fig, _ax, _im, _fontsize):
        cax = inset_axes(_ax, width="5%",   height="100%", loc='lower left',  bbox_to_anchor=(1.05, 0., 1, 1),
                         bbox_transform=_ax.transAxes, borderpad=0)
        _cbar = _fig.colorbar(_im, cax=cax)
        _cbar.ax.tick_params(labelsize=_fontsize)

    fig = plt.figure(figsize = [15, 7])
    axes = []
    gs = GridSpec(7, 15, figure = fig)
    axes.append([])
    axes[-1].append(fig.add_subplot(gs[1:5, 0:4]))
    axes.append([])
    axes[-1].append(fig.add_subplot(gs[0:3, 6:10]))
    axes[-1].append(fig.add_subplot(gs[0:3, 11:15]))
    axes[-1].append(fig.add_subplot(gs[4:7, 6:10]))
    axes[-1].append(fig.add_subplot(gs[4:7, 11:15]))

    ax = axes[0][0]
    plt.sca(ax)
    im = plt.imshow(eta.T)
    plt.gca().invert_yaxis() 
    plt.xlabel(r'$\ln{\alpha_{\mathrm{inf}}}$', fontsize = 18)
    plt.ylabel(r'$\ln{\alpha_{\mathrm{rem}}}$', fontsize = 18)
    ax.set_xticks(np.linspace(0, 30, 4), np.linspace(-0.3, 1.2, 4), fontsize = 12)
    ax.set_yticks(np.linspace(0, 30, 4), np.linspace(-0.3, 1.2, 4), fontsize = 12)
    set_colorbar(fig, axes[0][0], im, 12)
    ax.text(-0.1, 1.03, 'a', fontsize = 18, fontweight = 'bold', 
             horizontalalignment='right', verticalalignment='bottom', transform = ax.transAxes)

    locs = [0.55,0.2,0.2,0.55]
    ylabels = [r'$\partial R_0 / \partial \alpha_{\mathrm{inf}}$',
               r'$\partial R_0 / \partial \beta_{\mathrm{inf}}$',
               r'$\partial R_0 / \partial \alpha_{\mathrm{rem}}$',
               r'$\partial R_0 / \partial \beta_{\mathrm{rem}}$']
    colors = ['tab:red', 'tab:green', 'tab:blue', 'tab:orange']
    markers = ['1', '2', '+', 'x']
    title_idx = ['b', 'c', 'd', 'e']
    for i in range(4):
        ax = axes[1][i]
        plt.sca(ax)
        plt.plot(eta.reshape(31 * 31), prime_params[i].reshape(31 * 31), markers[i], color = colors[i])
        plt.xlabel(r'$\eta$', fontsize = 18)
        plt.ylabel(ylabels[i], fontsize = 18)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        inset_ax = axes[1][i].inset_axes([0.12, locs[i], 0.4, 0.4])
        #plt.sca(inset_ax)
        im = inset_ax.imshow(prime_params[i].T)
        inset_ax.invert_yaxis() 
        inset_ax.set_xlabel(r'$\ln{\alpha_{\mathrm{inf}}}$', fontsize = 12)
        inset_ax.set_ylabel(r'$\ln{\alpha_{\mathrm{rem}}}$', fontsize = 12)
        inset_ax.set_xticks(np.linspace(0, 30, 4), np.linspace(-0.3, 1.2, 4), fontsize = 9)
        inset_ax.set_yticks(np.linspace(0, 30, 4), np.linspace(-0.3, 1.2, 4), fontsize = 9)
        set_colorbar(fig, inset_ax, im, 9)
        ax.text(-0.1, 1.03, title_idx[i], fontsize = 18, fontweight = 'bold', 
                 horizontalalignment='right', verticalalignment='bottom', transform = ax.transAxes)
    plt.subplots_adjust(left=0.045, bottom=0.085, right=0.99, top=0.94, wspace=0.05, hspace=0.05)    
    return fig, axes


