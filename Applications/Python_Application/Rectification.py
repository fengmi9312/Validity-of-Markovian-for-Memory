# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 12:20:20 2023

@author: admin
"""


import numpy as np


######################################functions calulating the steady-state cumulative infected fraction##########################################
def calc_c(_r0, _c, _i_init, _contacts, _age_distribution):
    return 1 - (1 - _i_init) * np.exp(-_r0 / np.max(np.linalg.eigvals(_contacts * _age_distribution[None, :])) * (_contacts @ (_age_distribution * _c)))

def steady_c(_r0, _i_init, _contacts, _age_distribution):
    _tol = 1e-6
    _c = _i_init
    _c_prev = 0
    while True:
        _c_prev = _c
        _c = calc_c(_r0, _c, _i_init, _contacts, _age_distribution)
        if ((_c - _c_prev) ** 2).sum() ** 0.5 < _tol:
            break
    return _c, _c @ _age_distribution

############################################################################################################################


# the contacts matrix of the United States
usa_contacts = np.array([[10.10043284,  2.59381512,  1.62108103,  3.68725069,  2.01843086,
                    1.6522059 ,  1.10164538,  0.57036862],
                [ 2.59381512, 23.44951842,  3.52117847,  2.90479694,  4.17075415,
                    2.73156772,  0.91623497,  0.80602836],
                [ 1.62108103,  3.52117847, 11.41973385,  4.96882424,  3.96244646,
                    3.43683773,  1.0515248 ,  0.39471326],
                [ 3.68725069,  2.90479694,  4.96882424,  8.77711568,  5.4623625 ,
                    3.59454164,  1.62793598,  0.65586248],
                [ 2.01843086,  4.17075415,  3.96244646,  5.4623625 ,  7.44529254,
                    4.12432196,  1.34042595,  0.96830938],
                [ 1.6522059 ,  2.73156772,  3.43683773,  3.59454164,  4.12432196,
                    5.63607194,  1.6530207 ,  0.75739252],
                [ 1.10164538,  0.91623497,  1.0515248 ,  1.62793598,  1.34042595,
                    1.6530207 ,  2.83338601,  0.90734513],
                [ 0.57036862,  0.80602836,  0.39471326,  0.65586248,  0.96830938,
                    0.75739252,  0.90734513,  1.48919239]])

# the population distribution of the United States
usa_age_distribution = np.array([0.12000352, 0.12789141, 0.13925591, 0.13494838, 0.12189751, 0.12724997, 0.11627753, 0.11247577])



##################################### the parameters in this area can be customized by users ####################################################################
contacts = usa_contacts
age_distribution = usa_age_distribution
i_init = np.random.random(8) * 0.02
t_gen, t_rem, r0 = 5, 7, 2
###############################################################################################################################################################


#################################### the codes in this area generate results###########################################################################
eta = t_gen / t_rem                                                             # calculate eta
r0_real = r0 ** (eta ** 1.49)                                                   # calculate the real/rectified R_0
steady = steady_c(r0, i_init, contacts, age_distribution)                       # calculate the estimated steady state
steady_real = steady_c(r0_real, i_init, contacts, age_distribution)             # calculate the real/rectified steady state
######################################################################################################################################################


#################################### the codes in this area plot figures###########################################################################
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

fig = plt.figure(figsize=(12, 6))    
axes = [[]]
gs = GridSpec(6, 3, figure=fig)
axes[-1].append(fig.add_subplot(gs[0:4, 0:1]))
axes[-1].append(fig.add_subplot(gs[1:3, 1:2]))
axes[-1].append(fig.add_subplot(gs[1:3, 2:3]))
axes.append([])
for i in range(3):
    axes[-1].append(fig.add_subplot(gs[4:6, i:i+1]))

plt.sca(axes[0][0])
plt.imshow(contacts)
plt.gca().invert_yaxis()
plt.colorbar(shrink=0.65)
plt.sca(axes[0][1])
plt.bar(np.arange(len(age_distribution)), age_distribution, color = 'tab:green')
plt.sca(axes[0][2])
plt.bar(np.arange(len(i_init)), i_init, color = 'tab:orange')
plt.sca(axes[1][0])
plt.bar(['Estimated', 'Rectified'], [r0, r0_real], color=['tab:blue', 'tab:orange'], width = 0.6)
plt.sca(axes[1][1])
x = np.arange(len(age_distribution))  # the label locations
width = 0.4  # the width of the bars
labels = ['Estimated', 'Rectified']
for i, term in enumerate([steady[0], steady_real[0]]):
    offset = width * i
    rects = plt.bar(x + offset, term, width, label = labels[i])
plt.legend(loc='upper right')
plt.xticks(x+0.5*width, x)

plt.ylim(0,1)
plt.sca(axes[1][2])
plt.bar(['Estimated', 'Rectified'], [steady[1], steady_real[1]], color=['tab:blue', 'tab:orange'], width  = 0.6)

titles = [['Contact Matrix', 'Age Distribution', 'Initial Infected Fraction'], 
          [r'Basic Reproduction Number $R_0$', 'Group-level Cumulative Infected Fraction', 'Cumulative Infected Fraction']]
xlabels = [['Age Group Index', 'Age Group Index', 'Age Group Index'],[None, 'Age Group Index', None]]
ylabels = [['Age Group Index', 'Fraction', 'Fraction'],[r'$R_0$', 'Fraction', 'Fraction']]

for i in range(2):
    for j in range(3):
        plt.sca(axes[i][j])
        plt.title(titles[i][j])
        plt.xlabel(xlabels[i][j])
        plt.ylabel(ylabels[i][j])
plt.subplots_adjust(wspace=0.4, hspace=0.6)
#########################################################################################################################





















