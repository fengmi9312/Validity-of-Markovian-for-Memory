# -*- coding: utf-8 -*-
"""
Created on Sun Jun  5 13:46:14 2022

@author: 20481756
"""


import numpy as np
from scipy.optimize import fsolve, minimize, LinearConstraint, Bounds
from occur import calc_lambda, occur_mul, occur, gnr
import func
import scipy.special as sc
from copy import deepcopy

class sir_model:
    
    def __init__(self, i0, populations, contacts, ifrs, ylls, contacts_prop = None, tot_contacts = None, k = 1,
                 occur_inf = None, occur_rem = None, delay = 14, eta = 1, step = 1, time = 0, period_len = 1000, real_data = {}):
        
        self.__populations = np.array(populations)
        self.__contacts = {}
        self.__contacts['home'] = np.array(contacts['home'])
        self.__contacts['school'] = np.array(contacts['school'])
        self.__contacts['work'] = np.array(contacts['work'])
        self.__contacts['other_locations'] = np.array(contacts['other_locations'])
            
        if contacts_prop is None:
            self.__contacts_prop = {'home':np.array(1), 'school':np.array(1), 'work': np.array(1), 'other_locations': np.array(1)}
        else:
            self.__contacts_prop = {}
            self.__contacts_prop['home'] = np.array(contacts_prop['home'])
            self.__contacts_prop['school'] = np.array(contacts_prop['scholl'])
            self.__contacts_prop['work'] = np.array(contacts_prop['work'])
            self.__contacts_prop['other_locations'] = np.array(contacts_prop['other_locations'])
        
        self.__k = np.array(k)
        self.__group_amount = len(self.__populations)
        self.__tot_contacts = ((self.__contacts['home'] + self.__contacts['home'].T) 
                         * self.__contacts_prop['home']
                         + (self.__contacts['school'] + self.__contacts['school'].T) 
                         * self.__contacts_prop['school']
                         + (self.__contacts['work'] + self.__contacts['work'].T) 
                         * self.__contacts_prop['work']
                         + (self.__contacts['other_locations'] + self.__contacts['other_locations'].T) 
                         * self.__contacts_prop['other_locations'])
        self.__totmat = self.__tot_contacts * self.__populations[None, :] * self.__k  
        
        self.__i0 = np.array(i0)
        self.__ifrs = np.array(ifrs)
        self.__ylls = np.array(ylls)
        self.__delay = deepcopy(delay)
        self.__eta = deepcopy(eta)
        self.__inf = deepcopy(occur_inf)
        self.__rem = deepcopy(occur_rem)
        self.__step = deepcopy(step)
        self.__time0 = deepcopy(time)
        self.__period_len = deepcopy(period_len)
        
        self.__fitting_func_dict = func.srv_func_dict
        self.__real_data = deepcopy(real_data)
        if len(self.__real_data) > 0:
            for key in self.__real_data.keys():
                self.__real_data[key] = np.array(self.__real_data[key])
            if 'daily_infection' not in self.__real_data:
                self.__real_data['daily_infection'] = np.append(self.__real_data['confirmed'][0], np.diff(self.__real_data['confirmed']))
            if 'confirmed' not in self.__real_data:
                self.__real_data['confirmed'] = np.array(self.__real_data['daily_infection']).cumsum()
            if 'removal' not in self.__real_data:
                self.__real_data['removal'] = self.__real_data['recovery'] + self.__real_data['deaths']
            if 'daily_infection_arr' not in self.__real_data:
                self.__real_data['daily_infection_arr'] = np.append(self.__real_data['confirmed_arr'][0:1], 
                                                                    np.diff(self.__real_data['confirmed_arr'], axis = 0), axis = 0)
            self.__real_data['daily_infection'] = np.array(self.__real_data['daily_infection'])
            self.__real_data['daily_infection_arr'] = np.array(self.__real_data['daily_infection_arr'])
            
        self.__inf_fitting_info = {}
        self.__rem_fitting_info = {}
        self.__fitting_info = {}
        
        if self.__inf is not None and self.__rem is not None:
            self.__cum = occur_mul(self.__inf, self.__rem)
            self.return_to_init()
        
    def return_to_init(self):
        self.__time_line = [self.__time0,]
        self.__current_time_int = 0
        self.__time_int = 0
        
        self.__s = np.zeros((self.__period_len, self.__group_amount))
        self.__i_in = np.zeros((self.__period_len, self.__group_amount))
        self.__i = np.zeros((self.__period_len, self.__group_amount))
        self.__r = np.zeros((self.__period_len, self.__group_amount))
        self.__p = np.zeros((self.__period_len, self.__group_amount))
        self.__c = np.zeros((self.__period_len, self.__group_amount))
        self.__d = np.zeros((self.__period_len, self.__group_amount))
        self.__y = np.zeros((self.__period_len, self.__group_amount))
        
        self.__s[0] = 1 - self.__i0
        self.__i_in[0] = self.__i0
        self.__i[0] = self.__i0
        self.__r[0] = np.zeros(len(self.__i0))
        self.__p[0] = np.zeros(len(self.__i0))
        self.__c[0] = self.__i0
        self.__d[0] = self.__r[0] * self.__ifrs
        self.__y[0] = self.__r[0] * self.__ifrs * self.__ylls
        
        self.__j = np.zeros((self.__period_len, self.__group_amount))
        self.__j[0] = self.__calc_inf()
        
        self.__s_tot = np.zeros(self.__period_len)
        self.__i_in_tot = np.zeros(self.__period_len)
        self.__s_eff_tot = np.zeros(self.__period_len)
        self.__i_tot = np.zeros(self.__period_len)
        self.__r_tot = np.zeros(self.__period_len)
        self.__p_tot = np.zeros(self.__period_len)
        self.__c_tot = np.zeros(self.__period_len)
        self.__d_tot = np.zeros(self.__period_len)
        self.__y_tot = np.zeros(self.__period_len)
        self.__calc_tot()
        
        self.__lambda = None
        
        self.__s_eff = np.zeros((self.__period_len, self.__group_amount))
        self.__s_eff_tot = np.zeros(self.__period_len)
        self.__s_eff[0] = 1 - self.__i0
        self.__s_eff_tot[0] = self.__s_eff[self.__current_time_int] @ self.__populations
        self.__current_len = self.__period_len
        for i in range(self.__delay):
            self.spread_once()
     
    def __expand_len(self):
        self.__s = np.append(self.__s, np.zeros((self.__period_len, self.__group_amount)), axis = 0)
        self.__i_in = np.append(self.__i_in, np.zeros((self.__period_len, self.__group_amount)), axis = 0)
        self.__i = np.append(self.__i, np.zeros((self.__period_len, self.__group_amount)), axis = 0)
        self.__r = np.append(self.__r, np.zeros((self.__period_len, self.__group_amount)), axis = 0)
        self.__p = np.append(self.__p, np.zeros((self.__period_len, self.__group_amount)), axis = 0)
        self.__c = np.append(self.__c, np.zeros((self.__period_len, self.__group_amount)), axis = 0)
        self.__d = np.append(self.__d, np.zeros((self.__period_len, self.__group_amount)), axis = 0)
        self.__y = np.append(self.__y, np.zeros((self.__period_len, self.__group_amount)), axis = 0)
        self.__s_eff = np.append(self.__s_eff, np.zeros((self.__period_len, self.__group_amount)), axis = 0)
        
        self.__j = np.append(self.__j, np.zeros((self.__period_len, self.__group_amount)), axis = 0)
        
        self.__s_tot = np.append(self.__s_tot, np.zeros(self.__period_len))
        self.__i_in_tot = np.append(self.__i_in_tot, np.zeros(self.__period_len))
        self.__s_eff_tot = np.append(self.__s_eff_tot, np.zeros(self.__period_len))
        self.__i_tot = np.append(self.__i_tot, np.zeros(self.__period_len))
        self.__r_tot = np.append(self.__r_tot, np.zeros(self.__period_len))
        self.__p_tot = np.append(self.__p_tot, np.zeros(self.__period_len))
        self.__c_tot = np.append(self.__c_tot, np.zeros(self.__period_len))
        self.__d_tot = np.append(self.__d_tot, np.zeros(self.__period_len))
        self.__y_tot = np.append(self.__y_tot, np.zeros(self.__period_len))
        self.__s_eff_tot = np.append(self.__s_eff_tot, np.zeros(self.__period_len))
        
        self.__current_len = self.__current_len + self.__period_len
     
    
    def return_back(self, c_time_int, remain_v = False):
        if remain_v:
            m = 1
        else:
            m = 0
        
        if c_time_int + self.__delay + m == 0:
            self.return_to_init()
        else:
            self.__time_int = c_time_int + self.__delay + m - 1
            self.__current_time_int = c_time_int + m - 1
        
            if remain_v == False:
                self.spread_once()      
            
# the following functions calculate the spreading
###############################################################################            

    def __get_end(self, mark = True):
        if mark:
            return self.__time_int
        else:
            return self.__current_time_int
            
        
    def __get_cum_head(self, mark = True):
        tmp_head = self.__get_end(mark) - self.__cum.get_actl_len()
        if tmp_head < 0:
            return None
        else:
            return tmp_head
    
    
    def __get_cum_tail(self, mark = True):
        return min(self.__get_end(mark) + 1, self.__cum.get_actl_len())
          
          
    def __get_rem_head(self, mark = True):
        tmp_head = self.__get_end(mark) - self.__rem.get_actl_len()
        if tmp_head < 0:
            return None
        else:
            return tmp_head
    
    
    def __get_rem_tail(self, mark = True):
        return min(self.__get_end(mark) + 1, self.__rem.get_actl_len())
       
        
    def __calc_inf(self):
        return func.mulf(self.__totmat, self.__cum.get_rate()[:self.__get_cum_tail()] \
            @ np.array(self.__i_in)[self.__time_int:self.__get_cum_head():-1]) 
       
    
    def __calc_tot(self):
        self.__i_in_tot[self.__time_int] = self.__i_in[self.__time_int] @ self.__populations
        self.__s_tot[self.__time_int] = self.__s[self.__time_int] @ self.__populations
        self.__i_tot[self.__time_int] = self.__i[self.__time_int] @ self.__populations
        self.__r_tot[self.__time_int] = self.__r[self.__time_int] @ self.__populations
        self.__p_tot[self.__time_int] = self.__p[self.__time_int] @ self.__populations
        self.__c_tot[self.__time_int] = self.__c[self.__time_int] @ self.__populations
        self.__d_tot[self.__time_int] = self.__d[self.__time_int] @ self.__populations
        self.__y_tot[self.__time_int] = self.__y[self.__time_int] @ self.__populations
        if self.__time_int >= self.__delay:
            self.__s_eff_tot[self.__current_time_int] = self.__s_eff[self.__current_time_int] @ self.__populations
    
    
    def __get_eff(self):
        tmp = []
        for i in range(self.__group_amount):
            if self.__s_eff[self.__current_time_int][i] == 0:
                tmp.append(0)
            else:
                tmp.append(self.__s[self.__time_int][i] / self.__s_eff[self.__current_time_int][i])
        return np.array(tmp)
    
    
    def add_ordered_vaccine(self, vac, sttg):
        if self.__time_int == self.__current_len - 1:
            self.__expand_len()
        remaining_vac = vac
        if self.__time_int >= self.__delay:
            for i in sttg:
                if self.__s_eff[self.__current_time_int][i] > 0:
                    if self.__s_eff[self.__current_time_int][i] * self.__populations[i] > remaining_vac:
                        p_inc = remaining_vac * self.__eta * self.__s[self.__time_int][i] \
                        / self.__s_eff[self.__current_time_int][i] / self.__populations[i]
                        self.__p[self.__time_int][i] = self.__p[self.__time_int][i] + p_inc
                        self.__s[self.__time_int][i] = self.__s[self.__time_int][i] - p_inc
                        self.__s_eff[self.__current_time_int][i] = self.__s_eff[self.__current_time_int][i] - p_inc
                        remaining_vac = 0
                        break
                    else:
                        p_inc = self.__s[self.__time_int][i] * self.__eta
                        remaining_vac = remaining_vac - self.__s_eff[self.__current_time_int][i] * self.__populations[i]
                        self.__p[self.__time_int][i] = self.__p[self.__time_int][i] + p_inc
                        self.__s[self.__time_int][i] = self.__s[self.__time_int][i] - p_inc
                        self.__s_eff[self.__current_time_int][i] = self.__s_eff[self.__current_time_int][i] - p_inc
        return remaining_vac    
    
    
    def add_vaccine(self, vac_alloc = None):
        remaining_vac = 0
        if vac_alloc is not None:
            vac_alloc_tmp = deepcopy(vac_alloc)
            if self.__time_int >= self.__delay:
                for i in range(self.__group_amount):
                    if (vac_alloc_tmp[i] > self.__s_eff[self.__current_time_int][i]):
                        remaining_vac = remaining_vac + (vac_alloc_tmp[i] - self.__s_eff[self.__current_time_int][i]) * self.__populations[i]
                        vac_alloc_tmp[i] = self.__s_eff[self.__current_time_int][i]
                p_inc = vac_alloc_tmp * self.__eta * self.__get_eff()
                self.__p[self.__time_int] = self.__p[self.__time_int] + p_inc
                self.__s[self.__time_int] = self.__s[self.__time_int] - p_inc
                self.__s_eff[self.__current_time_int] = self.__s_eff[self.__current_time_int] - vac_alloc_tmp * self.__eta
        return remaining_vac
    
    
    def spread_once(self):
        if self.__time_int == self.__current_len - 1:
            self.__expand_len()
        
        i_in_tmp = self.__s[self.__time_int] * self.__j[self.__time_int]
        self.__i_in[self.__time_int + 1] = i_in_tmp
        self.__c[self.__time_int + 1] = self.__c[self.__time_int] + i_in_tmp
        self.__s[self.__time_int + 1] = self.__s[self.__time_int] - i_in_tmp
        self.__i[self.__time_int + 1] = self.__rem.get_srv()[:self.__get_rem_tail()] @ \
                                        self.__i_in[self.__time_int:self.__get_rem_head():-1]
        self.__r[self.__time_int + 1] = self.__c[self.__time_int + 1] - self.__i[self.__time_int + 1]
        self.__p[self.__time_int + 1] = self.__p[self.__time_int + 1]
        self.__d[self.__time_int + 1] = self.__r[self.__time_int + 1] * self.__ifrs
        self.__y[self.__time_int + 1] = self.__r[self.__time_int + 1] * self.__ifrs * self.__ylls
        if self.__time_int >= self.__delay:
            self.__s_eff[self.__current_time_int + 1] = self.__s_eff[self.__current_time_int] - \
                                                self.__s_eff[self.__current_time_int] * self.__j[self.__time_int - self.__delay]
            self.__current_time_int = self.__current_time_int + 1
            self.__time_line.append(self.__time_line[-1] + self.__step)
        self.__j[self.__time_int + 1] = self.__calc_inf()
        self.__time_int = self.__time_int + 1
        self.__calc_tot()
    
  
###############################################################################    

    
# the following functions are for the prediction of steady state    
#############################################################################
    
    def __calc_lambda(self):
        if self.__lambda is None:
            self.__lambda = calc_lambda(self.__inf,self.__rem)
        return self.__lambda

    def __prdt_func(self, x):
        return 1 - self.__s[0] * np.exp(-self.__totmat @ (self.__calc_lambda() * (x - self.__r[0]))) - x
    
    
    def solve_prdt(self, ctol = 1e-6):
        res = fsolve(self.__prdt_func, np.ones(self.__group_amount), xtol = ctol)
        return res
    
    
    def calc_prdt(self, **kwargs):
        ctol = kwargs.pop('ctol', 1e-6)
        lam = self.__calc_lambda()
        c_tmp = np.ones(self.__group_amount)
        while True:
            c_trans_arr = 1 - self.__s[0] * np.exp(-self.__totmat@(lam*(c_tmp-self.__r[0])))
            if np.all(abs(c_trans_arr-c_tmp) <= abs(c_tmp) * (np.ones(self.__group_amount)*ctol)):
                break
            else:
                c_tmp = c_trans_arr
        return c_trans_arr
    
    
    def __get_x_from_c(self, c, target):
        if target == 'c':
            return c @ self.__populations
        elif target == 'd':
            return c @ (self.__populations * self.__ifrs)
        elif target == 'y':
            return c @ (self.__populations * self.__ifrs * self.__ylls)
        else:
            return None
    
    
    def get_prdt_from_init(self,  **kwargs):
        target = kwargs.pop('target', 'c')
        method = kwargs.pop('method', 'solve')
        if method == 'calc':
            res = self.calc_prdt(**kwargs)
        elif method == 'solve':
            res = self.solve_prdt(**kwargs)
        else:
            pass
        return self.__get_x_from_c(res, target)  

###############################################################################


# the following functions return the details of the situations
###############################################################################

    def get_time_line(self):
        return self.__time_line
    
    def get_time(self):
        return self.__time_line[-1]
    
    def get_current_time_int(self):
        return self.__current_time_int
    
    
    def get_i_in(self):
        return self.__i_in[:self.__current_time_int + 1]
    
    def get_s(self):
        return self.__s[:self.__current_time_int + 1]
    
    def get_s_eff(self):
        return self.__s_eff[:self.__current_time_int + 1]
    
    def get_i(self):
        return self.__i[:self.__current_time_int + 1]
    
    def get_r(self):
        return self.__r[:self.__current_time_int + 1]
    
    def get_p(self):
        return self.__p[:self.__current_time_int + 1]
    
    def get_c(self):
        return self.__c[:self.__current_time_int + 1]
    
    def get_d(self):
        return self.__d[:self.__current_time_int + 1]
    
    def get_y(self):
        return self.__y[:self.__current_time_int + 1]
    
    
    def get_i_in_tot(self):
        return self.__i_in_tot[:self.__current_time_int + 1]
    
    def get_s_tot(self):
        return self.__s_tot[:self.__current_time_int + 1]
    
    def get_s_eff_tot(self):
        return self.__s_eff_tot[:self.__current_time_int + 1]
    
    def get_i_tot(self):
        return self.__i_tot[:self.__current_time_int + 1]
    
    def get_r_tot(self):
        return self.__r_tot[:self.__current_time_int + 1]
    
    def get_p_tot(self):
        return self.__p_tot[:self.__current_time_int + 1]
    
    def get_c_tot(self):
        return self.__c_tot[:self.__current_time_int + 1]
    
    def get_d_tot(self):
        return self.__d_tot[:self.__current_time_int + 1]
    
    def get_y_tot(self):
        return self.__y_tot[:self.__current_time_int + 1]
    
    def get_x_tot(self, target):
        if target == 'i_in':
            return self.get_i_in_tot()
        elif target == 's':
            return self.get_s_tot()
        elif target == 's_eff':
            return self.get_s_eff_tot()
        elif target == 'i':
            return self.get_i_tot()
        elif target == 'r':
            return self.get_r_tot()
        elif target == 'p':
            return self.get_p_tot()
        elif target == 'c':
            return self.get_c_tot()
        elif target == 'd':
            return self.get_d_tot()
        elif target == 'y':
            return self.get_y_tot()
        else:
            return None
        
    
    def getc_i_in(self):
        return self.__i_in[self.__current_time_int]
    
    def getc_s(self):
        return self.__s[self.__current_time_int]
    
    def getc_s_eff(self):
        return self.__s_eff[self.__current_time_int]
    
    def getc_i(self):
        return self.__i[self.__current_time_int]
    
    def getc_r(self):
        return self.__r[self.__current_time_int]
    
    def getc_p(self):
        return self.__p[self.__current_time_int]
    
    def getc_c(self):
        return self.__c[self.__current_time_int]
    
    def getc_d(self):
        return self.__d[self.__current_time_int]
    
    def getc_y(self):
        return self.__y[self.__current_time_int]
    
    
    def getc_i_in_tot(self):
        return self.__i_in_tot[self.__current_time_int]
    
    def getc_s_tot(self):
        return self.__s_tot[self.__current_time_int]
    
    def getc_s_eff_tot(self):
        return self.__s_eff_tot[self.__current_time_int]
    
    def getc_i_tot(self):
        return self.__i_tot[self.__current_time_int]
    
    def getc_r_tot(self):
        return self.__r_tot[self.__current_time_int]
    
    def getc_p_tot(self):
        return self.__p_tot[self.__current_time_int]
    
    def getc_c_tot(self):
        return self.__c_tot[self.__current_time_int]
    
    def getc_d_tot(self):
        return self.__d_tot[self.__current_time_int]
    
    def getc_y_tot(self):
        return self.__y_tot[self.__current_time_int]
    
    def getc_x_tot(self, target):
        if target == 'i_in':
            return self.__i_in_tot[self.__current_time_int]
        elif target == 's':
            return self.__s_tot[self.__current_time_int]
        elif target == 's_eff':
            return self.__s_eff_tot[self.__current_time_int]
        elif target == 'i':
            return self.__i_tot[self.__current_time_int]
        elif target == 'r':
            return self.__r_tot[self.__current_time_int]
        elif target == 'p':
            return self.__p_tot[self.__current_time_int]
        elif target == 'c':
            return self.__c_tot[self.__current_time_int]
        elif target == 'd':
            return self.__d_tot[self.__current_time_int]
        elif target == 'y':
            return self.__y_tot[self.__current_time_int]
        else:
            return None
    
    
    def get_populations(self):
        return self.__populations
    
    
###############################################################################    
    def get_k_from_steady(self, steady, **kwargs):
        prdt_method = kwargs.pop('prdt_method', 'calc')
        target = kwargs.pop('target', 'c')
        tol = kwargs.pop('tol',1e-6)
        
        def __get_steady_from_k(_k):
            self.__totmat = self.__tot_contacts * self.__populations[None, :] * _k   
            _res = self.get_prdt_from_init(target = target, method = prdt_method)
            self.__totmat = self.__tot_contacts * self.__populations[None, :] * self.__k   
            return _res
        
        res = None
        x0 = 1.0 / (self.__calc_lambda() * np.linalg.eig(self.__tot_contacts * self.__populations[None, :])[0].max())
        while True:
            res = func.binary_search(__get_steady_from_k, steady, x0, 2 * x0, lbound = 0, tol = tol)
            if res['success']:
                break
        self.__totmat = self.__tot_contacts * self.__populations[None, :] * self.__k
        return res['x']

    def adjust_i0(self):
        self.__i0 = np.array(self.__group_amount) * (self.__real_data['confirmed'][0] - self.__real_data['removal'][0])
        return self.__i0
        
    def fit_rem_from_data(self, init_params, **kwargs):
        srv_func = self.__fitting_func_dict[kwargs.pop('fitting_func')]
        _data_len = len(self.__real_data['daily_infection'])
        _removal_data = self.__real_data['removal'] - self.__real_data['removal'][0]
        def _loss_func(_params):
            _rem_srv = np.array([1 - srv_func(i * self.__step, *_params) for i in range(_data_len)])
            return ((np.convolve(_rem_srv, self.__real_data['daily_infection'], 'full')[:_data_len] - _removal_data)**2).sum()
        occur_length = kwargs.pop('occur_length')
        self.__rem_fitting_info['min_info'] = minimize(_loss_func, init_params, **kwargs)
        self.__rem = occur(gnr(srv_func, 'func', self.__rem_fitting_info['min_info'].x), 'srv', occur_length, self.__step)
        self.__lambda = None
        self.__lambda_arr = None
        return self.__rem_fitting_info['min_info']
        
    def fit_inf_from_data(self, init_params, **kwargs):
        srv_func = self.__fitting_func_dict[kwargs.pop('fitting_func')]
        _data_len = len(self.__real_data['daily_infection_arr'])
        occur_len = kwargs.pop('occur_length')
        _fitting_len = min(_data_len, occur_len)
        _i_in_data = self.__real_data['daily_infection_arr'][:-1].T
        _s_data = (1 - self.__real_data['confirmed_arr'][:-1]).T
        _c_data_target = self.__real_data['confirmed_arr'][1:].T
        _c_data0 = self.__real_data['confirmed_arr'][0:1].T
        def _loss_func(_params):
            _inf_rate = []
            for i in range(_fitting_len):
                _tmp_i = srv_func(i * self.__step, *_params)
                if _tmp_i == 0:
                    _inf_rate.append(1)
                else:
                    _inf_rate.append(1.0 - srv_func((i + 1) * self.__step, *_params) / _tmp_i)
            _cum_rate = np.array(_inf_rate) * self.__rem.get_srv()[:_fitting_len]
            _inf_tmp = np.array([np.convolve(_cum_rate, _i_in_data[i], 'full')[:_data_len - 1] for i in range(self.__group_amount)]) 
            _c_calc = ((1 - ((1 - _inf_tmp).reshape((self.__group_amount, _data_len - 1, 1)) ** 
                      self.__totmat.T.reshape((self.__group_amount, 1, self.__group_amount))).prod(axis = 0)).T 
                      * _s_data).cumsum(axis = 1) + _c_data0
            return ((self.__populations @ (_c_calc - _c_data_target))**2).sum()
        self.__inf_fitting_info['min_info'] = minimize(_loss_func, init_params, **kwargs)
        self.__inf = occur(gnr(srv_func, 'func', self.__inf_fitting_info['min_info'].x), 'srv', occur_len, self.__step)
        self.__cum = occur_mul(self.__inf, self.__rem)
        self.return_to_init()
        self.__lambda = self.__calc_lambda()
        return self.__inf_fitting_info['min_info']
    
    def fit_cum_from_data(self, init_params, **kwargs):
        srv_func_inf = self.__fitting_func_dict[kwargs.pop('fitting_func_inf')]
        srv_func_rem = self.__fitting_func_dict[kwargs.pop('fitting_func_rem')]
        _param_div = kwargs.pop('param_div')
        _data_len = len(self.__real_data['daily_infection_arr'])
        occur_len = kwargs.pop('occur_length')
        _fitting_len = min(_data_len, occur_len)
        _i_in_data = self.__real_data['daily_infection_arr'][:-1].T
        _c_data0 = self.__real_data['confirmed_arr'][0:1].T
        _s_data = (1 - self.__real_data['confirmed_arr'][:-1]).T
        _s_data_target = 1 - (self.__real_data['confirmed_arr'][1:] @ self.__populations)
        _r_data_target = self.__real_data['removal'][1:]
        def _loss_func(_params):
            _inf_rate = np.zeros(_fitting_len)
            _rem_srv = np.zeros(_fitting_len)
            for i in range(_fitting_len):
                _tmp_i = srv_func_inf(i * self.__step, *_params[:_param_div])
                _rem_srv[i] = srv_func_rem(i * self.__step, *_params[_param_div:])
                if _tmp_i == 0:
                    _inf_rate[i] = 1
                else:
                    _inf_rate[i] = 1.0 - srv_func_inf((i + 1) * self.__step, *_params[:_param_div]) / _tmp_i
            _cum_rate = _inf_rate * _rem_srv
            _inf_tmp = np.array([np.convolve(_cum_rate, _i_in_data[i], 'full')[:_data_len - 1] for i in range(self.__group_amount)]) 
            _s_calc = self.__populations @ (1 - (((1 - ((1 - _inf_tmp).reshape((self.__group_amount, _data_len - 1, 1)) ** 
                      self.__totmat.T.reshape((self.__group_amount, 1, self.__group_amount))).prod(axis = 0)).T 
                      * _s_data).cumsum(axis = 1) + _c_data0))
            _r_calc = np.convolve(1.0 - _rem_srv, self.__real_data['daily_infection'], 'full')[1:_data_len]
            return ((_s_calc - _s_data_target)**2).sum() + ((_r_calc - _r_data_target)**2).sum()
        self.__fitting_info['min_info'] = minimize(_loss_func, init_params, **kwargs)
        self.__rem = occur(gnr(srv_func_rem, 'func', self.__fitting_info['min_info'].x[:_param_div]), 'srv', occur_len, self.__step)
        self.__inf = occur(gnr(srv_func_inf, 'func', self.__fitting_info['min_info'].x[_param_div:]), 'srv', occur_len, self.__step)
        self.__cum = occur_mul(self.__inf, self.__rem)
        self.return_to_init()
        self.__lambda = self.__calc_lambda()
        return self.__fitting_info['min_info']
    
    
    def get_vac_x(self, vac_alloc, vac_freq, interval, consequent_step, target = 'c'):
        c_time_int = self.__current_time_int
        for i in range(vac_freq):
            self.add_vaccine(vac_alloc / vac_freq)
            for j in range(interval):
                self.spread_once()
        for i in range(consequent_step):
            self.spread_once()
        x = self.getc_x_tot(target)
        self.return_back(c_time_int)
        return x
    
    def set_i_in_init(self, i_in_data):
        self.__i0 = i_in_data[0]
        self.__time_line = [self.__time0,]
        self.__current_time_int = 0
        self.__time_int = 0
        
        self.__s = np.zeros((self.__period_len, self.__group_amount))
        self.__i_in = np.zeros((self.__period_len, self.__group_amount))
        self.__i = np.zeros((self.__period_len, self.__group_amount))
        self.__r = np.zeros((self.__period_len, self.__group_amount))
        self.__p = np.zeros((self.__period_len, self.__group_amount))
        self.__c = np.zeros((self.__period_len, self.__group_amount))
        self.__d = np.zeros((self.__period_len, self.__group_amount))
        self.__y = np.zeros((self.__period_len, self.__group_amount))
        
        self.__s[0] = 1 - self.__i0
        self.__i_in[0] = self.__i0
        self.__i[0] = self.__i0
        self.__r[0] = np.zeros(len(self.__i0))
        self.__p[0] = np.zeros(len(self.__i0))
        self.__c[0] = self.__i0
        self.__d[0] = self.__r[0] * self.__ifrs
        self.__y[0] = self.__r[0] * self.__ifrs * self.__ylls
        
        self.__j = np.zeros((self.__period_len, self.__group_amount))
        self.__j[0] = self.__calc_inf()
        
        self.__s_tot = np.zeros(self.__period_len)
        self.__i_in_tot = np.zeros(self.__period_len)
        self.__s_eff_tot = np.zeros(self.__period_len)
        self.__i_tot = np.zeros(self.__period_len)
        self.__r_tot = np.zeros(self.__period_len)
        self.__p_tot = np.zeros(self.__period_len)
        self.__c_tot = np.zeros(self.__period_len)
        self.__d_tot = np.zeros(self.__period_len)
        self.__y_tot = np.zeros(self.__period_len)
        self.__calc_tot()
        
        self.__lambda = None
        
        self.__s_eff = np.zeros((self.__period_len, self.__group_amount))
        self.__s_eff_tot = np.zeros(self.__period_len)
        self.__s_eff[0] = 1 - self.__i0
        self.__s_eff_tot[0] = self.__s_eff[self.__current_time_int] @ self.__populations
        self.__current_len = self.__period_len

        for _i_in in i_in_data[1:]:
            if self.__time_int == self.__current_len - 1:
                self.__expand_len()
            self.__i_in[self.__time_int + 1] = _i_in
            self.__c[self.__time_int + 1] = self.__c[self.__time_int] + _i_in
            self.__s[self.__time_int + 1] = self.__s[self.__time_int] - _i_in
            self.__i[self.__time_int + 1] = self.__rem.get_srv()[:self.__get_rem_tail()] @ \
                                            self.__i_in[self.__time_int:self.__get_rem_head():-1]
            self.__r[self.__time_int + 1] = self.__c[self.__time_int + 1] - self.__i[self.__time_int + 1]
            self.__p[self.__time_int + 1] = self.__p[self.__time_int + 1]
            self.__d[self.__time_int + 1] = self.__r[self.__time_int + 1] * self.__ifrs
            self.__y[self.__time_int + 1] = self.__r[self.__time_int + 1] * self.__ifrs * self.__ylls
            if self.__time_int >= self.__delay:
                self.__s_eff[self.__current_time_int + 1] = self.__s_eff[self.__current_time_int] - \
                                                    self.__s_eff[self.__current_time_int] * self.__j[self.__time_int - self.__delay]
                self.__current_time_int = self.__current_time_int + 1
                self.__time_line.append(self.__time_line[-1] + self.__step)
            self.__j[self.__time_int + 1] = self.__calc_inf()
            self.__time_int = self.__time_int + 1
            self.__calc_tot()
        for i in range(self.__delay):
            self.spread_once()

def calc_k(_alpha_inf, _beta_inf, _alpha_rem, _beta_rem, _occur_length, _steady_prdt, calc_params):
    _vgnr_inf = gnr(func.srv_weibull_scale, 'func', (_alpha_inf, _beta_inf))
    _vgnr_rem = gnr(func.srv_weibull_scale, 'func', (_alpha_rem, _beta_rem))
    _calc_params = deepcopy(calc_params)
    _calc_params['occur_inf'] = occur(_vgnr_inf, vgnr_type = 'srv', length = _occur_length, step = _calc_params['step'])
    _calc_params['occur_rem'] = occur(_vgnr_rem, vgnr_type = 'srv', length = _occur_length, step = _calc_params['step'])
    k_calc_obtain = sir_model(**_calc_params)
    _k = k_calc_obtain.get_k_from_steady(_steady_prdt, target = 'c', prdt_method = 'calc', tol = 1e-10)
    del k_calc_obtain
    return _k

def calc_general_k(_survivals_inf, _survivals_rem, _steady_prdt, calc_params):
    _vgnr_inf = gnr(_survivals_inf, 'arr', None)
    _vgnr_rem = gnr(_survivals_rem, 'arr', None)
    _calc_params = deepcopy(calc_params)
    _calc_params['occur_inf'] = occur(_vgnr_inf, vgnr_type = 'srv', length = len(_survivals_inf) - 1, step = _calc_params['step'])
    _calc_params['occur_rem'] = occur(_vgnr_rem, vgnr_type = 'srv', length = len(_survivals_rem) - 1, step = _calc_params['step'])
    k_calc_obtain = sir_model(**_calc_params)
    _k = k_calc_obtain.get_k_from_steady(_steady_prdt, target = 'c', prdt_method = 'calc', tol = 1e-10)
    del k_calc_obtain
    return _k


import metapopulation_simulation as ms
def get_fitting_data(_idx, _spreading_params, _k, _amount_seeds, 
             simu_structure_params, calc_params, params_fitting, 
             fitting_level= 0.1, occur_length = 4000):
    
    ################ initialize the simulation model ##########################
    _simu_structure_params = deepcopy(simu_structure_params)
    _simu_structure_params['k'] = _k
    simu_once = ms.simu(**_simu_structure_params)
    simu_once.set_generator_seed(_idx)
    simu_once.set_spreading_func('weibull', 'weibull')
    simu_once.set_total_spreading_params(_spreading_params[:2],
                                         _spreading_params[2:])
    simu_once.set_amount_seeds(_amount_seeds)
    ###########################################################################
    
    #################### get the simulation data ##############################
    _simu_once_data = {'s': [], 'i': [], 'r': [], 'w': [], 'c': [], 'd': [], 'y': []}
    _simu_once_arr_data = {'c': []}
    for _target in _simu_once_data.keys():
        _simu_once_data[_target].append(simu_once.get_x_amount(_target) / simu_once.get_node_amount())
    _simu_once_arr_data['c'].append(np.array(simu_once.get_x_amount_arr('c')) / np.array(simu_once.get_population_amounts()))
    while True:
        if simu_once.get_i_amount() == 0:
            break
        simu_once.spread_once()
        for _target in _simu_once_data.keys():
            _simu_once_data[_target].append(simu_once.get_x_amount(_target) / simu_once.get_node_amount())
        _simu_once_arr_data['c'].append(np.array(simu_once.get_x_amount_arr('c')) / np.array(simu_once.get_population_amounts()))
    for _target in _simu_once_data.keys():
        _simu_once_data[_target] = np.array(_simu_once_data[_target])
    _simu_once_arr_data['c'] = np.array(_simu_once_arr_data['c'])
    del simu_once
    ###########################################################################
    
    ##################### get the simulation data for fitting #################
    _fitting_data = {'confirmed': [], 'removal': [], 'confirmed_arr': []}
    _fitting_end_c = (_simu_once_data['c'][-1] - _simu_once_data['c'][0]) * fitting_level + _simu_once_data['c'][0]
    for _i in range(len(_simu_once_data['c'])):
        _fitting_data['confirmed'].append(_simu_once_data['c'][_i])
        _fitting_data['removal'].append(_simu_once_data['r'][_i])
        _fitting_data['confirmed_arr'].append(_simu_once_arr_data['c'][_i])
        if _fitting_data['confirmed'][-1] > _fitting_end_c:
            break
    _fitting_data['confirmed'] = np.array(_fitting_data['confirmed'])
    _fitting_data['removal'] = np.array(_fitting_data['removal'])
    _fitting_data['confirmed_arr'] = np.array(_fitting_data['confirmed_arr'])
    _fitting_length = len(_fitting_data['confirmed'])
    _i_in_data = np.append(_fitting_data['confirmed_arr'][0:1], 
                           np.diff(_fitting_data['confirmed_arr'], axis = 0), axis = 0)
    ###########################################################################
    
    ################## initialize fitting model ###############################
    _calc_params = deepcopy(calc_params)
    _calc_params['real_data'] = _fitting_data
    _calc_params['occur_inf'] = None
    _calc_params['occur_rem'] = None
    _calc_params['k'] = _k
    fitting_once = [sir_model(**_calc_params) for i in range(3)]
    _init_params_rem = [np.array([1.0, ]), np.array([1.0, 1.0]), np.array([1.0, 1.0])]
    _init_params_inf = [np.array([1.0, ]), np.array([1.0, 1.0]), np.array([1.0, 1.0])]
    _fitting_func = ['exponent', 'weibull_rate', 'gamma_rate']
    _bounds = [Bounds(np.ones(1) * 0.001, np.ones(1) * np.inf),
               Bounds(np.ones(2) * 0.001, np.ones(2) * np.inf),
               Bounds(np.ones(2) * 0.001, np.ones(2) * np.inf)]
    _rem_fitting_res = [None, None, None]
    _inf_fitting_res = [None, None, None]
    ###########################################################################
    
    ######################## fitting ##########################################
    for i in range(3):
        _init_params = _init_params_rem[i]
        while True:
            _rem_fitting_res[i] =\
            fitting_once[i].fit_rem_from_data(init_params = _init_params, 
                                              occur_length = occur_length,  
                                              fitting_func = _fitting_func[i], 
                                              bounds = _bounds[i], **params_fitting)
            if _rem_fitting_res[i].success:
                break
            else:
                _init_params = _init_params_rem[i] * np.random.uniform(0.7, 1.3) 
        _init_params = _init_params_inf[i]
        while True:     
            _inf_fitting_res[i] = \
            fitting_once[i].fit_inf_from_data(init_params = _init_params, 
                                              occur_length = occur_length, 
                                              fitting_func = _fitting_func[i], 
                                              bounds = _bounds[i],  **params_fitting)
            if _inf_fitting_res[i].success:
                break
            else:
                _init_params = _init_params_inf[i] * np.random.uniform(0.7, 1.3)
    ###########################################################################
    
    ###################### return data ###############################        
    _params = {'exponent_inf': _inf_fitting_res[0].x[0], 'exponent_rem': _rem_fitting_res[0].x[0], 
               'weibull_alpha_inf':_inf_fitting_res[1].x[0], 'weibull_beta_inf':1.0 / _inf_fitting_res[1].x[1], 
               'weibull_alpha_rem':_rem_fitting_res[1].x[0], 'weibull_beta_rem':1.0 / _rem_fitting_res[1].x[1], 
               'gamma_alpha_inf':_inf_fitting_res[2].x[0], 'gamma_beta_inf':1.0 / _inf_fitting_res[2].x[1], 
               'gamma_alpha_rem':_rem_fitting_res[2].x[0], 'gamma_beta_rem':1.0 / _rem_fitting_res[2].x[1],
               'fitting_length': _fitting_length}
    return {'params':_params, 'i_in_data': _i_in_data, 'simu_curves': _simu_once_data}



def get_general_fitting_data(_idx, _spreading_survivals, _k, _amount_seeds, 
             simu_structure_params, calc_params, params_fitting, 
             fitting_level= 0.1, occur_length = 4000):
    
    ################ initialize the simulation model ##########################
    _simu_structure_params = deepcopy(simu_structure_params)
    _simu_structure_params['k'] = _k
    simu_once = ms.simu(**_simu_structure_params)
    simu_once.set_generator_seed(_idx)
    simu_once.set_spreading_func('general', 'general')
    simu_once.set_total_spreading_params([_spreading_survivals[2],], [_spreading_survivals[2],])
    simu_once.set_total_spreading_survivals(_spreading_survivals[0], _spreading_survivals[1])
    simu_once.set_amount_seeds(_amount_seeds)
    ###########################################################################
    
    #################### get the simulation data ##############################
    _simu_once_data = {'s': [], 'i': [], 'r': [], 'w': [], 'c': [], 'd': [], 'y': []}
    _simu_once_arr_data = {'c': []}
    for _target in _simu_once_data.keys():
        _simu_once_data[_target].append(simu_once.get_x_amount(_target) / simu_once.get_node_amount())
    _simu_once_arr_data['c'].append(np.array(simu_once.get_x_amount_arr('c')) / np.array(simu_once.get_population_amounts()))
    while True:
        if simu_once.get_i_amount() == 0:
            break
        simu_once.spread_once()
        for _target in _simu_once_data.keys():
            _simu_once_data[_target].append(simu_once.get_x_amount(_target) / simu_once.get_node_amount())
        _simu_once_arr_data['c'].append(np.array(simu_once.get_x_amount_arr('c')) / np.array(simu_once.get_population_amounts()))
    for _target in _simu_once_data.keys():
        _simu_once_data[_target] = np.array(_simu_once_data[_target])
    _simu_once_arr_data['c'] = np.array(_simu_once_arr_data['c'])
    del simu_once
    ###########################################################################
    
    ##################### get the simulation data for fitting #################
    _fitting_data = {'confirmed': [], 'removal': [], 'confirmed_arr': []}
    _fitting_end_c = (_simu_once_data['c'][-1] - _simu_once_data['c'][0]) * fitting_level + _simu_once_data['c'][0]
    for _i in range(len(_simu_once_data['c'])):
        _fitting_data['confirmed'].append(_simu_once_data['c'][_i])
        _fitting_data['removal'].append(_simu_once_data['r'][_i])
        _fitting_data['confirmed_arr'].append(_simu_once_arr_data['c'][_i])
        if _fitting_data['confirmed'][-1] > _fitting_end_c:
            break
    _fitting_data['confirmed'] = np.array(_fitting_data['confirmed'])
    _fitting_data['removal'] = np.array(_fitting_data['removal'])
    _fitting_data['confirmed_arr'] = np.array(_fitting_data['confirmed_arr'])
    _fitting_length = len(_fitting_data['confirmed'])
    _i_in_data = np.append(_fitting_data['confirmed_arr'][0:1], 
                           np.diff(_fitting_data['confirmed_arr'], axis = 0), axis = 0)
    ###########################################################################
    
    ################## initialize fitting model ###############################
    _calc_params = deepcopy(calc_params)
    _calc_params['real_data'] = _fitting_data
    _calc_params['occur_inf'] = None
    _calc_params['occur_rem'] = None
    _calc_params['k'] = _k
    fitting_once = [sir_model(**_calc_params) for i in range(3)]
    _init_params_rem = [np.array([1.0, ]), np.array([1.0, 1.0]), np.array([1.0, 1.0])]
    _init_params_inf = [np.array([1.0, ]), np.array([1.0, 1.0]), np.array([1.0, 1.0])]
    _fitting_func = ['exponent', 'weibull_rate']
    _bounds = [Bounds(np.ones(1) * 0.001, np.ones(1) * np.inf),
               Bounds(np.ones(2) * 0.001, np.ones(2) * np.inf)]
    _rem_fitting_res = [None, None]
    _inf_fitting_res = [None, None]
    ###########################################################################
    
    ######################## fitting ##########################################
    for i in range(2):
        _init_params = _init_params_rem[i]
        while True:
            _rem_fitting_res[i] =\
            fitting_once[i].fit_rem_from_data(init_params = _init_params, 
                                              occur_length = occur_length,  
                                              fitting_func = _fitting_func[i], 
                                              bounds = _bounds[i], **params_fitting)
            if _rem_fitting_res[i].success:
                break
            else:
                _init_params = _init_params_rem[i] * np.random.uniform(0.7, 1.3) 
        _init_params = _init_params_inf[i]
        while True:     
            _inf_fitting_res[i] = \
            fitting_once[i].fit_inf_from_data(init_params = _init_params, 
                                              occur_length = occur_length, 
                                              fitting_func = _fitting_func[i], 
                                              bounds = _bounds[i],  **params_fitting)
            if _inf_fitting_res[i].success:
                break
            else:
                _init_params = _init_params_inf[i] * np.random.uniform(0.7, 1.3)
    ###########################################################################
    
    ###################### return data ###############################        
    _params = {'exponent_inf': _inf_fitting_res[0].x[0], 'exponent_rem': _rem_fitting_res[0].x[0], 
               'weibull_alpha_inf':_inf_fitting_res[1].x[0], 'weibull_beta_inf':1.0 / _inf_fitting_res[1].x[1], 
               'weibull_alpha_rem':_rem_fitting_res[1].x[0], 'weibull_beta_rem':1.0 / _rem_fitting_res[1].x[1], 
               'fitting_length': _fitting_length}
    return {'params':_params, 'i_in_data': _i_in_data, 'simu_curves': _simu_once_data}


def get_fitting_data_with_real(fitting_data, calc_params, params_fitting, occur_length = 4000):
    
    ################## initialize fitting model ###############################
    _calc_params = deepcopy(calc_params)
    _calc_params['real_data'] = fitting_data
    _calc_params['occur_inf'] = None
    _calc_params['occur_rem'] = None
    fitting_once = [sir_model(**_calc_params) for i in range(3)]
    _init_params_rem = [np.array([1.0, ]), np.array([1.0, 1.0]), np.array([1.0, 1.0])]
    _init_params_inf = [np.array([1.0, ]), np.array([1.0, 1.0]), np.array([1.0, 1.0])]
    _fitting_func = ['exponent', 'weibull_rate', 'gamma_rate']
    _bounds = [Bounds(np.ones(1) * 0.001, np.ones(1) * np.inf),
               Bounds(np.ones(2) * 0.001, np.ones(2) * np.inf),
               Bounds(np.ones(2) * 0.001, np.ones(2) * np.inf)]
    _rem_fitting_res = [None, None, None]
    _inf_fitting_res = [None, None, None]
    ###########################################################################
    
    ######################## fitting ##########################################
    for i in range(3):
        _init_params = _init_params_rem[i]
        while True:
            _rem_fitting_res[i] =\
            fitting_once[i].fit_rem_from_data(init_params = _init_params, 
                                              occur_length = occur_length,  
                                              fitting_func = _fitting_func[i], 
                                              bounds = _bounds[i], **params_fitting)
            if _rem_fitting_res[i].success:
                break
            else:
                _init_params = _init_params_rem[i] * np.random.uniform(0.7, 1.3) 
        _init_params = _init_params_inf[i]
        while True:     
            _inf_fitting_res[i] = \
            fitting_once[i].fit_inf_from_data(init_params = _init_params, 
                                              occur_length = occur_length, 
                                              fitting_func = _fitting_func[i], 
                                              bounds = _bounds[i],  **params_fitting)
            if _inf_fitting_res[i].success:
                break
            else:
                _init_params = _init_params_inf[i] * np.random.uniform(0.7, 1.3)
    ###########################################################################
    
    ###################### return data ###############################        
    _params = {'exponent_inf': _inf_fitting_res[0].x[0], 'exponent_rem': _rem_fitting_res[0].x[0], 
               'weibull_alpha_inf':_inf_fitting_res[1].x[0], 'weibull_beta_inf':1.0 / _inf_fitting_res[1].x[1], 
               'weibull_alpha_rem':_rem_fitting_res[1].x[0], 'weibull_beta_rem':1.0 / _rem_fitting_res[1].x[1], 
               'gamma_alpha_inf':_inf_fitting_res[2].x[0], 'gamma_beta_inf':1.0 / _inf_fitting_res[2].x[1], 
               'gamma_alpha_rem':_rem_fitting_res[2].x[0], 'gamma_beta_rem':1.0 / _rem_fitting_res[2].x[1]}
    return _params


def get_fitting_data_with_cum(_idx, _spreading_params, _k, _amount_seeds, 
             simu_structure_params, calc_params, params_fitting, 
             fitting_level= 0.1, occur_length = 4000):
    
    ################ initialize the simulation model ##########################
    _simu_structure_params = deepcopy(simu_structure_params)
    _simu_structure_params['k'] = _k
    simu_once = ms.simu(**_simu_structure_params)
    simu_once.set_generator_seed(_idx)
    simu_once.set_spreading_func('weibull', 'weibull')
    simu_once.set_total_spreading_params(_spreading_params[:2],
                                         _spreading_params[2:])
    simu_once.set_amount_seeds(_amount_seeds)
    ###########################################################################
    
    #################### get the simulation data ##############################
    _simu_once_data = {'s': [], 'i': [], 'r': [], 'w': [], 'c': [], 'd': [], 'y': []}
    _simu_once_arr_data = {'c': []}
    for _target in _simu_once_data.keys():
        _simu_once_data[_target].append(simu_once.get_x_amount(_target) / simu_once.get_node_amount())
    _simu_once_arr_data['c'].append(np.array(simu_once.get_x_amount_arr('c')) / np.array(simu_once.get_population_amounts()))
    while True:
        if simu_once.get_i_amount() == 0:
            break
        simu_once.spread_once()
        for _target in _simu_once_data.keys():
            _simu_once_data[_target].append(simu_once.get_x_amount(_target) / simu_once.get_node_amount())
        _simu_once_arr_data['c'].append(np.array(simu_once.get_x_amount_arr('c')) / np.array(simu_once.get_population_amounts()))
    for _target in _simu_once_data.keys():
        _simu_once_data[_target] = np.array(_simu_once_data[_target])
    _simu_once_arr_data['c'] = np.array(_simu_once_arr_data['c'])
    del simu_once
    ###########################################################################
    
    ##################### get the simulation data for fitting #################
    _fitting_data = {'confirmed': [], 'removal': [], 'confirmed_arr': []}
    _fitting_end_c = (_simu_once_data['c'][-1] - _simu_once_data['c'][0]) * fitting_level + _simu_once_data['c'][0]
    for _i in range(len(_simu_once_data['c'])):
        _fitting_data['confirmed'].append(_simu_once_data['c'][_i])
        _fitting_data['removal'].append(_simu_once_data['r'][_i])
        _fitting_data['confirmed_arr'].append(_simu_once_arr_data['c'][_i])
        if _fitting_data['confirmed'][-1] > _fitting_end_c:
            break
    _fitting_data['confirmed'] = np.array(_fitting_data['confirmed'])
    _fitting_data['removal'] = np.array(_fitting_data['removal'])
    _fitting_data['confirmed_arr'] = np.array(_fitting_data['confirmed_arr'])
    _i_in_data = np.append(_fitting_data['confirmed_arr'][0:1], 
                           np.diff(_fitting_data['confirmed_arr'], axis = 0), axis = 0)
    ###########################################################################
    
    ################## initialize fitting model ###############################
    _calc_params = deepcopy(calc_params)
    _calc_params['real_data'] = _fitting_data
    _calc_params['occur_inf'] = None
    _calc_params['occur_rem'] = None
    _calc_params['k'] = _k
    fitting_once = [sir_model(**_calc_params) for i in range(3)]
    _init_params_cum = [np.ones(2), np.ones(4), np.ones(4)]
    _fitting_func = ['exponent', 'weibull_rate', 'gamma_rate']
    _bounds = [Bounds(np.ones(2) * 0.01, np.ones(2) * np.inf),
               Bounds(np.ones(4) * 0.01, np.ones(4) * np.inf),
               Bounds(np.ones(4) * 0.01, np.ones(4) * np.inf)]
    _param_div = [1, 2, 2]
    _fitting_res = [None, None, None]
    ###########################################################################
    
    ######################## fitting ##########################################
    for i in range(3):
        _init_params = _init_params_cum[i]
        while True:
            _fitting_res[i] =\
            fitting_once[i].fit_cum_from_data(init_params = _init_params, 
                                              occur_length = occur_length,  
                                              fitting_func_inf = _fitting_func[i], 
                                              fitting_func_rem = _fitting_func[i], 
                                              param_div = _param_div[i],
                                              bounds = _bounds[i], **params_fitting)
            if _fitting_res[i].success:
                break
            else:
                _init_params = _init_params_cum[i] * np.random.uniform(0.7, 1.3) 
    ###########################################################################
    
    ###################### return data ###############################        
    _params = {'exponent_inf': _fitting_res[0].x[0], 'exponent_rem': _fitting_res[0].x[1], 
               'weibull_alpha_inf':_fitting_res[1].x[0], 'weibull_beta_inf':1.0 / _fitting_res[1].x[1], 
               'weibull_alpha_rem':_fitting_res[1].x[2], 'weibull_beta_rem':1.0 / _fitting_res[1].x[3], 
               'gamma_alpha_inf':_fitting_res[2].x[0], 'gamma_beta_inf':1.0 / _fitting_res[2].x[1], 
               'gamma_alpha_rem':_fitting_res[2].x[2], 'gamma_beta_rem':1.0 / _fitting_res[2].x[3]}
    return {'params':_params, 'i_in_data': _i_in_data, 'simu_curves': _simu_once_data}
    
def get_calc_data(_params, _i_in_data, _k, _simu_curves, calc_params, occur_length = 4000):
    fitting_once = []
    _srv_func = [func.srv_exponent, func.srv_weibull_scale, func.srv_gamma_scale]
    _inf_params = [(_params['exponent_inf'],),
                  (_params['weibull_alpha_inf'], _params['weibull_beta_inf']),
                  (_params['gamma_alpha_inf'], _params['gamma_beta_inf'])]
    _rem_params = [(_params['exponent_rem'],),
                  (_params['weibull_alpha_rem'], _params['weibull_beta_rem']),
                  (_params['gamma_alpha_rem'], _params['gamma_beta_rem'])]
    _rem_curves = {}
    _nmfuncs = ['exponent', 'weibull', 'gamma']
    for idx, _nmfunc in enumerate(_nmfuncs):
        _i_in = np.append(_simu_curves['c'][0], np.diff(_simu_curves['c']))
        _data_len = len(_i_in)
        _rem_srv = np.array([1 - _srv_func[idx](j * calc_params['step'], *_rem_params[idx]) for j in range(_data_len)])
        _rem_curves[_nmfunc] = np.convolve(_rem_srv, _i_in, 'full')[:_data_len] + _simu_curves['r'][0]
    for i in range(3):
        _calc_params = deepcopy(calc_params)
        vgnr_inf = gnr(_srv_func[i], 'func', _inf_params[i])
        vgnr_rem = gnr(_srv_func[i], 'func', _rem_params[i])
        _calc_params['occur_inf'] = occur(vgnr_inf, vgnr_type = 'srv', length = occur_length, step = _calc_params['step'])
        _calc_params['occur_rem'] = occur(vgnr_rem, vgnr_type = 'srv', length = occur_length, step = _calc_params['step'])
        _calc_params['k'] = _k
        fitting_once.append(sir_model(**_calc_params))
        fitting_once[i].set_i_in_init(_i_in_data)
    for j in range(len(_simu_curves['c']) - len(_i_in_data)):
        for i in range(3):
            fitting_once[i].spread_once()
    _calc_curves = {}
    _targets = ['s', 'i', 'r', 'c', 'd', 'y']
    for _target in _targets:
        for i, _nmfunc in enumerate(_nmfuncs):
            _calc_curves[_target +'_'+_nmfunc] = np.array(fitting_once[i].get_x_tot(_target))
    for _nmfunc in _nmfuncs:
        _calc_curves['w_' + _nmfunc] = _calc_curves['r_' + _nmfunc] - _calc_curves['d_' + _nmfunc]
    del fitting_once
    _targets = ['s', 'i', 'r', 'w', 'c', 'd', 'y']
    _res = {}
    _i_in_len = len(_i_in_data)
    for _target in _targets:
        _res[_target + '_tot_square'] = (_simu_curves[_target] ** 2).sum()
        _res[_target + '_tot'] = _simu_curves[_target].sum()
        _res[_target +'_tot_steady'] = _simu_curves[_target][-1]
        for _nmfunc in _nmfuncs:
            _res[_target +'_'+_nmfunc + '_square'] = ((_simu_curves[_target][_i_in_len:] - _calc_curves[_target +'_'+_nmfunc][_i_in_len:]) ** 2).sum()
            _res[_target +'_'+_nmfunc] = abs(_simu_curves[_target][_i_in_len:] - _calc_curves[_target +'_'+_nmfunc][_i_in_len:]).sum()
            _res[_target +'_'+_nmfunc+ '_steady'] = _calc_curves[_target +'_'+_nmfunc][-1]
    _rem_res = {}
    _rem_res['tot_square'] = (_simu_curves['r'] ** 2).sum()
    _rem_res['tot'] = _simu_curves['r'].sum()
    _rem_res['tot_steady'] = _simu_curves['r'][-1]
    for _nmfunc in _nmfuncs:
        _rem_res[_nmfunc + '_square'] = ((_simu_curves['r'] - _rem_curves[_nmfunc]) ** 2).sum()
        _rem_res[_nmfunc] = abs(_simu_curves['r'] - _rem_curves[_nmfunc]).sum()
        _rem_res[_nmfunc + '_steady'] = _rem_curves[_nmfunc][-1]
    return {'res': _res, 'calc_curves':_calc_curves, 'rem_curves': _rem_curves, 'rem_res': _rem_res}

def get_v_data(_idx, _simu_spreading_params, _fitting_spreading_params, spreading_len,
               _i_in_data, _k, _amount_seeds, 
               simu_structure_params, calc_params, 
               daily_vac_amount, vac_times = 7, vac_duration = 24, occur_length = 4000):
    
    def _get_alloc(s_arr, vac_amount, vac_group = None, adjust = True):
        if vac_group is None or len(vac_group) == 0:
            return None
        if adjust:
            vac_alloc = np.zeros(len(s_arr), dtype = np.int64)
        else:
            vac_alloc = np.zeros(len(s_arr))
        if s_arr[vac_group].sum() == 0:
            proportion = np.ones(len(vac_group)) / len(vac_group)
        else:
            proportion = s_arr[vac_group] / s_arr[vac_group].sum()
        if adjust:
            vac_alloc_tmp = (proportion * vac_amount).astype(np.int32)
            for i in np.argsort(proportion * vac_amount - vac_alloc_tmp):
                vac_alloc_tmp[i] += 1
                if vac_alloc_tmp.sum() == vac_amount:
                    break
            vac_alloc[vac_group] = vac_alloc_tmp
        else:
            vac_alloc[vac_group] = proportion * vac_amount
        return vac_alloc
    
    vac_groups = {'no_vaccine': None, 'under_20': np.arange(0, 2).tolist(), 
                  '20-49': np.arange(2, 5).tolist(), 
                  '20+': np.arange(2, 8).tolist(), '60+': np.arange(6, 8).tolist(), 
                  'all_ages': np.arange(0, 8).tolist()}
    vac_dates = [len(_i_in_data) + i * vac_duration for i in range(vac_times)]
    spreading_length = max(vac_dates[-1] + 30 * 24, spreading_len)
    
    _simu_structure_params = deepcopy(simu_structure_params)
    _simu_structure_params['k'] = _k
    simu_once = ms.simu(**_simu_structure_params)
    _node_amount = simu_once.get_node_amount()
    del simu_once
    simu_data = {'no_vaccine': None, 'under_20': None, '20-49': None, '20+': None,'60+': None, 'all_ages': None}
    for vac_group_key in vac_groups.keys():
        simu_data[vac_group_key] = {'s': [], 'i': [], 'r': [], 'w': [], 'c': [], 'd': [], 'y': []}
        _simu_structure_params = deepcopy(simu_structure_params)
        _simu_structure_params['k'] = _k
        simu_once = ms.simu(**_simu_structure_params)
        simu_once.set_generator_seed(_idx)
        simu_once.set_spreading_func('weibull', 'weibull')
        simu_once.set_total_spreading_params([_simu_spreading_params[0], _simu_spreading_params[1]],
                                             [_simu_spreading_params[2], _simu_spreading_params[3]])
        simu_once.set_amount_seeds(_amount_seeds)
        for _target in simu_data[vac_group_key].keys():
            simu_data[vac_group_key][_target].append(simu_once.get_x_amount(_target) / simu_once.get_node_amount())
        for i in range(spreading_length - 1):
            if i in vac_dates:
                vac_alloc_simu = _get_alloc(np.array(simu_once.get_s_amount_arr()), daily_vac_amount, vac_groups[vac_group_key])
                vac_groups_remain = None
                if vac_groups[vac_group_key] is not None:
                    vac_groups_remain = [gr for gr in np.arange(0,8).tolist() if gr not in vac_groups[vac_group_key]]
                if vac_alloc_simu is not None:
                    vac_remain = simu_once.add_vaccine(vac_alloc_simu)
                    if vac_remain > 0:
                        vac_alloc_simu_remain = _get_alloc(np.array(simu_once.get_s_amount_arr()), vac_remain, vac_groups_remain)
                        if vac_alloc_simu_remain is not None:
                            simu_once.add_vaccine(vac_alloc_simu_remain)
            simu_once.spread_once()
            for _target in simu_data[vac_group_key].keys():
                simu_data[vac_group_key][_target].append(simu_once.get_x_amount(_target) / simu_once.get_node_amount())
        for _target in simu_data[vac_group_key].keys():
            simu_data[vac_group_key][_target] = np.array(simu_data[vac_group_key][_target])
        del simu_once
    
    calc_data = {'no_vaccine': {}, 'under_20': {}, '20-49': {}, '20+': {}, '60+': {}, 'all_ages': {}}
    for vac_group_key in vac_groups.keys():
        calc_once = []
        _srv_func = [func.srv_exponent, func.srv_weibull_scale, func.srv_gamma_scale]
        _inf_params = [(_fitting_spreading_params['exponent_inf'],),
                      (_fitting_spreading_params['weibull_alpha_inf'], _fitting_spreading_params['weibull_beta_inf'])]
        _rem_params = [(_fitting_spreading_params['exponent_rem'],),
                      (_fitting_spreading_params['weibull_alpha_rem'], _fitting_spreading_params['weibull_beta_rem'])]
        for i in range(2):
            _calc_params = deepcopy(calc_params)
            vgnr_inf = gnr(_srv_func[i], 'func', _inf_params[i])
            vgnr_rem = gnr(_srv_func[i], 'func', _rem_params[i])
            _calc_params['occur_inf'] = occur(vgnr_inf, vgnr_type = 'srv', length = occur_length, step = _calc_params['step'])
            _calc_params['occur_rem'] = occur(vgnr_rem, vgnr_type = 'srv', length = occur_length, step = _calc_params['step'])
            _calc_params['k'] = _k
            calc_once.append(sir_model(**_calc_params))
            calc_once[i].set_i_in_init(_i_in_data)
            for j in range(len(_i_in_data) - 1, spreading_length - 1):
                if j in vac_dates:
                    vac_alloc_calc = _get_alloc(np.array(calc_once[i].getc_s_eff()), daily_vac_amount / _node_amount, vac_groups[vac_group_key], 
                                                adjust = False)
                    vac_groups_remain = None
                    if vac_groups[vac_group_key] is not None:
                        vac_groups_remain = [gr for gr in np.arange(0,8).tolist() if gr not in vac_groups[vac_group_key]]
                    if vac_alloc_calc is not None:
                        vac_remain = calc_once[i].add_vaccine(vac_alloc_calc / calc_once[i].get_populations())
                        if vac_remain > 0:
                            vac_alloc_calc_remain = _get_alloc(np.array(calc_once[i].getc_s_eff()), vac_remain, vac_groups_remain, adjust = False)
                            if vac_alloc_calc_remain is not None:
                                calc_once[i].add_vaccine(vac_alloc_calc_remain / calc_once[i].get_populations())
                calc_once[i].spread_once()
        for _target in ['s', 'i', 'r', 'c', 'd', 'y']:
            calc_data[vac_group_key]['exponent_' + _target] = calc_once[0].get_x_tot(_target)
            calc_data[vac_group_key]['weibull_' + _target] = calc_once[1].get_x_tot(_target)
        calc_data[vac_group_key]['exponent_w'] = calc_data[vac_group_key]['exponent_r'] - calc_data[vac_group_key]['exponent_d']
        calc_data[vac_group_key]['weibull_w'] = calc_data[vac_group_key]['weibull_r'] - calc_data[vac_group_key]['weibull_d']
        del calc_once
    
    return {'simu_data': simu_data, 'calc_data': calc_data}

def get_optimal_v_data(_idx, _simu_spreading_params, _fitting_spreading_params,
                       _i_in_data, _k, _amount_seeds, 
                       simu_structure_params, calc_params, 
                       daily_vac_amount, vac_times = 14, vac_duration = 24, consequent_step = 24 * 30, occur_length = 4000):
    
    _daily_vac_p = daily_vac_amount / simu_structure_params['node_amount']
    _simu_structure_params = deepcopy(simu_structure_params)
    _simu_structure_params['k'] = _k
    simu_once = ms.simu(**_simu_structure_params)
    _population_amounts = np.array(simu_once.get_population_amounts())
    del simu_once
    targets = ['c', 'd', 'y']
    simu_data = {}
    _srv_func = [func.srv_exponent, func.srv_weibull_scale, func.srv_gamma_scale]
    _inf_params = [(_fitting_spreading_params['exponent_inf'],),
                  (_fitting_spreading_params['weibull_alpha_inf'], _fitting_spreading_params['weibull_beta_inf']),
                  (_fitting_spreading_params['gamma_alpha_inf'], _fitting_spreading_params['gamma_beta_inf'])]
    _rem_params = [(_fitting_spreading_params['exponent_rem'],),
                  (_fitting_spreading_params['weibull_alpha_rem'], _fitting_spreading_params['weibull_beta_rem']),
                  (_fitting_spreading_params['gamma_alpha_rem'], _fitting_spreading_params['gamma_beta_rem'])]
    for i, way in enumerate(['exponent', 'weibull', 'gamma']):
        simu_data[way] = {}
        _calc_params = deepcopy(calc_params)
        vgnr_inf = gnr(_srv_func[i], 'func', _inf_params[i])
        vgnr_rem = gnr(_srv_func[i], 'func', _rem_params[i])
        _calc_params['occur_inf'] = occur(vgnr_inf, vgnr_type = 'srv', length = occur_length, step = _calc_params['step'])
        _calc_params['occur_rem'] = occur(vgnr_rem, vgnr_type = 'srv', length = occur_length, step = _calc_params['step'])
        _calc_params['k'] = _k
        calc_once = sir_model(**_calc_params)
        calc_once.set_i_in_init(_i_in_data)
        calc_once.spread_once()
        for target in targets:
            group_amount = len(calc_once.get_populations())
            bound = Bounds(np.zeros(group_amount), calc_once.getc_s_eff())
            linearconstraint = LinearConstraint(calc_once.get_populations(), 
                                                np.zeros(1), np.array([_daily_vac_p]))
            init_strategy = calc_once.getc_s_eff() * _daily_vac_p \
            / (calc_once.getc_s_eff() @ calc_once.get_populations())
            while True:
                _vac_res = minimize(calc_once.get_vac_x, init_strategy, 
                                    args = (vac_times, vac_duration, consequent_step, target), method='SLSQP',
                                    jac = None, options={'ftol': 1e-16, 'disp': True, 'maxiter':100000},
                                    constraints = linearconstraint, bounds = bound)
                if _vac_res.success == True:
                    break
                init_strategy = calc_once.getc_s_eff() * _daily_vac_p \
                / (calc_once.getc_s_eff() @ calc_once.get_populations()) * np.random.uniform(0.5, 1)
            _vac_amount_alloc = (_vac_res.x * _population_amounts).astype(np.int64)
            for j in np.argsort(_vac_res.x * _population_amounts - _vac_amount_alloc):
                if _vac_amount_alloc.sum() == daily_vac_amount:
                    break
                _vac_amount_alloc[j] += 1
            simu_once = ms.simu(**_simu_structure_params)
            simu_once.set_generator_seed(_idx)
            simu_once.set_total_spreading_params(_simu_spreading_params[0], _simu_spreading_params[1],
                                                 _simu_spreading_params[2], _simu_spreading_params[3])
            simu_once.set_amount_seeds(_amount_seeds)
            for j in range(len(_i_in_data)):
                simu_once.spread_once()
            for j in range(vac_times):
                simu_once.add_vaccine(_vac_amount_alloc)
                for _j in range(vac_duration):
                    simu_once.spread_once()
            simu_data[way][target] = [simu_once.get_x_amount(target) / simu_once.get_node_amount()]
            simu_once.spread_to_end()
            simu_data[way][target].append(simu_once.get_x_amount(target) / simu_once.get_node_amount())
            del simu_once
    return simu_data
###############################################################################












        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    
    
    
    