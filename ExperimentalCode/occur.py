# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 21:01:31 2020

@author: Tingting
"""

import numpy as np
from math import log

class gnr:
    
    def __init__(self, x, x_type, args = ()):
        self.__x = x
        self.__x_type = x_type
        self.__args = args

    def get_value(self, idx, step):
        if self.__x_type == 'arr':
            if idx >= len(self.__x):
                return 0
            else:
                return self.__x[idx]
        else:
            return self.__x(idx * step, *self.__args)


class occur:
    
    def __init__(self, vgnr, vgnr_type = 'dist' , length = 1000, step = 0.001, start = 0):
        if vgnr_type != 'srv' and vgnr_type != 'rate' and vgnr_type != 'dist':
            return None
        self.__values = {'srv': None, 'rate': None, 'dist': None}
        self.__step = step
        self.__vgnr = vgnr
        self.__vgnr_type = vgnr_type
        self.__start = start
        self.__len = length
        __length = self.__len
        if vgnr_type == 'srv':
            __length += 1
        __vals = []
        for i in range(__length):
            __val = vgnr.get_value(i + self.__start, self.__step)
            if __val >= 0 and __val <= 1:
                if vgnr_type == 'srv' and len(__vals) > 0:
                        if __val <= __vals[-1]:
                            __vals.append(__val)
                        else:
                            print('illegal values')
                            return None
                else:
                    __vals.append(__val)
            else:
                print('illegal values')
                return None
        self.__values[self.__vgnr_type] = np.array(__vals)
        if vgnr_type == 'srv':
            self.__values['rate'] = self.__from_srv_to_rate()
            self.__values['dist'] = self.__from_srv_to_dist()
        elif vgnr_type == 'rate':
            self.__values['srv'] = self.__from_rate_to_srv()
            self.__values['dist'] = self.__from_srv_to_dist()
        elif vgnr_type == 'dist':
            self.__values['srv'] = self.__from_dist_to_srv()
            self.__values['rate'] = self.__from_srv_to_rate()
        else:
            pass
        nonzeroidx = np.argwhere(self.__values['rate'] != 0)
        if len(nonzeroidx) == 0:
            self.__actl_len = 0
        else:
            self.__actl_len = nonzeroidx[-1][0] + 1
        
    def __from_rate_to_srv(self):
        __val_list = [1.0,]
        for __val in self.__values['rate']:
            __val_list.append(__val_list[-1] * (1 - __val))
        return np.array(__val_list)
        
    def __from_dist_to_srv(self):
        __val_list = [1.0,]
        for __val in self.__values['dist']:
            __val_list.append(__val_list[-1] - __val)
        return np.array(__val_list)
          
    def __from_srv_to_dist(self):
        __val_list = []
        for __val0, __val1 in zip(self.__values['srv'][:-1], self.__values['srv'][1:]):
            __val_list.append(__val0 - __val1)
        return np.array(__val_list)
        
    def __from_srv_to_rate(self):
        __val_list = []
        for __val0, __val1 in zip(self.__values['srv'][:-1], self.__values['srv'][1:]):
            if __val0 == 0:
                __val_list.append(0)
            else:
                __val_list.append(1.0 - __val1 / __val0)
        return np.array(__val_list)
    
    def get_srv(self):
        return self.__values['srv']
    
    def get_rate(self):
        return self.__values['rate']
    
    def get_dist(self):
        return self.__values['dist']
    
    def get_len(self):
        return self.__len
    
    def get_actl_len(self):
        return self.__actl_len
    
    def get_step(self):
        return self.__step
    
    def get_vgnr(self):
        return self.__vgnr
    
    def get_vgnr_type(self):
        return self.__vgnr_type
    
    def get_tline(self):
        return np.arange(self.__len) * self.__step
    
    def get_p(self):
        return self.__values['srv'][0] - self.__values['srv'][-1]
    
    def get_lambda(self):
        __p = self.get_p()
        if __p == 0:
            return np.inf
        else:
            return log(1.0 / (1.0 - self.get_p()))
        
    def cut_length(self, length):
        self.__values['srv'] = self.__values['srv'][:length + 1]
        self.__values['rate'] = self.__values['rate']
        self.__values['dist'] = self.__values['dist']
        self.__len = length
        nonzeroidx = np.argwhere(self.__values['rate'] != 0)
        if len(nonzeroidx) == 0:
            self.__actl_len = 0
        else:
            self.__actl_len = nonzeroidx[-1][0] + 1
        
    
def occur_mul(a, b, time_int = 0):
    if a.get_step() != b.get_step():
        print('The steps of a and b must be equal.')
        return None
    __length = max(a.get_actl_len(), b.get_actl_len())
    x = occur(a.get_vgnr(), a.get_vgnr_type(), __length, a.get_step(), time_int)
    y = occur(b.get_vgnr(), b.get_vgnr_type(), __length, b.get_step(), time_int)
    return occur(gnr(x.get_rate() * y.get_srv()[:-1], 'arr'), vgnr_type = 'rate', 
                 length = __length, step = a.get_step())

def occur_mul_x(a, b, time_int = 0):
    if a.get_step() != b.get_step():
        print('The steps of a and b must be equal.')
        return None
    __length = max(a.get_actl_len(), b.get_actl_len())
    x = occur(a.get_vgnr(), a.get_vgnr_type(), __length, a.get_step(), time_int)
    y = occur(b.get_vgnr(), b.get_vgnr_type(), __length, b.get_step(), time_int)
    return occur(gnr(x.get_dist() * y.get_srv()[:-1], 'arr'), vgnr_type = 'rate', 
                 length = __length, step = a.get_step())
        

def calc_lambda(a, b):
    if a.get_step() != b.get_step():
        print('The steps of a and b must be equal.')
        return None
    __length = max(a.get_actl_len(), b.get_actl_len()) * 2
    x = occur(a.get_vgnr(), a.get_vgnr_type(), __length, a.get_step())
    y = occur(b.get_vgnr(), b.get_vgnr_type(), __length, b.get_step())
    return np.log(1.0 / (1 - (x.get_rate() * y.get_srv()[:-1])).prod())


def calc_lambda_arr(a, b, time_int = 0):
    if a.get_step() != b.get_step():
        print('The steps of a and b must be equal.')
        return None
    __length = max(a.get_actl_len(), b.get_actl_len()) * 2
    x = occur(a.get_vgnr(), a.get_vgnr_type(), __length, a.get_step())
    y = occur(b.get_vgnr(), b.get_vgnr_type(), __length, b.get_step())
    res = []
    for i in range(min(x.get_actl_len(), y.get_actl_len())):
        __i_srv = y.get_srv()[i]
        res.append((1 - (x.get_rate()[i:] * y.get_srv()[i:-1] / __i_srv)).prod())
    res = np.array(res)
    res[res == 0] = res[res != 0].min()
    return np.log(1.0 / res)
    
        
        
        
        
        