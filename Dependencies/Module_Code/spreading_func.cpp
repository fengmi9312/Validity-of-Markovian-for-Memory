#include "spreading_func.h"
#include <stdexcept>

void dist_type::param(std::vector<double> const& _params)
{
	params = _params;
}

void dist_type::set_survivals(std::vector<double> const& _survivals)
{
	if (_survivals.size() <= 1) throw std::invalid_argument("the length of arguments must be larger than 1");
	std::vector<double>::const_iterator _r_it = _survivals.begin();
	++_r_it;
	std::vector<double>::const_iterator _l_it = _survivals.begin();
	while (_r_it != _survivals.end())
	{
		if ((*_r_it) > (*_l_it)) throw std::invalid_argument("the survivals must decrease over time");
		++_r_it;
		++_l_it;
	}
	survivals = _survivals;
}

double weibull_dist::operator()(std::default_random_engine& _generator)
{
	if (params.size() != 2) throw std::invalid_argument("the length of arguments must be 2");
	dist.param(std::weibull_distribution<double>::param_type(params[0], params[1]));
	return dist(_generator);
}

xburr_dist::xburr_dist() :dist(0, 1) {}

double xburr_dist::operator()(std::default_random_engine& _generator)
{
	if (params.size() != 3) throw std::invalid_argument("the length of arguments must be 3");
	return std::pow(std::pow(1.0 - dist(_generator), -params[1] / params[2]) - 1, 1.0 / params[0]);
}

xlog_logistic_dist::xlog_logistic_dist() :dist(0, 1) {}

double xlog_logistic_dist::operator()(std::default_random_engine& _generator)
{
	if (params.size() != 3) throw std::invalid_argument("the length of arguments must be 3");
	return params[0] * std::pow(std::pow(1.0 - dist(_generator), -params[2]) - 1.0, -1.0 / params[1]);
}

general_dist::general_dist() :dist(0, 1) {}


double general_dist::operator()(std::default_random_engine& _generator)
{
	if (params.size() != 1) throw std::invalid_argument("the length of arguments must be 3");
	double _step = (*(params.begin()));
	double _prob = dist(_generator);
	size_t _arr_head = 0;
	size_t _arr_tail = survivals.size() - 1;
	if (_prob > survivals[_arr_head]) return 0;
	else if (_prob < survivals[_arr_tail]) return survivals.size() * _step;
	else
	{
		size_t _arr_mid;
		while (_arr_head != _arr_tail - 1)
		{
			_arr_mid = (_arr_head + _arr_tail) / 2;
			if (_prob < survivals[_arr_mid]) _arr_head = _arr_mid;
			else _arr_tail = _arr_mid;
		}
		return ((_prob - survivals[_arr_tail]) / (survivals[_arr_tail] - survivals[_arr_head]) + _arr_tail) * _step;
	}
}

void hazard_type::param(std::vector<double> const& _params)
{
	params = _params;
}

void hazard_type::set_survivals(std::vector<double> const& _survivals)
{
	if (_survivals.size() <= 1) throw std::invalid_argument("the length of arguments must be larger than 1");
	std::vector<double>::const_iterator _r_it = _survivals.begin();
	++_r_it;
	std::vector<double>::const_iterator _l_it = _survivals.begin();
	while (_r_it != _survivals.end())
	{
		if ((*_r_it) > (*_l_it)) throw std::invalid_argument("the survivals must decrease over time");
		++_r_it;
		++_l_it;
	}
	survivals = _survivals;
}

double weibull_hazard::operator()(double _tau, double _actl_step)
{
	if (params.size() != 2) throw std::invalid_argument("the length of arguments must be 2");
	if (_actl_step == 0) return 0;
	return 1.0 - std::exp(std::pow(_tau / params[1], params[0]) - std::pow((_tau + _actl_step) / params[1], params[0]));
}

double xburr_hazard::operator()(double _tau, double _actl_step)
{
	if (params.size() != 3) throw std::invalid_argument("the length of arguments must be 3");
	if (_actl_step == 0) return 0;
	return 1.0 - std::pow((1.0 + std::pow(_tau + _actl_step, params[0])) / (1.0 + std::pow(_tau, params[0])), -params[1] * params[2]);
}

double xlog_logistic_hazard::operator()(double _tau, double _actl_step)
{
	if (params.size() != 3) throw std::invalid_argument("the length of arguments must be 3");
	if (_actl_step == 0) return 0;
	return 1.0 - std::pow((1 + std::pow(_tau / params[0], params[1])) / (1 + std::pow((_tau + _actl_step) / params[0], params[1])), params[2]);
}



double general_hazard::operator()(double _tau, double _actl_step)
{
	if (_tau < 0 || _actl_step < 0) throw std::invalid_argument("the parameters must be larger than 0");
	size_t _idx0 = size_t(_tau / _actl_step);
	size_t _idx1 = _idx0 + 1;
	double _survival0;
	double _survival1;
	if (_idx0 >= survivals.size() - 1) _survival0 = survivals[survivals.size() - 1];
	else
	{
		_survival0 = survivals[_idx0] - (survivals[_idx0] - survivals[_idx0 + 1]) * (_tau - _idx0 * _actl_step) / _actl_step;
	}
	if (_idx1 >= survivals.size() - 1) _survival1 = survivals[survivals.size() - 1];
	else
	{
		_survival1 = survivals[_idx1] - (survivals[_idx1] - survivals[_idx1 + 1]) * (_tau + _actl_step - _idx1 * _actl_step) / _actl_step;
	}
	if (_survival0 == 0) return 0;
	else return (_survival0 - _survival1) / _survival0;
}