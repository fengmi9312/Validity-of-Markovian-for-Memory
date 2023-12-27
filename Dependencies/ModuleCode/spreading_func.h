#pragma once
#include<vector>
#include<random>

typedef double (*hazard_func)(std::vector<double> const& _params_rem, double _tau, double _actl_step);

class dist_type
{
protected:
	std::vector<double> params;
	std::vector<double> survivals;
public:
	void param(std::vector<double> const& _params);
	void set_survivals(std::vector<double> const& _survivals);
	virtual double operator()(std::default_random_engine& _generator) = 0;
};

class weibull_dist : public dist_type
{
private:
	std::weibull_distribution<double> dist;
public:
	double operator()(std::default_random_engine& _generator);
};

class xburr_dist : public dist_type
{
private:
	std::uniform_real_distribution<double> dist;
public:
	xburr_dist();
	double operator()(std::default_random_engine& _generator);
};

class xlog_logistic_dist : public dist_type
{
private:
	std::uniform_real_distribution<double> dist;
public:
	xlog_logistic_dist();
	double operator()(std::default_random_engine& _generator);
};

class general_dist : public dist_type
{
private:
	std::uniform_real_distribution<double> dist;
public:
	general_dist();
	double operator()(std::default_random_engine& _generator);
};


class hazard_type
{
protected:
	std::vector<double> params;
	std::vector<double> survivals;
public:
	void param(std::vector<double> const& _params);
	void set_survivals(std::vector<double> const& _survivals);
	virtual double operator()(double _tau, double _actl_step) = 0;
};

class weibull_hazard : public hazard_type
{
public:
	double operator()(double _tau, double _actl_step);
};

class xburr_hazard : public hazard_type
{
public:
	double operator()(double _tau, double _actl_step);
};

class xlog_logistic_hazard : public hazard_type
{
public:
	double operator()(double _tau, double _actl_step);
};

class general_hazard : public hazard_type
{
public:
	double operator()(double _tau, double _actl_step);
};