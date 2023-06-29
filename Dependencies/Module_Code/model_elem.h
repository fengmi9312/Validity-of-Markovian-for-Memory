#pragma once
#include <map>
#include <list>
#include <random>

std::list<size_t>& diff_rand(size_t num, size_t _a, size_t _b, std::default_random_engine& _generator_ref);

class probability
{
private:
	std::uniform_real_distribution<double> distribution;
public:
	probability();
	bool operator()(double prob, std::default_random_engine& _generator_ref);
};

class node
{
private:
	char state;
	double trans_time;
	double rem_time;
	int remaining_ptime;
	std::vector<double> params_inf;
	std::vector<double> params_rem;
	void set_state(char _state, double _trans_time);
	void set_rem_time(double _rem_time);
	void set_remaining_ptime(int _remaining_ptime);
	void remaining_ptime_dec();
	void invalidate_rem_time();
	void invalidate_remaining_ptime();
	void set_spreading_params(std::vector<double> const& _params_inf, std::vector<double> const& _params_rem);
	void init();
	friend class metapopulation;
public:
	node();
	char get_state() const;
	double get_trans_time() const;
	double get_rem_time() const;
	int get_remaining_ptime() const;
	std::vector<double> get_params_inf() const;
	std::vector<double> get_params_rem() const;
};


class state_groups
{
private:
	std::map <char, std::list<node*>> state_group_map;
	std::vector<node*> all_nodes;
	state_groups(size_t _node_amount);
	friend class metapopulation;
public:
	~state_groups();
	size_t get_node_amount() const;
	size_t get_state_node_amount(char _state) const;
};


