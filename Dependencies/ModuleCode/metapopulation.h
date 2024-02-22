#pragma once
#include "model_elem.h"
#include "spreading_func.h"


class metapopulation
{
private:
	size_t const node_amount;
	size_t const group_amount;
	std::vector<double> populations;
	std::vector<std::vector<double>> contacts;
	std::vector<double> ifrs;
	std::vector<double> ylls;
	std::vector<size_t> population_amounts;
	double k_value;
	size_t delay;
	double eta;
	double const step;
	double time;
	std::vector<state_groups*> groups;
	size_t s_amount;
	std::vector<size_t> s_population_amounts;
	size_t i_amount;
	std::vector<size_t> i_population_amounts;
	size_t r_amount;
	std::vector<size_t> r_population_amounts;
	size_t w_amount;
	std::vector<size_t> w_population_amounts;
	size_t c_amount;
	std::vector<size_t> c_population_amounts;
	size_t d_amount;
	std::vector<size_t> d_population_amounts;
	double y_amount;
	std::vector<double> y_population_amounts;
	size_t v_amount;
	std::vector<size_t> v_population_amounts;
	size_t p_amount;
	std::vector<size_t> p_population_amounts;
	std::default_random_engine generator;
	probability prob;
	hazard_type* hazard_inf;
	dist_type* dist_rem;
public:
	metapopulation(size_t _node_amount, size_t _group_amount, std::vector<double>& _populations, std::vector<std::vector<double>>& _contacts, 
				   std::vector<double>& _ifrs, std::vector<double>& _ylls, double _k_value, double _step, double _time);
	~metapopulation();
	void set_vaccination_params(size_t _delay, double _eta);
	void set_spreading_func(hazard_type* _hazard_inf, dist_type* _dist_rem);
	void set_spreading_params(std::vector<double> const& _params_inf, std::vector<double> const& _params_rem);
	void set_spreading_params(std::vector<std::vector <std::vector<double>>> const& _params_inf, std::vector<std::vector <std::vector<double>>> const& _params_rem);
	void set_spreading_survivals(std::vector<double> const& _survivals_inf, std::vector<double> const& _survivals_rem);
	void init_all_states();
	void set_seeds(std::vector<size_t>& _init_amounts, double _time);
	void set_seeds(std::vector<double>& _init_probs, double _time);
	void set_seeds(std::vector<std::vector<size_t>>& _init_states, double _time);
	void adjust_rem_time();
	void reset_trans_time();
	void reset_rem_time();
	void set_generator_seed(int _seed);
	void spread_once();
	size_t add_vaccine(std::vector<size_t>& _allocation);
	size_t get_node_amount() const;
	size_t get_group_amount() const;
	double get_time() const;
	size_t get_s_amount() const;
	size_t get_i_amount() const;
	size_t get_r_amount() const;
	size_t get_w_amount() const;
	size_t get_c_amount() const;
	size_t get_d_amount() const;
	double get_y_amount() const;
	size_t get_v_amount() const;
	size_t get_p_amount() const;
	std::vector<size_t> get_s_amount_arr() const;
	std::vector<size_t> get_i_amount_arr() const;
	std::vector<size_t> get_r_amount_arr() const;
	std::vector<size_t> get_w_amount_arr() const;
	std::vector<size_t> get_c_amount_arr() const;
	std::vector<size_t> get_d_amount_arr() const;
	std::vector<double> get_y_amount_arr() const;
	std::vector<size_t> get_v_amount_arr() const;
	std::vector<size_t> get_p_amount_arr() const;
	size_t get_s_amount_idx(size_t _idx) const;
	size_t get_i_amount_idx(size_t _idx) const;
	size_t get_r_amount_idx(size_t _idx) const;
	size_t get_w_amount_idx(size_t _idx) const;
	size_t get_c_amount_idx(size_t _idx) const;
	size_t get_d_amount_idx(size_t _idx) const;
	double get_y_amount_idx(size_t _idx) const;
	size_t get_v_amount_idx(size_t _idx) const;
	size_t get_p_amount_idx(size_t _idx) const;
	double get_node_trans_time(size_t _group_idx, size_t _node_idx) const;
	double get_node_rem_time(size_t _group_idx, size_t _node_idx) const;
	char get_node_state(size_t _group_idx, size_t _node_idx) const;
	std::pair<std::vector<double>, std::vector<double>> const& get_node_spreading_params(size_t _group_idx, size_t _node_idx);
	size_t get_delay() const;
	double get_eta() const;
	double get_step() const;
	std::vector<double> const& get_populations() const;
	std::vector<std::vector<double>> const& get_contacts() const;
	std::vector<double> const& get_ifrs() const;
	std::vector<size_t> const& get_population_amounts() const;
};