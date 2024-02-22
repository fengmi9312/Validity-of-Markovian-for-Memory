#include"metapopulation.h"

metapopulation::metapopulation(size_t _node_amount, size_t _group_amount, std::vector<double>& _populations, std::vector<std::vector<double>>& _contacts,
	std::vector<double>& _ifrs, std::vector<double>& _ylls, double _k_value, double _step, double _time) :
	node_amount(_node_amount), group_amount(_group_amount), populations(_populations), contacts(_contacts), ifrs(_ifrs), ylls(_ylls), population_amounts(group_amount, 0),
	k_value(_k_value), delay(14), eta(1), step(_step), time(_time)
{
	double _population_sum = 0;
	size_t _remaining_node_amount = node_amount;
	for (std::vector<double>::const_iterator _it = populations.begin(); _it != populations.end(); ++_it)
	{
		_population_sum = _population_sum + (*_it);
	}
	for (size_t i = 0; i < group_amount; ++i)
	{
		populations[i] /= _population_sum;
		population_amounts[i] = (size_t)(populations[i] * node_amount);
		_remaining_node_amount -= population_amounts[i];
	}
	for (size_t i = 0; i < group_amount; ++i)
	{
		if (_remaining_node_amount == 0) break;
		--_remaining_node_amount;
		population_amounts[i] = population_amounts[i] + 1;
	}
	for (size_t i = 0; i < group_amount; ++i)
	{
		groups.push_back(new state_groups(population_amounts[i]));
	}
	s_amount = node_amount;
	s_population_amounts = population_amounts;
	i_amount = 0;
	i_population_amounts = std::vector<size_t>(group_amount, 0);
	r_amount = 0;
	r_population_amounts = std::vector<size_t>(group_amount, 0);
	w_amount = 0;
	w_population_amounts = std::vector<size_t>(group_amount, 0);
	c_amount = 0;
	c_population_amounts = std::vector<size_t>(group_amount, 0);
	d_amount = 0;
	d_population_amounts = std::vector<size_t>(group_amount, 0);
	y_amount = 0;
	y_population_amounts = std::vector<double>(group_amount, 0);
	v_amount = 0;
	v_population_amounts = std::vector<size_t>(group_amount, 0);
	p_amount = 0;
	p_population_amounts = std::vector<size_t>(group_amount, 0);
	hazard_inf = nullptr;
	dist_rem = nullptr;
}


metapopulation::~metapopulation()
{
	for (std::vector<state_groups*>::iterator _it = groups.begin(); _it != groups.end(); ++_it)
		delete* _it;
}


void metapopulation::set_vaccination_params(size_t _delay, double _eta)
{
	delay = _delay;
	eta = _eta;
}

void metapopulation::set_spreading_func(hazard_type* _hazard_inf, dist_type* _dist_rem)
{
	hazard_inf = _hazard_inf;
	dist_rem = _dist_rem;
}

void metapopulation::set_spreading_params(std::vector<double> const& _params_inf, std::vector<double> const& _params_rem)
{
	std::vector<char> keys = { 's', 'i', 'w', 'd', 'v', 'p' };
	for (std::vector<state_groups*>::iterator _it = groups.begin(); _it != groups.end(); ++_it)
	{
		for (std::vector<char>::const_iterator _key_it = keys.begin(); _key_it != keys.end(); ++_key_it)
		{
			for (std::list<node*>::iterator _it1 = (*_it)->state_group_map[*_key_it].begin(); _it1 != (*_it)->state_group_map[*_key_it].end(); ++_it1)
			{
				(*_it1)->set_spreading_params(_params_inf, _params_rem);
			}
		}
	}
}

void metapopulation::set_spreading_params(std::vector<std::vector <std::vector<double>>> const& _params_inf, std::vector<std::vector <std::vector<double>>> const& _params_rem)
{
	std::vector<state_groups*>::iterator _it_node_group = groups.begin();
	std::vector<std::vector <std::vector<double>>>::const_iterator _it_params_inf_group = _params_inf.begin();
	std::vector<std::vector <std::vector<double>>>::const_iterator _it_params_rem_group = _params_rem.begin();
	while (_it_node_group != groups.end() && _it_params_inf_group != _params_inf.end() && _it_params_rem_group != _params_rem.end())
	{
		std::vector<node*>::iterator _it_node = (*_it_node_group)->all_nodes.begin();
		std::vector <std::vector<double>>::const_iterator _it_params_inf = (*_it_params_inf_group).begin();
		std::vector <std::vector<double>>::const_iterator _it_params_rem = (*_it_params_rem_group).begin();
		while (_it_node != (*_it_node_group)->all_nodes.end() && _it_params_inf != (*_it_params_inf_group).end() && _it_params_rem != (*_it_params_rem_group).end())
		{
			(*_it_node)->set_spreading_params(*_it_params_inf, *_it_params_rem);
			++_it_node;
			++_it_params_inf;
			++_it_params_rem;
		}
		++_it_node_group;
		++_it_params_inf_group;
		++_it_params_rem_group;
	}
}

void metapopulation::set_spreading_survivals(std::vector<double> const& _params_inf, std::vector<double> const& _params_rem)
{
	hazard_inf->set_survivals(_params_inf);
	dist_rem->set_survivals(_params_rem);
}


void metapopulation::init_all_states()
{
	node* _node_ptr = nullptr;
	std::vector<char> keys = { 's', 'i', 'w', 'd', 'v', 'p' };
	std::list<node*>::iterator _current_node_it;
	std::list<node*>::iterator _next_node_it;
	for (std::vector<state_groups*>::iterator _it = groups.begin(); _it != groups.end(); ++_it)
	{
		for (std::vector<char>::const_iterator _key_it = keys.begin(); _key_it != keys.end(); ++_key_it)
		{
			_next_node_it = (*_it)->state_group_map[*_key_it].begin();
			while(_next_node_it != (*_it)->state_group_map[*_key_it].end())
			{
				_current_node_it = _next_node_it;
				++_next_node_it;
				(*_current_node_it)->init();
				if (*_key_it != 's')
				{
					_node_ptr = *_current_node_it;
					(*_it)->state_group_map[*_key_it].erase(_current_node_it);
					(*_it)->state_group_map[_node_ptr->get_state()].push_back(_node_ptr);
				}
			}
		}
	}
	s_amount = node_amount;
	s_population_amounts = population_amounts;
	i_amount = 0;
	i_population_amounts = std::vector<size_t>(group_amount, 0);
	r_amount = 0;
	r_population_amounts = std::vector<size_t>(group_amount, 0);
	w_amount = 0;
	w_population_amounts = std::vector<size_t>(group_amount, 0);
	c_amount = 0;
	c_population_amounts = std::vector<size_t>(group_amount, 0);
	d_amount = 0;
	d_population_amounts = std::vector<size_t>(group_amount, 0);
	y_amount = 0;
	y_population_amounts = std::vector<double>(group_amount, 0);
	v_amount = 0;
	v_population_amounts = std::vector<size_t>(group_amount, 0);
	p_amount = 0;
	p_population_amounts = std::vector<size_t>(group_amount, 0);
	time = 0;
}

void metapopulation::set_seeds(std::vector<std::vector<size_t>>& _init_states, double _time)
{
	init_all_states();
	time = _time;
	node* _node_ptr = nullptr;
	for (size_t i = 0; i < group_amount; ++i)
	{
		std::list<node*>::iterator _current_node_it;
		std::list<node*>::iterator _next_node_it = groups[i]->state_group_map['s'].begin();
		std::vector<size_t>::const_iterator _seed_it = _init_states[i].begin();
		while(_next_node_it != groups[i]->state_group_map['s'].end() && _seed_it != _init_states[i].end())
		{
			_current_node_it = _next_node_it;
			_node_ptr = *_current_node_it;
			++_next_node_it;
			if (*_seed_it == 1)
			{
				_node_ptr->set_state('i', time);
				dist_rem->param(_node_ptr->get_params_rem());
				_node_ptr->set_rem_time(time + (*dist_rem)(generator));
				_node_ptr->invalidate_remaining_ptime();
				groups[i]->state_group_map['s'].erase(_current_node_it);
				groups[i]->state_group_map['i'].push_back(_node_ptr);
				--s_amount;
				--s_population_amounts[i];
				++i_amount;
				++i_population_amounts[i];
				++c_amount;
				++c_population_amounts[i];
			}
			else if (*_seed_it == 2 || *_seed_it == 3)
			{
				if (*_seed_it == 2)
				{
					_node_ptr->set_state('w', time);
					_node_ptr->invalidate_rem_time();
					groups[i]->state_group_map['s'].erase(_current_node_it);
					groups[i]->state_group_map['w'].push_back(_node_ptr);
					++w_amount;
					++w_population_amounts[i];
				}
				else
				{
					_node_ptr->set_state('d', time);
					_node_ptr->invalidate_rem_time();
					groups[i]->state_group_map['s'].erase(_current_node_it);
					groups[i]->state_group_map['d'].push_back(_node_ptr);
					++d_amount;
					++d_population_amounts[i];
					y_amount = y_amount + ylls[i];
					y_population_amounts[i] = y_population_amounts[i] + ylls[i];
				}
				_node_ptr->invalidate_remaining_ptime();
				--s_amount;
				--s_population_amounts[i];
				++r_amount;
				++r_population_amounts[i];
				++c_amount;
				++c_population_amounts[i];
			}
			++_seed_it;
		}
	}
}

void metapopulation::adjust_rem_time()
{
	double _tau_rem;
	std::vector<char> keys = { 's', 'i', 'w', 'd', 'v', 'p' };
	for (std::vector<state_groups*>::iterator _group_it = groups.begin(); _group_it != groups.end(); ++_group_it)
	{
		for (std::vector<char>::const_iterator _key_it = keys.begin(); _key_it != keys.end(); ++_key_it)
		{
			for (std::list<node*>::iterator _node_it = (*_group_it)->state_group_map[*_key_it].begin(); _node_it != (*_group_it)->state_group_map[*_key_it].end(); ++_node_it)
			{
				if ((*_node_it)->get_state() == 'i')
				{
					do
					{
						dist_rem->param((*_node_it)->get_params_rem());
						_tau_rem = (*dist_rem)(generator);
					} while (_tau_rem <= time - (*_node_it)->get_trans_time());
					(*_node_it)->set_rem_time(_tau_rem + (*_node_it)->get_trans_time());
				}
			}
		}
	}
}

void metapopulation::reset_trans_time()
{
	std::vector<char> keys = { 's', 'i', 'w', 'd', 'v', 'p' };
	for (std::vector<state_groups*>::iterator _group_it = groups.begin(); _group_it != groups.end(); ++_group_it)
	{
		for (std::vector<char>::const_iterator _key_it = keys.begin(); _key_it != keys.end(); ++_key_it)
		{
			for (std::list<node*>::iterator _node_it = (*_group_it)->state_group_map[*_key_it].begin(); _node_it != (*_group_it)->state_group_map[*_key_it].end(); ++_node_it)
			{
				if ((*_node_it)->get_state() == 'i')
				{
					(*_node_it)->set_state('i', time);
				}
			}
		}
	}
}

void metapopulation::reset_rem_time()
{
	double _tau_rem;
	std::vector<char> keys = { 's', 'i', 'w', 'd', 'v', 'p' };
	for (std::vector<state_groups*>::iterator _group_it = groups.begin(); _group_it != groups.end(); ++_group_it)
	{
		for (std::vector<char>::const_iterator _key_it = keys.begin(); _key_it != keys.end(); ++_key_it)
		{
			for (std::list<node*>::iterator _node_it = (*_group_it)->state_group_map[*_key_it].begin(); _node_it != (*_group_it)->state_group_map[*_key_it].end(); ++_node_it)
			{
				if ((*_node_it)->get_state() == 'i')
				{
					dist_rem->param((*_node_it)->get_params_rem());
					_tau_rem = (*dist_rem)(generator);
					(*_node_it)->set_rem_time(_tau_rem + time);
				}
			}
		}
	}
}

void metapopulation::set_seeds(std::vector<size_t>& _init_amounts, double _time)
{
	std::vector<std::vector<size_t>> _init_states(group_amount);
	for (size_t i = 0; i < group_amount; ++i)
	{
		std::list<size_t> _inf_num_list = diff_rand(_init_amounts[i], 0, population_amounts[i], generator);
		_init_states[i] = std::vector<size_t>(population_amounts[i], 0);
		for (std::list<size_t>::const_iterator _it = _inf_num_list.begin(); _it != _inf_num_list.end(); _it++) { _init_states[i][*_it] = 1; }
	}
	set_seeds(_init_states, _time);
}

void metapopulation::set_seeds(std::vector<double>& _init_prob, double _time)
{
	std::vector<size_t> _seed_amounts(group_amount);
	for (size_t i = 0; i < group_amount; ++i)
	{
		_seed_amounts[i] = (size_t)(population_amounts[i] * _init_prob[i]);
	}
	set_seeds(_seed_amounts, _time);
}

void metapopulation::set_generator_seed(int _seed)
{
	generator.seed(_seed);
}

void metapopulation::spread_once()
{
	std::binomial_distribution<size_t> _distribution;
	std::vector<double> trans_inf(group_amount, 0);
	node* _node_ptr = nullptr;
	std::list<size_t> _inf_list;
	size_t _sv_num;
	double _inf_prob;
	double _actl_step;
	for (size_t i = 0; i < group_amount; ++i)
	{
		for (std::list<node*>::const_iterator _node_it = groups[i]->state_group_map['i'].begin(); _node_it != groups[i]->state_group_map['i'].end(); ++_node_it)
		{
			if ((*_node_it)->get_rem_time() - time < step) _actl_step = (*_node_it)->get_rem_time() - time;
			else _actl_step = step;
			hazard_inf->param((*_node_it)->get_params_inf());
			trans_inf[i] = trans_inf[i] + (*hazard_inf)(time - (*_node_it)->get_trans_time(), _actl_step);
		}
		trans_inf[i] = 1 - trans_inf[i] / groups[i]->get_node_amount();
	}
	std::vector<char> keys = { 's', 'i', 'w', 'd', 'v', 'p' };
	std::list<node*>::iterator _current_node_it;
	std::list<node*>::iterator _next_node_it;
	for (size_t i = 0; i < group_amount; ++i)
	{
		_inf_prob = 0;
		for (size_t j = 0; j < group_amount; ++j)
		{
			_inf_prob += (1 - pow(trans_inf[j], k_value * contacts[i][j])) * population_amounts[j] / node_amount;
		}
		_distribution.param(std::binomial_distribution<size_t>::param_type(s_population_amounts[i] + v_population_amounts[i], _inf_prob));
		_distribution.reset();
		_inf_list = diff_rand(_distribution(generator), 0, s_population_amounts[i] + v_population_amounts[i], generator);
		std::list<size_t>::const_iterator _list_it = _inf_list.begin();
		_sv_num = 0;
		for (std::vector<char>::const_iterator _key_it = keys.begin(); _key_it != keys.end(); ++_key_it)
		{
			if (*_key_it == 's' || *_key_it == 'v')
			{
				if (_list_it != _inf_list.end())
				{
					_next_node_it = groups[i]->state_group_map[*_key_it].begin();
					while (_next_node_it != groups[i]->state_group_map[*_key_it].end())
					{
						_current_node_it = _next_node_it;
						++_next_node_it;
						_node_ptr = *_current_node_it;
						if ((*_list_it) == _sv_num)
						{
							_node_ptr->set_state('i', time); 
							dist_rem->param(_node_ptr->get_params_rem());
							_node_ptr->set_rem_time(time + (*dist_rem)(generator));
							_node_ptr->invalidate_remaining_ptime();
							groups[i]->state_group_map[*_key_it].erase(_current_node_it);
							groups[i]->state_group_map['i'].push_back(_node_ptr);
							if (*_key_it == 's')
							{
								--s_amount;
								--s_population_amounts[i];
							}
							else
							{
								--v_amount;
								--v_population_amounts[i];
							}
							++i_amount;
							++i_population_amounts[i];
							++c_amount;
							++c_population_amounts[i];
							++_list_it;
							if (_list_it == _inf_list.end()) break;
						}
						++_sv_num;
					}
				}
				if (*_key_it == 'v')
				{
					_next_node_it = groups[i]->state_group_map[*_key_it].begin();
					while (_next_node_it != groups[i]->state_group_map[*_key_it].end())
					{
						_current_node_it = _next_node_it;
						++_next_node_it;
						_node_ptr = *_current_node_it;
						if (_node_ptr->get_remaining_ptime() == 0)
						{
							_node_ptr->set_state('p', time);
							_node_ptr->invalidate_remaining_ptime();
							_node_ptr->invalidate_rem_time();
							groups[i]->state_group_map[*_key_it].erase(_current_node_it);
							groups[i]->state_group_map['p'].push_back(_node_ptr);
							--v_amount;
							--v_population_amounts[i];
							++p_amount;
							++p_population_amounts[i];
						}
						else
						{
							_node_ptr->remaining_ptime_dec();
						}
					}
				}
			}
			else if (*_key_it == 'i')
			{
				_next_node_it = groups[i]->state_group_map[*_key_it].begin();
				while (_next_node_it != groups[i]->state_group_map[*_key_it].end())
				{
					_current_node_it = _next_node_it;
					++_next_node_it;
					_node_ptr = *_current_node_it;
					if (_node_ptr->get_rem_time() <= time + step)
					{
						if (prob(ifrs[i], generator))
						{
							_node_ptr->set_state('d', _node_ptr->get_rem_time());
							groups[i]->state_group_map[*_key_it].erase(_current_node_it);
							groups[i]->state_group_map['d'].push_back(_node_ptr);
							++d_amount;
							++d_population_amounts[i];
							y_amount = y_amount + ylls[i];
							y_population_amounts[i] = y_population_amounts[i] + ylls[i];
						}
						else
						{
							_node_ptr->set_state('w', _node_ptr->get_rem_time());
							groups[i]->state_group_map[*_key_it].erase(_current_node_it);
							groups[i]->state_group_map['w'].push_back(_node_ptr);
							++w_amount;
							++w_population_amounts[i];
						}
						_node_ptr->invalidate_rem_time();
						++r_amount;
						++r_population_amounts[i];
						--i_amount;
						--i_population_amounts[i];
					}
				}
			}
			else;
		}
	}
	time = time + step;
}

size_t metapopulation::add_vaccine(std::vector<size_t>& _allocation)
{
	size_t _vac_remaining = 0;
	size_t _vac_alloc = 0;
	node* _node_ptr = nullptr;
	std::list<size_t> _vac_list;
	size_t _s_num;
	std::list<node*>::iterator _current_node_it;
	std::list<node*>::iterator _next_node_it;
	for (size_t i = 0; i < group_amount; ++i)
	{
		_vac_alloc = _allocation[i];
		if (_vac_alloc > s_population_amounts[i])
		{
			_vac_remaining += (_vac_alloc - s_population_amounts[i]);
			_vac_alloc = s_population_amounts[i];
		}
		_vac_list = diff_rand(_vac_alloc, 0, s_population_amounts[i], generator);
		std::list<size_t>::const_iterator _list_it = _vac_list.begin();
		_s_num = 0;
		_next_node_it = groups[i]->state_group_map['s'].begin();
		while(_next_node_it != groups[i]->state_group_map['s'].end())
		{
			_current_node_it = _next_node_it;
			++_next_node_it;
			if (_list_it != _vac_list.end() && (*_list_it) == _s_num)
			{
				if (prob(eta, generator))
				{
					_node_ptr = *_current_node_it;
					_node_ptr->set_state('v', time);
					_node_ptr->set_remaining_ptime(int(delay));
					groups[i]->state_group_map['s'].erase(_current_node_it);
					groups[i]->state_group_map['v'].push_back(_node_ptr);
					--s_amount;
					--s_population_amounts[i];
					++v_amount;
					++v_population_amounts[i];
				}
				++_list_it;
			}
			++_s_num;
		}
	}
	return _vac_remaining;
}

size_t metapopulation::get_node_amount() const
{
	return node_amount;
}

size_t metapopulation::get_group_amount() const
{
	return group_amount;
}

double metapopulation::get_time() const
{
	return time;
}

size_t metapopulation::get_s_amount() const
{
	return s_amount;
}

size_t metapopulation::get_i_amount() const
{
	return i_amount;
}

size_t metapopulation::get_r_amount() const
{
	return r_amount;
}

size_t metapopulation::get_w_amount() const
{
	return w_amount;
}

size_t metapopulation::get_c_amount() const
{
	return c_amount;
}

size_t metapopulation::get_d_amount() const
{
	return d_amount;
}

double metapopulation::get_y_amount() const
{
	return y_amount;
}

size_t metapopulation::get_v_amount() const
{
	return v_amount;
}

size_t metapopulation::get_p_amount() const
{
	return p_amount;
}

std::vector<size_t> metapopulation::get_s_amount_arr() const
{
	return s_population_amounts;
}

std::vector<size_t> metapopulation::get_i_amount_arr() const
{
	return i_population_amounts;
}

std::vector<size_t> metapopulation::get_r_amount_arr() const
{
	return r_population_amounts;
}

std::vector<size_t> metapopulation::get_w_amount_arr() const
{
	return w_population_amounts;
}

std::vector<size_t> metapopulation::get_c_amount_arr() const
{
	return c_population_amounts;
}

std::vector<size_t> metapopulation::get_d_amount_arr() const
{
	return d_population_amounts;
}

std::vector<double> metapopulation::get_y_amount_arr() const
{
	return y_population_amounts;
}

std::vector<size_t> metapopulation::get_v_amount_arr() const
{
	return v_population_amounts;
}

std::vector<size_t> metapopulation::get_p_amount_arr() const
{
	return p_population_amounts;
}

size_t metapopulation::get_s_amount_idx(size_t _idx) const
{
	return s_population_amounts[_idx];
}

size_t metapopulation::get_i_amount_idx(size_t _idx) const
{
	return i_population_amounts[_idx];
}

size_t metapopulation::get_r_amount_idx(size_t _idx) const
{
	return r_population_amounts[_idx];
}

size_t metapopulation::get_w_amount_idx(size_t _idx) const
{
	return w_population_amounts[_idx];
}

size_t metapopulation::get_c_amount_idx(size_t _idx) const
{
	return c_population_amounts[_idx];
}

size_t metapopulation::get_d_amount_idx(size_t _idx) const
{
	return d_population_amounts[_idx];
}

double metapopulation::get_y_amount_idx(size_t _idx) const
{
	return y_population_amounts[_idx];
}

size_t metapopulation::get_v_amount_idx(size_t _idx) const
{
	return v_population_amounts[_idx];
}

size_t metapopulation::get_p_amount_idx(size_t _idx) const
{
	return p_population_amounts[_idx];
}

double metapopulation::get_node_trans_time(size_t _group_idx, size_t _node_idx) const
{
	return groups[_group_idx]->all_nodes[_node_idx]->get_trans_time();
}

double metapopulation::get_node_rem_time(size_t _group_idx, size_t _node_idx) const
{
	return groups[_group_idx]->all_nodes[_node_idx]->get_rem_time();
}

char metapopulation::get_node_state(size_t _group_idx, size_t _node_idx) const
{
	return groups[_group_idx]->all_nodes[_node_idx]->get_state();
}

std::pair<std::vector<double>, std::vector<double>> const& metapopulation::get_node_spreading_params(size_t _group_idx, size_t _node_idx)
{
	return *(new std::pair<std::vector<double>, std::vector<double>>{ groups[_group_idx]->all_nodes[_node_idx]->get_params_inf(),
																	  groups[_group_idx]->all_nodes[_node_idx]->get_params_rem()});
}

size_t metapopulation::get_delay() const
{
	return delay;
}

double metapopulation::get_eta() const
{
	return eta;
}

double metapopulation::get_step() const
{
	return step;
}

std::vector<double> const& metapopulation::get_populations() const
{
	return populations;
}

std::vector<std::vector<double>> const& metapopulation::get_contacts() const
{
	return contacts;
}

std::vector<double> const& metapopulation::get_ifrs() const
{
	return ifrs;
}

std::vector<size_t> const& metapopulation::get_population_amounts() const
{
	return population_amounts;
}