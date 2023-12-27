#include "model_elem.h"


std::list<size_t>& diff_rand(size_t _num, size_t _a, size_t _b, std::default_random_engine& _generator_ref)
{
	std::list<size_t>* _ptr = new std::list<size_t>;
	if ((_b <= _a) || (_num > (_b - _a))) return *_ptr;
	std::uniform_int_distribution<size_t> _distribution;
	size_t _rand_num;
	for (size_t i = 0; i < _num; ++i)
	{
		_distribution.param(std::uniform_int_distribution<size_t>::param_type(_a, _b - i - 1));
		_distribution.reset();
		_rand_num = _distribution(_generator_ref);
		std::list<size_t>::iterator _it = _ptr->begin();
		while (true)
		{
			if (_it == _ptr->end())
			{
				_ptr->insert(_it, _rand_num);
				break;
			}
			if (_rand_num < (*_it))
			{
				_ptr->insert(_it, _rand_num);
				break;
			}
			else ++_rand_num;
			++_it;
		}
	}
	return *_ptr;
}

probability::probability() : distribution(0, 1.0){}

bool probability::operator()(double _prob, std::default_random_engine& _generator_ref)
{
	if (_prob <= 0) return false;
	if (_prob >= 1) return true;
	while (_prob < 0.5)
	{
		_prob *= 2;
		distribution.reset();
		if (distribution(_generator_ref) < 0.5) return false;
	}
	distribution.reset();
	if (distribution(_generator_ref) < _prob) return true;
	return false;
}


node::node() : state('s'), trans_time(0), rem_time(-1), remaining_ptime(-1) {}


void node::init()
{
	state = 's';
	trans_time = 0;
	rem_time = -1;
	remaining_ptime = -1;
}

void node::set_spreading_params(std::vector<double> const& _params_inf, std::vector<double> const& _params_rem)
{
	params_inf = _params_inf;
	params_rem = _params_rem;
}

void node::set_rem_time(double _rem_time)
{
	rem_time = _rem_time;
}

void node::set_remaining_ptime(int _remaining_ptime)
{
	remaining_ptime = _remaining_ptime;
}

void node::remaining_ptime_dec()
{
	--remaining_ptime;
}

void node::invalidate_rem_time()
{
	rem_time = -1;
}

void node::invalidate_remaining_ptime()
{
	remaining_ptime = -1;
}


char node::get_state() const
{
	return state;
}

double node::get_trans_time() const
{
	return trans_time;
}

void node::set_state(char _state, double _trans_time)
{
	state = _state;
	trans_time = _trans_time;
}


double node::get_rem_time() const
{
	return rem_time;
}

int node::get_remaining_ptime() const
{
	return remaining_ptime;
}

std::vector<double> node::get_params_inf() const
{
	return params_inf;
}

std::vector<double> node::get_params_rem() const
{
	return params_rem;
}



state_groups::state_groups(size_t _node_amount)
{
	state_group_map['s'] = std::list<node*>();
	state_group_map['i'] = std::list<node*>();
	state_group_map['w'] = std::list<node*>();
	state_group_map['d'] = std::list<node*>();
	state_group_map['v'] = std::list<node*>();
	state_group_map['p'] = std::list<node*>();
	node* _node_ptr = nullptr;
	for (size_t i = 0; i < _node_amount; ++i)
	{
		_node_ptr = new node();
		state_group_map['s'].push_back(_node_ptr);
		all_nodes.push_back(_node_ptr);
	}
}

state_groups::~state_groups()
{
	std::vector<char> keys = {'s', 'i', 'w', 'd', 'v', 'p'};
	for (std::vector<char>::const_iterator _key_it = keys.begin(); _key_it != keys.end(); ++_key_it)
	{ 
		for (std::list<node*>::iterator _it = state_group_map[*_key_it].begin(); _it != state_group_map[*_key_it].end(); ++_it) 
			delete (*_it);
		state_group_map[*_key_it].clear();
	}
	all_nodes.clear();
}

size_t state_groups::get_node_amount() const
{
	return all_nodes.size();
}

size_t state_groups::get_state_node_amount(char _state) const
{
	if (_state == 's' || _state == 'i' || _state == 'w' ||  _state == 'd' ||  _state == 'v' || _state == 'p')
		state_group_map.find(_state)->second.size();
	return 0;
}

