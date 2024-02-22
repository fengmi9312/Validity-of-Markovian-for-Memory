#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "structmember.h"
#include "metapopulation.h"


std::vector<size_t>& PyObject_ToVectorSize_t_1D(PyObject* _pyobj)
{
	PyObject* _pyvec = PySequence_Fast(_pyobj, "argument must be iterable");
	if (!_pyvec) return *(new std::vector<size_t>);
	size_t _len = PySequence_Fast_GET_SIZE(_pyvec);
	std::vector<size_t>* _result = new std::vector<size_t>(_len);
	for (size_t i = 0; i < _len; ++i)
	{
		PyObject* _item = PySequence_Fast_GET_ITEM(_pyvec, i);
		if (!_item) {
			Py_DECREF(_pyvec);
			delete _result;
			return *(new std::vector<size_t>);
		}
		PyObject* _size_t_item = PyNumber_Long(_item);
		if (!_size_t_item) {
			Py_DECREF(_pyvec);
			delete _result;
			PyErr_SetString(PyExc_TypeError, "all items must be numbers");
			return *(new std::vector<size_t>);
		}
		(*_result)[i] = PyLong_AsSize_t(_size_t_item);
		Py_DECREF(_size_t_item);
	}
	Py_DECREF(_pyvec);
	return *_result;
}

std::vector<std::vector<size_t>>& PyObject_ToVectorSize_t_2D(PyObject* _pyobj)
{
	PyObject* _pyvec = PySequence_Fast(_pyobj, "argument must be iterable");
	if (!_pyvec) return *(new std::vector<std::vector<size_t>>);
	size_t _len = PySequence_Fast_GET_SIZE(_pyvec);
	std::vector<std::vector<size_t>>* _result = new std::vector<std::vector<size_t>>(_len);
	for (size_t i = 0; i < _len; ++i)
	{
		PyObject* _item = PySequence_Fast_GET_ITEM(_pyvec, i);
		if (!_item) {
			Py_DECREF(_pyvec);
			delete _result;
			return *(new std::vector<std::vector<size_t>>);
		}
		(*_result)[i] = PyObject_ToVectorSize_t_1D(_item);
	}
	Py_DECREF(_pyvec);
	return *_result;
}

std::vector<double>& PyObject_ToVectorDouble_1D(PyObject* _pyobj)
{
	PyObject* _pyvec = PySequence_Fast(_pyobj, "argument must be iterable");
	if (!_pyvec) return *(new std::vector<double>);
	size_t _len = PySequence_Fast_GET_SIZE(_pyvec);
	std::vector<double>* _result = new std::vector<double>(_len);
	for (size_t i = 0; i < _len; ++i)
	{
		PyObject* _item = PySequence_Fast_GET_ITEM(_pyvec, i);
		if (!_item) {
			Py_DECREF(_pyvec);
			delete _result;
			return *(new std::vector<double>);
		}
		PyObject* _double_item = PyNumber_Float(_item);
		if (!_double_item) {
			Py_DECREF(_pyvec);
			delete _result;
			PyErr_SetString(PyExc_TypeError, "all items must be numbers");
			return *(new std::vector<double>);
		}
		(*_result)[i] = PyFloat_AsDouble(_double_item);
		Py_DECREF(_double_item);
	}
	Py_DECREF(_pyvec);
	return *_result;
}

std::vector<std::vector<double>>& PyObject_ToVectorDouble_2D(PyObject* _pyobj)
{
	PyObject* _pyvec = PySequence_Fast(_pyobj, "argument must be iterable");
	if (!_pyvec) return *(new std::vector<std::vector<double>>);
	size_t _len = PySequence_Fast_GET_SIZE(_pyvec);
	std::vector<std::vector<double>>* _result = new std::vector<std::vector<double>>(_len);
	for (size_t i = 0; i < _len; ++i)
	{
		PyObject* _item = PySequence_Fast_GET_ITEM(_pyvec, i);
		if (!_item) {
			Py_DECREF(_pyvec);
			delete _result;
			return *(new std::vector<std::vector<double>>);
		}
		(*_result)[i] = PyObject_ToVectorDouble_1D(_item);
	}
	Py_DECREF(_pyvec);
	return *_result;
}

std::vector<std::vector<std::vector<double>>>& PyObject_ToVectorDouble_3D(PyObject* _pyobj)
{
	PyObject* _pyvec = PySequence_Fast(_pyobj, "argument must be iterable");
	if (!_pyvec) return *(new std::vector<std::vector<std::vector<double>>>);
	size_t _len = PySequence_Fast_GET_SIZE(_pyvec);
	std::vector<std::vector<std::vector<double>>>* _result = new std::vector<std::vector<std::vector<double>>>(_len);
	for (size_t i = 0; i < _len; ++i)
	{
		PyObject* _item = PySequence_Fast_GET_ITEM(_pyvec, i);
		if (!_item) {
			Py_DECREF(_pyvec);
			delete _result;
			return *(new std::vector<std::vector<std::vector<double>>>);
		}
		(*_result)[i] = PyObject_ToVectorDouble_2D(_item);
	}
	Py_DECREF(_pyvec);
	return *_result;
}




typedef struct {
	PyObject_HEAD
	metapopulation* simu_ptr;
} simuObject;

static void
simu_dealloc(simuObject* self)
{
	if (self->simu_ptr != nullptr)
		delete self->simu_ptr;
	Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject*
simu_new(PyTypeObject* type, PyObject* args, PyObject* kwds)
{
	simuObject* self;
	self = (simuObject*)type->tp_alloc(type, 0);
	if (self != NULL) {
		self->simu_ptr = nullptr;
	}
	return (PyObject*)self;
}

static int
simu_init(simuObject* self, PyObject* args, PyObject* kwds)
{
	static char* kwlist[] = { (char*)"node_amount", (char*)"group_amount", (char*)"populations", (char*)"contacts", (char*)"k", (char*)"ifrs", (char*)"ylls", (char*)"step", (char*)"time", NULL};
	size_t _node_amount;
	size_t _group_amount;
	PyObject* _population_array;
	PyObject* _contacts_matrix;
	double _k_value;
	PyObject* _ifr_array;
	PyObject* _yll_array;
	double _step = 0.01;
	double _time = 0;
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "KKOOdOO|dd", kwlist,
		&_node_amount, &_group_amount, &_population_array, &_contacts_matrix, &_k_value, &_ifr_array, &_yll_array, &_step, &_time)) 
	{ return -1; }
	self->simu_ptr = new metapopulation(_node_amount, _group_amount, PyObject_ToVectorDouble_1D(_population_array), PyObject_ToVectorDouble_2D(_contacts_matrix), 
										PyObject_ToVectorDouble_1D(_ifr_array), PyObject_ToVectorDouble_1D(_yll_array), _k_value, _step, _time);
	return 0;
}

static PyMemberDef simu_members[] = {
	{NULL}  /* Sentinel */
};


static PyObject* 
set_vaccination_params(simuObject* self, PyObject* args, PyObject* kwds)
{
	static char* kwlist[] = { (char*)"delay", (char*)"eta", NULL };
	size_t _delay = 14;
	double _eta = 1;
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "|Kd", kwlist, &_delay, &_eta)) {
		return NULL;
	}
	self->simu_ptr->set_vaccination_params(_delay, _eta);
	Py_INCREF(Py_None);
	return Py_None;
}

static PyObject*
set_spreading_func(simuObject* self, PyObject* args, PyObject* kwds)
{
	static char* kwlist[] = { (char*)"func_inf", (char*)"dist_rem", NULL};
	char* _func_inf_name_ptr;
	char* _dist_rem_name_ptr;
	hazard_type* _hazard_inf = nullptr;
	dist_type* _dist_rem = nullptr;
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "ss", kwlist, &_func_inf_name_ptr, &_dist_rem_name_ptr)) {
		return NULL;
	}
	std::string _hazfunc_inf_name(_func_inf_name_ptr);
	std::string _dist_rem_name(_dist_rem_name_ptr);
	if (_hazfunc_inf_name == "weibull") _hazard_inf = new weibull_hazard;
	else if (_hazfunc_inf_name == "xburr") _hazard_inf = new xburr_hazard;
	else if (_hazfunc_inf_name == "xlog_logistic") _hazard_inf = new xlog_logistic_hazard;
	else if (_hazfunc_inf_name == "general") _hazard_inf = new general_hazard;
	else;
	if (_dist_rem_name == "weibull") _dist_rem = new weibull_dist;
	else if (_dist_rem_name == "xburr") _dist_rem = new xburr_dist;
	else if (_dist_rem_name == "xlog_logistic") _dist_rem = new xlog_logistic_dist;
	else if (_dist_rem_name == "general") _dist_rem = new general_dist;
	else;
	self->simu_ptr->set_spreading_func(_hazard_inf, _dist_rem);
	Py_INCREF(Py_None);
	return Py_None;
}

static PyObject*
set_total_spreading_params(simuObject* self, PyObject* args, PyObject* kwds)
{
	static char* kwlist[] = { (char*)"params_inf", (char*)"params_rem", NULL };
	PyObject* _params_inf;
	PyObject* _params_rem;
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OO", kwlist, &_params_inf, &_params_rem)) {
		return NULL;
	}
	self->simu_ptr->set_spreading_params(PyObject_ToVectorDouble_1D(_params_inf), PyObject_ToVectorDouble_1D(_params_rem));
	Py_INCREF(Py_None);
	return Py_None;
}

static PyObject* 
set_detailed_spreading_params(simuObject* self, PyObject* args, PyObject* kwds)
{
	static char* kwlist[] = { (char*)"params_inf", (char*)"params_rem", NULL };
	PyObject* _params_inf;
	PyObject* _params_rem;
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OO", kwlist, &_params_inf, &_params_rem)) {
		return NULL;
	}
	self->simu_ptr->set_spreading_params(PyObject_ToVectorDouble_3D(_params_inf), PyObject_ToVectorDouble_3D(_params_rem));
	Py_INCREF(Py_None);
	return Py_None;
}

static PyObject*
set_total_spreading_survivals(simuObject* self, PyObject* args, PyObject* kwds)
{
	static char* kwlist[] = { (char*)"params_inf", (char*)"params_rem", NULL };
	PyObject* _survivals_inf;
	PyObject* _survivals_rem;
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OO", kwlist, &_survivals_inf, &_survivals_rem)) {
		return NULL;
	}
	self->simu_ptr->set_spreading_survivals(PyObject_ToVectorDouble_1D(_survivals_inf), PyObject_ToVectorDouble_1D(_survivals_rem));
	Py_INCREF(Py_None);
	return Py_None;
}

static PyObject*
init_all_states(simuObject* self, PyObject* Py_UNUSED(ignored))
{
	self->simu_ptr->init_all_states();
	Py_INCREF(Py_None);
	return Py_None;
}

static PyObject* 
set_prob_seeds(simuObject* self, PyObject* args, PyObject* kwds)
{
	static char* kwlist[] = { (char*)"init_prob", (char*)"time", NULL};
	PyObject* _init_prob_array;
	double _time = 0;
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|d", kwlist, &_init_prob_array, &_time)){
		return NULL;
	}
	self->simu_ptr->set_seeds(PyObject_ToVectorDouble_1D(_init_prob_array), _time);
	Py_INCREF(Py_None);
	return Py_None;
}

static PyObject* 
set_amount_seeds(simuObject* self, PyObject* args, PyObject* kwds)
{
	static char* kwlist[] = { (char*)"init_amount", (char*)"time", NULL };
	PyObject* _init_amount_array;
	double _time = 0;
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|d", kwlist, &_init_amount_array, &_time)) {
		return NULL;
	}
	self->simu_ptr->set_seeds(PyObject_ToVectorSize_t_1D(_init_amount_array), _time);
	Py_INCREF(Py_None);
	return Py_None;
}

static PyObject* 
set_state_seeds(simuObject* self, PyObject* args, PyObject* kwds)
{
	static char* kwlist[] = { (char*)"init_state", (char*)"time", NULL };
	PyObject* _init_state_array;
	double _time = 0;
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|d", kwlist, &_init_state_array, &_time)) {
		return NULL;
	}
	self->simu_ptr->set_seeds(PyObject_ToVectorSize_t_2D(_init_state_array), _time);
	Py_INCREF(Py_None);
	return Py_None;
}

static PyObject* 
set_generator_seed(simuObject* self, PyObject* args, PyObject* kwds)
{
	static char* kwlist[] = { (char*)"seed", NULL };
	int _seed;
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "i", kwlist, &_seed)) {
		return NULL;
	}
	self->simu_ptr->set_generator_seed(_seed);
	Py_INCREF(Py_None);
	return Py_None;
}

static PyObject* 
adjust_rem_time(simuObject* self, PyObject* Py_UNUSED(ignored))
{
	self->simu_ptr->adjust_rem_time();
	Py_INCREF(Py_None);
	return Py_None;
}

static PyObject*
reset_trans_time(simuObject* self, PyObject* Py_UNUSED(ignored))
{
	self->simu_ptr->reset_trans_time();
	Py_INCREF(Py_None);
	return Py_None;
}

static PyObject*
reset_rem_time(simuObject* self, PyObject* Py_UNUSED(ignored))
{
	self->simu_ptr->reset_rem_time();
	Py_INCREF(Py_None);
	return Py_None;
}

static PyObject* 
spread_once(simuObject* self, PyObject* Py_UNUSED(ignored))
{
	self->simu_ptr->spread_once();
	Py_INCREF(Py_None);
	return Py_None;
}

static PyObject*
spread(simuObject* self, PyObject* args, PyObject* kwds)
{
	static char* kwlist[] = { (char*)"times", NULL };
	size_t _times = 1;
	PyObject* PyList_s_ptr = PyList_New(0);
	PyObject* PyList_i_ptr = PyList_New(0);
	PyObject* PyList_r_ptr = PyList_New(0);
	PyObject* PyList_w_ptr = PyList_New(0);
	PyObject* PyList_c_ptr = PyList_New(0);
	PyObject* PyList_d_ptr = PyList_New(0);
	PyObject* PyList_y_ptr = PyList_New(0);
	PyObject* PyList_v_ptr = PyList_New(0);
	PyObject* PyList_p_ptr = PyList_New(0);
	PyObject* _item_s = PyLong_FromSize_t(self->simu_ptr->get_s_amount());
	PyList_Append(PyList_s_ptr, _item_s);
	Py_DECREF(_item_s);
	PyObject* _item_i = PyLong_FromSize_t(self->simu_ptr->get_i_amount());
	PyList_Append(PyList_i_ptr, _item_i);
	Py_DECREF(_item_i);
	PyObject* _item_r = PyLong_FromSize_t(self->simu_ptr->get_r_amount());
	PyList_Append(PyList_r_ptr, _item_r);
	Py_DECREF(_item_r);
	PyObject* _item_w = PyLong_FromSize_t(self->simu_ptr->get_w_amount());
	PyList_Append(PyList_w_ptr, _item_w);
	Py_DECREF(_item_w);
	PyObject* _item_c = PyLong_FromSize_t(self->simu_ptr->get_c_amount());
	PyList_Append(PyList_c_ptr, _item_c);
	Py_DECREF(_item_c);
	PyObject* _item_d = PyLong_FromSize_t(self->simu_ptr->get_d_amount());
	PyList_Append(PyList_d_ptr, _item_d);
	Py_DECREF(_item_d);
	PyObject* _item_y = PyFloat_FromDouble(self->simu_ptr->get_y_amount());
	PyList_Append(PyList_y_ptr, _item_y);
	Py_DECREF(_item_y);
	PyObject* _item_v = PyLong_FromSize_t(self->simu_ptr->get_v_amount());
	PyList_Append(PyList_v_ptr, _item_v);
	Py_DECREF(_item_v);
	PyObject* _item_p = PyLong_FromSize_t(self->simu_ptr->get_p_amount());
	PyList_Append(PyList_p_ptr, _item_p);
	Py_DECREF(_item_p);
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "|K", kwlist, &_times)) {
		return NULL;
	}
	for (size_t i = 0; i < _times; ++i)
	{
		self->simu_ptr->spread_once();
		_item_s = PyLong_FromSize_t(self->simu_ptr->get_s_amount());
		PyList_Append(PyList_s_ptr, _item_s);
		Py_DECREF(_item_s);
		_item_i = PyLong_FromSize_t(self->simu_ptr->get_i_amount());
		PyList_Append(PyList_i_ptr, _item_i);
		Py_DECREF(_item_i);
		_item_r = PyLong_FromSize_t(self->simu_ptr->get_r_amount());
		PyList_Append(PyList_r_ptr, _item_r);
		Py_DECREF(_item_r);
		_item_w = PyLong_FromSize_t(self->simu_ptr->get_w_amount());
		PyList_Append(PyList_w_ptr, _item_w);
		Py_DECREF(_item_w);
		_item_c = PyLong_FromSize_t(self->simu_ptr->get_c_amount());
		PyList_Append(PyList_c_ptr, _item_c);
		Py_DECREF(_item_c);
		_item_d = PyLong_FromSize_t(self->simu_ptr->get_d_amount());
		PyList_Append(PyList_d_ptr, _item_d);
		Py_DECREF(_item_d);
		_item_y = PyFloat_FromDouble(self->simu_ptr->get_y_amount());
		PyList_Append(PyList_y_ptr, _item_y);
		Py_DECREF(_item_y);
		_item_v = PyLong_FromSize_t(self->simu_ptr->get_v_amount());
		PyList_Append(PyList_v_ptr, _item_v);
		Py_DECREF(_item_v);
		_item_p = PyLong_FromSize_t(self->simu_ptr->get_p_amount());
		PyList_Append(PyList_p_ptr, _item_p);
		Py_DECREF(_item_p);
	}
	PyObject* _dict = PyDict_New();
	PyDict_SetItemString(_dict, "s", PyList_s_ptr);
	Py_DECREF(PyList_s_ptr);
	PyDict_SetItemString(_dict, "i", PyList_i_ptr);
	Py_DECREF(PyList_i_ptr);
	PyDict_SetItemString(_dict, "r", PyList_r_ptr);
	Py_DECREF(PyList_r_ptr);
	PyDict_SetItemString(_dict, "w", PyList_w_ptr);
	Py_DECREF(PyList_w_ptr);
	PyDict_SetItemString(_dict, "c", PyList_c_ptr);
	Py_DECREF(PyList_c_ptr);
	PyDict_SetItemString(_dict, "d", PyList_d_ptr);
	Py_DECREF(PyList_d_ptr);
	PyDict_SetItemString(_dict, "y", PyList_y_ptr);
	Py_DECREF(PyList_y_ptr);
	PyDict_SetItemString(_dict, "v", PyList_v_ptr);
	Py_DECREF(PyList_v_ptr);
	PyDict_SetItemString(_dict, "p", PyList_p_ptr);
	Py_DECREF(PyList_p_ptr);
	return _dict;
}

static PyObject*
spread_to_end(simuObject* self, PyObject* Py_UNUSED(ignored))
{
	PyObject* PyList_s_ptr = PyList_New(0);
	PyObject* PyList_i_ptr = PyList_New(0);
	PyObject* PyList_r_ptr = PyList_New(0);
	PyObject* PyList_w_ptr = PyList_New(0);
	PyObject* PyList_c_ptr = PyList_New(0);
	PyObject* PyList_d_ptr = PyList_New(0);
	PyObject* PyList_y_ptr = PyList_New(0);
	PyObject* PyList_v_ptr = PyList_New(0);
	PyObject* PyList_p_ptr = PyList_New(0);
	PyObject* _item_s = PyLong_FromSize_t(self->simu_ptr->get_s_amount());
	PyList_Append(PyList_s_ptr, _item_s);
	Py_DECREF(_item_s);
	PyObject* _item_i = PyLong_FromSize_t(self->simu_ptr->get_i_amount());
	PyList_Append(PyList_i_ptr, _item_i);
	Py_DECREF(_item_i);
	PyObject* _item_r = PyLong_FromSize_t(self->simu_ptr->get_r_amount());
	PyList_Append(PyList_r_ptr, _item_r);
	Py_DECREF(_item_r);
	PyObject* _item_w = PyLong_FromSize_t(self->simu_ptr->get_w_amount());
	PyList_Append(PyList_w_ptr, _item_w);
	Py_DECREF(_item_w);
	PyObject* _item_c = PyLong_FromSize_t(self->simu_ptr->get_c_amount());
	PyList_Append(PyList_c_ptr, _item_c);
	Py_DECREF(_item_c);
	PyObject* _item_d = PyLong_FromSize_t(self->simu_ptr->get_d_amount());
	PyList_Append(PyList_d_ptr, _item_d);
	Py_DECREF(_item_d);
	PyObject* _item_y = PyFloat_FromDouble(self->simu_ptr->get_y_amount());
	PyList_Append(PyList_y_ptr, _item_y);
	Py_DECREF(_item_y);
	PyObject* _item_v = PyLong_FromSize_t(self->simu_ptr->get_v_amount());
	PyList_Append(PyList_v_ptr, _item_v);
	Py_DECREF(_item_v);
	PyObject* _item_p = PyLong_FromSize_t(self->simu_ptr->get_p_amount());
	PyList_Append(PyList_p_ptr, _item_p);
	Py_DECREF(_item_p);
	while(self->simu_ptr->get_i_amount() != 0)
	{
		self->simu_ptr->spread_once();
		_item_s = PyLong_FromSize_t(self->simu_ptr->get_s_amount());
		PyList_Append(PyList_s_ptr, _item_s);
		Py_DECREF(_item_s);
		_item_i = PyLong_FromSize_t(self->simu_ptr->get_i_amount());
		PyList_Append(PyList_i_ptr, _item_i);
		Py_DECREF(_item_i);
		_item_r = PyLong_FromSize_t(self->simu_ptr->get_r_amount());
		PyList_Append(PyList_r_ptr, _item_r);
		Py_DECREF(_item_r);
		_item_w = PyLong_FromSize_t(self->simu_ptr->get_w_amount());
		PyList_Append(PyList_w_ptr, _item_w);
		Py_DECREF(_item_w);
		_item_c = PyLong_FromSize_t(self->simu_ptr->get_c_amount());
		PyList_Append(PyList_c_ptr, _item_c);
		Py_DECREF(_item_c);
		_item_d = PyLong_FromSize_t(self->simu_ptr->get_d_amount());
		PyList_Append(PyList_d_ptr, _item_d);
		Py_DECREF(_item_d);
		_item_y = PyFloat_FromDouble(self->simu_ptr->get_y_amount());
		PyList_Append(PyList_y_ptr, _item_y);
		Py_DECREF(_item_y);
		_item_v = PyLong_FromSize_t(self->simu_ptr->get_v_amount());
		PyList_Append(PyList_v_ptr, _item_v);
		Py_DECREF(_item_v);
		_item_p = PyLong_FromSize_t(self->simu_ptr->get_p_amount());
		PyList_Append(PyList_p_ptr, _item_p);
		Py_DECREF(_item_p);
	}
	PyObject* _dict = PyDict_New();
	PyDict_SetItemString(_dict, "s", PyList_s_ptr);
	Py_DECREF(PyList_s_ptr);
	PyDict_SetItemString(_dict, "i", PyList_i_ptr);
	Py_DECREF(PyList_i_ptr);
	PyDict_SetItemString(_dict, "r", PyList_r_ptr);
	Py_DECREF(PyList_r_ptr);
	PyDict_SetItemString(_dict, "w", PyList_w_ptr);
	Py_DECREF(PyList_w_ptr);
	PyDict_SetItemString(_dict, "c", PyList_c_ptr);
	Py_DECREF(PyList_c_ptr);
	PyDict_SetItemString(_dict, "d", PyList_d_ptr);
	Py_DECREF(PyList_d_ptr);
	PyDict_SetItemString(_dict, "y", PyList_y_ptr);
	Py_DECREF(PyList_y_ptr);
	PyDict_SetItemString(_dict, "v", PyList_v_ptr);
	Py_DECREF(PyList_v_ptr);
	PyDict_SetItemString(_dict, "p", PyList_p_ptr);
	Py_DECREF(PyList_p_ptr);
	return _dict;
}

static PyObject* 
add_vaccine(simuObject* self, PyObject* args, PyObject* kwds)
{
	static char* kwlist[] = { (char*)"allocation", NULL };
	PyObject* _allocation_array;
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "O", kwlist, &_allocation_array)) {
		return NULL;
	}
	PyObject* _res = PyLong_FromSize_t(self->simu_ptr->add_vaccine(PyObject_ToVectorSize_t_1D(_allocation_array)));
	return _res;
}

static PyObject* 
get_node_amount(simuObject* self, PyObject* Py_UNUSED(ignored))
{
	return PyLong_FromSize_t(self->simu_ptr->get_node_amount());
}

static PyObject* 
get_group_amount(simuObject* self, PyObject* Py_UNUSED(ignored))
{
	return PyLong_FromSize_t(self->simu_ptr->get_group_amount());
}

static PyObject* 
get_time(simuObject* self, PyObject* Py_UNUSED(ignored))
{
	return PyFloat_FromDouble(self->simu_ptr->get_time());
}

static PyObject* 
get_s_amount(simuObject* self, PyObject* Py_UNUSED(ignored))
{
	return PyLong_FromSize_t(self->simu_ptr->get_s_amount());
}

static PyObject* 
get_i_amount(simuObject* self, PyObject* Py_UNUSED(ignored))
{
	return PyLong_FromSize_t(self->simu_ptr->get_i_amount());
}

static PyObject* 
get_r_amount(simuObject* self, PyObject* Py_UNUSED(ignored))
{
	return PyLong_FromSize_t(self->simu_ptr->get_r_amount());
}

static PyObject* 
get_w_amount(simuObject* self, PyObject* Py_UNUSED(ignored))
{
	return PyLong_FromSize_t(self->simu_ptr->get_w_amount());
}

static PyObject* 
get_c_amount(simuObject* self, PyObject* Py_UNUSED(ignored))
{
	return PyLong_FromSize_t(self->simu_ptr->get_c_amount());
}

static PyObject* 
get_d_amount(simuObject* self, PyObject* Py_UNUSED(ignored))
{
	return PyLong_FromSize_t(self->simu_ptr->get_d_amount());
}

static PyObject* 
get_y_amount(simuObject* self, PyObject* Py_UNUSED(ignored))
{
	return PyFloat_FromDouble(self->simu_ptr->get_y_amount());
}

static PyObject* 
get_v_amount(simuObject* self, PyObject* Py_UNUSED(ignored))
{
	return PyLong_FromSize_t(self->simu_ptr->get_v_amount());
}

static PyObject* 
get_p_amount(simuObject* self, PyObject* Py_UNUSED(ignored))
{
	return PyLong_FromSize_t(self->simu_ptr->get_p_amount());
}

static PyObject*
get_s_amount_arr(simuObject* self, PyObject* Py_UNUSED(ignored))
{
	PyObject* PyList_ptr = PyList_New(0);
	PyObject* _item = NULL;
	std::vector<size_t> const& _s_amount_arr = self->simu_ptr->get_s_amount_arr();
	for (std::vector<size_t>::const_iterator _it = _s_amount_arr.begin(); _it != _s_amount_arr.end(); ++_it)
	{
		_item = PyLong_FromSize_t(*_it);
		PyList_Append(PyList_ptr, _item);
		Py_DECREF(_item);
	}
	return PyList_ptr;
}

static PyObject*
get_i_amount_arr(simuObject* self, PyObject* Py_UNUSED(ignored))
{
	PyObject* PyList_ptr = PyList_New(0);
	PyObject* _item = NULL;
	std::vector<size_t> const& _i_amount_arr = self->simu_ptr->get_i_amount_arr();
	for (std::vector<size_t>::const_iterator _it = _i_amount_arr.begin(); _it != _i_amount_arr.end(); ++_it)
	{
		_item = PyLong_FromSize_t(*_it);
		PyList_Append(PyList_ptr, _item);
		Py_DECREF(_item);
	}
	return PyList_ptr;
}


static PyObject*
get_r_amount_arr(simuObject* self, PyObject* Py_UNUSED(ignored))
{
	PyObject* PyList_ptr = PyList_New(0);
	PyObject* _item = NULL;
	std::vector<size_t> const& _r_amount_arr = self->simu_ptr->get_r_amount_arr();
	for (std::vector<size_t>::const_iterator _it = _r_amount_arr.begin(); _it != _r_amount_arr.end(); ++_it)
	{
		_item = PyLong_FromSize_t(*_it);
		PyList_Append(PyList_ptr, _item);
		Py_DECREF(_item);
	}
	return PyList_ptr;
}


static PyObject*
get_w_amount_arr(simuObject* self, PyObject* Py_UNUSED(ignored))
{
	PyObject* PyList_ptr = PyList_New(0);
	PyObject* _item = NULL;
	std::vector<size_t> const& _w_amount_arr = self->simu_ptr->get_w_amount_arr();
	for (std::vector<size_t>::const_iterator _it = _w_amount_arr.begin(); _it != _w_amount_arr.end(); ++_it)
	{
		_item = PyLong_FromSize_t(*_it);
		PyList_Append(PyList_ptr, _item);
		Py_DECREF(_item);
	}
	return PyList_ptr;
}

static PyObject*
get_c_amount_arr(simuObject* self, PyObject* Py_UNUSED(ignored))
{
	PyObject* PyList_ptr = PyList_New(0);
	PyObject* _item = NULL;
	std::vector<size_t> const& _c_amount_arr = self->simu_ptr->get_c_amount_arr();
	for (std::vector<size_t>::const_iterator _it = _c_amount_arr.begin(); _it != _c_amount_arr.end(); ++_it)
	{
		_item = PyLong_FromSize_t(*_it);
		PyList_Append(PyList_ptr, _item);
		Py_DECREF(_item);
	}
	return PyList_ptr;
}

static PyObject*
get_d_amount_arr(simuObject* self, PyObject* Py_UNUSED(ignored))
{
	PyObject* PyList_ptr = PyList_New(0);
	PyObject* _item = NULL;
	std::vector<size_t> const& _d_amount_arr = self->simu_ptr->get_d_amount_arr();
	for (std::vector<size_t>::const_iterator _it = _d_amount_arr.begin(); _it != _d_amount_arr.end(); ++_it)
	{
		_item = PyLong_FromSize_t(*_it);
		PyList_Append(PyList_ptr, _item);
		Py_DECREF(_item);
	}
	return PyList_ptr;
}

static PyObject*
get_y_amount_arr(simuObject* self, PyObject* Py_UNUSED(ignored))
{
	PyObject* PyList_ptr = PyList_New(0);
	PyObject* _item = NULL;
	std::vector<double> const& _y_amount_arr = self->simu_ptr->get_y_amount_arr();
	for (std::vector<double>::const_iterator _it = _y_amount_arr.begin(); _it != _y_amount_arr.end(); ++_it)
	{
		_item = PyFloat_FromDouble(*_it);
		PyList_Append(PyList_ptr, _item);
		Py_DECREF(_item);
	}
	return PyList_ptr;
}

static PyObject*
get_v_amount_arr(simuObject* self, PyObject* Py_UNUSED(ignored))
{
	PyObject* PyList_ptr = PyList_New(0);
	PyObject* _item = NULL;
	std::vector<size_t> const& _v_amount_arr = self->simu_ptr->get_v_amount_arr();
	for (std::vector<size_t>::const_iterator _it = _v_amount_arr.begin(); _it != _v_amount_arr.end(); ++_it)
	{
		_item = PyLong_FromSize_t(*_it);
		PyList_Append(PyList_ptr, _item);
		Py_DECREF(_item);
	}
	return PyList_ptr;
}


static PyObject*
get_p_amount_arr(simuObject* self, PyObject* Py_UNUSED(ignored))
{
	PyObject* PyList_ptr = PyList_New(0);
	PyObject* _item = NULL;
	std::vector<size_t> const& _p_amount_arr = self->simu_ptr->get_p_amount_arr();
	for (std::vector<size_t>::const_iterator _it = _p_amount_arr.begin(); _it != _p_amount_arr.end(); ++_it)
	{
		_item = PyLong_FromSize_t(*_it);
		PyList_Append(PyList_ptr, _item);
		Py_DECREF(_item);
	}
	return PyList_ptr;
}


static PyObject* 
get_s_amount_idx(simuObject* self, PyObject* args, PyObject* kwds)
{
	static char* kwlist[] = { (char*)"idx", NULL };
	size_t _idx;
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "K", kwlist, &_idx)) {
		return NULL;
	}
	return PyLong_FromSize_t(self->simu_ptr->get_s_amount_idx(_idx));
}

static PyObject* 
get_i_amount_idx(simuObject* self, PyObject* args, PyObject* kwds)
{
	static char* kwlist[] = { (char*)"idx", NULL };
	size_t _idx;
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "K", kwlist, &_idx)) {
		return NULL;
	}
	return PyLong_FromSize_t(self->simu_ptr->get_i_amount_idx(_idx));
}

static PyObject* 
get_r_amount_idx(simuObject* self, PyObject* args, PyObject* kwds)
{
	static char* kwlist[] = { (char*)"idx", NULL };
	size_t _idx;
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "K", kwlist, &_idx)) {
		return NULL;
	}
	return PyLong_FromSize_t(self->simu_ptr->get_r_amount_idx(_idx));
}

static PyObject* 
get_w_amount_idx(simuObject* self, PyObject* args, PyObject* kwds)
{
	static char* kwlist[] = { (char*)"idx", NULL };
	size_t _idx;
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "K", kwlist, &_idx)) {
		return NULL;
	}
	return PyLong_FromSize_t(self->simu_ptr->get_w_amount_idx(_idx));
}

static PyObject* 
get_c_amount_idx(simuObject* self, PyObject* args, PyObject* kwds)
{
	static char* kwlist[] = { (char*)"idx", NULL };
	size_t _idx;
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "K", kwlist, &_idx)) {
		return NULL;
	}
	return PyLong_FromSize_t(self->simu_ptr->get_c_amount_idx(_idx));
}

static PyObject* 
get_d_amount_idx(simuObject* self, PyObject* args, PyObject* kwds)
{
	static char* kwlist[] = { (char*)"idx", NULL };
	size_t _idx;
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "K", kwlist, &_idx)) {
		return NULL;
	}
	return PyLong_FromSize_t(self->simu_ptr->get_d_amount_idx(_idx));
}

static PyObject* 
get_y_amount_idx(simuObject* self, PyObject* args, PyObject* kwds)
{
	static char* kwlist[] = { (char*)"idx", NULL };
	size_t _idx;
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "K", kwlist, &_idx)) {
		return NULL;
	}
	return PyFloat_FromDouble(self->simu_ptr->get_y_amount_idx(_idx));
}

static PyObject* 
get_v_amount_idx(simuObject* self, PyObject* args, PyObject* kwds)
{
	static char* kwlist[] = { (char*)"idx", NULL };
	size_t _idx;
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "K", kwlist, &_idx)) {
		return NULL;
	}
	return PyLong_FromSize_t(self->simu_ptr->get_v_amount_idx(_idx));
}

static PyObject* 
get_p_amount_idx(simuObject* self, PyObject* args, PyObject* kwds)
{
	static char* kwlist[] = { (char*)"idx", NULL };
	size_t _idx;
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "K", kwlist, &_idx)) {
		return NULL;
	}
	return PyLong_FromSize_t(self->simu_ptr->get_p_amount_idx(_idx));
}

static PyObject* 
get_delay(simuObject* self, PyObject* Py_UNUSED(ignored))
{
	return PyLong_FromSize_t(self->simu_ptr->get_delay());
}

static PyObject* 
get_eta(simuObject* self, PyObject* Py_UNUSED(ignored))
{
	return PyFloat_FromDouble(self->simu_ptr->get_eta());
}

static PyObject* 
get_step(simuObject* self, PyObject* Py_UNUSED(ignored))
{
	return PyFloat_FromDouble(self->simu_ptr->get_step());
}

static PyObject* 
get_populations(simuObject* self, PyObject* Py_UNUSED(ignored))
{
	PyObject* PyList_ptr = PyList_New(0);
	PyObject* _item = NULL;
	std::vector<double> const& _populations = self->simu_ptr->get_populations();
	for (std::vector<double>::const_iterator _it = _populations.begin(); _it != _populations.end(); ++_it)
	{
		_item = PyFloat_FromDouble(*_it);
		PyList_Append(PyList_ptr, _item);
		Py_DECREF(_item);
	}
	return PyList_ptr;
}

static PyObject* 
get_population_amounts(simuObject* self, PyObject* Py_UNUSED(ignored))
{
	PyObject* PyList_ptr = PyList_New(0);
	PyObject* _item = NULL;
	std::vector<size_t> const& _population_amounts = self->simu_ptr->get_population_amounts();
	for (std::vector<size_t>::const_iterator _it = _population_amounts.begin(); _it != _population_amounts.end(); ++_it)
	{ 
		_item = PyLong_FromSize_t(*_it);
		PyList_Append(PyList_ptr, _item);
		Py_DECREF(_item);
	}
	return PyList_ptr;
}

static PyObject* 
get_ifrs(simuObject* self, PyObject* Py_UNUSED(ignored))
{
	PyObject* PyList_ptr = PyList_New(0);
	PyObject* _item = NULL;
	std::vector<double> const& _ifr = self->simu_ptr->get_ifrs();
	for (std::vector<double>::const_iterator _it = _ifr.begin(); _it != _ifr.end(); ++_it)
	{
		_item = PyFloat_FromDouble(*_it);
		PyList_Append(PyList_ptr, _item);
		Py_DECREF(_item);
	}
	return PyList_ptr;
}

static PyObject* 
get_contacts(simuObject* self, PyObject* Py_UNUSED(ignored))
{
	PyObject* PyListList_ptr = PyList_New(0);
	PyObject* PyList_ptr = NULL;
	PyObject* _item = NULL;
	std::vector<std::vector<double>> const& _contacts = self->simu_ptr->get_contacts();
	for (std::vector<std::vector<double>>::const_iterator _it = _contacts.begin(); _it != _contacts.end(); ++_it)
	{
		PyList_ptr = PyList_New(0);
		for (std::vector<double>::const_iterator __it = _it->begin(); __it != _it->end(); ++__it)
		{
			_item = PyFloat_FromDouble(*__it);
			PyList_Append(PyList_ptr, PyFloat_FromDouble(*__it));
			Py_DECREF(_item);
		}
		PyList_Append(PyListList_ptr, PyList_ptr);
		Py_DECREF(PyList_ptr);
	}
	return PyListList_ptr;
}

static PyObject* 
get_x_amount(simuObject* self, PyObject* args, PyObject* kwds)
{
	static char* kwlist[] = { (char*)"target", NULL };
	char* _target_ptr;
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "s", kwlist, &_target_ptr)) {
		return NULL;
	}
	std::string _target(_target_ptr, 1);
	if (_target == "s")
		return PyLong_FromSize_t(self->simu_ptr->get_s_amount());
	else if (_target == "i")
		return PyLong_FromSize_t(self->simu_ptr->get_i_amount());
	else if (_target == "r")
		return PyLong_FromSize_t(self->simu_ptr->get_r_amount());
	else if (_target == "w")
		return PyLong_FromSize_t(self->simu_ptr->get_w_amount());
	else if (_target == "c")
		return PyLong_FromSize_t(self->simu_ptr->get_c_amount());
	else if (_target == "d")
		return PyLong_FromSize_t(self->simu_ptr->get_d_amount());
	else if (_target == "y")
		return PyFloat_FromDouble(self->simu_ptr->get_y_amount());
	else if (_target == "v")
		return PyLong_FromSize_t(self->simu_ptr->get_v_amount());
	else if (_target == "p")
		return PyLong_FromSize_t(self->simu_ptr->get_p_amount());
	else
	{
		Py_INCREF(Py_None);
		return Py_None;
	}
}

static PyObject*
get_x_amount_arr(simuObject* self, PyObject* args, PyObject* kwds)
{
	static char* kwlist[] = { (char*)"target", NULL };
	char* _target_ptr;
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "s", kwlist, &_target_ptr)) {
		return NULL;
	}
	std::string _target(_target_ptr, 1);
	PyObject* PyList_ptr = PyList_New(0);
	PyObject* _item = NULL;
	if (_target == "s")
	{
		std::vector<size_t> const& _amount_arr = self->simu_ptr->get_s_amount_arr(); 
		for (std::vector<size_t>::const_iterator _it = _amount_arr.begin(); _it != _amount_arr.end(); ++_it)
		{
			_item = PyLong_FromSize_t(*_it);
			PyList_Append(PyList_ptr, _item);
			Py_DECREF(_item);
		}
		return PyList_ptr;
	}
	else if (_target == "i")
	{
		std::vector<size_t> const& _amount_arr = self->simu_ptr->get_i_amount_arr();
		for (std::vector<size_t>::const_iterator _it = _amount_arr.begin(); _it != _amount_arr.end(); ++_it)
		{
			_item = PyLong_FromSize_t(*_it);
			PyList_Append(PyList_ptr, _item);
			Py_DECREF(_item);
		}
		return PyList_ptr;
	}
	else if (_target == "r")
	{
		std::vector<size_t> const& _amount_arr = self->simu_ptr->get_r_amount_arr();
		for (std::vector<size_t>::const_iterator _it = _amount_arr.begin(); _it != _amount_arr.end(); ++_it)
		{
			_item = PyLong_FromSize_t(*_it);
			PyList_Append(PyList_ptr, _item);
			Py_DECREF(_item);
		}
		return PyList_ptr;
	}
	else if (_target == "w")
	{
		std::vector<size_t> const& _amount_arr = self->simu_ptr->get_w_amount_arr();
		for (std::vector<size_t>::const_iterator _it = _amount_arr.begin(); _it != _amount_arr.end(); ++_it)
		{
			_item = PyLong_FromSize_t(*_it);
			PyList_Append(PyList_ptr, _item);
			Py_DECREF(_item);
		}
		return PyList_ptr;
	}
	else if (_target == "c")
	{
		std::vector<size_t> const& _amount_arr = self->simu_ptr->get_c_amount_arr();
		for (std::vector<size_t>::const_iterator _it = _amount_arr.begin(); _it != _amount_arr.end(); ++_it)
		{
			_item = PyLong_FromSize_t(*_it);
			PyList_Append(PyList_ptr, _item);
			Py_DECREF(_item);
		}
		return PyList_ptr;
	}
	else if (_target == "d")
	{
		std::vector<size_t> const& _amount_arr = self->simu_ptr->get_d_amount_arr();
		for (std::vector<size_t>::const_iterator _it = _amount_arr.begin(); _it != _amount_arr.end(); ++_it)
		{
			_item = PyLong_FromSize_t(*_it);
			PyList_Append(PyList_ptr, _item);
			Py_DECREF(_item);
		}
		return PyList_ptr;
	}
	else if (_target == "y")
	{
		std::vector<double> const& _amount_arr = self->simu_ptr->get_y_amount_arr();
		for (std::vector<double>::const_iterator _it = _amount_arr.begin(); _it != _amount_arr.end(); ++_it)
		{
			_item = PyFloat_FromDouble(*_it);
			PyList_Append(PyList_ptr, _item);
			Py_DECREF(_item);
		}
		return PyList_ptr;
	}
	else if (_target == "v")
	{
		std::vector<size_t> const& _amount_arr = self->simu_ptr->get_v_amount_arr();
		for (std::vector<size_t>::const_iterator _it = _amount_arr.begin(); _it != _amount_arr.end(); ++_it)
		{
			_item = PyLong_FromSize_t(*_it);
			PyList_Append(PyList_ptr, _item);
			Py_DECREF(_item);
		}
		return PyList_ptr;
	}
	else if (_target == "p")
	{
		std::vector<size_t> const& _amount_arr = self->simu_ptr->get_p_amount_arr();
		for (std::vector<size_t>::const_iterator _it = _amount_arr.begin(); _it != _amount_arr.end(); ++_it)
		{
			_item = PyLong_FromSize_t(*_it);
			PyList_Append(PyList_ptr, _item);
			Py_DECREF(_item);
		}
		return PyList_ptr;
	}
	else
	{
		Py_INCREF(Py_None);
		return Py_None;
	}
}

static PyObject* 
get_x_amount_idx(simuObject* self, PyObject* args, PyObject* kwds)
{
	static char* kwlist[] = { (char*)"idx", (char*)"target", NULL };
	size_t _idx;
	char* _target_ptr;
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "KKs", kwlist, &_idx, &_target_ptr)) {
		return NULL;
	}
	std::string _target(_target_ptr, 1);
	if (_target == "s")
		return PyLong_FromSize_t(self->simu_ptr->get_s_amount_idx(_idx));
	else if (_target == "i")
		return PyLong_FromSize_t(self->simu_ptr->get_i_amount_idx(_idx));
	else if (_target == "r")
		return PyLong_FromSize_t(self->simu_ptr->get_r_amount_idx(_idx));
	else if (_target == "w")
		return PyLong_FromSize_t(self->simu_ptr->get_w_amount_idx(_idx));
	else if (_target == "c")
		return PyLong_FromSize_t(self->simu_ptr->get_c_amount_idx(_idx));
	else if (_target == "d")
		return PyLong_FromSize_t(self->simu_ptr->get_d_amount_idx(_idx));
	else if (_target == "y")
		return PyFloat_FromDouble(self->simu_ptr->get_y_amount_idx(_idx));
	else if (_target == "v")
		return PyLong_FromSize_t(self->simu_ptr->get_v_amount_idx(_idx));
	else if (_target == "p")
		return PyLong_FromSize_t(self->simu_ptr->get_p_amount_idx(_idx));
	else
	{
		Py_INCREF(Py_None);
		return Py_None;
	}
}

static PyObject* 
get_node_state(simuObject* self, PyObject* args, PyObject* kwds)
{
	static char* kwlist[] = { (char*)"group_idx", (char*)"node_idx", NULL };
	size_t _group_idx;
	size_t _node_idx;
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "KK", kwlist, &_group_idx, &_node_idx)) {
		return NULL;
	}
	char _state = self->simu_ptr->get_node_state(_group_idx, _node_idx);
	return PyUnicode_FromStringAndSize(&_state, 1);
}

static PyObject* 
get_node_trans_time(simuObject* self, PyObject* args, PyObject* kwds)
{
	static char* kwlist[] = { (char*)"group_idx", (char*)"node_idx", NULL };
	size_t _group_idx;
	size_t _node_idx;
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "KK", kwlist, &_group_idx, &_node_idx)) {
		return NULL;
	}
	return PyFloat_FromDouble(self->simu_ptr->get_node_trans_time(_group_idx, _node_idx));
}

static PyObject* 
get_node_rem_time(simuObject* self, PyObject* args, PyObject* kwds)
{
	static char* kwlist[] = { (char*)"group_idx", (char*)"node_idx", NULL };
	size_t _group_idx;
	size_t _node_idx;
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "KK", kwlist, &_group_idx, &_node_idx)) {
		return NULL;
	}
	return PyFloat_FromDouble(self->simu_ptr->get_node_rem_time(_group_idx, _node_idx));
}

static PyObject* 
get_node_spreading_params(simuObject* self, PyObject* args, PyObject* kwds)
{
	static char* kwlist[] = { (char*)"group_idx", (char*)"node_idx", NULL };
	size_t _group_idx;
	size_t _node_idx;
	if (!PyArg_ParseTupleAndKeywords(args, kwds, "KK", kwlist, &_group_idx, &_node_idx)) {
		return NULL;
	}
	std::pair<std::vector<double>, std::vector<double>> const& _node_spreading_params_ref = self->simu_ptr->get_node_spreading_params(_group_idx, _node_idx);
	PyObject* PyList_ptr = PyList_New(0);
	PyObject* _PyList_fist_ptr = PyList_New(0);
	PyObject* _PyList_second_ptr = PyList_New(0);
	PyObject* _item = NULL;
	for (std::vector<double>::const_iterator _it = _node_spreading_params_ref.first.begin(); _it != _node_spreading_params_ref.first.end(); ++_it)
	{
		_item = PyFloat_FromDouble(*_it);
		PyList_Append(_PyList_fist_ptr, _item);
		Py_DECREF(_item);
	}
	PyList_Append(PyList_ptr, _PyList_fist_ptr);
	Py_DECREF(_PyList_fist_ptr);
	for (std::vector<double>::const_iterator _it = _node_spreading_params_ref.second.begin(); _it != _node_spreading_params_ref.second.end(); ++_it)
	{
		_item = PyFloat_FromDouble(*_it);
		PyList_Append(_PyList_second_ptr, _item);
		Py_DECREF(_item);
	}
	PyList_Append(PyList_ptr, _PyList_second_ptr);
	Py_DECREF(_PyList_second_ptr);
	return PyList_ptr;
}

static PyMethodDef simu_methods[] = {
	{ "set_vaccination_params", (PyCFunction)set_vaccination_params, METH_VARARGS | METH_KEYWORDS, NULL },
	{ "set_spreading_func", (PyCFunction)set_spreading_func, METH_VARARGS | METH_KEYWORDS, NULL },
	{ "set_total_spreading_params", (PyCFunction)set_total_spreading_params, METH_VARARGS | METH_KEYWORDS, NULL },
	{ "set_detailed_spreading_params", (PyCFunction)set_detailed_spreading_params, METH_VARARGS | METH_KEYWORDS, NULL },
	{ "set_total_spreading_survivals", (PyCFunction)set_total_spreading_survivals, METH_VARARGS | METH_KEYWORDS, NULL },
	{ "init_all_states", (PyCFunction)init_all_states, METH_NOARGS, NULL },
	{ "set_prob_seeds", (PyCFunction)set_prob_seeds, METH_VARARGS | METH_KEYWORDS, NULL },
	{ "set_amount_seeds", (PyCFunction)set_amount_seeds, METH_VARARGS | METH_KEYWORDS, NULL },
	{ "set_state_seeds", (PyCFunction)set_state_seeds, METH_VARARGS | METH_KEYWORDS, NULL },
	{ "set_generator_seed", (PyCFunction)set_generator_seed, METH_VARARGS | METH_KEYWORDS, NULL },
	{ "adjust_rem_time", (PyCFunction)adjust_rem_time, METH_NOARGS, NULL },
	{ "reset_trans_time", (PyCFunction)reset_trans_time, METH_NOARGS, NULL },
	{ "reset_rem_time", (PyCFunction)reset_rem_time, METH_NOARGS, NULL },
	{ "spread_once", (PyCFunction)spread_once, METH_NOARGS, NULL },
	{ "spread", (PyCFunction)spread, METH_VARARGS | METH_KEYWORDS, NULL },
	{ "spread_to_end", (PyCFunction)spread_to_end, METH_NOARGS, NULL },
	{ "add_vaccine", (PyCFunction)add_vaccine, METH_VARARGS | METH_KEYWORDS, NULL },
	{ "get_node_amount", (PyCFunction)get_node_amount, METH_NOARGS, NULL },
	{ "get_group_amount", (PyCFunction)get_group_amount, METH_NOARGS, NULL },
	{ "get_time", (PyCFunction)get_time, METH_NOARGS, NULL },
	{ "get_s_amount", (PyCFunction)get_s_amount, METH_NOARGS, NULL },
	{ "get_i_amount", (PyCFunction)get_i_amount, METH_NOARGS, NULL },
	{ "get_r_amount", (PyCFunction)get_r_amount, METH_NOARGS, NULL },
	{ "get_w_amount", (PyCFunction)get_w_amount, METH_NOARGS, NULL },
	{ "get_c_amount", (PyCFunction)get_c_amount, METH_NOARGS, NULL },
	{ "get_d_amount", (PyCFunction)get_d_amount, METH_NOARGS, NULL },
	{ "get_y_amount", (PyCFunction)get_y_amount, METH_NOARGS, NULL },
	{ "get_v_amount", (PyCFunction)get_v_amount, METH_NOARGS, NULL },
	{ "get_p_amount", (PyCFunction)get_p_amount, METH_NOARGS, NULL },
	{ "get_x_amount", (PyCFunction)get_x_amount, METH_VARARGS | METH_KEYWORDS, NULL },
	{ "get_s_amount_arr", (PyCFunction)get_s_amount_arr, METH_NOARGS, NULL },
	{ "get_i_amount_arr", (PyCFunction)get_i_amount_arr, METH_NOARGS, NULL },
	{ "get_r_amount_arr", (PyCFunction)get_r_amount_arr, METH_NOARGS, NULL },
	{ "get_w_amount_arr", (PyCFunction)get_w_amount_arr, METH_NOARGS, NULL },
	{ "get_c_amount_arr", (PyCFunction)get_c_amount_arr, METH_NOARGS, NULL },
	{ "get_d_amount_arr", (PyCFunction)get_d_amount_arr, METH_NOARGS, NULL },
	{ "get_y_amount_arr", (PyCFunction)get_y_amount_arr, METH_NOARGS, NULL },
	{ "get_v_amount_arr", (PyCFunction)get_v_amount_arr, METH_NOARGS, NULL },
	{ "get_p_amount_arr", (PyCFunction)get_p_amount_arr, METH_NOARGS, NULL },
	{ "get_x_amount_arr", (PyCFunction)get_x_amount_arr, METH_VARARGS | METH_KEYWORDS, NULL },
	{ "get_s_amount_idx", (PyCFunction)get_s_amount_idx, METH_VARARGS | METH_KEYWORDS, NULL },
	{ "get_i_amount_idx", (PyCFunction)get_i_amount_idx, METH_VARARGS | METH_KEYWORDS, NULL },
	{ "get_r_amount_idx", (PyCFunction)get_r_amount_idx, METH_VARARGS | METH_KEYWORDS, NULL },
	{ "get_w_amount_idx", (PyCFunction)get_w_amount_idx, METH_VARARGS | METH_KEYWORDS, NULL },
	{ "get_c_amount_idx", (PyCFunction)get_c_amount_idx, METH_VARARGS | METH_KEYWORDS, NULL },
	{ "get_d_amount_idx", (PyCFunction)get_d_amount_idx, METH_VARARGS | METH_KEYWORDS, NULL },
	{ "get_y_amount_idx", (PyCFunction)get_y_amount_idx, METH_VARARGS | METH_KEYWORDS, NULL },
	{ "get_v_amount_idx", (PyCFunction)get_v_amount_idx, METH_VARARGS | METH_KEYWORDS, NULL },
	{ "get_p_amount_idx", (PyCFunction)get_p_amount_idx, METH_VARARGS | METH_KEYWORDS, NULL },
	{ "get_x_amount_idx", (PyCFunction)get_x_amount_idx, METH_VARARGS | METH_KEYWORDS, NULL },
	{ "get_delay", (PyCFunction)get_delay, METH_NOARGS, NULL },
	{ "get_eta", (PyCFunction)get_eta, METH_NOARGS, NULL },
	{ "get_step", (PyCFunction)get_step, METH_NOARGS, NULL },
	{ "get_populations", (PyCFunction)get_populations, METH_NOARGS, NULL },
	{ "get_population_amounts", (PyCFunction)get_population_amounts, METH_NOARGS, NULL },
	{ "get_ifrs", (PyCFunction)get_ifrs, METH_NOARGS, NULL },
	{ "get_contacts", (PyCFunction)get_contacts, METH_NOARGS, NULL },
	{ "get_x_amount", (PyCFunction)get_x_amount, METH_VARARGS | METH_KEYWORDS, NULL },
	{ "get_x_amount_idx", (PyCFunction)get_x_amount_idx, METH_VARARGS | METH_KEYWORDS, NULL },
	{ "get_node_state", (PyCFunction)get_node_state, METH_VARARGS | METH_KEYWORDS, NULL },
	{ "get_node_trans_time", (PyCFunction)get_node_trans_time, METH_VARARGS | METH_KEYWORDS, NULL },
	{ "get_node_rem_time", (PyCFunction)get_node_rem_time, METH_VARARGS | METH_KEYWORDS, NULL },
	{ "get_node_spreading_params", (PyCFunction)get_node_spreading_params, METH_VARARGS | METH_KEYWORDS, NULL },

	// Terminate the array with an object containing nulls.
	{ NULL }  /* Sentinel */
};

static PyTypeObject simuType = {
	PyVarObject_HEAD_INIT(NULL, 0)
	"metapopulation_simulation.simu",                                                       // tp_name
	sizeof(simuObject),                                                   // tp_basicsize
	0,                                                                      // tp_itemsize
	(destructor)simu_dealloc,                                             // tp_dealloc
	0,                                                                   // tp_vectorcall_offset  NULL
	0,                                                                      // tp_getattr
	0,                                                                      // tp_setattr
	0,                                                                      // tp_as_async
	0,                                                                      // tp_repr
	0,                                                                      // tp_as_number
	0,                                                                      // tp_as_sequence
	0,                                                                      // tp_as_mapping
	0,                                                                      // tp_hash
	0,                                                                      // tp_call
	0,                                                                      // tp_str
	0,                                                                      // tp_getattro
	0,                                                                      // tp_setattro
	0,                                                                      // tp_as_buffer
	Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,                               // tp_flags
	PyDoc_STR("metapopulation simulation objects"),                         // tp_doc
	0,																		// tp_traverse
	0,																		 // tp_clear
	0,																		// tp_richcompare
	0,                                                                      // tp_weaklistoffset
	0,                                                                      // tp_iter
	0,                                                                      // tp_iternext
	simu_methods,															 // tp_methods
	simu_members,															 // tp_members
	0,																		  // tp_getset
	0,                                                                      // tp_base
	0,                                                                     // tp_dict
	0,                                                                     // tp_descr_get
	0,                                                                     // tp_descr_set
	0,                                                                     // tp_dictoffset
	(initproc)simu_init,                                                 // tp_init
	0,                                                                     // tp_alloc
	simu_new,                                                            // tp_new
	0,                                                                     // tp_free
	0,                                                                     // tp_is_gc
	0,                                                                     // tp_bases
	0,                                                                     // tp_mro
	0,                                                                     // tp_cache
	0,                                                                     // tp_subclasses
	0,                                                                     // tp_weaklist
	0,                                                                     // tp_del
	0,                                                                     // tp_version_tag
	0,                                                                     // tp_finalize
};

static PyModuleDef metapopulation_simulation = {
	PyModuleDef_HEAD_INIT,
	"metapopulation_simulation",
	"The module that simulates the process of age-based metapopulation.",
	-1,
};

PyMODINIT_FUNC
PyInit_metapopulation_simulation(void)
{
	PyObject* m;
	if (PyType_Ready(&simuType) < 0)
		return NULL;

	m = PyModule_Create(&metapopulation_simulation);
	if (m == NULL)
		return NULL;

	Py_INCREF(&simuType);
	if (PyModule_AddObject(m, "simu", (PyObject*)&simuType) < 0) {
		Py_DECREF(&simuType);
		Py_DECREF(m);
		return NULL;
	}
	return m;
}