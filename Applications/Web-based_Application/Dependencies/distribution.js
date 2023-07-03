function range(start, end, step) {
    var ans = [];
    for (let i = start; i < end; i++) {
        ans.push(i * step);
    }
    return ans;
}

function range_calc(_range, _calc){
    let res = []
    for (let i = 0; i < _range.length; i++) {
            res.push(_calc(_range[i]));
        }
    return res
    }

function srv_to_dist(_srv, _step){
let res = [];
for (let i = 0; i < _srv.length - 1; i++) {
        res.push((_srv[i] - _srv[i + 1]) / _step);
    }
return res;
}

function srv_to_haz(_srv, _step){
let res = [];
for (let i = 0; i < _srv.length - 1; i++) {
    if (_srv[i] == 0){
        res.push(1 / _step);
    }  else {
        res.push(((_srv[i] - _srv[i + 1]) / _srv[i]) / _step);
    }
    }
return res;
}

function weibull_srv(_step, _len, _alpha, _beta){
let res = [];
let _tau = 0
for (let i = 0; i < _len + 1; i++) {
    _tau = i * _step
    res.push(math.exp(-math.pow(_tau / _beta, _alpha)));
    }
return res;
}

function weibull_dist(_step, _len, _alpha, _beta){
let res = [];
let _tau = 0
for (let i = 0; i < _len; i++) {
    _tau = i * _step
    res.push((_alpha / _beta) * math.pow(_tau /_beta, _alpha - 1) * math.exp(-math.pow(_tau/_beta, _alpha)));
    }
return res;
}

function weibull_mean(_alpha, _beta){
    return _beta * math.gamma(1 + 1 / _alpha);
}


function lognormal_srv(_step, _len, _alpha, _beta){
let res = [];
let _tau = 0
for (let i = 0; i < _len + 1; i++) {
    _tau = i * _step
    res.push(0.5 - 0.5 * math.erf((math.log(_tau) - _alpha)/ (_beta * math.sqrt(2))));
    }
return res;
}

function lognormal_dist(_step, _len, _alpha, _beta){
let res = [];
let _tau = 0
for (let i = 0; i < _len; i++) {
    _tau = i * _step
    res.push((1 / (_tau*_beta * math.sqrt(2 * Math.PI))) * math.exp(-math.square(math.log(_tau) - _alpha)/(2*math.square(_beta))));
    }
return res;
}

function lognormal_mean(_alpha, _beta){
    return math.exp(_alpha + math.square(_beta) / 2);
}

function gamma_dist(_step, _len, _alpha, _beta){
let res = [];
let _tau = 0
for (let i = 0; i < _len; i++) {
    _tau = i * _step
    res.push((1 / (math.gamma(_alpha) * math.pow(_beta, _alpha))) * math.pow(_tau, _alpha - 1) * math.exp(-_tau / _beta));
    }
return res;
}

function incompleteGamma(_s, _x) {
    const _epsilon = 1e-8; // Desired precision
    let _result = 0;
    let _term = 1;
    let _k = 0;

    while (_term > _epsilon) {
      let _numerator = math.pow(_x, _k);
      let _denominator = math.gamma(_s + _k + 1);
      _term = _numerator / _denominator;
      _result += _term;
      _k++;
    }
    return _result * math.exp(-_x) * math.pow(_x, _s);
}


function gamma_srv(_step, _len, _alpha, _beta){
    let res = [1];
    let _tau = 0
    for (let i = 0; i < _len; i++) {
        _tau = i * _step
        res.push(1 - incompleteGamma(_alpha, _tau / _beta));
        }
    return res;
}
  

function gamma_mean(_alpha, _beta){
    return _alpha * _beta;
}


function haz(_dist, _srv){
    let _res = []
    for (let i = 0; i < _dist.length; i++) {
        if (_srv[i] == 0) break;
        _res.push(_dist[i] / _srv[i]);
        }
    return _res;
}


function weibull_srv_func(_tau, _alpha, _beta){
    return math.exp(-math.pow(_tau / _beta, _alpha));
}

function gamma_srv_func(_tau, _alpha, _beta){
    return 1 - incompleteGamma(_alpha, _tau / _beta);
}

function lognormal_srv_func(_tau, _alpha, _beta){
    return 0.5 - 0.5 * math.erf((math.log(_tau) - _alpha)/ (_beta * math.sqrt(2)));
}


function generation_dist(_srv_func_inf, _srv_func_rem, _alpha_inf, _beta_inf, _alpha_rem, _beta_rem, _length, _step){
let _inf_epsilon = 1e-3;
let _tau = 0;
let _tot = 0;
let _res = [];
let _srv_inf_tmp;
let _srv_rem_tmp;
let _mean = 0;
let _i = 0;
while(true){
    _tau = _i * _step;
    _srv_inf_tmp = _srv_func_inf(_tau, _alpha_inf, _beta_inf);
    if (_i < _length){
        if (_srv_inf_tmp == 0) break;}
    else{
        if (_srv_inf_tmp < _inf_epsilon) break;
    }
    _srv_rem_tmp = _srv_func_rem(_tau, _alpha_rem, _beta_rem);
    _res.push((1 - _srv_func_inf(_tau + _step, _alpha_inf, _beta_inf) / _srv_inf_tmp) * _srv_rem_tmp)
    _tot += _res[_i];
    _mean += _res[_i] * _tau;
    ++_i;
}; 
for (let i = 0; i < _res.length; i++) {
    _res[i] /= _tot * _step;
}
_mean /= _tot;
return [_res, _mean, _tot];
}


function findDuration(f, _r_0, _t_gen, _step)
{
    let _left = -1;
    let _right = 1;
    while (f(_left, _r_0, _t_gen, _step) <= 0){
        _left *= 2;
    }
    while (f(_right, _r_0, _t_gen, _step) >= 0){
        _right *= 2;
    }
    return [_left, _right];
}

function findRoot(f, a, b, _r_0, _t_gen, _step, epsilon) {
    let f_a = f(a, _r_0, _t_gen, _step);
    let f_b = f(b, _r_0, _t_gen, _step);
    if (f_a * f_b >= 0) {
      throw new Error("f(a) and f(b) must have opposite signs.");
    }
    let c = a;
    let f_c = 0
    while ((b - a) >= epsilon) {
      c = (a + b) / 2;
      f_c = f(c, _r_0, _t_gen, _step);
      if (f_c === 0) {
        break;
      } else if (f_c < 0) {
        b = c;
      } else {
        a = c;
      }
    }
    return c;
  }

function growth_rate_func(_g, _r_0, _t_gen, _step){
    let _res = 0;
    let _tau = 0;
    for (let _i = 0; _i < _t_gen.length; ++_i){
        _tau = _i * _step;
        _res += math.exp(-_g*_tau) * _t_gen[_i] * _step;
    }
    return _res * _r_0 - 1;
}






function line_plot(_xArray, _yArray, _mean, xlabel, ylabel, plot_title, _div_id, line_color, _time_length) {
    let div_width = window.innerHeight * 0.4;
    let data = [{ x:_xArray, y:_yArray, mode:"lines", line: {color: line_color}}];
    let x_len = _time_length;
    let y_len = math.max(_yArray);
    let layout = {
        xaxis: { range: [- x_len / 25, x_len + x_len / 25], title: { text: xlabel, font: { size: div_width * 0.05 } } },
        yaxis: { range: [0, math.max(_yArray) + y_len / 10], title: { text: ylabel, font: { size: div_width * 0.05 } }},  
        title: { text: plot_title, font: { size: div_width * 0.05 } },
        margin: { t: div_width * 0.1, b: div_width * 0.1, l: div_width * 0.2, r: div_width * 0.2 },
        width: div_width,
        height: div_width * 0.8
    };
    Plotly.newPlot(_div_id, data, layout, { responsive: true, displayModeBar: false });
    if(_mean != false){
    update_layout = {shapes: [{type: 'line',xref: 'x',yref: 'paper',x0: _mean, y0: 0, x1: _mean, y1: 1,
        line: { color: 'gray', dash: 'dash', width: 1, },},],};
    Plotly.update(_div_id, {}, update_layout);
    }
}

function bar_plot(_val_gamma, _val_mu, _div_id) {
    let div_width = window.innerHeight * 0.4;
    let data = [ { category: '$\\gamma$', value: _val_gamma},{ category: '$\\mu$', value: _val_mu},];
    let _xValues = data.map(item => item.category);
    let _yValues = data.map(item => item.value);
    let trace = {x: _xValues, y: _yValues, type: 'bar', marker: {color: ['rgb(31, 119, 180)', 'rgb(255, 127, 14)'], bargap: 0.05 }};
    let layout = {
        title: { text: 'Markovian Parameters<br>in Transient-state Equivalence', font: { size: div_width * 0.05 } },
        xaxis: { title: { font: { size: div_width * 0.05 } } },
        yaxis: { title: { text: "Value", font: { size: div_width * 0.05 } } },
        margin: { t: div_width * 0.1, b: div_width * 0.1, l: div_width * 0.2, r: div_width * 0.2 },
        width: div_width,
        height: div_width * 0.8
    };
    Plotly.newPlot(_div_id, [trace], layout, { responsive: true, displayModeBar: false });
}


//const default_params = {'Weibull': ["1", "1"], 'Gamma': ["1", "1"], 'Log-normal': ["0", "0.1"]}
const default_ranges = {'Weibull': [["0.4", "4"], ["0.4", "4"]], 'Gamma': [["0.4", "4"], ["0.4", "4"]], 'Log-normal': [["-2", "2"], ["0.1", "1"]]}
const default_range_steps = {'Weibull': ["0.02", "0.02"], 'Gamma': ["0.02", "0.02"], 'Log-normal': ["0.02", "0.01"]}
const default_param_names = {'Weibull': ["alpha", "beta"], 
                            'Gamma': ["alpha", "beta"], 
                            'Log-normal': ["mu", "theta"]}


function update_sliders(process_name){  
    let func_name = document.getElementById(process_name + "_opt").value
    for (let j = 1; j < 3; j++){
        let p_b = document.getElementById(process_name + "_p"+ j +"_slider");
        let p_ele = document.getElementById("p_" + process_name +"_" + j);
        p_b.style.setProperty('--min',default_ranges[func_name][j-1][0]);
        p_b.style.setProperty('--max',default_ranges[func_name][j-1][1]);
        p_b.style.setProperty('--step',default_range_steps[func_name][j-1]);
        p_ele.setAttribute("min", default_ranges[func_name][j-1][0]); 
        p_ele.setAttribute("max", default_ranges[func_name][j-1][1]); 
        p_ele.setAttribute("step", default_range_steps[func_name][j-1]); 
        
        let _val = 0;
        if(p_ele.value > p_ele.max) _val = p_ele.max;
        else if (p_ele.value < p_ele.min) _val = p_ele.min;
        else _val = p_ele.value;
        p_ele.setAttribute("value", _val); 
        p_b.style.setProperty('--value',_val); 
        p_b.style.setProperty('--text-value', JSON.stringify((+_val).toLocaleString()));
    }   
}

update_sliders("inf");
update_sliders("rem");

var length = parseInt(document.getElementById("duration").value)
var step = parseFloat(document.getElementById("step").value)
var p_inf_1 = parseFloat(document.getElementById("p_inf_1").value)
var p_inf_2 = parseFloat(document.getElementById("p_inf_2").value)
var p_rem_1 = parseFloat(document.getElementById("p_rem_1").value)
var p_rem_2 = parseFloat(document.getElementById("p_rem_2").value)
var lam_max = parseFloat(document.getElementById("lam_max").value)

const func_forms = {'Weibull': [weibull_srv, weibull_dist], 'Gamma': [gamma_srv, gamma_dist], 'Log-normal': [lognormal_srv, lognormal_dist]}
const srv_funcs = {'Weibull': weibull_srv_func, 'Gamma': gamma_srv_func, 'Log-normal': lognormal_srv_func}
const mean_funcs = {'Weibull': weibull_mean, 'Gamma': gamma_mean, 'Log-normal': lognormal_mean}
var inf_func_form = func_forms[document.getElementById("inf_opt").value]
var rem_func_form = func_forms[document.getElementById("rem_opt").value]
var srv_func_inf = srv_funcs[document.getElementById("inf_opt").value]
var srv_func_rem = srv_funcs[document.getElementById("rem_opt").value]
var mean_func_inf = mean_funcs[document.getElementById("inf_opt").value]
var mean_func_rem = mean_funcs[document.getElementById("rem_opt").value]

var srv_inf = inf_func_form[0](step, length, p_inf_1, p_inf_2);
var dist_inf = inf_func_form[1](step, length, p_inf_1, p_inf_2);
var haz_inf = haz(dist_inf, srv_inf);
var mean_inf = mean_func_inf(p_inf_1, p_inf_2);
var tau_inf = range(0, haz_inf.length, step);
var srv_rem = rem_func_form[0](step, length, p_rem_1, p_rem_2);
var dist_rem = rem_func_form[1](step, length, p_rem_1, p_rem_2);
var haz_rem = haz(dist_rem, srv_rem);
var mean_inf = mean_func_rem(p_rem_1, p_rem_2);
var tau_rem = range(0, haz_rem.length, step);
var res_gen = generation_dist(srv_func_inf, srv_func_rem, p_inf_1, p_inf_2, p_rem_1, p_rem_2, length, step);
var tau_gen = range(0, res_gen[0].length, step);
var r_0 = res_gen[2] * lam_max;
var g_duration = findDuration(growth_rate_func, r_0, res_gen[0], step);
var growth_rate = findRoot(growth_rate_func, g_duration[0], g_duration[1], r_0, res_gen[0], step, 1e-6)
var val_gamma = growth_rate * r_0 / (lam_max * (r_0 - 1));
var val_mu = growth_rate / (r_0 - 1)

function change_calc(){
    document.getElementById("val_r_0").innerHTML="Basic Reproduction Number <i>R</i><sub>0</sub>: " + r_0.toFixed(2);
    document.getElementById("val_g").innerHTML="Growth Rate <i>g</i>: " + growth_rate.toFixed(2);
    bar_plot(val_gamma, val_mu, "markovian")
}


function plot_all_fig(){
    line_plot(tau_inf, srv_inf.slice(0,tau_inf.length), false, "$\\tau$", '$\\Psi_{\\mathrm{inf}}(\\tau)$', "$\\text{Infection Survival Function } \\Psi_{\\mathrm{inf}}(\\tau)$", "srv_inf", 'rgb(214,39,40)', length * step)
    line_plot(tau_inf, dist_inf.slice(0,tau_inf.length), mean_inf, "$\\tau$", '$\\psi_{\\mathrm{inf}}(\\tau)$', "$\\text{Infection Time Distribution } \\psi_{\\mathrm{inf}}(\\tau)$", "psi_inf", 'rgb(214,39,40)', length * step)
    line_plot(tau_inf, haz_inf, false, "$\\tau$", '$\\omega_{\\mathrm{inf}}(\\tau)$', "$\\text{Infection Hazard Function } \\omega_{\\mathrm{inf}}(\\tau)$", "haz_inf", 'rgb(214,39,40)', length * step)
    line_plot(tau_rem, srv_rem.slice(0,tau_rem.length), false, "$\\tau$", '$\\Psi_{\\mathrm{rem}}(\\tau)$', "$\\text{Removal Survival Function } \\Psi_{\\mathrm{rem}}(\\tau)$", "srv_rem", 'rgb(31, 119, 180)', length * step)
    line_plot(tau_rem, dist_rem.slice(0,tau_rem.length), mean_rem, "$\\tau$", '$\\psi_{\\mathrm{rem}}(\\tau)$', "$\\text{Removal Time Distribution } \\psi_{\\mathrm{rem}}(\\tau)$", "psi_rem", 'rgb(31, 119, 180)', length * step)
    line_plot(tau_rem, haz_rem, false, "$\\tau$", '$\\omega_{\\mathrm{rem}}(\\tau)$', "$\\text{Removal Hazard Function } \\omega_{\\mathrm{rem}}(\\tau)$", "haz_rem", 'rgb(31, 119, 180)', length * step)
    line_plot(tau_gen, res_gen[0], res_gen[1], "$\\tau$", '$\\psi_{\\mathrm{gen}}(\\tau)$', "$\\text{Generation Time Distribution } \\psi_{\\mathrm{gen}}(\\tau)$", "psi_gen", 'rgb(255, 127, 14)', length * step)
    bar_plot(val_gamma, val_mu, "markovian")
    document.getElementById("val_t").innerHTML="<i>T</i><sub>inf</sub>: " + mean_inf.toFixed(2) + ", " + "<i>T</i><sub>rem</sub>: " + mean_rem.toFixed(2) + ", " + "<i>T</i><sub>gen</sub>: " + res_gen[1].toFixed(2);
    change_calc();
}

function plot_inf_fig(){
    line_plot(tau_inf, srv_inf.slice(0,tau_inf.length), false, "$\\tau$", '$\\Psi_{\\mathrm{inf}}(\\tau)$', "$\\text{Infection Survival Function } \\Psi_{\\mathrm{inf}}(\\tau)$", "srv_inf", 'rgb(214,39,40)', length * step)
    line_plot(tau_inf, dist_inf.slice(0,tau_inf.length), mean_inf, "$\\tau$", '$\\psi_{\\mathrm{inf}}(\\tau)$', "$\\text{Infection Time Distribution } \\psi_{\\mathrm{inf}}(\\tau)$", "psi_inf", 'rgb(214,39,40)', length * step)
    line_plot(tau_inf, haz_inf, false, "$\\tau$", '$\\omega_{\\mathrm{inf}}(\\tau)$', "$\\text{Infection Hazard Function } \\omega_{\\mathrm{inf}}(\\tau)$", "haz_inf", 'rgb(214,39,40)', length * step)
    line_plot(tau_gen, res_gen[0], res_gen[1], "$\\tau$", '$\\psi_{\\mathrm{gen}}(\\tau)$', "$\\text{Generation Time Distribution } \\psi_{\\mathrm{gen}}(\\tau)$", "psi_gen", 'rgb(255, 127, 14)', length * step)
    bar_plot(val_gamma, val_mu, "markovian")
    document.getElementById("val_t").innerHTML="<i>T</i><sub>inf</sub>: " + mean_inf.toFixed(2) + ", " + "<i>T</i><sub>rem</sub>: " + mean_rem.toFixed(2) + ", " + "<i>T</i><sub>gen</sub>: " + res_gen[1].toFixed(2);
    change_calc();
}

function plot_rem_fig(){
    line_plot(tau_rem, srv_rem.slice(0,tau_rem.length), false, "$\\tau$", '$\\Psi_{\\mathrm{rem}}(\\tau)$', "$\\text{Removal Survival Function } \\Psi_{\\mathrm{rem}}(\\tau)$", "srv_rem", 'rgb(31, 119, 180)', length * step)
    line_plot(tau_rem, dist_rem.slice(0,tau_rem.length), mean_rem, "$\\tau$", '$\\psi_{\\mathrm{rem}}(\\tau)$', "$\\text{Removal Time Distribution } \\psi_{\\mathrm{rem}}(\\tau)$", "psi_rem", 'rgb(31, 119, 180)', length * step)
    line_plot(tau_rem, haz_rem, false, "$\\tau$", '$\\omega_{\\mathrm{rem}}(\\tau)$', "$\\text{Removal Hazard Function } \\omega_{\\mathrm{rem}}(\\tau)$", "haz_rem", 'rgb(31, 119, 180)', length * step)
    line_plot(tau_gen, res_gen[0], res_gen[1], "$\\tau$", '$\\psi_{\\mathrm{gen}}(\\tau)$', "$\\text{Generation Time Distribution } \\psi_{\\mathrm{gen}}(\\tau)$", "psi_gen", 'rgb(255, 127, 14)', length * step)
    bar_plot(val_gamma, val_mu, "markovian")
    document.getElementById("val_t").innerHTML="<i>T</i><sub>inf</sub>: " + mean_inf.toFixed(2) + ", " + "<i>T</i><sub>rem</sub>: " + mean_rem.toFixed(2) + ", " + "<i>T</i><sub>gen</sub>: " + res_gen[1].toFixed(2);
    change_calc();
}


function update_all_fig(){
    length = parseInt(document.getElementById("duration").value)
    step = parseFloat(document.getElementById("step").value)
    xArray = range(0, length, step);
    inf_func_form = func_forms[document.getElementById("inf_opt").value]
    rem_func_form = func_forms[document.getElementById("rem_opt").value]
    srv_func_inf = srv_funcs[document.getElementById("inf_opt").value]
    srv_func_rem = srv_funcs[document.getElementById("rem_opt").value]
    mean_func_inf = mean_funcs[document.getElementById("inf_opt").value]
    mean_func_rem = mean_funcs[document.getElementById("rem_opt").value]
    lam_max = parseFloat(document.getElementById("lam_max").value)
    srv_inf = inf_func_form[0](step, length, p_inf_1, p_inf_2);
    dist_inf = inf_func_form[1](step, length, p_inf_1, p_inf_2);
    haz_inf = haz(dist_inf, srv_inf);
    mean_inf = mean_func_inf(p_inf_1, p_inf_2);
    tau_inf = range(0, haz_inf.length, step);
    srv_rem = rem_func_form[0](step, length, p_rem_1, p_rem_2);
    dist_rem = rem_func_form[1](step, length, p_rem_1, p_rem_2);
    haz_rem = haz(dist_rem, srv_rem);
    mean_rem = mean_func_rem(p_rem_1, p_rem_2);
    tau_rem = range(0, haz_rem.length, step);
    res_gen = generation_dist(srv_func_inf, srv_func_rem, p_inf_1, p_inf_2, p_rem_1, p_rem_2, length, step);
    tau_gen = range(0, res_gen[0].length, step);
    plot_all_fig()
    update_lam_max()
}

function update_inf_fig(){
    p_inf_1 = parseFloat(document.getElementById("p_inf_1").value)
    p_inf_2 = parseFloat(document.getElementById("p_inf_2").value)
    inf_func_form = func_forms[document.getElementById("inf_opt").value]
    srv_func_inf = srv_funcs[document.getElementById("inf_opt").value]
    mean_func_inf = mean_funcs[document.getElementById("inf_opt").value]
    srv_inf = inf_func_form[0](step, length, p_inf_1, p_inf_2);
    dist_inf = inf_func_form[1](step, length, p_inf_1, p_inf_2);
    haz_inf = haz(dist_inf, srv_inf);
    mean_inf = mean_func_inf(p_inf_1, p_inf_2);
    tau_inf = range(0, haz_inf.length, step);
    res_gen = generation_dist(srv_func_inf, srv_func_rem, p_inf_1, p_inf_2, p_rem_1, p_rem_2, length, step);
    tau_gen = range(0, res_gen[0].length, step);
    plot_inf_fig()
    update_lam_max()
}

function update_rem_fig(){
    p_rem_1 = parseFloat(document.getElementById("p_rem_1").value)
    p_rem_2 = parseFloat(document.getElementById("p_rem_2").value)
    rem_func_form = func_forms[document.getElementById("rem_opt").value]
    srv_func_rem = srv_funcs[document.getElementById("rem_opt").value]
    mean_func_rem = mean_funcs[document.getElementById("rem_opt").value]
    srv_rem = rem_func_form[0](step, length, p_rem_1, p_rem_2);
    dist_rem = rem_func_form[1](step, length, p_rem_1, p_rem_2);
    haz_rem = haz(dist_rem, srv_rem);
    mean_rem = mean_func_rem(p_rem_1, p_rem_2);
    tau_rem = range(0, haz_rem.length, step);
    res_gen = generation_dist(srv_func_inf, srv_func_rem, p_inf_1, p_inf_2, p_rem_1, p_rem_2, length, step);
    tau_gen = range(0, res_gen[0].length, step);
    plot_rem_fig();
    update_lam_max()
}

function update_lam_max(){
    lam_max = parseFloat(document.getElementById("lam_max").value)
    r_0 = res_gen[2] * lam_max;
    g_duration = findDuration(growth_rate_func, r_0, res_gen[0], step);
    growth_rate = findRoot(growth_rate_func, g_duration[0], g_duration[1], r_0, res_gen[0], step, 1e-6)
    val_gamma = growth_rate * r_0 / (lam_max * (r_0 - 1));
    val_mu = growth_rate / (r_0 - 1)
    change_calc();
}

function update_inf_func(){
    update_sliders("inf");
    update_inf_fig();
}

function update_rem_func(){
    update_sliders("rem");
    update_rem_fig();
}
update_all_fig()

