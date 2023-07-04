
const usa_contacts = [[10.10043284,  2.59381512,  1.62108103,  3.68725069,  2.01843086,
                    1.6522059 ,  1.10164538,  0.57036862],
                [ 2.59381512, 23.44951842,  3.52117847,  2.90479694,  4.17075415,
                    2.73156772,  0.91623497,  0.80602836],
                [ 1.62108103,  3.52117847, 11.41973385,  4.96882424,  3.96244646,
                    3.43683773,  1.0515248 ,  0.39471326],
                [ 3.68725069,  2.90479694,  4.96882424,  8.77711568,  5.4623625 ,
                    3.59454164,  1.62793598,  0.65586248],
                [ 2.01843086,  4.17075415,  3.96244646,  5.4623625 ,  7.44529254,
                    4.12432196,  1.34042595,  0.96830938],
                [ 1.6522059 ,  2.73156772,  3.43683773,  3.59454164,  4.12432196,
                    5.63607194,  1.6530207 ,  0.75739252],
                [ 1.10164538,  0.91623497,  1.0515248 ,  1.62793598,  1.34042595,
                    1.6530207 ,  2.83338601,  0.90734513],
                [ 0.57036862,  0.80602836,  0.39471326,  0.65586248,  0.96830938,
                    0.75739252,  0.90734513,  1.48919239]]
const usa_population = [0.12000352, 0.12789141, 0.13925591, 0.13494838, 0.12189751, 0.12724997, 0.11627753, 0.11247577]
const init_i = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]



function calc_next_c(current_c, _r_0){
    return math.subtract(1, math.dotMultiply(math.subtract(1, calc_params[2]), math.map(math.dotMultiply(-_r_0/eigen_max,math.multiply(math.dotMultiply(calc_params[0], calc_params[1]), current_c)),math.exp)))
}

function calc_steady_c(_r_0){
    let current_c = calc_params[2];
    let last_c;
    do{
        last_c = current_c;
        current_c = calc_next_c(current_c, _r_0);
    }while(math.sqrt(math.sum(math.map(math.subtract(current_c, last_c), math.square))) > 0.000001);
    return [current_c, math.multiply(current_c, calc_params[1])]
}

function calc_real_r0(_r_0){
    return math.pow(_r_0, 1 / math.pow(t_gen / t_rem, -1.54706508))
}

var file_data = []
var calc_params = [usa_contacts, usa_population, init_i]
var eigen_max = math.max(math.eigs(math.dotMultiply(calc_params[0], calc_params[1])).values)
var r_0 = parseFloat(document.getElementById("r_0").value);
var t_gen = parseFloat(document.getElementById("t_gen").value);
var t_rem = parseFloat(document.getElementById("t_rem").value);

function get_r_0(){r_0 = parseFloat(document.getElementById("r_0").value);}
function get_t_gen(){t_gen = parseFloat(document.getElementById("t_gen").value);}
function get_t_rem(){t_rem = parseFloat(document.getElementById("t_rem").value);}
function get_default_params(){ 
    calc_params = [usa_contacts, usa_population, init_i];
}

function handleFile(event) {
    let file = event.target.files[0];
    let reader = new FileReader();
    file_data = []
    reader.onload = function (e) {
        let data = new Uint8Array(e.target.result);
        let workbook = XLSX.read(data, { type: "array" });
        for(let i = 0; i< 3; i++){
            if (i == 0){
                file_data.push(XLSX.utils.sheet_to_json(workbook.Sheets[workbook.SheetNames[i]], { header: 1 }));
            }
            else{
                file_data.push(XLSX.utils.sheet_to_json(workbook.Sheets[workbook.SheetNames[i]], { header: 1 }).map(row => row[0]));
            }
        }
    };
    reader.readAsArrayBuffer(file);
}

function get_file_data(){
    calc_params = file_data.slice()
}

function get_usa_data(){
    get_default_params();
    if (document.getElementById("identical_input").checked){
        calc_params[2] = []
        let _init = parseFloat(document.getElementById("identical_init_input").value)
        for(let i = 0; i < calc_params[1].length; i++){
            calc_params[2].push(_init);
        }
    }
    else if (document.getElementById("customized_input").checked){
        calc_params[2] = []
        for(let i = 0; i < calc_params[1].length; i++){
            calc_params[2].push(parseFloat(document.getElementById("customized_"+i).value));
        }
    }
    else if (document.getElementById("random_input").checked){
        calc_params[2] = []
        let _init = parseFloat(document.getElementById("random_init_input").value)
        for(let i = 0; i < calc_params[1].length; i++){
            calc_params[2].push(math.random(0, _init));
        }
    }
}








function plot_params() {
    let div_width = window.innerHeight * 0.4
    let trace = {z: calc_params[0], type: 'heatmap', colorscale: 'Viridis'};
    let layout = {
        title: { text: 'Contacts', font: { size: div_width * 0.05 } },
        xaxis: { title: { text: 'Age Group Index', font: { size: div_width * 0.04 } } },
        yaxis: { title: { text: 'Age Group Index', font: { size: div_width * 0.04 } } },
        margin: { t: div_width * 0.1, b: div_width * 0.15, l: div_width * 0.15, r: div_width * 0.15 },
        colorbar: { thickness: div_width * 0.01,},
        width: div_width ,
        height:div_width * 0.75
    };
    Plotly.newPlot('contact', [trace], layout, { responsive: true, displayModeBar: false });
    let bar_names = ['population', 'init_i']
    let bar_titles = ['Age Distribution', 'Initial Infected Fraction']
    let colors = ['rgb(44, 160, 44)', 'rgb(255, 127, 14)']
    for (let i = 1; i < 3; i++){
        div_width = window.innerHeight * 0.4
        trace = {x: range(0, calc_params[i].length, 1), y: calc_params[i], type: 'bar', marker: {color: colors[i - 1]}};
        layout = {
            title: { text: bar_titles[i - 1], font: { size: div_width * 0.05 } },
            xaxis: { title: { text: 'Age Group Index', font: { size: div_width * 0.05 } } },
            yaxis: { title: { text: 'fraction', font: { size: div_width * 0.05 } } },
            margin: { t: div_width * 0.1, b: div_width * 0.1, l: div_width * 0.2, r: div_width * 0.2 },
            width:div_width,
            height:div_width * 0.6
        };
        Plotly.newPlot(bar_names[i - 1], [trace], layout, { responsive: true, displayModeBar: false });
    }
}

function plot_r0(){
    get_r_0()
    get_t_gen();
    get_t_rem();
    let div_width = window.innerHeight * 0.4
    let real_r0 = calc_real_r0(r_0);
    let data = [ { category: 'Estimated', value: r_0},{ category: 'Rectified', value: real_r0},];
    let xValues = data.map(item => item.category);
    let yValues = data.map(item => item.value);
    let trace = {x: xValues, y: yValues, type: 'bar', marker: {color: ['rgb(31, 119, 180)', 'rgb(255, 127, 14)'], bargap: 0.05 }};
    let layout = {
        title: { text: 'Basic Reproduction Number <i>R</i><sub>0</sub>', font: { size: div_width * 0.05} },
        yaxis: { title: { text: '<i>R</i><sub>0</sub>', font: { size: div_width * 0.05 } } },
        margin: { t: div_width * 0.1, b: div_width * 0.1, l: div_width * 0.2, r: div_width * 0.2 },
        width: div_width,
        height: div_width * 0.6
    };
    Plotly.newPlot('r_0_rectified', [trace], layout, { responsive: true, displayModeBar: false });
}

function plot_c(){
    eigen_max = math.max(math.eigs(math.dotMultiply(calc_params[0], calc_params[1])).values);
    get_r_0()
    get_t_gen();
    get_t_rem();
    let div_width = window.innerHeight*0.4
    let _res = calc_steady_c(r_0);
    let _real_res = calc_steady_c(calc_real_r0(r_0));
    let trace1 = { x: range(0, _res[0].length, 1), y: _res[0], name: 'Estimated', type: 'bar'};
    let trace2 = { x: range(0, _real_res[0].length, 1), y: _real_res[0], name: 'Rectified', type: 'bar' };
    let data = [trace1, trace2];
    let layout = {
        barmode: 'group',
        title: { text: 'Group-level Cumulative Infected Fraction', font: { size: div_width * 0.04 } } ,
        xaxis: { title: { text: 'Age Group Index', font: { size: div_width * 0.05 } }},
        yaxis: { title: { text: 'Fraction', font: { size: div_width * 0.05 } } },
        legend: { font: { size: div_width * 0.027 }, x: 0.9, y: 1, orientation: 'v', bgcolor: 'rgba(0, 0, 0, 0)'},
        margin: { t: div_width * 0.1, b: div_width * 0.1, l: div_width * 0.2, r: div_width * 0.2 },
        width: div_width,
        height: div_width * 0.6
    };
    Plotly.newPlot('c_rectified', data, layout, { responsive: true, displayModeBar: false });
    div_width = window.innerHeight * 0.4
    data = [ { category: 'Estimated', value: _res[1]}, { category: 'Rectified', value: _real_res[1]}, ];
    let xValues = data.map(item => item.category);
    let yValues = data.map(item => item.value);
    let trace = { x: xValues, y: yValues, type: 'bar', marker: {color: ['rgb(31, 119, 180)', 'rgb(255, 127, 14)']}};
    layout = {
        title: { text: 'Cumulative Infected Fraction', font: { size: div_width * 0.05 } },
        yaxis: { title: { text: 'Fraction', font: { size: div_width * 0.05 } } },
        margin: { t: div_width * 0.1, b: div_width * 0.1, l: div_width * 0.2, r: div_width * 0.2 },
        width: div_width,
        height: div_width * 0.6};
    Plotly.newPlot('c_tot_rectified', [trace], layout, { responsive: true, displayModeBar: false });

}

get_usa_data()
plot_params()
plot_r0();
plot_c();
