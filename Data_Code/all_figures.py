# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 12:50:00 2023

@author: admin
"""

using_imported_data = False



if not using_imported_data:
    import import_experiment_data
    import generate_data
        
    ks = generate_data.get_ks()
        
    #import experiment data
    simu_data = import_experiment_data.import_simu_data()
    calc_data = import_experiment_data.import_calc_data()
    simu_data_m = import_experiment_data.import_simu_data_m()
    simu_data_nm = import_experiment_data.import_simu_data_nm()
    steady_simu_data = import_experiment_data.import_steady_simu_data()
    steady_calc_data = import_experiment_data.import_steady_calc_data()
    steady_simu_mar_data = import_experiment_data.import_steady_simu_mar_data()
    steady_calc_mar_data = import_experiment_data.import_steady_calc_mar_data()
    transient_data = import_experiment_data.import_transient_data()
    
    #generate data
    log_ratios = import_experiment_data.get_log_ratios(transient_data)
    figure_theories_data = generate_data.generate_figure_theories_data(simu_data, calc_data, simu_data_m, simu_data_nm, 
                                                                       steady_simu_data, steady_calc_data, steady_simu_mar_data, steady_calc_mar_data, transient_data)
    
    si_figure_percentile_data = generate_data.generate_si_figure_percentile_data(transient_data)
    si_figure_markovian_time_data = generate_data.generate_si_figure_markovian_time_data()
    si_figure_error_percentile_data = generate_data.generate_si_figure_error_percentile_data(transient_data, log_ratios)
    
    
    #delete imported data
    del simu_data
    del calc_data
    del simu_data_m
    del simu_data_nm
    del steady_simu_data
    del steady_calc_data
    del steady_simu_mar_data
    del steady_calc_mar_data
    del transient_data
    
    ############################################################################################################################################
    
    
    #import experiment data
    r0_g_data = import_experiment_data.import_r0_g_data()
    forecasting_curve_data = import_experiment_data.import_forecasting_curve_data()
    forecasting_data = import_experiment_data.import_forecasting_data()
    real_forecasting_data = import_experiment_data.import_real_forecasting_data()
    explain_data = import_experiment_data.import_explain_data()
    
    
    #generate data
    figure_analysis_data= generate_data.generate_figure_analysis_data(explain_data, forecasting_data, r0_g_data, log_ratios, ks)
    figure_forecasting_data, fitting_period_data= generate_data.generate_figure_forecasting_data(forecasting_curve_data, forecasting_data, real_forecasting_data, log_ratios)
    
    #delete imported data
    del forecasting_curve_data
    del forecasting_data
    del real_forecasting_data
    del explain_data
    
    ############################################################################################################################################
    
    #import experiment data
    vac_curve_data =  import_experiment_data.import_vac_curve_data()
    vac_data =  import_experiment_data.import_vac_data()
    real_vac_data =  import_experiment_data.import_real_vac_data()
    
    #generate data
    figure_prevention_data = generate_data.generate_figure_prevention_data(vac_curve_data, vac_data, real_vac_data, log_ratios)
    
    #delete imported data
    del vac_curve_data
    del vac_data
    del real_vac_data
    
    import pandas as pd
    
    def export_to_file(_data, _file_dir):
        _writer = pd.ExcelWriter(_file_dir)
        for _key in _data.keys():
            pd.DataFrame(_data[_key]).to_excel(_writer, sheet_name = _key)
        _writer.close()
    
    export_to_file(figure_theories_data, './generated_data/theories_data(Figure_2).xlsx')
    export_to_file(si_figure_percentile_data, './generated_data/percentile_data(Figure_S3).xlsx')
    export_to_file(si_figure_markovian_time_data, './generated_data/markovian_time_data(Figure_S1).xlsx')
    export_to_file(si_figure_error_percentile_data, './generated_data/error_percentile_data(Figure_S2).xlsx')
    export_to_file(figure_analysis_data, './generated_data/analysis_data(Figure_4).xlsx')
    export_to_file(fitting_period_data, './generated_data/fitting_period_data(Figure_4_and_5).xlsx')
    export_to_file(figure_forecasting_data, './generated_data/forecasting_data(Figure_5).xlsx')
    export_to_file(figure_prevention_data, './generated_data/prevention_data(Figure_6).xlsx')
    
    import figure_func
    
    figure_func.figure_theories(figure_theories_data)
    
    figure_func.figure_analysis(figure_analysis_data)
    
    figure_func.figure_forecasting(figure_forecasting_data, fitting_period_data)
    
    figure_func.figure_prevention(figure_prevention_data, fitting_period_data)
    
    figure_func.si_figure_percentile(si_figure_percentile_data)
    
    figure_func.si_figure_markovian_time(si_figure_markovian_time_data)
    
    figure_func.si_figure_error_percentile(si_figure_error_percentile_data)
    

else:
    
    def import_from_file(_file_dir):
        _data = pd.read_excel(_file_dir, index_col = 0, sheet_name = None)
        for _key in _data.keys():
            for _key1 in _data[_key].keys():
                _data[_key][_key1] = _data[_key][_key1].to_numpy()
        return _data
    
    figure_theories_data_imported = import_from_file('./generated_data/theories_data(Figure_2).xlsx')
    si_figure_percentile_data_imported = import_from_file('./generated_data/percentile_data(Figure_S3).xlsx')
    si_figure_markovian_time_data_imported = import_from_file('./generated_data/markovian_time_data(Figure_S1).xlsx')
    si_figure_error_percentile_data_imported = import_from_file('./generated_data/error_percentile_data(Figure_S2).xlsx')
    figure_analysis_data_imported = import_from_file('./generated_data/analysis_data(Figure_3).xlsx')
    figure_forecasting_data_imported = import_from_file('./generated_data/forecasting_data(Figure_4).xlsx')
    figure_prevention_data_imported = import_from_file('./generated_data/prevention_data(Figure_5).xlsx')
    
    import figure_func

    figure_func.figure_theories(figure_theories_data_imported)
    
    figure_func.figure_analysis(figure_analysis_data_imported)
    
    figure_func.figure_forecasting(figure_forecasting_data_imported)
    
    figure_func.figure_prevention(figure_prevention_data_imported)
    
    figure_func.si_figure_percentile(si_figure_percentile_data_imported)
    
    figure_func.si_figure_markovian_time(si_figure_markovian_time_data_imported)
    
    figure_func.si_figure_error_percentile(si_figure_error_percentile_data_imported)
    








