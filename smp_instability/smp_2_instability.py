#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
-----------------------------------------------------------
Script Name: smp_2_instability.py

Description:
    This script is meant as an "extended-use" example. Its a compilation of 
    use cases, which indicate how the different
    modules are combined and which input options are availible. 

Author: Bastian Bergfeld
Email: bbergfeld@gmx.net
Date Created: Fri May 2 08:59:48 2025
-----------------------------------------------------------
"""


from smp_instability.pre_processor import PreProcessor
from smp_instability.instability_modelling import R2015_point_instability
import smp_instability.post_processor as post_processor

from smp_instability.logging_class import LoggerConfig

S_params =  {"totallength": 4e4, "inclination": -38, "skierweight": 80, "slab_load": True}
PST_params = {"system": "pst-", "totallength": 4e3, "inclination": 0, "max_cracklength": 2e3, "num_da": 500}


#%%  run a single pnt file

pnt_file = 'D:\\Github_repos\\smp_instability\\smp_instability\\.pnt_data\\SMP1_FILE0003.pnt'
LoggerConfig.setup_logging(log_to_file = True, log_filename=pnt_file[:-4] + ".log")

profile = PreProcessor.run(pnt_file, save_2_pkl = False)

model = R2015_point_instability.run(profile,skier_stability_params=S_params,PST_params= PST_params, save_2_pkl = False)

post_processor.compute_logarithmic_sensitivity(model, "S_Reuter2015", "rc_Reuter2015")

plotter = post_processor.plotter_model(model)
plotter.density_vs_sth()


#%%  run all files in a folder (parallel computing)

root = 'D:\\Github_repos\\smp_instability\\smp_instability\\.pnt_data\\'
LoggerConfig.setup_logging(log_to_file= True, log_filename=root + "all_files.log")

_ = PreProcessor.run_multiple(root, save_2_pkl = True)
_ = R2015_point_instability.run_multiple(root,skier_stability_params=S_params,PST_params= PST_params, save_2_pkl = True)

#%% the cell above saved  model results. 
# Now we load the results an perform postprocessing.

models = post_processor.load_instability_results(root)

for meas in models:
    post_processor.compute_logarithmic_sensitivity(meas, "S_Reuter2015", "rc_Reuter2015")
    post_processor.compute_harmonic_mean(meas, "S_Reuter2015", "rc_Reuter2015")
    post_processor.compute_geometric_mean(meas, "S_Reuter2015", "rc_Reuter2015")
    post_processor.compute_logarithmic_mean(meas, "S_Reuter2015", "rc_Reuter2015")
    post_processor.compute_reziprocal_sum(meas, "S_Reuter2015", "rc_Reuter2015")
     


