#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
-----------------------------------------------------------
Script Name: <Prototyping.py>
Description: <Aims to develop new functionality->
Author: <Bastian Bergfeld>
Email: <bbergfeld@gmx.net>
Date Created: <Fri May 16 08:58:47 2025>
-----------------------------------------------------------
"""

from smp_instability.pre_processor import PreProcessor
from smp_instability.instability_modelling import R2015_point_instability
import smp_instability.post_processor as post_processor

from smp_instability.logging_class import LoggerConfig

#%%

pnt_file = 'D:\\Github_repos\\smp_instability\\prototyping\.pnt_data\\RHOSSA_stable.pnt'

LoggerConfig.setup_logging(log_to_file = True, log_filename=pnt_file[:-4] + ".log")

profile = PreProcessor.run(pnt_file, save_2_pkl = True)
model = R2015_point_instability.run(profile, save_2_pkl = True)
  
post_processor.robust_scaler(model)
post_processor.compute_harmonic_mean(model, "S_Reuter2015_scaled", "rc_Reuter2015_scaled")

#%%


plotter = post_processor.plotter_model(model)
fig,_,_ = plotter.density_vs_sth(ax2_metric='harmonic_mean')
fig.savefig(output_folder + "\\density_and_instability_profile_" + pnt_file.split("\\")[-1][:-4]+ ".png")

#%%

#%%
