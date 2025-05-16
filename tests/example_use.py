#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
-----------------------------------------------------------
Script Name: example_use.py

Description: 
    This module loads and used the processing pipelines implemented in pipelines.py
    It generates output that can be used in the dashboard.

Author: Bastian Bergfeld
Email: bbergfeld@gmx.net
Date Created: Fri May 2 08:59:48 2025
-----------------------------------------------------------
"""


from smp_instability.pipelines import plot_instability_profile,plot_hazard_map,plot_instability_bar
from smp_instability.logging_class import LoggerConfig

pnt_file = 'D:\\Github_repos\\smp_instability\\prototyping\.pnt_data\\SMP1_FILE0001.pnt'
output_folder = 'D:\\example_output_folder_for_smp_2_instability'

LoggerConfig.setup_logging(log_to_file = True, log_filename=pnt_file[:-4] + ".log")

plot_instability_profile(pnt_file, output_folder, force_computation = False)
plot_hazard_map(pnt_file, output_folder)
plot_instability_bar(pnt_file, output_folder)



