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


from smp_instability.pipelines import plot_smp_instability

pnt_file = 'D:\\Github_repos\\smp_instability\\prototyping\.pnt_data\\SMP1_FILE0003.pnt'
output_folder = 'D:\\example_output_folder_for_smp_2_instability'

plot_smp_instability(pnt_file, output_folder)


#%%
