#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
-----------------------------------------------------------
Script Name: <pipelines.py>

Description:
    This script implements  full processing pipelines for computing 
    an instability metric from SMP (Snow Micro Penetrometer) measurements, 
    including its visualization
    It orchestrates the steps of preprocessing, instability modeling, 
    and postprocessing by leveraging the corresponding modules.
    At the end, these pipelines are the modules to be used operationally.

Author: <Bastian Bergfeld>
Email: <bbergfeld@gmx.net>
Date Created: <Thu May  8 08:59:11 2025>
-----------------------------------------------------------
"""
from smp_instability.pre_processor import PreProcessor
from smp_instability.instability_modelling import R2015_point_instability
import smp_instability.post_processor as post_processor

from smp_instability.logging_class import LoggerConfig

#%%

def plot_smp_instability(pnt_file, output_folder):
    
    """
    pipline to attain an instability plot from a pnt file
    
    pnt_file (str): file pat to the input pnt file
    output_folder (str): directory to save resulting plots
    """
    LoggerConfig.setup_logging(log_to_file = True, log_filename=pnt_file[:-4] + ".log")
    profile = PreProcessor.run(pnt_file, save_2_pkl = True)
    model = R2015_point_instability.run(profile, save_2_pkl = True)
    post_processor.compute_logarithmic_sensitivity(model, "S_Reuter2015", "rc_Reuter2015")

    plotter = post_processor.plotter_model(model)
    fig,_,_ = plotter.density_vs_sth(ax2_metric='logarithmic_sensitivity')
    fig.savefig(output_folder + "\\density_vs_sth.png")
    return(True)