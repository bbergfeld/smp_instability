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
import os
import pickle
from smp_instability.pre_processor import PreProcessor
from smp_instability.instability_modelling import R2015_point_instability
import smp_instability.post_processor as post_processor

#%%
    
def load_or_compute_instability_model(pnt_file, force_computation = False):
    """
    Loads model-object from pickle if available, else computes it from pnt_file.

    Parameters:
        pnt_file (str): Path to the pnt file.
        compute_fn (callable): Function to compute the data if pickle is missing or invalid.
        *args, **kwargs: Arguments passed to compute_fn.

    Returns:
        object: The loaded or newly computed data.
    """
    pkl_path = pnt_file[:-4] + "_instability_instance.pkl"
    if force_computation == False:
        try:
            if os.path.isfile(pkl_path):
                with open(pkl_path, 'rb') as f:
                    return pickle.load(f)
        except Exception as e:
            print(f"Failed to load pickle: {e}")

    # If loading failed, compute and save
    profile = PreProcessor.run(pnt_file, save_2_pkl = True)
    model = R2015_point_instability.run(profile, save_2_pkl = True)
    return model
    
    
def plot_instability_profile(pnt_file, output_folder, force_computation = False):
    
    """
    pipline to attain an instability plot from a pnt file
    
    pnt_file (str): file path to the input pnt file
    output_folder (str): directory to save resulting plots
    force_computation (bool): if False, try to load pickled model, if True forcing a (re)computation

    """    
    model = load_or_compute_instability_model(pnt_file, force_computation)
   
    post_processor.robust_scaler(model)
    post_processor.compute_harmonic_mean(model, "S_Reuter2015_scaled", "rc_Reuter2015_scaled")
    
    plotter = post_processor.plotter_model(model)
    fig,_,_ = plotter.density_vs_sth(ax2_metric='harmonic_mean')
    fig.savefig(output_folder + "\\density_and_instability_profile_" + pnt_file.split("\\")[-1][:-4]+ ".png")
    return(True)


def plot_hazard_map(pnt_file, output_folder, force_computation = False):
    
    """
    pipline to attain a hazard plot from a pnt file
    
    pnt_file (str): file path to the input pnt file
    output_folder (str): directory to save resulting plots
    force_computation (bool): if False, try to load pickled model, if True forcing a (re)computation
    """
    model = load_or_compute_instability_model(pnt_file, force_computation)

    post_processor.robust_scaler(model)
    post_processor.compute_harmonic_mean(model, "S_Reuter2015_scaled", "rc_Reuter2015_scaled")

    plotter = post_processor.plotter_model(model)
    fig,_ = plotter.hazard_map()
    fig.savefig(output_folder + "\\hazard_map_" + pnt_file.split("\\")[-1][:-4]+ ".png")
    return(True)


def plot_instability_bar(pnt_file, output_folder, force_computation = False):
    
    """
    pipline to attain an instability bar from a pnt file
    
    pnt_file (str): file path to the input pnt file
    output_folder (str): directory to save resulting plots
    force_computation (bool): if False, try to load pickled model, if True forcing a (re)computation

    """
    model = load_or_compute_instability_model(pnt_file, force_computation)
    post_processor.robust_scaler(model)
    post_processor.compute_harmonic_mean(model, "S_Reuter2015_scaled", "rc_Reuter2015_scaled")

    plotter = post_processor.plotter_model(model)
    fig,_ = plotter.instability_bar()
    fig.savefig(output_folder + "\\instability_bar" + pnt_file.split("\\")[-1][:-4]+ ".png")
    return(True)