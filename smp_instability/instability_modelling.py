# -*- coding: utf-8 -*-
"""
#%% Routines for instability modelling
@author: Bergfeld Bastian
"""

import os
import concurrent.futures
import numpy as np
import pandas as pd
import logging
import pickle

import weac as weac
from smp_instability.logging_class import error_handling_decorator

    

class R2015_point_instability:
    """
    computes critical cut length and skier stability index from a given profile (provided by the preprocessor)
    the mechanical systems of the PST and "skier on snowpack" can be adjusted.

    Output: model instance (and dataframe of stability-profile is saved as pickle)
    """    
    def __init__(self, prof, skier_stability_params=None, PST_params=None):
        """
        Initializes the class with a layered snow profile.
        
        Parameters:
        - profile: pandas.dataframe, profile created with pre_processor.py
        
        - skier_stability_params: dict, to initialize weac:
            - totallength: float, lateral length of simulated snowpack in mm
            - skierweight: float, skier weight in kg
            - inclination: float, slope angle in degrees, negative clockwise
            - slab_load: bool, True or False matter if the stress, induced by the overlying slab, should be accounted for
        - PST_params: dict, to initialize weac:
            - system: string, slope normal beam ends '-pst', 'pst-' or vertial '-vpst', 'vpst-'; minus is cutting direction
            - totallength: float, lateral length of simulated snowpack in mm
            - inclination: float, slope angle in degrees, negative clockwise
            - max_cracklength: float, maximal crack length for which Energy release is computed in mm
            - num_da: float, number of infeniesimal crack increments within max_cracklength

        """

        if isinstance(prof, pd.DataFrame):  
                self.profile = prof        
        else:
            raise ValueError("Error: Provided profile is not a Dataframe.")
        self.file_path = self.profile.attrs["filepath"]
        self.skier_stability_params = skier_stability_params or {
            "totallength": 4e4, "inclination": -38, "skierweight": 80, "slab_load": True}
        self.PST_params = PST_params or {
            "system": "pst-", "totallength": 4e3, "inclination": -38, "max_cracklength": 2e3, "num_da": 500}
        self.stab = self.profile[["depthTop","thickness"]].copy(deep=True)
        self.stab.attrs = self.profile.attrs.copy()
        subset = {k: self.profile.attrs[k] for k in ['name', 'filepath'] if k in self.profile.attrs}
        subset["num_layers"] = len(self.profile)
        self.stab.attrs.update(subset)

    @error_handling_decorator
    def _get_slab_profile_for_weac(self, wl_id):
        """Loads the snow profile from a pickle file and computes additional parameters."""
        density = self.profile["CR2020_density"].values
        thickness = self.profile["thickness"].values
        slab_profile = np.column_stack((density, thickness))[:wl_id, :].tolist()
        return slab_profile

#---------- S computation --------------------------------    
    @error_handling_decorator
    def _compute_tau(self, weak_layer_id):
        """Computes shear stress along a given weak layer."""
        skier = weac.Layered(system='skier', layers=self._get_slab_profile_for_weac(weak_layer_id))
        seg_skier = skier.calc_segments(L=self.skier_stability_params["totallength"], m=self.skier_stability_params["skierweight"])['nocrack']
        C_skier = skier.assemble_and_solve(phi=self.skier_stability_params["inclination"], **seg_skier)
        xsl_skier, z_skier, xwl_skier = skier.rasterize_solution(C=C_skier, phi=self.skier_stability_params["inclination"], **seg_skier)
        x, tau = skier.get_weaklayer_shearstress(xwl_skier, z_skier, unit="kPa")
        return(x,tau)
    @error_handling_decorator    
    def _compute_max_tau(self, weak_layer_id):
        """Computes max shear stress for a given weak layer without the static stress induce by the slab"""
        x, tau = self._compute_tau(weak_layer_id)
        # Compute maximal shear stress in the weak layer max_tau_skier
        if self.skier_stability_params["slab_load"]: # take the full shear stress in the weak layer
            tau_skier = tau
        else: # take just the additional shear stress induced by the skier
            tau_skier = tau - self.profile["load_above"].iloc[weak_layer_id] / 1e3 * np.sin(np.deg2rad(self.skier_stability_params["inclination"]))
        return max(abs(tau_skier))
    @error_handling_decorator
    def compute_skier_stability_S(self):
        """Computes the stability ratio S_Reuter2015 and adds it to the DataFrame."""
        logging.info(f"   compute_skier_stability_S for {self.profile.attrs["name"]} {self.profile.shape[0]} layers")

        max_tau_skier = self.profile.index.to_series().apply(self._compute_max_tau)
        self.stab["max_shear_stress"] =  max_tau_skier
        self.stab["S_Reuter2015"] = self.profile["JS1999_sigma_macro"] * 1e3 / max_tau_skier

#---------- rc computation --------------------------------    
    @error_handling_decorator
    def _compute_rc_layer(self, weak_layer_id):
        """Computes rc for a given weak layer."""
        pst = weac.Layered(system=self.PST_params["system"], layers=self._get_slab_profile_for_weac(weak_layer_id))
    
        # Initialize outputs and crack lengths
        Gdif = np.zeros([3, self.PST_params["num_da"]])
        da = np.linspace(1e-6, self.PST_params["max_cracklength"], num=self.PST_params["num_da"])
        
        # Loop through crack lengths
        for i, a in enumerate(da):
            # Obtain lists of segment lengths, locations of foundations.
            seg_err = pst.calc_segments(L=self.PST_params["totallength"], a=a)
            # Assemble system and solve for free constants
            C1 = pst.assemble_and_solve(phi=self.PST_params["inclination"], **seg_err['crack']) 
            # Compute differential and incremental energy release rates
            Gdif[:, i] = pst.gdif(C1, self.PST_params["inclination"], **seg_err['crack'])
        
        w_f = self.profile["R2015_wf"][weak_layer_id]
        r_c = da[Gdif[0, :]*1000 > w_f].min()/10
        return(r_c)
        
    @error_handling_decorator
    def compute_critical_cut_length_rc(self):
        """Computes rc for a given weak layerand adds it to the DataFrame."""
        logging.info(f"   compute_critical_cut_length_rc for {self.profile.attrs["name"]} {self.profile.shape[0]} layers")
        self.stab["rc_Reuter2015"] = self.profile.index.to_series().apply(self._compute_rc_layer)

#---------- input and output --------------------------------  
    def save_dataframe(self):
        """Saves a DataFrame as a pickle file, modifying the filename."""
        try:
            new_filename = self.file_path[:-4] + "_instability_df.pkl"
            with open(new_filename, "wb") as f:
                pickle.dump(self.stab, f)
            logging.info(f"   Saved: {new_filename}")
        except: raise AssertionError("Error: saving instability dataframe failed")
            
    def save_instance(self):
        """Saves a class instance as a pickle file, modifying the filename."""
        try: 
            new_filename = self.file_path[:-4] + "_instability_instance.pkl"
            with open(new_filename, "wb") as f:
                pickle.dump(self, f)
            logging.info(f"   Saved: {new_filename}")
        except: raise AssertionError("Error: saving class instance failed")
    
    @classmethod
    def run(cls, profile, skier_stability_params=None, PST_params=None, save_2_pkl = False):
        """Creates an instance and runs all necessary computations."""
        instance = cls(profile, skier_stability_params, PST_params)
        instance.compute_skier_stability_S()
        instance.compute_critical_cut_length_rc()
        if save_2_pkl: 
            instance.save_dataframe()
            instance.save_instance()        
        return instance       

    @classmethod
    def run_multiple(cls, source, skier_stability_params=None, PST_params=None, save_2_pkl = False):
        """Runs instability modelling, saves results, and returns DataFrames.

        Checks if profile is:

        1. pd.DataFrame
        2. list of pd.DataFrames -> return as is
        
        3. .pkl file -> Returns the loaded profile.
        4. list of .pkl files -> Returns list of loaded profiles.
        
        5. folder -> Returns a list of profiles (all .pkl files in the folder).

        --> returns a list of R2015_point_instability instances (beside saving instability profile as pickle)
        """
        logging.info("starting R2015_point_instability...")
        
        def run_df_list(dfs, skier_stability_params, PST_params, save_2_pkl = False):
            # Multiple profiles, parallel processing
            results = []    
            with concurrent.futures.ThreadPoolExecutor() as executor:
                futures = {executor.submit(cls.run, i_prof, skier_stability_params, PST_params, save_2_pkl = False): i_prof for i_prof in dfs}
                for future in concurrent.futures.as_completed(futures):
                    try:
                        df_stab = future.result()  # Get DataFrame result
                        results.append(df_stab)  # Store DataFrame in results list
                    except Exception as e:
                        print(f"Error  {df_stab.file_path}: {e}")
            return(results)
                        
        # Case 1
        if isinstance(source, pd.DataFrame):  
                df_stab = cls.run(source, skier_stability_params, PST_params, save_2_pkl = False)
                stab_profiles = [df_stab]

        # Case 2 & 4
        elif isinstance(source, list): 
            # Case 2 (list of profiles)
            if isinstance(source[0], pd.DataFrame):
                stab_profiles = run_df_list(source, skier_stability_params, PST_params, save_2_pkl = False)
            # Case 4 (list of strings)
            elif isinstance(source[0], str):
                loaded_profiles = []
                for i_source in source:
                    if i_source.lower().endswith('.pkl'):
                        with open(i_source, 'rb') as f:
                            loaded_profiles.append(pickle.load(f)) 
                    else:
                        raise ValueError(f"Error: {i_source} is a file but not a .pkl file.")
                stab_profiles = run_df_list(loaded_profiles, skier_stability_params, PST_params, save_2_pkl = False)
            else: raise ValueError("Error: source is a list but not dataframe or filename")
                
        # Case 3
        elif os.path.isfile(source): 
            if source.lower().endswith('.pkl'):
                with open(source, 'rb') as f:
                    loaded_profile = pickle.load(f)
                df_stab = cls.run(loaded_profile, skier_stability_params, PST_params, save_2_pkl = False)
                stab_profiles = [df_stab]
            else:
                raise ValueError(f"Error: {source} is a file but not a .pkl file.")       

        # Case 5
        elif os.path.isdir(source):
            pkl_profiles = [os.path.join(source, f) for f in os.listdir(source) if f.lower().endswith('_clustered_profile.pkl')]
            if pkl_profiles:
                loaded_profiles = []
                for i_source in pkl_profiles:
                    if i_source.lower().endswith('.pkl'):
                        with open(i_source, 'rb') as f:
                            loaded_profiles.append(pickle.load(f))
                    else:
                        raise ValueError("Error: {i_source} is a file but not a .pkl file.")
                stab_profiles = run_df_list(loaded_profiles, skier_stability_params, PST_params, save_2_pkl = False)

            else: raise ValueError("No _clustered_profile.pkl files found in the folder.")
                
        else: raise ValueError("Not a valid source")
        
        logging.info("Instability modelling finished...")
        return(stab_profiles)  # Return list of DataFrames