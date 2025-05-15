
# -*- coding: utf-8 -*-
"""
Visualization and Metric Combination for Snowpack Instability Analysis
----------------------------------------------------------------------

Provides routines for:
- Loading and combining Reuter2015 snowpack instability model results
- Calculating composite metrics (e.g., harmonic, geometric, logarithmic means)
- Visualizing skier-induced stresses, density profiles, and depth-based instability metrics

Key Features:
- Supports batch loading of instability instances from files or directories
- Combines metrics into new indices for sensitivity analysis
- Plotting utilities for stress fields, dual-metric comparisons, and summary plots across profiles
- Encapsulated plotting tools via `plotter_model` class

Dependencies:
- numpy, pandas, matplotlib, pickle, logging
- Custom `weac` package for layered snowpack modeling and stress computation

Author: Bastian Bergfeld
"""
import logging
import pickle
import os
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import weac as weac


def load_instability_results(source):
    """Returns the loaded instability modelling instances as a list.
    Checks if sourse is:
    1. .pkl file -> Returns the loaded instability instance.
    2. list of .pkl files -> Returns list of instability instance.
    3. folder -> Returns a list of instability instance (all .pkl files in the folder).
    --> returns a list of R2015_point_instability instances
    """
    # Case 2
    if isinstance(source, list): 
        if isinstance(source[0], str):
            loaded_inst_obj = []
            for i_source in source:
                if i_source.lower().endswith('_instability_instance.pkl'):
                    with open(i_source, 'rb') as f:
                        loaded_inst_obj.append(pickle.load(f)) 
                else:
                    raise ValueError(f"Error: {i_source} is a file but not a .pkl file.")
        else: raise ValueError("Error: source is a list but not a correctfilename")
            
    # Case 1
    elif os.path.isfile(source): 
        if source.lower().endswith('_instability_instance.pkl'):
            with open(source, 'rb') as f:
                loaded_inst_obj = [pickle.load(f)]
        else:
            raise ValueError(f"Error: {source} is a file but not a .pkl file.")       

    # Case 3
    elif os.path.isdir(source):
        pkl_objs = [os.path.join(source, f) for f in os.listdir(source) if f.lower().endswith('_instability_instance.pkl')]
        if pkl_objs:
            loaded_inst_obj = []
            for i_source in pkl_objs:
                if i_source.lower().endswith('_instability_instance.pkl'):
                    with open(i_source, 'rb') as f:
                        loaded_inst_obj.append(pickle.load(f))
                else:
                    raise ValueError("Error: {i_source} is a file but not a .pkl file.")
        else: raise ValueError("No _clustered_profile.pkl files found in the folder.")
    else: raise ValueError("Not a valid source")
    
    logging.info("Instability objects loaded...")
    return(loaded_inst_obj)  # Return list of DataFrames






def compute_harmonic_mean(instability_instance, metric_1, metric_2):
    I1,I2 = instability_instance.stab[metric_1], instability_instance.stab[metric_2]
    instability_instance.stab["harmonic_mean"] = (2*I1*I2) / (I1+I2)
def compute_logarithmic_mean(instability_instance, metric_1, metric_2):
    I1,I2 = instability_instance.stab[metric_1], instability_instance.stab[metric_2]
    instability_instance.stab["logarithmic_mean"] = (I1-I2) / (np.log(I1)-np.log(I2))
def compute_geometric_mean(instability_instance, metric_1, metric_2):
    I1,I2 = instability_instance.stab[metric_1], instability_instance.stab[metric_2]
    instability_instance.stab["geometric_mean"] = np.sqrt(I1*I2)
def compute_reziprocal_sum(instability_instance, metric_1, metric_2):
    I1,I2 = instability_instance.stab[metric_1], instability_instance.stab[metric_2]
    instability_instance.stab["reziprocal_sum"] = 1/(1/I1+1/I2)
def compute_logarithmic_sensitivity(instability_instance, metric_1, metric_2):
    I1,I2 = instability_instance.stab[metric_1], instability_instance.stab[metric_2]
    instability_instance.stab["logarithmic_sensitivity"] = np.exp(np.log(I1)+np.log(I2))
def compute_euclidian_distance(instability_instance, metric_1, metric_2):
    I1,I2 = instability_instance.stab[metric_1], instability_instance.stab[metric_2]
    instability_instance.stab["euclidian_distance"] = np.sqrt(I1**2+I2**2)
def compute_weighted_geometric_mean(instability_instance, metric_1, metric_2):
    w1,w2 = 0.8, 0.2
    I1,I2 = instability_instance.stab[metric_1], instability_instance.stab[metric_2]
    scale_factor = np.sqrt(I1.mean() * I2.mean())
    instability_instance.stab["weighted_geometric_mean"] = (I1**w1) * (I2**w2) / scale_factor
    
def scale_metric(instability_instance, method='mean', postfix='_scaled', exclude_cols=['depthTop', 'thickness']):
    stats = {'mean': instability_instance.stab.mean,
             'median': instability_instance.stab.median,
             'min': instability_instance.stab.min,
             'max': instability_instance.stab.max }

    if method not in stats:
        raise ValueError(f"Unsupported method: {method}")

    for col in instability_instance.stab.columns.difference(exclude_cols):
        scale = stats[method]()[col]
        instability_instance.stab[col + postfix] = instability_instance.stab[col] / scale if scale != 0 else float('nan')
        

def robust_scaler(instability_instance, lower_quantile=0.1, upper_quantile=0.6, clip_quantile=0.9):
    """
    Skaliert alle metriken robust (aber ohne Zentrierung),
    sodass keine negativen Werte entstehen. Hohe Ausreißer werden geclippt.
    
    Args:
        model_instance.stab (pd.DataFrame): Der Input-DataFrame.
        lower_quantile (float): Unteres Quantil (z. B. 0.1).
        upper_quantile (float): Oberes Quantil (z. B. 0.9).
        clip_quantile (float): Clipping-Schwelle für hohe Ausreißer.
    
    Returns:
        adds scaled columns to the stab dataframe
    """
    df = instability_instance.stab
    columns_to_scale = [col for col in df.columns if col not in ['depthTop', 'Thickness']]

    for col in columns_to_scale:
        x = df[col].values

        # IQR ohne Zentrierung (kein Abzug von Median)
        q_low = np.quantile(x, lower_quantile)
        q_high = np.quantile(x, upper_quantile)
        iqr = q_high - q_low
        iqr = iqr if iqr != 0 else 1.0

        # Skalieren ohne Median-Abzug
        x_scaled = x / iqr

        # Clipping (falls gewünscht)
        upper_clip = np.quantile(x_scaled, clip_quantile)
        x_scaled = np.clip(x_scaled, None, upper_clip)

        df[f"{col}_scaled"] = x_scaled

    return df

def plot_averaged_from_profiles(models, metrics, window_size=2, overlap=50, how_to_average = "median", average_window = 10,save_fig=False,  log_x = None):
    """ models: list of instability instances
    metrics: list of metrics to plot
     ... various keywords....
    """
    if log_x == None:
         log_x = [False] * len(metrics)
    df_list = []
    for meas in models:
        df = meas.stab.merge(meas.profile)
        layering_thickness = window_size * overlap/100    
        new_depths = np.arange(0, df['thickness'].iloc[-1] + df['depthTop'].iloc[-1] , layering_thickness)
        df_resampled = pd.DataFrame({'depthTop': new_depths})
        df_resampled = df_resampled.merge(df, on='depthTop', how='left').ffill()
        df_list.append(df_resampled)
    df = pd.concat(df_list, ignore_index=True)
    df = df.sort_values(by="depthTop", ascending=True)
    fig, axes = plt.subplots(1, len(metrics), figsize=(len(metrics) * 4, 6), sharey=True)
    if isinstance(axes, np.ndarray):
        pass
    else: axes = [axes]
    for i, metric in enumerate(metrics):
        if how_to_average == "mean":
            df["temp_moving_avg"] = df[metric].rolling(window=average_window, center=True).mean()
        elif how_to_average == "median":
            df["temp_moving_avg"] = df[metric].rolling(window=average_window, center=True).median()
        else: print("use mean or median for how_to_average")
        axes[i].plot(df[metric], df["depthTop"], '-', alpha=0.5, label="range of "+str(len(models))+" measurments")  # Scatter + line
        axes[i].plot(df["temp_moving_avg"], df["depthTop"], 'r-', linewidth=2, label="Moving "+how_to_average+" ("+str(average_window)+"mm)")  # Smoothed line
        axes[i].set_xlabel(metric)
        axes[i].grid(True)
        if log_x[i]:
            axes[i].set_xscale("log")        
    axes[0].legend()
    axes[0].set_title(models[0].file_path[-12:-4]+"_to_"+models[-1].file_path[-12:-4])
    axes[0].set_ylabel("DepthTop (mm)")
    axes[0].invert_yaxis()  # Invert y-axis (depth increases downward)
    if save_fig:
        plt.savefig(models[0].file_path[:-4]+"_to_"+models[-1].file_path[-12:-4]+"_"+metric+"_with_depth.png")
    return(fig,axes)
    





class plotter_model:
    """Handles standardized plotting for instances of the instability class."""
    def __init__(self, instability_instance):
        self.model = instability_instance
    
    
    def skier_impact_2D(self, wl_id, field = "Sxx"):
        """
        Plots skier induced stress inside the snowpack.
    
        Parameters:
        - wl_id: id of the weak layer (determines slab and weak layer properties
        - field: {'u', 'w', 'Sxx', 'Txz', 'Szz', 'principal'}, optional
                 Field quantity for contour plot. Axial deformation 'u', vertical
                 deflection 'w', axial normal stress 'Sxx', shear stress 'Txz',
                 transverse normal stress 'Szz', or principal stresses 'principal'
        Returns:
        - ax1
        """
        slab_profile = self.model._get_slab_profile_for_weac(wl_id)
        skier = weac.Layered(system='skier', layers=slab_profile)
        seg_skier = skier.calc_segments(L=self.model.skier_stability_params["totallength"],
                                        a=0, m=self.model.skier_stability_params["skierweight"])['nocrack']
        C_skier = skier.assemble_and_solve(phi=self.model.skier_stability_params["inclination"], **seg_skier)
        xsl_skier, z_skier, xwl_skier = skier.rasterize_solution(C=C_skier, phi=self.model.skier_stability_params["inclination"], **seg_skier)
        weac.plot.deformed(skier, xsl=xsl_skier, xwl=xwl_skier, z=z_skier,
                   phi=self.model.skier_stability_params["inclination"], window=400, scale=100, aspect=2,
                   field=field)
    
        return()
    
    def skier_stresses_in_wl(self, wl_id):
        """
        Plots stress at the surface to a given weak layer and the strength of the weak layer.
    
        Parameters:
        - wl_id: id of the weak layer (determines slab and weak layer properties
        Returns:
        - ax1
        """
        slab_profile = self.model._get_slab_profile_for_weac(wl_id)
        skier = weac.Layered(system='skier', layers=slab_profile)
        seg_skier = skier.calc_segments(L=self.model.skier_stability_params["totallength"],
                                        a=0, m=self.model.skier_stability_params["skierweight"])['nocrack']
        C_skier = skier.assemble_and_solve(phi=self.model.skier_stability_params["inclination"], **seg_skier)
        xsl_skier, z_skier, xwl_skier = skier.rasterize_solution(C=C_skier, phi=self.model.skier_stability_params["inclination"], **seg_skier)
        weac.plot.stresses(skier, x=xwl_skier, z=z_skier, **seg_skier)
        print("weak layer strength:" + str(self.model.profile.JS1999_sigma_macro.iloc[wl_id]))
    
        return
    
    def density_vs_sth(self, ax1_metric = "CR2020_density", ax2_metric = "rc_Reuter2015", ax1=None, padding = 1.8):
        """
        Plots two metrics vs. depth for a given instability instance.
    
        Parameters:
        - model: instability_modelling instance
        - ax1_metric: standard is Density
        - ax2_metric: instability metric (e.g. rc_Reuter2015)
        - ax1: optional to plot in a given axis
        - padding: to seperate the left and right profile
        Returns:
        - ax1, ax2: Axes objects for further customization
        """
        def normalize(data, padding, percentile = 90):
            return(data / data.max(), np.percentile(data.dropna(), percentile) * padding/ data.max())
        df1 = self.model.profile
        df2 = self.model.stab
        df = pd.concat([df1, df2.drop(columns=df1.columns.intersection(df2.columns))], axis=1).copy()
    
        ax1_color = "lightblue"
        ax2_color = "red"
        df["depthBottom"] = df["depthTop"] + df["thickness"]
        ax1_data, ax1_xmax = normalize(df[ax1_metric], padding)
        ax2_data, ax2_xmax = normalize(df[ax2_metric], padding)
        
        if ax1 is None:
            fig, ax1 = plt.subplots(figsize=(6, 8))
            ax2 = ax1.twiny()
        else:
            fig = ax1.figure  # Use the figure from ax1 if provided
            if ax2:
                pass
            else: ax2 = ax1.twiny()
    
        # Plot as blocks
        for i, row in df.iterrows():
            ax1.fill_betweenx(
                [row["depthTop"], row["depthBottom"]],
                0, ax1_data.iloc[i],
                color=ax1_color, edgecolor='black', linewidth=1, alpha=0.6 )
            ax2.fill_betweenx(
                [row["depthTop"], row["depthBottom"]],
                0, ax2_data.iloc[i],
                color=ax2_color, edgecolor='black', linewidth=1, alpha=0.6 )
    
        # Labels & Formatting
        ax1.set_xlabel(ax1_metric, color=ax1_color)
        ax1.set_ylabel("Depth (mm)")
        ax1.invert_yaxis()
        ax1.grid(True, linestyle="--", alpha=0.2)
        ax1.set_xlim(0, ax1_xmax)
        ax1.set_ylim(df["depthBottom"].max()*1.05, -50)
        ax1.tick_params(axis='x', colors=ax1_color)
        ax1.grid(True, linestyle="--", alpha=0.2)
    
        ax2.set_xlabel(ax2_metric, color="red")
        ax2.set_xlim(0, ax2_xmax)  # Ensure full density range + padding
        ax2.invert_xaxis()  # Depth increases downward
        ax2.tick_params(axis='x', colors=ax2_color)
        # Manually create legend elements
        ax1_dummy_line = ax1.plot([], [], color=ax1_color, label=ax1_metric)
        ax1.legend(handles=ax1_dummy_line, loc = "upper left")
        ax2_dummy_line = ax2.plot([], [], color=ax2_color, label=ax2_metric)
        ax2.legend(handles=ax2_dummy_line, loc = "upper right")    
        return (fig, ax1, ax2)

    def plot_metrics(self, metrics, max_percentile = 0.9, max_limits = None, save_plot = False, fig = None, axes = None, color = "C0"):
            
        df = self.model.stab.merge(self.model.profile)
        profiles = df[metrics]
        
        if fig is None:
            fig, axes = plt.subplots(1, profiles.shape[1], figsize=(profiles.shape[1] * 4, 6), sharey=True)  
        
        # Plots für jedes Profil
        for idx, (ax, profile) in enumerate(zip(axes, profiles)):        
            # Compute depth_end directly in DataFrame
            df['depthEnd'] = df['depthTop'] + df['thickness']
            y_values = df[['depthTop', 'depthEnd']].to_numpy().flatten()
            x_values = df[profile].repeat(2).to_numpy()
            ax.plot(x_values, y_values, marker='o', linestyle='-',c=color)
            if max_limits is None:
                max_x = df[profile].quantile(max_percentile)
            else: max_x = max_limits[idx]
            ax.set_xlim(-max_x/20,max_x)
            ax.set_xlabel(profile)
            ax.grid(True)
            
        axes[0].set_ylabel("Depth (mm)")
        plt.gca().invert_yaxis()  # Tiefe soll von oben nach unten gehen
        axes[0].set_title(self.model.file_path  +" - "+ str(self.model.PST_params)   +" - "+ str(self.model.skier_stability_params),loc="left",)
        if save_plot:
            plt.savefig(self.model.file_path[:-4]+"_stability_metrics_with_depth.png")
        return (fig, axes)
    
    
    def hazard_map(self, x = "depthTop", y = "logarithmic_sensitivity", color = "CR2020_ssa", ax1=None):
        """
            to be written
        """
        if ax1 is None:
            fig, ax1 = plt.subplots(figsize=(6, 8))
            
        self.model.stab.plot.scatter( x=x,y=y,c=self.model.profile[color],ax = ax1,title="Hazard map", colorbar=False)
        sc = ax1.collections[0]
        cbar = plt.colorbar(sc, ax=ax1)
        cbar.set_label(color)  # <-- Set the colorbar title here
        ax1.set_ylabel("Instability ("+y+")")
        ax1.set_xlabel("Slab thickness (mm)")

        ax1.invert_yaxis()
        return (fig, ax1)    
        
        
    # @classmethod
    # def run(cls, instance):
    #     """Creates an instance and runs all necessary computations."""
    #     logging.info("starting PostProcessor...")
    #     compute_harmonic_mean(instance, "S_Reuter2015", "rc_Reuter2015")
    #     compute_logarithmic_mean(instance, "S_Reuter2015", "rc_Reuter2015")
    #     compute_geometric_mean(instance, "S_Reuter2015", "rc_Reuter2015")
    #     compute_reziprocal_sum(instance, "S_Reuter2015", "rc_Reuter2015")
    #     compute_logarithmic_sensitivity(instance, "S_Reuter2015", "rc_Reuter2015")
    #     plotter.density_vs_sth(instance, ax2_metric="logarithmic_sensitivity")  
    #     return instance       