# smp_instability

Welcome to **smp_instability** – a Python package to compute snowpack instability based on Snow Micro Penetrometer measurements.

## Key Modules

- **`smp_instability.pre_processor`**  
  
This module loads .pnt files measured with the Snow Micro Penetrometer. It processes these profiles to derive physical and mechanical properties relevant for snow instability analysis. 
    The key functionalities include:
    - Loading .pnt files and extracting profile data.
    - Computing microstructural features (e.g., density, SSA) using established parameterizations.
    - Segmenting the snowpack into layers via KMeans clustering based on weighted feature vectors.
    - Calculating mechanical properties per layer (e.g., fracture energy, compressive strength, elastic modulus).
    - Supporting batch processing for directories or lists of .pnt files.
    - Optionally visualizing and saving the processed profiles.
	
	
- **`smp_instability.instability_modelling`**  
	Implements routines for evaluating snowpack instability based on the method proposed by Reuter et al. (2015). 
	The core class `R2015_point_instability` analyzes preprocessed snow profiles to compute:

	- Skier-induced stability index (S_Reuter2015)
	- Critical crack length (rc_Reuter2015) using energy release modeling

	Features:
	- Works with individual profiles, lists of profiles, files, or folders
	- Customizable parameters for skier and PST (Propagation Saw Test) simulations
	- Parallelized batch processing
	- Integrated logging and error handling
	- Optional saving of results and model instances as pickle files

- **`smp_instability.post_processing`**  
	Visualization and Metric Combination for Snowpack Instability Analysis

## Purpose

This package implements pipelines to derive an instability assessment of a snowpack. It is based on mechanical modelling and fully automated.

## Installation

You can install **smp_instability** directly from GitHub using `pip`
(make sure git is installed in your active environment (conda install git):

```bash
pip install git+https://github.com/bbergfeld/smp_instability.git
```

## Examples

Example usage to derive instability from a pnt file:
```python
from smp_instability.pipelines import plot_instability_profile

pnt_file = 'filepath_to_pnt_file'
output_folder = 'outputfolder_to_save_plots'

plot_instability_profile(pnt_file, output_folder)


```



