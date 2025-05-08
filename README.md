# snow_models

Welcome to **snow_models** – a Python package designed for the easy and consistent application of snow- and avalanche-related parameterizations and models. This package provides several models and functions to support research and development in snow mechanics and avalanche forecasting.

## Key Modules

- **`snow_models.mechanical_params`**  
  Parametrizations and conversion functions for mechanical properties of snow.

- **`snow_models.crack_propagation`**  
  Models to estimate crack propagation speeds in weak snow layers.

- **`snow_models.wave_propagation`**  
  Models to compute elastic wave speeds for different types of waves.

## Purpose

This package implements published parameterizations to facilitate further research applications. It is designed to help researchers and practitioners integrate and use these models for snow and avalanche studies in an easy-to-use Python environment.

## Installation

You can install **smp_instability** directly from GitHub using `pip`
(make sure git is installed in your active environment (conda install git):

```bash
pip install git+https://github.com/bbergfeld/smp_instability.git
```

## Examples

Example usage to derive instability from a pnt file:
```python
from smp_instability.pipelines import plot_smp_instability

pnt_file = 'D:\\Github_repos\\smp_instability\\smp_instability\\.pnt_data\\SMP1_FILE0003.pnt'
output_folder = 'D:\\example_output_folder_for_smp_2_instability'

plot_smp_instability(pnt_file, output_folder)


```



