# Widening sexual size dimorphism in bees over a century of climate change
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17599627.svg)](https://doi.org/10.5281/zenodo.17599627)

# Abstract
Shrinking body size is a pervasive consequence of global change, with implications for biodiversity and ecosystem functioning. How size trends may differ among the sexes is poorly resolved, despite the destabilising potential of divergent responses. For bees, body size governs pollination efficiency and temperature sensitivity, yet evidence for size shifts is mixed and fragmented. We leverage natural history collections and community science measurements to disentangle temporal and spatial trends in body size across 136 bee species sampled over more than a century of climate change and land conversion in California (1900 – 2024). Male bees declined significantly in size since 1900, while females remained stable. The male size decline was driven by local precipitation, with smaller bees in drier years, accompanied by increasing precipitation variability over the study period. Across spatial gradients, males were also more sensitive than females to environmental variation, but with size variation driven by land use, not climate. This divergence in spatial and temporal drivers illustrates the risk of predicting temporal trends from spatial comparisons and underscores the utility of datasets spanning broad environmental gradients and long-time series. Further, these sex-specific responses raise the possibility of male sensitivity to resource scarcity and/or strong selection to maintain female body size. Together, these findings reveal male vulnerability and female resilience to environmental variability, with implications for pollination services and population stability under ongoing climate change.

# Repository Directory

### Code: Contains code for data analysis in R

**3D_heat_budget_analysis_github.R**: In this document, we analyze morphometric data (ITD, head width, costal vein length, dry mass, volume, surface area) to understand inter- and intra-specific predictors of body size.

## Data: Contains the raw morphometric and specimen data in two formats

**grid_filtered_CA_size_data.csv**: In this document, we present the size data used in the main analysis (ADD SCRIPT NAME.R). These data include only the specimens retained after spatiotemporal binning.

**before_filtering_CA_size_data.csv**: In this document, we present the size data used in the main analysis PRIOR to spatiotemporal binning, done in (insert script name).

**roberts_model_long.csv**: This document contains heat flux data from Roberts and Harrison, 1999, as well as longwave radiation and convective heat loss data estimated from 3D model data. These data are used to estimate heat budgets in the analysis code (3D_heat_budget_analysis_github.R).

**3D_apis_trait_format.csv**: In this document, we present the same raw measurement data as above, but in a standardized format that lends itself to functional trait data sharing (see this [functional trait data sharing repository](https://github.com/mostwald/Functional-trait-review) for descriptors of column headers).

Associated models and images can be found at the linked Zenodo repository.
