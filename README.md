# Realistic 3D morphology reshapes insect heat budgets
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17599627.svg)](https://doi.org/10.5281/zenodo.17599627)

# Abstract
Modeling insect heat exchange and predicting thermal responses depends on accurate representation of body size and shape. Still, most biophysical models approximate these complex forms using simplified geometric solids, whose relationships to real body forms have not been rigorously tested. Advances in surface modeling of small objects allow us to interrogate these assumptions by capturing the real 3D complexity of insect body forms. We used photogrammetry to construct 3D models of honey bee specimens and empirically measured body volume and surface area. Compared to empirical measurements, we found that traditional, geometric size estimation methods systematically underestimate body surface area and volume. We incorporated these error estimates into published heat budget data and found that these errors propagated non-linearly through the model, shifting the relative dominance of convective and radiative heat loss as temperature increases. These results suggest that body size and surface area assumptions can distort modeled heat fluxes, particularly under high temperatures, demonstrating that morphological simplifications can bias physiological inference. This work underscores the utility of empirical 3D morphology for refining biophysical models of insect thermoregulation.

# Repository Directory

### Code: Contains code for data analysis in R

**3D_heat_budget_analysis_github.R**: In this document, we analyze morphometric data (ITD, head width, costal vein length, dry mass, volume, surface area) to understand inter- and intra-specific predictors of body size.

## Data: Contains the raw morphometric and specimen data in two formats

**3D_heat_budget_measurement_data_wide_format.csv**: In this document, we present the raw measurement data used in the analysis (3D_heat_budget_analysis_github.R).

**roberts_model_long.csv**: This document contains heat flux data from Roberts and Harrison, 1999, as well as longwave radiation and convective heat loss data estimated from 3D model data. These data are used to estimate heat budgets in the analysis code (3D_heat_budget_analysis_github.R).

**3D_apis_trait_format.csv**: In this document, we present the same raw measurement data as above, but in a standardized format that lends itself to functional trait data sharing (see this [functional trait data sharing repository](https://github.com/mostwald/Functional-trait-review) for descriptors of column headers).

Associated models and images can be found at the linked Zenodo repository.
