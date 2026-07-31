# Widening sexual size dimorphism in bees over a century of climate change
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21672792.svg)](https://doi.org/10.5281/zenodo.21672792)
# Abstract
Shrinking body size is a pervasive consequence of global change, with implications for biodiversity and ecosystem functioning. Still, whether males and females differ in their vulnerability to shrinkage remains poorly resolved, despite the destabilising potential of sex-specific size trajectories. Further, disentangling the spatial and temporal drivers of size variation remains a persistent challenge. For bees, body size governs pollination efficiency and temperature sensitivity, yet evidence for size shifts is mixed and fragmented. We leverage natural history collections and community science measurements to investigate temporal and spatial trends in body size across 136 bee species sampled over more than a century of climate change and land conversion in California (1900 – 2024). Male bees declined significantly in size since 1900, while females remained stable. The male size decline was driven by local precipitation, with smaller bees in drier years, accompanied by increasing precipitation variability over the study period. Across spatial gradients, males were also more sensitive than females to environmental variation, but with size variation driven by land use, not climate. This divergence in spatial and temporal drivers illustrates the risk of predicting temporal trends from spatial comparisons and underscores the utility of datasets spanning broad environmental gradients and long time series. Further, these sex-specific responses raise the possibility of male sensitivity to resource scarcity and/or strong selection to maintain female body size. Together, these findings reveal male vulnerability and female resilience to environmental variability, with implications for pollination services and population stability under ongoing climate change.

# Repository Directory

### Code:

**CA_size_analysis.R** is the main script for the analysis presented in the paper. Uses grid_filtered_CA_size_data.csv as input.

**CA_size_figures.R** produces the figures in the paper. Run CA_size_analysis.R

**CA_size_spatiotemporal binning.R** is the script which takes all specimen data and filters them into 150x150km spatiotemporal bins, as described in the methods. Takes before_filtering_CA_size_data.csv as input and produces produces grid_filtered_CA_size_data.csv, which is used for the main analysis.


### Data:

**grid_filtered_CA_size_data.csv** includes the size data used in the main analysis (CA_size_analysis.R). These data include only the specimens retained after spatiotemporal binning.

**before_filtering_CA_size_data.csv** presents the size data used in the main analysis PRIOR to spatiotemporal binning, done in (CA_size_spatiotemporal binning.R).

**median_itd_global.csv** contains a global sample of >23K specimen measurements conducted by community scientists on Notes from Nature. These data include the California size dataset analysed in the paper as well as other measurements from other parts of the world, not used in the paper. Here, size data are presented as the median ITD measurement from all measurements for a given specimen, as recommended in Ostwald et al. 2025 (https://doi.org/10.1002/ece3.71665).

**raw_itd_global.csv** contains the global size dataset presented in median_itd_global.csv, but as raw measurements rather than specimen-level medians.

Specimen data are presented in standardised trait format described in: https://github.com/mostwald/Functional-trait-review.
