# Cross site comparison of ecosystem and plot scale methane fluxes across multiple sites
The R code used for processing and analyzing the EC and chamber methane (CH<sub>4</sub>) flux data used in the manuscript Määttä et al (2026) _A cross-site comparison of ecosystem- and plot-scale methane fluxes across multiple sites_ (in review at Biogeosciences). Prepint: https://doi.org/10.5194/egusphere-2025-5023

All code was created with R v4.3.3.

## Files

### Chamber data cleaning (Maatta_et_al_2026_Chamber_data_processing.R)
This file contains the R script for cleaning the raw chamber data sets. For the raw chamber data sets, please refer to the manuscript and the metadata of the published, temporally-aggregated multi-site datasets (available at https://doi.org/10.5281/zenodo.17312404) and contact the data providers.

### EC data cleaning
This file contains the R script for cleaning and subsetting the EC datasets obtained from the FLUXNET-CH<sub>4</sub> database (Delwiche et al. 2021, Knox et al. 2019; https://fluxnet.org/data/fluxnet-ch4-community-product/). 

### EC and chamber data combination
This file contains the R script for combining the EC and chamber data sets. First, the code combines the data and temporally aggregates the combined data per site, and then combines site-specific aggregations to multi-site temporal aggregations (half-hourly, hourly, daily, weekly, monthly, and annual). Finally, the code creates the published multi-site temporal aggregations (available at https://doi.org/10.5281/zenodo.17312404). Note: the unaggregated data sets are not part of the published data package but their processing is shown here nonetheless.

### Map (Maatta_et_al_2026_map.R)
This file contains the R script for the map in Fig. 1. 

### Statistical analyses (Maatta_et_al_2026_Statistical_analyses.R)
This file contains the R script used for all the statistical analyses and figures (except Fig. 1) in the manuscript. Please see the manuscript for more details on the methods and their rationale. When possible, the figures were created with colorblind-safe palettes (e.g., Okabe-Ito palette). 

The data sets used in this code can be downloaded at https://doi.org/10.5281/zenodo.17312404. 

## Contact
For further information, please contact Tiia Määttä (tiia.maatta@geo.uzh.ch)
