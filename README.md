# Bering Sea Marine Heatwave Consequences

**Project Leads:** Erin Fedewa, Mike Litzow and Emily Ryznar 

This repository contains the source code for data processing and analytical workflows for the following manuscript:

> **Title:** A Bering Sea marine heatwave triggers unprecedented crab hybridization, challenging fisheries management
> 
> **Authors:** Erin J. Fedewa*, Emily R. Ryznar, Connor Cleary, Jennifer L. Gardner, Shannon M. Hennessey, W. Christopher Long, Jon I. Richar, Leah S. Zacher and Michael A. Litzow
> 
> **Citation:** To be updated with DOI following acceptance

---

## Project Overview
* A 2018-2019 marine heatwave in the Bering Sea resulted in nearly ice-free winter conditions, leading to mass mortality of snow crab. Although the mass mortality event was the most immediate outcome of the heatwave, a number of emergent consequences have appeared in the following years. Here, we use a fisheries-independent trawl survey time series to evaluate competing hypotheses for the cascading impacts of a marine heatwave on Bering Sea crab stocks.

---

## Software Dependencies
* **Language environment:** R v4.3.2+ 
* **R Package dependencies:** Bering Sea crab abundance estimates are produced using the `crabpack` R package (Hennessey 2024). Dynamic structural equation models are fitted using the `dsem` R package (Thorson et al. 2024) 

---

## Data Availability
* **Processed Data:** Processed files required to run the final models and figures are included directly in the `output/` directory. 

---

## Repository Structure
To reproduce the manuscript results and figures, execute the pipeline scripts sequentially. 
```text

├── script/
│   ├── get_crab_abundance/         # Produce Bering Sea crab population abundance estimates
│   ├── get_sea_ice/                # Calculate Bering Sea ice concentration
│   ├── get_spatial_overlap/        # Calculate spatial overlap of mature Tanner crab and snow crab
│   ├── assess_lags/                # Assess cross-correlations of biologically-constrained lags
│   ├── dsem_tanner_model/          # Final Tanner crab DSEM model, diagnostics and simulation
│   ├── dsem_hybrid_model/          # Final hybrid crab DSEM model, diagnostics and simulation
│   ├── make_figure_1/              # Script for creating Figure 1
│   ├── make_figure_2/              # Script for creating Figure 2
│   ├── Additional_explorations/    # Exploratory script only, not used in final analysis
├── data/
│   ├── ERA5_ice_/                  # Sea ice data pulled from Climate Data Store API
│   ├── station_lookup.csv/         # EBS survey strata used for subsetting data 
├── output/
│   ├── crab_abundance/             # Output from get_crab_abundance; abundance time series of Tanner, snow and hybrid crab
│   ├── hybrid_final_dsem_summary/  # Output from dsem_hybrid_model; model estimates and variance
│   ├── hybrid_lags/                # Output from dsem_hybrid_model; model estimates across tested lags
│   ├── overlap_output/             # Output from get_spatial_overlap
│   ├── proportion_hybrid/          # Output from get_crab_abundance; % of 4 inch males that are hybrid and snow crab
│   ├── seaice_output/              # Output from get_sea_ice; average spring ice extent 
│   ├── tanner_final_dsem_summary/  # Output from dsem_tanner_model; model estimates and variance
│   ├── tanner_lags/                # Output from dsem_tanner_model; model estimates across tested lags
├── figures/                        # Figure 1 and Figure 2 output

```

---

## NOAA License
This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is provided on an ‘as is’ basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.

---

## Contact
For technical questions regarding the code, please open an issue in this repository or contact the corresponding author:
* **Name:** Erin Fedewa
* **Email:** erin.fedewa@noaa.gov

