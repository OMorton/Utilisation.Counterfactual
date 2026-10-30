# Historic and Ongoing Impacts of Human Utilization on Population Abundance

This repository contains the R code and analytical pipeline supporting:

> **Morton, O.***, McRae, L., O'Brien, D., & Johnson, T. F. — *"Historic and ongoing impacts of human utilization on population abundance" (2026)*

---

## 📌 Repository Overview

This project provides an end-to-end workflow to process Living Planet Database (LPD) records, fit Generalized Additive Mixed Models (GAMMs) accounting for auto-regressive structures, spatial correlation and co-occuring stressors, and extract chained temporal Average Treatment Effects on the Treated (ATT) denoting the impact of utilisation on population abundance trajectories, and reproduce the manuscript's figures.

---

## 📂 Project Structure & Script Summary

### 🛠️ Helper Functions
* **`Functions.R`**: Contains utility and helper functions used across the analysis to load and clean datasets, evaluate formula rank deficiency, and generate diagnostic and ATT plots.

---

### 1. Data Preparation (`Data.Preparation/`)
* **`01.PrepareLPD.R`**: Processes and filters raw LPD data. Excludes populations with unknown threat or utilization statuses and filters time-series by length criteria.
* **`02.Explore.R`**: Analyzes the spatial and taxonomic coverage of the processed LPD dataset and generates **Figure 1**.

---

### 2. Model Fitting (`Model.Fitting/`)
* **`03.Global.model.R`**: Fits an AR-naive GAMM to estimate the auto-regressive coefficient $\text{AR}(1)$, which is subsequently used to parameterize the final AR-informed model.
* **`04.Global.model.preds.R`**: Extracts posterior simulations from the global GAMM and calculates the chained Average Treatment Effect on the Treated (ATT).
* **`05.Disaggregated.models.R`**: Fits AR-informed GAMMs (using the pipeline from Script 03) across subsetted groups, including taxonomic groups, systems, geographic regions, and national levels.
* **`06.Disaggregated.model.preds.R`**: Draws simulated posteriors from the disaggregated GAMMs and calculates group-specific chained ATTs.

---

### 3. Visualization (`Plotting/`)
Processes the calculated ATTs and generates the figures presented in the paper:
* **`class.figures.plotting.R`**: Generates taxonomic class-level visualizations (**Figure 2**).
* **`region.figures.plotting.R`**: Generates geographic region-level visualizations (**Figure 3**).
* **`threat.figures.plotting.R`**: Generates threat-level and utilization visualizations (**Figure 4**).

---

The raw data were obtained from the LPD which can be requested at: www.livingplanetindex.org/data_portal (current data accessed 6th August 2026). All use of the dataset complied with the LPD’s Data Use Policy. The LPD’s Data Use Policy states that we cannot “pass the original datasets on to third parties and…direct all requests for use of the data back to LivingPlanetIndex@ioz.ac.uk” because of confidentiality concerns. 

---
