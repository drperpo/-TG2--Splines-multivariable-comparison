# STRATOS TG2- Splines project multivariable comparison


# Comparison of Multivariable Variable Selection Methods: MFP vs. Splines

## Project Overview

This repository contains code, data, and documentation related to the comparison of different statistical methods for multivariable model building, focusing on variable selection and functional form identification. The primary comparison is between the **Multivariable Fractional Polynomials (MFP)** approach and various **spline-based techniques**, including Multivariable Regression Splines (MVRS) and penalized splines within Generalized Additive Models (GAMs).

The first part relates to analysis presented in ISCB2022 on the PIMA Indians data. The newest part extends the analysis to bacteremia data (2024). 


## Repository Structure

The repository is organized as follows:

* **/Rscripts**: Contains R and R Markdown scripts for analyses.
    * `ISCB2022-analysis.Rmd`: Main analysis comparing MFP and spline methods ( using PIMA Indians data initially). Should have all the background code that was used to create the presentation. 
     * `mvrs.R`: Script  implementing / exploring the MVRS procedure.
It is  very crude and may contain errors. `
    * `MGCV_spline_var_selection_on_bact_data.R`: Specific analysis using `mgcv` (GAMs with penalized splines - TP, TS, PS basis) for variable selection on the Bacteremia dataset. *(This is  the most recent script).*
* **/data**: Contains the datasets used in the analyses. Includes data like:
    * PIMA Indians Diabetes (`diabetes.csv`)
    * Bacteremia study data (`Bacteremia_public_S2.csv`, `bacteremia-DataDictionary.csv`)
 
* **/docs**: Contains supplementary documents, presentations, and protocols.
    * `ISCB22-perperoglou-splines-presentation.pptx`: Presentation slides (likely from ISCB 2022).
    * `Draft protocol.docx`, `Aris-iscb-Suggestions-ws.docx`: Planning documents or drafts related to the project.
* `README.md`: This file, providing an overview of the repository.

## Methods Explored

* Multivariable Fractional Polynomials (MFP)
* Multivariable Regression Splines (MVRS)
* Generalized Additive Models (GAMs) with penalized splines for automatic variable/smoothness selection using `mgcv`:
    * Thin Plate splines (`bs="tp"`)
    * Tensor Product splines (`bs="ts"`)
    * P-splines (`bs="ps"`)
    * Cubic Regression splines (`bs="cr"`) / Natural Splines (`splines::ns`)

## Software & Packages

The analyses are primarily conducted in **R**. Key packages used include:

* `mgcv`
* `mfp`
* `splines`
* `glm` (base R)
* `dplyr`
* `tidyr`
* `ggplot2`

## How to Use

1.  Clone the repository.
2.  Ensure you have R and the necessary packages installed (see list above).
3.  Explore the scripts in the `/Rscripts` folder. The `.Rmd` file can be rendered, and the `.R` scripts can be run. Note that scripts might depend on data located in the `/data` folder.

---
