# README.md

The scripts in this repo can be used to run the analyses conducted in the Stagaman *et al*. 2024 article, {TITLE} published in *Communications Medicine*.

File descriptions:

1. `.Rprofile` creates a list of directories for reading in and saving data and plots
2. `_packages_and_sources.R` 
3. `_setup.R` loads packages, sources functions, set important variables
4. `01_variable_reduction_selection.R` Conducts the variable reduction step to identify the covariates to be used in all subsequent analyses.
5. `02_Analyses` all scripts used for statistical analyses
    1. `01_required_first.R` removes low abundance and low prevalence OTUs/KOs, rarefies counts, and generates data structures to be used in differential abundances and random forest analyses
    2. `02_taxonomic_alpha.R` conducts statistical analysis of alpha-diversity in relation to PD and the covariates of interest.
    3. `03_all_beta.R` conducts statistical analysis of beta-diversity in relation to PD and the covariates of interest.
    4. `04_random_forests` scripts used for random forest classfier models
        01. `04A_covariates_split_sources.R` RF model training and testing for saliva and stool separately using only covariate data (no OTU or KO abundances).
        02. `04B_covariates_combined_sources.R` RF model training and testing for saliva and stool combined using only covariate data (no OTU or KO abundances).
        03. `04C_taxon_abunds_split_types_split_sources.R` RF model training and testing for saliva and stool separately using OTU, species, and genus abundances separately.
        04. `04D_taxon_abunds_split_types_combined_sources.R` RF model training and testing for saliva and stool combined using OTU, species, and genus abundances separately.
        05. `04E_taxon_abunds_aggregated_split_sources.R` RF model training and testing for saliva and stool separately using OTU, species, and genus abundances combined.
        06. `04F_taxon_abunds_aggregated_combined_sources.R` RF model training and testing for saliva and stool combined using OTU, species, and genus abundances combined.
        07. `04G_taxon_abunds_covariates_split_types_split_sources.R` RF model training and testing for saliva and stool separately using OTU, species, and genus abundances separately, including covariates.
        08. `04H_taxon_abunds_covariates_split_types_combined_sources.R` RF model training and testing for saliva and stool combined using OTU, species, and genus abundances separately, including covariates.
        09. `04I_taxon_abunds_covariates_aggregated_split_sources.R` RF model training and testing for saliva and stool separately using OTU, species, and genus abundances combined, including covariates.
        10. `04J_taxon_abunds_covariates_aggregated_combined_sources.R` RF model training and testing for saliva and stool combined using OTU, species, and genus abundances combined, including covariates.
        11. `04K_function_abunds_split_types_split_sources.R` RF model training and testing for saliva and stool separately using KO, module, and pathway abundances separately
        12. `04L_function_abunds_covars_split_types_split_sources.R` RF model training and testing for saliva and stool separately using KO, module, and pathway abundances separately, including covariates.
    5. `05_LDAs` scripts used for differential abundance analsyes (a.k.a linear discriminate analyses)
        1. `05A_taxon_ANCOMBC2_wCovariates.R` Differential OTU abundance analysis using ANCOMBC2, including covariates.
        2. `05B_taxon_ANCOMBC2_PDonly.R` Differential OTU abundance analysis using ANCOMBC2, without covariates.
        3. `05C_taxon_ALDEx2_wCovariates.R` Differential OTU abundance analysis using ALDEx2, including covariates.
        4. `05D_taxon_ALDEx2_PDonly.R` Differential OTU abundance analysis using ALDEx2, without covariates.
        5. `05E_function_ANCOMBC2_wCovariates.R` Differential KO abundance analysis using ANCOMBC2, including covariates.
        6. `05F_function_ANCOMBC2_PDonly.R` Differential KO abundance analysis using ANCOMBC2, without covariates.
        7. `05G_function_ALDEx2_wCovariates.R` Differential KO abundance analysis using ALDEx2, including covariates.
        8. `05H_function_ALDEx2_PDonly.R` Differential KO abundance analysis using ALDEx2, without covariates.
    6. `06_networks` scripts used to generate feature-feature (OTUs/KOs) networks
        1. `06A_between_saliva_stool_associations.R` Generates OTU-OTU and KO-KO co-abundance networks (utilizing SpeicEasi) between saliva and stool microbiomes.
        2. `06B_within_sample_associations.R` Generates OTU-OTU and KO-KO co-abundance networks (utilizing SpeicEasi) within saliva and stool microbiomes.
    7. `99_required_last.R` merges data from across different scripts in prepartion for plotting and further analysis.
6. `Helpers` scripts for custom functions and model specifications used across scripts
    1. `functions.R` custom functions
    2. `model_specs.R` model specifications
