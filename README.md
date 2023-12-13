# BONDpipeline

This repository contains the scripts for the BOND pipeline and scripts for subsequent analysis and visualisation.

1) The BOND pipeline was developed for pre-processing of EEG power and connectivity and is based on the MADE pipeline. 

This pipeline is based on the Maryland Analysis of Developmental EEG (MADE) pipeline with small adjustments to make it suitable for the Longitudinal European Autism Project (LEAP) dataset and analyses of EEG connectivity. 
To run this pipeline, the dependencies of the MADE pipeline will need to be downloaded as well. We refer to the MADE paper for the code: Debnath, R., Buzzell, G. A., Morales, S., Bowers, M. E., Leach, S. C., & Fox, N. A. (2020). The Maryland analysis of developmental EEG (MADE) pipeline. Psychophysiology, 57(6), 1–13. https://doi.org/10.1111/psyp.13580 

EEG power and connectivity metrics from the BOND pipeline were compared with metrics calculated after pre-processing with a manual pipeline, the MADE pipeline, and the Harvard Automated Processing Pipeline for EEG. 

2) The scripts 'LEAPpipelines_xxx' were used for subsequent analyses. LEAPpipelines_1_EEGmetrics.m calculates the EEG metrics (power and connectivity) for trials across all conditions, the social video, and the non-social video. LEAP_2_BetweenPipelineComparisons.m calculates the intra-class correlation between the pipelines across the delta, theta, alpha, and beta canonical frequency bands for all trials and condition differences (social - non-social video). LEAPpipelines_3_WithinPipelineComparisons.m selects alternating trials for all conditions, and for social and non-social conditions only. It then calculates the correlations between the split-half datasets as measure of internal consistency.

3) The folder 'Visualisation' contains scritps used to create the figures in the manuscript.

(written by Rianne Haartsen, PhD., Dec-2023)
