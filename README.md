## Unsupervised Anomaly Detection in PVP Signals with Hidden Markov Models

- This repository contains the codes of the paper 'Unsupervised Anomaly Detection in Peripheral Venous Pressure Signals with Hidden Markov Models'. 
- This paper is going to be published in [Biomedical Signal Processing and Control](https://www.journals.elsevier.com/biomedical-signal-processing-and-control) (currently under second review).
- (Till now) We have decided _not_ to publish the dataset used in this study.
- Please send an email to `mahayat@uark.edu` or `wuj@uark.edu` for any further queries.

## Description

- `MLEestimate.R` - MLE estimates of DLM parameters using `dlm` package of R
- `Kalman_Filter.m` - Kalman filter on DLM using estimates
- `HMM_GMM.m` - HMM based state estimation assuming Gaussian distribution
