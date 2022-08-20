# fNIRS-RSFC

## Signal Analysis in fNIRS RSFC

Functional near-infrared spectroscopy (fNIRS) is widely used as a non-invasive neuroimaging modality. FNIRS resting state functional connectivity (RSFC) can be measured using the correlation-based method via Pearson's correlation coefficient. 

Assumptions are imposed on the fNIRS signals when performing correlation-based connectivity analyses. Pearson's correlation coefficients and the subsequent statistical tests assume stationarity and whiteness of signals. Violation of any one of these assumptions may invalidate either the sample correlation estimate, or the statistical significance results.

The following two aspects of fNIRS signal properties are considered, 1) autocorrelated fNIRS signals violate the whiteness assumption, 2) non-stationary components in fNIRS signals violate the stationarity assumption. Signal analysis solutions that correct the statistical tests are proposed based on the considerations.  

### Correct for coloured frequency spectra in fNIRS signals

*Abstract: FNIRS is a non-invasive neuroimaging modality for monitoring brain oxygenation levels. Resting-state functional connectivity (FC) methods applied to fNIRS data assume that the signals are white. If this assumption is violated, statistical tests are invalid and connectivity may be artificially recorded. Temporal filtering is known to cause a violation of the whiteness assumption, and corrections to statistical tests of connectivity have been proposed for fMRI data. However, such corrections assume that unfiltered signals are white. It is known that fNIRS signals are typically not white prior to preprocessing. We propose a correction to statistical tests of functional connectivity, that accounts for both non-white fNIRS data and the further effects of filtering.*

- M. Wang, C. Davey and L. Johnston, "Correction of induced functional connectivity in filtered resting state fNIRS data," *The 27th Annual Meeting of the Organization for Human Brain Mapping (OHBM)*, Virtual, 2021. [[Full paper](https://github.com/mengmwang/fNIRS-RSFC/blob/main/Wang2021_1.pdf), [Poster](https://github.com/mengmwang/fNIRS-RSFC/blob/main/Wang2021_Poster_1.pdf)]

### Correct for non-stationary components in fNIRS signals

*Abstract: Functional Near-Infrared Spectroscopy (fNIRS) signals are contaminated by various sources of noise, including physiological noise and motion artefacts, both of which can change over time. We empirically establish the non-stationarity of fNIRS signals, identifying the distribution of the time-varying signal power. Non-stationarity of fNIRS signals violates the stationarity assumption imposed by resting-state functional connectivity (RSFC) measures, invalidating statistical significance tests. We propose a simple correction for the time-varying signal power of fNIRS signals that restores the integrity of RSFC analyses. The correction is drawn from a similar method proposed for non-stationary functional Magnetic Resonance Imaging (fMRI) data, that was analytically established to restore the voracity of RSFC tests.*

- M. Wang, C. Davey and L. Johnston, "Correction for time-varying signal power in fNIRS connectivity analyses," *Society of fNIRS Virtual Conference 2021 (fNIRS2021)*, Virtual, 2021. [[Full paper](https://github.com/mengmwang/fNIRS-RSFC/blob/main/Wang2021_2.pdf), [Poster](https://github.com/mengmwang/fNIRS-RSFC/blob/main/Wang2021_Poster_2.pdf)]

### MATLAB Implementation

This repository contains the MATLAB code for the signal and statistical correction methods to correct for the non-white and non-stationary fNIRS signals in fNIRS RSFC (simulations). 

**Code Demo for Simulated Signals**

`nw_sim.m` - A statistical correction method to correct for the effects for non-white signals in sample correlation distributions

`ns_sim.m` - A correction method to correct for the effects for multiplicative non-stationary signals in sample correlation distributions

**Note** 

The code files run in MATLAB, simply type the file names `nw_sim` and `ns_sim` in the MATLAB command. Results will be generated in various figures. 

The code files in this repository contain simulation signals and results only. The correction methods can be applied in experimental fNIRS data. 

This is part of my PhD research. Please get in touch if you are interested to know more! :wink:
