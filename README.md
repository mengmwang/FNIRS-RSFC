# fNIRS-RSFC

## Signal Analysis in fNIRS RSFC

Functional near-infrared spectroscopy (fNIRS) is widely used as a non-invasive neuroimaging modality. FNIRS resting state functional connectivity (RSFC) can be measured using the correlation-based method via Pearson's correlation coefficient. 

Assumptions are imposed on the fNIRS signals when performing correlation-based connectivity analyses. Pearson's correlation coefficients and the subsequent statistical tests assume stationarity and whiteness of signals. Violation of any one of these assumptions may invalidate either the sample correlation estimate, or the statistical significance results.

The following two aspects of fNIRS signal properties are considered, 1) autocorrelated fNIRS signals violate the whiteness assumption, 2) non-stationary components in fNIRS signals violate the stationarity assumption. Signal analysis solutions that correct the statistical tests are proposed based on the considerations.  

### Correct for coloured frequency spectra in fNIRS signals

*Abstract: *

### Correct for non-stationary components in fNIRS signals




### Code files

**Correction methods to correct for the assumption violations**

Code for simulation signals 

The methods can be applied in experimental fnirs data

Check [my Google site](https://sites.google.com/view/mengmengwang/research) for detailed methods and publications

`nw_sim.m` : a statistical correction method to correct for the effects for non-white signals in sample correlation distributions

`ns_sim.m` : a correction method to correct for the effects for multiplicative non-stationary signals in sample correlation distributions

*Note: The code files run in Matlab, include generate simulated signals, apply the proposed corrections and generate results*
