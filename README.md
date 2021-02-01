# T2-DWI
This repository provides functions for performing three different model fits to combined T2- and diffusion-weighted (T2-DWI) data acquired using two echo times (TEs) and two b-values:
1. Two-component model
2. Bi-exponential model
3. Mono-exponential apparent diffusion coefficient (ADC)

The method was developed at the MR Cancer group at the Norwegian University of Science and Technology (NTNU) in Trondheim, Norway. https://www.ntnu.edu/isb/mr-cancer

# Note
The provided code was developed for research and is not approved for clinical use.

# How to cite T2-DWI
In case of using or referring to T2-DWI, please cite it as:

Syversen IF, Elschot M, Bathen TF, Goa PE. Two-component model of prostate tissue using hybrid multidimensional T2 and diffusion-weighted imaging. (abstract 2419) In: Proceedings of the 2020 ISMRM & SMRT Virtual Conference & Exhibition, 2020.

# How to use T2-DWI
This is a MATLAB® function, and the function was written and tested using MATLAB® R2019b.

### Two-component model:
Use the calculateTwoComponent.m function, which uses the function hybridfit_TCmodel.m to perform the actual model fit.

### Bi-exponential model:
Use the calculateBiExponential.m function, which uses the function hybridfit_BE.m to perform the actual model fit.

### Mono-exponential ADC:
Use the calculateADC.m function.

Note that for all three models, there are certain functions that have to be implemented by the user itself, e.g. for loading of image files.

# Contact us:
Feel free to contact us: ingrid.f.syversen@ntnu.no
