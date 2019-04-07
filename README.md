# AxialOpt
AxialOpt is a code for the preliminary design and optimization of axial turbines. The output of AxialOpt can be used in system-level analyses (such as a power cycle optimization) to estimate the efficiency or footprint of axial turbines for a given set of thermodynamic specifications. In addition, the information provided by AxialOpt can be used as the starting point for the aerodynamic design of the turbine blades using more advanced flow models based on CFD.

The models and optimization methodology of AxialOpt are documented in a peer-reviewed, [open-access publication](#1) and the source code is also stored in a [Zenodo repository](#2).

## Features

* The axial turbine model is composed of three sub-models that are used as building blocks:
  1. A meanline model to describe the flow in each cascade
  2. An empirical loss model to evaluate the entropy generation in each cascade
  3. A [one-dimensional flow diffuser model](https://github.com/RoberAgro/AnnularDiffuser1D) to compute the exit kinetic energy recovery
* The model is formulated for axial turbines with any number of stages
* The model is formulated to use arbitrary equations of state to compute thermodynamic properties:
  1. The current version uses the [REFPROP v9.1](#3) fluid library
  2. Other fluid libraries and look-up tables may be implemented in the future
* The loss model is formulated in a general way to use:
  1. Any set of empirical correlations to compute the losses:
      1. [Ainley and Mathieson](#4) (implemented)
      2. [Dunhan and Came](#5) (implemented)
      3. [Kacker and Okapuu](#6)(implemented)
      4. [Craig and Cox](#7) (will be implemented soon)
      5. Other loss model contributions are welcome
  2. Different definitions for the loss coefficient:
      1. Stagnation pressure loss coefficient
      2. Enthalpy loss coefficient (also known as kinetic energy loss coefficient)
      3. Entropy loss coefficient
* The preliminary design is formulated as a constrained optimization problem
  1. This allows explore the design space in a systematic way and account for technical constraints
  2. It is straighforward to modify the objective function and constraints depending on the problem
  3. There are several gradient-based algorithms available to solve the optimization problem, including:
      1. Sequential Quadratic Programming (SQP)
      2. Interior Point (IP)
* The output of AxialOpt can be saved as:
  1. Text file containing a summary of the solution of the optimization problem
  2. Text files containing the geometry of each cascade and variables at each station
  3. Figures with the T-s and h-s diagrams of the expansion
  4. Figures with the the velocity triangles
  5. Figures with the axial-radial and axial-tangential views of the turbine
  6. Figures with the breakdown of the losses within the turbine
  
_Note: AxialOpt is not suitable to estimate the performance of an existing design under different conditions. The extension of the code to compute the performance at off-design conditions is underway._



## Requisites
AxialOpt was implemented in [MATLAB R2018a](https://nl.mathworks.com/) and requires a [REFPROP v9.1](#3) installation. The folder [link_refprop_matlab](link_refprop_matlab) contains a instructions about how to link REFPROP with MATLAB.

AxialOpt has the option to use the [export_fig library](https://github.com/altmany/export_fig) to produce publication-quality figures. Using this library requires _ghostcript_ and _pdftops_. See the installation instructions in the [original repository](https://github.com/altmany/export_fig).



## Examples
The folder [AxialOpt_examples](AxialOpt_examples) contains five examples commented in detail to get started with the code:
  * A supercritical Carbon dioxide turbine
  * A organic Rankine cycle (ORC) turbine using R125 as working fluid
  * A organic Rankine cycle (ORC) turbine using hexane as working fluid
  * An industrial gas turbine
  * A high-pressure steam turbine
  
These examples show the capabilities AxialOpt and you can use them as a template to start your own projects.


## License
AxialOpt is licensed under the terms of the MIT license. See the [license file](LICENSE.md) for more information.


## Contact information
AxialOpt was developed by PhD candidate [Roberto Agromayor](https://www.ntnu.edu/employees/roberto.agromayor) and Associate Professor [Lars Olof Nord](https://www.ntnu.edu/employees/lars.nord) at the [Process and Power Research Group](https://www.ntnu.edu/ept/process-power#/view/about) of the [Norwegian University of Science and Technology (NTNU)](https://www.ntnu.no/)

Please drop us an email to [roberto.agromayor@ntnu.no](mailto:roberto.agromayor@ntnu.no) if you have questions about the code or you have a bug to report. We would also love to hear about your experiences with AxialOpt in general.



## References
<a name="1"></a>
R. Agromayor and L. O. Nord, Preliminary Design and Optimization of Axial Turbines Accounting for Diffuser Performance, International Journal of Turbomachinery, Propulsion and Power (submitted).

[![DOI Preliminary_Design_and_Optimization_of_Axial_Turbines_Accounting_for_Diffuser_Performance](https://img.shields.io/badge/DOI-Preliminary_Design_and_Optimization_of_Axial_Turbines_Accounting_for_Diffuser_Performance-blue.svg)](https://www.google.com/)

<a name="2"></a>
R. Agromayor, and L. O. Nord, AxialOpt - A Mean-Line Model for the Design and Optimization of Axial Turbines, Zenodo repository, 2019

[![DOI AxialOpt_Zenodo_repository](https://img.shields.io/badge/DOI-AxialOpt_Zenodo_repository-blue.svg)](http://doi.10.5281/zenodo.2616406)


<a name="3"></a>
E. W. Lemmon, M. L. Huber, and M. O. McLinden, NIST Standard Reference Database 23: Reference Fluid Thermodynamic and Transport Properties (REFPROP) version 9.1, National Institute of Standards and Technology, 2013.

[[![DOI REFPROP_fluid_library](https://img.shields.io/badge/DOI-REFPROP_fluid_library-blue.svg)](https://dx.doi.org/10.18434/T4JS3C)]



<a name="4"></a>
Ainley, D.G.; Mathieson, C.R. A Method of Performance Estimation for Axial-Flow Turbines. Aeronautical Research Council Reports and Memoranda 1951, 2974, 1–30.

[![URL Ainley_and_Mathieson_loss_model_(1951)](https://img.shields.io/badge/DOI-Ainley_and_Mathieson_loss_model_(1951)-blue.svg)]

<a name="5"></a>
Dunham, J.; Came, P.M. Improvements to the Ainley-MathiesonMethod of Turbine Performance Prediction. Journal of Engineering for Power 1970, 92, 252–256.

[![DOI Dunham_and_Came_loss_model_(1970)](https://img.shields.io/badge/DOI-Dunham_and_Came_loss_model_(1970)-blue.svg)]


<a name="6"></a>
Kacker, S.C.; Okapuu, U. A Mean Line Prediction Method for Axial Flow Turbine Efficiency. Journal of Engineering for Power 1982, 104, 111–119.

[![DOI Kacker_and_Okapuu_loss_model_(1982)](https://img.shields.io/badge/DOI-Kacker_and_Okapuu_loss_model_(1982)-blue.svg)]


<a name="7"></a>
Craig, H.R.M.; Cox, H.J.A. Performance Estimation of Axial Flow Turbines. Proceedings of the Institution of Mechanical Engineers 1971, 185, 407–424.

[![DOI Craig_and_Cox_loss_model_(1970)](https://img.shields.io/badge/DOI-Craig_and_Cox_loss_model_(1970)-blue.svg)](https://doi.org/10.1243/PIME_PROC_1970_185_048_02)

