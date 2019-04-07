# AxialOpt
AxialOpt is a code for the preliminary design and optimization of axial turbines. The output of AxialOpt can be used in system-level analyses (such as a power cycle optimization) to estimate the efficiency or footprint of axial turbines for a given set of thermodynamic specifications. In addition, the information provided by AxialOpt can be used as the starting point for the aerodynamic design of the turbine blades using more advanced flow models based on CFD.

The models and optimization methodology of AxialOpt are documented in a peer-reviewed, open-access publication [1](#1) and the source code is also stored in a Zenodo repository [2].

## Features

* The axial turbine model is composed of three sub-models that are used as building blocks:
  1. A meanline model to describe the flow in each cascade.
  2. An empirical loss model to evaluate the entropy generation in each cascade.
  3. A one-dimensional flow model to account for the influence of the diffuser.
* The model is formulated for axial turbines with any number of stages.
* The model is formulated to use arbitrary equations of state to compute thermodynamic properties:
  1. The current version uses the REFPROP v9.1 fluid library.
  2. Other fluid libraries and look-up tables may be implemented in the future.
* The loss model is formulated in a general way to use:
  1. Any set of empirical correlations to compute the losses:
      1. Ainley and Mathieson (implemented).
      2. Dunhan and Came (implemented).
      3. Kacker and Okapuu (implemented).
      4. Craig and Cox (will be implemented soon).
      5. Other loss model contributions are welcome.
  2. Different definitions for the loss coefficient:
      1. Stagnation pressure loss coefficient.
      2. Enthalpy loss coefficient (also known as kinetic energy loss coefficient).
      3. Entropy loss coefficient.
* The preliminary design is formulated as a constrained optimization problem.
  1. This allows explore the design space in a systematic way and account for technical constraints.
  2. It is straighforward to modify the objective function and constraints depending on the problem.
  3. There are several gradient-based algorithms available to solve the optimization problem, including:
      1. Sequential Quadratic Programming (SQP).
      2. Interior Point (IP).
* The output of AxialOpt can be saved as:
  1. Text file containing a summary of the solution of the optimization problem.
  2. Text files containing the geometry of each cascade and variables at each station.
  3. Figures with the T-s and h-s diagrams of the expansion.
  4. Figures with the the velocity triangles.
  5. Figures with the axial-radial and axial-tangential views of the turbine.
  6. Figures with the breakdown of the losses within the turbine.

_Note: AxialOpt is not suitable to estimate the performance of an existing design under different conditions.
  The extension of the code to compute the performance at off-design conditions is underway._



## Requisites
AxialOpt was implemented in MATLAB R2018a [3] and requires a REFPROP v9.1 installation [4]. The folder [refprop-matlab](AnnularDiffuser1D/refprop-matlab) contains a short guide about how to link REFPROP with MATLAB.

AxialOpt has the option to use the export_fig library to produce publication-quality figures. Using the export_fig library requires _ghostcript_ and _pdftops_. See the installation instructions in the original repository.



## Examples
The folder [AxialOpt_examples](AxialOpt_examples)  AxialOpt_examples contains five commented examples to get started with the code:
  * A supercritical Carbon dioxide turbine.
  * A organic Rankine cycle (ORC) turbine using R125 as working fluid.
  * A organic Rankine cycle (ORC) turbine using hexane as working fluid.
  * An industrial gas turbine.
  * A high-pressure steam turbine.
These examples show the capabilities AxialOpt and you can use them as a template to start your own projects!

## Contributions

    _Note: Feel free to contact if you have implemented another loss model and would like to contribute to this repository_  


## License
AxialOpt is licensed under the terms of the MIT license. See the [license file](LICENSE.md) for more information.


## Contact information
AxialOpt was developed by PhD candidate [Roberto Agromayor](https://www.ntnu.edu/employees/roberto.agromayor) and Associate Professor [Lars Olof Nord](https://www.ntnu.edu/employees/lars.nord) at the [Thermal Energy Research Group](https://www.ntnu.edu/ept/thermal-energy1
) of the [Norwegian University of Science and Technology (NTNU)](https://www.ntnu.no/)

Contact the email address [roberto.agromayor@ntnu.no](mailto:roberto.agromayor@ntnu.no) for inquiries about the code.


Please feel free to submit any bugs or feature requests. We'd also love to hear about your experiences with autograd in general. Drop us an email!

Contributions and bug reports are welcome!



## References
<a name="1"></a>[1] R. Agromayor and L. O. Nord, Preliminary Design and Optimization of Axial Turbines Accounting for Diffuser Performance, International Journal of Turbomachinery, Propulsion and Power (submitted).

[![DOI](https://img.shields.io/badge/DOI-Diffuser_paper_DOI-blue.svg)](https://www.google.com) (not ready yet)


[2] R. Agromayor, and L. O. Nord, AxialOpt - A Mean-Line Model for the Design and Optimization of Axial Turbines, Zenodo repository, 2019

[![DOI](https://zenodo.org/badge/178391900.svg)](https://zenodo.org/badge/latestdoi/178391900)


[3] The MathWorks Inc., MATLAB version R2018a, 2018.

[![URL](https://img.shields.io/badge/URL-https://nl.mathworks.com/-blue.svg)](https://nl.mathworks.com/)


[4] E. W. Lemmon, M. L. Huber, and M. O. McLinden, NIST Standard Reference Database 23: Reference Fluid Thermodynamic and Transport Properties (REFPROP) version 9.1, National Institute of Standards and Technology, 2013.

[![DOI](https://img.shields.io/badge/DOI-https://dx.doi.org/10.18434/T4JS3C-blue.svg)](https://dx.doi.org/10.18434/T4JS3C)



