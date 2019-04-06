# AxialOpt
AxialOpt contains a mean-line model and optimization methodology for the preliminary design of axial turbines.



mean-line model for the design and optimization of axial turbines with any number of stages allows to use arbitrary equations of state and loss models. The model is documented in a peer-reviewed, open-access publication [1] and the source code is also stored in a Zenodo repository [2].


## Features
* The preliminary design is formulated as a constrained optimization problem.
* The optimization problem can be solved using various gradient based algorithms, including:
  * Sequential Quadratic Programming (SQP)
  * Interior Point (IP)

* The code was formulated so that it is straighforward to modify the objective function and constraints depending on the problem.

* The turbine model is formulated to use arbitrary equations of state to compute thermodynamic properties.
  * The current version uses the REFPROP v9.1 fluid library
  * Other fluid libraries and look-up tables may be implemented in the future

* The code can be used to design turbine with any number of stages

* The turbine model accounts for the influence of the diffuser using a one-dimensional flow model

* The turbine model is formulated in a general way so that it is possible to use any empirical loss model
   * Ainley and Mathieson (implemented)
   * Dunhan and Came (implemented)
   * Kacker and Okapuu (implemented)
   * Craig and Cox (will be implemented soon)

_Note: Feel free to contact us if you implement a new loss model and would like to add it to this repository._

* There are several definitions available for the loss coefficient
  * Stagnation pressure loss coefficient
  * Enthalpy loss coefficient (also known as kinetic energy loss coefficient)
  * Entropy loss coefficient
  
  
* The code is not suitable to estimate the performance of an existing design (off-design simulations)
_Note: The extension of the code to compute the performance at off-design conditions is underway_

* The results of the optimization are saved in 


* Plotting features:
  * Ts diagram of the expansion (with saturation line and isobars)
  * Ts diagram of the expansion (with saturation line and isobars)
  * Two options to plot the velocity triangles of each stage
  * Axial-radial view of the turbine (with or without diffuser)
  * Axial-tangential view of the turbine
  * Breakdown of the losses (profile, secondary, tip clearance, trailing edge, diffuser friction, and kinetic energy)


Different possibilities to save the figures


* Saving features
  * The results of the optimization are saved in the form of .txt files and also as a .mat file containing a MATLAB structure with the
  * Overall parameters
  * Geometry of the cascades
  * Thermodynamic states
  * Velocity triangles
  * Information about the optimization problem solution (independent variables and constraints)


## Examples
The folder xxxx contains five examples that illustrate the capabilities of the code, including

  * Gas turbine
  * Steam turbine
  * Organic Rankine cycle turbines (saturated and transcritical)
  * Supercritical Carbon dioxide turbine

These examples illustrate that AxialOpt can be used to optimize turbines for a wide variety of applications and they

can be used as a template for your own project


## Requisites
AxialOpt was implemented in MATLAB R2018a [3] and requires a REFPROP v9.1 installation [4] (see instructions).

Saving the figures using the export_fig library requitres pdf tops and ghoscript. See installation instructions in the export_fig repository

## Instructions
* Just run the script [`AxialOpt_demo.m`](AxialOpt/Examples/AxialOpt_demo.m) to get started!
* The folder [refprop-matlab](AnnularDiffuser1D/refprop-matlab) contains a short guide about how to link REFPROP with MATLAB.


## License
AxialOpt is licensed under the terms of the MIT license. See the [license file](LICENSE.md) for more information.


## Contact information
AxialOpt was developed by PhD candidate [Roberto Agromayor](https://www.ntnu.edu/employees/roberto.agromayor) and Associate Professor [Lars Olof Nord](https://www.ntnu.edu/employees/lars.nord) at the [Thermal Energy Research Group](https://www.ntnu.edu/ept/thermal-energy1
) of the [Norwegian University of Science and Technology (NTNU)](https://www.ntnu.no/)

Contact the email address [roberto.agromayor@ntnu.no](mailto:roberto.agromayor@ntnu.no) for inquiries about the code.


## References
[1] R. Agromayor and L. O. Nord, Preliminary Design and Optimization of Axial Turbines Accounting for Diffuser Performance, International Journal of Turbomachinery, Propulsion and Power (submitted).

[![DOI](https://img.shields.io/badge/DOI-Diffuser_paper_DOI-blue.svg)](https://www.google.com) (not ready yet)


[2] R. Agromayor, and L. O. Nord, AxialOpt - A Mean-Line Model for the Design and Optimization of Axial Turbines, Zenodo repository, 2019

[![DOI](https://zenodo.org/badge/178391900.svg)](https://zenodo.org/badge/latestdoi/178391900)


[3] The MathWorks Inc., MATLAB version R2018a, 2018.

[![URL](https://img.shields.io/badge/URL-https://nl.mathworks.com/-blue.svg)](https://nl.mathworks.com/)


[4] E. W. Lemmon, M. L. Huber, and M. O. McLinden, NIST Standard Reference Database 23: Reference Fluid Thermodynamic and Transport Properties (REFPROP) version 9.1, National Institute of Standards and Technology, 2013.

[![DOI](https://img.shields.io/badge/DOI-https://dx.doi.org/10.18434/T4JS3C-blue.svg)](https://dx.doi.org/10.18434/T4JS3C)



