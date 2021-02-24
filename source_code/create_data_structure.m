function turbine_data = create_data_structure(n_stages)

% Create a data structure to store the turbine variables

% This function is not really necessary because MATLAB does not require to
% declare variables before using them
% I created this function as a reference to describe the variables that are
% used within the code and the units of each variable

% Variables associated with cascades
var_cascade = struct('type',      [], ...                                  % Specify if the cascade is a stator or rotor
                     'r_m',       [], ...                                  % Mean radius of the cascade (m)
                     'c',         [], ...                                  % Chord (m)
                     'b',         [], ...                                  % Axial chord (m)
                     'stagger',   [], ...                                  % Stagger angle (rad)
                     'H',         [], ...                                  % Mean blade height (m)
                     'delta_fl',  [], ...                                  % Flaring angle (rad)
                     's',         [], ...                                  % Spacing between blades (also known as pitch) (m)
                     'o',         [], ...                                  % Opening (m)
                     't_max',     [], ...                                  % Maximum thickness (m)
                     't_te',      [], ...                                  % Trailing edge thickness (m)
                     't_cl',      [], ...                                  % Tip clearance thickness (m)
                     'cs',        [], ...                                  % Spacing between cascades (m)
                     'z',         [], ...                                  % Number of blades (-)
                     'Re',        [], ...                                  % Relative Reynolds number of the cascade (-)
                     'theta_in',  [], ...                                  % Metal angle (rad)
                     'theta_out', [], ...                                  % Metal angle (rad)
                     'deflection',[], ...                                  % Flow deflection (rad)
                     'camber',    [], ...                                  % Blade camber angle (rad)
                     'incid',     [], ...                                  % Incidence angle (rad)
                     'delta',     [], ...                                  % Deviation angle (rad)
                     'eta_drop',  [], ...                                  % Fraction of the isentropic efficiency drop caused by the cascade
                     'Y_guess',   [], ...                                  % Loss coefficient computed from the degrees of freedom (-)
                     'Y_loss',    [], ...                                  % Loss coefficient computed from the loss correlations (-)
                     'Y_error',   [], ...                                  % Error in the loss coefficient (-)
                     'Y_p',       [], ...                                  % Profile loss coefficient computed from the loss correlations (-)
                     'Y_s',       [], ...                                  % Secondary loss coefficient computed from the loss correlations (-)
                     'Y_cl',      [], ...                                  % Tip clearance loss coefficient computed from the loss correlations (-)
                     'Y_te',      [], ...                                  % Trailing edge loss coefficient computed from the loss correlations (-)
                     'dh_s',      [], ...                                  % Difference between the outlet static enthalpy and the local isentropic enthalpy (J/kg)
                     'fraction',  []);                                     % Share of each loss mechanism to the total loss coefficient

                 
% Variables associated with the planes at the inlet and outlet of cascades
var_plane = struct('T',      [], ...                                       % Static temperature (K)
                   'T0',     [], ...                                       % Stagnation temperature (K)
                   'T0rel',  [], ...                                       % Relative stagnation temperature (K)
                   'p',      [], ...                                       % Static pressure (kPa)
                   'p0',     [], ...                                       % Stagnation pressure (kPa)
                   'p0rel',  [], ...                                       % Relative stagnation pressure (kPa)
                   'h',      [], ...                                       % Enthalpy (J/kg)
                   'h0',     [], ...                                       % Stagnation enthalpy (J/kg)
                   'h0rel',  [], ...                                       % Relative stagnation enthalpy (J/kg)
                   'd',      [], ...                                       % Density (kg/m3)
                   's',      [], ...                                       % Entropy (J/kg.K)
                   'Z',      [], ...                                       % Compressibility factor (-)
                   'a',      [], ...                                       % Speed of sound (m/s)
                   'mu',     [], ...                                       % Dynamic viscosity (Pa.s)
                   'v',      [], ...                                       % Absolute velocity (m/s)
                   'v_t',    [], ...                                       % Tangential component of the absolute velocity (m/s)
                   'v_m',    [], ...                                       % Meridional component of the absolute velocity (m/s)
                   'w',      [], ...                                       % Relative velocity (m/s)
                   'w_t',    [], ...                                       % Tangential component of the relative velocity (m/s)
                   'w_m',    [], ...                                       % Meridional component of the relative velocity (m/s)
                   'u',      [], ...                                       % Peripherial velocity (m/s)
                   'alpha',  [], ...                                       % Absolute flow angle (rad)
                   'beta',   [], ...                                       % Relative flow angle (rad)
                   'Ma',     [], ...                                       % Absolute Mach number (-)
                   'Ma_rel', [], ...                                       % Relative Mach number (-)
                   'A',      [], ...                                       % Annulus area (m2)
                   'H',      [], ...                                       % Blade height (m)
                   'r_m',    [], ...                                       % Mean radius (m)
                   'r_h',    [], ...                                       % Hub radius (m)
                   'r_t',    [], ...                                       % Tip radius (m)
                   'r_ht',   []);                                          % Hub to tip radii ratio (-)

               
% Variables associated with a turbine stage (stator and rotor)
var_stage = struct('R',       [], ...                                      % Degree of reaction
                   'phi_in',  [], ...                                      % Flow coefficient at the rotor inlet
                   'phi_out', [], ...                                      % Flow coefficient at the rotor outlet
                   'psi',     []);                                         % Work coefficient
               
               
% Variable associated with the diffuser
diffuser = struct('area_ratio',        [], ...                             % Area ratio
                  'axial_length',      [], ...                             % Diffuser axial length
                  'Cp_compressible',   [], ...                             % Pressure recovery factor
                  'Cp_incompressible', [], ...                             % Pressure recovery factor for incompressible flow
                  'Cp_ideal',          [], ...                             % Ideal pressure recovery factor for inviscid, incompressible flowç
                  'eta_drop_friction', [], ...                             % Total to static efficiency drop due to friction in the diffuser walls
                  'eta_drop_kinetic',  [], ...                             % Total to static efficiency drop due to discharge kinetic energy 
                  'Cf',                [], ...                             % Skin friction coefficient (-)
                  'phi',               [], ...                             % Mean cant angle
                  'div',               [], ...                             % Divergence semiangle
                  'phi_1',             [], ...                             % Inner wall angle
                  'phi_2',             [], ...                             % Outer wall angle
                  'r_1_in',            [], ...                             % Hub radius at the inlet
                  'r_2_in',            [], ...                             % Casing radius at the inlet
                  'r_1_out',           [], ...                             % Hub radius at the outlet
                  'r_2_out',           [], ...                             % Casing radius at the outlet
                  'b_in',              [], ...                             % Channel height at the inlet
                  'b_out',             [], ...                             % Channer height at the outlet
                  'A_in',              [], ...                             % Flow area at the inlet
                  'A_out',             [], ...                             % Flow area at the outlet
                  'dh_s',              [], ...                             % Difference between the outlet static enthalpy and the local isentropic enthalpy (J/kg)
                  'ODE_solution',      []);                                % Structure containing the solution of the ODE system along the diffuser
                    
                                       
% Variables associated with the turbine as a whole
overall = struct('fluid',                [], ...                           % Working fluid
                 'n_stages',             [], ...                           % Number of stages
                 'n_cascades',           [], ...                           % Number of cascades
                 'p0_in',                [], ...                           % Stagnation pressure at the inlet of the turbine (kPa)
                 'T0_in',                [], ...                           % Stagnation temperature at the inlet of the turbine (K)
                 'p_out',                [], ...                           % Static pressure at the outlet of the turbine (kPa)
                 'T_out',                [], ...                           % Static temperature at the outlet of the turbine (K)
                 'pressure_ratio',       [], ...                           % Pressure ratio (stagnation to static)
                 'volume_ratio',         [], ...                           % Volume ratio (stagnation to static)
                 'mass_flow',            [], ...                           % Mass flow rate (kg/s)
                 'vol_flow_in',          [], ...                           % Volumetric flow rate at the inlet of the turbine (m3/s)
                 'vol_flow_out',         [], ...                           % Volumetric flow rate at the outlet of the turbine (m3/s)
                 'eta_ts',               [], ...                           % Total to static efficiency
                 'eta_tt',               [], ...                           % Total to total efficiency
                 'power',                [], ...                           % Power output (W)
                 'angular_speed',        [], ...                           % Angular speed (rad/s)
                 'revs',                 [], ...                           % Angular speed (rpm)
                 'torque',               [], ...                           % Mechanical torque (Nm)             
                 'mean_radius',          [], ...                           % Turbine mean radius (m)
                 'mean_diameter',        [], ...                           % Turbine mean diameter (m)
                 'axial_length',         [], ...                           % Axial length (m)
                 'volume',               [], ...                           % Turbine volume (m3)
                 'specific_speed',       [], ...                           % Specific speed (-)
                 'specific_diameter',    [], ...                           % Specific diameter (-)
                 'thrust',               [], ...                           % Axial thrust (N)
                 'Ma_rel_max',           [], ...                           % Maximum relative Mach number within the turbine
                 'isentropic_velocity',  [], ...                           % Isentropic velocity (also known as spouting velocity) (m/s)
                 'isentropic_mach',      []);                              % Isentropic Mach number (based on the speed of sound at the inlet)
              
             
% Thermodynamic properties of the working fluid                        
fluid_properties = struct('fluid',     [],   ...                           % Name of the working fluid
                          'T_crit',    [],   ...                           % Critical temperature of the working fluid
                          'p_crit',    [],   ...                           % Critical pressure of the working fluid
                          'T_trip',    [],   ...                           % Triple temperature of the working fluid
                          'p_trip',    [],   ...                           % Triple pressure of the working fluid
                          'MW',        []);                                % Molecular mass of the working fluid
             
                      
% Variables associated with the optimization problem
optimization = struct('x',                   [],   ...                     % Vector of degrees of freedom
                      'lb',                  [],   ...                     % Vector of upper bounds
                      'ub',                  [],   ...                     % Vector of lower bounds
                      'c',                   [],   ...                     % Vector of inequality constraints
                      'c_eq',                []);                          % Vector of equality constraints
           
% Create the turbine data structure
n_cascades = 2*n_stages;
cascade(1:n_cascades) = var_cascade;
plane(1:1+2*n_stages) = var_plane;       
stage(1:n_stages)     = var_stage;
turbine_data = struct('cascade',            cascade,             ...
                      'plane',              plane,               ...
                      'stage',              stage,               ...
                      'diffuser',           diffuser,            ...
                      'overall',            overall,             ...
                      'fluid_properties',   fluid_properties,    ...
                      'optimization',       optimization);
             

end