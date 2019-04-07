%% Initial definitions
% Choose the loss system
% Ainley-Mathieson -> 'AM'
% Dunham-Came      -> 'DC'
% Kacker Okapuu    -> 'KO'
loss_system = 'KO';

% Choose the definition of the loss coefficient
% See Denton (1993) for a discussion about the different definitions
% Stagnation pressure loss coefficient -> 'p0'
% Energy loss coefficient              -> 'h'
% Entropy loss coefficient             -> 's'
loss_coefficient = 'p0';

% Objective function
% 'total-to-static' or 'total-to-total'
obj_function = 'total-to-static';


%% Define turbine model input fixed parameters
% Number of stages
n_stages = 4;                                                              % Any integer number
n_cascades = 2*n_stages;                                                   % Each stage has one stator cascade and one rotor cascade

% Thermodynamic boundaries (obtained from a system analysis)
% fluid = 'air.mix';                                                         % Name used by REFPROP
% The code also works with mixtures but the computational time required by
% REFPROP is too high. This may be addressed in the future if look-up
% tables are implemented to compute the thermodynamic properties
fluid = 'nitrogen';                                                        % Name used by REFPROP
T0_in = 1533;                                                              % Inlet stagnation temperature in (K) (alternatively provide the stagnation enthalpy if the expansion starts in the two-phase region)
p0_in = 14.6*100;                                                          % Inlet stagnation pressure in (kPa)
h0_in = refpropm('h','T',T0_in,'p',p0_in,fluid);                           % Inlet stagnation enthalpy in (J/kg)   
p_out = 1.00*100;                                                          % Outlet static pressure in (Pa)

% Define the mass flow rate
mass_flow = 435.0;                                                        % Mass flow rate (kg/s)


% Alternatively it is possible to provide the isentropic power output
% isentropic_power = 250e6;                                                  % Isentropic power in (W)
% s_in = refpropm('s','T',T0_in,'p',p0_in,fluid);                            % Inlet entropy in (K/kg K)
% h_out_s = refpropm('h','p',p_out,'s',s_in,fluid);                          % Outlet static isentropic enthalpy (J/kg)
% mass_flow = power_isentropic/(h0_in-h_out_s);                              % Mass flow rate (kg/s)

% Other parameters
t_cl = 5e-4;                                                               % Tip clearance height manufacturing limit (m). Values between 0.2-0.5 mm are reasonable                  
angle_in = 0.00*pi/180;                                                    % First stator inlet angle (rad). Using axial entry stators (algle_in = 0.00) is reasonable                                      
phi = 30*pi/180;                                                           % Diffuser mean cant angle (rad)
div = 5*pi/180;                                                            % Diffuser divergence semiangle (rad)
AR = 3.000;                                                                % Diffuser area ratio (set a value close to one for the no-diffuser case)
Cf = 0.010;                                                                % Skin friction coefficient. Dubitsky–Japikse (2008) suggest 0.010 as a reasonable estimate

                       
%% Define the independent variables and optimization bounds
% 1) Specify a an initial guess for the degrees of freedom
% 2) Specify the lower and upper bounds for the degrees of freedom
% Specific speed
w_s = 1.00;                                                                % Easier to guess than the angular speed
w_s_min = 0.1;                                                             % The range 0.1-10.0 covers almost all cases
w_s_max = 10.0;                                                            % This range can be extended if desired

% Specific diameter
d_s = 2/sqrt(n_stages)/w_s;                                                % Easier to guess than the actual diameter. 2/sqrt(n_stages)/w_s
d_s_min = 0.1;                                                             % The range 0.1-10.0 covers almost all cases
d_s_max = 10.0;                                                            % This range can be extended if desired

% Reduced velocity at the inlet of the first stator
vel_in = 0.20;
vel_in_min = 0.001;                                                        % The reduced velocity is positive
vel_in_max = 0.500;                                                        % The reduced velocity is lower than 1.00

% Reduced relative velocity at the outlet of each cascade
vel_out(1:2:n_cascades-1) = 1/sqrt(2*n_stages);                            % Reduced relative velocity at the outlet of the stators. 1/sqrt(2*n_stages); is a good initial guess
vel_out(2:2:n_cascades)   = 1/sqrt(2*n_stages);                            % Reduced relative velocity at the outlet of the rotors. 1/sqrt(2*n_stages); is a good initial guess
vel_out_min(1:n_cascades) = 0.05;                                          % The reduced velocity is positive
vel_out_max(1:n_cascades) = 1.25;                                          % The reduced velocity is lower than 1.00

% Relative angle at the outlet of each cascade
ang_out(1:2:n_cascades-1) = +70/180*pi;                                    % Relative angle at the outlet of each stator
ang_out(2:2:n_cascades)   = -70/180*pi;                                    % Relative angle at the outlet of each rotor
ang_out_min(1:2:n_cascades-1) = +40/180*pi;                                % The low limit of the Ainley-Mathieson profile loss correlation for the relative angle at the outlet of the stator is +40 deg (tricky sign convention)
ang_out_min(2:2:n_cascades)   = -80/180*pi;                                % The low limit of the Ainley-Mathieson profile loss correlation for the relative angle at the outlet of the rotor is -80 deg (tricky sign convention)
ang_out_max(1:2:n_cascades-1) = +80/180*pi;                                % The high limit of the Ainley-Mathieson profile loss correlation for the relative angle at the outlet of the stator is +80 deg (tricky sign convention)
ang_out_max(2:2:n_cascades)   = -40/180*pi;                                % The high limit of the Ainley-Mathieson profile loss correlation for the relative angle at the outlet of the rotor is -40 deg (tricky sign convention)

% Aspect ratio
r_Hc(1:n_cascades) = 1.30;                                                 % Blade height to actual chord aspect ratio (not axial chord)
r_Hc_min(1:n_cascades) = 1.000;                                            % Saravanamuttoo advises against aspect ratios lower than 1.00
r_Hc_max(1:n_cascades) = 2.000;                                            % Saravanamuttoo indicates that high values may leed to vibration problems. Values around 3.00-4.00 are safe

% Pitch to chord ratio
r_sc(1:n_cascades) = 0.75;                                                 % Spacing to chord ratio (inverse of solidity in US terminology)
r_sc_min(1:n_cascades) = 0.30;                                             % The low limit of the Ainley-Mathieson profile loss correlation for the pitch to chord ratio is 0.30                              
r_sc_max(1:n_cascades) = 1.10;                                             % The high limit of the Ainley-Mathieson profile loss correlation for the pitch to chord ratio is 1.10                             

% Entropy at the outlet of each cascade
% Entropy is a degree of freedom but it is constrained by the loss model 
% Assume reasonable value for the total-to-static efficiency and compute
% the entropy at the outlet of each cascade assuming a linear distribution
% from the inlet to the outlet of the turbine
eta = 0.85;                                                                % Total-to-static isentropic efficiecy guess
s_in  = refpropm('s','h',h0_in,'p',p0_in,fluid);
h_out_s = refpropm('h','p',p_out,'s',s_in,fluid);
h_out = h0_in-eta*(h0_in-h_out_s);
s = refpropm('s','p',p_out,'h',h_out,fluid);
delta_s = (s-s_in)/n_cascades;
s_out  = (s_in+delta_s:delta_s:s)/s_in;                                    % Entropy at the outlet of each cascade (J/kg K). Linear increase of entropy along the turbine
s_out_min(1:n_cascades) = 0.99;                                            % The entropy is cannot decrease along the turbine
eta = 0.60;                                                                % Minimum isentropic efficiency used to give an upper bound for the entropy at the outlet of each cascade (s_out_max)
h_out = h0_in-eta*(h0_in-h_out_s);
s_out_max(1:n_cascades) = refpropm('s','p',p_out,'h',h_out,fluid)/s_in;    % Conservative upper bound for the entropy



%% Design constraints
% 1) Specify the minimum and maximum values (use brackets ignore the bound)
% 2) Specify whether to apply the constraint or not
% 3) Specify a reference value to scale the problem 

% Flaring angle
constraints.flare_angle.min = -10*pi/180;                                  % Ainley and Mathieson (March 1951) recommend a maximum of 12.5 deg
constraints.flare_angle.max = +10*pi/180;                                  % Ainley and Mathieson (March 1951) recommend a maximum of 12.5 deg
constraints.flare_angle.applied = 'yes';
constraints.flare_angle.ref = 1.00;

% Hub-to-tip ratio
constraints.r_ht.min = 0.600;                                              % Kacker-Okapuu correlation holds for values as low as 0.500. Saravanamuttoo discusses ratios as low as 0.50.
constraints.r_ht.max = 0.900;                                              % High values imply large secondary losses
constraints.r_ht.applied = 'yes';
constraints.r_ht.ref = 1.00;

% Degree of reaction 
constraints.reaction.min = 0.100;                                          % Avoid pressure increase in rotor cascades  (R>0)
constraints.reaction.max = 0.900;                                          % Avoid pressure increase in stator cascades (R<1)
constraints.reaction.applied = 'yes';
constraints.reaction.ref = 1.00;

% Relative Mach number
constraints.Ma_rel.min = [];
constraints.Ma_rel.max = 5.000;                                            % Set an arbitrary limit for the maximum relative Mach number if desired
constraints.Ma_rel.applied = 'yes';
constraints.Ma_rel.ref = 1.00;

% Diffuser inlet Mach number
constraints.Ma_diffuser.min = [];
constraints.Ma_diffuser.max = 0.95;                                        % Constraint the diffuser to be subsonic (prevent the diffuser to behave as a supersonic nozzle)
constraints.Ma_diffuser.applied = 'yes';
constraints.Ma_diffuser.ref = 1.00;

% Constrain inlet angles (avoid too low deflection)
constraints.beta_in_stator.min = [];
constraints.beta_in_stator.max = +15/180*pi;                               % Kacker and Okapuu (1982) propose a correlation with a maximum value of +30 deg for stator cascades
constraints.beta_in_stator.applied = 'yes';
constraints.beta_in_stator.ref = 1.00;

% Constrain inlet angles (avoid too low deflection)
constraints.beta_in_rotor.min = -15/180*pi;                                % Kacker and Okapuu (1982) propose a correlation with a minimum value of -30 deg for rotor cascades 
constraints.beta_in_rotor.max = [];
constraints.beta_in_rotor.applied = 'yes';
constraints.beta_in_rotor.ref = 1.00;

% Minimum height constraint
constraints.height.min = 0.01;                                             % Insert any value in [m]
constraints.height.max = [];
constraints.height.applied = 'no';
constraints.height.ref = 0.01;
                                                                    
% Minimum chord constraint
constraints.chord.min = 0.01;                                              % Insert any value in [m]
constraints.chord.max = [];
constraints.chord.applied = 'no';
constraints.chord.ref = 0.01;

% Angular speed constraint
constraints.RPM.min = 3000;                                                % Insert any value in [rpm]
constraints.RPM.max = 3000;                                                % 3000 or 3600 for syncronous speed in Europe or USA, respectively
constraints.RPM.applied = 'yes';
constraints.RPM.ref = 10000;

% Feel free to add additional constraints making small modifications to:
% 1) The constraints. structure defined in this section of the code
% 2) The function AxialOpt_computation_constraints() of the source code


%% Working fluid properties
% Critical properties
% These are necessary to avoid flash computations when the pressure or
% temperature are supercritical
p_crit = refpropm('p','C',0,' ',0,fluid);                                  % Critical pressure
T_crit = refpropm('T','C',0,' ',0,fluid);                                  % Critical temperature

% Pressure and temperature limits
% These limits are used by the function p_hs_flash() in case the ordinary
% enthalpy-entropy function call to REFPROP fails
p_min = 0.5*p_out;                                                         % Minimum pressure
p_max = 2.0*p0_in;                                                         % Maximum pressure
p_0 = 101.325;                                                             % Ambient pressure

% Molecular weight
MW = refpropm('M','T',T0_in,'p',p0_in,fluid);                                  % Molecular weight


%% Organize the turbine parameters in a structure
% Thermodynamic input from the cycle optimization    
thermodynamic_input = struct('fluid',       fluid,       ...
                             'mass_flow',   mass_flow,   ...
                             'h0_in',       h0_in,       ...
                             'p0_in',       p0_in,       ...
                             'p_out',       p_out);
         
% Thermodynamic properties of the working fluid                        
fluid_properties = struct('fluid',     fluid,    ...
                          'T_crit',    T_crit,   ...
                          'p_crit',    p_crit,   ...
                          'p_min',     p_min,    ...
                          'p_max',     p_max,    ...
                          'MW',        MW);
                         
% Turbine design specifications                         
design_input = struct('n_stages',         n_stages,         ...
                      'n_cascades',       n_cascades,       ...
                      'loss_system',      loss_system,      ...
                      'loss_coefficient', loss_coefficient, ...
                      'obj_function',     obj_function,     ...
                      't_cl',             t_cl,             ...
                      'angle_in',         angle_in,         ...
                      'phi',              phi,              ...
                      'div',              div,              ...
                      'AR',               AR,               ...
                      'Cf',               Cf);
       
% Lower bounds for the degrees of freedom                 
lower_bounds = struct('w_s',       w_s_min,       ...
                      'd_s',       d_s_min,       ...
                      'vel_in',    vel_in_min,    ...
                      'vel_out',   vel_out_min,   ...
                      'ang_out',   ang_out_min,   ...
                      'r_Hc',      r_Hc_min,      ...
                      'r_sc',      r_sc_min,      ...
                      's_out',     s_out_min);

% Upper bounds for the degrees of freedom     
upper_bounds = struct('w_s',       w_s_max,       ...
                      'd_s',       d_s_max,       ...
                      'vel_in',    vel_in_max,    ...
                      'vel_out',   vel_out_max,   ...
                      'ang_out',   ang_out_max,   ...
                      'r_Hc',      r_Hc_max,      ...
                      'r_sc',      r_sc_max,      ...
                      's_out',     s_out_max);
                  

% Organize the previous structures in a compact way     
parameters_turbine = struct('thermodynamic_input', thermodynamic_input, ...
                            'fluid_properties',    fluid_properties,    ...
                            'design_input',        design_input,        ...
                            'lower_bounds',        lower_bounds,        ...
                            'upper_bounds',        upper_bounds,        ...
                            'constraints',         constraints);
                    
                        
                                        
                        
%% Create the optimization_problem structure for fmincon
optimization_problem = struct('lb',        [], ...                         % Vector of lower bounds
                              'ub',        [], ...                         % Vector of upper bounds
                              'Aeq',       [], ...                         % Matrix for linear equality constraints
                              'beq',       [], ...                         % Vector for linear equality constraints
                              'Aineq',     [], ...                         % Matrix for linear inequality constraints
                              'bineq',     [], ...                         % Vector for linear inequality constraints
                              'nonlcon',   [], ...                         % Nonlinear constraints function
                              'objective', [], ...                         % Objective function
                              'x0',        [], ...                         % Initial guess for the degrees of freedom
                              'options',   [], ...                         % Optimization options structure
                              'solver', 'fmincon');

% Define the initial guess vector
x0_turb = [w_s,       ...
           d_s,       ...
           vel_in,    ...
           vel_out,   ...
           ang_out,   ...
           r_Hc,      ...
           r_sc,      ...
           s_out]';

optimization_problem.x0 = x0_turb;

% Define the vector of lower bounds
lb_turb = [w_s_min,       ...
           d_s_min,       ...
           vel_in_min,    ...
           vel_out_min,   ...
           ang_out_min,   ...
           r_Hc_min,      ...
           r_sc_min,      ...
           s_out_min]'; 

optimization_problem.lb = lb_turb;
                       
% Define the vector of upper bounds                  
ub_turb = [w_s_max,       ...
           d_s_max,       ...
           vel_in_max,    ...
           vel_out_max,   ...
           ang_out_max,   ...
           r_Hc_max,      ...
           r_sc_max,      ...
           s_out_max]';
       
optimization_problem.ub = ub_turb;
    

% Define the options using the optimoptions function
% 'interior-point', 'active-set', and 'sqp' algorithms
% 'active-set' is my favourite because it is usually a bit faster
% 'spq' can be more robust but takes a bit longer to converge
% My experience with 'interior-point' is not so good
algorithm = 'sqp';


% Parallel computing
% Might fail in some cases (probably relate to h-s flash function calls)
% true or false
parallel = false;


% Create options structure
optimization_problem.options = optimoptions(@fmincon,                   ...
                               'Display', 'iter-detailed',              ...
                               'Algorithm', algorithm   ,               ...
                               'StepTolerance', 1e-06,                  ...
                               'ConstraintTolerance', 1e-06,            ...
                               'OptimalityTolerance', 1e-06,            ...
                               'FunctionTolerance', 1e-6,               ...
                               'MaxFunctionEvaluations', 5000,          ...
                               'MaxIterations', 1000,                   ...
                               'FiniteDifferenceStepSize', sqrt(eps),   ...
                               'UseParallel', parallel,                 ...
                               'PlotFcn', {@optimplotfval,              ...
                                           @optimplotconstrviolation,   ...
                                           @optimplotfirstorderopt,     ...
                                           @optimplotstepsize});
              
                       
% Description of the possible outputs of the optimization function
% Exitflag =  1 -> Success. First-order optimality measure was less than options and constraints are not violated.
% Exitflag =  2 -> Success. Step size was less than options and constraints are not violated.
% Exitflag =  0 -> Unsuccess. Number of iterations or function evaluations exceeded option
% Exitflag = -1 -> Unsuccess. Solver stopped by an output function
% Exitflag = -2 -> Unsuccess. No feasible point was found
% Only active set algorithm
% Exitflag =  4 -> Success. Magnitude of the search direction was less than 2*options.StepTolerance and maximum constraint violation was less than options.ConstraintTolerance.
% Exitflag =  5 -> Success. Magnitude of directional derivative in search direction was less than 2*options.OptimalityTolerance and maximum constraint violation was less than options.ConstraintTolerance.

