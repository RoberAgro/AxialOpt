function turbine_data = evaluate_turbine_model(x,fixed_parameters)

% This function contains the mathematical model of the axial turbine

% 1) Load data and do preliminary computations
% 2) Evaluate the cascade model (once for each cascade)
% 3) Evaluate the loss model (once for each cascade)
% 4) Evaluate the diffuser model (at the exit of the last rotor)
% 5) Store data and do miscellaneous computations


%% Load parameters
angle_in   = fixed_parameters.design_input.angle_in;
t_cl       = fixed_parameters.design_input.t_cl;
fluid      = fixed_parameters.thermodynamic_input.fluid;
mass_flow  = fixed_parameters.thermodynamic_input.mass_flow;
h0_in      = fixed_parameters.thermodynamic_input.h0_in;
p0_in      = fixed_parameters.thermodynamic_input.p0_in;
p_out      = fixed_parameters.thermodynamic_input.p_out;
n_stages   = fixed_parameters.design_input.n_stages;
n_cascades = 2*n_stages;


%% Create the turbine_data structure
% All the computation results are stored in this structure
% See the function create_data_structure() for more information
turbine_data = create_data_structure(fixed_parameters.design_input.n_stages);


%% Convert the vector of independent variables to physical variables
% Isentropic velocity (also known as spouting velocity)
s_in = prop_calculation('S','H',h0_in,'P',p0_in,fluid);
h_out_s = prop_calculation('H','P',p_out,'S',s_in,fluid);
v_0 = sqrt(2*(h0_in-h_out_s));

% Isentropic Mach number at the inlet
a_in = compute_speed_of_sound(p0_in,'H',h0_in,fluid);
Ma_0 = v_0/a_in;

% Rename and scale the degrees of freedom
w_s     = x(1);
d_s     = x(2);
vel_in  = x(3)*v_0;
vel_out = x(4+0*n_cascades:4+1*n_cascades-1)*v_0;
ang_out = x(4+1*n_cascades:4+2*n_cascades-1);
r_Hc    = x(4+2*n_cascades:4+3*n_cascades-1);
r_sc    = x(4+3*n_cascades:4+4*n_cascades-1);
r_to    = x(4+4*n_cascades:4+5*n_cascades-1);
s_out   = x(4+5*n_cascades:4+6*n_cascades-1)*s_in;


%% Compute the mean radius and the angular speed
d_out_s = prop_calculation('D','P',p_out,'S',s_in,fluid);
radius = d_s/2*(mass_flow/d_out_s)^(1/2)/(h0_in-h_out_s)^(1/4);
omega = w_s*(h0_in-h_out_s)^(3/4)/(mass_flow/d_out_s)^(1/2);


%% Store some parameters in the turbine data structure for further use
% Overall (global) parameters
turbine_data.overall.fluid                = fluid;                         % Working fluid
turbine_data.overall.n_stages             = n_stages;                      % Number of stages
turbine_data.overall.n_cascades           = n_cascades;                    % Number of cascades
turbine_data.overall.mass_flow            = mass_flow;                     % Mass flow rate
turbine_data.overall.h0_in                = h0_in;                         % Stagnation enthalpy at the inlet
turbine_data.overall.p0_in                = p0_in;                         % Stagnation pressure at the inlet
turbine_data.overall.p_out                = p_out;                         % Static pressure at the outlet
turbine_data.overall.s_in                 = s_in;                          % Entropy at the inlet
turbine_data.overall.h_out_s              = h_out_s;                       % Enthalpy at the outlet for a isentropic expansion
turbine_data.overall.angular_speed        = omega;                         % Angular speed
turbine_data.overall.mean_radius          = radius;                        % Mean radius
turbine_data.overall.mean_diameter        = 2*radius;                      % Mean diameter
turbine_data.overall.specific_speed       = w_s;                           % Specific speed
turbine_data.overall.specific_diameter    = d_s;                           % Speficic diameter
turbine_data.overall.specific_blade_speed = w_s*d_s/2;                     % Specific blade velocity
turbine_data.overall.t_cl                 = t_cl;                          % Manufacturing limit for the clearance gap

% Working fluid properties
turbine_data.fluid_properties.fluid = fluid;
turbine_data.fluid_properties.MW = prop_calculation('MOLARMASS',fluid);
turbine_data.fluid_properties.p_trip = prop_calculation('P_TRIPLE',fluid);
turbine_data.fluid_properties.T_trip = prop_calculation('T_TRIPLE',fluid);
turbine_data.fluid_properties.p_crit = prop_calculation('P_CRITICAL',fluid);
turbine_data.fluid_properties.T_crit = prop_calculation('T_CRITICAL',fluid);
% turbine_data.fluid_properties.p_min = [];
% turbine_data.fluid_properties.p_max = [];
turbine_data.fluid_properties.p_min = 0.5*p_out;
turbine_data.fluid_properties.p_max = 1.5*p0_in;


%% Computate the main variables of each cascade
for k = 1:n_cascades
    
    % Define the degrees of freedom of the current cascade
    x_cascade = [vel_out(k), ...
                 ang_out(k), ...
                 r_Hc(k), ...
                 r_sc(k), ...
                 r_to(k), ...
                 s_out(k)]';

    % Compute the velocity triangles, thermodynamic properties and
    % geometric parameters of the current cascade
    turbine_data = evaluate_cascade_model(k, ...                           % Current cascade number. Stator (k is odd) or rotor (k is even)
                                          vel_in, ...                      % Absolute velocity at the inlet of the first stator (degree of freedom)
                                          angle_in, ...                    % Absolute angle at the inlet of the first stator (fixed parameter)
                                          x_cascade, ...                   % Degrees of freedom of the current cascade
                                          turbine_data);                   % Turbine data structure where all the variables are stored

    % Compute the loss coefficients of the current cascade
    turbine_data = evaluate_loss_model(k,turbine_data,fixed_parameters);
    
end


%% Compute the design parameters of each stage
% This section of the code is not critical
for k = 1:n_stages
    
    % Import variables from the cascade computation
    v_m_2 = turbine_data.plane(3+4*(k-1)).v_m;
    v_m_3 = turbine_data.plane(4+4*(k-1)).v_m;
    u_2 = turbine_data.plane(3+4*(k-1)).u;
    u_3 = turbine_data.plane(4+4*(k-1)).u;
    h1 = turbine_data.plane(1+4*(k-1)).h;
    h2 = turbine_data.plane(3+4*(k-1)).h;
    h3 = turbine_data.plane(4+4*(k-1)).h;
    h01 = turbine_data.plane(1+4*(k-1)).h0;
    h03 = turbine_data.plane(4+4*(k-1)).h0;
    
    % Compute the design parameters
    phi_in = v_m_2/u_2;                                                    % Flow coefficient at the inlet of the rotor
    phi_out = v_m_3/u_3;                                                   % Flow coefficient at the outlet of the rotor
    psi = (h01-h03)/u_3^2;                                                 % Work coefficient of the stage
    R = (h2-h3)/(h1-h3);                                                   % Degree of reaction of the stage
    
    % Save the design parameters in the turbine_data structure
    turbine_data.stage(k).phi_in = phi_in;
    turbine_data.stage(k).phi_out = phi_out;
    turbine_data.stage(k).psi = psi;
    turbine_data.stage(k).R = R;    
    
end
 

%% Compute the flow in the diffuser
if strcmp(fixed_parameters.design_input.diffuser_model, 'no') ~= 1
    turbine_data = evaluate_diffuser_model(turbine_data,fixed_parameters);
    turbine_data.diffuser.diffuser_model = fixed_parameters.design_input.diffuser_model;
else
    turbine_data.diffuser.diffuser_model = fixed_parameters.design_input.diffuser_model;
end


%% Compute the turbine global parameters
% Total-to-static and total-to-total efficiencies
h_out = turbine_data.plane(end).h;
h0_out = turbine_data.plane(end).h0;
vel_out = turbine_data.plane(end).v;
eta_ts = (h0_in-h0_out)/(h0_in-h_out_s);
eta_tt = (h0_in-h0_out)/(h0_in-h_out_s-vel_out^2/2);

% Pressure ratio (total to static)
PR = p0_in/p_out;

% Volume ratio across the turbine (total to static)
d_in = prop_calculation('D','P',p0_in,'H',h0_in,fluid);
d_out = turbine_data.plane(end).d;
VR = d_in/d_out;

% Volumetric flows at inlet and outlet
vol_flow_in = mass_flow/d_in;
vol_flow_out = mass_flow/d_out;

% Power extracted from the fluid
W_dot = mass_flow*(h0_in-h0_out);

% Revolutions per minute
N = omega*60/(2*pi);

% Torque done by the fluid on the rotor
T = W_dot/omega;

% Compute the axial length of the turbine
% Axial length of the blades plus the cascade spacing length (except for
% the cascade spacing length of the last rotor)
L = sum([turbine_data.cascade.b])+sum([turbine_data.cascade.cs])-turbine_data.cascade(end).cs;

% Volume of the turbine casing
% Computed as the volume of a cone frustum
radius   = turbine_data.overall.mean_radius;
H_in  = turbine_data.plane(1).H;
H_out = turbine_data.plane(end).H;
R1 = radius+H_in/2;
R2 = radius+H_out/2;
Vol = 1/3*pi*L*(R1^2+R1*R2+R2^2);

% Axial thrust done by the rotor on the bearings
% Computed from a control volume analysis (not yet implemented)
thrust = [];

% Store overall (global) parameters in the turbine data structure
turbine_data.overall.T_out                  = turbine_data.plane(end).T;
turbine_data.overall.T0_in                  = turbine_data.plane(1).T0;
turbine_data.overall.pressure_ratio         = PR;
turbine_data.overall.volume_ratio           = VR;
turbine_data.overall.vol_flow_in            = vol_flow_in;
turbine_data.overall.vol_flow_out           = vol_flow_out;
turbine_data.overall.eta_ts                 = eta_ts;
turbine_data.overall.eta_tt                 = eta_tt;
turbine_data.overall.power                  = W_dot;
turbine_data.overall.revs                   = N;
turbine_data.overall.torque                 = T;
% turbine_data.overall.Cp_compressible        = turbine_data.diffuser.Cp_compressible;
% turbine_data.overall.area_ratio             = turbine_data.diffuser.area_ratio;
turbine_data.overall.axial_length           = L;
turbine_data.overall.volume                 = Vol;
turbine_data.overall.thrust                 = thrust;
turbine_data.overall.Ma_rel_max             = max([turbine_data.plane.Ma_rel]);
turbine_data.overall.isentropic_velocity    = v_0;
turbine_data.overall.isentropic_mach        = Ma_0;


%% Breakdown of losses in the turbine (loss_matrix)
% Compute correction factor "a" such that sum(dh_s)*a = dh_ss
% dh_s is h_out-h_out_s in each cascade
% dh_ss is h_out-h_out_ss for the whole expansion

if strcmp(fixed_parameters.design_input.diffuser_model, 'no') ~= 1
    
    dh_s = 0;
    for k = 1:n_cascades
        dh_s = dh_s + turbine_data.cascade(k).dh_s;
    end
    dh_s = dh_s + turbine_data.diffuser.dh_s;
    dh_ss = h_out-h_out_s;
    a = dh_ss/dh_s;

    % Loss breakdown in each cascade
    % Split among profile, secondary, clearance, and trailing edge losses
    loss_matrix = zeros(n_cascades+1,6);
    for k = 1:n_cascades
        dh_s = turbine_data.cascade(k).dh_s;
        eta_drop = a*dh_s/(h0_in-h_out_s);
        fraction = [turbine_data.cascade(k).fraction];
        loss_matrix(k,1:4) = fraction*eta_drop;    
        turbine_data.cascade(k).eta_drop = eta_drop;
        turbine_data.overall.loss_matrix = loss_matrix;
    end

    % Loss breakdown in the diffuser
    eta_drop_friction = a*turbine_data.diffuser.dh_s/(h0_in-h_out_s);
    eta_drop_kinetic = (1/2*vel_out^2)/(h0_in-h_out_s);
    loss_matrix(n_cascades+1,5) = eta_drop_friction;
    loss_matrix(n_cascades+1,6) = eta_drop_kinetic;
    turbine_data.overall.loss_matrix = loss_matrix;
    turbine_data.diffuser.eta_drop_friction = eta_drop_friction;
    turbine_data.diffuser.eta_drop_kinetic = eta_drop_kinetic;
    
else
    
    dh_s = 0;
    for k = 1:n_cascades
        dh_s = dh_s + turbine_data.cascade(k).dh_s;
    end
    dh_ss = h_out-h_out_s;
    a = dh_ss/dh_s;

    % Loss breakdown in each cascade
    % Split among profile, secondary, clearance, and trailing edge losses
    loss_matrix = zeros(n_cascades+1,5);
    for k = 1:n_cascades
        dh_s = turbine_data.cascade(k).dh_s;
        eta_drop = a*dh_s/(h0_in-h_out_s);
        fraction = [turbine_data.cascade(k).fraction];
        loss_matrix(k,1:4) = fraction*eta_drop;    
        turbine_data.cascade(k).eta_drop = eta_drop;
        turbine_data.overall.loss_matrix = loss_matrix;
    end

    % Loss breakdown at the exit
    eta_drop_kinetic = (1/2*vel_out^2)/(h0_in-h_out_s);
    loss_matrix(n_cascades+1,5) = eta_drop_kinetic;
    turbine_data.overall.loss_matrix = loss_matrix;
    
end

% Check the loss breakdown is consistent: eta_ts +sum(delta_eta) = 1
% disp(1-sum(sum(loss_matrix))-eta_ts)


%% Prepare a structure to summarize the optimization problem variables
% Initialize cells
header    = {};
name      = {};
satisfied = {};
active    = {};

% Define the header names
header(1) = {'Variable name'};
header(2) = {'Lower bound'};
header(3) = {'Value'};
header(4) = {'Upper bound'};
header(5) = {'Satisfied'};
header(6) = {'Active'};


% Define variable names
k = 1;
name(k) = {'Specific speed'};    k = k + 1;
name(k) = {'Specific diameter'}; k = k + 1;
name(k) = {'Inlet velocity'};    k = k + 1;
name(k:k+n_cascades-1) = {'Outlet velocity'};        k = k+n_cascades;
name(k:k+n_cascades-1) = {'Exit flow angle'};        k = k+n_cascades;
name(k:k+n_cascades-1) = {'Aspect ratio'};           k = k+n_cascades;
name(k:k+n_cascades-1) = {'Pitch-chord ratio'};      k = k+n_cascades;
name(k:k+n_cascades-1) = {'Trailing-opening ratio'}; k = k+n_cascades;
name(k:k+n_cascades-1) = {'Outlet entropy'};

% Define the lower bounds
value_min = [fixed_parameters.lower_bounds.w_s,     ...
             fixed_parameters.lower_bounds.d_s,     ...
             fixed_parameters.lower_bounds.vel_in,  ...
             fixed_parameters.lower_bounds.vel_out, ...
             fixed_parameters.lower_bounds.ang_out, ...
             fixed_parameters.lower_bounds.r_Hc,    ...
             fixed_parameters.lower_bounds.r_sc,    ...
             fixed_parameters.lower_bounds.r_to,    ...
             fixed_parameters.lower_bounds.s_out];
          
% Define the upper bounds
value_max = [fixed_parameters.upper_bounds.w_s,     ...
             fixed_parameters.upper_bounds.d_s,     ...
             fixed_parameters.upper_bounds.vel_in,  ...
             fixed_parameters.upper_bounds.vel_out, ...
             fixed_parameters.upper_bounds.ang_out, ...
             fixed_parameters.upper_bounds.r_Hc,    ...
             fixed_parameters.upper_bounds.r_sc,    ...
             fixed_parameters.upper_bounds.r_to,    ...
             fixed_parameters.upper_bounds.s_out];
                      
% Check if the bounds are satisfied
N_var = length(value_min);
for k = 1:N_var
    if x(k)-value_min(k) >= 0 && value_max(k) - x(k) >= 0
        satisfied(k) = {'yes'};
    else
        satisfied(k) = {'no'};
    end
end

% Check if the bounds are active
for k = 1:N_var
    if abs(x(k)-value_min(k)) <= 1e-3 || abs(value_max(k) - x(k)) <= 1e-3
        active(k) = {'yes'};
    else
        active(k) = {'no'};
    end
end

% Store cells as a structure
bounds_summary.header    = header;
bounds_summary.name      = name;
bounds_summary.value     = num2cell(x);
bounds_summary.value_min = num2cell(value_min);
bounds_summary.value_max = num2cell(value_max);
bounds_summary.active    = active;
bounds_summary.satisfied = satisfied;

% Store the optimization information into turbine_data structure
turbine_data.optimization.x = x;
turbine_data.optimization.lb = value_min;
turbine_data.optimization.ub = value_max;
turbine_data.optimization.bounds_summary = bounds_summary;


%% Evaluate and store the nonlinear constraints
turbine_data = evaluate_optimization_constraints(turbine_data,fixed_parameters);


end








