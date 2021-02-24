function turbine_data = evaluate_diffuser_model(turbine_data,fixed_parameters)

% This function contains the diffuser flow model

% If diffuser_model = '1D_flow'
% Solve the ODE system describing the flow in a vaneless annular channel
% 1) steady
% 2) compressible
% 3) axisymmetric flow
% 4) with heat transfer and friction
%

%% Input parameters
% Design parameters
Cf = fixed_parameters.design_input.Cf;          % Skin friction coefficient
phi = fixed_parameters.design_input.phi;        % Mean wall cant angle
div = fixed_parameters.design_input.div;        % Wall divergence semi-angle
AR_target = fixed_parameters.design_input.AR;   % Objective area ratio
diffuser_model = fixed_parameters.design_input.diffuser_model;

% Outlet plane index
k = 4*fixed_parameters.design_input.n_stages;

% Inlet state
fluid = turbine_data.overall.fluid;
p_min = turbine_data.fluid_properties.p_min;
p_max = turbine_data.fluid_properties.p_max;
d_in = turbine_data.plane(k).d;
p_in = turbine_data.plane(k).p;
p0_in = turbine_data.plane(k).p0;
h_in = turbine_data.plane(k).h;
h0_in = turbine_data.plane(k).h0;
s_in = turbine_data.plane(k).s;
mass_flow = turbine_data.overall.mass_flow;

d0_in = prop_calculation('D','P',p0_in,'S',s_in,fluid);

% Inlet velocity
v_m_in = turbine_data.plane(k).v_m;
v_t_in = turbine_data.plane(k).v_t;

% Inlet geometry
r_in = turbine_data.plane(k).r_m;
H_in = turbine_data.plane(k).H;
b_in = H_in/cos(phi);
x_in = 0.00;


%% Compute the geometry of the diffuser
% Solve quadratic equation to find the exit meridional coordinate
a0 = -r_in*b_in*(AR_target-1);
a1 = 2*tan(div)*r_in + sin(phi)*b_in;
a2 = 2*tan(div)*sin(phi);
m_out = roots([a2, a1, a0]);    % Find the two solutions
m_out = m_out(m_out>0);         % Pick the positive solution

    function error = compute_energy_balance_error(density)
        v_t = (r_in/r_out)*v_t_in;
        v_m = mass_flow/(2*pi*r_in*r_out*density);
        h = prop_calculation('H','D',density,'S',s_in,fluid);
        error = h0_in - (h + v_t^2/2 + v_m^2/2);
    end

%% Evaluate the diffuser model
if strcmp(diffuser_model, 'isentropic')
    
    % Compute exit geometry
    r_out = r_fun(r_in,phi,m_out);
    x_out = x_fun(x_in,phi,m_out);
    b_out = b_fun(b_in,div,m_out);
    
    % Compute the exit density
    options = optimset('TolX',1e-6);
    d_out = fzero(@compute_energy_balance_error, [0.75*d_in, 1.25*d0_in], options);

    % Compute the exit state
    p_out = prop_calculation('P','D',d_out,'S',s_in,fluid);
    h_out = prop_calculation('H','D',d_out,'S',s_in,fluid);
    v_out = sqrt(2*(h0_in - h_out));
    v_t_out = r_in/r_out*v_t_in;
    v_m_out = sqrt(v_out^2 - v_t_out^2);

    % Organize information
    v_m = [v_m_in v_m_out];
    v_t = [v_t_in v_t_out];
    d = [d_in d_out];
    p = [p_in p_out];
    h = [h_in h_in];
    h0 = [h0_in h0_in];
    s = [s_in s_in];
    s_gen = [s_in s_in];
    m = [0 m_out];


elseif strcmp(diffuser_model, '1D')
    
    % Define the initial conditions
    U0 = [v_m_in v_t_in d_in p_in s_in];

    % % Use the options RelTol and AbsTol to set the integration tolerance
    % options = odeset('Events',@(m,U)area_ratio(m,U,AR_target,phi,div,r_in,b_in), ...
    %                  'RelTol',1e-4, ...
    %                  'AbsTol',1e-4);
    % 
    % % Integrate the ode system using ode45
    % % Use m and U and variables and the rest of inputs as extra parameters
    % % Integrate between 0 and infinity (the AR is the stopping criterion)
    % [m,U] = ode45(@(m,U)ode_diffuser(m,U,phi,div,r_in,b_in,x_in,Cf,fluid),[0,inf],U0,options);

    % Use the options RelTol and AbsTol to set the integration tolerance
    options = odeset('RelTol',1e-4, 'AbsTol',1e-4);
    
    % Integrate the ode system using ode45
    % Use m and U and variables and the rest of inputs as extra parameters
    % Integrate between 0 and infinity (the AR is the stopping criterion)
    [m,U] = ode45(@(m,U)ode_diffuser(m,U,phi,div,r_in,b_in,x_in,Cf,fluid),[0,m_out],U0,options);

    % Rename the solution to the other calculations
    v_m = U(:,1);
    v_t = U(:,2);
    d = U(:,3);
    p = U(:,4);
    s_gen = U(:,5);

    % Compute the entropy and stagnation enthalpy using the equations of state
    [N, ~] = size(s_gen);
    s = zeros(N,1);
    h  = zeros(N,1);
    h0 = zeros(N,1);
    for i = 1:size(s)
        s(i) = prop_calculation('S','P',p(i),'D',d(i),fluid);
        h(i) = prop_calculation('H','P',p(i),'D',d(i),fluid);
        h0(i) = h(i) +(v_m(i)^2 + v_t(i)^2)/2;
    end

else
    error("Choose a valid option for diffuser_model: '1D', 'isentropic', or 'no'")
end


% Compute the geometry of the diffuser
r = r_fun(r_in,phi,m);
x = x_fun(x_in,phi,m);
b = b_fun(b_in,div,m);
A = 2*pi*r.*b;
AR = (b.*r)/(b_in*r_in);
x_outer = x - b/2.*sin(phi);
x_inner = x + b/2.*sin(phi);
r_outer = r + b/2.*cos(phi);
r_inner = r - b/2.*cos(phi);
phi_2 = phi + div;
phi_1 = phi - div;

%% Compute the variables at the outlet of the diffuser
% Velocities
v_m_out = v_m(end);
v_t_out = v_t(end);
v_out = sqrt(v_m_out^2+v_t_out^2);
alpha_out = atan(v_t_out/v_m_out);
w_m_out = v_m_out;
w_t_out = v_t_out;
w_out = v_out;
beta_out = atan(w_t_out/w_m_out);

% Thermodynamic states
d_out = d(end);
p_out = p(end);
T_out = prop_calculation('T','P',p_out,'D',d_out,fluid);
h_out = prop_calculation('H','P',p_out,'D',d_out,fluid);
s_out = prop_calculation('S','P',p_out,'D',d_out,fluid);
Z_out = prop_calculation('Z','P',p_out,'D',d_out,fluid);
h_out_s = prop_calculation('H','P',p_out,'S',s_in,fluid);
a_out = compute_speed_of_sound(p_out,'H',h_out,fluid);
mu_out = compute_viscosity(p_out,'H',h_out,fluid);
h0_out = h_out+v_out^2/2;
h0rel_out = h0_out;
p0_out = p_hs_flash(h0_out,s_out,fluid,p_min,p_max);
p0rel_out = p0_out;
T0_out = prop_calculation('T','P',p0_out,'H',h0_out,fluid);
T0rel_out = T0_out;
Ma_out = v_out/a_out;
Ma_rel_out = Ma_out;


%% Compute the pressure recovery factor
% Compressible definition (general)
Cp_compressible = (p-p_in)/(p0_in-p_in);

% Incmpressible definition
v_in = sqrt(v_m_in^2+v_t_in^2);
Cp_incompressible = (p-p_in)/(1/2*d_in*v_in^2);

% Ideal pressure recovery coefficient
tan_a = v_t_in(1)/v_m_in(1);
Cp_ideal = 1 - (r_in./r).^2.*((b_in./b).^2 + tan_a^2)/(1 + tan_a^2);


%% Store the computed variables in the turbine_data structure
% Save the variables at the outlet plane
k = k+1;
turbine_data.plane(k).v      = v_out;
turbine_data.plane(k).v_t    = v_t_out;
turbine_data.plane(k).v_m    = v_m_out;
turbine_data.plane(k).w      = w_out;
turbine_data.plane(k).w_t    = w_t_out;
turbine_data.plane(k).w_m    = w_m_out;
turbine_data.plane(k).u      = 0;
turbine_data.plane(k).alpha  = alpha_out;
turbine_data.plane(k).beta   = beta_out;
turbine_data.plane(k).T      = T_out;
turbine_data.plane(k).T0     = T0_out;
turbine_data.plane(k).T0rel  = T0rel_out;
turbine_data.plane(k).p      = p_out;
turbine_data.plane(k).p0     = p0_out;
turbine_data.plane(k).p0rel  = p0rel_out;
turbine_data.plane(k).d      = d_out;
turbine_data.plane(k).h      = h_out;
turbine_data.plane(k).h0     = h0_out;
turbine_data.plane(k).h0rel  = h0rel_out;
turbine_data.plane(k).s      = s_out;
turbine_data.plane(k).Z      = Z_out;
turbine_data.plane(k).a      = a_out;
turbine_data.plane(k).mu     = mu_out;
turbine_data.plane(k).Ma     = Ma_out;
turbine_data.plane(k).Ma_rel = Ma_rel_out;
turbine_data.plane(k).A      = A(end);
turbine_data.plane(k).H      = b(end);
turbine_data.plane(k).r_m    = r(end);
turbine_data.plane(k).r_h    = r_inner(end);
turbine_data.plane(k).r_t    = r_outer(end);
turbine_data.plane(k).r_ht   = r_inner(end)/r_outer(end);

% Save the main diffuser parameters
turbine_data.diffuser.area_ratio = AR(end);
turbine_data.diffuser.axial_length = x(end);
turbine_data.diffuser.phi = phi;
turbine_data.diffuser.div = div;
turbine_data.diffuser.phi_1 = phi_1;
turbine_data.diffuser.phi_2 = phi_2;
turbine_data.diffuser.Cf = Cf;
turbine_data.diffuser.Cp_compressible = Cp_compressible(end);
turbine_data.diffuser.Cp_incompressible = Cp_incompressible(end);
turbine_data.diffuser.Cp_ideal = Cp_ideal(end);
turbine_data.diffuser.r_1_in = r_inner(1);
turbine_data.diffuser.r_2_in = r_outer(1);
turbine_data.diffuser.r_1_out = r_inner(end);
turbine_data.diffuser.r_2_out = r_outer(end);
turbine_data.diffuser.b_in = b(1);
turbine_data.diffuser.b_out = b(end);
turbine_data.diffuser.A_in = A(1);
turbine_data.diffuser.A_out = A(end);
turbine_data.diffuser.dh_s = h_out-h_out_s;

% Store parameter distribution
turbine_data.diffuser.ODE_solution.AR = AR;
turbine_data.diffuser.ODE_solution.Cp = Cp_compressible;
turbine_data.diffuser.ODE_solution.v_m = v_m;
turbine_data.diffuser.ODE_solution.v_t = v_t;
turbine_data.diffuser.ODE_solution.p = p;
turbine_data.diffuser.ODE_solution.d = d;
turbine_data.diffuser.ODE_solution.s_gen = s_gen;
turbine_data.diffuser.ODE_solution.s = s;
turbine_data.diffuser.ODE_solution.h = h;
turbine_data.diffuser.ODE_solution.h0 = h0;
turbine_data.diffuser.ODE_solution.m = m;
turbine_data.diffuser.ODE_solution.r = r;
turbine_data.diffuser.ODE_solution.x = x;
turbine_data.diffuser.ODE_solution.b = b;
turbine_data.diffuser.ODE_solution.x_outer = x_outer;
turbine_data.diffuser.ODE_solution.x_inner = x_inner;
turbine_data.diffuser.ODE_solution.r_outer = r_outer;
turbine_data.diffuser.ODE_solution.r_inner = r_inner;

end


function dUdm = ode_diffuser(m,U,phi,div,r_in,b_in,x_in,Cf,fluid)

% Rename variables
v_m = U(1);
v_t = U(2);
d   = U(3);
p   = U(4);
alpha = atan(v_t/v_m);
v = sqrt(v_m^2+v_t^2);

% Local geometry
r = r_fun(r_in,phi,m);              % Radius as a function of m
x = x_fun(x_in,phi,m);              % Axial distance as a function of m
b = b_fun(b_in,div,m);              % Channel width as a function of m

% Derivative of the area change (forward finite differences)
delta = 1e-4;
diff_br = (b_fun(b_in,div,m+delta)*r_fun(r_in,phi,m+delta) - b*r)/delta;

% Derivative of internal energy with respect to pressure (constant density)
% e1 = prop_calculation('U','P',p-delta,'D',d,fluid);
% e2 = prop_calculation('U','P',p+delta,'D',d,fluid);
% dedp_d = (e2 - e1)/(2*delta);
dedp_d = prop_calculation('d(U)/d(P)|D','P',p,'D',d,fluid);

% % Ideal gas limit (check)
% cp = prop_calculation('CPMASS','P',p,'D',d,fluid);
% cv = prop_calculation('CVMASS','P',p,'D',d,fluid);
% dedp_d_ideal = 1/(d*(cp/cv-1));

% Speed of sound (avoid computations in the two phase region)
a = compute_speed_of_sound(p,'D',d,fluid);

% Stress at the wall
tau_w = Cf*d*v^2/2;        % Skin friction coefficient

% Heat flux at the wall
q_w = 0;                   % Adiabatic wall

% Coefficient matrix A
A = [d           0          v_m       0;
     d*v_m       0            0       1;
     0       d*v_m            0       0;
     0           0   -d*v_m*a^2   d*v_m];

% Source term vector
S = zeros(4,1);
S(1) = -d*v_m/(b*r)*diff_br;
S(2) = +d*v_t*v_t/r*sin(phi) - 2*tau_w/b*cos(alpha);
S(3) = -d*v_t*v_m/r*sin(phi) - 2*tau_w/b*sin(alpha);
S(4) = 2*(tau_w*v + q_w)/b/dedp_d;

% Obtain the slope of the solution by Gaussian elimination
dUdm = A\S;

% Check entropy generation
T = prop_calculation('T','D',d,'P',p,fluid);
sigma = 2/b*(tau_w*v);
dUdm(5) = sigma/(d*v_m)/T;      % ds/dm

end


function r = r_fun(r_in,phi,m)
r = r_in + sin(phi)*m;
end


function x = x_fun(x_in,phi,m)
x = x_in + cos(phi)*m;
end


function b = b_fun(b_in,div,m)
b = b_in + 2*tan(div)*m;
end


function [AR_check,isterminal,direction] = area_ratio(m,~,AR_prescribed,phi,div,r_in,b_in)

% Geometry
r = r_fun(r_in,phi,m);              % Radius as a function of m
b = b_fun(b_in,div,m);              % Channel width as a function of m
AR_current = (b*r)/(b_in*r_in);     % Current area ratio

% Stopping criterion
AR_check = AR_prescribed - AR_current;
isterminal = 1;   % stop the integration
direction = 0;    % negative direction

end



%% Extensions of the code to make it more general:
% Use general functions for the geometry as input (instead of linear funcs)
% If the geometry is provided as general functions the local angle phi has
% to be computed by differentiation (finite differences)
% Perhaps prescribe the diffuser geometry using NURBS curves

% If the parametrization variable is not the meridional coordinate (m) then
% it is necessary to integrate the arclength to compute m

% Provide an arbitrary variation of the skin friction coefficient or use an
% empirical correlation to compute it

% Implement the heat transfer model in the code
% Perhaps use the Chilton-Colburn analogy to compute the heat transfer
% coefficient. It is also necessary to compute the stagnation temperature
% and to prescribe a wall temperature distribution



%% Compute the diffuser efficiency drop due to friction (deprecated)
% h01 = turbine_data.overall.h0_in;                                          % Stagnation enthalpy at the inlet of the turbine
% s1  = turbine_data.overall.s_in;                                           % Entropy at the inlet of the turbine
% h2s = turbine_data.overall.h_out_s;                                        % Static enthalpy at the outlet of the turbine for an isentropic expansion
% h_in_ss = prop_calculation('H','P',p_in,'S',s1,fluid);                             % Static enthalpy at the inlet of the current cascade computed using the entropy at the inlet of the turbine
% h_out_ss = prop_calculation('H','P',p_out,'S',s1,fluid);                           % Static enthalpy at the outlet of the current cascade computed using the entropy at the inlet of the turbine
% eta_drop_friction = ((h_in_ss-h_out_ss)-(h_in-h_out))/(h01-h2s);           % Efficiency drop due to skin friction (total-to-static efficiency)
% eta_drop_kinetic = (v_out^2/2)/(h01-h2s);                                  % Efficiency drop due to discharge kinetic energy (total-to-static efficiency)


