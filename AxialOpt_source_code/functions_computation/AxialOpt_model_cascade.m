function turbine_data = AxialOpt_model_cascade(k,vel_in,angle_in,x_cascade,turbine_data)
%% Evaluate the flow in a turbine cascade
% Author: Roberto Agromayor

% 1) Velocity triangles (in and out)
% 2) Thermodynamic conditions (in and out)
% 3) Cascade geometry


%% Determine if the cascade is an stator or rotor
if mod(k,2) == 0
    type = 'rotor';             % k is odd for stator cascades
else
    type = 'stator';            % k is even for rotor cascades
end


%% Rename the degrees of freedom
vel_out   = x_cascade(1);
angle_out = x_cascade(2);
r_Hc      = x_cascade(3);
r_sc      = x_cascade(4);
s_out     = x_cascade(5);


%% Rename the fixed variables for convenience
omega = turbine_data.overall.angular_speed;
radius = turbine_data.overall.mean_radius;
m = turbine_data.overall.mass_flow;
fluid = turbine_data.overall.fluid;
p_crit = turbine_data.fluid_properties.p_crit;
p_min = turbine_data.fluid_properties.p_min;
p_max = turbine_data.fluid_properties.p_max;


%% Velocity triangle at the inlet
% The inlet velocity triangle is defined by the inlet absolute velocity
if k == 1
    v_in = vel_in;                               % Degree of freedom
    alpha_in = angle_in;                         % Fixed parameter
elseif k > 1
    v_in     = turbine_data.plane(2*k-2).v;      % Previous cascade outlet
    alpha_in = turbine_data.plane(2*k-2).alpha;  % Previous cascade outlet
end

% Define the blade speed for stator and rotor cascades
if strcmp(type,'stator') == 1
    u_in = 0.00;
elseif strcmp(type,'rotor') == 1
    u_in = omega*radius;
else
    error('The cascade must be a stator or a rotor')
end

% Compute the inlet velocity triangle
[u_in,v_in,v_t_in,v_m_in,w_in,w_t_in,w_m_in,alpha_in,beta_in] = ...
    AxialOpt_compute_triangle_in(u_in,v_in,alpha_in);


%% Velocity triangle at the outlet
% The outlet velocity triangle is defined by the outlet relative velocity
w_out = vel_out;         % Degree of freedom
beta_out = angle_out;    % Degree of freedom

% Define the blade speed for stator and rotor cascades
if strcmp(type,'stator') == 1
    u_out = 0.00;
elseif strcmp(type,'rotor') == 1
    u_out = omega*radius;
else
    error('The cascade must be a stator or a rotor')
end

% Compute the outlet velocity triangle
[u_out,v_out,v_t_out,v_m_out,w_out,w_t_out,w_m_out,alpha_out,beta_out] =...
    AxialOpt_compute_triangle_out(u_out,w_out,beta_out);


%% Compute the thermodynamic properties
% Thermodynamic properties at the inlet of the cascade
if k == 1
    
    % Thermodynamic properties at the inlet of the first stator
    h0_in = turbine_data.overall.h0_in;
    p0_in = turbine_data.overall.p0_in;
    T0_in = refpropm('T','p',p0_in,'h',h0_in,fluid);
    s_in  = refpropm('s','p',p0_in,'h',h0_in,fluid);
    h_in  = h0_in-v_in^2/2;
    p_in  = p_hs_flash(h_in,s_in,fluid,p_min,p_max);
    T_in  = refpropm('T','p',p_in,'s',s_in,fluid);
    d_in = refpropm('d','p',p_in,'h',h_in,fluid);
    Z_in = refpropm('z','p',p_in,'h',h_in,fluid);

    q_in = quality('h',h_in,fluid,p_in,p_crit);
    % Compute the other thermodynamic properties avoiding computations in 
    % the two phase region that REFPROP is not able to handle
    if q_in > 1.00
        % Ordinary computations if the fluid is in the gas phase
        a_in = refpropm('a','p',p_in,'h',h_in,fluid);
        mu_in = refpropm('V','p',p_in,'h',h_in,fluid);
    elseif q_in <= 1.00
        % Use saturated vapor properties is the fluid is in the 2-phase region
        a_in = refpropm('a','p',p_in,'q',1.00,fluid);
        mu_in = refpropm('V','p',p_in,'q',1.00,fluid);
    end
    
    h0rel_in = h0_in;
    p0rel_in = p0_in;
    T0rel_in = T0_in;
    
elseif k > 1
    
    % Thermodynamic properties at the inlet of the current cascade
    % Static properties are equal at the inlet of cascade k and and the
    % outlet of cascade k-1
    T_in     = turbine_data.plane(2*k-2).T;
    p_in     = turbine_data.plane(2*k-2).p;
    d_in     = turbine_data.plane(2*k-2).d;
    h_in     = turbine_data.plane(2*k-2).h;
    s_in     = turbine_data.plane(2*k-2).s;
    Z_in     = turbine_data.plane(2*k-2).Z;
    a_in     = turbine_data.plane(2*k-2).a;
    mu_in    = turbine_data.plane(2*k-2).mu;
    
    % Absolute stagnation properties are equal at the inlet of cascade k 
    % and the outlet of cascade k-1    
    h0_in    = turbine_data.plane(2*k-2).h0;
    p0_in    = turbine_data.plane(2*k-2).p0;
    T0_in    = turbine_data.plane(2*k-2).T0;
    
    % Relative stagnation properties change from the inlet of cascade k 
    % to the outlet of cascade k-1    
    h0rel_in = h_in+w_in^2/2;
    p0rel_in = p_hs_flash(h0rel_in,s_in,fluid,p_min,p_max);
    T0rel_in = refpropm('T','p',p0rel_in,'h',h0rel_in,fluid);

end

% Thermodynamic properties at the outlet of the cascade
% Use the conservation of rothalpy (or relative staganation enthalpy) to
% compute the enthalpy at the outlet. Then use the entropy at the outlet
% (one of the degrees of freedom) to determine the thermodynamic state
h0rel_out = h0rel_in;   % This equation works for both rotor and stator
h_out = h0rel_out-w_out^2/2;
h0_out = h_out+v_out^2/2;
p0_out = p_hs_flash(h0_out,s_out,fluid,p_min,p_max);
T0_out = refpropm('T','p',p0_out,'h',h0_out,fluid);
p0rel_out = p_hs_flash(h0rel_out,s_out,fluid,p_min,p_max);
T0rel_out = refpropm('T','p',p0rel_out,'h',h0rel_out,fluid);
p_out = p_hs_flash(h_out,s_out,fluid,p_min,p_max);  % s_out is a DoF
T_out = refpropm('T','p',p_out,'h',h_out,fluid);
d_out = refpropm('d','p',p_out,'h',h_out,fluid);
Z_out = refpropm('z','p',p_out,'h',h_out,fluid);
h_out_s = refpropm('h','p',p_out,'s',s_in,fluid);

% Compute the other thermodynamic properties avoiding computations in the
% two phase region that REFPROP is not able to handle
q_out = quality('h',h_out,fluid,p_out,p_crit);
if q_out > 1.00
    % Ordinary computations if the fluid is in the gas phase
    a_out = refpropm('a','p',p_out,'h',h_out,fluid);
    mu_out = refpropm('V','p',p_out,'h',h_out,fluid);
elseif q_out <= 1.00
    % Use saturated vapor properties is the fluid is in the 2-phase region
    a_out = refpropm('a','p',p_out,'q',1.00,fluid);
    mu_out = refpropm('V','p',p_out,'q',1.00,fluid);
end



%% Compute the Mach numbers
% Mach number at the inlet of the cascade
Ma_in = v_in/a_in;
Ma_rel_in = w_in/a_in;

% Mach number at the outlet of the cascade
Ma_out = v_out/a_out;
Ma_rel_out = w_out/a_out;


%% Compute the geometry of the stage
% Geometry at the inlet of the cascade
A_in = m/(d_in*v_m_in);
H_in = A_in/(2*pi*radius);
r_h_in = radius-H_in/2;
r_t_in = radius+H_in/2;
r_ht_in = r_h_in/r_t_in;

% Geometry at the outlet of the cascade
A_out = m/(d_out*v_m_out);
H_out = A_out/(2*pi*radius);
r_h_out = radius-H_out/2;
r_t_out = radius+H_out/2;
r_ht_out = r_h_out/r_t_out;


%% Compute the cascade geometry
% Compute incidence, deviation and metal angle
% These lines could be the staring point for off-design computations
incid = 0.00;
delta = 0.00;
theta_in = beta_in-incid;         
theta_out = beta_out-delta;
camber     = theta_out-theta_in;
deflection = beta_out-beta_in;

% Compute other geometric parameters
H        = (H_in+H_out)/2;                                                 % Mean height
c        = (1/r_Hc)*H;                                                     % Chord
s        = r_sc*c;                                                         % Spacing
o        = s*cos(beta_out);                                                % Opening
t_max    = geometry_correlation_thickness(abs(theta_in-theta_out),c);      % Maximum thickness
t_te     = o/10;                                                           % Trailing edge thickness (Macchi 1981 use 0.1. Seems reasonable based on Kacker-Okapuu correlation. Ainley Mathieson proposed 0.02*s)
t_cl     = turbine_data.overall.t_cl;                                      % Tip and shroud clearance (Luca Da Lio 2016 suggest 0.5 mm or is it 0.2 mm?)
stagger  = (theta_in+theta_out)/2;                                         % Stagger angle (circular arc blades)
b        = c*cos(stagger);                                                 % Axial chord
delta_fl = atan((H_out-H_in)/(2*b));                                       % Flaring angle
cs       = 0.40*b;                                                         % Axial spacing after the cascade
z        = 2*pi*radius/s;                                                  % Number of blades
Re       = d_out*w_out*c/mu_out;                                           % Reynolds number (evaluated at the outlet of the stage)
% stagger  = 0.90*geometry_correlation_stagger(beta_in,beta_out,type);     % Stagger angle correlation proposed by Okapuu
% cs       = max(1.5*o,0.3*b);
% t_te     = max([0.02*s,L_min]);                                          % Trailing edge thickness
% t_cl     = max(H/100,t_cl);                                        % Tip and shroud clearance (no reference yet)


% Some comments:
%{

The opening is computed with an approximate relation (can be found in
Dixon or other texbooks. Other authors propose correlationss that would
give more approximate values. These correlations could be added in the
future. 1) Ainley and Mathieson 2) Aungier

The trailing edge thickness is computed with the rule o=0.02s proposed in
several papers including the one from Kacker and Okappu

Saravanamutto recommends ss=(0.25-0.50)b for the spacing after each
cascade. G. Persico recommends to use an spacing of the same order of
magnitude as the opening ss~o (personal conversation). This later approach
has more physical meaning because it is related to the lenght scale of the
trailing edge vortices and their dissipation

The stagger angle was computed assuming circular arc blades
It could also be computed according to the correlation propossed by Okappu
stagger  = geometry_stagger(theta_in,theta_out,type);
This approach did not work so well for me in the past

Tip and shroud clearance is computed according to the recommendations of
some papers by Lazzaretto et al. (2016)

%}


%% Store computed variables
% Save the variables at the inlet plane
turbine_data.plane(2*k-1).v      = v_in;
turbine_data.plane(2*k-1).v_t    = v_t_in;
turbine_data.plane(2*k-1).v_m    = v_m_in;
turbine_data.plane(2*k-1).w      = w_in;
turbine_data.plane(2*k-1).w_t    = w_t_in;
turbine_data.plane(2*k-1).w_m    = w_m_in;
turbine_data.plane(2*k-1).u      = u_in;
turbine_data.plane(2*k-1).alpha  = alpha_in;
turbine_data.plane(2*k-1).beta   = beta_in;
turbine_data.plane(2*k-1).T      = T_in;
turbine_data.plane(2*k-1).T0     = T0_in;
turbine_data.plane(2*k-1).T0rel  = T0rel_in;
turbine_data.plane(2*k-1).p      = p_in;
turbine_data.plane(2*k-1).p0     = p0_in;
turbine_data.plane(2*k-1).p0rel  = p0rel_in;
turbine_data.plane(2*k-1).d      = d_in;
turbine_data.plane(2*k-1).h      = h_in;
turbine_data.plane(2*k-1).h0     = h0_in;
turbine_data.plane(2*k-1).h0rel  = h0rel_in;
turbine_data.plane(2*k-1).s      = s_in;
turbine_data.plane(2*k-1).Z      = Z_in;
turbine_data.plane(2*k-1).a      = a_in;
turbine_data.plane(2*k-1).mu     = mu_in;
turbine_data.plane(2*k-1).Ma     = Ma_in;
turbine_data.plane(2*k-1).Ma_rel = Ma_rel_in;
turbine_data.plane(2*k-1).A      = A_in;
turbine_data.plane(2*k-1).H      = H_in;
turbine_data.plane(2*k-1).r_m    = radius;
turbine_data.plane(2*k-1).r_h    = r_h_in;
turbine_data.plane(2*k-1).r_t    = r_t_in;
turbine_data.plane(2*k-1).r_ht   = r_ht_in;

% Save the variables at the outlet plane
turbine_data.plane(2*k).v      = v_out;
turbine_data.plane(2*k).v_t    = v_t_out;
turbine_data.plane(2*k).v_m    = v_m_out;
turbine_data.plane(2*k).w      = w_out;
turbine_data.plane(2*k).w_t    = w_t_out;
turbine_data.plane(2*k).w_m    = w_m_out;
turbine_data.plane(2*k).u      = u_out;
turbine_data.plane(2*k).alpha  = alpha_out;
turbine_data.plane(2*k).beta   = beta_out;
turbine_data.plane(2*k).T      = T_out;
turbine_data.plane(2*k).T0     = T0_out;
turbine_data.plane(2*k).T0rel  = T0rel_out;
turbine_data.plane(2*k).p      = p_out;
turbine_data.plane(2*k).p0     = p0_out;
turbine_data.plane(2*k).p0rel  = p0rel_out;
turbine_data.plane(2*k).d      = d_out;
turbine_data.plane(2*k).h      = h_out;
turbine_data.plane(2*k).h0     = h0_out;
turbine_data.plane(2*k).h0rel  = h0rel_out;
turbine_data.plane(2*k).s      = s_out;
turbine_data.plane(2*k).Z      = Z_out;
turbine_data.plane(2*k).a      = a_out;
turbine_data.plane(2*k).mu     = mu_out;
turbine_data.plane(2*k).Ma     = Ma_out;
turbine_data.plane(2*k).Ma_rel = Ma_rel_out;
turbine_data.plane(2*k).A      = A_out;
turbine_data.plane(2*k).H      = H_out;
turbine_data.plane(2*k).r_m    = radius;
turbine_data.plane(2*k).r_h    = r_h_out;
turbine_data.plane(2*k).r_t    = r_t_out;
turbine_data.plane(2*k).r_ht   = r_ht_out;

% Save the variables of the cascade
turbine_data.cascade(k).type      = type;
turbine_data.cascade(k).r_m       = radius;
turbine_data.cascade(k).c         = c;
turbine_data.cascade(k).b         = b;
turbine_data.cascade(k).stagger   = stagger;
turbine_data.cascade(k).H         = H;
turbine_data.cascade(k).delta_fl  = delta_fl;
turbine_data.cascade(k).s         = s;
turbine_data.cascade(k).o         = o;
turbine_data.cascade(k).t_max     = t_max;
turbine_data.cascade(k).t_te      = t_te;
turbine_data.cascade(k).t_cl      = t_cl;
turbine_data.cascade(k).z         = z;
turbine_data.cascade(k).cs        = cs;
turbine_data.cascade(k).Re        = Re;
turbine_data.cascade(k).theta_in  = theta_in;
turbine_data.cascade(k).theta_out = theta_out;
turbine_data.cascade(k).incid     = incid;
turbine_data.cascade(k).delta     = delta;
turbine_data.cascade(k).camber    = camber;
turbine_data.cascade(k).deflection = deflection;
turbine_data.cascade(k).dh_s      = h_out-h_out_s;


end




%% Compute the efficiency drop in the current cascade (deprecated)
% % This computation may not seem clear, see handwritten notes
% h01 = turbine_data.overall.h0_in;                                          % Stagnation enthalpy at the inlet of the turbine
% s1  = turbine_data.overall.s_in;                                           % Entropy at the inlet of the turbine
% h2s = turbine_data.overall.h_out_s;                                        % Static enthalpy at the outlet of the turbine for an isentropic expansion
% h_in_ss = refpropm('h','p',p_in,'s',s1,fluid);                             % Static enthalpy at the inlet of the current cascade computed using the entropy at the inlet of the turbine
% h_out_ss = refpropm('h','p',p_out,'s',s1,fluid);                           % Static enthalpy at the outlet of the current cascade computed using the entropy at the inlet of the turbine
% eta_drop = ((h_in_ss-h_out_ss)-(h_in-h_out))/(h01-h2s);                    % Efficiency drop in the current cascade (total-to-static efficiency)





