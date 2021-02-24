function turbine_data = evaluate_optimization_constraints(turbine_data,fixed_parameters)

% This function computes the equality and inequality constraint violation

%% Compute the vector of inequality constraints
% MATLAB uses constraints of the type c_ineq<0 
% Change sign for c_ineq>0 constraints

% Load variables
values_ineq.flare_angle    = [turbine_data.cascade.delta_fl]';             % Vector of flaring angles
values_ineq.r_ht           = [turbine_data.plane(1:end-1).r_ht]';          % Vector of hub to tip radii ratios (do not consider the section at the outlet of the diffuser)
values_ineq.beta_in_stator = [turbine_data.plane(1:4:end-3).beta]';        % Vector of relative angles at the inlet of the stator cascades
values_ineq.beta_in_rotor  = [turbine_data.plane(3:4:end-1).beta]' ;       % Vector of relative angles at the inlet of the rotor cascades
values_ineq.reaction       = [turbine_data.stage.R]';                      % Vector with the degree of reaction of each stage
values_ineq.Ma_rel         = [turbine_data.plane.Ma_rel]';                 % Vector of relative Mach numbers
values_ineq.height         = [turbine_data.cascade.H]';                    % Vector of cascade mean blade height
values_ineq.chord          = [turbine_data.cascade.c]';                    % Vector of cascade chords
values_ineq.thickness      = [turbine_data.cascade.t_max]';                % Vector of cascade maximum blade thickness
values_ineq.thickness_te   = [turbine_data.cascade.t_te]';                 % Vector of cascade blade trailing edge thickness
values_ineq.RPM            = turbine_data.overall.revs;                    % Angular speed in [rpm]
Ma_diffuser_in             = turbine_data.plane(end-1).Ma;                 % Mach number at the inlet of the diffuser
alpha_diffuser_in          = turbine_data.plane(end-1).alpha;              % Absolute flow angle at the inlet of the diffuser
values_ineq.Ma_diffuser    = Ma_diffuser_in*cos(alpha_diffuser_in);        % Meridional Mach number at the inlet of the diffuser

% Rename variables to make code more readable
constraints = fixed_parameters.constraints;

% Get constraint names
constraint_names = fieldnames(constraints);

% Compute the vector of inequality constraints
c_ineq = [];
for k = 1:length(constraint_names)
    
    % Determine wether to apply the constraint or not
    if strcmp(constraints.(constraint_names{k}).applied,'yes') == 1
        
        % Get the value of the current variable
        value = values_ineq.(constraint_names{k});
        
        % Get the limits of the current variable
        value_min = constraints.(constraint_names{k}).min;
        value_max = constraints.(constraint_names{k}).max;
        
        % Get the reference value of the current variable
        value_ref = constraints.(constraint_names{k}).ref;
        
        % Apply a lower limit to the current variable if relevant
        if isempty(value_max) ~= 1
            c_ineq = [c_ineq; (value-value_max)/value_ref];
        end
        
        % Apply an upper limit to the current variable if relevant
        if isempty(value_min) ~= 1
            c_ineq = [c_ineq; -(value-value_min)/value_ref];
        end
        
    end
    
end


%% Compute the vector of inequality constraints
% MATLAB uses constraints of the type c_eq=0 

% Load variables
p_out_target = fixed_parameters.thermodynamic_input.p_out;                 % Target for the static outlet pressure
p_out_comp   = turbine_data.plane(end).p;                                  % Computed static outlet pressure
p_out_error  = (p_out_comp-p_out_target)/p_out_target;                     % Relative error in the static outlet pressure
Y_error      = [turbine_data.cascade.Y_error]';                            % Vector of pressure loss coefficient errors

% Compute the vector of inequality constraints
c_eq = [p_out_error; Y_error];

% Store the vectors of equality and inequality constraints
turbine_data.optimization.c = c_ineq;
turbine_data.optimization.c_eq = c_eq;


%% Prepare a structure to summarize the constraints
% Initialize cells
header    = {};
name      = {};
value     = {};
value_min = {};
value_max = {};
applied   = {};
satisfied = {};
active    = {};

% Define the header names
header(1) = {'Variable name'};
header(2) = {'Lower bound'};
header(3) = {'Value'};
header(4) = {'Upper bound'};
header(5) = {'Applied'};
header(6) = {'Satisfied'};
header(7) = {'Active'};

% Add the equality constraints manually
% Outlet pressure error
kk = 1;
name(kk)      = {'Exit pressure error'};
value(kk)     = {p_out_error};
value_min(kk) = {0.00};
value_max(kk) = {0.00};
applied(kk)   = {'yes'};
active(kk)    = {'yes'};
if abs(p_out_error) <= 1e-3
    satisfied(kk) = {'yes'};
else
    satisfied(kk) = {'no'};
end
kk = kk + 1;

% Loss coefficient error
for i = 1:length(Y_error)
    name(kk)      = {'Loss coeff. error'};
    value(kk)     = {Y_error(i)};
    value_min(kk) = {0.00};
    value_max(kk) = {0.00};
    applied(kk)   = {'yes'};
    active(kk)    = {'yes'};
    if abs(Y_error(i)) <= 1e-3
        satisfied(kk) = {'yes'};
    else
        satisfied(kk) = {'no'};
    end
    kk = kk + 1;
end

% Loop through all the variables with inqeuality constraints
for k = 1:length(constraint_names)
    
    % Get the value of the current variable
    value_ = values_ineq.(constraint_names{k});
    
    % Get the limits of the current variable
    value_min_ = constraints.(constraint_names{k}).min;
    value_max_ = constraints.(constraint_names{k}).max;
    
    % Get the reference value of the current variable
    value_ref = constraints.(constraint_names{k}).ref;
    
    % Number of elements of the array containing the current variable
    N_var = length(value_);
    
    % Convert constraint names to full names
    constraint_name = constraint_names{k};
    
    if strcmp(constraint_name, 'flare_angle')
        constraint_name = 'Flaring angle';
    end 
    
    if strcmp(constraint_name, 'r_ht')
        constraint_name = 'Hub-to-tip ratio';
    end 
    
    if strcmp(constraint_name, 'reaction')
        constraint_name = 'Degree of reaction';
    end
    
    if strcmp(constraint_name, 'Ma_rel')
        constraint_name = 'Relative Mach number';
    end
    
    if strcmp(constraint_name, 'Ma_diffuser')
        constraint_name = 'Outlet Mach number';
    end
    
    if strcmp(constraint_name, 'beta_in_stator')
        constraint_name = 'Stator inlet angle';
    end
    
    if strcmp(constraint_name, 'beta_in_rotor')
        constraint_name = 'Rotor inlet angle';
    end
    
    if strcmp(constraint_name, 'height')
        constraint_name = 'Blade height';
    end
    
    if strcmp(constraint_name, 'chord')
        constraint_name = 'Blade chord';
    end
    
    if strcmp(constraint_name, 'RPM')
        constraint_name = 'Angular speed (RPM)';
    end
    
    % Save the information of the current variable
    for i = 1:N_var
        
        % Save constraint information
        name(kk)      = {constraint_name};
        value(kk)     = {value_(i)};
        value_min(kk) = {value_min_};
        value_max(kk) = {value_max_};
        applied(kk)   = {constraints.(constraint_names{k}).applied};
        
        % No lower limit (the limit is minus infinity)
        if isempty(value_min_) == 1
            value_min__ = -inf;
        else
            value_min__ = value_min_;
        end
        
        % No upper limit (the limit is plus infinity)
        if isempty(value_max_) == 1
            value_max__ = inf;
        else 
            value_max__ = value_max_;
        end
        
        % Check if the constraints are satisfied (arbitrary tolerance)
        if strcmp(constraints.(constraint_names{k}).applied,'yes') == 1
            if (value_(i)-value_min__)/value_ref >= -1e-3 && (value_max__-value_(i))/value_ref >= -1e-3
                satisfied(kk) = {'yes'};
            else
                satisfied(kk) = {'no'};
            end
        else
            satisfied(kk) = {[]};
        end
        
        % Check if the constraints are active (arbitrary tolerance)
        if strcmp(constraints.(constraint_names{k}).applied,'yes') == 1
            if abs((value_(i)-value_min__)/value_ref) <= 1e-3 || abs((value_max__-value_(i))/value_ref) <= 1e-3
                active(kk) = {'yes'};
            else
                active(kk) = {'no'};
            end
        else
            active(kk) = {[]};
        end
        
        % Update looping veriable
        kk = kk + 1;
        
    end

end


% Store cells as a structure
constraint_summary.header    = header;
constraint_summary.name      = name;
constraint_summary.value     = value;
constraint_summary.value_min = value_min;
constraint_summary.value_max = value_max;
constraint_summary.applied   = applied;
constraint_summary.active    = active;
constraint_summary.satisfied = satisfied;
turbine_data.optimization.constraint_summary = constraint_summary;


end