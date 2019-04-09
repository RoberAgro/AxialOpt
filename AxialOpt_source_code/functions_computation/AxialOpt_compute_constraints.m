function turbine_data = AxialOpt_compute_constraints(turbine_data,parameters)
%% Compute the vector of inequality constraints
% Author: Roberto Agromayor

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
constraints = parameters.constraints;

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
p_out_target = parameters.thermodynamic_input.p_out;                       % Target for the static outlet pressure
p_out_comp   = turbine_data.plane(end).p;                                  % Computed static outlet pressure
p_out_error  = (p_out_comp-p_out_target)/p_out_target;                     % Relative error in the static outlet pressure
Y_error      = [turbine_data.cascade.Y_error]';                            % Vector of pressure loss coefficient errors

% Compute the vector of inequality constraints
c_eq = [p_out_error; Y_error];


%% Store the vectors of equality and inequality constraints
% Store the constraints
turbine_data.optimization.c = c_ineq;
turbine_data.optimization.c_eq = c_eq;



%% Prepare an array of structures to summarize the constraints
% The information of constraint_summary is later printed on a .txt file
% Initialize the array of structures
constraint_summary = struct('variable',        [], ...
                            'min_value',       [], ...
                            'value',           [], ...
                            'max_value',       [], ...
                            'applied',         [], ...
                            'satisfied',       [], ...
                            'active',          []);
                        
% Add the equality constraints manually
% Outlet pressure error
kk = 1;
constraint_summary(kk).variable  = 'p_out_error';
constraint_summary(kk).min_value = 0.00;
constraint_summary(kk).value     = p_out_error;
constraint_summary(kk).max_value = 0.00;
constraint_summary(kk).applied   = 'yes';
constraint_summary(kk).active    = 'yes';
if abs(p_out_error) <= 1e-3
    constraint_summary(kk).satisfied = 'yes';
else
    constraint_summary(kk).satisfied = 'no';
end
kk = kk + 1;

% Loss coefficient error
for i = 1:length(Y_error)
    constraint_summary(kk).variable  = 'Y_error';
    constraint_summary(kk).min_value = 0.00;
    constraint_summary(kk).value     = Y_error(i);
    constraint_summary(kk).max_value = 0.00;
    constraint_summary(kk).applied   = 'yes';
    constraint_summary(kk).active    = 'yes';
    if abs(Y_error(i)) <= 1e-3
        constraint_summary(kk).satisfied = 'yes';
    else
        constraint_summary(kk).satisfied = 'no';
    end
    kk = kk + 1;
end

% Loop through all the variables with inqeuality constraints
for k = 1:length(constraint_names)
    
    % Get the value of the current variable
    value = values_ineq.(constraint_names{k});
    
    % Get the limits of the current variable
    value_min = constraints.(constraint_names{k}).min;
    value_max = constraints.(constraint_names{k}).max;
    
    % Get the reference value of the current variable
    value_ref = constraints.(constraint_names{k}).ref;
    
    % Number of elements of the array containing the current variable
    N_var = length(value);
    
    % Save the information of the current variable
    for i = 1:N_var
        
        % Save constraint information
        constraint_summary(kk).variable  = constraint_names{k};
        constraint_summary(kk).min_value = value_min;
        constraint_summary(kk).value     = value(i);
        constraint_summary(kk).max_value = value_max;
        constraint_summary(kk).applied   = constraints.(constraint_names{k}).applied;
        
        % Use uxiliary variables in case the min/max limit does not exist
        % These are required for the logical statements used to determine
        % if the constraint is satisfied or active
        
        % No lower limit (the limit is minus infinity)
        if isempty(value_min) == 1
            min_value = -inf;
        else
            min_value = value_min;
        end
        
        % No upper limit (the limit is plus infinity)
        if isempty(value_max) == 1
            max_value = inf;
        else
            max_value = value_max;
        end
        
        % Check if the constraints are satisfied (arbitrary tolerance)
        if strcmp(constraints.(constraint_names{k}).applied,'yes') == 1
            if (value(i)-min_value)/value_ref >= -1e-3 && (max_value-value(i))/value_ref >= -1e-3
                constraint_summary(kk).satisfied = 'yes';
            else
                constraint_summary(kk).satisfied = 'no';
            end
        end
        
        % Check if the constraints are active (arbitrary tolerance)
        if strcmp(constraints.(constraint_names{k}).applied,'yes') == 1
            if abs((value(i)-min_value)/value_ref) <= 1e-3 || abs((max_value-value(i))/value_ref) <= 1e-3
                constraint_summary(kk).active = 'yes';
            else
                constraint_summary(kk).active = 'no';
            end
        end
        
        % Update looping veriable
        kk = kk + 1;
        
    end


end

% constraint_summary is used to print the optimization results in .txt files
turbine_data.optimization.constraint_summary = constraint_summary;


end