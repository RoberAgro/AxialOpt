function [f, c, c_eq] = AxialOpt_compute_objfun_constr(x,parameters)
%% Evalutate the turbine model
% Author: Roberto Agromayor


%% Call the turbine model function for the current set of design variables
turbine_data = AxialOpt_model_turbine(x,parameters);                         


%% Return the objective function to the optimization algorithm
if strcmp(parameters.design_input.obj_function,'total-to-static') == 1
    
    % Total-to-static efficiency
    eta_ts = turbine_data.overall.eta_ts;
    f = -eta_ts;
    
elseif strcmp(parameters.design_input.obj_function,'total-to-total') == 1
    
    % Total-to-total efficiency
    eta_tt = turbine_data.overall.eta_tt;
    f = -eta_tt;
    
else
    
    error('The objective function must the the total-to-static or the total-to-total isentropic efficiency')
    
end


%% Return the nonlinear constraints to the optimization algorithm
% Inequality constraints
c = turbine_data.optimization.c;


% Equality constraints
c_eq = turbine_data.optimization.c_eq;


end