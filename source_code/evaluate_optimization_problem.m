function [turbine_data, f, c, c_eq] = evaluate_optimization_problem(x,parameters)

% Evaluate the turbine mean-line model
turbine_data = evaluate_turbine_model(x,parameters);                         

% Extract the objective function value
if strcmp(parameters.design_input.obj_function,'total-to-static') == 1
    eta_ts = turbine_data.overall.eta_ts;
    f = -eta_ts;
elseif strcmp(parameters.design_input.obj_function,'total-to-total') == 1
    eta_tt = turbine_data.overall.eta_tt;
    f = -eta_tt;
else
    error('The objective function must the the total-to-static or the total-to-total isentropic efficiency')
end

% Extract the vector of inequality constraints
c = turbine_data.optimization.c;

% Extract the vector of  equality constraints
c_eq = turbine_data.optimization.c_eq;

end