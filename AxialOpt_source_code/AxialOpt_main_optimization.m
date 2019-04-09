function [x_opt,f_opt,exitflag,output,solutions]  = AxialOpt_main_optimization(parameters,optimization_problem,multistart,N)
%% Solve the optimization problem
% Author: Roberto Agromayor

% The optimization algorithm (fmincon) is called in this function
% This function uses nested functions to have a temporal storage of the 
% objective function and constraints to avoid the evaluation of turbine
% model each time the objective function or the constraints are required
% This trick halves the number of function evaluations

x_last   = []; % Degrees of freedom in the last computation
f_bis    = []; % Objective function in storage
c_bis    = []; % Nonlinear inequality constraints in storage
c_eq_bis = []; % Nonlinear equality constraints in storage

% Objective function and nonlinear constraints functions (nested below)
optimization_problem.objective = @(x) AxialOpt_objective_function(x,parameters);
optimization_problem.nonlcon   = @(x) AxialOpt_constraints(x,parameters);

% Use fmincon to optimize the turbine
if strcmp(multistart,'no') == 1
    [x_opt,f_opt,exitflag,output] = fmincon(optimization_problem);
    solutions = [];
elseif strcmp(multistart,'yes') == 1
    ms = MultiStart('Display','iter','StartPointsToRun','bounds');
    % Display must be 'off', 'on', 'iter', or 'final'.
    [x_opt,f_opt,exitflag,output,solutions] = run(ms,optimization_problem,N);
else
    error('Select ''yes'' or ''no'' for the variable mustistart')
end


    function f = AxialOpt_objective_function(x,parameters)
        if ~isequal(x,x_last)      % Check if computation is necessary
           [f_bis,c_bis,c_eq_bis] = AxialOpt_compute_objfun_constr(x,parameters);
            x_last = x;
        end
        f = f_bis;        
    end

    function [c, c_eq] = AxialOpt_constraints(x,parameters)
        if ~isequal(x,x_last)      % Check if computation is necessary
            [f_bis,c_bis,c_eq_bis] = AxialOpt_compute_objfun_constr(x,parameters);
            x_last = x;
        end
        c = c_bis;
        c_eq = c_eq_bis;          
    end


end