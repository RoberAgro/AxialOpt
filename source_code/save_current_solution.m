function stop = save_current_solution(x,~,~,fixed_parameters)

% Use a persistent variable to keep track of the number of iterations
persistent iter
if isempty(iter)
    iter = 0;
end
iter = iter+1;

% Evaluate the cycle model for the current vector of independent variables
turbine_data = evaluate_turbine_model(x,fixed_parameters);

% Save the current solution as a MATLAB data structure
save(fullfile(fixed_parameters.results_path,[fixed_parameters.project_name, '_iter', num2str(iter), '.mat']),'turbine_data')

% Return a false stop flag
stop = false;

end