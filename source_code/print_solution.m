function [] = print_solution(turbine_data)

% Print bounds
bounds_summary = turbine_data.optimization.bounds_summary;
header    = bounds_summary.header;
name      = bounds_summary.name;
value     = bounds_summary.value;
value_min = bounds_summary.value_min;
value_max = bounds_summary.value_max;
satisfied = bounds_summary.satisfied;
active    = bounds_summary.active;
fprintf('\n')
fprintf('|-----------------------------------------------------------------------------------------------|\n')
fprintf('|-------------------------- Values and bounds of the design variables --------------------------|\n')
fprintf('|-----------------------------------------------------------------------------------------------|\n')
fprintf('%26s %15s %12s %15s %12s %10s \n', header{1}, header{2}, header{3}, header{4}, header{5}, header{6});
for i = 1:size(name,2)
    fprintf('%26s %15.4f %12.4f %15.4f %12s %10s \n', name{i}, value_min{i}, value{i}, value_max{i}, satisfied{i}, active{i});
end
fprintf('|-----------------------------------------------------------------------------------------------|\n')

% Print constraints
constraint_summary = turbine_data.optimization.constraint_summary;
header    = constraint_summary.header;
name      = constraint_summary.name;
value     = constraint_summary.value;
value_min = constraint_summary.value_min;
value_max = constraint_summary.value_max;
applied   = constraint_summary.applied;
satisfied = constraint_summary.satisfied;
active    = constraint_summary.active;
fprintf('\n')
fprintf('|-------------------------------------------------------------------------------------------------------------|\n')
fprintf('|-------------------------------- Values and limits of the problem constraints -------------------------------|\n')
fprintf('|-------------------------------------------------------------------------------------------------------------|\n')
fprintf('%26s %15s %12s %15s %10s %12s %10s \n', header{1}, header{2}, header{3}, header{4}, header{5}, header{6}, header{7});
for i = 1:size(name,2)
    if strcmp(name{i}, 'Angular speed (RPM)')
        fprintf('%26s %15.1f %12.1f %15.1f %10s %12s %10s \n', name{i}, value_min{i}, value{i}, value_max{i}, applied{i}, satisfied{i}, active{i});
    else
        fprintf('%26s %15.4f %12.4f %15.4f %10s %12s %10s \n', name{i}, value_min{i}, value{i}, value_max{i}, applied{i}, satisfied{i}, active{i});
    end
end
fprintf('|-------------------------------------------------------------------------------------------------------------|\n')

end

