function [] = AxialOpt_print_results(turbine_data,my_filename,my_path)
%% Print the results of the optimization in .txt files
% Author: Roberto Agromayor

% Print global variables
file_name = fopen(fullfile(my_path,[my_filename,'_overall_variables.txt']), 'w');
field_name = fieldnames(turbine_data.overall);
for k = 1:length(field_name)-1
    my_var = turbine_data.overall.(field_name{k}); % Use to check data type
    if ischar(my_var) == 1
        fprintf(file_name,'%s\t%s\n',field_name{k},my_var);
    elseif isinteger(my_var) == 1
        fprintf(file_name,'%s\t%i\n',field_name{k},my_var);
    else
        fprintf(file_name,'%s\t%.6e\n',field_name{k},turbine_data.overall.(field_name{k}));
    end
end
fclose(file_name);

% Print diffuser variables
file_name = fopen(fullfile(my_path,[my_filename,'_diffuser_variables.txt']), 'w');
field_name = fieldnames(turbine_data.diffuser);
for k = 1:length(field_name)-1
    my_var = turbine_data.diffuser.(field_name{k}); % Use to check data type
    if ischar(my_var) == 1
        fprintf(file_name,'%s\t%s\n',field_name{k},my_var);
    elseif isinteger(my_var) == 1
        fprintf(file_name,'%s\t%i\n',field_name{k},my_var);
    else
        fprintf(file_name,'%s\t%.6e\n',field_name{k},turbine_data.diffuser.(field_name{k}));
    end
end
fclose(file_name);


% Print cascade variables
n_cascades = turbine_data.overall.n_cascades;
file_name = fopen(fullfile(my_path,[my_filename,'_cascade_variables.txt']), 'w');
field_name = fieldnames(turbine_data.cascade);
for k = 1:length(field_name)-1
    my_var = turbine_data.cascade(:).(field_name{k});  % Use to check data type
    my_format = '%s';
    for i = 1:n_cascades
        if ischar(my_var) == 1
            my_format = strcat(my_format,'\t%s');
        elseif isinteger(my_var) == 1
            my_format = strcat(my_format,'\t%i');
        else
            my_format = strcat(my_format,'\t%.6e');
        end
    end
    my_format = strcat(my_format,'\n');
    fprintf(file_name, my_format,field_name{k},turbine_data.cascade(:).(field_name{k}));
end
fclose(file_name);


% Print plane variables
file_name = fopen(fullfile(my_path,[my_filename,'_plane_variables.txt']), 'w');
field_name = fieldnames(turbine_data.plane);
for k = 1:length(field_name)-1
    my_var = turbine_data.plane(:).(field_name{k});  % Use to check data type
    my_format = '%s';
    for i = 1:2*n_cascades+1
        if ischar(my_var) == 1
            my_format = strcat(my_format,'\t%s');
        elseif isinteger(my_var) == 1
            my_format = strcat(my_format,'\t%i');
        else
            my_format = strcat(my_format,'\t%.6e');
        end
    end
    my_format = strcat(my_format,'\n');        
    fprintf(file_name,my_format,field_name{k},turbine_data.plane(:).(field_name{k}));
end
fclose(file_name);


% Print optimization results
exitflag = turbine_data.optimization.exitflag;
output = turbine_data.optimization.output;
file_name = fopen(fullfile(my_path,[my_filename,'_optimization.txt']), 'w');
fprintf(file_name, 'Independent variables\n');
fprintf(file_name, ['variable','\t','lower_bound','\t','value','\t','upper_bound','\t','satisfied','\t','active','\n']);
for k = 1:length(turbine_data.optimization.check_bounds)
    fprintf(file_name, ['%s','\t','%.6e','\t','%.6e','\t','%.6e','\t','%s','\t','%s','\n'], ...
            turbine_data.optimization.check_bounds(k).variable, ...
            turbine_data.optimization.check_bounds(k).lower_bound, ...
            turbine_data.optimization.check_bounds(k).value, ...
            turbine_data.optimization.check_bounds(k).upper_bound, ...
            turbine_data.optimization.check_bounds(k).satisfied, ...
            turbine_data.optimization.check_bounds(k).active);    
end
fprintf(file_name,'\n\n');
fprintf(file_name, 'Nonlinear constraints\n');
fprintf(file_name, ['variable','\t','min_value','\t','value','\t','max_value','\t','applied','\t','satisfied','\t','active','\n']);
for k = 1:length(turbine_data.optimization.constraint_summary)
    fprintf(file_name, ['%s','\t','%.6e','\t','%.6e','\t','%.6e','\t','%s','\t','%s','\t','%s','\n'], ...
            turbine_data.optimization.constraint_summary(k).variable,    ...
            turbine_data.optimization.constraint_summary(k).min_value,   ...
            turbine_data.optimization.constraint_summary(k).value,       ...
            turbine_data.optimization.constraint_summary(k).max_value,   ...
            turbine_data.optimization.constraint_summary(k).applied,     ...
            turbine_data.optimization.constraint_summary(k).satisfied,   ...
            turbine_data.optimization.constraint_summary(k).active);     
end
fprintf(file_name,'\n\n');
fprintf(file_name,'Optimization algorithm output\n');
fprintf(file_name,['exitflag','\t','iterations','\t','function_count','\t','con_violation','\t','stepsize','\n']);
fprintf(file_name,'%i\t%i\t%i\t%.6e\t%.6e', ...
        exitflag, ...
        output.iterations, ...
        output.funcCount, ...
        output.constrviolation, ...
        output.stepsize);
fprintf(file_name,'\n\n');
fclose(file_name);


end
