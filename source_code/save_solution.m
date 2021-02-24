function [] = save_solution(turbine_data,my_filename,my_path)

% Print the results in .txt files

% Print global variables
file_name = fopen(fullfile(my_path,[my_filename,'_solution_overall.txt']), 'w');
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

% Print cascade variables
n_cascades = turbine_data.overall.n_cascades;
file_name = fopen(fullfile(my_path,[my_filename,'_solution_cascades.txt']), 'w');
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
    fprintf(file_name,my_format,field_name{k},turbine_data.cascade(:).(field_name{k}));
end
fclose(file_name);

% Print plane variables
file_name = fopen(fullfile(my_path,[my_filename,'_solution_stations.txt']), 'w');
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

% Print optimization problem output
file_name = fopen(fullfile(my_path,[my_filename,'_solution_optimization.txt']), 'w');
bounds_summary = turbine_data.optimization.bounds_summary;
header    = bounds_summary.header;
name      = bounds_summary.name;
value     = bounds_summary.value;
value_min = bounds_summary.value_min;
value_max = bounds_summary.value_max;
satisfied = bounds_summary.satisfied;
active    = bounds_summary.active;
fprintf(file_name,'|-----------------------------------------------------------------------------------------------|\n');
fprintf(file_name,'|-------------------------- Values and bounds of the design variables --------------------------|\n');
fprintf(file_name,'|-----------------------------------------------------------------------------------------------|\n');
fprintf(file_name,'%25s %15s %12s %15s %12s %10s \n', header{1}, header{2}, header{3}, header{4}, header{5}, header{6});
for i = 1:size(name,2)
    fprintf(file_name,'%25s %15.4f %12.4f %15.4f %12s %10s \n', name{i}, value_min{i}, value{i}, value_max{i}, satisfied{i}, active{i});
end
fprintf(file_name,'|-----------------------------------------------------------------------------------------------|\n');

constraint_summary = turbine_data.optimization.constraint_summary;
header    = constraint_summary.header;
name      = constraint_summary.name;
value     = constraint_summary.value;
value_min = constraint_summary.value_min;
value_max = constraint_summary.value_max;
applied   = constraint_summary.applied;
satisfied = constraint_summary.satisfied;
active    = constraint_summary.active;
fprintf(file_name,'\n');
fprintf(file_name,'|-------------------------------------------------------------------------------------------------------------|\n');
fprintf(file_name,'|-------------------------------------------- Problem constraints --------------------------------------------|\n');
fprintf(file_name,'|-------------------------------------------------------------------------------------------------------------|\n');
fprintf(file_name,'%25s %15s %12s %15s %10s %12s %10s \n', header{1}, header{2}, header{3}, header{4}, header{5}, header{6}, header{7});
for i = 1:size(name,2)
    if strcmp(name{i}, 'Angular speed (RPM)')
        fprintf(file_name,'%25s %15.1f %12.1f %15.1f %10s %12s %10s \n', name{i}, value_min{i}, value{i}, value_max{i}, applied{i}, satisfied{i}, active{i});
    else
        fprintf(file_name,'%25s %15.4f %12.4f %15.4f %10s %12s %10s \n', name{i}, value_min{i}, value{i}, value_max{i}, applied{i}, satisfied{i}, active{i});
    end
end
fprintf(file_name,'|-------------------------------------------------------------------------------------------------------------|\n');
fclose(file_name);


end
