%% AxialOpt - Optimization of axial turbines
% Author: Roberto Agromayor
% Date: 05/04/2019
% Description: This program is capable of optimizing axial turbines of any
% number of stages using the thermodynamic boundary conditions obtained
% from a system analysis (for instance, a power cycle optimization)


%% Initialize the program
% Clear all variables and close all figures
clearvars
close all
clc

% Add the path to the AxialOpt_source directory
% The function genpath() is very convenient
addpath(genpath('../../AxialOpt_source'))

% Call the script to define the design parameters and optimization problem
input_parameters

% Create a directory to store the results (you can personalize the name)
my_filename = ['results_example_', fluid];
my_results_path = fullfile(pwd, my_filename);
if exist(my_results_path, 'dir') ~= 7    % This is a MATLAB's convention
    mkdir(my_results_path)
end

% Create a directory to store the figures (you can personalize the name)
my_filename = ['figures_example_', fluid];
my_figures_path = fullfile(pwd, my_filename);
if exist(my_figures_path, 'dir') ~= 7    % This is a MATLAB's convention
    mkdir(my_figures_path)
end

% Be tidy and clear the worskpace now that we stored everything (optional)
clearvars -except optimization_problem parameters_turbine my_filename my_figures_path my_results_path    


                      
%% Choose what figures should be plotted
choose_plots = struct('diagram_hs',            'yes', ...                  % Enthalpy-entropy diagram of the expansion
                      'diagram_Ts',            'yes', ...                  % Temperature-entropy diagram of the expansion
                      'axial_tangential',      'yes', ...                  % View of the blade cascades
                      'axial_radial',          'yes', ...                  % View of the meridional plane
                      'axial_radial_diffuser', 'yes', ...                  % View of the meridional plane including the diffuser
                      'triangles_A',           'yes', ...                  % Velocity triangles with blade velocity common origin
                      'triangles_B',           'yes', ...                  % Velocity triangles with flow velocity common origin
                      'loss_breakdown',        'yes', ...                  % Breakdown of the efficiency loss
                      'save',                  'no');                    % Choose wether to save the figure or not. Options:
                                                                           % 1) 'full' saves the figure in vector and in high quality raster format using the export_fig library (long computational time)
                                                                           % 2) 'lite' saves the figure in vector format using the default MATLAB export function (short computational time)
                                                                           % 3) 'no' do not save the figures
                                                                           
% Using saving option full requires: 1) ghoscript and 2) pdftops
%  http://www.ghostscript.com
%  http://xpdfreader.com

                                                                           
%% Load a previous solution as initial guess
% % Comment these line to optimize a new turbine design from scratch
% load(fullfile(my_results_path,'turbine_data.mat'))
% optimization_problem.x0 = turbine_data.optimization.x;


%% Plot the initial guess
% close all
% turbine_data = AxialOpt_computation(optimization_problem.x0,parameters_turbine);
% AxialOpt_plot(turbine_data,choose_plots,my_filename,my_figures_path);


%% Solve the optimization problem
% Run a single optimization or a multistart optimization
%(this feature of the code is not mature)
multistart = 'no';
N_multistart = 15;
tic
[x_opt,f_opt,exitflag,output,solutions] = AxialOpt_optimization(parameters_turbine,optimization_problem,multistart,N_multistart);
toc


%% Plot the optimum solution
close all
turbine_data = AxialOpt_computation(x_opt,parameters_turbine);
AxialOpt_plot(turbine_data,choose_plots,my_filename,my_figures_path);


%% Print the optimization results in .txt files
turbine_data = AxialOpt_computation(x_opt,parameters_turbine);
turbine_data.optimization.exitflag = exitflag;
turbine_data.optimization.output = output;
AxialOpt_print_results(turbine_data,my_filename,my_results_path)


%% Save the turbine_data structure
save(fullfile(my_results_path,[my_filename,'.mat']),'turbine_data')
save(fullfile(my_results_path,['turbine_data','.mat']),'turbine_data')


        