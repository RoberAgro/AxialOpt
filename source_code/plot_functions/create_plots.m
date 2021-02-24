function turbine_data = create_plots(x,fixed_parameters,filename_suffix)

% Load setting for beautiful plots
set_plot_options()

% Evaluate the optimization problem
if nargin == 2
    filename = fixed_parameters.project_name;
elseif nargin == 3
    filename = [fixed_parameters.project_name, '_', filename_suffix];
else
    error('The number of arguments must be 2 or 3')
end
filepath = fixed_parameters.results_path;
turbine_data = evaluate_optimization_problem(x,fixed_parameters);

% Choose whether to save the plots or not
% 'choose_plots' is an structure that contains what diagrams to draw
save = fixed_parameters.choose_plots.save;

% Choose whether to plot the saturation line or not
plot_satline = fixed_parameters.choose_plots.plot_satline;

% Plot the diagrams specified in the choose_plots variable
if strcmp(fixed_parameters.choose_plots.diagram_hs,'yes')
    plot_diagram_hs(turbine_data,plot_satline,filename,filepath,save)
end

if strcmp(fixed_parameters.choose_plots.diagram_Ts,'yes')
    plot_diagram_Ts(turbine_data,plot_satline,filename,filepath,save)
end

if strcmp(fixed_parameters.choose_plots.axial_tangential,'yes')
    plot_view_axial_tangential(turbine_data,filename,filepath,save)
end

if strcmp(fixed_parameters.choose_plots.axial_radial,'yes')
    diffuser_model = fixed_parameters.design_input.diffuser_model;
    if strcmp(diffuser_model,'no')
        plot_view_axial_radial(turbine_data,filename,filepath,save)
    elseif strcmp(diffuser_model,'1D') || strcmp(diffuser_model,'isentropic')
        plot_view_axial_radial_diffuser(turbine_data,filename,filepath,save)
    else
        error("The option diffuser_model must be 'no', '1D' or 'isentropic'")
    end
end

if strcmp(fixed_parameters.choose_plots.triangles_A,'yes')
    plot_velocity_triangles_A(turbine_data,filename,filepath,save)
end

if strcmp(fixed_parameters.choose_plots.triangles_B,'yes')
    plot_velocity_triangles_B(turbine_data,filename,filepath,save)
end

if strcmp(fixed_parameters.choose_plots.loss_breakdown,'yes')
    plot_loss_breakdown(turbine_data,filename,filepath,save)
end


end
