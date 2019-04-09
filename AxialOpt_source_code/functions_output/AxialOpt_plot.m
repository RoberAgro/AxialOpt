function [] = AxialOpt_plot(turbine_data,choose_plots,my_filename,my_path)
%% Plot the graphs specified in the "choose_plots" structure
% Author: Roberto Agromayor

% Run the figure settings function
AxialOpt_plot_settings()

% Choose whether to save the plots or not
save = choose_plots.save;

% Plot the diagrams specified in the choose_plots variable
if strcmp(choose_plots.diagram_hs,'yes')
    AxialOpt_plot_diagram_hs(turbine_data,my_filename,my_path,save)
end

if strcmp(choose_plots.diagram_Ts,'yes')
    AxialOpt_plot_diagram_Ts(turbine_data,my_filename,my_path,save)
end

if strcmp(choose_plots.axial_tangential,'yes')
    AxialOpt_plot_view_axial_tangential(turbine_data,my_filename,my_path,save)
end

if strcmp(choose_plots.axial_radial,'yes')
    AxialOpt_plot_view_axial_radial(turbine_data,my_filename,my_path,save)
end

if strcmp(choose_plots.axial_radial_diffuser,'yes')
    AxialOpt_plot_view_axial_radial_diffuser(turbine_data,my_filename,my_path,save)
end

if strcmp(choose_plots.triangles_A,'yes')
    AxialOpt_plot_velocity_triangles_A(turbine_data,my_filename,my_path,save)
end

if strcmp(choose_plots.triangles_B,'yes')
    AxialOpt_plot_velocity_triangles_B(turbine_data,my_filename,my_path,save)
end

if strcmp(choose_plots.loss_breakdown,'yes')
    AxialOpt_plot_loss_breakdown(turbine_data,my_filename,my_path,save)
end


end
