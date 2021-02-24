function [] = plot_diagram_Ts(turbine_data,plot_satline,my_filename,my_path,save)

% Plot the T-s diagram for the expansion along the turbine

% Load variables
p_turbine = [turbine_data.plane.p];
T0_turbine = [turbine_data.plane.T0];
T_turbine = [turbine_data.plane.T];
h_turbine = [turbine_data.plane.h];
h0_turbine = [turbine_data.plane.h0];
s_turbine = [turbine_data.plane.s];
p_out = turbine_data.overall.p_out;
fluid = turbine_data.overall.fluid;
T_crit = turbine_data.fluid_properties.T_crit;
p_crit = turbine_data.fluid_properties.p_crit;
T_trip = turbine_data.fluid_properties.T_trip;

% Create the figure
fig = figure(); ax_fig = gca; 
hold on; box on;
pbaspect([1.15 1 1])
xlabel({' ';'$s$ -- Entropy (kJ/kg$\,$K)'});
ylabel({'$T$ -- Temperature ($^{\circ}$C)';' '});
ax_fig.YAxis.TickLabelFormat = '%.0f';
ax_fig.XAxis.TickLabelFormat = '%.2f';

% Define axis limits
delta_T = max(T_turbine) - min(T_turbine);
delta_s = max(s_turbine) - min(s_turbine);
T_minplot = min(T_turbine) - delta_T/4;
T_maxplot = max(T_turbine) + delta_T/4;
s_minplot = min(s_turbine) - delta_s/1;
s_maxplot = max(s_turbine) + delta_s/1;
ax_limits = [s_minplot/1000 s_maxplot/1000 T_minplot-273.15 T_maxplot-273.15];
axis(ax_limits)

% Plot the saturation line
if strcmp(plot_satline,'yes') == 1
    [T_sat, s_sat] = compute_sat_line(fluid,T_trip,T_crit,'T','S',200);
    plot(s_sat/1000,T_sat-273.15,'k','LineWidth',0.75)
    T_minplot = min(T_minplot, min(T_sat));
    T_maxplot = 1.1*max(T_maxplot, max(T_sat));
    s_minplot = min(s_minplot, min(s_sat));
%     s_maxplot = max(s_maxplot, max(s_sat));
    s_maxplot = s_maxplot + (max(s_sat)-min(s_sat))/2;
    ax_limits = [s_minplot/1000 s_maxplot/1000 T_minplot-273.15 T_maxplot-273.15];
    axis(ax_limits)
end

% Plot the isobars corresponding to the pressures within the turbine
for i = 1:length(p_turbine)
    [~, s_isobar, T_isobar] = compute_isobar(p_turbine(i),fluid,200);
    plot(s_isobar/1000, T_isobar-273.15, 'k')
end

% Plot the isobars corresponding to the exhaust pressure
[~, s_isobar, T_isobar] = compute_isobar(p_out,fluid,200);
plot(s_isobar/1000, T_isobar-273.15, 'b')

% Plot the critical point
s_crit = prop_calculation('S','T',T_crit,'P',p_crit,fluid);
plot(s_crit/1000,T_crit-273.15,'ko','MarkerFaceColor','w','MarkerSize',2.5)

% Plot the static and stagnation states
plot(s_turbine/1000,T_turbine-273.15,'bo-','LineWidth',0.5,'MarkerSize',2.5,'MarkerFaceColor','w')
plot(s_turbine/1000,T0_turbine-273.15,'ro-','LineWidth',0.5,'MarkerSize',2.5,'MarkerFaceColor','w')

% Save the figure
name = 'Ts_diagram';
if save == 1
    saveas(fig,fullfile(my_path,[my_filename,'_',name,'.pdf']),'pdf')
elseif save == 2
%     saveas(fig,fullfile(my_path,[my_filename,'_',name]),'fig')
%     export_fig(fig,fullfile(my_path,[my_filename,'_',name]),'-png','-r1000')
%     export_fig(fig,fullfile(my_path,[my_filename,'_',name]),'-eps','-painters')
    export_fig(fig,fullfile(my_path,[my_filename,'_',name]),'-pdf','-painters')
elseif save ~= 0
    error('Choose a valid saving option')
end

end

