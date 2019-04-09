function [] = AxialOpt_plot_diagram_Ts(turbine_data,my_filename,my_path,save)
%% Plot the h-s diagram for the expansion along the turbine
% Author: Roberto Agromayor

% The h-s diagram was chosen to represent the expansion because it is the
% most relevant for turbomachinery (enthalpy is related to work exchange
% and entropy is related to irrversibility


%% Initialize the problem
% User input parameters
p_turbine = [turbine_data.plane.p];
T0_turbine = [turbine_data.plane.T0];
T_turbine = [turbine_data.plane.T];
h_turbine = [turbine_data.plane.h];
h0_turbine = [turbine_data.plane.h0];
s_turbine = [turbine_data.plane.s];
p_out = turbine_data.overall.p_out;

% Define fluid properties
fluid = turbine_data.overall.fluid;                                        % Set the composition of the working fluid
T_crit = turbine_data.fluid_properties.T_crit;
p_crit = turbine_data.fluid_properties.p_crit;
p_min = turbine_data.fluid_properties.p_min;
p_max = turbine_data.fluid_properties.p_max;
T_min = 0.9*min(T_turbine);
T_max = max(1.1*T_crit, 1.1*max(T_turbine));


% Compute saturation line
try
    % Preliminary definitions
    N_sat = 50;                                                            % Number of points
    T_trip = refpropm('T','R',0,' ',0,fluid);                              % Triple temperature
    p_trip = refpropm('p','R',0,' ',0,fluid);                              % Triple temperature
    T_min = T_trip;
    p_min = p_trip;
    
    % Computation of the saturation line
    % Do not reach the critical point exactly (REFPROP may fail)
    try
        [T_sat, s_sat] = sat_line(fluid,T_trip,T_crit-0.1,'T','s',N_sat);
    catch
        [T_sat, s_sat] = sat_line(fluid,T_trip,T_crit-1.0,'T','s',N_sat);
    end
    saturation_line = 'yes';
    
    % Compute the location of the critical point
    s_crit = refpropm('s','T',0.999*T_crit,'p',p_crit,fluid);
    
catch
    % Avoid saturation line computation in case it fails
    s_sat = [];
    T_sat = [];
    s_crit = [];
    T_crit = [];
    saturation_line = 'no';
    disp('Skipping saturation line computation in the h-s diagram')
    
end


%% Enthalpy-entropy diagram
fig = figure; ax_fig = gca;
hold on; box on; axis square;
xlabel({' ';'$s$ -- Entropy (kJ/kg$\cdot$K)'});
ylabel({'$T$ -- Temperature ($^{\circ}$C)';' '});
ax_fig.YAxis.TickLabelFormat = '%.0f';
ax_fig.XAxis.TickLabelFormat = '%.2f';

% % Axis limits
% T_minplot = 0.9*min(T_turbine);
% T_maxplot = 1.10*max(T0_turbine);
% s_minplot = 0.9*min(s_turbine);
% s_maxplot = 1.10*max(s_turbine);
% ax_limits = [s_minplot/1000 s_maxplot/1000 T_minplot-273.15 T_maxplot-273.15];
% axis(ax_limits)


%% Plot the diagram
% Number of isobar points
N_isobar = 100;

% Internal stage isobars
for i = 1:length(p_turbine)
[~, s_isobar, T_isobar] = hs_isobar(p_turbine(i),p_crit,T_crit,p_min,T_min,T_max,fluid,N_isobar,saturation_line);
plot(s_isobar/1000, T_isobar-273.15, 'k')
end

% Exhaust pressure isobar
[~, s_isobar, T_isobar] = hs_isobar(p_out,p_crit,T_crit,p_min,T_min,T_max,fluid,N_isobar,saturation_line);
plot(s_isobar/1000, T_isobar-273.15, 'b')

% Plot the saturation line
plot(s_sat/1000,T_sat-273.15,'k','LineWidth',0.75)

% Plot the critical point
plot(s_crit/1000,T_crit-273.15,'ko','MarkerFaceColor','w','MarkerSize',2.5)

% Plot the static and stagnation isobars
plot(s_turbine/1000,T_turbine-273.15,'bo-','LineWidth',0.5,'MarkerSize',2.5,'MarkerFaceColor','w')
plot(s_turbine/1000,T0_turbine-273.15,'ro-','LineWidth',0.5,'MarkerSize',2.5,'MarkerFaceColor','w')

% Save the figure
name = '_diagram_Ts';
if strcmp(save,'lite') == 1
    saveas(fig,fullfile(my_path,[my_filename,name,'.pdf']),'pdf')
elseif strcmp(save,'full') == 1
%     saveas(fig,fullfile(my_path,[my_filename,name]),'fig')
    export_fig(fig,fullfile(my_path,[my_filename,name]),'-png','-r1000')
    export_fig(fig,fullfile(my_path,[my_filename,name]),'-pdf','-painters')
%     export_fig(fig,fullfile(my_path,[my_filename,name]),'-eps','-painters')
elseif strcmp(save,'no') ~= 1
    error('Choose a valid saving option')
end



end