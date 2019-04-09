function [] = AxialOpt_plot_loss_breakdown(turbine_data,my_filename,my_path,save)
%% Loss breakdown in each cascade
% Author: Roberto Agromayor

% Plot the total-to-static efficiency drop in each cascade, including what
% fraction corresponds to profile, secondary, tip clearance, trailing edge
% and kinetic energy losses

% Load the required variables
loss_matrix = turbine_data.overall.loss_matrix;
eta_ts = turbine_data.overall.eta_ts;
n_stages = turbine_data.overall.n_stages;

% Prepare the figure
fig = figure; ax_fig = gca;
title({['Turbine total-to-static efficiency: $\eta_{ts} = ',num2str(eta_ts,'%0.3f'),'$']; ' '})
hold on; box on;
pbaspect([1.5 1 1])
ylabel({' '; '$\Delta \eta_{ts}$ -- Total-to-static effciency drop';' '});
ax_fig.XAxis.TickLabelFormat = '%0.2f';
ax_fig.YAxis.TickLabelFormat = '%0.3f';
ax_fig.XTick = 1:1:50;

% Prepare the axis labels
name = cell(2*n_stages,1);
for k = 1:n_stages
    name{2*k-1} = ['S$_{',num2str(k,'%0.0f'),'}$'];
    name{2*k+0} = ['R$_{',num2str(k,'%0.0f'),'}$'];
end
name{end+1} = 'Diffuser';
name{end+1} = ' ';  % Add one white label at the end
set(gca,'xticklabel',name)

% Plot the losses
bar(loss_matrix,0.5,'stacked')
legend('Profile loss','Secondary loss','Tip clearance loss','Trailing edge loss','Skin friction loss','Kinetic energy loss','Location','EastOutside')
% colormap(gray)

% Save the figure
name = '_loss_breakdown';
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