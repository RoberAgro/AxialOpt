function [] = plot_loss_breakdown(turbine_data,my_filename,my_path,save)

% Plot the loss distribution in each cascade

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
ylabel({' '; '$\Delta \eta_{ts}$ -- Efficiency drop';' '});
ax_fig.XAxis.TickLabelFormat = '%0.2f';
ax_fig.YAxis.TickLabelFormat = '%0.3f';
ax_fig.XTick = 1:1:50;

% Prepare the axis labels
diffuser_model = turbine_data.diffuser.diffuser_model;
if strcmp(diffuser_model,'1D') || strcmp(diffuser_model,'isentropic')
    
    name = cell(2*n_stages,1);
    for k = 1:n_stages
        name{2*k-1} = ['S$_{',num2str(k,'%0.0f'),'}$'];
        name{2*k+0} = ['R$_{',num2str(k,'%0.0f'),'}$'];
    end
    name{end+1} = 'Diffuser';
    name{end+1} = ' ';  % Add one white label at the end
    set(gca,'xticklabel',name)
    b = bar(loss_matrix,0.5,'stacked','FaceColor','flat');
    legend('Profile loss','Secondary loss','Tip clearance loss','Trailing edge loss','Skin friction loss','Kinetic energy loss','Location','EastOutside')
    
elseif strcmp(diffuser_model,'no')
        
    name = cell(2*n_stages,1);
    for k = 1:n_stages
        name{2*k-1} = ['S$_{',num2str(k,'%0.0f'),'}$'];
        name{2*k+0} = ['R$_{',num2str(k,'%0.0f'),'}$'];
    end
    name{end+1} = 'Exit';
    name{end+1} = ' ';  % Add one white label at the end
    set(gca,'xticklabel',name)
    b = bar(loss_matrix,0.5,'stacked','FaceColor','flat');
    legend('Profile loss','Secondary loss','Tip clearance loss','Trailing edge loss','Kinetic energy loss','Location','EastOutside')
   
else
    error("The option diffuser_model must be 'no', '1D' or 'isentropic'")
end

% Adjust colormap
map = colormap('gray');
map = 0.30 + 0.60*map;
colormap(map);
for k = 1:size(loss_matrix,2)
    b(k).CData = k;
end

% Save the figure
name = 'loss_breakdown';
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
