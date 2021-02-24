function [] = plot_view_axial_tangential(turbine_data,my_filename,my_path,save)

% Axial-Tangential view of the turbine

% The information of a mean-line turbine optimization is not enough to
% define the detailed shape of the blades
% Arbitrary blade profiles based on NURBS curves were used for the plotting
% The shape of the blades plotted with this function only has aethetic
% purposes and it should not be taken as a realistic blade

% Load the required variables
s          = [turbine_data.cascade.s];
cs         = [turbine_data.cascade.cs];
c          = [turbine_data.cascade.c];
b          = [turbine_data.cascade.b];
stagger    = [turbine_data.cascade.stagger];
theta_in   = [turbine_data.cascade.theta_in];
theta_out  = [turbine_data.cascade.theta_out];
t_max      = [turbine_data.cascade.t_max];
t_te       = [turbine_data.cascade.t_te];
L          = turbine_data.overall.axial_length;
n_stages   = turbine_data.overall.n_stages;
n_cascades = 2*n_stages;
reaction   = [turbine_data.stage.R];

% Number of blades plotted for the last rotor row
N_blades   = 5;

% Prepare the digure
fig = figure;
hold on; axis image; axis off;

% Define the axes limits
axis([-0.40*L 1.40*L 0.15*s(end) 0.80*(N_blades-1)*s(end)])

% Plot a bounding box
plot([0 0 L L 0], ...
     [0.15 0.80*(N_blades-1) 0.80*(N_blades-1) 0.15 0.15]*s(end), ...
      'k','LineWidth',0.50)

  
x = 0;
for j = 1:n_cascades
    
    % Define the cascade planes
    cascade_inlet  = x;
    cascade_outlet = cascade_inlet+c(j)*cos(stagger(j));
    
    % Space to the inlet of the next cascade
    x = cascade_outlet+cs(j);
    
    % Plot the cascade planes
    plot([cascade_inlet cascade_inlet],[-1 1]*s(end)*10,'k','LineWidth',0.50)
    plot([cascade_outlet cascade_outlet],[-1 1]*s(end)*10,'k','LineWidth',0.50)
    
end


x_0 = 0;
for j = 1:n_cascades
    
    % Plots based on the 2D parametrization that I developed for ParaBlade
    % Parametrization based on 8 engineering parameters
    % 1) Blade chord
    % 2) Stagger angle
    % 3) Leading edge metal angle
    % 4) Trailing edge metal angle
    % 5) Leading edge radius
    % 6) Trailing edge radius
    % 7) Leading edge wedge semiangle
    % 8) Trailing edge wedge semiangle
    
    
    if reaction(ceil(j/2)) < 0.25 && mod(j,2) == 0
        % Prepare blade_parameters structure
        % Parameters for an impulse blade (just aesthetics for rotor)
        blade_parameters.chord = c(j);
        blade_parameters.stagger = stagger(j);
        blade_parameters.theta_in = theta_in(j);
        blade_parameters.theta_out = theta_out(j);
        blade_parameters.wedge_in = 10.0*pi/180;
        blade_parameters.wedge_out = 10.0*pi/180;
        blade_parameters.radius_in = 1.0*t_te(j);
        blade_parameters.radius_out = 1.0*t_te(j);
        blade_parameters.d1 = 0.60;
        blade_parameters.d2 = 0.40;
        blade_parameters.d3 = 0.40;
        blade_parameters.d4 = 0.60;
                
    else
        % Prepare blade_parameters structure
        % Parameters for a reaction blade (just aesthetics)
        blade_parameters.chord = c(j);
        blade_parameters.stagger = stagger(j);
        blade_parameters.theta_in = theta_in(j);
        blade_parameters.theta_out = theta_out(j);
        blade_parameters.wedge_in = 10*pi/180;
        blade_parameters.wedge_out = 7.5*pi/180;
        blade_parameters.radius_in = 0.50*t_max(j);
        blade_parameters.radius_out = 1.0*t_te(j);
        blade_parameters.d1 = 0.30;
        blade_parameters.d2 = 0.30;
        blade_parameters.d3 = 0.30;
        blade_parameters.d4 = 0.30;
        
    end
    
    
    % Get the coordinates of the blade
    [x,y] = plot_blade_profile(blade_parameters);

    % Plot the blades for the current cascade
    for k = 1:N_blades*ceil(s(end)/s(j))
        plot(x+x_0,y+s(j)*(k-1),'k','LineWidth',0.50)
    end
    
    % Space to the inlet of the next cascade
    x_0 = x_0+b(j)+cs(j);
    
end


% Save the figure
name = 'view_axial_tangential';
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



% % Deprecated blade parametrization based on
% % 1) Circular arc blades
% % 2) NACA thickness distributions
% for j = 1:n_cascades
%     
%     % The blades are parametriced as circular arcs
%     % The parametrization is different for stator and rotor blades (pi rad)
%     % A thickness distribution is imposed normal to the mean line to
%     % determine the shape of the blade
%     % The thickness distribution corresponds to that of NACA 4-digits
%     % airfoils. The equation is obtained from the NACA report of 1935
%     % See handwritten notes for more information about blade definition
%     
%     % Define the blade mean lines
%     if strcmp(type{j},'stator') == 1
%         phi_in = theta_in(j)+3*pi/2;
%         phi_out = theta_out(j)+3*pi/2;
%     elseif strcmp(type{j},'rotor') == 1
%         phi_in = theta_in(j)+pi/2;
%         phi_out = theta_out(j)+pi/2;
%     end
%     
%     R = c(j)/sqrt((2*(1-cos(phi_out-phi_in))));
%     phi = linspace(phi_in,phi_out,500);
%     x = R*cos(phi)-R*cos(phi_in);
%     y = R*sin(phi)-R*sin(phi_in);
%     
%     % Define the blade thickess distribution
%     z = linspace(0,1,500);
%     t = (t_max(j)/0.20)*(0.2969*z.^0.5-0.1260*z-0.3516*z.^2+0.2843*z.^3-0.1036*z.^4);
%     x_u = x+t.*cos(phi);
%     x_l = x-t.*cos(phi);
%     y_u = y+t.*sin(phi);
%     y_l = y-t.*sin(phi);
% 
%     % Plot the blades for the current cascade
%     for k = 1:N_blades*ceil(s(end)/s(j))
%         plot(x_u+x_0,y_u+s(j)*(k-1),'k','LineWidth',0.50)
%         plot(x_l+x_0,y_l+s(j)*(k-1),'k','LineWidth',0.50)
%         % plot(x+x_0,y+s(j)*(k-1),'g','LineWidth',0.50)
%     end
%     
%     % Space to the inlet of the next cascade
%     x_0 = x_0+b(j)+cs(j);
%     
% end