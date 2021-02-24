function [] = plot_view_axial_radial_diffuser(turbine_data,my_filename,my_path,save)

% Axial-Radial view of the turbine including the exhaust diffuser

% The code for this plot can be difficult to understand
% The geometry is defined stage by stage (cascade by cascade is simpler)
% The function 'hatchfill' is used to get different hatched patters for the
% static and moving parts of the turbine
% Set symmetry = -1 to see the whole turbine
% Set symmetry = 1 to see onlyy half of the turbine

% Load the required variables
b          = [turbine_data.cascade.b];
cs         = [turbine_data.cascade.cs];
H          = [turbine_data.plane.H];
t_cl       = [turbine_data.cascade.t_cl];
r_m        = turbine_data.overall.mean_radius;
L_turb     = turbine_data.overall.axial_length;
n_stages   = turbine_data.overall.n_stages;

% Diffuser
L_diff = turbine_data.diffuser.axial_length;
r_1_in = turbine_data.diffuser.r_1_in;
r_2_in = turbine_data.diffuser.r_2_in;
r_1_out = turbine_data.diffuser.r_1_out;
r_2_out = turbine_data.diffuser.r_2_out;



% Define the symmetry parameters
symmetry = -1;

% Define my metal color
my_grey = 0.7*[1 1 1];

% Prepare the figure
fig = figure;
hold on; axis image; axis off;

% Plot the cascades
x_0 = 0;
for j = 1:n_stages
    
    % Prepare the geometry
    stator_x = [        0          0   b(2*j-1)    b(2*j-1)         0];
    stator_y = [-H(4*j-3)   H(4*j-3)   H(4*j-2)   -H(4*j-2) -H(4*j-3)]/2+r_m;
    rotor_x  = [b(2*j-1)+cs(2*j-1)   b(2*j-1)+cs(2*j-1)   b(2*j-1)+cs(2*j-1)+b(2*j)   b(2*j-1)+cs(2*j-1)+b(2*j)   b(2*j-1)+cs(2*j-1)];
    rotor_y  = [-H(4*j-1)                      H(4*j-1)                      H(4*j)                     -H(4*j)            -H(4*j-1)]/2+r_m;
    
    % Plot the geometry
    plot(stator_x+x_0,stator_y,'k','LineWidth',0.75);
    plot(rotor_x+x_0,rotor_y,'k','LineWidth',0.75);
    plot(stator_x+x_0,symmetry*stator_y,'k','LineWidth',0.75);
    plot(rotor_x+x_0,symmetry*rotor_y,'k','LineWidth',0.75);
    x_0 = x_0+(b(2*j-1)+cs(2*j-1)+b(2*j)+cs(2*j));
    
end


% Plot the disks and the first stator
x_0 = 0;
dy = 0.15*(r_2_in-r_1_in);  % Good thickness unit
cc = 0.25;                         % Shaft gap fraction of the radius
dd = 0.20;                         % Shaft fraction of the radius
for j = 1:n_stages
    
    % Plot the rotor disks
    rotor_cl_x = [b(2*j-1)+cs(2*j-1)   b(2*j-1)+cs(2*j-1)+b(2*j)   b(2*j-1)+cs(2*j-1)+b(2*j)   b(2*j-1)+cs(2*j-1)   b(2*j-1)+cs(2*j-1)];
    rotor_cl_y = [-H(4*j-1)/2+r_m                  -H(4*j)/2+r_m                           0                    0      -H(4*j-1)/2+r_m];
    patch(rotor_cl_x+x_0,rotor_cl_y,my_grey,'LineWidth',0.50);
    patch(rotor_cl_x+x_0,symmetry*rotor_cl_y,my_grey,'LineWidth',0.50);
    
    % Plot the first stator
    if j == 1
        first_stator_x = [          -3*dy                 0          b(2*j-1)             b(2*j-1)          b(2*j-1)-dy   b(2*j-1)-dy    -3*dy       -3*dy   b(2*j-1)-2*dy        b(2*j-1)-2*dy                -3*dy             -3*dy];
        first_stator_y = [-H(4*j-3)/2+r_m   -H(4*j-3)/2+r_m   -H(4*j-2)/2+r_m   -H(4*j-2)/2+r_m-dy   -H(4*j-3)/2+r_m-dy        cc*r_m   cc*r_m   cc*r_m+dy       cc*r_m+dy   -H(4*j-3)/2+r_m-dy   -H(4*j-3)/2+r_m-dy   -H(4*j-3)/2+r_m];
        patch(first_stator_x,first_stator_y,my_grey,'LineWidth',0.50);
        patch(first_stator_x,symmetry*first_stator_y,my_grey,'LineWidth',0.50);
    end

    x_0 = x_0+(b(2*j-1)+cs(2*j-1)+b(2*j)+cs(2*j));
   
end


% Plot the hub contour
x_0 = 0;
for j = 1:n_stages
    hub_x(1+4*(j-1):4*j) = x_0 + [                     0                 b(2*j-1)     b(2*j-1)+cs(2*j-1)   b(2*j-1)+cs(2*j-1)+b(2*j)];
    hub_y(1+4*(j-1):4*j) = r_m - [H(4*j-3)/2+t_cl(2*j-1)   H(4*j-2)/2+t_cl(2*j-1)             H(4*j-1)/2                    H(4*j)/2]; 
    x_0 = x_0+(b(2*j-1)+cs(2*j-1)+b(2*j)+cs(2*j));
end
hub_x = hub_x(3:end);
hub_y = hub_y(3:end);
hub_x = [hub_x hub_x(end:-1:1)    hub_x(1)];
hub_y = [hub_y hub_y(end:-1:1)-dy hub_y(1)];
patch(hub_x,hub_y,my_grey,'LineWidth',0.50);
patch(hub_x,hub_y*symmetry,my_grey,'LineWidth',0.50);


% Plot the shaft
t = linspace(0,2*pi,100); w = 0.025*r_m;
curvy_line_x = w*sin(t);
curvy_line_y = linspace(-dd,dd,100)*r_m;
x_shaft = [-5*dy+curvy_line_x  L_turb+7*dy+curvy_line_x(end:-1:1) -5*dy];
y_shaft = [curvy_line_y  curvy_line_y(end:-1:1) -dd*r_m];
patch(x_shaft,y_shaft,my_grey,'LineWidth',0.50);

t = linspace(0,2*pi,100); w = 0.025*r_m;
curvy_line_x = w*sin(t+pi);
curvy_line_y = linspace(-dd,dd,100)*r_m;
x_shaft = [-5*dy+curvy_line_x  L_turb+7*dy+curvy_line_x(end:-1:1) -5*dy];
y_shaft = [curvy_line_y  curvy_line_y(end:-1:1) -dd*r_m];
patch(x_shaft,y_shaft,my_grey,'LineWidth',0.50);


% Plot the casing
r_2_in = r_2_in + t_cl(2*n_stages);
x_0 = 0;
for j = 1:n_stages
    casing_x(1+4*(j-1):4*j) = x_0 + [         0     b(2*j-1)     b(2*j-1)+cs(2*j-1)   b(2*j-1)+cs(2*j-1)+b(2*j)];
    casing_y(1+4*(j-1):4*j) = r_m + [H(4*j-3)/2   H(4*j-2)/2   H(4*j-1)/2+t_cl(2*j)          H(4*j)/2+t_cl(2*j)];
    x_0 = x_0+(b(2*j-1)+cs(2*j-1)+b(2*j)+cs(2*j));
end


casing_x = [-3*dy          casing_x   L_turb+2*dy   L_turb+L_diff     L_turb+L_diff        L_turb+2*dy             L_turb              -3*dy              -3*dy];
casing_y = [ casing_y(1)   casing_y        r_2_in         r_2_out        r_2_out+dy        r_2_in+1*dy        r_2_in+1*dy        r_2_in+1*dy        casing_y(1)];
patch(casing_x,casing_y,my_grey,'LineWidth',0.50);
patch(casing_x,symmetry*casing_y,my_grey,'LineWidth',0.50);


% Plot the outlet annulus
L_start = L_turb+dy;
% bearing_x = [     2*dy    2*dy    4*dy   4*dy   dy          dy];
% bearing_y = [r_1_in-dy   bb+dy   bb+dy     bb   bb   r_1_in-dy];
bearing_x = [     2*dy        2*dy        4*dy       4*dy       dy          dy];
bearing_y = [r_1_in-dy   cc*r_m+dy   cc*r_m+dy     cc*r_m   cc*r_m   r_1_in-dy];
outlet_x = [     0     2*dy    L_diff       L_diff        2*dy   bearing_x           0        0] + L_start;
outlet_y = [r_1_in   r_1_in   r_1_out   r_1_out-dy   r_1_in-dy   bearing_y   r_1_in-dy   r_1_in];
patch(outlet_x,outlet_y,my_grey,'LineWidth',0.50);
patch(outlet_x,symmetry*outlet_y,my_grey,'LineWidth',0.50);


% % % Try to draw a scale (it needs manual tweaking)
% % scale = 0.1;  % 10 cm
% % x_0 = 0;
% % for i = 1:10    
% %     plot([x_0 x_0 x_0+scale/10 x_0+scale/10 x_0]-0.20*L,[0 scale/50 scale/50 0 0]-2*r_m,'k')
% %     x_0 = x_0+scale/10;    
% % end
% % x_0 = 0;
% % for i = 1:5
% %     patch([x_0 x_0 x_0+scale/10 x_0+scale/10 x_0]-0.20*L,[0 scale/50 scale/50 0 0]-2*r_m,'k')
% %     x_0 = x_0+scale/5;
% % end
% % text(-0.20*L,-2.15*r_m,'0.00')
% % text(-0.20*L+scale*1.0,-2.15*r_m,'0.10  (m)')
% % plot([-0.20*L -0.20*L],[-2.0*r_m -2.0*r_m-scale/30],'k')
% % plot([-0.20*L+scale*1.0 -0.20*L+scale*1.0],[-2.0*r_m -2.0*r_m-scale/30],'k')


% Save the figure
name = 'view_axial_radial_diffuser';
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

