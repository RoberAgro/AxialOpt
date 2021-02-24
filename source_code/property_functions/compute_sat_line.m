function [X_sat, Y_sat] = compute_sat_line(fluid,T_min,T_max,x,y,N)

% Compute the liquid and the vapor saturation lines
[X_liq, Y_liq] = compute_liq_line(fluid,T_min,T_max,x,y,ceil(N/2));
[X_vap, Y_vap] = compute_vap_line(fluid,T_min,T_max,x,y,floor(N/2));

% Export the results
X_sat = [X_liq; X_vap];    % Create a single vector with property x
Y_sat = [Y_liq; Y_vap];    % Create a single vector with property y

end


function [X_vap, Y_vap] = compute_vap_line(fluid,T_min,T_max,x,y,N_vap)

% Clustering parameter (>1)
beta = 1.0010;               

% Vapor saturation line
z_vap = linspace(0,1,N_vap)';
T_vap = T_max+cluster_func(z_vap,beta)*(T_min-T_max);
X_vap = zeros(N_vap,1);                 % Pre-allocate space
Y_vap = zeros(N_vap,1);                 % Pre-allocate space
for i = 1:N_vap    
    X_vap(i) = prop_calculation(x,'T',T_vap(i),'Q',1,fluid);    
    Y_vap(i) = prop_calculation(y,'T',T_vap(i),'Q',1,fluid);
end

end


function [X_liq, Y_liq] = compute_liq_line(fluid,T_min,T_max,x,y,N_liq)

% Clustering parameter (>1)
beta = 1.0010;

% Liquid saturation line
z_liq = linspace(0,1,N_liq)';
T_liq = T_max+cluster_func(z_liq,beta)*(T_min-T_max);
T_liq(1:end) = T_liq(end:-1:1);         % Reverse the vector (clustering)
X_liq = zeros(N_liq,1);                 % Pre-allocate space
Y_liq = zeros(N_liq,1);                 % Pre-allocate space
for i = 1:N_liq
    X_liq(i) = prop_calculation(x,'T',T_liq(i),'Q',0,fluid);
    Y_liq(i) = prop_calculation(y,'T',T_liq(i),'Q',0,fluid);    
end

end

