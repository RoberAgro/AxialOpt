function [h_isobar, s_isobar, T_isobar] = compute_isobar(p_value,fluid,N)

% Compute isobar curves in the h-s diagram and T-s diagrams

% Initialize variables
h_isobar = zeros(N,1);     % Pre-allocate space
s_isobar = zeros(N,1);     % Pre-allocate space

% Define temperature limits
T_min  = prop_calculation('T_TRIPLE',fluid)+5;
T_max  = prop_calculation('T_MAX',fluid)-5;
T_crit = prop_calculation('T_CRITICAL',fluid);
p_crit = prop_calculation('P_CRITICAL',fluid);

% Define temperature points

% Compute enthalpy and entropy

if p_value > p_crit
    
    % Subcritical isobar
    beta = 1.0010;
    z1 = linspace(0,1,ceil(N/2))';
    T1 = T_crit+cluster_func(z1,beta)*(T_min-T_crit);
    z2 = linspace(0,1,floor(N/2))';
    T2 = T_crit+cluster_func(z2,beta)*(T_max-T_crit);
    T_isobar = [T1(end:-1:1); T2];
    for i = 1:length(T_isobar)
        h_isobar(i) = prop_calculation('H','T',T_isobar(i),'P',p_value,fluid);
        s_isobar(i) = prop_calculation('S','T',T_isobar(i),'P',p_value,fluid);
    end

else
    
    % Supercritical isobar (start from saturated vapor)
    beta = 1.0010;
    z1 = linspace(0,1,N)';
    T_isobar = prop_calculation('T','P',p_value,'Q',1,fluid)+0.1;
    T_isobar = T_isobar+cluster_func(z1,beta)*(T_max-T_isobar);
    for i = 1:length(T_isobar)
        h_isobar(i) = prop_calculation('H','T',T_isobar(i),'P',p_value,fluid);
        s_isobar(i) = prop_calculation('S','T',T_isobar(i),'P',p_value,fluid);
    end
end


% Note that temperature is used as a parameter to get the hs points
% Useful because s increases monotically with T
% Other choice could fail with complex substances like MM
% See the h-s diagram of MM to see why it would fail

end