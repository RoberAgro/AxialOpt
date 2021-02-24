function [p,s] = ps_grid(fluid,N_p,N_s,p_max,T_min,T_max)
%% Compute a grid of thermodynamic properties in the p-s plane
% Author: Roberto Agromayor
% This plane is convenient for Rankine power cycles


%% Initialize the problem variables
[T_crit, s_crit] = refpropm('Ts','C',0,' ',0,fluid);     % Critical point
p = zeros(N_p,N_s);                                      % Pressure
s = zeros(N_p,N_s);                                      % Entropy

%% Subcritical pressures
Np1 = round(N_p/2);
Np2 = N_p-Np1;
beta = 1.018;
z = linspace(0,1,Np1)';
T1 = 0.999*T_crit+cluster_func(z,beta)*(T_min-0.999*T_crit);
T1(1:end) = T1(end:-1:1);  % Reverse the vector

% Initialize more variables
s1 = zeros(Np1,1);
p1 = zeros(Np1,1);
p2 = zeros(Np2,1);

% Compute the vapor saturation line in the T-s diagram and move a small
% distance into the vapor region in the direction normal to the vapor
% saturation line
delta = 1; % Distance normal to the saturation line
for i = 1:Np1
    s1(i) = refpropm('s','T',T1(i),'q',1.0,fluid);
    xi    = dsdT_sat(T1(i),fluid);
    dT    = -xi/sqrt((xi^2+1));
    ds    = 1/sqrt((xi^2+1));
    T1(i) = T1(i)+delta*dT;         % Overwrite the values of T1
    s1(i) = s1(i)+delta*ds;         % Overwrite the values of s1
    p1(i) = refpropm('p','T',T1(i),'s',s1(i),fluid);
end

%% Supercritical pressures
beta = 1.20;
z = linspace(0,1,Np2)';
tau = refpropm('T','p',p_max,'s',s_crit,fluid);    % T at s_crit and p_max
T2 = T1(end)+cluster_func(z,beta)*(tau-T1(end));

% Smooth distribution of entropy from s1(end) to s_crit
beta = 1.500;
s2 = s_crit+cluster_func(z,beta)*(s1(end)-s_crit);
s2(1:end) = s2(end:-1:1);

for i = 1:Np2
    p2(i) = refpropm('p','T',T2(i),'s',s2(i),fluid);
end

%% Get the vectors with the minimum and maximum entropies
% Combine the subcritical and supercritical results
P = [p1; p2];
s_min = [s1; s2];
s_max = zeros(N_p,1);
for i = 1:N_p
    s_max(i) = refpropm('s','T',T_max,'p',P(i),fluid);
end

%% Define the grid
% Cluster the entropy points close to the saturation curve
z = linspace(0,1,N_s);
beta = 1.1;
for i = 1:N_p
    s(i,:) = s_min(i)+cluster_func(z,beta)*(s_max(i)-s_min(i));
end

for j = 1:N_s
    p(:,j) = P;
end


end