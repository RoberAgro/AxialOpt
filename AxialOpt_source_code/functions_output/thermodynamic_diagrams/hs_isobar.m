function [h_isobar, s_isobar, T_isobar] = hs_isobar(p,p_crit,T_crit,p_min,T_min,T_max,fluid,N,sat_line)
%% Compute the isobars in the h-s diagram or in the T-s diagram
% Author: Roberto Agromayor

% Initialize variables
Np = length(p);             % Number of isobars
h_isobar = zeros(N,Np);     % Pre-allocate space
s_isobar = zeros(N,Np);     % Pre-allocate space

for j = 1:Np
 
    if strcmp(sat_line,'yes') == 1
        
        if p(j) <= p_min
            T_isobar = linspace(T_min+0.01,T_max,N);
            for i = 1:length(T_isobar)
                h_isobar(i,j) = refpropm('h','T',T_isobar(i),'p',p(j),fluid);
                s_isobar(i,j) = refpropm('s','T',T_isobar(i),'p',p(j),fluid);
            end
        end
        
        if p(j) < p_crit && p(j) > p_min
            T_isobar = linspace(0.01+refpropm('T','p',p(j),'q',1,fluid),T_max,N);
            for i = 1:length(T_isobar)
                h_isobar(i,j) = refpropm('h','T',T_isobar(i),'p',p(j),fluid);
                s_isobar(i,j) = refpropm('s','T',T_isobar(i),'p',p(j),fluid);
            end
        end
        
        if p(j) >= p_crit
            T_isobar = linspace(T_crit,T_max,N);
            for i = 1:length(T_isobar)
                h_isobar(i,j) = refpropm('h','T',T_isobar(i),'p',p(j),fluid);
                s_isobar(i,j) = refpropm('s','T',T_isobar(i),'p',p(j),fluid);
            end
        end
        
    else
        
        T_isobar = linspace(T_min+10.01,T_max,N);
        for i = 1:length(T_isobar)
            h_isobar(i,j) = refpropm('h','T',T_isobar(i),'p',p(j),fluid);
            s_isobar(i,j) = refpropm('s','T',T_isobar(i),'p',p(j),fluid);
        end
        
        
    end
    

end

% Note that temperature is used as a parameter to get the hs points
% Useful because s increases monotically with T
% Other choice could fail with complex substances like MM
% See the hs diagram of MM to see why it would fail

end