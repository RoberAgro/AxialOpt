function p = p_hs_flash(h,s,fluid,p_min,p_max)
%% Function to do h-s flash computations in REFPROP
% Author: Roberto Agromayor

try
    % Refprop default h-s flash function
    % It is faster, but not so robust
    % It might fail for some substances and states
    p = refpropm('p','h',h,'s',s,fluid);
    
catch  
    % This alternative uses the maximum and minimum expected pressures to
    % help the convergence of the flash calculation
    % It is more robust than the refprop h-s flash computation but it is 
    % more computationally expensive (factor of 10)
    options = optimset('TolX',1e-5);
    rho_error = @(h,s,fluid,p) refpropm('d','p',p,'h',h,fluid)-refpropm('d','p',p,'s',s,fluid);
    p = fzero(@(p) rho_error(h,s,fluid,p), [p_min p_max],options);
    % disp('Using long h-s flash calculations')

end
