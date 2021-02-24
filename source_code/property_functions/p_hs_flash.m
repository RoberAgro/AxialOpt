function p = p_hs_flash(h,s,fluid,p_min,p_max)

% p = prop_calculation('P','H',h,'S',s,fluid);

try
    % Refprop default h-s flash function
    % It is faster, but not so robust
    % It might fail for some substances and states
    p = prop_calculation('P','H',h,'S',s,fluid);
    
catch  
    % This alternative uses the maximum and minimum expected pressures to
    % help the convergence of the flash calculation
    % It is more robust than the refprop h-s flash computation but it is 
    % more computationally expensive (factor of 10)
    options = optimset('TolX',1e-5);
    rho_error = @(h,s,fluid,p) prop_calculation('D','P',p,'H',h,fluid)-prop_calculation('D','P',p,'S',s,fluid);
    p = fzero(@(p) rho_error(h,s,fluid,p), [p_min p_max],options);
%     disp('Using long h-s flash calculations')

end
