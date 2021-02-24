function q = compute_quality(p,prop2,prop2_value,fluid)

% Give an output for quality even if the pressure is supercritical
p_crit = prop_calculation('P_CRITICAL',fluid);
if p < p_crit
    q = prop_calculation('Q','P',p,prop2,prop2_value,fluid);
elseif p >= p_crit
    q = 1.1;
else
    error('Ooops, something went wrong in quality computation');
end

% % Coolprop gives a negative value when the fluid is superheated
% % This is a dirty fix that should be avoided
% if q < 0
%     q = 1.1;

end