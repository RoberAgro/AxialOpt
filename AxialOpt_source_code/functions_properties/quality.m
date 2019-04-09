function q = quality(prop2,prop2_value,fluid,p,p_crit)
%% Compute the vapor quality of a thermodynamic state
% Author: Roberto Agromayor

% Give an output for quality even if the pressure is supercritical
if p < p_crit
    q = refpropm('q','p',p,prop2,prop2_value,fluid);
elseif p >= p_crit
    q = 1.1;
else
    error('Ooops, something went wrong');
end

end