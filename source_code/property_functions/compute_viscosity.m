function v = compute_viscosity(p,prop2,prop2_value,fluid)

% Avoid viscosity computation in the two-phase region
q = compute_quality(p,prop2,prop2_value,fluid);
if (q >= 0) && (q <= 1)
    v = prop_calculation('V','P',p,'Q',1,fluid);
else
    v = prop_calculation('V','P',p,prop2,prop2_value,fluid);
end

end