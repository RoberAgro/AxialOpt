function a = compute_speed_of_sound(p,prop2,prop2_value,fluid)

% Avoid speed of sound computation in the two-phase region
q = compute_quality(p,prop2,prop2_value,fluid);
if (q >= 0) && (q <= 1)
    a = prop_calculation('A','P',p,'Q',1,fluid);
else
    a = prop_calculation('A','P',p,prop2,prop2_value,fluid);
end

end