function t_max = evaluate_thickness_correlation(deflection,c)

% Compute the blade thickness using the formula proposed by Kacker (1982)

% This function is used to compute the blade thickness as a function of the
% real chord and the deflection of the flow
deflection = deflection*180/pi;

if deflection <= 40
    t_max = 0.15*c;
elseif deflection > 40 && deflection < 120
    t_max = (0.15 + (0.25-0.15)/(120-40)*(deflection-40))*c;
elseif deflection >= 120
    t_max = 0.25*c;
end

end

