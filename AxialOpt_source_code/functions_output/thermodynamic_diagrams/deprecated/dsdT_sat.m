function xi = dsdT_sat(T,fluid)
% Slope of the vapor saturation line in the T-s diagram
% Author: Roberto Agromayor
dT = T/1000;
s2 = refpropm('s','T',T+dT/2,'q',1,fluid);
s1 = refpropm('s','T',T-dT/2,'q',1,fluid);
xi = (s2-s1)/dT;  % Greek letter xi
end

