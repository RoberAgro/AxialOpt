function [Y,Y_p,Y_s,Y_cl,Y_te] = loss_model_DC(cascade,plane_in,plane_out)
%% Compute the loss coefficient using the Dunham-Came loss model
% Author: Roberto Agromayor


%% Load the variables of the current cascade
type       = cascade.type;
Re         = cascade.Re;
angle_in   = plane_in.beta;
angle_out  = plane_out.beta;
c          = cascade.c;
s          = cascade.s;
H          = cascade.H;
t_cl       = cascade.t_cl;
t_max      = cascade.t_max;
t_te       = cascade.t_te;
Ma_rel_out = plane_out.Ma_rel;


%% Compute the loss coefficients
% Mach number correction factor
f_Ma = mach_correction(Ma_rel_out);
% f_Ma = 1;     % Uncomment this line to ignore the supersonic penalty

% Reynolds number correction factor
f_Re = reynolds_correction(Re);

% Profile loss coefficient
Y_p = profile_loss(angle_in,angle_out,c,s,t_max);
Y_p = Y_p*f_Ma*f_Re;

% Secondary loss coefficient
Y_s = secondary_loss(angle_in,angle_out,H,c);
Y_s = Y_s*f_Re;

% Clearance loss coefficient
Y_cl = clearance_loss(type,angle_in,angle_out,H,c,t_cl);

% Trailing edge loss coefficient
y_te = trailing_loss(t_te,s);
Y_te = (Y_p+Y_s+Y_cl)*(y_te-1);

% Overall loss coefficient
Y = Y_p+Y_s+Y_cl+Y_te;


end


function Y_p = profile_loss(angle_in,angle_out,c,s,t_max)

% Compute the profile loss coefficient
% The absolute value of the outlet angle is used as input
% This is meaningless for stator cascades as angle_out>0
% However, for rotor cascades angle_out<0 and the absolute value is used to
% compute Y_p_reaction and Y_p_impulse according to Aungier correlation
% These formulas are valid for 40<abs(angle_out)<80
% Extrapolating outside of this limits might give completely wrong results
% If the optimization algorithm has upper and lower bounds for the outlet 
% angle there is no need to worry about this problem

Y_p_reaction  = nozzle_blades(s/c,abs(angle_out));
Y_p_impulse = impulse_blades(s/c,abs(angle_out));

% This formula works for both stator and rotor cascades
% The logic for this formula is quite tricky and I had to put a lot of
% thought into it. See my handwritten notes for information
% There is no need to change the sign of any angle in this formula
Y_p = Y_p_reaction-abs(angle_in/angle_out)*(angle_in/angle_out)*(Y_p_impulse-Y_p_reaction);

% Limit the extrapolation of the profile loss to avoid negative values for
% blade profiles with little deflection
% Low limit to 80% of the axial entry nozzle profile loss
% This value is completely arbitrary
Y_p = max(Y_p,0.80*Y_p_reaction);

% Avoid unphysical effect on the thickness by defining the variable aa
aa = max(0,-angle_in/angle_out);
Y_p = Y_p*((t_max/c)/0.20)^aa;

end

function Y_p_reaction = nozzle_blades(r_sc,angle_out)

% Use Aungier correlation to compute the pressure loss coefficient
phi = 90-angle_out*180/pi;
r_sc_min = (0.46+phi/77).*(phi < 30) + (0.614+phi/130).*(phi >= 30);
X = r_sc-r_sc_min;
A = (0.025+(27-phi)/530).*(phi < 27) + (0.025+(27-phi)/3085).*(phi >= 27);
B = 0.1583-phi/1640;
C = 0.08*((phi/30).^2-1);
n = 1+phi/30;
Y_p_reaction = (A+B.*X.^2+C.*X.^3).*(phi < 30) + (A+B.*abs(X).^n).*(phi >= 30);

end

function Y_p_impulse = impulse_blades(r_sc,angle_out)

% Use Aungier correlation to compute the pressure loss coefficient
phi = 90-angle_out*180/pi;
r_sc_min = 0.224+1.575*(phi/90)-(phi/90).^2;
X = r_sc-r_sc_min;
A = 0.242-phi/151+(phi/127).^2;
B = (0.30+(30-phi)/50).*(phi < 30) + (0.30+(30-phi)/275).*(phi >=30);
C = 0.88-phi/42.4+(phi/72.8).^2;
Y_p_impulse = A+B.*X.^2-C.*X.^3;

end

function f_Ma = mach_correction(Ma_rel_out)

f_Ma = 1+60*(Ma_rel_out-1).^2.*(Ma_rel_out > 1);

end


function f_Re = reynolds_correction(Re)

f_Re = (Re/2e5).^(-0.20);

end


function Y_s = secondary_loss(angle_in,angle_out,H,c)

% Compute the secondary loss coefficient
angle_m = atan((tan(angle_in)+tan(angle_out))/2);
Z       = 4*(tan(angle_in)-tan(angle_out))^2*cos(angle_out)^2/cos(angle_m);
Y_s     = 0.0334*(c/H)*cos(angle_out)/cos(angle_in)*Z;

end


function Y_cl = clearance_loss(type,angle_in,angle_out,H,c,t_cl)

if strcmp(type,'stator') == 1
    B = 0.00;                                                              % Empirical parameter for the stator
elseif strcmp(type,'rotor') == 1                               
    B = 0.37;                                                              % Empirical parameter for the rotor (choose 0.37 for shrouded blades)
else
    error('Specify the type of cascade: ''rotor'' or ''stator''')
end

angle_m = atan((tan(angle_in)+tan(angle_out))/2);
Z       = 4*(tan(angle_in)-tan(angle_out))^2*cos(angle_out)^2/cos(angle_m);
Y_cl    = B*(c/H)*(t_cl/H)^0.78*Z;

end


function y_te = trailing_loss(t_te,s)

% Import correlation data
r_te_s_data = [0.000 0.020 0.040 0.080 0.120];
y_te_data   = [0.914 1.000 1.105 1.375 1.680];

% Interpolate the value of y_te
y_te = interp1(r_te_s_data,y_te_data,t_te/s,'linear','extrap');

end


