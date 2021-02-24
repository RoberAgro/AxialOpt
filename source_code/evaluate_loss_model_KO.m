function [Y,Y_p,Y_s,Y_cl,Y_te] = evaluate_loss_model_KO(cascade,plane_in,plane_out)

% Compute the loss coefficient using the Kacker-Okapuu loss model

%% Load the variables of the current cascade
type       = cascade.type;
Re         = cascade.Re;
angle_in   = plane_in.beta;
angle_out  = plane_out.beta;
c          = cascade.c;
b          = cascade.b;
o          = cascade.o;
s          = cascade.s;
H          = cascade.H;
t_cl       = cascade.t_cl;
t_max      = cascade.t_max;
t_te       = cascade.t_te;
Ma_rel_in  = plane_in.Ma_rel;
Ma_rel_out = plane_out.Ma_rel;
r_ht_in    = plane_in.r_ht;
p_in       = plane_in.p;
p0rel_in   = plane_in.p0rel;
p_out      = plane_out.p;
p0rel_out  = plane_out.p0rel;


%% Compute the loss coefficients
% Mach number correction factor
f_Ma = mach_correction(Ma_rel_out);
% f_Ma = 1;     % Uncomment this line to ignore the supersonic penalty

% Reynolds number correction factor
f_Re = reynolds_correction(Re);
% f_Re = 1;

% Profile loss coefficient
Y_p = profile_loss(type,angle_in,angle_out,c,s,t_max, ...
      Ma_rel_in,Ma_rel_out,r_ht_in,p_in,p0rel_in,p_out,p0rel_out);

% Corrected profile loss coeffcient
Y_p = Y_p*f_Ma*f_Re;

% Secondary loss coefficient
Y_s  = secondary_loss(angle_in,angle_out,Ma_rel_in,Ma_rel_out,H,c,b);

% Clearance loss coefficient
Y_cl = clearance_loss(type,angle_in,angle_out,H,c,t_cl);

% Trailing edge loss coefficient
Y_te = trailing_loss(angle_in,angle_out,t_te,o);

% Overall loss coefficient
Y = Y_p+Y_s+Y_cl+Y_te;


end


function Y_p = profile_loss(type,angle_in,angle_out,c,s,t_max,Ma_rel_in,Ma_rel_out,r_ht_in,p_in,p0rel_in,p_out,p0rel_out)

% Inlet shock loss
a = max(0,f_hub(r_ht_in,type)*Ma_rel_in-0.40);
Y_shock = 0.75*(a)^1.75*(r_ht_in)*(p0rel_in-p_in)/(p0rel_out-p_out);
% Avoid unphysical results if r_ht_in becomes negative during the
% optimization iterations
Y_shock = max(0,Y_shock);

% Compressible flow correction factor
% Limit excessively low values (it might be a problem during optimization)
% The limit to 0.1 is quite arbitrary, but it worked for me
Kp = max(0.1,K_p(Ma_rel_in,Ma_rel_out));

% Compute the profile loss coefficient
% The absolute value of the outlet angle is used as input
% This is meaningless for stator cascades as angle_out>0
% However, for rotor cascades angle_out<0 and the absolute value is used to
% compute Y_p_reaction and Y_p_impulse according to Aungier correlation
% These formulas are valid for 40<abs(angle_out)<80
% Extrapolating outside of this limits might give completely wrong results
% If the optimization algorithm has upper and lower bounds for the outlet 
% angle there is no need to worry about this problem
% angle_out_bis keeps the 40deg-losses for outlet angles lower than 40deg
angle_out_bis = max(abs(angle_out),40*pi/180);
Y_p_reaction  = nozzle_blades(s/c,angle_out_bis);
Y_p_impulse = impulse_blades(s/c,angle_out_bis);

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
Y_p = 0.914*(2/3*Y_p*Kp+Y_shock);

end


function Y_s = secondary_loss(angle_in,angle_out,Ma_rel_in,Ma_rel_out,H,c,b)

% Compute the secondary loss coefficient
% Limit excessively low values (it might be a problem during optimization)
% The limit to 0.1 is quite arbitrary, but it worked for me
Ks      = max(0.10,K_s(Ma_rel_in,Ma_rel_out,H/b));   
angle_m = atan((tan(angle_in)+tan(angle_out))/2);
Z       = 4*(tan(angle_in)-tan(angle_out))^2*cos(angle_out)^2/cos(angle_m);
Y_s     = 1.2*Ks*0.0334*f_AR(H/c)*cos(angle_out)/cos(angle_in)*Z;

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


function Y_te = trailing_loss(angle_in,angle_out,t_te,o)
% Reacting blading
r_to_data_reaction = [0.000 0.200 0.400];
phi_data_reaction  = [0.000 0.045 0.150];

% Impulse blading
r_to_data_impulse = [0.000 0.200 0.400];
phi_data_impulse  = [0.000 0.025 0.075];

% Numerical trick to avoid too big r_to's
r_to = min(0.400,t_te/o);

% Interpolate data
d_phi2_reaction = interp1(r_to_data_reaction,phi_data_reaction,r_to,'linear','extrap');
d_phi2_impulse  = interp1(r_to_data_impulse ,phi_data_impulse ,r_to,'linear','extrap');

d_phi2 = d_phi2_reaction-abs(angle_in/angle_out)*(angle_in/angle_out)*(d_phi2_impulse-d_phi2_reaction);

% Limit the extrapolation of the trailing edge loss
d_phi2 = max(d_phi2,d_phi2_impulse/2);
Y_te = 1/(1-d_phi2)-1;

end


function f = f_hub(r_ht,type)

if r_ht < 0.50
    r_ht = 0.50;        % Numerical trick to prevent extrapolation
end

% Stator curve
r_ht_data_S = [0.50 0.60 0.70 0.80 0.90 1.00];
f_data_S    = [1.40 1.18 1.05 1.00 1.00 1.00];

% Rotor curve
r_ht_data_R = [0.50 0.60 0.70 0.80 0.90 1.00];
f_data_R    = [2.15 1.70 1.35 1.12 1.00 1.00];

if strcmp(type,'stator') == 1
    f = interp1(r_ht_data_S,f_data_S,r_ht,'linear','extrap');
elseif strcmp(type,'rotor') == 1
    f = interp1(r_ht_data_R,f_data_R,r_ht,'linear','extrap');
else
    error('Specify the type of cascade: ''rotor'' or ''stator''')
end

end


function f = f_AR(x)

f = (1-0.25*sqrt(2-x))/x.*(x < 2) + 1/x.*(x >= 2);

end


function f_Ma = mach_correction(Ma_rel_out)

f_Ma = 1+60*(Ma_rel_out-1).^2.*(Ma_rel_out > 1);
% f_Ma = 1;

end


function f_Re = reynolds_correction(Re)

f_Re = (Re/2e5).^(-0.40).*(Re < 2e5) + 1*(Re >= 2e5 & Re <= 1e6) + (Re/1e6).^(-0.20).*(Re > 1e6);

end


function Y_p_reaction = nozzle_blades(r_sc,angle_out)

% Use Aungier correlation to compute the pressure loss coefficient
% This correlation is a formula that reproduces the figures from the Ainley
% and Mathieson original figures

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
% This correlation is a formula that reproduces the figures from the Ainley
% and Mathieson original figures

phi = 90-angle_out*180/pi;
r_sc_min = 0.224+1.575*(phi/90)-(phi/90).^2;
X = r_sc-r_sc_min;
A = 0.242-phi/151+(phi/127).^2;
B = (0.30+(30-phi)/50).*(phi < 30) + (0.30+(30-phi)/275).*(phi >=30);
C = 0.88-phi/42.4+(phi/72.8).^2;
Y_p_impulse = A+B.*X.^2-C.*X.^3;

end

function Kp = K_p(Ma_in,Ma_out)
Kp = 1-K_2(Ma_in/Ma_out)*(1-K_1(Ma_out));
% Kp = 1;
end


function Ks = K_s(Ma_in,Ma_out,r_Hb)
Ks = 1-K_3(1/r_Hb)*(1-K_p(Ma_in,Ma_out));
% Ks = 1;
end


function K1 = K_1(x)
K1 = 1*(x < 0.20) + (1-1.25*(x-0.20)).*(x > 0.20 & x < 1.00);
end


function K2 = K_2(x)
K2 = x^2;
end


function K3 = K_3(x)
K3 = x^2;
end




