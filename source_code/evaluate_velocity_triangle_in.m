function [u,v,v_t,v_m,w,w_t,w_m,alpha,beta] = evaluate_velocity_triangle_in(u,v,alpha)

% Compute the velocity triangle at the inlet of the cascade
% Input:
% 1) Blade speed
% 2) Absolute velocity
% 3) Absolute flow angle

% Absolute velocities
v_t = v*sin(alpha);
v_m = v*cos(alpha);

% Relative velocities
w_t = v_t-u;
w_m = v_m;
w = sqrt(w_t^2+w_m^2);

% Relative flow angle
beta = atan(w_t/w_m);

end