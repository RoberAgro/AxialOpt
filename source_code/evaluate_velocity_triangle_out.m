function [u,v,v_t,v_m,w,w_t,w_m,alpha,beta] = evaluate_velocity_triangle_out(u,w,beta)

% Compute the velocity triangle at the outlet of the cascade
% Input:
% 1) Blade speed
% 2) Relative velocity
% 3) Relative flow angle

% Relative velocities
w_t = w*sin(beta);
w_m = w*cos(beta);

% Absolute velocities
v_t = w_t+u;
v_m = w_m;
v = sqrt(v_t^2+v_m^2);

% Absolute flow angle
alpha = atan(v_t/v_m);

end