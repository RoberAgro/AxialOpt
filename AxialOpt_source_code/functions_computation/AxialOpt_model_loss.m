function turbine_data = AxialOpt_model_loss(k,turbine_data,parameters)
%% Evaluate the loss coefficient error
% Author: Roberto Agromayor

% Compute the loss coefficient
% 1) According to its definition
% 2) Using an empirical loss model
% The difference between loss coefficients will be driven to zero using an
% equality constraint for the optimization


%% Compute the loss coefficient using its definition
% Load the variables required for the computation
p0rel_in         = turbine_data.plane(2*k-1).p0rel;
p0rel_out        = turbine_data.plane(2*k).p0rel;
p_out            = turbine_data.plane(2*k).p;
s_in             = turbine_data.plane(2*k-1).s;
h_out            = turbine_data.plane(2*k).h;
s_out            = turbine_data.plane(2*k).s;
T_out            = turbine_data.plane(2*k).T;
w_out            = turbine_data.plane(2*k).w;
fluid            = turbine_data.overall.fluid;
h_out_s          = refpropm('h','p',p_out,'s',s_in,fluid);
loss_system      = parameters.design_input.loss_system;
loss_coefficient = parameters.design_input.loss_coefficient;

% Choose the definition of the loss coefficient
if strcmp(loss_coefficient,'p0') == 1
    % Stagnation pressure drop coefficient
    Y_guess = (p0rel_in-p0rel_out)/(p0rel_out-p_out);
elseif strcmp(loss_coefficient,'h') == 1   
    % Energy loss coefficient
    Y_guess = (h_out-h_out_s)/(w_out^2/2);
elseif strcmp(loss_coefficient,'s') == 1
    % Entropy loss coefficient
    Y_guess = (s_out-s_in)*T_out/(w_out^2/2);
else
    error('Choose a valid definition for the loss coefficient: ''p0'', ''h'', or ''s''')
end


%% Compute the loss coefficient from the correlations
cascade   = turbine_data.cascade(k);
plane_in  = turbine_data.plane(2*k-1);
plane_out = turbine_data.plane(2*k);

if strcmp(loss_system,'AM') == 1
    % Ainley-Mathieson loss system
    [Y_loss,Y_p,Y_s,Y_cl,Y_te] = loss_model_AM(cascade,plane_in,plane_out);
elseif strcmp(loss_system,'DC') == 1   
    % Dunham-Came loss system
    [Y_loss,Y_p,Y_s,Y_cl,Y_te] = loss_model_DC(cascade,plane_in,plane_out);
elseif strcmp(loss_system,'KO') == 1
    % Kacker-Okapuu loss system
    [Y_loss,Y_p,Y_s,Y_cl,Y_te] = loss_model_KO(cascade,plane_in,plane_out);
else
    error('Choose a valid loss system: ''AM'', ''DC'', or ''KO''')
end

Y_error = Y_guess-Y_loss;
fraction = [Y_p Y_s Y_cl Y_te]/Y_loss;


%% Store the computed variables in the turbine_data structure
turbine_data.cascade(k).Y_guess  = Y_guess;
turbine_data.cascade(k).Y_loss   = Y_loss;
turbine_data.cascade(k).Y_error  = Y_error;
turbine_data.cascade(k).Y_p      = Y_p;
turbine_data.cascade(k).Y_s      = Y_s;
turbine_data.cascade(k).Y_cl     = Y_cl;
turbine_data.cascade(k).Y_te     = Y_te;
turbine_data.cascade(k).fraction = fraction;


end