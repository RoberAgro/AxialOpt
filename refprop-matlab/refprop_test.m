%% Test to check if REFPROP and MATLAB were correctly linked
% Author: Roberto Agromayor
% Contact: roberto.agromayor@ntnu.no
% Date: 15/08/2018

clear all
close all
clc

% Substance
fluid = 'water';                % Define the fluid as a string

% Pressure
pressure = 101.315;             % Ambient pressure in kPa

% Vapor quality
quality = 1;                    % Saturated steam


% Saturation temperature of water at ambient temperature (in Kelvin)
T_kelvin = refpropm('T','p',pressure,'q',quality,fluid);

% Saturation temperature of water at ambient temperature (in Celsius)
T_celsius = T_kelvin - 273.15;

% Print the result
disp(['The saturation temperature of water at ambient pressure is ',num2str(T_celsius),' Celsius degrees'])
disp('If you can read this message MATLAB and REFPROP are correctly linked')