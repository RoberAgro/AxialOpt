function stagger = evaluate_stagger_correlation(angle_in,angle_out,type)

% Compute the stagger angle using the formula proposed by Kacker (1982)

% Import the inlet and outlet angles and the thickness to chord ratio
if strcmp(type,'stator') == 1
    angle_in  = -angle_in;                                                 % alpha_1 (change sign)
    angle_out = +angle_out;                                                % alpha_2
elseif strcmp(type,'rotor') == 1
    angle_in  = +angle_in;                                                 % beta_2
    angle_out = -angle_out;                                                % beta_3 (change sign)
else
    error('Specify the type of cascade: ''rotor'' or ''stator''')
end

% Convert to degrees
angle_in = 180/pi*angle_in;
angle_out = 180/pi*angle_out;

if angle_in < -30
    angle_in = -30;
end

if angle_in > 70
    angle_in = 70;
end

if angle_out < 40
    angle_out = 40;
end

if angle_out > 80
    angle_out = 80;
end

% Sample the figure from Kacker and Okappu
angle_in_data  = [-30 -10 10 30 50 70];
angle_out_data = [80 70 60  50  40];
stagger_data   = [69 60 52  45  38
                  68 56 46  37  28
                  67 52 38  27  19
                  64 45 27  16   6
                  62 38 17   4  -8
                  58 29  6 -10 -22]';

% Interpolate data from the figure
[angle_in_array, angle_out_array] = meshgrid(angle_in_data,angle_out_data);
stagger_interp = @(x,y) interp2(angle_in_array,angle_out_array,stagger_data,x,y,'linear');
stagger = stagger_interp(angle_in,angle_out);

% Correct the sign for the case of rotor
if strcmp(type,'stator') == 1
    stagger = stagger*pi/180;   
elseif strcmp(type,'rotor') == 1
    stagger = -stagger*pi/180; % change sign
end


end

