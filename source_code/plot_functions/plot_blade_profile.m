function [x,y] = plot_blade_profile(blade_parameters)

% Compute the shape of the blades using a 2D blade parametrization
% The parametrization is based on NURBS curves and 8 engineering parameters

% The set of engineering parameters:
% 1) Blade chord
% 2) Stagger angle
% 3) Leading edge camber semiangle
% 4) Leading edge wedge semiangle
% 5) Leading edge radius (thickness)
% 6) Trailing edge camber angle
% 7) Trailing edge wedge angle
% 8) Trailing edge radius (thickness)

% The control points are used to impose these engineering parameters on the
% blade shape and to be able to control with mode detail (additions control
% points) the upper and lower surfaces, if necessary

% The pressure distribution on the upper and lower surfaces is:
% Upper surface: Suction (rotor blades) - Pressure (stator blades)
% Lower surface: Pressure (rotor blades) - Suction (stator blades)

% The blade shape is defined connecting four points with NURBS curves
% Point 1 is the end of the leading edge on the upper surface
% Point 2 is the end of the leading edge on the lower surface
% Point 3 is the beginning of the trailing edge on the lower surface
% Point 4 is the end of the trailing edge on the upper surface

% The blade shape and its first and second derivatives are computed in
% order to check that the surface is smooth (G1 and G2)

% The function also calculates other engineering parameters that depend on
% the main variables such as the channel area distribution, the location
% and width of the throat, the uncovered turning, and the mean curvature of
% the uncovered fraction of the upper surface.


%% Load the blade parameters
chord      = blade_parameters.chord;       % Blade chord
stagger    = blade_parameters.stagger;     % Stagger angle
theta_in   = blade_parameters.theta_in;    % Inlet camber angle
theta_out  = blade_parameters.theta_out;   % Outlet camber angle
wedge_in   = blade_parameters.wedge_in;    % Leading edge wedge semiangle
wedge_out  = blade_parameters.wedge_out;   % Trailing edge wedge semiangle
radius_in  = blade_parameters.radius_in;   % Leading edge radius
radius_out = blade_parameters.radius_out;  % Trailing edge radius
d1         = blade_parameters.d1;          % Distance of the control point
d2         = blade_parameters.d2;          % Distance of the control point
d3         = blade_parameters.d3;          % Distance of the control point
d4         = blade_parameters.d4;          % Distance of the control point


%% Compute the location of the 4 connecting points
% Point 1
a1 = theta_in + wedge_in;
x1 = radius_in*(1-sin(a1));
y1 = radius_in*cos(a1);

% Point 2
a2 = theta_in - wedge_in;
x2 = radius_in*(1+sin(a2));
y2 = -radius_in*cos(a2);

% Point 3
a3 = theta_out + wedge_out;
x3 = chord*cos(stagger) - radius_out*(1-sin(a3));
y3 = chord*sin(stagger) - radius_out*cos(a3);

% Point 4
a4 = theta_out - wedge_out;
x4 = chord*cos(stagger) - radius_out*(1+sin(a4));
y4 = chord*sin(stagger) + radius_out*cos(a4);


%% Compute the circular arc connecting 1 and 2
% The circular arc is represented by a second order rational Bezier curve
% Compute the control points
P0 = [x1 y1]';                           % First control point
P2 = [x2 y2]';                           % Third control point
n = [-(P2(2)-P0(2)); (P2(1)-P0(1))];     % Normal direction
n = n/norm(n);                           % Unitary normal vector
a = wedge_in;                            % Circular arc semiangle
R = radius_in;                           % Circular arc radius
P1 = (P0+P2)/2 - R*cos(a)/tan(a)*n;      % Second control point (see notes)
P = [P0 P1 P2];                          % Matrix of control points

% Compute the variable n - number of control points minus one
[~,nn] = size(P);
n = nn-1;

% Weight of the control points
W = ones(n+1,1);    % Unitary weight
W(2) = sin(a);      % Special weight for the second control point

% Define the parameter vector
u = linspace(0,1,100);

% Define the order of the NURBS
% The order of the NURBS is p (order of the basis polynomials
% Linear (p = 1), Quadratic (p = 1), Cubic (p = 1), etc.
% Set p = n (number of control points one) to obtain a Bezier curve 
p = 2; 

% Definition of the knot vector (clamped NURBS)
% p+1 zeros, n-p equispaced points between 0 and 1, and p+1 ones
% m+1 points where m = n+p+1
U = [zeros(1,p),linspace(0,1,n-p+2),ones(1,p)];

% Compute the NURBS curve
r12   =  NURBS(P,p,U,u,W);


%% Compute the curve connecting 2 and 3
% Compute the control points
P0 = [x2 y2]';
P1 = [x2 y2]' + d2*chord*[cos(a2) sin(a2)]';
P2 = [x3 y3]' - d3*chord*[cos(a3) sin(a3)]';
P3 = [x3 y3]';
P = [P0 P1 P2 P3];           % Matrix of control points

% Compute the variable n - number of control points minus one
n = 3;

% Weight of the control points
W = ones(n+1,1);    % Unitary weight

% Define the parameter vector
u = linspace(0,1,100);

% Define the order of the NURBS
% The order of the NURBS is p (order of the basis polynomials
% Linear (p = 1), Quadratic (p = 1), Cubic (p = 1), etc.
% Set p = n (number of control points one) to obtain a Bezier curve 
p = 3; 

% Definition of the knot vector (clamped NURBS)
% p+1 zeros, n-p equispaced points between 0 and 1, and p+1 ones
% m+1 points where m = n+p+1
U = [zeros(1,p),linspace(0,1,n-p+2),ones(1,p)];

% Compute the NURBS curve
r23 = NURBS(P,p,U,u,W);


%% Compute the circular arc connecting 3 and 4
% The circular arc is represented by a second order rational Bezier curve
% Compute the control points
P0 = [x3 y3]';                           % First control point
P2 = [x4 y4]';                           % Third control point
n = [-(P2(2)-P0(2)); (P2(1)-P0(1))];     % Normal direction
n = n/norm(n);                           % Unitary normal vector
a = wedge_out;                           % Circular arc semiangle
R = radius_out;                          % Circular arc radius
P1 = (P0+P2)/2 - R*cos(a)/tan(a)*n;      % Second control point (see notes)
P =[P0 P1 P2];                           % Matrix of control points

% Compute the variable n - number of control points minus one
[~,nn] = size(P);
n = nn-1;

% Weight of the control points
W = ones(n+1,1);    % Unitary weight
W(2) = sin(a);      % Special weight for the second control point

% Define the parameter vector
u = linspace(0,1,100);

% Define the order of the NURBS
% The order of the NURBS is p (order of the basis polynomials
% Linear (p = 1), Quadratic (p = 1), Cubic (p = 1), etc.
% Set p = n (number of control points one) to obtain a Bezier curve 
p = 2; 

% Definition of the knot vector (clamped NURBS)
% p+1 zeros, n-p equispaced points between 0 and 1, and p+1 ones
% m+1 points where m = n+p+1
U = [zeros(1,p),linspace(0,1,n-p+2),ones(1,p)];

% Compute the NURBS curve
r34 = NURBS(P,p,U,u,W);


%% Compute the curve connecting 4 and 1
% Compute the variable n - number of control points minus one
% [~,nn] = size(P);
% n = nn-1;
n = 3;

% Weight of the control points
W = ones(n+1,1);    % Unitary weight

% Define the parameter vector
u = linspace(0,1,100);

% Define the order of the NURBS
% The order of the NURBS is p (order of the basis polynomials
% Linear (p = 1), Quadratic (p = 1), Cubic (p = 1), etc.
% Set p = n (number of control points one) to obtain a Bezier curve 
p = 3; 

% Definition of the knot vector (clamped NURBS)
% p+1 zeros, n-p equispaced points between 0 and 1, and p+1 ones
% m+1 points where m = n+p+1
U = [zeros(1,p),linspace(0,1,n-p+2),ones(1,p)];

% Compute the control points
P0 = [x4 y4]';
P1 = [x4 y4]' - d1*chord*[cos(a4) sin(a4)]';
P2 = [x1 y1]' + d4*chord*[cos(a1) sin(a1)]';
P3 = [x1 y1]';
P =[P0 P1 P2 P3];           % Matrix of control points

% Compute the NURBS curve
r41 = NURBS(P,p,U,u,W);


%% Combine the blade segments and prepare output
r = [r12; r23; r34; r41];
x = r(:,1);
y = r(:,2);


end


function [r,R,N] = NURBS(P,p,U,u,W)
% Information
% Author: Roberto Agromayor
% Date: 15/09/2018
% Description:
% This function computes the clamped NURBS of order p defined by the n+1
% control points contained in matrix P and weight vector W.
% Matrix P is given by:
% P = [p_0 p1_ p2_ ... p_n] (n+1 columns)
% Where p_i are the position vectors of the control points.
% p_i are given as column vectors with any number of dimensions.
% U is the vector of knot points given as an input
% The curve r(u) is parametrized by the variable u given as an input vector

% Force u to be a column vector
if iscolumn(u) ~= 1
    u = u(:);
end

% Force w to be a column vector
if iscolumn(W) ~= 1
    W = W(:);
end

% Preliminary computations
Nt = length(u);         % Number of points where the b-spline is evaluated
[~,nn] = size(P);       % Number of control points
n = nn-1;               % Variable n in the notes   
N = zeros(Nt,n);        % Initialize the array of basis polynomials

% Number of basis polynomials in the current step of the recursion
m = n+p+1;

% First step of the recursion formula (p = 0)
% The non-strict inequality (<=) is kept in the second term to avoid
% problems in the last control point
% This does not affect the rest of the curve
for i = 1:m
    N(:,i) = 0+1.*(u >= U(i)).*(u <= U(i+1));
end

% Second and next steps of the recursion formula (p = 1, 2, ...)
p = 1;
while m > n+1    % Advance the recursion until there are n+1 polynomials
    
    m = m-1;                       % Update the number of polynomials
    N_bis = zeros(Nt,m);           % Initialize the dummy variable
    
    for i = 1:m        
        n1 = (u-U(i))/(U(i+p)-U(i)).*N(:,i);
        n2 = (U(i+p+1)-u)/(U(i+p+1)-U(i+1)).*N(:,i+1);
        n1(isnan(n1)) = 0;         % Trick to avoid NaN (inf*zero)
        n2(isnan(n2)) = 0;         % Trick to avoid NaN (inf*zero)
        N_bis(:,i) = n1 + n2;      % Recursion formula
    end
    
    N = N_bis;                     % Update the array of basis polynomials
    p = p+1;                       % Update the order of the polynomials
    
end

% Compute the NURBS curve given by r(u)
% The summation is performed exploiting matrix multiplication
% Convert the W vector into a Nt x (n+1) matrix
W = ones(Nt,1)*W';   
R = (N.*W)./sum(N.*W,2);
r = R*P';

end


