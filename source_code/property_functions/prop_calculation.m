function prop = prop_calculation(varargin)

% CoolProp function call
prop = py.CoolProp.CoolProp.PropsSI(varargin{:});

% REFPROP function call
% prop = refpropm(varargin{:});

end