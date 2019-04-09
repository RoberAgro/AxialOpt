function prop = property_grid(fluid,p,s,name)
% Author: Roberto Agromayor
% Initialize the problem
[N_p,N_s] = size(p);
prop = zeros(N_p,N_s);

% Compute the property grid
for i = 1:N_p
    for j = 1:N_s
        
        if strcmp(name,'Gamma') == 1
            prop(i,j) = fundamental_derivative(p(i,j),s(i,j),fluid);
        else
            prop(i,j) = refpropm(name,'p',p(i,j),'s',s(i,j),fluid);
        end
        
    end
end