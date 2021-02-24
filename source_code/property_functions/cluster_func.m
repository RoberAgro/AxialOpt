function f = cluster_func(z,beta)

% The cluster function f(z) increases the point density when z -> 0
% The amount of clustering can be controled with the parameter beta>1

f = 1 + beta*(1-((beta+1)/(beta-1)).^(1-z))./(1+((beta+1)/(beta-1)).^(1-z));

% If you want more points when z -> 1, then define the inverse function and
% reverse the resulting vector

end