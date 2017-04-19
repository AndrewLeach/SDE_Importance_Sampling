function [g] = g_eval(x,parameters)
if(strcmp(parameters.g,'normal'))
g_mu = parameters.g_mu;
g_sigma = parameters.g_sigma;

y = g_sigma\(x-g_mu);
g = (y'*y)/2;
end

end
