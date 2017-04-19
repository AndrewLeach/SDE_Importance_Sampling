function [Dg] = Dg_eval(x,parameters)
if(strcmp(parameters.g,'normal'))
g_mu = parameters.g_mu;
g_sigma = parameters.g_sigma;

Dg = (g_sigma*(g_sigma'))\(x-g_mu);
Dg = Dg';
end
end
