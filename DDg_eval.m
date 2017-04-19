function [DDg] = DDg_eval(x,parameters)
if(strcmp(parameters.g,'normal'))
g_sigma = parameters.g_sigma;

DDg = inv(g_sigma*(g_sigma'));
end
end
