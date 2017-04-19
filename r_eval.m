function [r] = r_eval(X0,X,parameters)

dt = parameters.dt;
sigma = parameters.sigma;
D = parameters.D;
N = length(X)/D;

r = zeros(N*D+1,1);

block = 1:D;
r(block) = sqrt(dt/2)*(sigma\((X(block)-X0)/dt-f_eval(X0,parameters)));

for n = 2:N
    block = (n-1)*D+1:n*D;
    r(block) = sqrt(dt/2)*(sigma\((X(block)-X(block-D))/dt-f_eval(X(block-D),parameters)));
end

if(strcmp(parameters.g,'normal'))
    r(N*D+1) = sqrt(g_eval(X(block),parameters));
end
end