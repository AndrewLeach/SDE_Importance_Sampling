function [output,logW,logO,X] = Standard_Sampler(parameters)

output = struct();

N = parameters.N;
D = parameters.D;
dt = parameters.dt;
sigma = parameters.sigma;
ep = parameters.ep;
    
X = zeros(N*D,1);    
Xt = parameters.X0;

etas = normrnd(zeros(N*D,1),ones(N*D,1));
for n = 1:N
    block = (n-1)*D + 1 : n*D;
    Xt = Xt + dt*f_eval(Xt,parameters);
    eta = etas(block);
    Xt = Xt + sqrt(dt*ep)*sigma*eta;
    X(block) = Xt;
end

logW = 0;
logO = -g_eval(Xt,parameters)/ep;

end