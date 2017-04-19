function [output,logW,logO,X] = Linear_Map_Sampler(phi,H,parameters)

output = struct();

N = parameters.N;
D = parameters.D;
dt = parameters.dt;
ep = parameters.ep;
sigma = parameters.sigma;
X0 = parameters.X0;

if(isempty(H))
    H = Hessian_eval([],[],X0,phi,parameters);
end
U = chol(H,'upper');
detSH =  (det(dt*(sigma*sigma'))^N)*det(H);

block = (N-1)*D+1:N*D;
eta = normrnd(zeros(N*D,1),ones(N*D,1));
Y = (U/sqrt(ep))\eta;
X = Y + phi;
rx = r_eval(X0,X,parameters);
gx = g_eval(X(block),parameters);    
Rx = rx'*rx-gx;
Gx = (X-phi)'*H*(X-phi)/2;    

logW = -log(detSH)/2 - (Rx-Gx)/ep;
logO = -gx/ep;

end