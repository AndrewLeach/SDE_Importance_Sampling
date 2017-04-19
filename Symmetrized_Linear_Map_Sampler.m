function [output,logWp,logWm,logOp,logOm,Xp,Xm] = Symmetrized_Linear_Map_Sampler(phi,H,parameters)

output = struct();

N = parameters.N;
D = parameters.D;
dt = parameters.dt;
ep = parameters.ep;
sigma = parameters.sigma;
X0 = parameters.X0;

% Find the zero viscosity optimal path
if(isempty(H))
    H = Hessian_eval([],[],X0,phi,parameters);
end
U = chol(H,'upper');
detSH =  (det(dt*(sigma*sigma'))^N)*det(H);

eta = normrnd(zeros(N*D,1),ones(N*D,1));
Y = (U/sqrt(ep))\eta;
block = (N-1)*D + 1:N*D;
for mi = [-1 1]
    X =  mi*Y + phi;
    rx = r_eval(X0,X,parameters);
    gx = g_eval(X(block),parameters);
    Rx = rx'*rx-gx;
    Gx = (X-phi)'*H*(X-phi)/2;

    if(mi == 1)
        logWp = -log(detSH)/2 - (Rx-Gx)/ep;
        logOp = -gx/ep;
        Xp = X;
    elseif(mi == -1)
        logWm = -log(detSH)/2 - (Rx-Gx)/ep;
        logOm = -gx/ep;
        Xm = X;
    end
end

end