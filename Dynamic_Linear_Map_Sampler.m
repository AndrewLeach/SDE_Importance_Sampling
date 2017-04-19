function [output,logW,logO,X] = Dynamic_Linear_Map_Sampler(phi,parameters)

output = struct();

N = parameters.N;
D = parameters.D;
dt = parameters.dt;
sigma = parameters.sigma;
ep = parameters.ep;

Xt = parameters.X0;    
X = zeros(N*D,1);    
Z = 0;
replaced_time = [];

etas = normrnd(zeros(N*D,1),ones(N*D,1));
for n = 1:N
    block = (n-1)*D + 1 : n*D;
    
    % Find the zero viscosity optimal path
    [val,phi,replaced] = phi_eval(1,phi,Xt,1e00,1e00,parameters);
    if(replaced==1)
        replaced_time = [replaced_time,n];
    end
    H = Hessian_eval([],[],Xt,phi,parameters);
    if(n<N)
        Hh = H(D+1:end,D+1:end);
        H1 = H(D+1:end,1:D);
        H11 = H(1:D,1:D);
        
        Ht = H11 - H1'*(Hh\H1);
        phit = phi(1:D);
        phih = phi(D+1:end);
    else
        Ht = H;
        phit = phi;
    end
    
    U = chol(Ht,'upper');    
    eta = etas(block);
    Xtp = Xt;
    Xt = (U/sqrt(ep))\eta + phit;
    psit = Xtp + dt*f_eval(Xtp,parameters);
    Sit = inv(dt*(sigma*sigma'));
    
    Rxt = (Xt-psit)'*Sit*(Xt-psit)/2;
    Gxt = (Xt-phit)'*Ht*(Xt-phit)/2;
    detSHt = det(Sit\Ht);
    Z = Z - (1/2)*log(detSHt) - (Rxt-Gxt)/ep;
    X(block) = Xt;
    
    if(n<N)
        %phi = phih;
        phi = phih - (Hh\H1)*(Xt-phit);
    end
end
output.replaced_time = replaced_time;
logW = Z;
logO = -g_eval(Xt,parameters)/ep;
end