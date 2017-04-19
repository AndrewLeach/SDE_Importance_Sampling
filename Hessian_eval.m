function [H] = Hessian_eval(Dr,r,X0,X,parameters)

dt = parameters.dt;
sigma = parameters.sigma;
D = parameters.D;
N = length(X)/D;

if(isempty(Dr) || isempty(r))
    r = r_eval(X0,X,parameters);
    Dr = Dr_eval(X0,X,parameters);
end

H = 2*(Dr'*Dr);

for n = 2:N
    block = (n-1)*D+1:n*D;
    DDf = DDf_eval(X(block-D),parameters);
    for d = 1:D
        H(block-D,block-D) = H(block-D,block-D) + 2*r((n-1)*D + d)*(-sqrt(dt/2))*(sigma\squeeze(DDf(d,:,:)));
    end
end

if(strcmp(parameters.g,'normal'))
    block = (N-1)*D+1:N*D;
    g = g_eval(X(block),parameters);
    Dg = Dg_eval(X(block),parameters);
    DDg = DDg_eval(X(block),parameters);
    
    if(g>0)
        H(block,block) = H(block,block) + 2*r(N*D+1)*(DDg/(2*sqrt(g)) - (Dg'*Dg)/(4*(g^(3/2))));
    else
        H(block,block) = zeros(D,D);
        display('zero g');
    end
end

end