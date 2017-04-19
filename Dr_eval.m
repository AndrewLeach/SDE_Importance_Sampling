function [Dr] = Dr_eval(X0,X,parameters)

dt = parameters.dt;
sigma = parameters.sigma;
D = parameters.D;
N = length(X)/D;

Dr = spalloc(N*D+1,N*D,3*(D^2)*N + D);

block = 1:D;
Dr(block,block) = sqrt(dt/2)*(sigma\(eye(D)/dt));

for n = 2:N
    block = (n-1)*D+1:n*D;
    Dr(block,block) = sqrt(dt/2)*(sigma\(eye(D)/dt));
    Dr(block,block-D) = sqrt(dt/2)*(sigma\((-eye(D)/dt-Df_eval(X(block-D),parameters))));
end

if(strcmp(parameters.g,'normal'))
    g = g_eval(X(block),parameters);
    Dg = Dg_eval(X(block),parameters);
    if(g>0)
        Dr(N*D+1,block) = Dg/(2*sqrt(g));
    else
        Dr(N*D+1,block) = zeros(1,D);
        display('zero g');
    end
end

end