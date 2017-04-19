function[] = Compare_Samplers_Pole_Reversal(batch,sampler,solver,phi_init,X0,sigma,eps,M,N,dt,D,mu,nu,gamma,g,g_mu,g_sigma)

%profile memory on;

S_CLOCK = clock;
rng(mod(ceil(S_CLOCK(6)*cputime*99999)*str2num(batch),2^31),'twister');

tic;
parameters = struct('solver',solver,'X0',X0,'sigma',sigma,'ep',0,'M',M,'N',N,...
                    'dt',dt,'D',D,'mu',mu,'nu',nu,'gamma',gamma,...
                    'g',g,'g_mu',g_mu,'g_sigma',g_sigma);
M = parameters.M;
L = length(eps);
Q = zeros(L,1);
logW = zeros(M,L);
logO = zeros(M,L);
logOW = zeros(M,L);
logWb = zeros(2*M,L);
logOb = zeros(2*M,L);
X = zeros(D*N,M,L);
Xp = zeros(D*N,M,L);
Xm = zeros(D*N,M,L);

if(~strcmp(sampler,'NONE'))
    if(isempty(phi_init))
        phi_init = zeros(N*D,1);
        [value,phi_init,replaced] = phi_eval(0,phi_init,X0,1e02,1e00,parameters);
        value
    end
    H_init = Hessian_eval([],[],X0,phi_init,parameters);
end
    
for l = 1:L
    parameters.ep = eps(l);
    for m = 1:M
        if(strcmp(sampler,'NONE'))
            [output,logW(m,l),logO(m,l),X(:,m,l)] = Standard_Sampler(parameters);
         elseif(strcmp(sampler,'DLM'))
            [output,logW(m,l),logO(m,l),X(:,m,l)] = Dynamic_Linear_Map_Sampler(phi_init,parameters);
        elseif(strcmp(sampler,'SDLM'))
            [output,logWb(m,l),logWb(M+m,l),logOb(m,l),logOb(M+m,l),Xp(:,m,l),Xm(:,m,l)] = Symmetrized_Dynamic_Linear_Map_Sampler(phi_init,parameters);
        elseif(strcmp(sampler,'LM'))
            [output,logW(m,l),logO(m,l),X(:,m,l)] = Linear_Map_Sampler(phi_init,H_init,parameters);
        elseif(strcmp(sampler,'SLM'))
            [output,logWb(m,l),logWb(M+m,l),logOb(m,l),logOb(M+m,l),Xp(:,m,l),Xm(:,m,l)] = Symmetrized_Linear_Map_Sampler(phi_init,H_init,parameters);
        end
        outputs(m,l) = output;
    end
    if(strcmp(sampler,'SDLM')||strcmp(sampler,'SLM'))
        logOW(:,l) = logOb(1:M,l) + logWb(1:M,l) + log(1+exp(logOb(M+1:2*M,l) + logWb(M+1:2*M,l)...
            - logOb(1:M,l) - logWb(1:M,l))) - log(2);
        logW(:,l) = logWb(1:M,l) + log(1+exp(logWb(M+1:2*M,l)-logWb(1:M,l)))-log(2);
        logO(:,l) = logOb(1:M,l) + logWb(1:M,l) + log(1+exp(logOb(M+1:2*M,l) + logWb(M+1:2*M,l)...
            - logOb(1:M,l) - logWb(1:M,l))) - logW(:,l) - log(2);

        preR = (logOW(:,l)) - max(logOW(:,l));
        R = mean(exp(preR).^2)/mean(exp(preR)).^2;
    else
        preR = (logW(:,l) + logO(:,l))-max(logW(:,l) + logO(:,l));
        R = mean(exp(preR).^2)/mean(exp(preR)).^2;
    end
    Q(l) = R-1;
end

elapsed_time = toc
output_file_name = strcat(sampler,'_',solver,'_',batch);
save(output_file_name);

end