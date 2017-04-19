function [ value_best, phi_best, replaced ] = phi_eval(init,phi_init,Xt,num_init,std_noise,parameters )

    min_tol = 1e-04;
    N = parameters.N;
    D = parameters.D;
    dt = parameters.dt;
    n = N-length(phi_init)/D+1;
    replaced = 0;
    
    if(init==1)
        % Find the zero viscosity optimal path
        f=@(x)(opt_fun_full(Xt,x,parameters));
        options = optimoptions(@fminunc,'GradObj','on','Hessian','on','Algorithm','trust-region',...
        'Diagnostics','off',...
        'Display','off',...
        'MaxFunEvals',1e3,...
        'MaxIter',400,...
        'TolX',1e-15,...
        'TolFun',1e-15...
        );
        [phi_best,value_best,exitflag,output] = fminunc(f,phi_init,options);
    end
    
    if(num_init~=0)
        xi = normrnd(zeros(D*(N-n+1),num_init),ones(D*(N-n+1),num_init));
        for m = 1:num_init
            phi = zeros((N-n+1)*D,1);
            phi(1:D) = Xt + dt*f_eval(Xt,parameters);
            for k = 2:(N-n+1)
                block = (k-1)*D + 1 : k*D;
                phi(block) = phi(block-D) + dt*f_eval(phi(block-D),parameters) + sqrt(dt)*std_noise*xi(block,m);
            end

            % Find the zero viscosity optimal path

            f=@(x)(opt_fun_full(Xt,x,parameters));
            options = optimoptions(@fminunc,'GradObj','on','Hessian','on','Algorithm','trust-region',...
            'Diagnostics','off',...
            'Display','off',...
            'MaxFunEvals',1e3,...
            'MaxIter',400,...
            'TolX',1e-15,...
            'TolFun',1e-15...
            );
            [phi,fval,exitflag,output] = fminunc(f,phi,options);

            if((init==0)&&(m==1))
                phi_best = phi;
                value_best = fval;
            end

            if(value_best*(1-min_tol)>fval)
                    % Update with the better phi
                    value_best = fval;
                    phi_best = phi;
                    replaced = 1;
            end        
        end
    end
end

