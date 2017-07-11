%% phi_eval.m

%% This file is part of SDE_Importance_Sampling.  This code
%% accompanies the paper "Symmetrized importance samplers
%% for stochastic differential equations" by Andrew Leach,
%% Kevin K Lin, and Matthias Morzfeld.  A preprint version
%% of the paper can be found at
%% http://arxiv.org/abs/1707.02695 .

%% Copyright (C) 2017 by Andrew Leach <imaleach@gmail.com>.

%% This program is free software; you can redistribute
%% it and/or modify it under the terms of the GNU
%% General Public License as published by the Free
%% Software Foundation; either version 2 of the
%% License, or (at your option) any later version.
%% 
%% This program is distributed in the hope that it
%% will be useful, but WITHOUT ANY WARRANTY; without
%% even the implied warranty of MERCHANTABILITY or
%% FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%% General Public License for more details.
%% 
%% You should have received a copy of the GNU General
%% Public License along with this program; if not, write
%% to the Free Software Foundation, Inc., 51 Franklin
%% Street, Fifth Floor, Boston, MA 02110-1301 USA.

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

