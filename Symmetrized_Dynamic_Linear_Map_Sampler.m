%% Symmetrized_Dynamic_Linear_Map_Sampler.m

%% This file is part of SDE_Importance_Sampling.  This code
%% accompanies the paper "Symmetrized importance samplers
%% for stochastic differential equations" by Andrew Leach,
%% Kevin K Lin, and Matthias Morzfeld.

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

function [output,logWp,logWm,logOp,logOm,Xp,Xm] = Symmetrized_Dynamic_Linear_Map_Sampler(phi,parameters)

output = struct();

M = parameters.M;
N = parameters.N;
D = parameters.D;
dt = parameters.dt;
sigma = parameters.sigma;
ep = parameters.ep;

Xp = zeros(N*D,1);
Xm = zeros(N*D,1);
replaced_time = [];

etas = normrnd(zeros(N*D,1),ones(N*D,1));
phistart = phi;
for mi = [-1 1]
    Xt = parameters.X0;    
    X = zeros(N*D,1);
    Z = 0;
    phi = phistart;
    for n = 1:N
        block = (n-1)*D + 1 : n*D;

        % Find the zero viscosity optimal path
        [val,phi,replaced] = phi_eval(1,phi,Xt,1e00,1e00,parameters);
        if(replaced==1)
            replaced_time = [replaced_time,[n;mi]];
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
        eta = mi*etas(block);
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
        if(mi == 1)
            Xp(block) = X(block);
        elseif(mi == -1)
            Xm(block) = X(block);
        end
    end
    if(mi == 1)
        logWp = Z;
        logOp = -g_eval(Xt,parameters)/ep;
    elseif(mi == -1)
        logWm = Z;
        logOm = -g_eval(Xt,parameters)/ep;
    end
end
output.replaced_time = replaced_time;
end