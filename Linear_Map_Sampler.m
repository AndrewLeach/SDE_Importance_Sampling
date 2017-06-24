%% Linear_Map_Sampler.m

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