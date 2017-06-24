%% Standard_Sampler.m

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

function [output,logW,logO,X] = Standard_Sampler(parameters)

output = struct();

N = parameters.N;
D = parameters.D;
dt = parameters.dt;
sigma = parameters.sigma;
ep = parameters.ep;
    
X = zeros(N*D,1);    
Xt = parameters.X0;

etas = normrnd(zeros(N*D,1),ones(N*D,1));
for n = 1:N
    block = (n-1)*D + 1 : n*D;
    Xt = Xt + dt*f_eval(Xt,parameters);
    eta = etas(block);
    Xt = Xt + sqrt(dt*ep)*sigma*eta;
    X(block) = Xt;
end

logW = 0;
logO = -g_eval(Xt,parameters)/ep;

end