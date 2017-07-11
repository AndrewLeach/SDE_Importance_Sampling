%% Dr_eval.m

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