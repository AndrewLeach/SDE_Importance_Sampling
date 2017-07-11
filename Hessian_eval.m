%% Hessian_eval.m

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