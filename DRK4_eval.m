%% DRK4_eval.m

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

function [Df] = DRK4_eval(x,parameters)

dt = parameters.dt;
D = length(x);

k1 = b_eval(x,parameters);
k2 = b_eval(x + dt*(k1/2),parameters);
k3 = b_eval(x + dt*(k2/2),parameters);

Dk1 = Db_eval(x,parameters);
Dk2 = Db_eval(x + dt*(k1/2),parameters)*(eye(D) + dt*(Dk1/2));
Dk3 = Db_eval(x + dt*(k2/2),parameters)*(eye(D) + dt*(Dk2/2));
Dk4 = Db_eval(x + dt*(k3),parameters)*(eye(D) + dt*(Dk3));

Df = Dk1/6 + Dk2/3 + Dk3/3 + Dk4/6;

end