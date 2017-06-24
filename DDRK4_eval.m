%% DDRK4_eval.m

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

function [DDf] = DDRK4_eval(x,parameters)

dt = parameters.dt;
D = length(x);

k1 = b_eval(x,parameters);
k2 = b_eval(x + dt*(k1/2),parameters);
k3 = b_eval(x + dt*(k2/2),parameters);

Dk1 = Db_eval(x,parameters);
Dk2 = Db_eval(x + dt*(k1/2),parameters)*(eye(D) + dt*(Dk1/2));
Dk3 = Db_eval(x + dt*(k2/2),parameters)*(eye(D) + dt*(Dk2/2));

% Dimension 1 carries the different dimensions of the vector field
% Dimensions 2 and 3 carry the derivatives wrt each of those dimensions.
% DDk(d,:,:) is the Hessian of the dth dimension of the vector field
DDk1 = DDb_eval(x,parameters);

Dbk1 = Db_eval(x+dt*(k1/2),parameters);
DDbk1 = DDb_eval(x+dt*(k1/2),parameters);
DDk2 = zeros(D,D,D);
Id = eye(D);
for d = 1:D
    for p = 1:D
        for l = 1:D
            for m = 1:D
                DDk2(d,p,l) = DDk2(d,p,l) + (dt/2)*Dbk1(d,m)*DDk1(m,p,l);
                for n = 1:D
                    DDk2(d,p,l) = DDk2(d,p,l) + DDbk1(d,m,n)*(Id(m,l) + (dt/2)*Dk1(m,l))*...
                        (Id(n,p) + (dt/2)*Dk1(n,p));
                end
            end
        end
    end
end

Dbk2 = Db_eval(x+dt*(k2/2),parameters);
DDbk2 = DDb_eval(x+dt*(k2/2),parameters);
DDk3 = zeros(D,D,D);
Id = eye(D);
for d = 1:D
    for p = 1:D
        for l = 1:D
            for m = 1:D
                DDk3(d,p,l) = DDk3(d,p,l) + (dt/2)*Dbk2(d,m)*DDk2(m,p,l);
                for n = 1:D
                    DDk3(d,p,l) = DDk3(d,p,l) + DDbk2(d,m,n)*(Id(m,l) + (dt/2)*Dk2(m,l))*...
                        (Id(n,p) + (dt/2)*Dk2(n,p));
                end
            end
        end
    end
end

Dbk3 = Db_eval(x+dt*k3,parameters);
DDbk3 = DDb_eval(x+dt*k3,parameters);
DDk4 = zeros(D,D,D);
Id = eye(D);
for d = 1:D
    for p = 1:D
        for l = 1:D
            for m = 1:D
                DDk4(d,p,l) = DDk4(d,p,l) + dt*Dbk3(d,m)*DDk3(m,p,l);
                for n = 1:D
                    DDk4(d,p,l) = DDk4(d,p,l) + DDbk3(d,m,n)*(Id(m,l) + dt*Dk3(m,l))*...
                        (Id(n,p) + dt*Dk3(n,p));
                end
            end
        end
    end
end

DDf = DDk1/6 + DDk2/3 + DDk3/3 + DDk4/6;

end