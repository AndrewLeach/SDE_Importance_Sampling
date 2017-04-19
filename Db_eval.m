function [Db] = Db_eval(x,parameters)

mu = parameters.mu;
nu = parameters.nu;

Db = [   mu, -x(3),  -x(2);
       x(3),   -nu,   x(1);
       x(2),  x(1),     -1];

end
