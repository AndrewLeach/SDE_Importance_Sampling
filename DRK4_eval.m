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