function [f] = RK4_eval(x,parameters)

dt = parameters.dt;

k1 = b_eval(x,parameters);
k2 = b_eval(x + dt*(k1/2),parameters);
k3 = b_eval(x + dt*(k2/2),parameters);
k4 = b_eval(x + dt*(k3),parameters);

f = k1/6 + k2/3 + k3/3 + k4/6;

end