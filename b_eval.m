function [b] = b_eval(x,parameters)

mu = parameters.mu;
nu = parameters.nu;
gamma = parameters.gamma;

b = [       mu*x(1) - x(3)*x(2);
           -nu*x(2) + x(3)*x(1);
     gamma -   x(3) + x(1)*x(2)];

end
