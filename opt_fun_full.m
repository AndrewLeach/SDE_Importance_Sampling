function [ f, Grad, Hessian ] = opt_fun_full(X0,X,parameters )

D = parameters.D;
N = length(X)/D;

r = r_eval(X0,X,parameters);
Dr = Dr_eval(X0,X,parameters);
f = r'*r;
Grad = Gradient_eval(Dr,r,X0,X,parameters);
Hessian = Hessian_eval(Dr,r,X0,X,parameters);

end

