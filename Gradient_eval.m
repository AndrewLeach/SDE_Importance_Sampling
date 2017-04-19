function [Grad] = Gradient_eval(Dr,r,X0,X,parameters)

Grad = 2*(r')*Dr;

end