function [f] = f_eval(x,parameters)

if(strcmp(parameters.solver,'EM'))
    f = b_eval(x,parameters);
elseif(strcmp(parameters.solver,'RK4'))
    f = RK4_eval(x,parameters);
end

end
