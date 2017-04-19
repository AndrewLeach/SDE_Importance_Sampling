function [DDf] = DDf_eval(x,parameters)

if(strcmp(parameters.solver,'EM'))
    DDf = DDb_eval(x,parameters);
elseif(strcmp(parameters.solver,'RK4'))
    DDf = DDRK4_eval(x,parameters);
end

end
