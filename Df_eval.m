function [Df] = Df_eval(x,parameters)

if(strcmp(parameters.solver,'EM'))
    Df = Db_eval(x,parameters);
elseif(strcmp(parameters.solver,'RK4'))
    Df = DRK4_eval(x,parameters);
end

end
