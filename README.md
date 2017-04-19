# SDE_Importance_Sampling

Test_Samplers will run and save a batch (of the chaotic pole reversal ex.) as a .mat file

## input:
batch number (unique identifier for naming, not the number of batches, takes a string)  
sampler ('NONE','LM','SLM','DLM','SDLM')  
solver (use 'RK4')  
eps (the epsilons to test, ex., 10.^(-7:-3) )  
batch size (number of samples, ex., 30 or 120 )  
T ( units of time, use 10 )  
dt ( size of time step, use 1/10)  

## output:
Q  (relative variance for each epsilon)  
X  (samples drawn)  
logW (log weights)  
logO (log observations)

## example:
Test_Samplers('1','LM','RK4',10.^(-7:-2),30,10,1/10)
