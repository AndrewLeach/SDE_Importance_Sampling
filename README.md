This package is SDE_Importance_Sampling.  This code
accompanies the paper "Symmetrized importance samplers for
stochastic differential equations" by Andrew Leach, Kevin K
Lin, and Matthias Morzfeld.

Copyright (C) 2017 by Andrew Leach <imaleach@gmail.com>.

This program is free software; you can redistribute
it and/or modify it under the terms of the GNU
General Public License as published by the Free
Software Foundation; either version 2 of the
License, or (at your option) any later version.

This program is distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General
Public License along with this program; if not, write
to the Free Software Foundation, Inc., 51 Franklin
Street, Fifth Floor, Boston, MA 02110-1301 USA.


# DynamicLinearMap

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
