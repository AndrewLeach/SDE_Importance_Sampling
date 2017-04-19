function [] = Test_Samplers(batch_post,sampler_post,solver_post,eps_post,M_post,T_post,dt_post)

load('Initialization_T_10');

batch = batch_post;
sampler = sampler_post;
solver = solver_post;
eps = eps_post;
M = M_post;
T = T_post;
dt = dt_post;

parameters.mu = 0.119;
mu = parameters.mu;
parameters.nu = 0.1;
nu = parameters.nu;
parameters.gamma = 0.9;
gamma = parameters.gamma;
parameters.D = 3;
D = parameters.D;
parameters.dt = dt;

parameters.sigma = 1e00*eye(3);
sigma = parameters.sigma;
parameters.M = M;

N = T/parameters.dt;
parameters.N = N;


start_idx = 1;
end_idx = N;


fp1 = [sqrt(nu+gamma*sqrt(nu/mu));-sqrt(mu+gamma*sqrt(mu/nu));-sqrt(mu*nu)];
fp2 = [-sqrt(nu+gamma*sqrt(nu/mu));sqrt(mu+gamma*sqrt(mu/nu));-sqrt(mu*nu)];

%X0 = fp1;
%XT = fp2;
parameters.X0 = X0;
parameters.XT = XT;

%parameters.X0 = X0;
parameters.g = 'normal';
g = parameters.g;
parameters.g_sigma = 1e00*eye(3);

g_sigma = parameters.g_sigma;
parameters.g_mu = XT;
g_mu = parameters.g_mu;

Compare_Samplers_Pole_Reversal(batch,sampler,solver,phi_init,X0,sigma,eps,M,N,dt,D,mu,nu,gamma,g,g_mu,g_sigma)
end