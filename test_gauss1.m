n = 1000;
v = zeros(n,1);
v(floor(0.250*n):floor(0.280*n)) = 5e-4; % eV

L = 1e-5;
k = 5e6 * pi;
d_mu = 0.1;
d_omega = 0.03;

t_max = 4e-9;

sim = EVSimulation(v, t_max, d_mu, d_omega, L, k);
sim.plot = true;