n = 1000;
v = zeros(n,1);
v(floor(0.250*n):floor(0.280*n)) = 5e-4;

L = 1e-5;

k = 1e8 * pi;

sim = EVSimulation(v, 2e-9, 0.1, 0.03, L, k);
sim.plot = true;