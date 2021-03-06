n = 1000;
v = zeros(n,1);
v(floor(0.250*n):floor(0.280*n)) = 5e-4;

L = 1e-5;
k = 5e7 * pi;

sim = EVSimulation(v, 1e-9, 0.1, 0.03, L, k);
sim.plot = true;