n = 1000;
u = (1:n)';
v = zeros(n,1);
v(floor(0.250*n):floor(0.280*n)) = 5e-4;
%v(350:450) = 0.004;

L = 1e-5;

%[e,E] = get_hamiltonian_eigenvectors(v, L);

k = 5e6 * pi;
d_mu = 0.1;
d_omega = 0.03;

%d_orig = create_gauss_distr(n, L, d_mu, d_omega, k);
%d = fit_distribution(E, d_orig, L);

t_max = 4e-9;

%plot_slider(e, E, d, 1e-8*0.001*v/(norm(v,inf) * norm(d,inf)^2), t_max);

sim = EVSimulation(v, t_max, d_mu, d_omega, L, k);
sim.plot = true;
H = sim.H;