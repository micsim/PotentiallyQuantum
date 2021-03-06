n = 1000;
u = (1:n)';
v = zeros(n,1);
%v(floor(0.250*n):floor(0.280*n)) = 5e-4;
%v(350:450) = 0.004;

L = 1e-5;
v(floor(n/2)) = n/L;

%[e,E] = get_hamiltonian_eigenvectors(v, L);


k = 20e6 * pi;
%d_orig = create_gauss_distr(n, L, 0.1, 0.03, k);
%d = fit_distribution(E, d_orig, L);

%plot_slider(e, E, d, 1e-8*0.001*v/(norm(v,inf) * norm(d,inf)^2), 6e-9);
sim = EVSimulation(v, 6e-9, 0.1, 0.03, L, k);
sim.plot = true;