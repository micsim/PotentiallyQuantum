n = 1000;
u = (1:n)';
v = zeros(n,1);
v(floor(0.250*n):floor(0.280*n)) = 5e-4;
%v(350:450) = 0.004;

L = 1e-5;

%[e,E] = get_hamiltonian_eigenvectors(v, L);

%K = dft_eigenvectors(E);
% Fourier transform each eigenvector to get the eigenvectors in k space.

k = 5e7 * pi;
%d_orig = create_gauss_distr(n, L, 0.1, 0.03, k);
%d = fit_distribution(E, d_orig, L);

%plot_slider(e, E, d, 1e-8*0.001*v/(norm(v,inf) * norm(d,inf)^2), 1e-9, K);

simulation = EVSimulation(v, 1e-9, 0.1, 0.03, L, k);
simulation.plot = true;