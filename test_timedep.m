n = 300;
L = 1e-5;
u = (1:n)';
v = zeros(n,1);
v(floor(0.250*n):floor(0.280*n)) = 3.5e-5;

n_t = 1000;
dt = 3e-12;

%d_orig = create_gauss_distr(n, L, 0.1, 0.03, 1e7 * pi);

v_mat = zeros(n, n_t);
for i = 1:n_t
    v_mat(:,i) = v;
end

%psi = cn_evolution(v_mat, d_orig, L, dt);

%plot_timedep_slider(psi, n_t, dt, 4e0*0.001*v_mat/(norm(v,inf) * norm(d_orig,inf)^2));

%[e,E] = get_hamiltonian_eigenvectors(v, L);
%d = fit_distribution(E, d_orig, L);
%plot_slider(e, E, d, 4e0*0.001*v/(norm(v,inf) * norm(d_orig,inf)^2), n_t*dt);

simulation = StepwiseSimulation(v_mat, n_t*dt, 0.1, 0.03, L, 1e7 * pi);
simulation.plot = true;