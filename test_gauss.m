u = (1:500)';
v = zeros(500,1);
%v(350:450) = 0.004;

[e,E] = get_hamiltonian_eigenvectors(v, 3);

d = fit_distribution(E, normpdf(u, 200, 50));

evolve_eigenvectors(e,E,0.01);

plot_distr(E, d, v);
