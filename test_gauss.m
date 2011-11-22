n = 500;
u = (1:n)';
v = zeros(500,1);
%v(350:450) = 0.004;

[e,E] = get_hamiltonian_eigenvectors(v, 3);


k = 1;

d_orig = normpdf(u, 200, 50);
for j = 1:n
    d_orig(j) = exp(1i*k*j) * d_orig(j);
end

d = fit_distribution(E, d_orig);

%evolve_eigenvectors(e,E,0.01);

plot_distr(E, d, v);
