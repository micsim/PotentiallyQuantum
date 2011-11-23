n = 500;
u = (1:n)';
v = zeros(n,1);
v(350:450) = 0.275;
%v(350:450) = 0.004;

[e,E] = get_hamiltonian_eigenvectors(v, 3);


k = 1;

d_orig = normpdf(u, 150, 50);
for j = 1:n
    d_orig(j) = exp(1i*k*j) * d_orig(j);
end

d = fit_distribution(E, d_orig);

plot_slider(e, E, d, 0.1*v, 1000);
