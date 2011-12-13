n = 1000;
u = (1:n)';
v = zeros(n,1);
v(floor(0.250*n):floor(0.280*n)) = 8e-6;
%v(350:450) = 0.004;

L = 1e-5;

[e,E] = get_hamiltonian_eigenvectors(v, L);


k = 5e6 * pi;

d_orig = sqrt(normpdf(u, 0.1*n, 0.03*n));
for j = 1:n
    pos = j/n * L;
    d_orig(j) = exp(1i*k*pos) * d_orig(j);
end

d = fit_distribution(E, d_orig, L);

%disp('This should be 1:');
%disp(sum(abs(E*d).^2));

plot_slider(e, E, d, 1e-8*0.001*v/(norm(v,inf) * norm(d,inf)^2), 5e-9);
%plot_probabilities(E,e,1);
