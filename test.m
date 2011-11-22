v = sin((1:100)' * 0.001).^2;

[e,E] = get_hamiltonian_eigenvectors(v, 3);
plot_probabilities(E,e,5);
