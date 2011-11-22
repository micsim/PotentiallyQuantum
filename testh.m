v = ones(1000,1);
for(i=1:600)
v(i,1)=0;
end

[e,E] = get_hamiltonian_eigenvectors(v, 100);
plot_probabilities(E,e,5);
