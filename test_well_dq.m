n = 1000;
u = (1:n)';

v = 1e12*ones(n,1);
v(1/4 * n:3/4 * n) = 0;

L = 1e-9; % [m]
[e,E] = get_hamiltonian_eigenvectors(v, L,@(H) eig(full(H)),@get_kinetic_hamiltonian_dq);

%groundstate energy for an infinite well:
h_bar = 1.054571506e-34; % [J*s]
joule_to_eV = 6.24150974e18; % [eV/J]
electron_mass = 9.1093829140e-31; % [kg]

genergy = (1:10)'.^2 * joule_to_eV * h_bar^2 * pi^2 / (2*electron_mass * (L/2)^2);
both = [genergy, e(1:10)];

fprintf('groundstate energy for an infinite well:\nby formula | the calculated eigenvalue:\n');
disp(both);

fprintf('relative error:\n');
disp((genergy(1:10) - e(1:10))./genergy(1:10));

plot_amplitudes(E,e,5, v);

