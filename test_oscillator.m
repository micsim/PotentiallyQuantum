n = 2000;
u = (1:n)';
L = 1e-5;

omega = 5e8; % [rad/s]

electron_mass = 9.1093829140e-31; % [kg]
joule_to_eV = 6.24150974e18; % [eV/J]

v = (((u/n - 0.5).*L).^2).*(0.5*electron_mass*omega^2*joule_to_eV); % potential of harmonic oscillator
%v(350:450) = 0.004;



[e,E,H] = get_hamiltonian_eigenvectors(v, L, @(H) eigs(H,100,0));
% To calculate all eigenvectors delete the last argument.


k = 5e6 * pi;

d_orig = sqrt(normpdf(u, 0.5*n, 0.03*n));
for j = 1:n
    pos = j/n * L;
    d_orig(j) = exp(1i*k*pos) * d_orig(j);
end

d = fit_distribution(E, d_orig, L);

%calculate energy of particle
totalenergy = (E*d)'* H * (E*d);
fprintf('total energy of particle: %d\n',totalenergy);

%disp('This should be 1:');
%disp(sum(abs(E*d).^2));
plot_slider(e, E, d, 1e-8*0.001*v/(norm(v,inf) * norm(d,inf)^2), 1e-7);
%plot_probabilities(E,e,1);

%lowest eigenenergies of harmonic oscillator:
h_bar = 1.054571506e-34; % [J*s]

genergy = ((0:99)'+0.5).*(joule_to_eV * h_bar * omega);
both = [genergy, e(1:100)];

fprintf('lowest eigenenergies for a harmonic oscillator:\nby formula | the calculated eigenvalue:\n');
disp(both);

fprintf('relative error:\n');
disp((genergy(1:100) - e(1:100))./genergy(1:100));
