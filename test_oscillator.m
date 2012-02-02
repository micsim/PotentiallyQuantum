n = 2000;
u = (1:n)';
L = 1e-5;

omega = 5e8; % [rad/s]

electron_mass = 9.1093829140e-31; % [kg]
joule_to_eV = 6.24150974e18; % [eV/J]

v = (((u/n - 0.5).*L).^2).*(0.5*electron_mass*omega^2*joule_to_eV); % potential of harmonic oscillator
%v(350:450) = 0.004;



%[e,E,H] = get_hamiltonian_eigenvectors(v, L, @(H) eigs(H,100,0));
% To calculate all eigenvectors delete the last argument.


%K = dft_eigenvectors(E);
% Fourier transform each eigenvector to get the eigenvectors in k space.


k = 5e6 * pi;
%d_orig = create_gauss_distr(n, L, 0.5, 0.03, k);
%d = fit_distribution(E, d_orig, L);

sim = EVSimulation(v, 1e-7, 0.5, 0.03, L, k);
sim.eig_function = @(H) eigs(H,100,0);
sim.distribute();
e = sim.eigenvalues;
E = sim.eigenvectors;
K = sim.k_eigenvectors;
d = sim.d_fitted;
H = sim.H;

%calculate energy of particle
totalenergy = (E*d)'* H * (E*d);
fprintf('total energy of particle: %d eV\n',totalenergy);

%plot_slider(e, E, d, 1e-8*0.001*v/(norm(v,inf) * norm(d,inf)^2), 1e-7, K);

%lowest eigenenergies of harmonic oscillator, calculated analytically:
h_bar = 1.054571506e-34; % [J*s]

genergy = ((0:55)'+0.5).*(joule_to_eV * h_bar * omega);
both = [genergy, e(1:56)];

fprintf('lowest eigenenergies for a harmonic oscillator:\nby formula | the calculated eigenvalue:\n');
disp(both);

fprintf('relative error:\n');
disp((genergy(1:56) - e(1:56))./genergy(1:56));

sim.plot = true;
