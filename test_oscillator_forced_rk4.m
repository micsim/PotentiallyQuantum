n = 600;
u = (1:n)';
L = 2e-5;
n_t = 4000;

omega = 1e9; % [rad/s]
dt = 1*(2*pi/(omega*n_t));%3e-11;


electron_mass = 9.1093829140e-31; % [kg]
joule_to_eV = 6.24150974e18; % [eV/J]

omega_pot = 0.8*omega;  %angular frequency of the harmonically oscillating potential
ampl_pot = 0.03;        %amplitude of the potential oscillation, relative to L

v_mat = zeros(n, n_t);
for i = 1:n_t
    v_mat(:,i) = (((u/n - 0.5+ampl_pot*sin(omega_pot*i*dt)).*L).^2).*(0.5*electron_mass*omega^2*joule_to_eV); % potential of harmonic oscillator;
end

%[e,E,H] = get_hamiltonian_eigenvectors(v, L, @(H) eigs(H,100,0));
% To calculate all eigenvectors delete the last argument.


%K = zeros(n,1);
% Fourier transformation disabled for testing.


k = 0e6 * pi;

d_orig = sqrt(normpdf(u, (0.5)*n, 0.03*n));
for j = 1:n
    pos = j/n * L;
    d_orig(j) = exp(1i*k*pos) * d_orig(j);
end

psi = rk4_evolution(v_mat, d_orig, L, dt);

plot_timedep_slider(psi, n_t, dt, 4e0*0.001*v_mat/(norm(v_mat(:,1),inf) * norm(d_orig,inf)^2));

exp_values = zeros(n_t,1);

for l = 1:n_t
    exp_values(l,1) = u'*(abs(psi(:,l)).^2)/(norm(psi(:,l)))^2;
end
plot(exp_values);
