%test case for a time dependent potential:
%Forced harmonic oscillation, no damping,
%drive frequency = 1.2 * eigenfrequency

n = 300;
u = (1:n)';
L = 1e-5;
n_t = 4000;
dt = 4*(2*pi/(omega*n_t));%3e-11;

omega = 1.2e9; % [rad/s]

electron_mass = 9.1093829140e-31; % [kg]
joule_to_eV = 6.24150974e18; % [eV/J]

omega_pot = 1*omega;  %angular frequency of the harmonically oscillating potential
ampl_pot = 0.03;        %amplitude of the potential oscillation, relative to L

v_mat = zeros(n, n_t);
for i = 1:n_t
    v_mat(:,i) = (((u/n - 0.5+ampl_pot*sin(omega_pot*i*dt)).*L).^2).*(0.5*electron_mass*omega^2*joule_to_eV); % potential of harmonic oscillator;
end


k = 0; %particle has no impact at t=0

d_orig = sqrt(normpdf(u, (0.5+ampl_pot)*n, 0.03*n)); 
for j = 1:n
    pos = j/n * L;
    d_orig(j) = exp(1i*k*pos) * d_orig(j);
end

psi = cn_evolution(v_mat, d_orig, L, dt);

plot_timedep_slider(psi, n_t, dt, 4e0*0.001*v_mat/(norm(v,inf) * norm(d_orig,inf)^2));

exp_values = zeros(n_t,1);

for l = 1:n_t
    exp_values(l,1) = u'*(abs(psi(:,l)).^2)/(norm(psi(:,l)))^2;
end
plot(exp_values);
