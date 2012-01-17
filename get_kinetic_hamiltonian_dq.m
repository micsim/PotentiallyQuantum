function T = get_kinetic_hamiltonian_dq(n, L, mass)
    h_bar = 1.054571506e-34; % [J*s]
    joule_to_eV = 6.24150974e18; % [eV/J]
    electron_mass = 9.1093829140e-31; % [kg]

    % factor is   h^2 ("h-bar")/ 2m(dx^2)
    %           = h^2 ("h-bar") (n-1)^2 / (2m L^2)
    % where dx = L/(n-1) is the step-width of position-space.
    factor = joule_to_eV * h_bar^2 * (n-1)^2 / (2*mass*electron_mass * L^2); % [eV]

    % Build the hamiltonian:
    %  First the kinetic energy:
    %    The diagonal:
    T = factor * 2* eye(n);
    % An diagonal matrix of all the same values (that value times the identity).

    %    The other elements:
    for i = 2:n
        T(i,i-1) = -factor;
        T(i-1,i) = -factor;
    end
    
    % Add border values for periodic boundary conditions
    T(1,n) = -factor;
    T(n,1) = -factor;
end
