function T = get_kinetic_hamiltonian(n, L, mass)
    assert(mod(n,2) == 0, 'n must be even!');

    if nargin < 3
        mass = 1;
    end

    h_bar = 1.054571506e-34; % [J*s]
    joule_to_eV = 6.24150974e18; % [eV/J]
    electron_mass = 9.1093829140e-31; % [kg]

    % factor is   h^2 ("h-bar") K^2  / 2m
    %           = h^2 ("h-bar") pi^2 / 2m(dx^2)
    %           = h^2 ("h-bar") pi^2 (n-1)^2 / (2m L^2)
    factor = joule_to_eV * h_bar^2 * pi^2 * (n-1)^2 / (2*mass*electron_mass * L^2); % [eV]

    % Build the hamiltonian:
    %  First the kinetic energy:
    %    The diagonal:
    T = factor * 1/3 * (1 + 2/(n^2)) * eye(n);
    % An diagonal matrix of all the same values (that value times the identity).

    %    The other elements:
    for i = 1:n
        for j = 1:n
            if i ~= j
                T(i,j) = factor * 2/(n^2) * (-1)^(j-i) / (sin(pi*(j-i)/n))^2;
            end
        end
    end
end
