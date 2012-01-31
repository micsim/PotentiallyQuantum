classdef EVSimulation < Simulation
    %EVSIMULATION Simulation with an invariant set of eigenvectors and an
    %timeindependet potential.
    
    properties
        eig_function = @(H) eig(full(H));
    end
    properties (SetAccess = protected)
        eigenvalues = [];
        eigenvectors = [];
        k_eigenvectors = [];
        H = [];
        % The hamiltonian.
        
        d_fitted = [];
    end
    
    methods
        function o = EVSimulation(varargin)
            o = o@Simulation(varargin{:});
        end
        
        % Calculate T, H, eigenvalues, eigenvectors:
        function internal_recompute(o)
            assert(size(o.V,2) == 1, 'Potential V must be a column vector!');

            % The hamiltonian itself:
            o.H = o.internal_recompute_get_hamiltonian(o.V);

            % Find eigenvalues and eigenvectors:
            [eve,eva] = o.eig_function(o.H);

            % normalize:
            eve = sqrt(o.n/o.L)*eve;

            % The eigenvalues as a vector:
            eva = diag(eva);
            % Vector of the values on the diagonal.

            % Sort eigenvalues and store order:
            [o.eigenvalues, o.eigenvectors] = sort_eigenvectors(eva, eve);
            % TODO Optimize: Not needed when using eigs as eig_function!
            
            % calculate k eigenvectors:
            o.k_eigenvectors = fft(o.eigenvectors);
        end
        
        function internal_redistribute(o)
            if isempty(o.eigenvectors) || isempty(o.eigenvalues)...
                || isempty(o.k_eigenvectors)
                o.internal_recompute();
            end

            distribution = o.internal_redistribute_momentum();
            
            o.d_fitted = fit_distribution(o.eigenvectors, distribution/norm(distribution), o.L);
        end
                        
        function internal_replot(o)
            if isempty(o.eigenvalues) || isempty(o.eigenvectors)...
                    || isempty(o.k_eigenvectors) || isempty(o.d_fitted)
                o.internal_redistribute();
            end
            
            value = get(o.slider, 'Value')*o.t_max;
            E_new = evolve_eigenvectors(o.eigenvalues, o.eigenvectors,...
                                        value);

            set(0,'CurrentFigure',o.x_plot);
            o.plot_distr(E_new * o.d_fitted, o.V_factor * o.V);

            k_new = evolve_eigenvectors(o.eigenvalues, o.k_eigenvectors, value);

            set(0,'CurrentFigure',o.k_plot);
            o.plot_distr(k_new * o.d_fitted);
        end
    end
end

