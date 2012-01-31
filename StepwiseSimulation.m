classdef StepwiseSimulation < Simulation
    %STEPWISESIMULATION Simulation with discrete timestepts and an
    %timedependent potential.
    
    properties (SetAccess = protected)
        psi = [];
    end
    
    properties (Dependent = true, SetAccess = protected)
        n_t;
        % Number of points the time is discretised in.
    end
    
    properties (Dependent = true)
        dt;
        % Time between two steps.
    end
    
    methods
        function o = StepwiseSimulation(varargin)
            o = o@Simulation(varargin{:});
        end

        function set.dt(o, value)
            o.t_max = value*n_t;
            o.redistribute();
            % TODO can that be deleted?
        end
        function out = get.dt(o)
            out = o.t_max / o.n_t;
        end
                
        function out = get.n_t(o)
            out = size(o.V,2) + 1;
        end
        
        function internal_recompute(o)
            o.internal_recompute_kinetic_hamiltonian();
        end
        
        function internal_redistribute(o)
            if isempty(o.T)
                o.internal_recompute();
            end
            
            distribution = o.internal_redistribute_momentum();
            
            joule_to_eV = 6.24150974e18; % [eV/J]
            h_bar = 1.054571506e-34; % [J*s]
            h_bar = h_bar * joule_to_eV; % covert to [eV*s]

            assert(size(distribution, 1) == o.n, 'd must have the same dimension as the potential in one point in time.');
            assert(size(distribution, 2) == 1, 'd must be a column vector!');

            o.psi = zeros(o.n, o.n_t);
            o.psi(:,1) = distribution;

            dt_local  = o.dt;
            n_local   = o.n;
            n_t_local = o.n_t;
            % Local copies for efficiency.
            for i = 1:n_t_local-1
                H = diag(o.V(:,i)) + o.T;
                H2 = 1i/h_bar * H * dt_local / 2;

                o.psi(:, i+1) = ((eye(o.n) + H2)\(eye(n_local) - H2)) * o.psi(:,i);
                % DOES NOT WORK: TODO why?
                %tmp = (eye(n_local) + 2*H2) * o.psi(:,i);
                %o.psi(:, i+1) = tmp/norm(tmp);
            end
        end
        
        function internal_replot(o)
            if isempty(o.T) || isempty(o.psi)
                o.internal_redistribute();
            end
            
            value = floor(get(o.slider, 'Value')*(o.n_t - 1)) + 1;
            
            set(0,'CurrentFigure',o.x_plot);
            o.plot_distr(o.psi(:,value), o.V_factor * o.V(:,max(1,value-1)));

            set(0,'CurrentFigure',o.k_plot);
            o.plot_distr(fft(o.psi(:,value)));
        end
                
        function t_max_was_reset(o)
            o.redistribute();
        end
    end
end

