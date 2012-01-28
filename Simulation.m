classdef Simulation < handle
    %SIMULATION Parent class for EVSimulation, StepwiseSimulation.
    %       Simulate a particle in 1D with periodic boundary conditions
    
    properties
        V = [];
        % Potential timeindependent -> n x 1       Matrix
        %           timedependent   -> n x n_t - 1 Matrix
        L = 1e-5;
        % Length of simulated space.
        k = 0;
        % k to start with.
        d = [];
        % The probability distribution of the particle at the starting
        % time.
        d_mu = 0.5;
        d_omega = 0.1;
        % If d_mu and d_omega are not empty, compute d from them in
        % reevolve as a normal distribution.
        
        plot;
        % Wether the class is mainting a plot at the moment.
    end
    
    properties (Dependent = true, SetAccess = private)
        n;
        % Number of points the space is discretised in.
    end
    
    methods
        function set.V(o, value)
            o.set_property_recompute('V', value);
        end
        function set.L(o, value)
            o.set_property_recompute('L', value);
        end
        function set.k(o, value)
            o.set_property_reevolve('k', value);
        end
        
        function set.d(o, value)
            if(o.d ~= value)
                if isempty(value)
                    if isempty(o.d_mu)
                        o.d_mu = 0.5;
                    end
                    if isempty(o.d_omega)
                        o.d_omega = 0.1;
                    end
                else
                    o.d_mu = [];
                    o.d_omega = [];
                end
                % Eigther d is nonempty or d_mu and d_omega are.
                
                o.set_property_reevolve('d', value);
            end
        end
        function set.d_mu(o, value)
            o.d_mu = value;
            
            if ~isempty(o.d_mu) && ~isempty(o.d_omega)
                o.reevolve();
            end
        end
        function set.d_omega(o, value)
            o.d_omega = value;
            
            if ~isempty(o.d_mu) && ~isempty(o.d_omega)
                o.reevolve();
            end
        end
        
        function set.plot(o, value)
            o.plot = value;
            
            if value
                o.replot();
            end
        end
        
        function out = get.n(o)
            out = size(o.V,1);
        end
        
    end
    
    methods (Access = private)
        % Redo all computations.
        function set_property_recompute(o, name, value)
            o.(name) = value;
            
            o.recompute();
        end
        
        % Redo the evolution including changes from the starting k
        function set_property_reevolve(o, name, value)
            o.(name) = value;
            
            o.reevolve();
        end
        
        function recompute(o)
            % TODO
            o.replot();
        end
        function reevolve(o)
            % TODO
            o.replot();
        end
        
        function replot(o)
            % TODO
            if o.plot
                o.replot();
            end
        end
    end
end

