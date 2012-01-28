classdef Simulation
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
        
        plot;
        % Wether the class is mainting a plot at the moment.
    end
    
    properties (Dependent = true, SetAccess = private)
        n;
        % Number of points the space is discretised in.
    end
    
    methods
        function set.plot(o, value)
            o.set_property_replot('plot', value);
        end
        
        function out = get.n(o)
            out = size(o.V,1);
        end
        
    end
    
    methods (Access = private)
        function set_property_replot(o, name, value)
            o.(name) = value;
            
            if o.plot
                o.replot();
            end
        end
        
        function replot(o)
            % TODO
        end
    end
end

