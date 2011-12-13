% plots the time evolution of the position expectation
% value of the particle represented by the distribution.
% e = eigenvalues
% E = eigenvectors
% d = distribution
% max_time = max of plotted time in seconds
% num_steps = number of steps plotted between 0 and max_time
function plot_pos_exp_val(e, E, d, max_time, num_steps)

    n = size(E,1);
    assert(size(d,2) == 1, 'distribution must be a column vector.');   
    
    timestep = max_time/num_steps;
    
    %time values for x-axis
    time = (0:num_steps)';
    time = time.*timestep;
    
    %pos values for y-axis
    exp_val_pos = zeros(num_steps,1);
    for t=0:num_steps 
        
        %first version calculated the index of the position with
        %the highest probability value:
        %[a,ind] = max(abs(E*d));
        %pos(t+1) = ind;
        
        %second version calculates the expectation value:
        x = (1:n)';
        wavefct = E*d;
        exp_val_pos(t+1) = sum(x.*(abs(wavefct).^2))/sum(abs(wavefct).^2);
        E = evolve_eigenvectors(e, E, timestep);
    end
    
    %plot expectation value of position over time:
    plot(time,exp_val_pos);
    
end
