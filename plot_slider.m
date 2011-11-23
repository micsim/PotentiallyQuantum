function plot_slider(e, E, d, v, t_max)
    if nargin < 5
        t_max = 2*size(v,1);
    end

    fig = figure();
    slider = uicontrol(fig, 'Style', 'Slider', 'Units', 'normalized',...
        'Position', [0,-0.95,1,1]);
    warning('off','MATLAB:addlistener:invalidEventName');
    addlistener(slider,'ActionEvent',@(~,~) slider_moved());

    plot_distr(E, d, v);

    function slider_moved()
        value = get(slider, 'Value');
        E_new = evolve_eigenvectors(e, E, value*t_max);

        plot_distr(E_new, d, v);
    end
end
