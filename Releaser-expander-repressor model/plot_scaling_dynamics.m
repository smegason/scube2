% runs a model for two different lengths of the patterning
% field, then plots ligand, scube, and scube expression dynamics

% plot relative or absolute x-axis
scale_axis = 1;
t_end = 1440;

load('variables/uniform_degradation');
paramsModel.a_S = 1e-02;
paramsModel.c_S = 3e-03;

%%% generate model
handle = @models.HS;

paramsSys.t_min = t_end;
paramsSys.max_t = t_end;
paramsSys.t_step = 2;

lengths = [200, 300];
phis = [1, 1];

indices_to_plot = {1, 2, 'p'};
output_names = {'HH.gif', 'Scube.gif', 'Scube_expression.gif'};
y_names = {'Scaled ligand concentration', 'Scaled scube concentration', 'Scube expression'};


for i = 1 : 3
    make_anim(indices_to_plot{i}, paramsModel, paramsSys, handle, phis, lengths, scale_axis, y_names{i}, output_names{i});
end
    

function make_anim(index_to_plot, paramsModel, paramsSys, handle, phis, lengths, scale_axis, y_name, output_name)
    disp('Running model for short field...')
    [A_short, t_short] = get_dynamics(index_to_plot, paramsModel, paramsSys, handle, phis(1), lengths(1));
    disp('Running model for long field...')
    [A_long, t_long] = get_dynamics(index_to_plot, paramsModel, paramsSys, handle, phis(2), lengths(2));

    % as time, pick the longest 
    if length(t_long) > length(t_short)
        t = t_long;
    else
        t = t_short;
    end

    scaling = metrics.domains(A_long, A_short);
    disp(['Scaling score: ' num2str(scaling)]);

    % disp(['short end: ' num2str(A_short(end, end))]);
    % disp(['long end: ' num2str(A_long(end, end))]);

    disp('Generating animation...')
    plotfunc = @(i) plot_data(paramsSys, lengths, A_short, A_long, scale_axis, y_name, t, i);
    utils.animate(length(t), plotfunc, output_name);
end


function [output, t] = get_dynamics(index, paramsModel, paramsSys, handle, phi, L)
    paramsSys.L = L;
    paramsModel.a_H = phi;
    model = handle(paramsModel, paramsSys);
    [t, y] = run_model(paramsSys, model);
    
    A = zeros(length(t), paramsSys.N);
    for i = 1 : length(t)
        A(i, :) = model.get_perception(squeeze(y(i, :, :)));
    end

    if strcmp(index(1), 'p')
        output = A;
    else
        output = squeeze(y(:, index, :));
    end
end


function plot_data(paramsSys, lengths, A_short, A_long, scale_axis, y_name, t, i);
    % plots the ith timestepfor all data sets in the cell array 
    % data at the ith timstep
   
    % set up axes
    y_min = min(min(A_short, [], 'all'), min(A_long, [], 'all'));
    y_max = max(max(A_short, [], 'all'), max(A_long, [], 'all'));
    a = 0.05;  % how much extra space to show
    
    y_min = max(y_min, -7);  % so it doesn't break when on a log scale
    % y_max = 2;
    
    % y_axis = [y_min - a*(y_max - y_min), y_max + a*(y_max - y_min)];
    y_axis = [y_min, y_max + a*(y_max - y_min)];
    x_axis = lengths(2) / paramsSys.N * (1 : paramsSys.N);
    
    clf;
    hold on;
    
    gradient_short = squeeze(A_short(i, :));
    gradient_long = squeeze(A_long(i, :));
    
    if scale_axis
        plot(x_axis / lengths(2), gradient_short);
        plot(x_axis / lengths(2), gradient_long);
        xlabel('Relative position x/L');
    else
        x_axis_short = 1 : ceil(lengths(1) / lengths(2)) * paramsSys.N;
        interped = interp1((1 : paramsSys.N) * (lengths(1) / paramsSys.N), ...
            gradient_short, x_axis_short);
        plot(x_axis_short, interped);
        
        plot(x_axis, gradient_long);
        xlabel('Position (um)');
    end
    
    ylim(y_axis);
    ylabel(y_name);
    legend(['L = ' num2str(lengths(1))], ['L = ' num2str(lengths(2))]);
    title(['$$t = ' num2str(t(i)) '$$'], 'interpreter', 'latex')
    
    hold off;
end
