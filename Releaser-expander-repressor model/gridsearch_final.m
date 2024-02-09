% conducts gridsearch analysis; generates the scaling metric for each
% gradient on a grid of the given parameters
    
ranges = struct();
ranges.K_H = [1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1];
ranges.a_S = [1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1];
ranges.D_HS_rel = [1, 3, 10, 30];
ranges.Phi_HS_rel = [1, 3, 10, 30];
ranges.c_S = [0, 0.001, 0.01, 0.1];
ranges.n_H = [1, 3, 10];
ranges.D_S = [10, 100, 1000];

% turn parameter ranges into grid
param_names = fieldnames(ranges);
param_ranges = cell(length(param_names), 1);
for i=1:length(param_names)
    param_ranges{i} = ranges.(param_names{i});
end 

grid_cell = cell(length(param_names), 1);
[grid_cell{:}] = ndgrid(param_ranges{:});

grid = [];
for i = 1 : length(grid_cell)
    for j = 1 : numel(grid_cell{i})
        grid(j).(param_names{i}) = grid_cell{i}(j); %#ok<SAGROW>
    end
end

% set up some things to display
ppm = utils.ParforProgressbar(numel(grid) *9);
disp('Running models...')
tic
ticBytes(gcp)


%% Model 1: First order degradation; constant flux 

lengths = [200, 300];
phis = [1, 1];
load('variables/uniform_degradation', 'paramsModel', 'paramsSys');

paramsSys.max_t = 1440;
paramsSys.t_min = 1440;

handle = @models.HS;

% loop through parameter grid
base_paramsModel = paramsModel;
output = zeros(numel(grid), 3);
parfor i = 1 : numel(grid)
    % set up parameters and model
    this_paramsModel = utils.merge_structs(base_paramsModel, grid(i));
    model = handle(this_paramsModel, paramsSys);
    
    output(i, :) = get_scaling(this_paramsModel, paramsSys, handle, lengths, phis)
    
    ppm.increment();
end

% put output into a table
params_T = struct2table(grid);
output_T = array2table(output);
output_T.Properties.VariableNames = {'scaling_L2', 'scaling_integrated', 'scaling_domains'}; 

T = [params_T, output_T];
writetable(T, 'uniform_degradation_gridsearch.xlsx');
save('uniform_degradation_gridsearch.mat', 'T', 'ranges');


%% Model 2: Second order degradation, low production regime

lengths = [200, 300];
phis = [1, 1] * 0.1;
load('variables/self_enhanced', 'paramsModel', 'paramsSys');

paramsSys.max_t = 1440;
paramsSys.t_min = 1440;

handle = @models.HS_self_enhanced;

% loop through parameter grid
base_paramsModel = paramsModel;
output = zeros(numel(grid), 3);
parfor i = 1 : numel(grid)
    % set up parameters and model
    this_paramsModel = utils.merge_structs(base_paramsModel, grid(i));
    model = handle(this_paramsModel, paramsSys);
    
    output(i, :) = get_scaling(this_paramsModel, paramsSys, handle, lengths, phis)
    
    ppm.increment();
end

% put output into a table
params_T = struct2table(grid);
output_T = array2table(output);
output_T.Properties.VariableNames = {'scaling_L2', 'scaling_integrated', 'scaling_domains'}; 

T = [params_T, output_T];
writetable(T, 'self_enhanced_low_gridsearch.xlsx');
save('self_enhanced_low_gridsearch.mat', 'T', 'ranges');

%% Model 3: Second order degradation; high production regime 

lengths = [200, 300];
phis = [1, 1] * 30;
load('variables/self_enhanced', 'paramsModel', 'paramsSys');

paramsSys.max_t = 1440;
paramsSys.t_min = 1440;

handle = @models.HS_self_enhanced;

% loop through parameter grid
base_paramsModel = paramsModel;
output = zeros(numel(grid), 3);
parfor i = 1 : numel(grid)
    % set up parameters and model
    this_paramsModel = utils.merge_structs(base_paramsModel, grid(i));
    model = handle(this_paramsModel, paramsSys);
    
    output(i, :) = get_scaling(this_paramsModel, paramsSys, handle, lengths, phis)
    
    ppm.increment();
end

% put output into a table
params_T = struct2table(grid);
output_T = array2table(output);
output_T.Properties.VariableNames = {'scaling_L2', 'scaling_integrated', 'scaling_domains'}; 

T = [params_T, output_T];
writetable(T, 'self_enhanced_high_gridsearch.xlsx');
save('self_enhanced_high_gridsearch.mat', 'T', 'ranges');


%% Model 4: Flux scales with L

lengths = [200, 300];
phis = lengths / 300;
load('variables/uniform_degradation', 'paramsModel', 'paramsSys');

paramsSys.max_t = 1440;
paramsSys.t_min = 1440;

handle = @models.HS;

% loop through parameter grid
base_paramsModel = paramsModel;
output = zeros(numel(grid), 3);
parfor i = 1 : numel(grid)
    % set up parameters and model
    this_paramsModel = utils.merge_structs(base_paramsModel, grid(i));
    model = handle(this_paramsModel, paramsSys);
    
    output(i, :) = get_scaling(this_paramsModel, paramsSys, handle, lengths, phis)
    
    ppm.increment();
end

% put output into a table
params_T = struct2table(grid);
output_T = array2table(output);
output_T.Properties.VariableNames = {'scaling_L2', 'scaling_integrated', 'scaling_domains'}; 

T = [params_T, output_T];
writetable(T, 'flux_scaling_gridsearch.xlsx');
save('flux_scaling_gridsearch.mat', 'T', 'ranges');


%% Model 5: First order degradation; heterozygous HH knockout

lengths = [300, 300];
phis = [0.5, 1];
load('variables/uniform_degradation', 'paramsModel', 'paramsSys');

paramsSys.max_t = 1440;
paramsSys.t_min = 1440;

handle = @models.HS;

% loop through parameter grid
base_paramsModel = paramsModel;
output = zeros(numel(grid), 3);
parfor i = 1 : numel(grid)
    % set up parameters and model
    this_paramsModel = utils.merge_structs(base_paramsModel, grid(i));
    model = handle(this_paramsModel, paramsSys);
    
    output(i, :) = get_scaling(this_paramsModel, paramsSys, handle, lengths, phis)
    
    ppm.increment();
end

% put output into a table
params_T = struct2table(grid);
output_T = array2table(output);
output_T.Properties.VariableNames = {'scaling_L2', 'scaling_integrated', 'scaling_domains'}; 

T = [params_T, output_T];
writetable(T, 'het_KO_gridsearch.xlsx');
save('het_KO_gridsearch.mat', 'T', 'ranges');

%% Model 6: Modify degradation

lengths = [200, 300];
phis = [1, 1];
load('variables/uniform_degradation', 'paramsModel', 'paramsSys');
paramsModel.D_HS_rel = 1;
paramsModel.c_HS_rel = 1;

paramsSys.max_t = 1440;
paramsSys.t_min = 1440;

handle = @models.HS;

% change the ranges...
ranges = struct();
ranges.K_H = [1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1];
ranges.a_S = [1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1];
ranges.c_HS_rel = [0.001, 0.03, 0.1, 0.3, 1, 3, 10, 30];
ranges.Phi_HS_rel = [0.001, 0.03, 0.1, 0.3, 1, 3, 10, 30];
ranges.c_S = [0, 0.001, 0.01, 0.1];
ranges.n_H = [1, 3, 10];
ranges.D_S = [10, 100, 1000];

% turn parameter ranges into grid
param_names = fieldnames(ranges);
param_ranges = cell(length(param_names), 1);
for i=1:length(param_names)
    param_ranges{i} = ranges.(param_names{i});
end 

grid_cell = cell(length(param_names), 1);
[grid_cell{:}] = ndgrid(param_ranges{:});

grid = [];
for i = 1 : length(grid_cell)
    for j = 1 : numel(grid_cell{i})
        grid(j).(param_names{i}) = grid_cell{i}(j); %#ok<SAGROW>
    end
end

% loop through parameter grid
base_paramsModel = paramsModel;
output = zeros(numel(grid), 3);
parfor i = 1 : numel(grid)
    % set up parameters and model
    this_paramsModel = utils.merge_structs(base_paramsModel, grid(i));
    model = handle(this_paramsModel, paramsSys);
    
    output(i, :) = get_scaling(this_paramsModel, paramsSys, handle, lengths, phis)
    
    ppm.increment();
end

% put output into a table
params_T = struct2table(grid);
output_T = array2table(output);
output_T.Properties.VariableNames = {'scaling_L2', 'scaling_integrated', 'scaling_domains'}; 

T = [params_T, output_T];
writetable(T, 'degradation_gridsearch.xlsx');
save('degradation_gridsearch.mat', 'T', 'ranges');

%%

disp('Finished!')
toc
tocBytes(gcp)
delete(ppm);


function eta = get_scaling(paramsModel, paramsSys, handle, lengths, phis)
    % runs model for the two given lengths and computes the scaling metric
    % between the two steady state gradients
    
    paramsSys.L = lengths(1);
    paramsModel.a_H = phis(1);
    model = handle(paramsModel, paramsSys);
    [~, y_short] = run_model(paramsSys, model);
    A_short = squeeze(y_short(:, 1, :));
    
    paramsSys.L = lengths(2);
    paramsModel.a_H = phis(2);
    model = handle(paramsModel, paramsSys);
    [~, y_long] = run_model(paramsSys, model);
    A_long = squeeze(y_long(:, 1, :));

    eta1 = metrics.L2_area(A_long, A_short);
    eta2 = metrics.integrated_area(A_long, A_short, 120, paramsSys.t_step, paramsSys.max_t);
    eta3 = metrics.domains(A_long, A_short, 3);
    
    eta = [eta1, eta2, eta3];
end
