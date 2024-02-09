% generates a plot of the scaling factor on the relative change 
% in the diffusion constant and the release rate

%% uniform degradation
% D_HSs = 10 .^ (0 : 0.1 : 2);
% Phi_HSs = 10 .^ (0 : 0.1 : 2);
% lengths = [200, 300];
% phis = [1, 1];
% load('variables\uniform_degradation.mat');
% handle = @models.HS;
% metric = 'profiles';


%% self enhanced; low production
% D_HSs = 10 .^ (0 : 0.1 : 2);
% Phi_HSs = 10 .^ (0 : 0.1 : 2);
% lengths = [200, 300];
% phis = [1, 1] * 0.1;
% load('variables\self_enhanced.mat');
% handle = @models.HS_self_enhanced;
% metric = 'domains';


%% self enhanced; high production
% D_HSs = 10 .^ (0 : 0.1 : 2);
% Phi_HSs = 10 .^ (0 : 0.1 : 2);
% lengths = [200, 300];
% phis = [1, 1] * 30;
% load('variables\self_enhanced.mat');
% handle = @models.HS_self_enhanced;
% metric = 'domains';

%% self enhanced; high production, MSE scaling metric
D_HSs = 10 .^ (0 : 0.1 : 2);
Phi_HSs = 10 .^ (0 : 0.1 : 2);
lengths = [200, 300];
phis = [1, 1] * 30;
load('variables\self_enhanced.mat');
handle = @models.HS_self_enhanced;
metric = 'L2';


%% first order degradation, phi scales with L
% D_HSs = 10 .^ (0 : 0.1 : 2);
% Phi_HSs = 10 .^ (0 : 0.1 : 2);
% lengths = [200, 300];
% phis = lengths / 300;
% load('variables\uniform_degradation.mat');
% handle = @models.HS;
% metric = 'domains';


%% heterozygous knockout
% D_HSs = 10 .^ (0 : 0.1 : 2);
% Phi_HSs = 10 .^ (0 : 0.1 : 2);
% lengths = [300, 300];
% phis = [0.5, 1];
% load('variables\self_enhanced.mat');
% handle = @models.HS_self_enhanced;
% metric = 'domains';


%%
paramsSys.t_min = 1440;
paramsSys.max_t = 1440;

base_paramsModel = paramsModel;

num_iterations = numel(Phi_HSs) * numel(D_HSs);
data = zeros(1, num_iterations);

% for parfor loop later on
coords = zeros(numel(Phi_HSs), numel(D_HSs), 2);
for i = 1 : numel(Phi_HSs)
    for j = 1 : numel(D_HSs)
        coords(i, j, :) = [i, j];
    end
end
coords = reshape(coords, num_iterations, 2);

% set up progress bar
ppm = utils.ParforProgressbar(num_iterations);
disp('Running models...');
tic

parfor i = 1 : num_iterations
    % generate parameter set
    x = coords(i, 1);
    y = coords(i, 2);
    paramsModel = base_paramsModel;
    paramsModel.Phi_HS_rel = Phi_HSs(x);
    paramsModel.D_HS_rel = D_HSs(y);
        
    % run models to get scaling parameter
    data(i) = get_scaling(paramsModel, paramsSys, handle, phis, lengths, metric);
    
    ppm.increment();
end

data = reshape(data, numel(Phi_HSs), numel(D_HSs));

delete(ppm);
disp('Finished!');
toc

% generate heatmap of eta on D/A
h = imagesc([1,100], [1,100], data);
xlabel('Relative change in diffusion constant D_{HS} / D_H');   
ylabel('Relative change in release rate \kappa_{HS} / \kappa_H');
cb = colorbar();
ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';
ax.YDir = 'normal';
ax.ColorScale = 'log';
cb.Colormap = turbo;
ax.Colormap = turbo;
ylabel(cb, 'Scaling (fractional shift in domain boundaries)')

disp(' ');
disp('Best parameter set: ');
best_scaling = min(min(data));
[i, j] = find(data == best_scaling);
disp(['D_HS = ' num2str(D_HSs(j))]);
disp(['Phi_HS = ' num2str(Phi_HSs(i))]);
disp(['Scaling metric: ' num2str(best_scaling)]);
% caxis([0.03, 0.16]);


function eta = get_scaling(paramsModel, paramsSys, handle, phis, lengths, metric)
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

    if strcmp('domains', metric)
        eta = metrics.domains(A_short, A_long);
    else
        eta = metrics.L2_area(A_short, A_long);
    end
end
        
