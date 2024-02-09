%%% takes in gridsearch results and outputs a variety of 
%%% visualizations and metrics

%% uniform degradation dots
% metric = 'scaling_L2';
% load('gridsearch_results/uniform_degradation_gridsearch.mat', 'T', 'ranges');
% param_indices = 1 : 7;
% labels = {'K_H', '\alpha_S', 'D_{HS, rel}', '\kappa_{HS, rel}', ...
%     '\gamma_S', 'n_H', 'D_S'};


%% self-enhanced dots; low production
% metric = 'scaling_domains';
% load('gridsearch_results/self_enhanced_low_gridsearch.mat', 'T', 'ranges');
% param_indices = 3 : 4;
% labels = {'', '', 'D_{HS, rel}', '\kappa_{HS, rel}'};


%% self-enhanced dots; high production
% metric = 'scaling_domains';
% load('gridsearch_results/self_enhanced_high_gridsearch.mat', 'T', 'ranges');
% param_indices = 3 : 4;
% labels = {'', '', 'D_{HS, rel}', '\kappa_{HS, rel}'};


%% self-enhanced dots; high production, MSE error metric
% metric = 'scaling_L2';
% load('gridsearch_results/self_enhanced_high_gridsearch.mat', 'T', 'ranges');
% param_indices = 3 : 4;
% labels = {'', '', 'D_{HS, rel}', '\kappa_{HS, rel}'};


%% flux with phi dots
% metric = 'scaling_domains';
% load('gridsearch_results/flux_scaling_gridsearch.mat', 'T', 'ranges');
% param_indices = 3 : 4;
% labels = {'', '', 'D_{HS, rel}', '\kappa_{HS, rel}'};


%% scube changes degradation
% metric = 'scaling_L2';
% load('gridsearch_results/degradation_gridsearch.mat', 'T', 'ranges');
% param_indices = 3 : 4;
% labels = {'', '', '\gamma_{HS, rel}', '\kappa_{HS, rel}'};


%% heterzygous knockout
metric = 'scaling_domains';
load('gridsearch_results/2_5_het_HH_KO.mat', 'T', 'ranges');
param_indices = 3 : 4;
labels = {'', '', 'D_{HS, rel}', '\kappa_{HS, rel}'};


%%
params = fieldnames(ranges);
output_tables = struct();

for i = param_indices
    figure(i)
    
    % param, scaling pairs
    scaling = T.(metric);
    param_vals = T.(params{i});
    
    xs = scaling;
    
    % replace "param_vals" with an index for each value of the parameter
    param_indices = ones(size(xs));
    uniques = unique(param_vals);
    for j = 1 : numel(uniques)
        param_indices(param_vals == uniques(j)) = j;
    end
    
    dy = 10;
    jitter = 4;
    
    ys = dy * param_indices + (rand(size(param_indices)) * 2 - 1) * jitter;

    scatter(xs, ys, 0.3)
    
    yticks(unique(param_indices) * dy);
    yticklabels(unique(param_vals));
    
    xlabel('scaling');
    ax = gca;
    ax.XScale = 'log';
    xlabel('Scaling (fractional shift in domain boundaries)');
    ylabel(labels{i});
    % xline(0.1);
end


