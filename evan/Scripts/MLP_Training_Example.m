%% MLP Example
load('80train_20test.mat');
%loads data: training_in test_in desired_train desired_test 
desired_test = desired_test';

ncomp=2500;
weight_ini = 0.01;
num_hid_nodes = 25;
num_out_nodes = length(unique(desired_train)); % # of classes
num_iter=4000;
save_times=unique(ceil(logspace(0,log10(num_iter-1e-8),500)));
stop_point=-1; %-1 = don't stop (controlled by annealing break logic instead)

%% Learning Rate (n) Vector
% Generalized logistic curve
% in log(iterations) time-scale
vvo=1;QQo=.5;AAo=.4e-4;KKo=2e-4;BBo=3;MMo=0;
x=linspace(-1,1,length(save_times));
nvec=AAo+(KKo-AAo)./(1+QQo*exp(-BBo*(x-MMo))).^(1/vvo);

% Summary parameters (inputs) for generation internally
%nlogistic=[vvo,QQo,AAo,KKo,BBo,MMo];

% Use nvec as input for precomputed curve
nvec=interp1(save_times,nvec,1:num_iter);

%% Annealing parameters    
K=[0 0 0.3 1e-5 .01 .95 1];
%% Run
[weights1, out_weights1, train_out, test_out,AUC,MSE,Annealing_Log] = ...
    MLP_PCA_variable_learning_MSE_differential_annealing_v6(...
        training_in, test_in, desired_test, desired_train,... 
        weight_ini, num_hid_nodes, num_out_nodes, ...
        nvec, num_iter,save_times,stop_point,K,'AUC',0);