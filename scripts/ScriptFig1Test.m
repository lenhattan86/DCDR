%% Reset & clear memory
if ~exist('DONT_RESET')
    clear; clear global variable; close all; 
    clc;
end
%% classpath
addpath('functions');
%% Load parameter

% Debuging parameters
global DEBUG PROGRESS_BAR;
global MAX_ITER GRADIENT_ESTIMATE
DEBUG = 0;
PROGRESS_BAR = 1;
if(exist('DISABLE_PROGRESS_BAR'))
    if DISABLE_PROGRESS_BAR
        PROGRESS_BAR = 0;
    end
end

RUN_TEST_CASES = 0;
TEST_CASE = 7; % indicate the test case to be tested.
VERIFY_RESULTS = 0;

SAVE   = 1;
SKIP_PLOTTING = 0;
if ~exist('GENERATE_DATA')
    GENERATE_DATA = 1;
end


% FD for low dimensional problem, SP for high dimensional problem
GRADIENT_ESTIMATE = 'CVX' % CVX
% GRADIENT_ESTIMATE = 'FD' % FD
% GRADIENT_ESTIMATE = 'SP' % SP 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE_DATA = 0; 
GEN_FILE = 'scriptConvergence.mat';
extraName = ['Fig1_maxiter_' int2str(MAX_ITER)]
MC_samples    = 100; % Number of MC   Samples.
MAX_ITER      = 30;
RENEWABLE_TYPE = 'wind';
% RENEWABLE_TYPE = 'solar';
W_RATIO = 0.95;
beta = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation parameters
LoadParameters
alp = 0.602; [ A, a] = getStepSizeSetting(1, 30, 2000, 100, alp) %( Start, End, stepStart, stepEnd, alpha )
STEP_SIZE_SETTINGS =   [A   a  alp   0.1 0.101]

if MC_samples > 300
    RUN_OFFICAL_SIM = 1;
else
    RUN_OFFICAL_SIM = 0;
end

if RUN_OFFICAL_SIM==1
    SKIP_PLOTTING  = 0;
    RUN_TEST_CASES = 0;
end

%% Start simulation
tic
w_r_pred = w_mean;
p_r_pred = p_r_mean;
L_r_pred = ones(J,1)*L_mean;
[q_l_pred, m, lambda, a, b, c, d, isSuccessful] = ...
    pred_procurement(p_l, w_r_pred, p_r_pred, L_r_pred , mu, M, beta ,pi_ij);
[exp_cost_pred, rt_energy_cost_pred, queuing_delay_cost_pred, network_delay_cost_pred] = ...
    expected_cost( p_l, q_l_pred, w_rSamples, p_rSamples, L_rSamples, mu , M , beta , pi_ij);

%% Stochastic Gradient Estimate based Algorithm (SGEA)
q_l_pred
q_l_0 = q_l_pred;
[expect_lt_cost, q_l_SGEA, grad_norm, gradMatrix, objValues] = ...
    long_term_procurement(q_l_0, p_l, w_rSamples, p_rSamples, L_rSamples, mu ,M , beta ,pi_ij, STEP_SIZE_SETTINGS);
[exp_cost_SGEA, rt_energy_cost_SGEA, queuing_delay_cost_SGEA, network_delay_cost_SGEA] = ...
    expected_cost( p_l, q_l_SGEA, w_rSamples, p_rSamples, L_rSamples, mu , M , beta , pi_ij);
q_l_SGEA
%% VI. Plot Results
expect_SGEA = [(exp_cost_SGEA - rt_energy_cost_SGEA - queuing_delay_cost_SGEA - network_delay_cost_SGEA) ...
    rt_energy_cost_SGEA queuing_delay_cost_SGEA, network_delay_cost_SGEA];   

if ~SKIP_PLOTTING     
    figure;
    xIteration = 1:MAX_ITER;
    % plot figure of gradients.
    plot(xIteration, gradMatrix');
    grid on;
    % plot figure of objective values
    figure
    plot(xIteration, objValues);    
    grid on;
end

%% V. Saving Results
if ~exist('extraName')
    extraName='';
end

if or(SAVE,RUN_OFFICAL_SIM)
    formatOut = 'mmdd_hhMM';
    SIM_TIMESTAMP_STOP = datetime('now'); 
    savefile = ['results/' extraName '_' datestr(SIM_TIMESTAMP_STOP,formatOut) '.mat']
    save(savefile);
else
    save(['results/' extraName '.mat']);
end
toc