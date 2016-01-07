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
MC_samples    = 2000; % Number of MC   Samples.
MAX_ITER      = 150;
RENEWABLE_TYPE = 'wind';
% RENEWABLE_TYPE = 'solar';
W_RATIO = 0.35;
beta = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation parameters
LoadParameters
alp = 0.902; [ A, a] = getStepSizeSetting(1, 50, 5000, 100, alp) %( Start, End, stepStart, stepEnd, alpha )
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
% NoLT+GLB-RT
q_l_NoLT = zeros(N,1);
[exp_cost_NoLT, rt_energy_cost_NoLT, queuing_delay_cost_NoLT, network_delay_cost_NoLT] ...
    = expected_cost( p_l, q_l_NoLT, w_rSamples, p_rSamples, L_rSamples, mu ,M , beta ,pi_ij);
%% Stochastic Gradient Estimate based Algorithm (SGEA)
% q_l_0 = q_l_pred;
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