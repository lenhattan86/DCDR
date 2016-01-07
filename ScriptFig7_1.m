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
%GENERATE_DATA = 0;
GEN_FILE = 'Fig7_1.mat';
extraName = 'Fig7_1';
MC_samples    = 2000; % Number of MC   Samples.
MAX_ITER      = 50;
RENEWABLE_TYPE = 'wind';
% RENEWABLE_TYPE = 'solar';
W_RATIO = 0.35;
beta = 1;
P_STD = 0.15;
L_STD = P_STD;
W_STD = P_STD;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation parameters
LoadParameters % load default parameters

if MC_samples > 400
    RUN_OFFICAL_SIM = 1;
else
    RUN_OFFICAL_SIM = 0;
end

if RUN_OFFICAL_SIM==1
    DEBUG = 0;
    SKIP_PLOTTING  = 1;
    RUN_TEST_CASES = 0;
end

%% Start simulation
tic
%% Plot the objective function with 2-dim q_l
% q_l_max = 1000;
% d_q = 100;
% [ bestCost_plot3d, q_l_best_plot3d, X, Y, Z ] = plot3d_EPLT(d_q, q_l_max, p_l, w_rSamples, p_rSamples, L_rSamples);
% figure(2);
% surfc(X,Y,Z);

%% Other methods
% NoLT+NoGLB-RT
q_l_NoLT_NoGLB = zeros(N,1);
[lt_cost_NoLT_NoGLB, rt_energy_cost_NoLT_NoGLB, queuing_delay_cost_NoLT_NoGLB, network_delay_cost_NoLT_NoGLB] ...
    = fixedLT_noGLB(p_l, q_l_NoLT_NoGLB, w_rSamples, p_rSamples, L_rSamples, mu, M, beta, pi_ij);
%%
% fixedLT+NoGLB-RT
q_l_fixedLTNoGLB = fixedLT*LOAD_SCALE*M;
[exp_cost_fixedLTNoGLB, rt_energy_cost_fixedLTNoGLB, queuing_delay_cost_fixedLTNoGLB, network_delay_cost_fixedLTNoGLB] ...
    = fixedLT_noGLB(p_l, q_l_fixedLTNoGLB, w_rSamples, p_rSamples, L_rSamples, mu, M, beta, pi_ij);


%%
% NoLT+GLB-RT
q_l_NoLT = zeros(N,1);
[exp_cost_NoLT, rt_energy_cost_NoLT, queuing_delay_cost_NoLT, network_delay_cost_NoLT] ...
    = expected_cost( p_l, q_l_NoLT, w_rSamples, p_rSamples, L_rSamples, mu ,M , beta ,pi_ij);
%%
% fixedLT+GLB-RT
q_l_fixedLT = fixedLT*LOAD_SCALE*M;
[exp_cost, rt_energy_cost_fixedLT, queuing_delay_cost_fixedLT, network_delay_cost_fixedLT] ...
    = expected_cost( p_l, q_l_fixedLT, w_rSamples, p_rSamples, L_rSamples, mu ,M , beta ,pi_ij);
exp_cost_fixedLT = p_l'*q_l_fixedLT; 
%%
% Prediction mean based method (PA)
w_r_pred = w_mean;
p_r_pred = p_r_mean;
L_r_pred = ones(J,1)*L_mean;
[q_l_pred, m, lambda, a, b, c, d, isSuccessful] = ...
    pred_procurement(p_l, w_r_pred, p_r_pred, L_r_pred , mu, M, beta ,pi_ij);
[exp_cost_pred, rt_energy_cost_pred, queuing_delay_cost_pred, network_delay_cost_pred] = ...
    expected_cost( p_l, q_l_pred, w_rSamples, p_rSamples, L_rSamples, mu , M , beta , pi_ij);
%%
% Deterministic optimization method (BA)
[ exp_offline_cost, offline_worst_cost, offline_best_cost ] = ...
    expected_offline_cost( p_l, q_l, w_rSamples, p_rSamples, L_rSamples ,mu, M ,beta, pi_ij );

%% Stochastic Gradient Estimate based Algorithm (SGEA)
close all;
q_l_0 = q_l_pred;
% q_l_SGEA = q_l_pred; grad_norm = 0;
[expect_lt_cost, q_l_SGEA, grad_norm, gradMatrix, objValues ] = ...
    long_term_procurement(q_l_0, p_l, w_rSamples, p_rSamples, L_rSamples, mu ,M , beta ,pi_ij, STEP_SIZE_SETTINGS);

[exp_cost_SGEA, rt_energy_cost_SGEA, queuing_delay_cost_SGEA, network_delay_cost_SGEA] = ...
    expected_cost( p_l, q_l_SGEA, w_rSamples, p_rSamples, L_rSamples, mu , M , beta , pi_ij);


%% V. Verify output
if RUN_TEST_CASES
    VERIFY_RESULTS = 1;    
    %RunTestCase
end

%% VI. Plot Results

% NoLT_NoGLB = [lt_cost_NoLT_NoGLB, rt_energy_cost_NoLT_NoGLB, queuing_delay_cost_NoLT_NoGLB, network_delay_cost_NoLT_NoGLB];
NoLT_NoGLB_ignoreQDelay = [lt_cost_NoLT_NoGLB, rt_energy_cost_NoLT_NoGLB, queuing_delay_cost_NoLT_NoGLB, network_delay_cost_NoLT_NoGLB];
% NoLT_NoGLB = [0, 0, 0, 0];

fixedLT_NoGLB = [exp_cost_fixedLTNoGLB, rt_energy_cost_fixedLTNoGLB, ...
        queuing_delay_cost_fixedLTNoGLB, network_delay_cost_fixedLTNoGLB];

NoLT = [(exp_cost_NoLT - rt_energy_cost_NoLT - queuing_delay_cost_NoLT - network_delay_cost_NoLT) ...
            rt_energy_cost_NoLT queuing_delay_cost_NoLT network_delay_cost_NoLT];    
        
fixedLT = [exp_cost_fixedLT, rt_energy_cost_fixedLT, queuing_delay_cost_fixedLT, network_delay_cost_fixedLT]        

lt_energy_cost_pred = exp_cost_pred - rt_energy_cost_pred - queuing_delay_cost_pred - network_delay_cost_pred;
pred_alg = [lt_energy_cost_pred rt_energy_cost_pred queuing_delay_cost_pred network_delay_cost_pred];

expect_SGEA = [(exp_cost_SGEA - rt_energy_cost_SGEA - queuing_delay_cost_SGEA - network_delay_cost_SGEA) ...
    rt_energy_cost_SGEA queuing_delay_cost_SGEA, network_delay_cost_SGEA];   

if ~SKIP_PLOTTING     
    figure;
    y = [NoLT_NoGLB_ignoreQDelay; fixedLT_NoGLB; NoLT; fixedLT; pred_alg; expect_SGEA; exp_offline_cost];
    bar_chart = bar( y, 'stacked');
    AX=legend(bar_chart, { 'long-term energy cost', 'real-time energy cost', ...
            'queueing delay cost', 'network delay cost'});
    set(gca,'XTickLabel',{'NoLT NoGLB', 'fixedLT NoGLB','NoLT', 'fixed LT', 'PA', 'SGEA', 'BA'})
    titleStr = ['\beta=' num2str(beta) 'grad-norm=' num2str(grad_norm) ...
            ' E[p^r_{std}]=' num2str(mean(p_r_std)) ' L_{std}=' num2str(L_std)  ...
            ' E[L^r]=' num2str(L_mean) ' w_{std}=' num2str(w_std)  ' E[w^r]=' num2str(w_mean) ];
%     title(titleStr);
    grid on;
end

%% V. Saving Results
performance_gap = (exp_cost_pred - sum(exp_offline_cost));
if ~exist('extraName')
    extraName='';
end
extraName = [extraName '_gap' num2str(round(performance_gap/sum(exp_offline_cost) *100))];

if or(SAVE,RUN_OFFICAL_SIM)
    formatOut = 'mmdd_hhMM';
    SIM_TIMESTAMP_STOP = datetime('now'); 
    savefile = ['results/comparisons/' extraName '_' datestr(SIM_TIMESTAMP_STOP,formatOut) '.mat']
    save(savefile);
else
    save(['results/comparisons/' extraName '.mat']);
end
toc