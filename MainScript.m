%% Reset & clear memory
clear all; clc; close all; clear global variable;
%% classpath
addpath('functions');
tic

%% Load parameter
% Debuging parameters
RUN_OFFICAL_SIM = 0;
global DEBUG PROGRESS_BAR;
DEBUG = 0;
PROGRESS_BAR = 1;
RUN_TEST_CASES = 1;
TEST_CASE = 1; % indicate the test case to be tested.
VERIFY_RESULTS = 0;

SAVE   = 0;
SKIP_PLOTTING = 0;
GENERATE_DATA = 1;

if RUN_OFFICAL_SIM==1
    DEBUG = 0;
    RUN_TEST_CASES = 0;
end

% Simulation parameters
LoadParameters % load default parameters

%% Start simulation

%% Long term procurement

dqmin = norm(partial_q*ones(N,1));
[ expect_cost, q_l_best, grad_norm] = long_term_procurement(q_l_0, p_l, w_rSamples, p_rSamples, L_rSamples)

%% Other methods

% real_time only (RT_only)
q_l_RT_only = zeros(N,1);
exp_cost_RT_only = expected_cost( p_l, q_l_RT_only, w_rSamples, p_rSamples, L_rSamples)

%% V. Verify output
if RUN_TEST_CASES
    VERIFY_RESULTS = 1;    
    RunTestCase
end

%% V. Plot Results
if ~SKIP_PLOTTING
end

%% VI. Saving Results
if or(SAVE,RUN_OFFICAL_SIM)
    SIM_TIMESTAMP_STOP = clock;    	
    savefile = ['results/' strrep(strcat(num2str(SIM_TIMESTAMP_STOP(1:3))),' ','') ...
            '_' strrep(strcat(num2str(SIM_TIMESTAMP_STOP(4:5))),' ','')  '_simResult.mat']
    save(savefile);
else
    save(strcat('results/simResult.mat'));    
end