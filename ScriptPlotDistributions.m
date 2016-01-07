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

DEBUG = 1;
PROGRESS_BAR = 1;
if(exist('DISABLE_PROGRESS_BAR'))
    if DISABLE_PROGRESS_BAR
        PROGRESS_BAR = 0;
    end
end

RUN_TEST_CASES = 0;
TEST_CASE = 7; % indicate the test case to be tested.

VERIFY_RESULTS = 0;

SAVE   = 0;
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
GEN_FILE = 'Fig6_1.mat';
extraName = 'Fig6_1';
MC_samples    = 10000; % Number of MC   Samples.
MAX_ITER      = 50;
RENEWABLE_TYPE = 'wind';
% RENEWABLE_TYPE = 'solar';
W_RATIO = 0.35;
beta = 1;
distName = 'real';
% distName = 'normal';
% distName = 'tlocationscale';
% distName = 'generalized extreme value';
% distName = 'uniform';
P_STD = 0.35;
L_STD = P_STD;
W_STD = P_STD;
legendStr = {'real','normal','GEV'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation parameters
distName = 'real';
LoadParameters % load default parameters
figure(1)
[probabilities,xf] = ksdensity(w_rSamples(1,:));
plot(xf,probabilities);
hold on;
figure(2)
[probabilities,xf] = ksdensity(p_rSamples(1,:));
plot(xf,probabilities);
hold on;
figure(3)
[probabilities,xf] = ksdensity(L_rSamples(1,:));
plot(xf,probabilities);
hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
distName = 'normal';
LoadParameters % load default parameters
figure(1)
[probabilities,xf] = ksdensity(w_rSamples(1,:));
plot(xf,probabilities);
hold on;
figure(2)
[probabilities,xf] = ksdensity(p_rSamples(1,:));
plot(xf,probabilities);
hold on;
figure(3)
[probabilities,xf] = ksdensity(L_rSamples(1,:));
plot(xf,probabilities);
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
distName = 'generalized extreme value';
LoadParameters % load default parameters
figure(1)
[probabilities,xf] = ksdensity(w_rSamples(1,:));
plot(xf,probabilities);
hold on;
legend(legendStr);
figure(2)
[probabilities,xf] = ksdensity(p_rSamples(1,:));
plot(xf,probabilities);
hold on;
legend(legendStr);
figure(3)
[probabilities,xf] = ksdensity(L_rSamples(1,:));
plot(xf,probabilities);
hold on;
legend(legendStr);