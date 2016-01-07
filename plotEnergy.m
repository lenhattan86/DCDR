% energy plot
clear; close all; clc;
%% classpath
addpath('functions');
%% load data
opfficialRunNo = 6;

%path = 'C:\Users\NhatTan\Dropbox\Research-Tan\SharedManuscript\figs\';
% figPath = 'C:\Users\NhatTan\Dropbox\Research-Tan\SharedManuscript\figs\';
figPath = 'results/figures/';

fontAxis = 14;
fontTitle = 14;
fontLegend = 14;
LineWidth = 2;
FontSize = 14;

brown = [0.635294139385223 0.0784313753247261 0.184313729405403];


fileFig_2 = 'results/comparisons/Fig2_M__gap3_1128_1946.mat';

PAColor = 'x-.b';
PAString = 'PA';
SGAColor = 'o-r';
SGAString = 'SGA';
BAColor = '>-.k';
BAString = 'PA';

TotalCostColor = '>-.k';
TotalCostString = 'Total cost';
EnergyCostColor = 'o-r';
EnergyCostString = 'energy cost';
LTCostColor = 'x-.b';
LTCostString = 'long-term energy cost';



disp('Figure 2b: energy consumption');
load(fileFig_2)


%% hold(axes1,'on');
w_renew = sum(mean(w_rSamples,2));
%%nLTnGLB
energy_nLTnGLB = [sum(M)-w_renew w_renew];
%% fLTnGLB
energy_fLTnGLB = [sum(M)-w_renew w_renew];
%% nLT
energy_nLT =  computeEnergyConsumption(p_l, q_l_NoLT, w_rSamples, p_rSamples, L_rSamples, mu ,M , beta ,pi_ij);
%% fLT
energy_fLT =  computeEnergyConsumption(p_l, q_l_fixedLT, w_rSamples, p_rSamples, L_rSamples, mu ,M , beta ,pi_ij);
%% PA
energyPA = computeEnergyConsumption(p_l, q_l_pred, w_rSamples, p_rSamples, L_rSamples, mu ,M , beta ,pi_ij);
%% SGEA
energySGEA = computeEnergyConsumption(p_l, q_l_SGEA, w_rSamples, p_rSamples, L_rSamples, mu ,M , beta ,pi_ij);
%% Offline
[ exp_offline_cost, offline_worst_cost, offline_best_cost, energyOffline ] = ...
    expected_offline_cost( p_l, q_l, w_rSamples, p_rSamples, L_rSamples ,mu, M ,beta, pi_ij );
%%
figure1 = figure;
axes1 = axes('Parent',figure1);
% y = [NoLT; fixedLT; pred_alg; expect_SGEA; exp_offline_cost]; lengendStr = {'NoLT', 'fixed long-term', 'PA', 'SGA', 'OA'};
y = [energy_nLTnGLB; energy_fLTnGLB; energy_nLT;  energy_fLT; energyPA; energySGEA; energyOffline];  
xLabels = {'nLTnGLB', 'fLTnGLB', 'nLT', 'fLT', 'PA', 'SGA', 'OA'};
legendStr = {'non-renewable energy cost', 'renewable energy'};
bar_chart = bar(y, 'stacked','barwidth',0.5);
set(bar_chart,{'FaceColor'},{'g';brown}); 
axis([0.5 7.5 0 10.5e6]);
set(gca,'fontsize',fontAxis);
legend(legendStr,'Location','northeast','FontSize',fontLegend/2,'Orientation','vertical');
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 4.0 3.0]);
set(gca,'XTickLabel',xLabels,'FontSize',fontAxis/2);
xlabel('algorithms','FontSize',fontAxis/2);
ylabel('cost','FontSize',fontAxis/2);
print ('-depsc', [figPath 'energy_consumption.eps']);