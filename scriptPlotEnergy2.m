energy_nLTnGLB = [5.1108  0.6035]*10^5;
%% fLTnGLB
energy_fLTnGLB = [5.1108 0.6035]*10^5;
%% nLT
energy_nLT =  [2.2658 0.6035]*10^5;
%% fLT
energy_fLT =  [2.2954  0.6035]*10^5;
%% PA
energyPA = [2.4353 0.5998]*10^5;
%% SGEA
energySGEA = [2.5332 0.5970]*10^5;
%% Offline
energyOffline= [2.4800 0.6035]*10^5;
%%
brown = [0.635294139385223 0.0784313753247261 0.184313729405403];
%path = 'C:\Users\NhatTan\Dropbox\Research-Tan\SharedManuscript\figs\';
% figPath = 'C:\Users\NhatTan\Dropbox\Research-Tan\SharedManuscript\figs\';
figPath = 'results/figures/';

fontAxis = 14;
fontTitle = 14;
fontLegend = 14;
LineWidth = 2;
FontSize = 14;

figure1 = figure;
axes1 = axes('Parent',figure1);
% y = [NoLT; fixedLT; pred_alg; expect_SGEA; exp_offline_cost]; lengendStr = {'NoLT', 'fixed long-term', 'PA', 'SGA', 'OA'};
y = [energy_nLTnGLB; energy_fLTnGLB; energy_nLT;  energy_fLT; energyPA; energySGEA; energyOffline];  
xLabels = {'nLTnGLB', 'fLTnGLB', 'nLT', 'fLT', 'PA', 'SGA', 'OA'};
legendStr = {'non-renewable energy', 'renewable energy'};
bar_chart = bar(y, 'stacked','barwidth',0.5);
set(bar_chart,{'FaceColor'},{'g';brown}); 
% axis([0.5 7.5 0 10.5e6]);
set(gca,'fontsize',fontAxis);
legend(legendStr,'Location','northeast','FontSize',fontLegend,'Orientation','vertical');
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 4.0 3.0]);
set(gca,'XTickLabel',xLabels,'FontSize',fontAxis);
xlabel('algorithms','FontSize',fontAxis);
ylabel('cost','FontSize',fontAxis);
print ('-depsc', [figPath 'energy_consumption.eps']);