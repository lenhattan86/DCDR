clear; close all; clc;
%% classpath
addpath('functions');
%% load data
opfficialRunNo = 6;

%path = 'C:\Users\NhatTan\Dropbox\Research-Tan\SharedManuscript\figs\';
%figPath = 'C:\Users\NhatTan\Dropbox\Research-Tan\SharedManuscript\figs\';
figPath = 'results/figures/';

fontAxis = 14;
fontTitle = 14;
fontLegend = 14;
LineWidth = 2;
FontSize = 14;

brown = [0.635294139385223 0.0784313753247261 0.184313729405403];


switch opfficialRunNo
    case 6
        % Nov. 28
        disp('Results on Nov. 28');
        fileFig_1 = 'results/Fig1_maxiter_30_1128_2323.mat';

        fileFig_2 = 'results/comparisons/Fig2_M__gap3_1128_1946.mat';

        fileFig_5a = ['results/comparisons/Fig5a_5_gap2_1128_2011.mat  ' ...
                     ;'results/comparisons/Fig2_M__gap3_1128_1946.mat  ' ...
 
                     ;'results/comparisons/Fig5a_65_gap5_1128_2011.mat ' ...
                     ;'results/comparisons/Fig5a_95_gap6_1128_2131.mat ' ... %Fig5a_95_gap6_1125_1439 
                    ];

        fileFig_5b = ['results/comparisons/Fig5b_5_gap2_1128_0942.mat ' ...
                     ;'results/comparisons/Fig5b_35_gap3_1128_0934.mat' ...
                     ;'results/comparisons/Fig5b_65_gap3_1128_0921.mat' ...
                     ;'results/comparisons/Fig5b_95_gap4_1128_1031.mat' ... 
                    ]; 

        fileFig_6a = ['results/comparisons/Fig6_1_gap2_1128_1602.mat ' ...  % normal  
                     ;'results/comparisons/Fig6_5_gap1_1128_0719.mat ' ... % real
                     ;'results/comparisons/Fig6_3_gap1_1128_0708.mat ' ... % GEV
%                      ;'results/comparisons/Fig6_2_gap2_1128_0722.mat ' ... % tlocationscale
%                      ;'results/comparisons/Fig6_4_gap5_1128_0711.mat ' ... % uniform
                     
                    ];
        fileFig_7 = ['results/comparisons/Fig7_1_gap0_1128_1133.mat ' ...     
                     ;'results/comparisons/Fig6_5_gap1_1128_0719.mat ' ...
                     ;'results/comparisons/Fig7_2_gap3_1128_1127.mat ' ...
                     ;'results/comparisons/Fig7_3_gap4_1128_1127.mat ' ...                     
                    ];        

   case 1
        % Nov. 19
        fileFig_1 = 'results/converge_maxiter_150_1119_1723.mat';
        fileFig_2 = 'results/comparisons/Fig2_gap21_1120_0905.mat';
        fileFig_5a = ['results/comparisons/Fig5a_5_gap8_1119_1500.mat  ' ...
                     ;'results/comparisons/Fig5a_20_gap11_1119_1511.mat' ...
                     ;'results/comparisons/Fig2_gap16_1119_1519.mat    ' ...
                     ;'results/comparisons/Fig5a_50_gap21_1119_1532.mat' ...
                     ;'results/comparisons/Fig5a_60_gap25_1119_1512.mat' ...
                     ;'results/comparisons/Fig5a_80_gap30_1119_1524.mat' ...
                    ];
   case 2
      % Nov. 20
    fileFig_1 = 'results/converge_maxiter_150_1119_1723.mat';
    fileFig_2 = 'results/comparisons/Fig2_gap21_1120_0905.mat';
    fileFig_5a = ['results/comparisons/Fig5a_5_gap8_1119_1500.mat  ' ...
                 ;'results/comparisons/Fig5a_65_gap34_1120_0304.mat' ...
                 ;'results/comparisons/Fig5a_95_gap49_1120_0433.mat' ...
                ];
    fileFig_5b = ['results/comparisons/Fig5b_5_gap6_1120_0354.mat  ' ...
                 ;'results/comparisons/Fig5b_35_gap12_1120_0802.mat' ...
    %              ;'results/comparisons/Fig5b_65_gap17_1120_0808.mat' ...
                 ;'results/comparisons/Fig5b_95_gap36_1120_0502.mat' ...
                ]; 
    case 3
        % Nov. 22
        disp('Results on Nov. 22');
        fileFig_1 = 'results/Fig1_maxiter_200_1123_0038.mat';
        fileFig_2 = 'results/comparisons/Fig2_gap10_1122_0225.mat';
        fileFig_5a = ['results/comparisons/Fig5a_5_gap5_1121_2204.mat  ' ...
                     ;'results/comparisons/Fig2_gap10_1122_0225.mat    ' ...
                     ;'results/comparisons/Fig5a_65_gap19_1121_2159.mat' ...
                     ;'results/comparisons/Fig5a_95_gap26_1121_2158.mat' ...
                    ];
        fileFig_5b = ['results/comparisons/Fig5b_5_gap5_1121_2155.mat  ' ...
                     ;'results/comparisons/Fig5b_35_gap8_1122_0242.mat ' ...
                     ;'results/comparisons/Fig5b_65_gap12_1122_0242.mat' ...
                     ;'results/comparisons/Fig5b_95_gap18_1121_2157.mat' ...
                    ]; 
                
        fileFig_6a = ['results/comparisons/Fig6_1_gap12_1122_0247.mat' ...                     
                     ;'results/comparisons/Fig6_3_gap11_1122_0246.mat' ...
                     ;'results/comparisons/Fig6_2_gap12_1122_0245.mat' ...
                     ;'results/comparisons/Fig2_gap10_1122_0225.mat  ' ...
                    ]; 
    case 4
        % Nov. 25
        disp('Results on Nov. 25');
        fileFig_1 = 'results/Fig1_maxiter_200_1123_0038.mat';
        
        fileFig_2 = 'results/comparisons/Fig2_M__gap3_1125_1443.mat';
        
        fileFig_5a = ['results/comparisons/Fig5a_5_gap2_1125_1446.mat  ' ...
                     ;'results/comparisons/Fig2_M__gap3_1125_1443.mat  ' ...
                     ;'results/comparisons/Fig5a_65_gap5_1125_1445.mat ' ...
                     ;'results/comparisons/Fig5b_95_gap5_1125_1126.mat ' ...
                    ];
                
        fileFig_5b = ['results/comparisons/Fig5b_5_gap2_1125_1128.mat ' ...
                     ;'results/comparisons/Fig5b_35_gap3_1125_1126.mat' ...
                     ;'results/comparisons/Fig5b_65_gap4_1125_1120.mat' ...
                     ;'results/comparisons/Fig5a_95_gap6_1125_1439.mat' ...
                    ]; 
                
        fileFig_6a = ['results/comparisons/Fig6_1_gap3_1125_1441.mat ' ...                     
                     ;'results/comparisons/Fig6_3_gap2_1125_1440.mat ' ...
                     ;'results/comparisons/Fig6_2_gap3_1125_1444.mat ' ...
                     ;'results/comparisons/Fig2_M__gap3_1125_1443.mat' ...
                     ;'results/comparisons/Fig6_4_gap5_1125_1444.mat ' ... % uniform
                    ];  
    case 5
        % Nov. 26
        disp('Results on Nov. 26');
        fileFig_1 = 'results/Fig1_maxiter_50_1127_1406.mat';
        
        fileFig_2 = 'results/comparisons/Fig2_M__gap3_1125_1443.mat';
        
        fileFig_5a = ['results/comparisons/Fig5a_5_gap2_1125_1446.mat  ' ...
                     ;'results/comparisons/Fig2_M__gap3_1125_1443.mat  ' ...
                     ;'results/comparisons/Fig5b_65_gap4_1126_1225.mat ' ...
                     ;'results/comparisons/Fig5b_95_gap5_1126_1226.mat ' ...
                    ];
                
        fileFig_5b = ['results/comparisons/Fig5b_5_gap2_1125_1128.mat ' ...
                     ;'results/comparisons/Fig5b_35_gap3_1125_1126.mat' ...
                     ;'results/comparisons/Fig5b_65_gap4_1125_1120.mat' ...
                     ;'results/comparisons/Fig5a_95_gap6_1125_1439.mat' ...
                    ]; 
                
        fileFig_6a = ['results/comparisons/Fig6_1_gap3_1125_1441.mat ' ...                     
                     ;'results/comparisons/Fig6_3_gap2_1125_1440.mat ' ...
                     ;'results/comparisons/Fig6_2_gap3_1125_1444.mat ' ...
                     ;'results/comparisons/Fig2_M__gap3_1125_1443.mat' ...
                     ;'results/comparisons/Fig6_4_gap5_1125_1444.mat ' ... % uniform
                    ]; 
  
   otherwise
      disp('no result found');
end

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

%% Figure 1: convergence
figure1 = figure;
axes1 = axes('Parent',figure1);
% hold(axes1,'on');

load(fileFig_1);
titleStr = 'Convergence';
xIteration = 1:MAX_ITER;
if MAX_ITER> 300
    semilogx(xIteration, gradMatrix');
else
    plot(xIteration, gradMatrix(1:4,:)', 'linewidth',LineWidth);
end
% legendStr = {'\nabla_{queueing^l_1}', '\nabla_{queueing^l_2}', '\nabla_{queueing^l_3}', '\nabla_{queueing^l_4}', '\nabla_{queueing^l_5}', '\nabla_{queueing^l_6}', '\nabla_{queueing^l_7}', '\nabla_{queueing^l_8}', '\nabla_{queueing^l_9}', '\nabla_{queueing^l_{10}}'};
% legendStr = {'\nabla_{queueing^l_1}', '\nabla_{queueing^l_2}', '\nabla_{queueing^l_3}', '\nabla_{queueing^l_4}', '\nabla_{queueing^l_5}'};
legendStr = {'\nabla_{q^l_1}', '\nabla_{q^l_2}', '\nabla_{q^l_3}', '\nabla_{q^l_4}'};
ylim([-max(max(abs(gradMatrix))) max(max(abs(gradMatrix)))]);
% axis([1 10 0 3.2]);
set(gca,'fontsize',fontAxis);
legend(legendStr,'Location','northeast','FontSize',fontLegend); %,'Orientation','horizontal'
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 4.0 3.0]);
xlabel('iterations','FontSize',fontAxis);
ylabel('gradient','FontSize',fontAxis);
print ('-depsc', [figPath 'grad_converge.eps']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot figure of objective values
figure1 = figure;
axes1 = axes('Parent',figure1);
% hold(axes1,'on');

if MAX_ITER> 300
    semilogx(xIteration, objValues);
else
    plot(xIteration, objValues, 'linewidth', LineWidth);
end
legendStr = {'long-term objective'};
% axis([1 10 0 3.2]);
set(gca,'fontsize',fontAxis);
legend(legendStr,'Location','northeast','FontSize',fontLegend,'Orientation','horizontal');
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 4.0 3.0]);
xlabel('iterations','FontSize',fontAxis);
ylabel('cost','FontSize',fontAxis);
print ('-depsc', [figPath 'obj_converge.eps']);

%% Figure 2: Cost comparison
disp('Figure 2: Cost comparison');
titleStr = 'Cost comparison';
load(fileFig_2)

figure1 = figure;
axes1 = axes('Parent',figure1);
% hold(axes1,'on');

% y = [NoLT; fixedLT; pred_alg; expect_SGEA; exp_offline_cost]; lengendStr = {'NoLT', 'fixed long-term', 'PA', 'SGA', 'OA'};
y = [NoLT_NoGLB_ignoreQDelay; fixedLT_NoGLB; NoLT;  fixedLT; pred_alg; expect_SGEA; exp_offline_cost];  
% y = [NoLT_NoGLB_ignoreQDelay(1:2); fixedLT_NoGLB(1:2); NoLT(1:2);  fixedLT(1:2); pred_alg(1:2); expect_SGEA(1:2); exp_offline_cost(1:2)];  
xLabels = {'nLTnGLB', 'fLTnGLB', 'nLT', 'fLT', 'PA', 'SGA', 'OA'};
legendStr = {'long-term energy cost', 'real-time energy cost', ...
        'queueing delay cost', 'network delay cost'};
bar_chart = bar(y, 'stacked','barwidth',0.5);
set(bar_chart,{'FaceColor'},{'r';'b';'g';brown}); 
axis([0.5 7.5 0 10.5e6]);
set(gca,'fontsize',fontAxis);
legend(legendStr,'Location','northeast','FontSize',fontLegend/2,'Orientation','vertical');
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 4.0 3.0]);
set(gca,'XTickLabel',xLabels,'FontSize',fontAxis/2);
xlabel('algorithms','FontSize',fontAxis/2);
ylabel('cost','FontSize',fontAxis/2);
print ('-depsc', [figPath 'cost_comparison.eps']);
 


%% Figure 5: Show impact of percentage of renewable energy \alpha in data centers on cost performance. 
% The impact of different types of renewable energy (different error distributions) are also interesting.

% Figure 5b solar
disp('Figure 5b - solar');
titleStr = 'solar impact';

fileLen = length(fileFig_5b(:,1));
legendStr = {'long-term energy cost', 'real-time energy cost', ...
        'queueing delay cost', 'network delay cost'};
y = 0;
Y = zeros(fileLen,2,4);
xAxisStr = cell(fileLen,1);
for iFile = 1:fileLen
    load(fileFig_5b(iFile,:)); 
    scales(iFile) = W_RATIO*100;
%     Y(iFile,1,:) = NoLT;
%     Y(iFile,2,:) = fixedLT;
    Y(iFile,1,:) = pred_alg;
    Y(iFile,2,:) = expect_SGEA;
    xAxisStr{iFile} = ['PA&SGA ' num2str(scales(iFile)) '%'];
    
    ltPercentPA(iFile) = sum(q_l_pred)/(sum(M))*100;
    ltPercentSGA(iFile) = sum(q_l_SGEA)/(sum(M))*100;
%     ltPercentPA(iFile) = sum(q_l_pred)/(sum(L_mean./mu))*100;
%     ltPercentSGA(iFile) = sum(q_l_SGEA)/(sum(L_mean./mu))*100;
end

% groupLabels = xAxisStr;     % set labels
% plotBarStackGroups(Y, groupLabels, legendStr); % plot groups of stacked bars
% %title(titleStr);

% Cost curve
figure1 = figure;
axes1 = axes('Parent',figure1);

plot(scales,sum(Y(:,2,1:4),3),TotalCostColor, 'LineWidth', LineWidth);
hold on;
plot(scales,sum(Y(:,2,1:2),3), EnergyCostColor,'LineWidth', LineWidth);
hold on;
plot(scales,Y(:,2,1), LTCostColor,'LineWidth', LineWidth);

legendStr = {TotalCostString, EnergyCostString, LTCostString}
axis([0 100 0 6e6]);
set(gca,'fontsize',fontAxis);
legend(legendStr,'Position',[0.57,0.5,0.25,0.1],'FontSize',fontLegend,'Orientation','vertical'); % [left, bottom, width, height]
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 4.0 3.0]);
xlabel('ratio (%)','FontSize',fontAxis);
ylabel('cost','FontSize',fontAxis);
print ('-depsc', [figPath 'solar_impact.eps']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% percentages of long-term energy
figure1 = figure;
axes1 = axes('Parent',figure1);

plot(scales,ltPercentPA, PAColor, 'LineWidth', LineWidth);
hold on;
plot(scales,ltPercentSGA, SGAColor, 'LineWidth', LineWidth);

legendStr ={PAString,SGAString};
axis([0 100 0 100]);
set(gca,'fontsize',fontAxis);
legend(legendStr,'Location','northeast','FontSize',fontLegend,'Orientation','vertical'); 
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 4.0 3.0]);
xlabel('ratio (%)','FontSize',fontAxis);
ylabel('long-term procurement (%)','FontSize',fontAxis);
print ('-depsc', [figPath 'solar_impact2.eps']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 5a - wind
disp('Figure 5a - wind');
titleStr = 'wind impact';
fileLen = length(fileFig_5a(:,1));
% fileLen = max(fileLen, scales);

y = 0;
Y = zeros(fileLen,2,4);
xAxisStr = cell(fileLen,1);
for iFile = 1:fileLen
    load(fileFig_5a(iFile,:)); 
    scales(iFile) = W_RATIO*100;
%     Y(iFile,1,:) = NoLT;
%     Y(iFile,2,:) = fixedLT;
    Y(iFile,1,:) = pred_alg;
    Y(iFile,2,:) = expect_SGEA;
    xAxisStr{iFile} = ['PA&SGA ' num2str(scales(iFile)) '%'];
    
    ltPercentPA(iFile) = sum(q_l_pred)/(sum(M))*100;
    ltPercentSGA(iFile) = sum(q_l_SGEA)/(sum(M))*100;
    ltPercentSGA(iFile) = sum(q_l_SGEA)/(sum(M))*100;
%     ltPercentPA(iFile) = sum(q_l_pred)/(sum(L_mean./mu))*100;
%     ltPercentSGA(iFile) = sum(q_l_SGEA)/(sum(L_mean./mu))*100;
end

% groupLabels = xAxisStr;     % set labels
% plotBarStackGroups(Y, groupLabels, legendStr); % plot groups of stacked bars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cost curve
figure1 = figure;
legendStr = {TotalCostString, EnergyCostString, LTCostString}
% plot(scales,sum(Y(:,1,1:4),3),TotalCostColor, 'LineWidth', LineWidth);
% hold on;
plot(scales,sum(Y(:,2,1:4),3),TotalCostColor, 'LineWidth', LineWidth);
hold on;
plot(scales,sum(Y(:,2,1:2),3),EnergyCostColor, 'LineWidth', LineWidth);
hold on;
plot(scales,Y(:,2,1),LTCostColor, 'LineWidth', LineWidth);

axis([0 100 0 6e6]);
set(gca,'fontsize',fontAxis);
legend(legendStr,'Position',[0.57,0.5,0.25,0.1],'FontSize',fontLegend,'Orientation','vertical'); % [left, bottom, width, height]
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 4.0 3.0]);
xlabel('ratio (%)','FontSize',fontAxis);
ylabel('cost','FontSize',fontAxis);
print ('-depsc', [figPath 'wind_impact.eps']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% percentages of long-term energy
figure1 = figure;
axes1 = axes('Parent',figure1);

plot(scales,ltPercentPA, PAColor, 'LineWidth', LineWidth);
hold on;
plot(scales,ltPercentSGA, SGAColor, 'LineWidth', LineWidth);
ylim([0 100]);
%title(titleStr);
legendStr ={PAString,SGAString};
axis([0 100 0 100]);
set(gca,'fontsize',fontAxis);
legend(legendStr,'Location','northeast','FontSize',fontLegend,'Orientation','vertical'); 
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 4.0 3.0]);
xlabel('ratio (%)','FontSize',fontAxis);
ylabel('long-term procurement (%)','FontSize',fontAxis);
print ('-depsc', [figPath 'wind_impact2.eps']);


%% Figure 6
disp('Figure 6a - Distribution');
titleStr = 'Prediction impact';
fileLen = length(fileFig_6a(:,1));
legendStr = {'long-term energy cost', 'real-time energy cost', ...
        'queueing delay cost', 'network delay cost'};
% fileLen = max(fileLen, scales);
y = 0;
Y = zeros(fileLen,4);
xAxisStr = cell(fileLen,1);
for iFile = 1:fileLen
    load(fileFig_6a(iFile,:)); 
%     Y(iFile,1,:) = NoLT;
%     Y(iFile,2,:) = fixedLT;
%     Y(iFile,1,:) = pred_alg;
    Y(iFile,:) = expect_SGEA;
    if strcmp(distName,'generalized extreme value')
        distName = 'GEV';
    end
    xAxisStr{iFile} = distName;
end

% Plot distritions

% Plot performance comparisons
% 
% groupLabels = xAxisStr;     % set labels
% plotBarStackGroups(Y, groupLabels, legendStr); % plot groups of stacked bars

figure1 = figure;
axes1 = axes('Parent',figure1);

y = sum(Y,2);
% y =  Y(:,2,:);
barChart = bar(y, 0.4); legend('PA','SGA'); ylabel('cost');
set(barChart,{'FaceColor'},{'b'}); 
% [im_hatch,colorlist] = applyhatch_pluscolor(gcf,'\-x.',0,0,[],150);
% imwrite(im_hatch,'im_hatch.png','png');
legendStr = {'SGA'};
xLabelStr = {'normal', 'real' , 'GEV', 'tlocationscale', 'uniform'};
set(gca,'XTickLabel', xLabelStr,'FontSize',fontAxis);
axis([0.5 3.5 0 11e6]);
set(gca,'fontsize',fontAxis);
legend(legendStr,'Location','northwest','FontSize',fontLegend); %,'Orientation','horizontal'
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 4.0 3.5]);
xlabel('distributions','FontSize',fontAxis);
ylabel('cost','FontSize',fontAxis);
print ('-depsc', [figPath 'distribution_impact.eps']);
%% Figure 6
disp('Figure 6b - Standard Deviation (RMSE)');
titleStr = 'Prediction impact';
fileLen = length(fileFig_7(:,1));
legendStr = {'long-term energy cost', 'real-time energy cost', ...
        'queueing delay cost', 'network delay cost'};
% fileLen = max(fileLen, scales);
y = 0;
Y = zeros(fileLen,4);
X = zeros(fileLen,4);
xAxisStr = cell(fileLen,1);
scales = zeros(fileLen,1);
for iFile = 1:fileLen
    load(fileFig_7(iFile,:)); 
    X(iFile,:) = pred_alg;
    Y(iFile,:) = expect_SGEA;
    xAxisStr{iFile} = [num2str(round(P_STD*100)) '%'];
    scales(iFile)  = round(P_STD*100);
end

figure1 = figure;
axes1 = axes('Parent',figure1);

% plot(scales, sum(X,2) );
% hold on;
% plot(scales, sum(Y,2) );
% hold on;
% ylim([0 inf]);
barChart = bar(Y,'stacked','barwidth',0.5); xlim([0 5]);
set(barChart,{'FaceColor'},{'r';'b';'g';brown}); 
% barChart = bar(X,'stacked','barwidth',0.5); xlim([0 5]);
xlim([0 6]);
ylabel('cost'); xlabel('RMSE');
% legend('PA','SGA');
set(gca,'XTickLabel', xAxisStr,'FontSize',fontAxis);
axis([0.5 4.5 0 11e6]);
set(gca,'fontsize',fontAxis);
% columnlegend(2, legendStr, 'Location','northwest');
legend(legendStr,'Location','northwest','FontSize',fontLegend);  %,'Orientation','horizontal'
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 4.0 3.5]);
xlabel('RMSE','FontSize',fontAxis);
ylabel('cost','FontSize',fontAxis);
print ('-depsc', [figPath 'rmse_impact.eps']);

close all;