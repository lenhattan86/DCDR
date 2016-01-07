%% Trade off
close all;

beta_array = [0.4 1 1.4 2];
no_long_term =  [3.0644e+07 6.7061e+07 9.2653e+07 1.2941e+08];
pred_based =    [2.7630e+07 6.3335e+07 8.8630e+07 1.2535e+08];
offline_based = [2.6940e+07 6.2749e+07 8.8066e+07 1.2479e+08];
worstOffline =  [4.4118e+07 8.9941e+07 1.3125e+08 1.8674e+08];
bestOffline =   [1.7592e+07 4.1227e+07 5.9241e+07 8.3901e+07];
long_term =     [2.7626e+07 6.3335e+07 8.8629e+07 1.2535e+08];

proposedStr = 'Proposed Algorithm';
withoutLongTermStr = 'Without long-term procurement';

no_long_termLineStyle = '-r';
pred_basedLineStype = '>--';
offlineLineStype = 'o-k';
long_termLineStyle = 'x-.b';

worstOfflineLineStyle = '-.';
bestOfflineLineStyle = ':';


figure(1);


plot(beta_array, no_long_term, no_long_termLineStyle, 'LineWidth',1);
hold on;
plot(beta_array, pred_based, pred_basedLineStype,'LineWidth',1);
hold on;
plot(beta_array, offline_based, offlineLineStype,'LineWidth',1);
hold on;
plot(beta_array, long_term, long_termLineStyle,'LineWidth',1);
hold on;
plot(beta_array, worstOffline, worstOfflineLineStyle,'LineWidth',1);
hold on;
plot(beta_array, bestOffline, bestOfflineLineStyle,'LineWidth',1);
hold on;


xlabel('\beta');
ylabel('cost') % left y-axis
title('Operation cost vs. \beta');
legend(withoutLongTermStr, 'Prediction based', 'Offline mean',  proposedStr, 'worst offline', 'best offline');
grid on;
%%

figure(2);

beta_array = [1 1.4 2];
no_long_term =  [4.9371e+04  7.0948e+04 8.1633e+04];
pred_based =    [4.9034e+04 7.0235e+04 7.9612e+04];
long_term =     [4.9036e+04 7.0262e+04 7.9710e+04];

plot(beta_array, no_long_term, no_long_termLineStyle, 'LineWidth',1);
hold on;
plot(beta_array, pred_based, pred_basedLineStype,'LineWidth',1);
hold on;
plot(beta_array, long_term, long_termLineStyle,'LineWidth',1);


xlabel('\beta');
ylabel('cost') % left y-axis
title('Delay cost vs. \beta');
legend(withoutLongTermStr, 'Prediction based', proposedStr);
grid on;

%%

figure(3);

beta_array = [1 1.4 2];
no_long_term =  [7.9952e+04 7.7185e+04 9.9665e+04];
pred_based =    [6.5834e+04 6.6999e+04 8.1094e+04];
long_term =     [6.5825e+04 6.6972e+04 8.0995e+04];

plot(beta_array, no_long_term, no_long_termLineStyle, 'LineWidth',1);
hold on;
plot(beta_array, pred_based, pred_basedLineStype,'LineWidth',1);
hold on;
plot(beta_array, long_term, long_termLineStyle,'LineWidth',1);

xlabel('\beta');
ylabel('cost') % left y-axis
title('Energy cost vs. \beta');
legend(withoutLongTermStr, 'Prediction based', 'Offline mean',  proposedStr, 'worst offline', 'best offline');
grid on;