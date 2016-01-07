% Vary all uncertain parameters to find the best settings.
clear; clear global variable; close all; 
clc;
DONT_RESET = 0;
GENERATE_DATA = 1;
DISABLE_PROGRESS_BAR = 1;

p_r_std_array    = [0 5 10 15 20 25 30]./100;
w_std_array      = [0 5 10 15 20 25 30]./100;
L_std_array      = [0 5 10 15 20 25 30]./100;
%p_r_SCALE_array  = [1 1.5 2 2.5 3 3.5];

cnt = 0;
total_run = length(p_r_std_array)*length(w_std_array)*length(L_std_array);
for iP=1:length(p_r_std_array)
    GENERATE_p_r = 1;
    p_r_std = p_r_std_array(iP);
    for iW=1:length(w_std_array)
        GENERATE_w_r = 1;
        w_std = w_std_array(iW);
        for iL=1:length(L_std_array)
            GENERATE_L_r = 1;
            L_std = L_std_array(iL);
            extraName = ['_' num2str(p_r_std) '_' num2str(w_std) '_' num2str(L_std)];
            ScriptMethodComparison
            GENERATE_DATA = 0;
            cnt = cnt + 1;
            progressbar(cnt/total_run);
        end
    end
end