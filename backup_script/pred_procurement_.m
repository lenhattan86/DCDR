function [lt_energy_cost, rl_energy_cost, queueing_delay_cost, network_delay_cost, isSuccessful] ...
    = pred_procurement(p_l, w_r_pred_set, p_r_pred_set, L_r_pred_set,  w_rSamples, p_rSamples, L_rSamples  )
%PREDICTION_PROCUREMENT Summary of this function goes here
%   provide the long term procurement plan based on predition based
%   information
    MC_samples = length(w_r_pred_set(1,:));
    lt_energy_cost = 0;
    rl_energy_cost = 0;
    queueing_delay_cost = 0;
    network_delay_cost = 0;   
    isSuccessful = 1;
    for i = 1:MC_samples
        [q_l_pred,  m, lambda, lt_temp, rt_temp, queueing_temp, network_temp,  flag] ...
            = offline_procurement(p_l, w_r_pred_set(:,i), p_r_pred_set(:,i), L_r_pred_set(:,i) );
        if flag == 0
            isSuccessful = 0;
        end        
        % compute the cost with q_l_pred
        [optVal, rt_temp, queueing_temp, network_temp, m, lambda, isSuccessful] = ...
                    geo_load_balance(q_l_pred, w_rSamples(:,i), p_rSamples(:,i), L_rSamples(:,i));
                
        lt_energy_cost = lt_energy_cost + lt_temp;
        rl_energy_cost = rl_energy_cost + rt_temp;
        queueing_delay_cost = queueing_delay_cost + queueing_temp;
        network_delay_cost = network_delay_cost + network_temp;
        
        if flag == 0
            isSuccessful = 0;
        end
    end
    lt_energy_cost = lt_energy_cost/MC_samples;
    rl_energy_cost = rl_energy_cost/MC_samples;
    queueing_delay_cost = queueing_delay_cost/MC_samples;
    network_delay_cost = network_delay_cost/MC_samples;   
end

