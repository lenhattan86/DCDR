function [total_cost, rt_energy_cost, queueing_delay_cost,network_delay_cost, q_r_sum, q_renew] = ...
    expected_cost( p_l, q_l, w_rSamples, p_rSamples, L_rSamples, mu ,M , beta ,pi_ij )
%EXPECTED_COST Summary of this function goes here
%   Detailed explanation goes here
    global DEBUG     
    total_cost =  0;
    rt_energy_cost = 0;
    queueing_delay_cost = 0;
    network_delay_cost = 0;
    MC_samples = length(w_rSamples(1,:));
    q_r_sum = 0;
    q_renew = 0;
    for idx = 1:MC_samples
        w_r = w_rSamples(:,idx);
        p_r = p_rSamples(:,idx);
        L_r = L_rSamples(:,idx);      
        [optVal, opt_rt_energy_cost, opt_queueing_delay_cost, opt_network_delay_cost , m, lambda, isSuccessful] = ...                    
            geo_load_balance(q_l, w_r, p_r, L_r, mu ,M , beta , pi_ij);
        if isSuccessful<0           
            error('Cannot solve this geo_load_balance for q_l');            
        end
        % Compute the mean of cost.
        total_cost = total_cost + optVal; 
        rt_energy_cost = rt_energy_cost + opt_rt_energy_cost;
        queueing_delay_cost = queueing_delay_cost + opt_queueing_delay_cost;
        network_delay_cost = network_delay_cost + opt_network_delay_cost;
        q_r = pos(m-q_l-w_r);
        q_r_sum = q_r_sum + sum(q_r);
        q_renew = q_renew + sum(pos(m-q_l-q_r));
    end
    total_cost = total_cost/(MC_samples)  + p_l'*q_l;
    rt_energy_cost = rt_energy_cost/(MC_samples);
    queueing_delay_cost = queueing_delay_cost/(MC_samples);
    network_delay_cost = network_delay_cost/(MC_samples);
    q_r_sum = q_r_sum/MC_samples;
    q_renew = q_renew/MC_samples;
end