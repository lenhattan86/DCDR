function [glb_cost, exp_q_r] = expected_glb( q_l, w_rSamples, p_rSamples, L_rSamples, mu ,M , beta ,pi_ij )
%EXPECTED_COST Summary of this function goes here
%   Detailed explanation goes here
    global DEBUG
    global MC_samples
    
    total_cost =  0;
    exp_q_r = zeros(N,1);
    for idx = 1:MC_samples      
        w_r = w_rSamples(:,idx);
        p_r = p_rSamples(:,idx);
        L_r = L_rSamples(:,idx);
        [optVal, rt_energy_cost, queuing_delay_cost, network_delay_cost, m, lambda, isSuccessful] = ...                    
            geo_load_balance(q_l, w_r, p_r, L_r, mu ,M , beta ,pi_ij);
        if ~isSuccessful
            if DEBUG
                error('Cannot solve this geo_load_balance given q_l');            
            end
        end
        % Compute the mean of cost. 
        total_cost = total_cost + optVal;
        exp_q_r = exp_q_r + pos(m-w_r-q_l);
    end 
    glb_cost = total_cost/MC_samples;
    exp_q_r = exp_q_r/MC_samples;
end

