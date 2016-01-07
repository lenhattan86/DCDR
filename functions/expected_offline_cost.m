function [ exp_offline_cost, offline_worst_cost, offline_best_cost, energy ] = expected_offline_cost( p_l, q_l, w_rSamples, p_rSamples, L_rSamples ,mu, M ,beta, pi_ij)
%EXPECTED_OFFLINE_COST Summary of this function goes here
%   Detailed explanation goes here

    global DEBUG

    offline_worst_cost = -inf*ones(1,4);
    offline_best_cost = inf*ones(1,4);
     energy = zeros(2,1);
    lt_e_cost = 0;
    rt_e_cost = 0;
    queueing_delay_cost = 0;
    network_delay_cost = 0;
    MC_samples = length(w_rSamples(1,:));
    N = length(p_l);
    J = length(L_rSamples(:,1));
    
    for idx = 1:MC_samples      
        w_r = w_rSamples(:,idx);
        p_r = p_rSamples(:,idx);
        L_r = L_rSamples(:,idx);
        [ q_l_offline,  m, lambda, lt_e_cost_temp, rt_e_cost_temp, queueing_delay_cost_temp, network_delay_cost_temp,  isSuccessful ] ...
                = offline_procurement(p_l, w_r, p_r, L_r , mu, M, beta ,pi_ij);
            
        lt_e_cost = lt_e_cost + lt_e_cost_temp;
        rt_e_cost = rt_e_cost + rt_e_cost_temp;
        queueing_delay_cost = queueing_delay_cost + queueing_delay_cost_temp;
        network_delay_cost = network_delay_cost + network_delay_cost_temp;

        optVal = lt_e_cost_temp + rt_e_cost_temp  + queueing_delay_cost_temp + network_delay_cost_temp;

        % update worst & best cases        
        if optVal > sum(offline_worst_cost)
            offline_worst_cost = [lt_e_cost_temp  rt_e_cost_temp queueing_delay_cost_temp network_delay_cost_temp];
        end
        if optVal < sum(offline_best_cost)
            offline_best_cost =  [lt_e_cost_temp rt_e_cost_temp queueing_delay_cost_temp network_delay_cost_temp];
        end
        q_r = pos(m-q_l_offline-w_r);
        energy(1) = energy(1) + sum(q_l_offline) + sum(q_r);
        energy(2) = energy(2) + sum(pos(m-q_r-q_l_offline));
    end    
    exp_offline_cost = [lt_e_cost/(MC_samples) rt_e_cost/(MC_samples) queueing_delay_cost/(MC_samples) network_delay_cost/(MC_samples)];
    energy = energy./MC_samples;
end



