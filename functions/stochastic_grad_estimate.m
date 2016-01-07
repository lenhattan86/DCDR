function [gradientVector] = ...
stochastic_grad_estimate(p_l, q_l, partial_q, w_rSamples, p_rSamples, L_rSamples)
%STOCHASTIC_GRAD_DESCENT Summary of this function goes here
%   Estimate the stochastic gradient using finite difference estimator 
% [book: Handbook of simulation optimization -chapter 5] 

% TODO: improve the speed of gradient estimation as N times
    global DEBUG
    
    MC_samples = length(w_rSamples(1,:));
    
    % Foreach sample do steps 1 & 2
    gradientVector = zeros(N,1);     
    for i = 1:N
        q_l_right = q_l;
        q_l_left  = q_l;
        gradTemp = 0;
        for idx = 1:MC_samples      
            if q_l(i) > 0
                q_l_left(i) = q_l(i) - partial_q;  
            end

            q_l_right(i) = q_l(i) + partial_q;
            
            gap = q_l_right(i) - q_l_left(i);
            
            [optVal , rt_energy_cost, queuing_delay_cost, network_delay_cost, m, lambda, isSuccessful] = ...                    
                geo_load_balance(q_l_left, w_rSamples(:,idx), p_rSamples(:,idx), L_rSamples(:,idx), mu ,M , beta ,pi_ij);
            if isSuccessful<0 
                idx
                q_l_left
                error('Cannot solve this geo_load_balance for q_l_left');                     
            end
            F_left = optVal + p_l'*q_l_left;            
            [optVal , rt_energy_cost, queuing_delay_cost, network_delay_cost, m, lambda, isSuccessful] = ...
                geo_load_balance( q_l_right, w_rSamples(:,idx), p_rSamples(:,idx), L_rSamples(:,idx), mu ,M , beta ,pi_ij);
            if isSuccessful<0
                idx
                q_l_right
                error('Cannot solve this geo_load_balance for q_l_right'); 
            end
            F_right = optVal + p_l'*q_l_right;
            % Compute the mean of gradient. 
            gradTemp = gradTemp + (F_right - F_left)./gap; 
        end
        gradientVector(i) = gradTemp/(totalSamples);
    end
end