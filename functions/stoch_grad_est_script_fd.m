grad = zeros(N,1);    
for idx = 1:MC_samples     
    q_l_left  = q_l_best;
    [optVal, rt_energy_cost, queuing_delay_cost, network_delay_cost,m_left, lambda_left, isSuccessful] = ...                    
        geo_load_balance(q_l_left, w_rSamples(:,idx), p_rSamples(:,idx), L_rSamples(:,idx), mu ,M , beta ,pi_ij);
    if isSuccessful<=0 
        idx
        q_l_left
        error('Cannot solve this geo_load_balance for q_l_left');                     
    end
    F_left = optVal + p_l'*q_l_left;  
    objValues(iter) = objValues(iter) + F_left;
    
    for i = 1:N   
        q_l_right = q_l_best;
        q_l_right(i) = q_l_best(i) + ck;        
        gap = q_l_right(i) - q_l_left(i);          
        
        [optVal, rt_energy_cost, queuing_delay_cost, network_delay_cost, m_right, lambda_right, isSuccessful] = ...
            geo_load_balance(q_l_right, w_rSamples(:,idx), p_rSamples(:,idx), L_rSamples(:,idx), mu ,M , beta ,pi_ij);
        if isSuccessful<=0
            idx
            q_l_right
            error('Cannot solve this geo_load_balance for q_l_right'); 
        end
        F_right = optVal + p_l'*q_l_right;
        
        % Compute the mean of gradient. 
        grad(i) = grad(i) + (F_right - F_left)./gap;
        
        if PROGRESS_BAR
            progressbarCount = progressbarCount+1;
            progressbar(progressbarCount/(MAX_ITER*N*MC_samples));
        end
    end    
end
grad = grad./MC_samples;
objValues(iter) = objValues(iter)/MC_samples;