grad = zeros(N,1);  

for idx = 1:MC_samples     
    [optVal, rt_energy_cost, queuing_delay_cost, network_delay_cost,m_left, lambda_left, isSuccessful, q_lMultipliers] = ...                    
    geo_load_balance(q_l_best, w_rSamples(:,idx), p_rSamples(:,idx), L_rSamples(:,idx), mu ,M , beta ,pi_ij);    
    gradTemp = 0;
    grad = grad - q_lMultipliers;
    objValues(iter) = objValues(iter) + optVal + p_l'*q_l_best;
    if PROGRESS_BAR
        progressbarCount = progressbarCount+1;
        progressbar(progressbarCount/(MAX_ITER*MC_samples));
    end
end
grad = grad./MC_samples + p_l;
objValues(iter) = objValues(iter)/MC_samples;