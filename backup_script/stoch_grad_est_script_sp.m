MC_samples = length(w_rSamples(1,:));
grad = zeros(N,1);

DELTA = [ones(1,N/2) (-1)*ones(1,N/2)];
%delta=2*round(rand(N,1))-1;
delta=DELTA(randperm(N))';

q_l_plus   = q_l_best + ck*delta;
q_l_minus  = q_l_best - ck*delta;
gradTemp = 0;

% yplus=loss(thetaplus);
% yminus=loss(thetaminus);
for idx = 1:MC_samples      
    [optVal, rt_energy_cost, queuing_delay_cost, network_delay_cost, m_left, lambda_left, isSuccessful] = ...                    
        geo_load_balance(q_l_minus, w_rSamples(:,idx), p_rSamples(:,idx), L_rSamples(:,idx));
    if isSuccessful<0 
        idx
        q_l_minus
        error('Cannot solve this geo_load_balance for q_l_minus');                     
    end
    F_left = optVal + p_l'*q_l_minus;    

    [optVal, rt_energy_cost, queuing_delay_cost, network_delay_cost, m_right, lambda_right, isSuccessful] = ...
        geo_load_balance(q_l_plus, w_rSamples(:,idx), p_rSamples(:,idx), L_rSamples(:,idx));
    if isSuccessful<0
        idx
        q_l_plus
        error('Cannot solve this geo_load_balance for q_l_plus'); 
    end
    F_right = optVal + p_l'*q_l_plus;

    % Compute the mean of gradient. 
    gradTemp = gradTemp + (F_right - F_left)./(2*ck*delta);

    if PROGRESS_BAR
        progressbarCount = progressbarCount+1;
        progressbar(progressbarCount/(MAX_ITER*MC_samples));
    end
end
% ghat=(yplus-yminus)./(2*ck*delta);
% theta=theta-ak*ghat;
grad = gradTemp./MC_samples;