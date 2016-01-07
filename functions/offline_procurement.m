function [q_l_offline,  m, lambda, lt_energy_cost, rl_energy_cost, queueing_delay_cost, network_delay_cost,  isSuccessful] ...
    = offline_procurement(p_l, w_r, p_r, L_r, mu, M, beta ,pi_ij)
%OFFLINE_PROCUREMENT Summary of this function goes here
%   Detailed explanation goes here
    global DEBUG PROGRESS_BAR QUEUEING_DELAY_SCALE
   
    % prediction based offline optimal solution
    % long term + GLB-RT
    N = length(p_l);
    J = length(L_r);
    cvx_begin quiet
        variables q_l_offline(N) m(N) lambda(N,J);
        minimize(p_l'* q_l_offline ...
            + p_r'*pos(m-w_r-q_l_offline) ...                        
            + beta * QUEUEING_DELAY_SCALE * sum(sum(lambda,2) + quad_over_lin(sum(lambda,2), mu.*m-sum(lambda,2),0)) ...         
            + beta * (sum(sum(lambda.*pi_ij))))        
        subject to
            lambda >= 0;
            sum(lambda,2) <= m.*mu;       
            m >= 0;
            m <= M;        
            sum(lambda,1)' == L_r;
            q_l_offline >= 0;
    cvx_end
    isSuccessful = 1;
    if (~(or(strcmp(cvx_status,'Solved'),strcmp(cvx_status,'Inaccurate/Solved'))))        
        cvx_status        
%         optVal = (p_r'*pos(m-q_l-w_r) ...                        
%             + beta * sum(inv_pos(mu .* lambda_i - inv_pos(m))) ...         
%             + beta * (sum(sum(lambda.*pi_ij))));
        isSuccessful = 0;
    end
    lt_energy_cost = p_l'* q_l_offline;
    rl_energy_cost = p_r'*pos(m-q_l_offline-w_r);
    queueing_delay_cost = beta* QUEUEING_DELAY_SCALE *sum(sum(lambda,2) + quad_over_lin(sum(lambda,2), mu.*m-sum(lambda,2),0));
    network_delay_cost = beta*(sum(sum(lambda.*pi_ij)));
end

