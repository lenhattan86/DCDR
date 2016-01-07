function [optVal, rt_energy_cost, queuing_delay_cost, network_delay_cost, m, lambda, isSuccessful, q_lMultipliers] = ...
            geo_load_balance(q_l, w_r, p_r, L_r, mu , M , beta , pi_ij)
% GEO_LOAD_BALANCE This function is to compute the optimal solution for
% geo-graphical load balancing problem.
% Requirement:  CVX
% Input: 
% q_l: is the long term energy procurement
% w_r: is the renewable energy generated inside data centers
% p_r: real-time electricity prices at data centers
% L_r: real-time workload from sources J
% mu: service rate
% beta: the importance of delay
% pi_ij: newwork delay
%
% Output: 
% m: number of active servers.
% lambda: distributed workload. 

    global DEBUG PROGRESS_BAR QUEUEING_DELAY_SCALE 

    N = length(w_r);
    J = length(L_r);

    %             + beta * (sum(sum(lambda.*pi_ij))) ...
    % sum(lambda,2) + quad_over_lin(sum(lambda,2), mu.*m-sum(lambda,2),0) <= 2;
    cvx_begin quiet
        variables q_r(N) m(N)  lambda(N,J)
        dual variable q_lMultipliers
        minimize(p_r'*q_r)        
        subject to
            lambda >= 0;
            sum(lambda,2) <= m.*mu;       
            m >= 0;
            m <= M;        
            sum(lambda,1)' == L_r;
            q_r >= 0;
            q_lMultipliers: -q_r+m-w_r-q_l <= 0;
            sum(lambda,2) + quad_over_lin(sum(lambda,2), mu.*m-sum(lambda,2),0) <= 10/6*sum(lambda,2);
    cvx_end
    optVal = cvx_optval;    
    rt_energy_cost = p_r'*pos(m-q_l-w_r);
    queuing_delay_cost = beta * QUEUEING_DELAY_SCALE * sum(sum(lambda,2) + quad_over_lin(sum(lambda,2), mu.*m-sum(lambda,2),0));
    network_delay_cost = beta * (sum(sum(lambda.*pi_ij)));
    isSuccessful = 1;
    if (~(or(strcmp(cvx_status,'Solved'),strcmp(cvx_status,'Inaccurate/Solved'))))        
        cvx_status        
        sum(L_r)
        sum(mu.*M)
        error('This optimization is unsolved');
        isSuccessful = -1;
    end    
        
%     utilizationRate =  mean(sum(lambda,2)./m); 
%     utilizatedServer = mean(m./M);
%     if utilizationRate < 0.3 || utilizationRate > 0.56
%         utilizationRate        
% %         error('Need to scale network delay and queing delay');
%     end
end