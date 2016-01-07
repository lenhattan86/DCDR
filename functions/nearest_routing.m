function [rt_energy_cost, queuing_delay_cost, network_delay_cost, m, lambda, isSuccessful] = ...
                    nearest_routing(q_l, w_r, p_r, L_r, mu ,M , beta ,pi_ij)
    global DEBUG PROGRESS_BAR QUEUEING_DELAY_SCALE UTIL_RATE
    
	J = length(L_r);
    N = length(w_r);
    lambda = zeros(N,J);
    % solve m and lambda without queueing delay cost.
    % find the closest data centers.
    closestDC = zeros(J,1);
    for j = 1:1:J
        minDist = inf;
        for i = 1:1:N            
            if minDist > pi_ij(i,j)
                minDist = pi_ij(i,j);
                closestDC(j) = i;            
            end
        end
        lambda(closestDC(j),j) = L_r(j);
    end  
    Lambda_i = sum(lambda,2);
    m = (1/UTIL_RATE)*Lambda_i./mu;
    m = max(m,M); % turn on all server
    if(min(m) == 0)
        error('One of data center is not active');
    end
    
    % compute energy cost + delay cost
    rt_energy_cost = p_r'*pos(m-q_l-w_r);
    queuing_delay_cost = beta * QUEUEING_DELAY_SCALE * sum(sum(lambda,2) + quad_over_lin(sum(lambda,2), mu.*m-sum(lambda,2),0));
%     queuing_delay_cost = beta * sum(Lambda_i./(mu-Lambda_i./m));
    network_delay_cost = beta * (sum(sum(lambda.*pi_ij)));
    isSuccessful = 1;
end
