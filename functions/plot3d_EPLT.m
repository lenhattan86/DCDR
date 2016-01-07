function [ bestCost, q_l_best, X, Y, Z ] = plot3d_EPLT(d_q, q_l_max, p_l, w_rSamples, p_rSamples, L_rSamples)
    
    bestCost = inf;
    q_l_best = zeros(2,1);
    q_l_vector = zeros(2,1);
    X = 0:d_q:q_l_max;
    Y = X; 
    Z = ones(length(X));
    idx = 1;
    for q_l_1 = 0:d_q:q_l_max
        q_l_vector(1) =  q_l_1; 
        idy = 1;
        for q_l_2 = 0:d_q:q_l_max
            q_l_vector(2) = q_l_2; 
            [total_cost, rt_energy_cost, delay_cost] = expected_cost( p_l, q_l_vector, w_rSamples, p_rSamples, L_rSamples );
            Z(idx, idy) = total_cost;
            if bestCost > total_cost
                bestCost = total_cost;
                q_l_best = q_l_vector;
            end    
            idy = idy + 1;
        end
        idx = idx + 1;    
    end
end