function [energy] = computeEnergyConsumption(p_l, q_l, w_rSamples, p_rSamples, L_rSamples, mu ,M , beta ,pi_ij)
    [total_cost, rt_energy_cost, queueing_delay_cost,network_delay_cost, q_r_sum, q_renew] = ...
        expected_cost( p_l, q_l, w_rSamples, p_rSamples, L_rSamples, mu ,M , beta ,pi_ij );
    energy = zeros(2,1);
    energy(1) = sum(q_l)+q_r_sum;
    energy(2) = q_renew;
end

