function [q_l_pred, m, lambda, lt_energy_cost, rl_energy_cost, queueing_delay_cost, network_delay_cost, isSuccessful] ...
    = pred_procurement(p_l, w_r_pred, p_r_pred, L_r_pred, mu, M, beta ,pi_ij)
%PREDICTION_PROCUREMENT Summary of this function goes here
%   provide the long term procurement plan based on predition based
%   information
    [q_l_pred,  m, lambda, lt_energy_cost, rl_energy_cost, queueing_delay_cost, network_delay_cost,  isSuccessful] ...
        = offline_procurement(p_l, w_r_pred, p_r_pred, L_r_pred, mu, M, beta ,pi_ij);
end

