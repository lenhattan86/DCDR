function [lt_energy_cost_fixedLT_NoGLB, rt_energy_cost_fixedLT_NoGLB, queuing_delay_cost_fixedLT_NoGLB, network_delay_cost_fixedLT_NoGLB] ...
    = fixedLT_noGLB( p_l, q_l, w_rSamples, p_rSamples, L_rSamples, mu ,M , beta ,pi_ij)

    MC_samples = length(w_rSamples(1,:));   
    
    lt_energy_cost_fixedLT_NoGLB = 0;
    rt_energy_cost_fixedLT_NoGLB = 0;    
    queuing_delay_cost_fixedLT_NoGLB = 0;
    network_delay_cost_fixedLT_NoGLB = 0;
    
    for i=1:MC_samples
        [rt_energy_cost, queuing_delay_cost, network_delay_cost, m, lambda, isSuccessful] = ...
                    nearest_routing(q_l, w_rSamples(:,i), p_rSamples(:,i), L_rSamples(:,i), mu ,M , beta ,pi_ij);   
                
        lt_energy_cost_fixedLT_NoGLB = lt_energy_cost_fixedLT_NoGLB + p_l'*q_l;        
        rt_energy_cost_fixedLT_NoGLB = rt_energy_cost_fixedLT_NoGLB + rt_energy_cost;    
        queuing_delay_cost_fixedLT_NoGLB = queuing_delay_cost_fixedLT_NoGLB + queuing_delay_cost;
        network_delay_cost_fixedLT_NoGLB = network_delay_cost_fixedLT_NoGLB + network_delay_cost;
        
    end
    lt_energy_cost_fixedLT_NoGLB = lt_energy_cost_fixedLT_NoGLB/MC_samples;
    rt_energy_cost_fixedLT_NoGLB = rt_energy_cost_fixedLT_NoGLB/MC_samples;    
    queuing_delay_cost_fixedLT_NoGLB = queuing_delay_cost_fixedLT_NoGLB/MC_samples;
    network_delay_cost_fixedLT_NoGLB = network_delay_cost_fixedLT_NoGLB/MC_samples;
end

