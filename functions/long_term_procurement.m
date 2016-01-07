function [ expect_cost, q_l_best, grad_norm, gradMatrix, objValues] = ...
    long_term_procurement(q_l_0, p_l, w_rSamples, p_rSamples, L_rSamples, mu ,M , beta ,pi_ij, STEP_SIZE_SETTINGS)
%LONG_TERM_PROCUREMENT Summary of this function goes here
%   Detailed explanation goes here

    global DEBUG PROGRESS_BAR    
    global MAX_ITER GRADIENT_ESTIMATE 
    MC_samples = length(w_rSamples(1,:));
    
    A = STEP_SIZE_SETTINGS(1);
    a = STEP_SIZE_SETTINGS(2); 
    alp = STEP_SIZE_SETTINGS(3); 
    
    C =STEP_SIZE_SETTINGS(4);    
    gamma = STEP_SIZE_SETTINGS(5);
            
    N = length(w_rSamples(:,1));
    J = length(L_rSamples(:,1));

    % initialize q_l = {0}, 
    q_l_best = q_l_0;
    iter = 1;
    stepIter = ones(N,1);
    dq = inf;
    grad_norm = inf;   
%     expect_cost = inf;
    if PROGRESS_BAR
        progressbar 
        progressbarCount = 0;
    end
    grad = zeros(N,1);
    prev_grad = 0*grad;
    step_size = ones(N,1);
    gradMatrix = inf*ones(N,MAX_ITER);
    objValues  = zeros(MAX_ITER,1);    
    while iter <= MAX_ITER        
        % calculate gradient:
        %[grad] = stochastic_grad_estimate(p_l, q_l_best, partial_q, ...
        %    w_rSamples, p_rSamples, L_rSamples);
        
        % Call script instead of  function "stochastic_grad_estimate"
        for i = 1:N
            if prev_grad(i)*grad(i)<=0
                step_size(i) = a./((stepIter(i)+A)^alp);  
                stepIter(i) = stepIter(i) + 1;
            end
        end        
        if strcmp(GRADIENT_ESTIMATE, 'CVX')
            prev_grad=grad;        
            stoch_grad_est_lagrange
            gradMatrix(:,iter) = grad; 
        elseif strcmp(GRADIENT_ESTIMATE, 'FD')        
            ck=C./(iter^gamma);
            prev_grad=grad;            
            stoch_grad_est_script_fd % Finite Difference
            gradMatrix(:,iter) = grad; 
        elseif strcmp(GRADIENT_ESTIMATE, 'SP')            
            stoch_grad_est_script_sp % Simutaneous Pertubation
        end
        
        grad_norm = norm(grad);        
        % take step:
        q_l_new = q_l_best - step_size.*grad;             
        q_l_best = max(q_l_new,0); % project q_l_best on the feasible domain.         
        
        if DEBUG            
            if N > 1
                figure(2);
                plot((1:iter),gradMatrix(:,(1:iter)));
                xlim([0 MAX_ITER+30]);
                legendStr = {'\nabla_{q^l_1}', '\nabla_{q^l_2}', '\nabla_{q^l_3}', '\nabla_{q^l_4}', '\nabla_{q^l_5}', '\nabla_{q^l_6}', '\nabla_{q^l_7}', '\nabla_{q^l_8}', '\nabla_{q^l_9}', '\nabla_{q^l_{10}}'};
                legend(legendStr);   
            
                figure(3);
                scatter((iter),objValues(iter));
                xlim([0 MAX_ITER+1]);
                hold on;
                drawnow
                q_l_best
                grad
            else
                grad
                q_l_best
            end            

        end
        iter = iter + 1; 
    end  
    
            % update expected total cost
    [expect_temp, rt_energy_cost_temp, delay_cost_temp] = expected_cost( p_l, q_l_best, ...
        w_rSamples, p_rSamples, L_rSamples, mu ,M , beta ,pi_ij);
    %         cost_gap = expect_cost - expect_temp;
    expect_cost = expect_temp;
    
    if PROGRESS_BAR
        progressbar(1);
    end
end

