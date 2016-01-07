%% TEST CASE 1
if TEST_CASE == 1
    % p_l > p_r, w_r = 0, 
    % p_r, w_r, L_r are constants
    % So, expected q_l must be 0;
    disp('p_l > p_r, w_r = 0 ; p_r, w_r, L_r are constants');
    expeted_q_l = zeros(N,1);    
    if VERIFY_RESULTS
        % expected output
        [optVal m lambda isSuccessful] = ...
                geo_load_balance(q_l_best, w_rSamples(:,1), p_rSamples(:,1), L_rSamples(:,1))
%         expectect_cost = expected_cost(p_l, q_l_best, w_rSamples, p_rSamples, L_rSamples);    
        if or(q_l_best ~= expeted_q_l, ~isSuccessful)
            display('FAILED');
            q_l_best
            expeted_q_l    
            l
        else
            display('PASSED');
            q_l_best
        end
        %expected_optimal_cost
        VERIFY_RESULTS = 0;     
    else
        N = 2; % Number of Data centers
        J = 4; % Number of sources
        M = 1000 * ones (N,1); % Capacity: Number of servers in a data center.        
        pi_ij = ones(N,J); % network delay
        mu = 10*ones(N,1); % service rate    

        q_l = zeros(N, 1); % long term procurement   
        q_l_0 = 10*ones(N, 1);
        
        p_l = 5*ones(N, 1); % long term prices    

        p_rSamples = 3*ones(N, MC_samples); % real-time prices    
        L_rSamples = 20*ones(J, MC_samples); % Amount of workload at source
        w_rSamples = zeros(N, MC_samples); % Renewable energy generation.    
        VERIFY_RESULTS = 1;
    end    
%% TEST CASE 2    
elseif TEST_CASE == 2
    % p_l > p_r,
    % p_r, w_r, L_r are random variables
    % So, expected q_l must be 0;
    disp('p_l > p_r; p_r, w_r, L_r are random variables ');
    expeted_q_l = zeros(N,1);

    if VERIFY_RESULTS
        [optVal m lambda isSuccessful] = ...
            geo_load_balance(q_l_best, w_rSamples(:,1), p_rSamples(:,1), L_rSamples(:,1))
        % expected output
        if q_l_best == expeted_q_l
            display('PASSED');
            q_l_best
        else
            display('FAILED');
            q_l_best
            expeted_q_l
        end
        %expected_optimal_cost
    else
        N = 2; % Number of Data centers
        J = 4; % Number of sources
        M = 1000 * ones (N,1); % Capacity: Number of servers in a data center.
        
        pi_ij = ones(N,J); % network delay
        mu = 10*ones(N,1); % service rate
        mu_min = 100;
        mu_max = 200;

        q_l = zeros(N, 1); % long term procurement   
        q_l_0 = 10*ones(N, 1);

        p_l = zeros(N, 1); % long term prices
        p_l_max = 5;
        p_l_min = 4;
        

        p_rSamples = zeros(N, MC_samples); % real-time prices
        p_max = 3;
        p_min = 4;
        L_rSamples = zeros(J, MC_samples); % Amount of workload at source
        L_max = 20;
        L_min = 15;
        w_rSamples = zeros(N, MC_samples); % Renewable energy generation.
        w_max = 10;
        w_min = 5;

        p_l = (p_l_max-p_l_min)*rand(N,1) + p_l_min;
        for i = 1:N
            % renewable
            w_rSamples(i,:) = (w_max-w_min)*rand(1,MC_samples) + w_min;
            % prices
            p_rSamples(i,:) = (p_max-p_min)*rand(1,MC_samples) + p_min;
        end
        for j = 1:J
            % workload demand at source j
            L_rSamples(j,:) = (L_max-L_min)*rand(1,MC_samples) + L_min;
        end    
        pi_ij = rand(N,J);
        mu = (mu_max-mu_min)*rand(N,1) + mu_min;
        VERIFY_RESULTS = 1;
    end    
%% TEST CASE 3    
elseif TEST_CASE == 3
    % p_l < p_r, 
    % p_r, w_r, L_r are constants
    % So, expected q_l must be maximum, and q_r = 0
    disp('p_l < p_r; p_r, w_r, L_r are constants');
    expeted_q_r = zeros(N,1);
    
    if VERIFY_RESULTS
        % expected output
        [optVal m lambda isSuccessful] = ...
                geo_load_balance(q_l_best, w_rSamples(:,1), p_rSamples(:,1), L_rSamples(:,1))
        [glb_cost, exp_q_r] = expected_glb( q_l_best, w_rSamples, p_rSamples, L_rSamples );
        if exp_q_r == expeted_q_r
            display('PASSED');
            exp_q_r
            glb_cost
        else
            display('FAILED');
            exp_q_r
            glb_cost
        end
        %expected_optimal_cost
    else
        N = 2; % Number of Data centers
        J = 4; % Number of sources
        M = 100 * ones (N,1); % Capacity: Number of servers in a data center.
        
        pi_ij = ones(N,J); % network delay
        mu = 100*ones(N,1); % service rate       

        q_l = zeros(N, 1); % long term procurement   
        q_l_0 = 10*ones(N, 1);       
        
        q_r = 10*ones(N, 1); % real-time procurement

        p_l = zeros(N, 1); % long term prices
        p_l_max = 8;
        p_l_min = 6; 

        p_rSamples = 3*ones(N, MC_samples); % real-time prices    
        L_rSamples = 20*ones(J, MC_samples); % Amount of workload at source
        w_rSamples = zeros(N, MC_samples); % Renewable energy generation.    
        
        VERIFY_RESULTS = 1;
    end    

%% TEST CASE 4 

elseif TEST_CASE == 4
    % E[p_l] < E[p_r]
    % p_r, w_r = 0, L_r are random         
    display('avg[p_l] < E[p_r] & p_r, w_r = 0, L_r are random');
    if VERIFY_RESULTS
        % expected output
        q_l_best
        [glb_cost, exp_q_r] = expected_glb( q_l_best, w_rSamples, p_rSamples, L_rSamples )
        total_cost = glb_cost + p_l'*q_l_best
%         if exp_q_r == expeted_q_r
%            display('PASSED');
%            exp_q_r
%             glb_cost
%         else
%             display('FAILED');
%         end
        %expected_optimal_cost
    else
        N = 5; % Number of Data centers
        J = 20; % Number of sources
        M = 4000 * ones (N,1); % Capacity: Number of servers in a data center.
        
        pi_ij = ones(N,J); % network delay
        mu = 10*ones(N,1); % service rate

        q_l = zeros(N, 1); % long term procurement   
        q_l_0 = 0*ones(N, 1);       

        p_l = zeros(N, 1); % long term prices
        p_l_max = 8;
        p_l_min = 3;

        p_rSamples = zeros(N, MC_samples); % real-time prices
        p_max = 11;
        p_min = 5;
        L_rSamples = zeros(J, MC_samples); % Amount of workload at source
        L_max = 100;
        L_min = 50;
        w_rSamples = zeros(N, MC_samples); % Renewable energy generation.
        w_max = 0;
        w_min = 0;

        p_l = (p_l_max-p_l_min)*rand(N,1) + p_l_min;
        for i = 1:N
            % renewable
            w_rSamples(i,:) = (w_max-w_min)*rand(1,MC_samples) + w_min;
            % prices
            p_rSamples(i,:) = (p_max-p_min)*rand(1,MC_samples) + p_min;
        end
        for j = 1:J
            % workload demand at source j
            L_rSamples(j,:) = (L_max-L_min)*rand(1,MC_samples) + L_min;
        end    
        pi_ij = 2*rand(N,J);
        VERIFY_RESULTS = 1;           
    end  
    
elseif TEST_CASE == 5
    % p_l < p_r & p_r, w_r, L_r are random     
    display('p_l < p_r & p_r, w_r = 0, L_r are random');
    
    if VERIFY_RESULTS
        % expected output
        q_l_best
        [glb_cost, exp_q_r] = expected_glb( q_l_best, w_rSamples, p_rSamples, L_rSamples )
        total_cost = glb_cost + p_l'*q_l_best
        w_r = w_rSamples(:,1); p_r = p_rSamples(:,1); L_r = L_rSamples(:,1);
        [optVal m lambda isSuccessful] = ...
                geo_load_balance(q_l_best, w_r, p_r, L_r)
%         if exp_q_r == expeted_q_r
%            display('PASSED');
%            exp_q_r
%             glb_cost
%         else
%             display('FAILED');
%         end
        %expected_optimal_cost
    else
        N = 2; % Number of Data centers
        J = 4; % Number of sources
        M = 1000 * ones (N,1); % Capacity: Number of servers in a data center.
        
        pi_ij = ones(N,J); % network delay
        mu = 100*ones(N,1); % service rate

        q_l = zeros(N, 1); % long term procurement   
        q_l_0 = 0*ones(N, 1);       

        p_l = zeros(N, 1); % long term prices
        p_l_max = 5;
        p_l_min = 3;

        p_rSamples = zeros(N, MC_samples); % real-time prices
        p_max = 20;
        p_min = 10;
        L_rSamples = zeros(J, MC_samples); % Amount of workload at source
        L_max = 100;
        L_min = 50;
        w_rSamples = zeros(N, MC_samples); % Renewable energy generation.
        w_max = 0;
        w_min = 0;

        p_l = (p_l_max-p_l_min)*rand(N,1) + p_l_min;
        for i = 1:N
            % renewable
            w_rSamples(i,:) = (w_max-w_min)*rand(1,MC_samples) + w_min;
            % prices
            p_rSamples(i,:) = (p_max-p_min)*rand(1,MC_samples) + p_min;
        end
        for j = 1:J
            % workload demand at source j
            L_rSamples(j,:) = (L_max-L_min)*rand(1,MC_samples) + L_min;
        end    
        pi_ij = 2*rand(N,J);
        VERIFY_RESULTS = 1;           
    end  
elseif TEST_CASE == 6
    % normal cases with small number of data centers 
    disp('normal cases with small number of data centers ');
    
    if VERIFY_RESULTS
        % expected output
        total_cost = expected_cost(p_l, q_l_best, w_rSamples, p_rSamples, L_rSamples)
        q_l_best     
       
%         if exp_q_r == expeted_q_r
%            display('PASSED');
%            exp_q_r
%             glb_cost
%         else
%             display('FAILED');
%         end
        %expected_optimal_cost
    else
        N = 2; % Number of Data centers
        J = 4; % Number of sources
        M = 1000 * ones (N,1); % Capacity: Number of servers in a data center.
        
        pi_ij = ones(N,J); % network delay
        mu = 1*ones(N,1); % service rate
        mu_min = 1;
        mu_max = 2;

        q_l = zeros(N, 1); % long term procurement   
        q_l_0 = 0*ones(N, 1);       

        p_l = zeros(N, 1); % long term prices
        p_l_max = 10;
        p_l_min = 5;

        p_rSamples = zeros(N, MC_samples); % real-time prices
        p_max = 20;
        p_min = 5;
        L_rSamples = zeros(J, MC_samples); % Amount of workload at source
        L_max = 350;
        L_min = 100;
        w_rSamples = zeros(N, MC_samples); % Renewable energy generation.
        w_max = 0;
        w_min = 0;

        p_l = (p_l_max-p_l_min)*rand(N,1) + p_l_min;
        for i = 1:N
            % renewable
            w_rSamples(i,:) = (w_max-w_min)*rand(1,MC_samples) + w_min;
            % prices
            p_rSamples(i,:) = (p_max-p_min)*rand(1,MC_samples) + p_min;
        end
        for j = 1:J
            % workload demand at source j
            L_rSamples(j,:) = (L_max-L_min)*rand(1,MC_samples) + L_min;
        end    
        pi_ij = rand(N,J);
        mu = (mu_max-mu_min)*rand(N,1) + mu_min;
        VERIFY_RESULTS = 1;           
    end 
%%    
elseif TEST_CASE == 7
    % E[p_l] < E[p_r]
    % p_r, w_r, L_r are random      
    
    if VERIFY_RESULTS
        % expected output
        q_l_best
        [glb_cost, exp_q_r] = expected_glb( q_l_best, w_rSamples, p_rSamples, L_rSamples )
        total_cost = glb_cost + p_l'*q_l_best
        w_r = w_rSamples(:,1); p_r = p_rSamples(:,1); L_r = L_rSamples(:,1);
        [optVal m lambda isSuccessful] = ...
                geo_load_balance(q_l_best, w_r, p_r, L_r)
%         if exp_q_r == expeted_q_r
%            display('PASSED');
%            exp_q_r
%             glb_cost
%         else
%             display('FAILED');
%         end
        %expected_optimal_cost
    else
        N = 20; % Number of Data centers
        J = 40; % Number of sources
        M = 2000 * ones (N,1); % Capacity: Number of servers in a data center.
        
        pi_ij = ones(N,J); % network delay
        mu = 100*ones(N,1); % service rate
        mu_min = 100;
        mu_max = 200;

        q_l = zeros(N, 1); % long term procurement   
        q_l_0 = 0*ones(N, 1);       

        p_l = zeros(N, 1); % long term prices
        p_l_max = 8;
        p_l_min = 3;

        p_rSamples = zeros(N, MC_samples); % real-time prices
        p_max = 12;
        p_min = 3;
        L_rSamples = zeros(J, MC_samples); % Amount of workload at source
        L_max = 500;
        L_min = 50;
        w_rSamples = zeros(N, MC_samples); % Renewable energy generation.
        w_max = 5;
        w_min = 1;

        p_l = (p_l_max-p_l_min)*rand(N,1) + p_l_min;
        for i = 1:N
            % renewable
            w_rSamples(i,:) = (w_max-w_min)*rand(1,MC_samples) + w_min;
            % prices
            p_rSamples(i,:) = (p_max-p_min)*rand(1,MC_samples) + p_min;
        end
        for j = 1:J
            % workload demand at source j
            L_rSamples(j,:) = (L_max-L_min)*rand(1,MC_samples) + L_min;
        end    
        pi_ij = 2*rand(N,J);
        mu = (mu_max-mu_min)*rand(N,1) + mu_min;
        VERIFY_RESULTS = 1;           
    end
end