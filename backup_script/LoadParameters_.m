

%% Simulation parameters
global MAX_ITER partial_q STEP_SIZE CK GRADIENT_ESTIMATE
global N J mu M beta pi_ij 
global MC_samples;

MC_samples    = 30; % Number of MC   Samples.
MAX_ITER      = 30;

% FD for low dimensional problem, SP for high dimensional problem
GRADIENT_ESTIMATE = 'FD'; % SP or FD

partial_q  = 1;

beta = 0.01;
DELAY_SCALE = 0.5;
p_r_SCALE = 2.5;
PRICE_SCALE = 1;

if RUN_TEST_CASES
    disp(['Run test case: ', num2str(TEST_CASE)]);
    %DEBUG=1;
    RunTestCase
else
    
    datacenter_location = [37.390073 122.081852;
    47.237054	119.852629;
    45.600287	121.143165;
    39.791127	89.647865;
    31.515337	82.850018;
    37.554376	77.433815;
    29.784641	95.364304;
    27.976211	82.455368;
    35.920474	81.540928;
    33.165145	80.002899];

    state_location = [32 87; 65 150; 34 111; 35 92; 37 120; 39 105; 41.5 72.5; ...
        39 75.5; 39 77; 28 81; 33 83; 21 158; 44 114; 40 89; 40 86; 42 93; 38 98; ...
        38 85; 31 92; 45 69; 39 76.5; 42.5 72; 43 84; 46 94.5; 33 90; 38 92; 47 110; ...
        41 99; 39 116; 43.5 71.5; 40 74.5; 34 106; 43 76; 35 79; 47 100; 40 83; ...
        35 97; 44 120; 41 78; 41.5 71.5; 34 81; 44 100; 36 86; 31 99; 40 112; ...
        44 73; 38 79; 47 121; 39 81; 44 89; 43 107];

    N = 10; % Number of Data centers
    J = 20; % Number of sources
    M = 10000 * ones (N,1); % Capacity: Number of servers in a data center.  
    
     pi_ij = ones(N,J); % network delay
%     pi_ij_avg = 10;
%     pi_ij_std  = 2;
%     pi_ij = max(pi_ij_std.*randn(N,J) + pi_ij_avg, 0);
    for i = 1:1:N
        for j = 1:1:J
            pi_ij(i,j)=sqrt((datacenter_location(i,1)-state_location(j,1))^2+(datacenter_location(i,2)-state_location(j,2))^2);
        end
    end   
    pi_ij = pi_ij * DELAY_SCALE;
    
    mu = 100*ones(N,1);
%     mu_avg = 100;
%     mu_std = 20;
    %mu = max(mu_std.*randn(N,1) + mu_avg, 0);
    

    q_l = zeros(N, 1); % long term procurement
    q_l_max = 10;  
    q_l_0 = q_l;
    q_l_pred = q_l;
    
    %p_l = 10*randn(N,1) + 100;
    p_l = [92  101   102   105   115   102    93    95    92   103]';
%     p_l = [87  93   108   112]';
    p_l = PRICE_SCALE*p_l;
    
    p_rSamples = zeros(N, MC_samples); % real-time prices
    p_r_mean = p_r_SCALE*p_l;
    if ~exist('p_r_std')
        %p_r_std =  0.3 * mean(p_r_mean);
        p_r_std = 75;% max 75/2.5*100
    end   
    
    L_mean = 400000; % max L_mean should be less than M*mu*N/J
    if ~exist('L_std')
        L_std =  120000; %max: 130000/400000
    end
    
    w_mean = 3000;
    if ~exist('w_std')
        w_std = 1100; % max: 1100/3000
    end
    
    %% prediction errors        
    w_r_pred_set = zeros(N,MC_samples);
    w_error_std = 20/100 * w_mean;
    
    p_r_pred_set = zeros(N,MC_samples);
    p_r_error_std = 20/100 * p_r_mean;    
    
    L_r_pred_set = zeros(J,MC_samples);
    L_error_std =  20/100 * L_mean;     
    %% II. Generate data or Load from files
    if GENERATE_DATA
        for i = 1:N
            % renewable
            w_rSamples(i,:) = max(w_std.*randn(1,MC_samples) + w_mean,0);
            %w_rSamples(i,:) = w_std.*randn(MC_samples,1) + w_mean;
            if(min(w_rSamples(i,:)) < 0)
                error('Please decrease w_std');
            end
            w_r_pred_set(i,:)  = max(w_error_std.*randn(1,MC_samples) + w_rSamples(i,:), 0);
            
            % prices
            p_rSamples(i,:) = max(p_r_std.*randn(1,MC_samples) + p_r_mean(i),0);
            if(min(p_rSamples(i,:)) < 0)
                error('Please decrease p_r_std');
            end
            p_r_pred_set(i,:)  = max(p_r_error_std(i).*randn(1,MC_samples) + p_rSamples(i,:),0);
        end
        for j = 1:J
            % workload demand at source j
            L_rSamples(j,:) = max(L_std.*randn(1,MC_samples) + L_mean,0);
            if(min(L_rSamples(j,:)) < 0)
                error('Please decrease L_std');
            end
            L_r_pred_set(j,:)  = max(w_error_std.*randn(1,MC_samples) + L_rSamples(j,:),0); 
        end
        
        save('data/gen/uncentainty.mat');
    else
        load('data/gen/uncentainty.mat');
        for i = 1:N            
            % renewable
            if exist('GENERATE_w_r')
                w_rSamples(i,:) = max(w_std.*randn(1,MC_samples) + w_mean,0);
                if(min(w_rSamples(i,:)) < 0)
                    error('Please decrease w_std');
                end
            end
            % prices
            if exist('GENERATE_p_r')
                p_rSamples(i,:) = max(p_r_std.*randn(1,MC_samples) + p_r_mean(i),0);
                if(min(p_rSamples(i,:)) < 0)
                    error('Please decrease p_r_std');
                end
            end
        end
        for j = 1:J
            % workload demand at source j
            if exist('GENERATE_L_r')
                L_rSamples(j,:) = max(L_std.*randn(1,MC_samples) + L_mean,0);
                if(min(L_rSamples(j,:)) < 0)
                    error('Please decrease L_std');
                end
            end
        end
    end   
end