warning('on','all');
alp = 0.602; [ A, a] = getStepSizeSetting(1, 30, 2500, 100, alp) %( Start, End, stepStart, stepEnd, alpha )
STEP_SIZE_SETTINGS =   [A   a  alp   0.1 0.101]
% step_size = a./((iter+A)^alp);
% ck = C./(iter^gamma);
if ~exist('N')
    N = 10; % Number of Data centers
    J = 40; % Number of sources
end

if (~exist('P_STD'))
    P_STD = 0.53;
end
if (~exist('L_STD'))
    L_STD = 0.6;
end
if strcmp(RENEWABLE_TYPE, 'wind')
    if (~exist('W_STD'))
        W_STD = 0.97;
    end
else
    if (~exist('W_STD'))
        W_STD = 0.54;
    end    
end


if (~exist('distName'))
    distName = 'real';
end
% distName = 'normal';
% distName = 'tlocationscale';
% distName = 'generalized extreme value';



PRICE_SCALE = 2.3;
p_l_RATIO = 1/2.5;

LOAD_SCALE = 0.35; % if this value is large, the queueing delay for nearest routing becomes inf
if (~exist('fixedLT'))
    fixedLT = 0.5;
end
NETWORK_DELAY_SCALE = 2;
global QUEUEING_DELAY_SCALE  UTIL_RATE
% QUEUEING_DELAY_SCALE = 10;
QUEUEING_DELAY_SCALE = 1;
UTIL_RATE = 0.5;

%% Simulation parameters
partial_q  = 1;

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

    state_location = [32 87; 65 150; 34 111; 35 92; 37 120; ...
                    39 105; 41.5 72.5; 39 75.5; 39 77; 28 81; ...
                    33 83; 21 158; 44 114; 40 89; 40 86; ...
                    34 81; 44 100; 36 86; 31 99; 40 112; ...
                    44 73; 38 79; 47 121; 39 81; 44 89; 43 107; ...
                    42 93; 38 98; 38 85; 31 92; 45 69; ...
                    39 76.5; 42.5 72; 43 84; 46 94.5; 33 90; ...
                    38 92; 47 110; 41 99; 39 116; 43.5 71.5; ...
                    40 74.5; 34 106; 43 76; 35 79; 47 100; ...
                    40 83; 35 97; 44 120; 41 78; 41.5 71.5; ...
                    ];
%     N = length(datacenter_location);
%     J = length(state_location);
    closestDC = zeros(J,1);
    pi_ij = ones(N,J); % network delay
    for j = 1:1:J
        minDist = inf;
        for i = 1:1:N
            pi_ij(i,j)=sqrt((datacenter_location(i,1)-state_location(j,1))^2+(datacenter_location(i,2)-state_location(j,2))^2);
            if minDist > pi_ij(i,j)
                minDist = pi_ij(i,j);
                closestDC(j) = i;            
            end
        end
    end
    dcLoads = hist(closestDC,N);
    L_mean = 5000;    
    % setup capacity
    MU = 1;
    M = zeros (N,1); % Capacity: Number of servers in a data center.  
    for i=1:N
        M(i) = dcLoads(i)*L_mean/LOAD_SCALE/MU;               
    end

    
    pi_ij = pi_ij * NETWORK_DELAY_SCALE;
    
    mu = MU*ones(N,1);

    q_l = zeros(N, 1); % long term procurement
    q_l_max = 10;  
    q_l_0 = q_l;
    q_l_pred = q_l;
    
    % tx.houston, miso.cin, ca.np15, ohio, aep, nj, ny, chi, bos
%     p_r_mean = [18.57 4.65 26.83  31.15 32.16 36.04 50.66 30.92 26.7 66.81'; % midnight 1
%     p_r_mean = [45.06 48.13 24.74 31.96 32.16 32.7 34.72 88.85 31.77 74.66]'; % noon 12
%     p_r_mean = [45.06 48.13 24.74 31.96 32.16 32.7 34.72 46.85 31.77 41.66]'; % adjust.
    p_r_mean = [10.41 3.73 5.87 7.48 5.86 6.67 6.44 8.60 6.03 5.49]'; % from GLB.
    if N == 1
        p_r_mean = mean(p_r_mean);
    end
    p_r_mean = PRICE_SCALE*p_r_mean;
    p_r_mean = p_r_mean(1:N);
    p_l = p_l_RATIO*p_r_mean;
    p_rSamples = zeros(N, MC_samples); 

    if ~exist('p_r_std')
        %p_r_std =  0.3 * mean(p_r_mean);
        p_r_std = P_STD*p_r_mean; % percent
    end   
    
    if ~exist('L_std')
        L_std =  L_STD*L_mean;% percent 
    end
    
    w_mean = W_RATIO*LOAD_SCALE*M;
    if ~exist('w_std')
        w_std = W_STD*w_mean; % percent 
    end
    
    if strcmp(distName, 'normal') || strcmp(distName, 'uniform')
        L_rMax = 1.95*L_mean; L_rMin = 0.05*L_mean;
        p_rMax = 1.95*p_r_mean; p_rMin = 0.05*p_r_mean;
        w_rMax = 1.99*w_mean; w_rMin = 0.01*w_mean;
    else
        L_rMax = 2.5*L_mean; L_rMin   = 0.05*L_mean;
        p_rMax = 2.5*p_r_mean; p_rMin = 0.05*p_r_mean;
        w_rMax = 2.5*w_mean; w_rMin   = 0.01*w_mean;
    end
        
    %% II. Generate data or Load from files
    w_rSamples = zeros(N,MC_samples);
    p_rSamples = zeros(N,MC_samples);
    L_rSamples = zeros(J,MC_samples);
    if GENERATE_DATA        
        for i = 1:N
            % renewable
            if strcmp(distName,'real')
                w_rSamples(i,:) = genRandValues(RENEWABLE_TYPE, MC_samples, w_mean(i), w_std(i), w_rMin(i), w_rMax(i));
            else
                w_rSamples(i,:) = genRandValues(distName, MC_samples, w_mean(i), w_std(i), w_rMin(i), w_rMax(i));
            end
%             w_rSamples(i,1:MC_samples*0.8) = 0;
%             w_rSamples(i,MC_samples*0.8+1:MC_samples) = 5*max(w_std*w_mean.*randn(1,0.2*MC_samples) + w_mean,0);
            if(min(w_rSamples(i,:)) < 0)
                warning('Please decrease w_std');
                w_rSamples(i,:) = max(w_rSamples(i,:), 0);
            end            
            % prices
            if strcmp(distName,'real')
                p_rSamples(i,:) = genRandValues('price', MC_samples, p_r_mean(i), p_r_std(i), p_rMin(i), p_rMax(i));
            else
                p_rSamples(i,:) = genRandValues(distName, MC_samples, p_r_mean(i), p_r_std(i), p_rMin(i), p_rMax(i));
            end            
            if(min(p_rSamples(i,:)) < 0)
                warning('Please decrease p_r_std');
                p_rSamples(i,:) = max(p_rSamples(i,:), 0);
            end            
        end
        for j = 1:J
            % workload demand at source j
            if strcmp(distName,'real')
                L_rSamples(j,:) = genRandValues('workload', MC_samples, L_mean, L_std, L_rMin, L_rMax);                
            else
                L_rSamples(j,:) = genRandValues(distName, MC_samples, L_mean, L_std, L_rMin, L_rMax);
            end            
            if(min(L_rSamples(j,:)) < 0)
                warning('Please decrease L_std');
                L_rSamples(j,:) = max(L_rSamples(j,:), 0);
            end                
        end  
        
    else
        if ~exist('GEN_FILE')
            load(['data/gen/uncentainty.mat']);
        else
            load(['data/gen/' GEN_FILE]);
        end
        
        for i = 1:N            
            % renewable
            if exist('GENERATE_w_r')
                if strcmp(distName,'real')
                    w_rSamples(i,:) = genRandValues(RENEWABLE_TYPE, MC_samples, w_mean(i), w_std(i), w_rMin(i), w_rMax(i));
                else
                    w_rSamples(i,:) = genRandValues(distName, MC_samples, w_mean(i), w_std(i), w_rMin(i), w_rMax(i));
                end
                if(min(w_rSamples(i,:)) < 0)
                    warning('Please decrease w_std');
                    w_rSamples(i,:) = max(w_rSamples(i,:), 0);
                end 
            end
            % prices
            if exist('GENERATE_p_r')
                if strcmp(distName,'real')
                    p_rSamples(i,:) = genRandValues('price', MC_samples, p_r_mean(i), p_r_std(i), p_rMin(i), p_rMax(i));
                else
                    p_rSamples(i,:) = genRandValues(distName, MC_samples, p_r_mean(i), p_r_std(i), p_rMin(i), p_rMax(i));
                end 
                if(min(p_rSamples(i,:)) < 0)
                    warning('Please decrease p_r_std');
                    p_rSamples(i,:) = max(p_rSamples(i,:), 0);
                end
            end
        end
        for j = 1:J
            % workload demand at source j
            if exist('GENERATE_L_r')
                if strcmp(distName,'real')
                    L_rSamples(j,:) = genRandValues('workload', MC_samples, L_mean, L_std, L_rMin, L_rMax);
                else
                    L_rSamples(j,:) = genRandValues(distName, MC_samples, L_mean, L_std, L_rMin, L_rMax);
                end
                if(min(L_rSamples(j,:)) <= 0)
                    warning('Please decrease L_std');
                    L_rSamples(j,:) = max(L_rSamples(j,:), L_rMin);
                end 
            end
        end
    end   
    
    % Adjust L_rSamples to prevent infeasible for GLB-RT & Nearest Routing
%     for iMC = 1:MC_samples
%         reserverGap = 0.01;        
%         loadGap = sum(L_rSamples(:,iMC)) - sum(M.*mu);        
%         if loadGap  > reserverGap*sum(M.*mu)/J
%             L_rSamples(:,iMC) = L_rSamples(:,iMC) - ceil(loadGap/J) - reserverGap*sum(M.*mu)/J;
%         end
%     end
    realL_STD = std(L_rSamples(1,:))/mean(L_rSamples(1,:))  
    realw_STD = std(w_rSamples(1,:))/mean(w_rSamples(1,:))  
    realp_STD = std(p_rSamples(1,:))/mean(p_rSamples(1,:))
%     mean(L_rSamples(1,:))
%     mean(w_rSamples(1,:))
%     mean(p_rSamples(1,:))
%     hist(w_rSamples(1,:),MC_samples)
%     figure
%     hist(p_rSamples(1,:),MC_samples)
%     figure
%     hist(L_rSamples(1,:),MC_samples)    
    if ~exist('GEN_FILE')
        save('data/gen/uncentainty.mat', 'w_rSamples', 'p_rSamples', 'L_rSamples');
    else
        save(['data/gen/' GEN_FILE], 'w_rSamples', 'p_rSamples', 'L_rSamples');
    end    
end