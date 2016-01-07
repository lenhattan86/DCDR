has_quadprog = exist( 'quadprog' );
has_quadprog = has_quadprog == 2 | has_quadprog == 3;
has_linprog  = exist( 'linprog' );
has_linprog  = has_linprog == 2 | has_linprog == 3;
rnstate = randn( 'state' ); randn( 'state', 1 );
s_quiet = cvx_quiet(true);
s_pause = cvx_pause(false);
cvx_clear;
clc;
clear;

RESHAPING = 1;
RT = 1; %renewable weight, no renewable if RT=1;
solar = 0.9;
renewable = 0.1;
T = 11;
REVERSE = 1; %using reverse price
RANDOM = 0; %using random price
ROUND = 510; %number of rounds per hour
NON_CONST_GAMMA = 0;
ERROR = 0.001;
REUSE = 0 % 1 means reuse the results of last timestep for initial value
%gamma = min(M.^4./(M+init_value).^3/b);
gamma = 10;

a = 10;
b = 51;
I = ones(a,1);
J = ones(b,1);
U = ones(1,a);
beta = 1;
L = load('hot-load1s-norm');
N = size(L);
for i=1:1:N(1)
    L(i,1) = L(i,1)/10000000-12862303208;
end
duration = floor(L(N(1),1)/3600);
hour = 1;
count = 0;
sum_t = 0;


for i=1:1:duration
    time(i) = i;
    average_load(i) = 0;
end
for i=1:1:N(1)
    if L(i,1) <= 3600*hour
        count = count + 1;
        sum_t = sum_t + L(i,2);
    else
        average_load(hour) = sum_t/count;
        hour = hour + 1;
        count = 0;
        sum_t = 0;
    end
end


%plot(t,average_load)

%traffic reshaping
if RESHAPING
    duration = 48;
	average_load(1) = (average_load(48) + average_load(2))/2;
	average_load(25) = (average_load(24) + average_load(26))/2;
	average_load(49) = (average_load(48) + average_load(50))/2;
end

geography = [930 179 1536 585 8614 1354 931 213 163 4921 2261 271 330 3137 1498 767 700 960 875 352 1443 1710 2334 1407 523 1344 241 455 654 399 2252 441 4596 2015 162 2702 762 1031 3038 278 976 197 1333 5040 623 181 2049 1091 369 1556 132];

for i=3:1:duration
    hourly_load_hawaii(i) = average_load(i-2);
end
hourly_load_hawaii(1) = average_load(duration-1);
hourly_load_hawaii(2) = average_load(duration);

for i=1:1:duration
    hourly_load_pacific(i) = average_load(i);
end

for i=1:1:(duration-1)
    hourly_load_mountain(i) = average_load(i+1);
end
hourly_load_mountain(duration) = average_load(1);

for i=1:1:(duration-2)
    hourly_load_central(i) = average_load(i+2);
end
hourly_load_central(duration) = average_load(2);
hourly_load_central(duration-1) = average_load(1);

for i=1:1:(duration-3)
    hourly_load_eastern(i) = average_load(i+3);
end
hourly_load_eastern(duration) = average_load(3);
hourly_load_eastern(duration-1) = average_load(2);
hourly_load_eastern(duration-2) = average_load(1);

time_zone = [2 -2 1 2 0 1 3 3 3 3 3 -2 1 2 3 2 2 3 2 3 3 3 3 2 2 2 1 2 0 3 3 1 3 3 2 3 2 0 3 3 3 2 2 2 1 3 3 0 3 2 1];

for i=1:1:51
    if time_zone(i) == -2
        hourly_load(i,:) = geography(i)*hourly_load_hawaii(:);
    elseif time_zone(i) == 0
        hourly_load(i,:) = geography(i)*hourly_load_pacific(:);
    elseif time_zone(i) == 1
        hourly_load(i,:) = geography(i)*hourly_load_mountain(:);
    elseif time_zone(i) == 2
        hourly_load(i,:) = geography(i)*hourly_load_central(:);
    elseif time_zone(i) == 3
        hourly_load(i,:) = geography(i)*hourly_load_eastern(:);
    end
end

max_load = max(sum(hourly_load));
M = floor(max_load/9)*[5 1 2 1 2 3 1 1 1 1]';

price = [10.41 3.73 5.87 7.48 5.86 6.67 6.44 8.6 6.03 5.49];

%price = [1 1 1 1 1 1 1 1 1 1];

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

state_location = [32 87; 65 150; 34 111; 35 92; 37 120; 39 105; 41.5 72.5; 39 75.5; 39 77; 28 81; 33 83; 21 158; 44 114; 40 89; 40 86; 42 93; 38 98; 38 85; 31 92; 45 69; 39 76.5; 42.5 72; 43 84; 46 94.5; 33 90; 38 92; 47 110; 41 99; 39 116; 43.5 71.5; 40 74.5; 34 106; 43 76; 35 79; 47 100; 40 83; 35 97; 44 120; 41 78; 41.5 71.5; 34 81; 44 100; 36 86; 31 99; 40 112; 44 73; 38 79; 47 121; 39 81; 44 89; 43 107];

for i = 1:1:10
    for j = 1:1:51
        delay(i,j)=sqrt((datacenter_location(i,1)-state_location(j,1))^2+(datacenter_location(i,2)-state_location(j,2))^2);
    end
end

for rt=1:1:RT
    renewable(rt) = 0.01 * (rt - 1);
    
    %intial point
    for i=1:1:a
        for j=1:1:b
            lambda(i,j)=hourly_load(j,11)/a;
        end
    end
    
    m = min([M,lambda*J./U'.*(I./sqrt(price)'+I)]')';
    m_record_initial = m;
    init_value = beta*(price*m) + sum(sum(delay.*lambda)) - sum(m) + sum(quad_over_lin([m';zeros(1,a)],m'-J'*lambda'));
    %least_value = 

    
    
    lip = (init_value./(M+init_value)).*M.*U';

    for t=11:1:T
        % cvx version
        time(t) = t;
        
        %change between times
        if(t>1 && REUSE)
            for i=1:1:a
                for j=1:1:b
                    lambda(i,j)=lambda(i,j)*hourly_load(j,t)/hourly_load(j,t-1);
                end
            end
        end
     
        %optimal value
        cvx_begin
            variables m_opt(a) lambda_opt(a,b);
            minimize( beta*(price*m_opt) + sum(sum(delay.*lambda_opt)) - sum(m_opt) + sum(quad_over_lin([m_opt';zeros(1,a)],m_opt'-J'*lambda_opt')) );
            subject to
                m_opt >= 0;
                lambda_opt >= 0;    
                lambda_opt'*I == hourly_load(:,t);
                m_opt <= M;       
        cvx_end
        m_opt_record(t,rt,:) = m_opt;
        lambda_opt_record(t,rt,:,:) = lambda_opt;
        optimal_record(t,rt) = cvx_optval;
               
        %Distributed Descent
        for round=1:1:ROUND
            round
            temp = zeros(a,b);
            slot = (t-11)*ROUND + round;
            s(slot) = slot;
            for j=1:1:b
                %j
                for i = 1:1:a
                    if lambda(i,:)*J <= sqrt(beta*price(i))/(sqrt(beta*price(i))+1)*M(i)*U(i)
                        F(i,j) = (delay(i,j)+(sqrt(beta*price(i))+1)^2/U(i));
                        F1(i,j) = F(i,j);
                    else
                        F(i,j) = (delay(i,j)+U(i)/(U(i)-lambda(i,:)*J/M(i))^2);
                        F1(i,j) = F(i,j);
                    end
                end
                count = 0;
                sum1 = 0;
                Omega = zeros(1,a);
                for i = 1:1:a
                    if lambda(i,j) > hourly_load(j,t)*ERROR
                        count = count + 1;
                        sum1 = sum1 + F(i,j);
                        F(i,j) = inf;
                        Omega(i) = 1;
                    end
                end
                [FJ,FI] = sort(F(:,j));
                for i = 1:1:a-count
                    if FJ(i) < sum1/count
                        Omega(FI(i)) = 1;
                        sum1 = sum1 + FJ(i);
                        count = count + 1;
                    end
                end
                theta = sum1/count; 
                gamma_curr = gamma;
                for i = 1:1:a
                    if lambda(i,j) > hourly_load(j,t)*ERROR & (F1(i,j) > theta)
                        gamma_curr = min(gamma_curr, lambda(i,j)/(F1(i,j)- theta));
                    end
                end
                for i=1:1:a
                    if Omega(i) == 1
                        lambda(i,j) = lambda(i,j) - gamma_curr * (F1(i,j) - theta);
                    end
                end
            end
            
            
            m = min([M,lambda*J./U'.*(I./sqrt(beta*price')+I)]')';
            m_record(:,t,rt,slot) = m;
            
            
            
            objective(slot,rt) = beta*(price*m) + sum(sum(delay.*lambda)) - sum(m) + sum(quad_over_lin([m';zeros(1,a)],m'-J'*lambda'));              
            %objective(t,rt,round,j) = beta*(P(t,:,rt)*m) + sum(sum(delay.*lambda)) - sum(m) + sum(quad_over_lin([m';zeros(1,a)],m'-J'*lambda'));
            optimal_t(slot,rt) = optimal_record(t,rt);
            lambda_record(:,:,slot,rt) = lambda;
            sum(lambda_record(:,1,slot,rt))
            %lambda_record(:,:,t,rt,round,j) = lambda;
            %brown_consumption_optimal(t,rt) = brown_ratio(t,:,rt)*m;           
        end
    end   
end

fid = fopen('alg1.txt','r');
[objective_rr, count] = fread(fid,'double');
fclose(fid);

fid = fopen('alg2.txt','r');
[objective_gp, count] = fread(fid,'double');
fclose(fid);

s = 1:1:510;
figure;
plot(s, objective_rr(:,1), 'k-', s, objective_gp(:,1),'b-', s, objective(:,1), 'r-', s, optimal_t(:,1), 'k--');
legend('Algorithm 1','Algorithm 2','Algorithm 3','Optimal');

%title(sprintf('total cost for PMR=%.2f (original PMR=%.2f)',newPMR,originalPMR));
xlabel('time(slot)');
ylabel('cost');
ymax=max((objective_rr));
ymin = 0;
ylim([0,ymax]);
xlim([0,510]);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 3.6 2.8]);
print ('-depsc', '../figs/convergence.eps');
eps2pdf('../figs/convergence.eps','/usr/local/bin/gs');

cc=hsv(12);
pp = 1;
figure;
for i = 1:1:a
    hold on;
    plot(1:slot, reshape(lambda_record(i,pp,:,rt),1,slot),'color',cc(i,:));
end
box;
xlim([1,510]);
ylim([0,100]);
set (gcf, 'PaperUnits', 'inches', 'PaperPosition', [0.1 0 3.6 2.8]);
print ('-depsc', '../figs/stability3.eps');
eps2pdf('../figs/stability3.eps','/usr/local/bin/gs');


save('gd.mat');