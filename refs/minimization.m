% x = [m(N) lambda(NxJ)]
%global N J mu M beta pi_ij
TEST = 1;
TEST_CASE = 4;

%global N J mu M beta pi_ij
addpath('../functions');
addpath('../classes');    
addpath('../');  
if TEST
    clear global variable; clc; close all;
    DEBUG = 1;
    beta = 2;
    N = 2;
    J = 4;
    L_r = 40*ones(J,1);
    w_r = 0*ones(N,1);
    p_r = 5*rand(N,1)+1;
    q_l = 0*rand(N,1);
    M = 5*ones(N,1);
    mu = 10*ones(N,1);
    pi_ij = 2.5*rand(N,J)+1;          
elseif TEST_CASE > 0
    %x0=[lambda m];
    clear global variable; clc; close all;
    DEBUG = 1;
    LoadTestCase
end 

%lambda_avg = sum(L_r)/(N*J);
%x0 = 0.001*ones(N,J+1); % m + lambda

fun = @(x)p_r'*pos(x(1:N,J+1)-w_r-q_l) ...
    + beta*(sum(sum(x(:,1:J).*pi_ij))) ...
    + beta*sum(x(:,1:J)*ones(J,1).*1./(mu - x(:,1:J)*ones(J,1).*inv_pos(x(1:N,J+1))));
%test simple function vs cvx: 
%     fun = @(x)p_r'*pos(x(1:N,J+1)-q_l-w_r) + beta*(sum(sum(x(:,1:J).*pi_ij))) + beta * sum(inv_pos(mu - inv_pos(x(1:N,J+1))))

lb = zeros(N,J+1); % x>=lb
ub = ones(N,J+1); ub(:,J+1) = M; ub(:,1:J) = inf; %x<=ub
A = [];
b = [];

Aeq = zeros(J,N*(J+1)); 
for idx =1:J
    Aeq(idx,(idx-1)*N+1:idx*N) = ones(N,1)';
end
beq = L_r;
%Aeq = []; beq=[];
% round(Aeq*reshape(x,N*(J+1),1)) == beq

nonlcon = @(x)glb_nonlin_constraint(x,N,J,mu);
if DEBUG
    display = 'iter';
else
    display = 'off';
end

options = optimoptions('fmincon','Display',display,'Algorithm','interior-point', ...
    'MaxIter', 10000, 'MaxFunEvals', 100000, ...
    'TolFun',1.0e-02,'TolX',1.0e-03,...
    'DiffMinChange',1);

xInit = ones(N,J+1); % m + lambda
f = zeros(size(xInit));
x = linprog(f, A, b, Aeq, beq, lb, ub);
f_zero = @(x)0
[x,optVal,isSuccessful] = fmincon(f_zero, xInit, A, b, Aeq, beq, lb, ub, nonlcon, options);

%[x,fval,exitflag] = fmincon(fun, x0, A, b, Aeq, beq, lb, ub, nonlcon);
[x,optVal,isSuccessful] = fmincon(fun, x0, A, b, Aeq, beq, lb, ub, nonlcon,options);
%x=x0;
%x_vector=reshape(x,N*(J+1),1);
m = x(:,J+1);
lambda=x(:,1:J);