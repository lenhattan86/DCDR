
TEST = 1;
TEST_CASE = 0;

global N J mu M beta pi_ij
if TEST
    clear global variable; clc; close all;
    DEBUG = 1;
    beta = 1/2;
    N = 2;
    J = 4;
    L_r = 40*ones(J,1);
    w_r = 0*ones(N,1);
    p_r = 5*rand(N,1)+1;
    q_l = 0*rand(N,1);
    M = 5*ones(N,1);
    mu = 10*ones(N,1);
    pi_ij = 2.5*rand(N,J)+1;          
end 


%[x cvx_optval delay] = hetero_opt(x0, L_r, mu, p_l, delay_cost, beta, pi_ij, M, w_r)

[N J] = size(pi_ij);
w = 1; %size(w_r,1);
w_r = reshape(w_r,1,N);
q_l = reshape(q_l,1,N);
delay_cost = beta*ones(N,1);
x0 = zeros(1,N);
%trick, duplicating mu to match the shape of ld_TSJ so that we can use ./
mumu=repmat(reshape(mu,1,N,1),[w,1,J]);
dd=repmat(reshape(pi_ij,1,N,J),[w,1,1])...
    .* repmat(reshape(delay_cost,1,N,1),[w,1,J]);
capcap=repmat(reshape(M,1,N,1),[w,1,1]);

cvx_precision default;
cvx_begin 
    variables x(w,N) ld_TSJ(w,S,J);
    ld=sum(ld_TSJ./mumu,3);
    size(sum(ld_TSJ,1));
    if w > 1
       minimize(sum((ld+quad_over_lin(ld,x-ld,0))*delay_cost + pos(x-w_r)*p_r)...
          + sum(sum(sum(ld_TSJ.*dd))) ...  %adding propagation dealy
	      + sum( pos(x(2:w,:)-x(1:(w-1),:))*beta) + pos(x(1,:)-x0)*beta);
    else
       minimize(sum((ld+quad_over_lin(ld,x-ld,0))*delay_cost + pos(x-w_r-q_l)*p_r)...
          + sum(sum(sum(ld_TSJ.*dd))) ...  %adding propagation dealy
	      + 0);
    end
    subject to
	%squeeze(sum(ld_TSJ,2))==L_r;	% fails for w=1
	sum(ld_TSJ,2)==reshape(L_r, [w, 1, J]);
        ld_TSJ>=0;
        x<=capcap;
        %x>=capcap/100;
        x>=0;
cvx_end
delay = (sum(sum(sum(ld_TSJ.*dd)))+sum((ld+quad_over_lin(ld,x-ld,0))*delay_cost))/(sum(sum(sum(ld_TSJ))));
if (strcmp(cvx_status,'Solved')==0 && strcmp(cvx_status,'Inaccurate/Solved')==0)
    cvx_status
end

